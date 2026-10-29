/*
 * E4 stationarity-windows measurement - ESP32-S3 + TMC2209.
 *
 * Measurement only: homing and end-to-end warm-up must already be complete.
 * The stage must be physically parked at x_ref=61.35 mm and the IDS recording
 * before the operator presses 'p' to arm the coordinate, then 'r' to run.
 *
 * STEP is produced by the ESP32 RMT peripheral with the same buffered
 * trapezoidal-profile implementation used by the proven 2/5/8 mm/s timing
 * test. No UART access or serial output occurs inside a pulse train. MSCNT and
 * TMC status are read only after each record window.
 *
 * Wiring matches the ISO and RMT bandwidth-test sketches:
 *   STEP 7, DIR 8, ENN 9, TMC UART TX 17/RX 18, end switches 10/11.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/rmt_tx.h"
#include "e4_config.h"

// Pins.
constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t SW_NEG_PIN = 10;
constexpr uint8_t SW_POS_PIN = 11;
constexpr int SWITCH_ACTIVE_LEVEL = LOW;

// Driver and geometry.
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr double USTEP_PER_MM =
    (double)CFG_STEPS_PER_REV * (double)CFG_MRES / CFG_LEAD_MM;
constexpr uint16_t STEP_HIGH_US = 3;
constexpr uint32_t MAX_RMT_PERIOD_US = 30000;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);
rmt_channel_handle_t rmtChannel = nullptr;
rmt_encoder_handle_t rmtEncoder = nullptr;
rmt_symbol_word_t rmtBuffers[CFG_RMT_RING][CFG_RMT_CHUNK];
size_t rmtBufferIndex = 0;

enum State : uint8_t {
  ST_IDLE, ST_HOME, ST_INIT, ST_DATUM, ST_W1, ST_REF, ST_W2, ST_W3, ST_W4,
  ST_W5, ST_W6, ST_W7, ST_END_REF, ST_PARK, ST_DONE, ST_ABORT
};
const char *STATE_NAME[] = {
  "IDLE", "HOME", "INIT", "DATUM", "W1", "REF", "W2", "W3", "W4",
  "W5", "W6", "W7", "END_REF", "PARK", "DONE", "ABORT"
};

State state = ST_IDLE;
int32_t positionU = CFG_XREF_U;
bool driverReady = false;
bool rmtReady = false;
bool positionArmed = false;
bool abortRequested = false;
bool unsafeFault = false;
bool running = false;
bool autoRun = false;  // 'g': operator gates continue by themselves
bool timeLimited = false;
bool lastDirPositive = true;
bool directionKnown = false;
uint32_t sequenceStartedMs = 0;
uint32_t measuredStops = 0;
uint32_t recordedStops = 0;
uint8_t windowsComplete = 0;

inline double mmOf(int32_t u) { return (double)u / USTEP_PER_MM; }

bool switchPressed(uint8_t pin) {
  if (digitalRead(pin) != SWITCH_ACTIVE_LEVEL) return false;
  delayMicroseconds(300);
  return digitalRead(pin) == SWITCH_ACTIVE_LEVEL;
}

/* Only the switch in the direction of travel can stop a move.
 *
 * Testing both refuses every move while either is pressed, which leaves no
 * way to drive clear of one without reflashing. That is the defect found on
 * the 2/5/8 mm/s run, where it blocked the homing back-offs and the first
 * warm-up move. These sequences never approach a switch - the directional
 * soft-limit guard below stops them first - so this only matters for getting
 * out of a fault, but it costs nothing to get right.
 */
bool switchAhead(bool positive) {
  return switchPressed(positive ? SW_POS_PIN : SW_NEG_PIN);
}

void setState(State next) {
  state = next;
  Serial.printf("# EVENT,%lu,%s,state_enter,x_usteps=%ld\n",
                (unsigned long)millis(), STATE_NAME[state], (long)positionU);
}

void syncWrite(bool high) {
  if (CFG_SYNC_GPIO >= 0) digitalWrite((uint8_t)CFG_SYNC_GPIO, high ? HIGH : LOW);
}

bool pollAbort() {
  while (Serial.available()) {
    const char c = (char)Serial.read();
    if (c == 3) abortRequested = true;
  }
  return abortRequested;
}

bool waitMs(uint32_t durationMs) {
  const uint32_t until = millis() + durationMs;
  while ((int32_t)(millis() - until) < 0) {
    if (pollAbort()) return false;
    delay(2);
  }
  return true;
}

// RMT pulse generator copied from the 2/5/8 mm/s timing-test command path.
bool beginRmt() {
  rmt_tx_channel_config_t channel = {};
  channel.gpio_num = (gpio_num_t)STEP_PIN;
  channel.clk_src = RMT_CLK_SRC_DEFAULT;
  channel.resolution_hz = CFG_RMT_HZ;
  channel.mem_block_symbols = 64;
  channel.trans_queue_depth = CFG_RMT_QUEUE_DEPTH;
  channel.flags.invert_out = false;
  channel.flags.with_dma = false;
  if (rmt_new_tx_channel(&channel, &rmtChannel) != ESP_OK) return false;
  rmt_copy_encoder_config_t encoder = {};
  if (rmt_new_copy_encoder(&encoder, &rmtEncoder) != ESP_OK) return false;
  return rmt_enable(rmtChannel) == ESP_OK;
}

inline void makeRmtSymbol(rmt_symbol_word_t &symbol, uint32_t periodUs) {
  if (periodUs <= STEP_HIGH_US) periodUs = STEP_HIGH_US + 1;
  if (periodUs > MAX_RMT_PERIOD_US) periodUs = MAX_RMT_PERIOD_US;
  symbol.level0 = 1;
  symbol.duration0 = STEP_HIGH_US;
  symbol.level1 = 0;
  symbol.duration1 = periodUs - STEP_HIGH_US;
}

void setDirection(bool positive) {
  if (directionKnown && positive == lastDirPositive) return;
  digitalWrite(DIR_PIN, positive ? HIGH : LOW);
  delayMicroseconds(20);
  lastDirPositive = positive;
  directionKnown = true;
}

uint32_t emitPulses(uint32_t pulses, bool positive,
                    double vMmS = CFG_V_MM_S,
                    uint32_t chunk = CFG_RMT_CHUNK,
                    bool seeking = false) {
  if (!pulses) return 0;
  setDirection(positive);
  const double vMax = vMmS * USTEP_PER_MM;
  const double accel = CFG_A_MM_S2 * USTEP_PER_MM;
  uint32_t nAccel = (uint32_t)(vMax * vMax / (2.0 * accel));
  if (2 * nAccel > pulses) nAccel = pulses / 2;

  uint32_t queued = 0;
  while (queued < pulses) {
    if (pollAbort()) break;
    // Homing seeks the switch on purpose; anywhere else it is a fault.
    if (switchAhead(positive)) {
      if (!seeking) unsafeFault = true;
      break;
    }
    const uint32_t count = min(chunk, pulses - queued);
    const int32_t next = positionU +
        (positive ? (int32_t)count : -(int32_t)count);
    if (!seeking &&
        ((positive && next > CFG_SOFT_MAX_U) || (!positive && next < CFG_SOFT_MIN_U))) {
      unsafeFault = true;
      break;
    }
    rmt_symbol_word_t *buffer = rmtBuffers[rmtBufferIndex];
    rmtBufferIndex = (rmtBufferIndex + 1) % CFG_RMT_RING;
    for (uint32_t k = 0; k < count; ++k) {
      const uint32_t i = queued + k;
      double rate;
      if (i < nAccel) rate = sqrt(2.0 * accel * (double)(i + 1));
      else if (i >= pulses - nAccel) rate = sqrt(2.0 * accel * (double)(pulses - i));
      else rate = vMax;
      if (rate < 1.0) rate = 1.0;
      makeRmtSymbol(buffer[k], (uint32_t)llround(1000000.0 / rate));
    }
    rmt_transmit_config_t tx = {};
    tx.loop_count = 0;
    if (rmt_transmit(rmtChannel, rmtEncoder, buffer,
                     count * sizeof(rmt_symbol_word_t), &tx) != ESP_OK) {
      unsafeFault = true;
      break;
    }
    queued += count;
    positionU += positive ? (int32_t)count : -(int32_t)count;
  }
  if (rmt_tx_wait_all_done(rmtChannel, -1) != ESP_OK) unsafeFault = true;
  return queued;
}

bool moveTo(int32_t targetU, const char *label, uint32_t *tStart = nullptr,
            uint32_t *tArrive = nullptr) {
  if (!positionArmed || targetU < CFG_SOFT_MIN_U || targetU > CFG_SOFT_MAX_U) {
    unsafeFault = true;
    Serial.printf("# ERR,reason=target_or_position,target_usteps=%ld\n", (long)targetU);
    return false;
  }
  const int32_t delta = targetU - positionU;
  const bool positive = delta >= 0;
  const uint32_t wanted = (uint32_t)(positive ? delta : -delta);
  const uint32_t started = millis();
  if (tStart) *tStart = started;
  Serial.printf("# EVENT,%lu,%s,move_start,label=%s,target_usteps=%ld\n",
                (unsigned long)started, STATE_NAME[state], label, (long)targetU);
  const uint32_t emitted = emitPulses(wanted, positive);
  const uint32_t arrived = millis();
  if (tArrive) *tArrive = arrived;
  Serial.printf("# EVENT,%lu,%s,arrive,label=%s,target_usteps=%ld,emitted=%lu\n",
                (unsigned long)arrived, STATE_NAME[state], label, (long)targetU,
                (unsigned long)emitted);
  if (emitted != wanted) {
    Serial.printf("# ERR,reason=incomplete_move,wanted=%lu,emitted=%lu\n",
                  (unsigned long)wanted, (unsigned long)emitted);
    return false;
  }
  return true;
}

/* Homing, as in the step_timing_test routine that homes the stage before the
 * warm-up: seek each switch at 2 mm/s, then three 0.15 mm/s touches that must
 * agree within 0.02 mm. The datum is the negative switch's trigger, which is
 * the frame the generated targets and soft limits are expressed in.
 */
constexpr double HOME_SEEK_MM_S = 2.0;
constexpr double HOME_SLOW_MM_S = 0.15;
constexpr double HOME_BACKOFF_MM = 2.0;
constexpr uint8_t HOME_TOUCHES = 3;
constexpr double HOME_TOL_MM = 0.02;
constexpr double HOME_MAX_TRAVEL_MM = 200.0;
constexpr double END_MARGIN_MM = 2.0;
constexpr uint32_t HOME_CHUNK = 16;  // finer stop granularity on a switch

inline int32_t ustepsOf(double mm) { return (int32_t)llround(mm * USTEP_PER_MM); }

// A relative homing move: the switch ahead ends it, the soft limits do not.
void jog(double mm, double vMmS) {
  const int32_t delta = ustepsOf(mm);
  emitPulses((uint32_t)(delta >= 0 ? delta : -delta), delta >= 0, vMmS,
             HOME_CHUNK, true);
}

bool homeEnd(bool positive, int32_t *trigger) {
  const double sign = positive ? 1.0 : -1.0;
  if (switchAhead(positive)) jog(-sign * HOME_BACKOFF_MM, HOME_SLOW_MM_S);
  jog(sign * HOME_MAX_TRAVEL_MM, HOME_SEEK_MM_S);
  if (abortRequested || !switchAhead(positive)) return false;
  int32_t lo = 0, hi = 0;
  for (uint8_t t = 0; t < HOME_TOUCHES; ++t) {
    jog(-sign * HOME_BACKOFF_MM, HOME_SLOW_MM_S);
    jog(sign * 3.0 * HOME_BACKOFF_MM, HOME_SLOW_MM_S);
    if (abortRequested || !switchAhead(positive)) return false;
    if (t == 0 || positionU < lo) lo = positionU;
    if (t == 0 || positionU > hi) hi = positionU;
  }
  Serial.printf("# HOME_END,%s,spread_mm=%.4f\n", positive ? "pos" : "neg",
                mmOf(hi - lo));
  if (mmOf(hi - lo) > HOME_TOL_MM) return false;
  *trigger = positive ? lo : hi;
  return true;
}

/* Home both ends, check the travel holds the soft limits, then centre at x_ref
 * from below - the same approach as every reference visit. Leaves the
 * coordinate armed, replacing the manual 'p' precondition. */
bool homeAndCentre() {
  positionArmed = false;
  abortRequested = unsafeFault = false;
  setState(ST_HOME);
  digitalWrite(EN_PIN, LOW);
  int32_t negT = 0, posT = 0;
  if (!homeEnd(false, &negT)) {
    Serial.println("# ERR,reason=home_neg");
    return false;
  }
  jog(HOME_BACKOFF_MM, HOME_SLOW_MM_S);
  if (!homeEnd(true, &posT)) {
    Serial.println("# ERR,reason=home_pos");
    return false;
  }
  positionU -= negT;
  const int32_t travel = posT - negT;
  Serial.printf("# HOMED,travel_mm=%.4f,soft_mm=%.3f..%.3f\n", mmOf(travel),
                mmOf(CFG_SOFT_MIN_U), mmOf(CFG_SOFT_MAX_U));
  if (CFG_SOFT_MIN_U < ustepsOf(END_MARGIN_MM) ||
      CFG_SOFT_MAX_U > travel - ustepsOf(END_MARGIN_MM)) {
    Serial.println("# ERR,reason=travel_short_for_soft_limits");
    return false;
  }
  positionArmed = true;
  if (!moveTo(CFG_REF_PRE_U, "centre_pre") || !moveTo(CFG_XREF_U, "centre_xref"))
    return false;
  Serial.printf("# CENTRED,x_usteps=%ld,x_mm=%.6f\n", (long)positionU,
                mmOf(positionU));
  return true;
}

bool driverHealthy() {
  const uint32_t status = driver.DRV_STATUS();
  const bool fault = status & ((1UL << 1) | (1UL << 2) | (1UL << 3) |
                               (1UL << 4) | (1UL << 5));
  if (status & ((1UL << 0) | (1UL << 6) | (1UL << 7)))
    Serial.printf("# WARN,driver_status=0x%08lX\n", (unsigned long)status);
  if (fault) {
    Serial.printf("# ERR,reason=driver_fault,status=0x%08lX\n", (unsigned long)status);
    unsafeFault = true;
  }
  return !fault;
}

bool configureDriver() {
  digitalWrite(EN_PIN, HIGH);
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(100);
  driver.begin();
  driver.toff(5);
  driver.blank_time(24);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.rms_current(CFG_RUN_CURRENT_MA);
  driver.ihold(driver.irun());
  driver.iholddelay(CFG_IHOLDDELAY);
  driver.TPOWERDOWN(CFG_TPOWERDOWN);
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);
  driver.semin(0);
  driver.TCOOLTHRS(0);
  driver.pwm_autoscale(false);
  uint8_t code = 8;
  for (uint16_t value = CFG_MRES; value > 1; value >>= 1) --code;
  const uint32_t mask = 0x0F000000UL;
  driver.CHOPCONF((driver.CHOPCONF() & ~mask) | ((uint32_t)code << 24));
  uint8_t connection = 2;
  for (uint8_t attempt = 0; attempt < 5 && connection; ++attempt) {
    connection = driver.test_connection();
    if (connection) delay(50);
  }
  if (connection) return false;
  const uint32_t gconf = driver.GCONF();
  const uint32_t chop = driver.CHOPCONF();
  const uint32_t status = driver.DRV_STATUS();
  const bool readback = (gconf & (1UL << 2)) && !(status & (1UL << 30)) &&
                        !(chop & (1UL << 28)) && (((chop & mask) >> 24) == code) &&
                        driver.ihold() == driver.irun();
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,status=0x%08lX,"
                "mres=%u,irun=%u,ihold=%u,tpowerdown=%u,iholddelay=%u\n",
                (unsigned long)gconf, (unsigned long)chop, (unsigned long)status,
                CFG_MRES, driver.irun(), driver.ihold(), CFG_TPOWERDOWN,
                CFG_IHOLDDELAY);
  if (!readback) return false;
  digitalWrite(EN_PIN, LOW);
  return true;
}

bool recordStop(const char *window, uint8_t passNumber, int8_t direction,
                uint16_t index, int32_t targetU) {
  uint32_t tMove = 0, tArrive = 0;
  if (!moveTo(targetU, window, &tMove, &tArrive)) return false;
  if (!waitMs(CFG_SETTLE_MS)) return false;
  const uint32_t tRecord0 = millis();
  syncWrite(true);
  Serial.printf("# EVENT,%lu,%s,record_start,pass=%u,dir=%d,index=%u,target_usteps=%ld\n",
                (unsigned long)tRecord0, STATE_NAME[state], passNumber,
                direction, index, (long)targetU);
  const bool completed = waitMs(CFG_RECORD_MS);
  syncWrite(false);
  const uint32_t tRecord1 = millis();
  if (!completed) return false;
  Serial.printf("# EVENT,%lu,%s,record_end,pass=%u,dir=%d,index=%u,target_usteps=%ld\n",
                (unsigned long)tRecord1, STATE_NAME[state], passNumber,
                direction, index, (long)targetU);
  const uint32_t mscnt = driver.MSCNT();  // after the record window, never in motion
  const uint16_t phaseBin = (uint16_t)((targetU % CFG_SCREW_REV_U) / CFG_MOTOR_CYCLE_U);
  Serial.printf("# STOP,%s,%s,%u,%d,%u,%ld,%.6f,%u,%lu,%lu,%lu,%lu,%lu\n",
                STATE_NAME[state], window, passNumber, direction, index,
                (long)targetU, mmOf(targetU), phaseBin,
                (unsigned long)tMove, (unsigned long)tArrive,
                (unsigned long)tRecord0, (unsigned long)tRecord1,
                (unsigned long)mscnt);
  ++recordedStops;
  if (window[0] == 'W') ++measuredStops;
  return driverHealthy();
}

bool visitReference(State asState, uint8_t passNumber, uint8_t index) {
  setState(asState);
  return moveTo(CFG_REF_PRE_U, "ref_pre") &&
         recordStop("REF", passNumber, +1, index, CFG_XREF_U);
}

bool runWindow(State windowState, const char *name, const int32_t *targets,
               uint8_t passNumber) {
  setState(windowState);
  if (!moveTo(targets[0] - CFG_OVERRUN_U, "lower_pre")) return false;
  for (uint16_t i = 0; i < CFG_TARGETS_PER_WINDOW; ++i)
    if (!recordStop(name, passNumber, +1, i + 1, targets[i])) return false;
  if (!moveTo(targets[CFG_TARGETS_PER_WINDOW - 1] + CFG_OVERRUN_U,
              "upper_overrun")) return false;
  for (int i = CFG_TARGETS_PER_WINDOW - 1; i >= 0; --i)
    if (!recordStop(name, passNumber, -1, i + 1, targets[i])) return false;
  ++windowsComplete;
  return true;
}

bool runDatum() {
  for (uint8_t visit = 1; visit <= CFG_DATUM_VISITS; ++visit)
    if (!visitReference(ST_DATUM, 0, visit)) return false;
  return true;
}

uint32_t estimateMoveMs(int32_t fromU, int32_t toU) {
  const double distance = fabs((double)(toU - fromU) / USTEP_PER_MM);
  const double rampDistance = CFG_V_MM_S * CFG_V_MM_S / CFG_A_MM_S2;
  const double seconds = distance <= rampDistance
      ? 2.0 * sqrt(distance / CFG_A_MM_S2)
      : 2.0 * CFG_V_MM_S / CFG_A_MM_S2 +
        (distance - rampDistance) / CFG_V_MM_S;
  return (uint32_t)ceil(seconds * 1000.0);
}

uint32_t estimateWindowMs(const int32_t *targets, int32_t startU, int32_t *endU) {
  uint32_t total = estimateMoveMs(startU, targets[0] - CFG_OVERRUN_U);
  int32_t at = targets[0] - CFG_OVERRUN_U;
  for (uint16_t i = 0; i < CFG_TARGETS_PER_WINDOW; ++i) {
    total += estimateMoveMs(at, targets[i]) + CFG_SETTLE_MS + CFG_RECORD_MS;
    at = targets[i];
  }
  const int32_t overrun = targets[CFG_TARGETS_PER_WINDOW - 1] + CFG_OVERRUN_U;
  total += estimateMoveMs(at, overrun);
  at = overrun;
  for (int i = CFG_TARGETS_PER_WINDOW - 1; i >= 0; --i) {
    total += estimateMoveMs(at, targets[i]) + CFG_SETTLE_MS + CFG_RECORD_MS;
    at = targets[i];
  }
  *endU = at;
  return total;
}

uint32_t estimateReferenceMs(int32_t startU, uint8_t visits, int32_t *endU) {
  uint32_t total = 0;
  int32_t at = startU;
  for (uint8_t i = 0; i < visits; ++i) {
    total += estimateMoveMs(at, CFG_REF_PRE_U);
    total += estimateMoveMs(CFG_REF_PRE_U, CFG_XREF_U) + CFG_SETTLE_MS + CFG_RECORD_MS;
    at = CFG_XREF_U;
  }
  *endU = at;
  return total;
}

uint32_t estimateRemainingAfterCheckpoint() {
  uint32_t total = 0;
  int32_t at = CFG_XREF_U;
  for (uint8_t w = CFG_WATCHDOG_AFTER_WINDOW; w < CFG_WINDOW_COUNT; ++w) {
    total += estimateWindowMs(CFG_WINDOWS[w], at, &at);
    if (w + 1 < CFG_WINDOW_COUNT)
      total += estimateReferenceMs(at, 1, &at);
  }
  total += estimateReferenceMs(at, CFG_END_REF_VISITS, &at);
  return total + 10000UL;
}

bool timeWatchdog() {
  const uint32_t elapsed = millis() - sequenceStartedMs;
  const uint32_t remaining = estimateRemainingAfterCheckpoint();
  const uint32_t projected = elapsed + remaining;
  Serial.printf("# TIME_WATCHDOG,elapsed_ms=%lu,remaining_ms=%lu,projected_ms=%lu,limit_ms=%lu\n",
                (unsigned long)elapsed, (unsigned long)remaining,
                (unsigned long)projected, (unsigned long)CFG_WATCHDOG_MS);
  if (projected <= CFG_WATCHDOG_MS) return true;
  timeLimited = true;
  Serial.println("# WARN,time_limit_projected,skip_remaining_windows=1,stationarity_complete=0");
  return false;
}

void printHeader() {
  Serial.printf("# SESSION,experiment=E4_stationarity_windows,config=%s,windows=%u,"
                "targets_per_window=%u,grid_usteps=%u,time_limit_ms=%lu,motion=RMT\n",
                CFG_HASH, CFG_WINDOW_COUNT, CFG_TARGETS_PER_WINDOW, CFG_GRID_SPACING_U,
                (unsigned long)CFG_TIME_LIMIT_MS);
  Serial.printf("# SETTINGS,lead_mm=%.3f,steps_rev=%u,mres=%u,chopper=SpreadCycle,"
                "run_ma=%u,hold_ma=%u,tpowerdown=%u,iholddelay=%u,v_mm_s=%.3f,"
                "a_mm_s2=%.3f,settle_ms=%lu,record_ms=%lu,overrun_usteps=%ld,"
                "xref_usteps=%ld,rmt_hz=%lu,rmt_chunk=%u\n",
                CFG_LEAD_MM, CFG_STEPS_PER_REV, CFG_MRES, CFG_RUN_CURRENT_MA,
                CFG_HOLD_CURRENT_MA, CFG_TPOWERDOWN, CFG_IHOLDDELAY, CFG_V_MM_S,
                CFG_A_MM_S2, (unsigned long)CFG_SETTLE_MS,
                (unsigned long)CFG_RECORD_MS, (long)CFG_OVERRUN_U,
                (long)CFG_XREF_U, (unsigned long)CFG_RMT_HZ, CFG_RMT_CHUNK);
  Serial.printf("# SETTINGS2,nominal_stop_ms=%lu,aquadb_clock_ns=%u,"
                "soft_min_usteps=%ld,soft_max_usteps=%ld,ref_pre_usteps=%ld,"
                "screw_rev_usteps=%u,ref_drift_warning_um=%.3f,sync_gpio=%d\n",
                (unsigned long)CFG_NOMINAL_STOP_MS,
                CFG_AQUADB_CLOCK_NS, (long)CFG_SOFT_MIN_U, (long)CFG_SOFT_MAX_U,
                (long)CFG_REF_PRE_U, CFG_SCREW_REV_U,
                CFG_REF_DRIFT_WARNING_UM, CFG_SYNC_GPIO);
  for (uint8_t w = 0; w < CFG_WINDOW_COUNT; ++w)
    for (uint16_t i = 0; i < CFG_TARGETS_PER_WINDOW; ++i)
      Serial.printf("# TARGET,W%u,%u,%ld,%.6f,%u\n", w + 1, i + 1,
                    (long)CFG_WINDOWS[w][i], mmOf(CFG_WINDOWS[w][i]),
                    (unsigned)((CFG_WINDOWS[w][i] % CFG_SCREW_REV_U) / CFG_MOTOR_CYCLE_U));
}

void parkIfSafe() {
  setState(ST_PARK);
  if (positionArmed && !unsafeFault) {
    abortRequested = false;
    moveTo(CFG_REF_PRE_U, "park_pre");
    moveTo(CFG_XREF_U, "park_xref");
  } else {
    digitalWrite(EN_PIN, HIGH);
    positionArmed = false;
    Serial.println("# PARK_SKIPPED,position_not_safe=1,motor_disabled=1");
  }
}

void runE4() {
  if (!driverReady || !rmtReady || !positionArmed || running) {
    Serial.printf("# ERR,reason=not_ready,driver=%u,rmt=%u,armed=%u,running=%u\n",
                  driverReady, rmtReady, positionArmed, running);
    return;
  }
  running = true;
  abortRequested = unsafeFault = timeLimited = false;
  measuredStops = 0;
  recordedStops = 0;
  windowsComplete = 0;
  sequenceStartedMs = millis();
  setState(ST_INIT);
  printHeader();

  bool ok = driverHealthy() && !switchPressed(SW_NEG_PIN) && !switchPressed(SW_POS_PIN);
  if (ok) ok = runDatum();

  const State windowStates[CFG_WINDOW_COUNT] = {
    ST_W1, ST_W2, ST_W3, ST_W4, ST_W5, ST_W6, ST_W7
  };
  uint8_t refIndex = 0;
  for (uint8_t w = 0; ok && w < CFG_WINDOW_COUNT; ++w) {
    char windowName[3];
    snprintf(windowName, sizeof(windowName), "W%u", unsigned(w + 1));
    ok = runWindow(windowStates[w], windowName, CFG_WINDOWS[w], 1);
    if (ok && (w + 1) < CFG_WINDOW_COUNT)
      ok = visitReference(ST_REF, 1, ++refIndex);
    if (ok && (w + 1) == CFG_WATCHDOG_AFTER_WINDOW && !timeWatchdog())
      break;
  }

  if (ok || timeLimited) {
    for (uint8_t visit = 1; visit <= CFG_END_REF_VISITS; ++visit)
      if (!visitReference(ST_END_REF, 1, ++refIndex)) {
        ok = false;
        break;
      }
    Serial.printf("# ENV_SNAPSHOT_REQUEST,%lu\n", (unsigned long)millis());
  }

  if (!ok && !timeLimited) setState(ST_ABORT);
  parkIfSafe();
  const bool complete = ok && !timeLimited &&
                        windowsComplete == CFG_WINDOW_COUNT &&
                        measuredStops == CFG_ANALYSIS_STOP_COUNT &&
                        recordedStops == CFG_TOTAL_RECORD_STOP_COUNT;
  setState(complete ? ST_DONE : ST_ABORT);
  Serial.printf("# SUMMARY,elapsed_ms=%lu,analysis_stops=%lu,recorded_stops=%lu,windows_complete=%u,"
                "complete=%u,time_limited=%u,aborted=%u,unsafe=%u,x_usteps=%ld\n",
                (unsigned long)(millis() - sequenceStartedMs),
                (unsigned long)measuredStops, (unsigned long)recordedStops,
                unsigned(windowsComplete), complete,
                timeLimited, (!ok && !timeLimited), unsafeFault, (long)positionU);
  running = false;
}

void printMenu() {
  Serial.println("\n--- E4 stationarity windows ---");
  Serial.printf("driver=%u rmt=%u armed_at_xref=%u position_mm=%.6f\n",
                driverReady, rmtReady, positionArmed, mmOf(positionU));
  Serial.println("p  accept precondition and arm at x_ref");
  Serial.println("r  run E4 measurement    t  print target summary");
  Serial.println("s  status                0  release motor");
  Serial.println("h  home both ends, centre at x_ref and arm");
  Serial.println("g  home, centre, then run without operator gates");
  Serial.println("Ctrl+C requests a controlled abort");
  Serial.print("> ");
}

void setup() {
  pinMode(EN_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(DIR_PIN, HIGH);
  pinMode(SW_NEG_PIN, INPUT_PULLUP);
  pinMode(SW_POS_PIN, INPUT_PULLUP);
  if (CFG_SYNC_GPIO >= 0) {
    pinMode((uint8_t)CFG_SYNC_GPIO, OUTPUT);
    syncWrite(false);
  }
  Serial.begin(115200);
  const uint32_t started = millis();
  while (!Serial && millis() - started < 3000) delay(10);
  delay(200);
  Serial.printf("# BOOT,fw=e4_stationarity_windows,config=%s\n", CFG_HASH);
  rmtReady = beginRmt();
  driverReady = configureDriver();
  Serial.printf("# READY,driver=%u,rmt=%u,precondition_armed=0\n",
                driverReady, rmtReady);
  Serial.println("# PRECONDITION_REQUIRED,homed=1,warm=1,xref_mm=61.350,ids_recording=1");
  printMenu();
}

void loop() {
  if (!Serial.available()) { delay(5); return; }
  const char key = (char)Serial.read();
  if (key == '\r' || key == '\n') return;
  if (key == 3) {
    abortRequested = true;
    Serial.println("# ABORT_REQUEST");
    return;
  }
  switch (key) {
    case 'p':
      if (!driverReady || !rmtReady || switchPressed(SW_NEG_PIN) || switchPressed(SW_POS_PIN)) {
        Serial.println("# ERR,reason=cannot_arm");
      } else {
        digitalWrite(EN_PIN, LOW);
        positionU = CFG_XREF_U;
        positionArmed = true;
        abortRequested = false;
        unsafeFault = false;
        Serial.printf("# PRECONDITION_ARMED,x_usteps=%ld,x_mm=%.6f\n",
                      (long)positionU, mmOf(positionU));
      }
      break;
    case 'r': runE4(); break;
    case 'h': homeAndCentre(); break;
    case 'g':
      if (homeAndCentre()) {
        autoRun = true;
        runE4();
        autoRun = false;
      }
      break;
    case 's':
      Serial.printf("# STATUS,state=%s,armed=%u,x_usteps=%ld,x_mm=%.6f,"
                    "neg=%u,pos=%u,driver_ok=%u\n",
                    STATE_NAME[state], positionArmed, (long)positionU, mmOf(positionU),
                    switchPressed(SW_NEG_PIN), switchPressed(SW_POS_PIN),
                    driverHealthy());
      break;
    case 't':
      Serial.printf("# PLAN,windows=%u,targets_per_window=%u,analysis_stops=%lu,grid_usteps=%u\n",
                    CFG_WINDOW_COUNT, CFG_TARGETS_PER_WINDOW,
                    (unsigned long)CFG_ANALYSIS_STOP_COUNT, CFG_GRID_SPACING_U);
      for (uint8_t w = 0; w < CFG_WINDOW_COUNT; ++w)
        Serial.printf("# WINDOW,index=%u,start_mm=%.3f,centre_mm=%.3f,end_mm=%.3f\n",
                      unsigned(w + 1), mmOf(CFG_WINDOWS[w][0]),
                      mmOf(CFG_WINDOWS[w][CFG_TARGETS_PER_WINDOW / 2]),
                      mmOf(CFG_WINDOWS[w][CFG_TARGETS_PER_WINDOW - 1]));
      break;
    case '0':
      digitalWrite(EN_PIN, HIGH);
      positionArmed = false;
      Serial.println("# MOTOR_RELEASED,coordinate_invalid=1");
      break;
    case '?': break;
    default: Serial.println("# ERR,reason=unknown_command"); break;
  }
  printMenu();
}

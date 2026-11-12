/*
 * F0b_pt2 - reversal vs microstep resolution, and speed.  ESP32-S3 + TMC2209.
 *
 * The host-generated f0b_pt2_config.h holds every integer command in
 * execution order.  Nothing is computed here but the step timing.  This is
 * F0b's motion engine with the four changes procedure section 0 asks for.
 *
 *   1  ONE INTERNAL UNIT.  positionU is always in 1/256 counts, whatever MRES
 *      is in force (1 mm = 25600).  A segment's travel is converted to driver
 *      steps by dividing by uStepPerDriverStep = 1 << mres_code, and a
 *      segment whose travel does not divide is a fault, not a rounding.
 *
 *   2  MRES SWITCH SEGMENT, at rest only.  CHOPCONF.MRES is written over
 *      UART, CHOPCONF and MSCNT are read back before and after, both are
 *      logged, and the run aborts if MSCNT moved or the field did not take.
 *      MSCNT counts the microstep table position 0..1023 at every MRES, so a
 *      change of MRES must leave it alone; if it does not, the coil current
 *      vector moved and the comparison between levels is void.
 *
 *   3  WINDOWS TIMED FROM THE LAST STEP EDGE.  lastEdgeUs is stamped when the
 *      final pulse of a motion segment has actually been transmitted.  A
 *      DWELL ends at lastEdgeUs + end_offset_ms, a MARK at
 *      lastEdgeUs + end_offset_ms, with the offsets generated as
 *      settle and settle + record.  The settle therefore ends 2.000 s after
 *      the stage stopped moving and the mark 2.500 s after it, whatever the
 *      scheduler did in between, and the analysis can read 0.3-0.5 s and
 *      2.0-2.5 s from one dwell.
 *
 *   4  Telemetry buffered in RAM, printed after SEQ_DONE, as before.
 *
 * Two step paths, unchanged from F0b: at or above ~33 Hz the buffered RMT
 * path, below it one pulse per micros() deadline.  Block V's slowest visit is
 * 0.01 mm/s, which at 1/256 is 256 Hz, so every measurement segment in this
 * session uses the RMT path; the slow path is retained for homing and for
 * the fallback configurations.
 *
 * Before flashing, perform the E0 warm-up and the bench check that the RMT
 * path produces 51.2 kHz cleanly (2 mm/s at 1/256).  Then 'g' to home,
 * centre and run, or 'p'/'r' when already homed and parked at 60.000 mm.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/rmt_tx.h"
#include "f0b_pt2_config.h"

constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t SW_NEG_PIN = 10;
constexpr uint8_t SW_POS_PIN = 11;
constexpr int SWITCH_ACTIVE_LEVEL = LOW;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t STEP_HIGH_US = 3;

constexpr double USTEP_PER_MM = 256.0 * CFG_STEPS_PER_REV / CFG_LEAD_MM;
constexpr double ACCEL_U_S2 = CFG_A_MM_S2 * USTEP_PER_MM;
constexpr double TRANSIT_RATE_U_S = CFG_V_MM_S * USTEP_PER_MM;
constexpr double RMT_MIN_RATE_HZ = 1000000.0 / (double)CFG_RMT_MAX_PERIOD_US;
constexpr uint8_t HOMING_MRES_CODE = 0;          // 1/256 while homing

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);
rmt_channel_handle_t rmtChannel = nullptr;
rmt_encoder_handle_t rmtEncoder = nullptr;
rmt_symbol_word_t rmtBuffers[CFG_RMT_RING][CFG_RMT_CHUNK];
size_t rmtBufferIndex = 0;

struct SegRecord {
  uint16_t row;
  uint32_t startMs;
  uint32_t endMs;
  uint32_t lastEdgeUs;      // when the final pulse of this segment went out
  uint16_t mscntStart;
  uint16_t mscntEnd;
  uint32_t drvStatus;
  uint32_t aux;             // SETMRES: CHOPCONF read back after the write
};

SegRecord records[CFG_SEG_COUNT];
uint16_t recordCount = 0;
int32_t positionU = CFG_PARK_U;
uint32_t uStepPerDriverStep = 1;     // 1/256 counts per driver step
uint8_t currentMresCode = HOMING_MRES_CODE;
uint32_t lastEdgeUs = 0;
uint16_t currentBlock = 0;
bool driverReady = false;
bool rmtReady = false;
bool positionArmed = false;
bool abortRequested = false;
bool unsafeFault = false;
bool running = false;
bool timeLimited = false;
bool mresFault = false;
bool lastDirPositive = true;
bool directionKnown = false;
uint32_t sequenceStartedMs = 0;
uint16_t homeNegMscnt = 0;
uint16_t homePosMscnt = 0;
uint16_t centreMscnt = 0;
bool homePhaseValid = false;

inline double mmOf(int32_t u) { return (double)u / USTEP_PER_MM; }
inline int32_t ustepsOf(double mm) { return (int32_t)llround(mm * USTEP_PER_MM); }
inline uint16_t mresOf(uint8_t code) { return (uint16_t)(256u >> code); }

bool switchPressed(uint8_t pin) {
  if (digitalRead(pin) != SWITCH_ACTIVE_LEVEL) return false;
  delayMicroseconds(300);
  return digitalRead(pin) == SWITCH_ACTIVE_LEVEL;
}

bool switchAhead(bool positive) {
  return switchPressed(positive ? SW_POS_PIN : SW_NEG_PIN);
}

void syncWrite(bool high) {
  if (CFG_SYNC_GPIO >= 0) digitalWrite((uint8_t)CFG_SYNC_GPIO, high ? HIGH : LOW);
}

bool pollAbort() {
  if (Serial.available() && Serial.read() == 3) {
    abortRequested = true;
    return true;
  }
  return abortRequested || unsafeFault;
}

bool waitUntilUs(uint32_t deadlineUs) {
  while ((int32_t)(deadlineUs - micros()) > 0) {
    if (pollAbort()) return false;
    const int32_t remaining = (int32_t)(deadlineUs - micros());
    if (remaining > 4000) delay(1);
  }
  return true;
}

bool waitMs(uint32_t durationMs) {
  return waitUntilUs(micros() + durationMs * 1000UL);
}

bool beginRmt() {
  rmt_tx_channel_config_t channel = {};
  channel.gpio_num = (gpio_num_t)STEP_PIN;
  channel.clk_src = RMT_CLK_SRC_DEFAULT;
  channel.resolution_hz = CFG_RMT_HZ;
  channel.mem_block_symbols = 64;
  channel.trans_queue_depth = 4;
  if (rmt_new_tx_channel(&channel, &rmtChannel) != ESP_OK) return false;
  rmt_copy_encoder_config_t encoder = {};
  if (rmt_new_copy_encoder(&encoder, &rmtEncoder) != ESP_OK) return false;
  return rmt_enable(rmtChannel) == ESP_OK;
}

inline void makeRmtSymbol(rmt_symbol_word_t &symbol, uint32_t periodUs) {
  if (periodUs < STEP_HIGH_US + 1) periodUs = STEP_HIGH_US + 1;
  if (periodUs > CFG_RMT_MAX_PERIOD_US) periodUs = CFG_RMT_MAX_PERIOD_US;
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

// One safety gate for every pulse-emitting path.  count is in driver steps.
bool stepAllowed(bool positive, uint32_t count, bool seeking) {
  if (switchAhead(positive)) {
    if (!seeking) unsafeFault = true;
    return false;
  }
  const int32_t travel = (int32_t)(count * uStepPerDriverStep);
  const int32_t next = positionU + (positive ? travel : -travel);
  // Homing leaves the stage on the positive switch, outside the soft window,
  // so only a step that goes further out is a fault.
  if (!seeking && ((positive && next > CFG_SOFT_MAX_U) ||
                   (!positive && next < CFG_SOFT_MIN_U))) {
    unsafeFault = true;
    return false;
  }
  return true;
}

bool emitSinglePulse(bool positive, bool seeking) {
  if (!stepAllowed(positive, 1, seeking)) return false;
  rmt_symbol_word_t symbol;
  makeRmtSymbol(symbol, STEP_HIGH_US + 1);
  rmt_transmit_config_t tx = {};
  tx.loop_count = 0;
  if (rmt_transmit(rmtChannel, rmtEncoder, &symbol, sizeof(symbol), &tx) != ESP_OK ||
      rmt_tx_wait_all_done(rmtChannel, -1) != ESP_OK) {
    unsafeFault = true;
    return false;
  }
  lastEdgeUs = micros();
  positionU += positive ? (int32_t)uStepPerDriverStep
                        : -(int32_t)uStepPerDriverStep;
  return true;
}

// Slow path: one pulse per micros() deadline, below RMT_MIN_RATE_HZ.
uint32_t emitSlowRate(uint32_t pulses, bool positive, double rateHz) {
  setDirection(positive);
  const uint32_t periodUs = (uint32_t)llround(1000000.0 / rateHz);
  uint32_t deadline = micros();
  uint32_t done = 0;
  while (done < pulses) {
    if (pollAbort()) break;
    if (!emitSinglePulse(positive, false)) break;
    ++done;
    deadline += periodUs;
    if (done < pulses && !waitUntilUs(deadline)) break;
  }
  return done;
}

// Buffered RMT path: trapezoidal profile peaking at rateHz, in driver steps.
uint32_t emitProfile(uint32_t pulses, bool positive, double rateHz,
                     uint32_t chunk = CFG_RMT_CHUNK, bool seeking = false) {
  if (!pulses) return 0;
  setDirection(positive);
  // Acceleration is a property of the stage, so it is expressed in counts and
  // converted here: the same mm/s^2 at every MRES.
  const double accelSteps = ACCEL_U_S2 / (double)uStepPerDriverStep;
  uint32_t nAccel = (uint32_t)(rateHz * rateHz / (2.0 * accelSteps));
  if (2 * nAccel > pulses) nAccel = pulses / 2;
  uint32_t queued = 0;
  while (queued < pulses) {
    if (pollAbort()) break;
    const uint32_t count = min(chunk, pulses - queued);
    if (!stepAllowed(positive, count, seeking)) break;
    rmt_symbol_word_t *buffer = rmtBuffers[rmtBufferIndex];
    rmtBufferIndex = (rmtBufferIndex + 1) % CFG_RMT_RING;
    for (uint32_t k = 0; k < count; ++k) {
      const uint32_t i = queued + k;
      double rate;
      if (i < nAccel) rate = sqrt(2.0 * accelSteps * (double)(i + 1));
      else if (i >= pulses - nAccel)
        rate = sqrt(2.0 * accelSteps * (double)(pulses - i));
      else rate = rateHz;
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
    positionU += positive ? (int32_t)(count * uStepPerDriverStep)
                          : -(int32_t)(count * uStepPerDriverStep);
  }
  if (rmt_tx_wait_all_done(rmtChannel, -1) != ESP_OK) unsafeFault = true;
  lastEdgeUs = micros();
  return queued;
}

uint32_t emitAtRate(uint32_t pulses, bool positive, double rateHz) {
  if (!pulses) return 0;
  if (rateHz >= RMT_MIN_RATE_HZ) return emitProfile(pulses, positive, rateHz);
  return emitSlowRate(pulses, positive, rateHz);
}

// rateUperS is in 1/256 counts per second: the same number at every MRES.
bool travelTo(int32_t targetU, double rateUperS) {
  if (!positionArmed || targetU < CFG_SOFT_MIN_U || targetU > CFG_SOFT_MAX_U) {
    unsafeFault = true;
    return false;
  }
  const int32_t delta = targetU - positionU;
  if (!delta) return true;
  // Change 1: the travel must be a whole number of steps at this MRES.
  if (delta % (int32_t)uStepPerDriverStep) {
    Serial.printf("# ERR,reason=not_divisible,delta_u=%ld,mres=1/%u\n",
                  (long)delta, mresOf(currentMresCode));
    unsafeFault = true;
    return false;
  }
  const bool positive = delta > 0;
  const uint32_t wanted =
      (uint32_t)((positive ? delta : -delta) / (int32_t)uStepPerDriverStep);
  const double rateSteps = rateUperS / (double)uStepPerDriverStep;
  const uint32_t emitted = emitAtRate(wanted, positive, rateSteps);
  return emitted == wanted && !unsafeFault && !abortRequested;
}

// Change 2: MRES is written at rest, and MSCNT must not move across it.
bool applyMres(uint8_t code, uint16_t *before, uint16_t *after,
               uint32_t *chopAfter) {
  *before = (uint16_t)driver.MSCNT();
  const uint32_t chop = driver.CHOPCONF();
  driver.CHOPCONF((chop & ~0x0F000000UL) | ((uint32_t)code << 24));
  delay(5);
  *chopAfter = driver.CHOPCONF();
  *after = (uint16_t)driver.MSCNT();
  const bool took = ((*chopAfter >> 24) & 0x0FUL) == code;
  // At a coarser MRES the TMC2209 only holds MSCNT at half-step positions
  // (step/2 + k*step), and snaps to the nearest one when the field changes:
  // 433 -> 440 at 1/16.  No MSCNT is on both the 1/16 and 1/32 grids, so an
  // unchanged MSCNT cannot be required.  Accept that snap, nothing larger;
  // MSCNT before/after is logged for every switch.  Going finer snaps too
  // (421 -> 420 at 1/128 -> 1/256), so the bound is half the coarser step.
  const uint8_t prev = currentMresCode > 8 ? 0 : currentMresCode;
  const int16_t step = (int16_t)(1 << code);
  const int16_t coarser = (int16_t)(1 << max(code, prev));
  int16_t moved = (int16_t)((*after - *before + 1536) % 1024) - 512;
  const bool onGrid = code == 0 || (*after % step) == step / 2;
  const bool still = onGrid && abs(moved) <= coarser / 2;
  if (took && still) {
    currentMresCode = code;
    uStepPerDriverStep = 1UL << code;
    return true;
  }
  mresFault = true;
  unsafeFault = true;
  return false;
}

// Homing is the proven E0/E1-E4 implementation, run at 1/256 throughout.
constexpr double HOME_SEEK_MM_S = 2.0;
constexpr double HOME_SLOW_MM_S = 0.15;
constexpr double HOME_BACKOFF_MM = 2.0;
constexpr uint8_t HOME_TOUCHES = 3;
constexpr double HOME_TOL_MM = 0.02;
constexpr double HOME_MAX_TRAVEL_MM = 200.0;
constexpr double END_MARGIN_MM = 2.0;
constexpr uint32_t HOME_CHUNK = 16;

void jog(double mm, double vMmS) {
  const int32_t delta = ustepsOf(fabs(mm));
  if (!delta) return;
  emitProfile((uint32_t)(delta / (int32_t)uStepPerDriverStep), mm > 0,
              vMmS * USTEP_PER_MM / (double)uStepPerDriverStep,
              HOME_CHUNK, true);
}

bool homeEnd(bool positive, int32_t *trigger) {
  // One seeking profile that the switch cuts short, as in F0b.  Issuing it
  // as separate 16-step profiles ramps every chunk from rest and seeks at
  // ~0.1 mm/s instead of HOME_SEEK_MM_S.
  if (!switchPressed(positive ? SW_POS_PIN : SW_NEG_PIN))
    jog(positive ? HOME_MAX_TRAVEL_MM : -HOME_MAX_TRAVEL_MM, HOME_SEEK_MM_S);
  if (abortRequested || !switchPressed(positive ? SW_POS_PIN : SW_NEG_PIN))
    return false;
  int32_t sum = 0;
  for (uint8_t touch = 0; touch < HOME_TOUCHES; ++touch) {
    jog(positive ? -HOME_BACKOFF_MM : HOME_BACKOFF_MM, HOME_SEEK_MM_S);
    int32_t crept = 0;
    while (!switchPressed(positive ? SW_POS_PIN : SW_NEG_PIN)) {
      if (pollAbort()) return false;
      emitProfile(1, positive,
                  HOME_SLOW_MM_S * USTEP_PER_MM / (double)uStepPerDriverStep,
                  1, true);
      crept += (int32_t)uStepPerDriverStep;
      if (crept > ustepsOf(2 * HOME_BACKOFF_MM)) return false;
    }
    sum += crept;
  }
  *trigger = sum / HOME_TOUCHES;
  return true;
}

bool homeAndCentre() {
  int32_t negTrigger = 0, posTrigger = 0;
  positionArmed = false;
  if (!homeEnd(false, &negTrigger)) return false;
  homeNegMscnt = (uint16_t)driver.MSCNT();
  const int32_t crossStart = positionU;
  if (!switchPressed(SW_POS_PIN)) jog(HOME_MAX_TRAVEL_MM, HOME_SEEK_MM_S);
  if (abortRequested || !switchPressed(SW_POS_PIN)) return false;
  const int32_t span = positionU - crossStart;
  if (!homeEnd(true, &posTrigger)) return false;
  homePosMscnt = (uint16_t)driver.MSCNT();
  Serial.printf("# HOMED,travel_mm=%.4f,soft_mm=%.3f..%.3f\n",
                mmOf(span), mmOf(CFG_SOFT_MIN_U), mmOf(CFG_SOFT_MAX_U));
  positionU = span;
  positionArmed = true;
  if (!travelTo(CFG_PARK_U, HOME_SEEK_MM_S * USTEP_PER_MM)) {
    Serial.println("# ERR,reason=centre_failed");
    return false;
  }
  centreMscnt = (uint16_t)driver.MSCNT();
  homePhaseValid = true;
  Serial.printf("# CENTRED,x_usteps=%ld,x_mm=%.6f,MSCNT=%u,"
                "home_neg_MSCNT=%u,home_pos_MSCNT=%u\n",
                (long)positionU, mmOf(positionU), centreMscnt,
                homeNegMscnt, homePosMscnt);
  return true;
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
  uint16_t before = 0, after = 0;
  uint32_t chop = 0;
  currentMresCode = 0xFF;                     // force the write
  if (!applyMres(HOMING_MRES_CODE, &before, &after, &chop)) {
    Serial.println("# ERR,reason=mres_init_failed");
    return false;
  }
  const uint32_t status = driver.DRV_STATUS();
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,status=0x%08lX,"
                "mres=%u,irun=%u,ihold=%u,tpowerdown=%u,iholddelay=%u\n",
                (unsigned long)driver.GCONF(), (unsigned long)chop,
                (unsigned long)status, mresOf(currentMresCode),
                driver.irun(), driver.ihold(), CFG_TPOWERDOWN, CFG_IHOLDDELAY);
  digitalWrite(EN_PIN, LOW);
  return (status & 0x0000003FUL) == 0;
}

bool driverStatusFault(uint32_t status) {
  return status & ((1UL << 1) | (1UL << 2) | (1UL << 3) |
                   (1UL << 4) | (1UL << 5));
}

bool executeSegment(const Pt2Seg &seg) {
  if (recordCount >= CFG_SEG_COUNT) {
    unsafeFault = true;
    return false;
  }
  SegRecord &record = records[recordCount];
  record.row = seg.row;
  record.mscntStart = (uint16_t)driver.MSCNT();
  record.startMs = millis();
  record.aux = 0;
  bool ok = true;
  switch (seg.type) {
    case PT2_MOVE:
      ok = travelTo(seg.u_target, TRANSIT_RATE_U_S);
      break;
    case PT2_RAMP:
      ok = travelTo(seg.u_target, (double)seg.rate_mcps / 1000.0);
      break;
    case PT2_DWELL:
      // Change 3: the window ends at a fixed offset from the last step edge.
      ok = waitUntilUs(lastEdgeUs + seg.end_offset_ms * 1000UL);
      break;
    case PT2_MARK:
      syncWrite(true);
      ok = waitUntilUs(lastEdgeUs + seg.end_offset_ms * 1000UL);
      syncWrite(false);
      break;
    case PT2_SETMRES: {
      uint16_t before = 0, after = 0;
      uint32_t chop = 0;
      ok = applyMres(seg.mres_code, &before, &after, &chop);
      record.aux = chop;
      break;
    }
    default:
      unsafeFault = true;
      ok = false;
      break;
  }
  record.endMs = millis();
  record.lastEdgeUs = lastEdgeUs;
  record.mscntEnd = (uint16_t)driver.MSCNT();
  record.drvStatus = driver.DRV_STATUS();
  ++recordCount;
  if (driverStatusFault(record.drvStatus)) {
    unsafeFault = true;
    return false;
  }
  return ok;
}

const char *typeName(uint8_t type) {
  switch (type) {
    case PT2_MOVE: return "move";
    case PT2_RAMP: return "ramp";
    case PT2_DWELL: return "dwell";
    case PT2_MARK: return "mark";
    case PT2_SETMRES: return "setmres";
    default: return "?";
  }
}

void printHeader() {
  Serial.printf("# SESSION,experiment=%s,procedure=pt2,config_sha256=%s,"
                "table_sha256=%s,segments=%u,marks=%u,motion=RMT+timed,"
                "unit=1/256_counts,windows_from=t_last_edge\n",
                CFG_EXPERIMENT, CFG_HASH, CFG_TABLE_HASH, CFG_SEG_COUNT,
                CFG_MARK_COUNT);
  Serial.printf("# SETTINGS,lead_mm=%.3f,steps_rev=%u,start_mres=%u,"
                "chopper=SpreadCycle,run_ma=%u,hold_ma=%u,v_mm_s=%.3f,"
                "a_mm_s2=%.3f,settle_ms=%lu,record_ms=%lu,ustep_per_mm=%.1f,"
                "rmt_min_rate_hz=%.2f\n",
                CFG_LEAD_MM, CFG_STEPS_PER_REV, CFG_START_MRES,
                CFG_RUN_CURRENT_MA, CFG_HOLD_CURRENT_MA, CFG_V_MM_S,
                CFG_A_MM_S2, (unsigned long)CFG_SETTLE_MS,
                (unsigned long)CFG_RECORD_MS, USTEP_PER_MM, RMT_MIN_RATE_HZ);
  Serial.printf("# HOME_PHASE,valid=%u,home_neg_MSCNT=%u,home_pos_MSCNT=%u,"
                "centre_MSCNT=%u,x_usteps=%ld,x_mm=%.6f\n",
                homePhaseValid, homeNegMscnt, homePosMscnt, centreMscnt,
                (long)positionU, mmOf(positionU));
}

void dumpTiming() {
  for (uint16_t i = 0; i < recordCount; ++i) {
    const SegRecord &record = records[i];
    const Pt2Seg &seg = CFG_SEGS[record.row];
    Serial.printf("# SEG,%u,%s,%s,%u,%s,%ld,%lu,%lu,%u,%u,%lu,%lu,%lu,"
                  "0x%08lX,0x%08lX\n",
                  seg.row, CFG_BLOCK_NAMES[seg.block], typeName(seg.type),
                  mresOf(seg.mres_code), CFG_SIDE_NAMES[seg.side],
                  (long)seg.u_target, (unsigned long)seg.rate_mcps,
                  (unsigned long)seg.end_offset_ms,
                  record.mscntStart, record.mscntEnd,
                  (unsigned long)record.startMs, (unsigned long)record.endMs,
                  (unsigned long)record.lastEdgeUs,
                  (unsigned long)record.drvStatus, (unsigned long)record.aux);
    if ((i & 0x1F) == 0x1F) delay(2);
  }
}

void runSession() {
  if (!driverReady || !rmtReady || !positionArmed || running) {
    Serial.printf("# ERR,reason=not_ready,driver=%u,rmt=%u,armed=%u,running=%u\n",
                  driverReady, rmtReady, positionArmed, running);
    return;
  }
  running = true;
  abortRequested = unsafeFault = timeLimited = mresFault = false;
  recordCount = 0;
  sequenceStartedMs = millis();
  lastEdgeUs = micros();
  centreMscnt = (uint16_t)driver.MSCNT();
  printHeader();
  Serial.println("# SEQ_START,uart_quiet=1,timing_buffer=RAM");
  bool ok = true;
  for (uint16_t i = 0; ok && i < CFG_SEG_COUNT; ++i) {
    if (millis() - sequenceStartedMs >= CFG_WATCHDOG_MS) {
      timeLimited = true;
      ok = false;
      break;
    }
    currentBlock = CFG_SEGS[i].block;
    ok = executeSegment(CFG_SEGS[i]);
    if (abortRequested || unsafeFault) ok = false;
  }
  const uint32_t elapsed = millis() - sequenceStartedMs;
  const bool complete = ok && recordCount == CFG_SEG_COUNT;
  Serial.printf("# SEQ_DONE,elapsed_ms=%lu,planned_ms=%lu,segments=%u,"
                "complete=%u,aborted=%u,unsafe=%u,time_limited=%u,"
                "mres_fault=%u,block=%s\n",
                (unsigned long)elapsed, (unsigned long)CFG_PLAN_DURATION_MS,
                recordCount, complete, abortRequested, unsafeFault,
                timeLimited, mresFault, CFG_BLOCK_NAMES[currentBlock]);
  dumpTiming();
  Serial.printf("# SUMMARY,config_sha256=%s,table_sha256=%s,home_neg_MSCNT=%u,"
                "final_mres=%u,final_x_usteps=%ld,final_x_mm=%.6f,"
                "records_dumped=%u\n",
                CFG_HASH, CFG_TABLE_HASH, homeNegMscnt,
                mresOf(currentMresCode), (long)positionU, mmOf(positionU),
                recordCount);
  if (!complete || unsafeFault) {
    digitalWrite(EN_PIN, HIGH);
    positionArmed = false;
  }
  running = false;
}

// A bench check for the one thing that can stop this session: whether the
// RMT path delivers the rate block M needs at the finest resolution.  The
// driver output is disabled for the duration and the same number of pulses
// is sent back, so the stage does not move and MSCNT returns to where it
// started; the test reports whether it did.  Safe to run unhomed.
void rateSelfTest() {
  if (running) {
    Serial.println("# RATE_TEST,result=busy");
    return;
  }
  const uint8_t saved = currentMresCode;
  uint16_t before = 0, after = 0;
  uint32_t chop = 0;
  if (!applyMres(0, &before, &after, &chop)) {
    Serial.println("# RATE_TEST,result=mres_switch_failed");
    return;
  }
  const int32_t savedPosition = positionU;
  const uint16_t mscntBefore = (uint16_t)driver.MSCNT();
  digitalWrite(EN_PIN, HIGH);                  // no motion while we measure
  const double wanted = 2.0 * USTEP_PER_MM;    // 2 mm/s at 1/256 = 51.2 kHz
  const uint32_t pulses = (uint32_t)llround(wanted * 0.25);
  const uint32_t t0 = micros();
  const uint32_t done = emitProfile(pulses, true, wanted);
  const uint32_t t1 = micros();
  emitProfile(done, false, wanted);            // return MSCNT to where it was
  digitalWrite(EN_PIN, LOW);
  positionU = savedPosition;
  const uint16_t mscntAfter = (uint16_t)driver.MSCNT();
  Serial.printf("# RATE_TEST,target_hz=%.0f,pulses=%lu,done=%lu,"
                "elapsed_us=%lu,achieved_hz=%.0f,mscnt_before=%u,"
                "mscnt_after=%u,mscnt_restored=%u,moved=0\n",
                wanted, (unsigned long)pulses, (unsigned long)done,
                (unsigned long)(t1 - t0),
                done ? done * 1e6 / (double)(t1 - t0) : 0.0,
                mscntBefore, mscntAfter, mscntBefore == mscntAfter);
  applyMres(saved, &before, &after, &chop);
  unsafeFault = false;                         // the test is not a run fault
}

void printMenu() {
  Serial.println();
  Serial.println("--- F0b_pt2 reversal vs MRES and speed ---");
  Serial.printf("driver=%u rmt=%u armed=%u position_mm=%.6f mres=1/%u\n",
                driverReady, rmtReady, positionArmed, mmOf(positionU),
                mresOf(currentMresCode));
  Serial.println("g  home, centre, run        h  home and centre only");
  Serial.println("p  arm manually at x0        r  run generated table");
  Serial.println("t  table/hash summary        s  status");
  Serial.println("c  51.2 kHz rate self-test   0  release motor");
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
  Serial.begin(921600);
  delay(300);
  Serial.printf("# BOOT,fw=f0b_pt2,config_sha256=%s\n", CFG_HASH);
  rmtReady = beginRmt();
  driverReady = configureDriver();
  Serial.printf("# READY,driver=%u,rmt=%u,precondition_armed=0\n",
                driverReady, rmtReady);
  Serial.println("# PRECONDITION_REQUIRED,E0_warmup_complete=1,"
                 "IDS_recording=1,scope_record_min=60,rate_bench_check=1");
  printMenu();
}

void loop() {
  if (!Serial.available()) return;
  const int key = Serial.read();
  switch (key) {
    case 'g':
      if (homeAndCentre()) {
        positionArmed = true;
        runSession();
      }
      break;
    case 'h':
      if (homeAndCentre()) positionArmed = true;
      break;
    case 'p':
      if (!driverReady || !rmtReady || switchPressed(SW_NEG_PIN) ||
          switchPressed(SW_POS_PIN)) {
        Serial.println("# ERR,reason=cannot_arm");
      } else {
        digitalWrite(EN_PIN, LOW);
        positionU = CFG_PARK_U;
        positionArmed = true;
        abortRequested = unsafeFault = false;
        Serial.printf("# ARMED,x_usteps=%ld\n", (long)positionU);
      }
      break;
    case 'r':
      runSession();
      break;
    case 'c':
      rateSelfTest();
      break;
    case 't':
      Serial.printf("# PLAN,segments=%u,marks=%u,planned_ms=%lu,"
                    "watchdog_ms=%lu,x0_usteps=%ld,config_sha256=%s,"
                    "table_sha256=%s\n",
                    CFG_SEG_COUNT, CFG_MARK_COUNT,
                    (unsigned long)CFG_PLAN_DURATION_MS,
                    (unsigned long)CFG_WATCHDOG_MS, (long)CFG_X0_U,
                    CFG_HASH, CFG_TABLE_HASH);
      break;
    case 's':
      Serial.printf("# STATUS,x_usteps=%ld,x_mm=%.6f,mres=1/%u,"
                    "mscnt=%u,drv=0x%08lX\n",
                    (long)positionU, mmOf(positionU),
                    mresOf(currentMresCode), (uint16_t)driver.MSCNT(),
                    (unsigned long)driver.DRV_STATUS());
      break;
    case '0':
      digitalWrite(EN_PIN, HIGH);
      positionArmed = false;
      Serial.println("# MOTOR_RELEASED,coordinate_invalid=1");
      break;
    default:
      return;
  }
  printMenu();
}

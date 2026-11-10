/*
 * F0b friction scoping - ESP32-S3 + TMC2209.
 *
 * The host-generated f0b_config.h contains every integer command in execution
 * order.  Nothing is computed here but the step timing.  Five segment types
 * cover the whole session:
 *
 *   MOVE   trapezoidal transit at CFG_V_MM_S (1 mm/s)
 *   RAMP   trapezoidal profile peaking at the segment's own rate.  Below
 *          ~1 kHz the acceleration phase is shorter than one microstep, so a
 *          measurement ramp is constant-rate in practice; at 19.2 and
 *          25.6 kHz it is a real trapezoid, which is what stops the stage
 *          being asked to start at 0.75 or 1 mm/s from standstill.
 *   STEP   n single microsteps with a dwell after each (the L ladder's
 *          1, 2 and 5 ustep strokes)
 *   DWELL  hold in place, motor energised at hold current
 *   MARK   the recorded window of a reference stop
 *
 * Two step paths.  Rates at or above ~33 Hz go through the buffered RMT path
 * proven in E1-E5.  Below that the period exceeds what one RMT symbol can
 * express (CFG_RMT_MAX_PERIOD_US), so pulses are emitted one at a time
 * against a micros() deadline.  Block S runs at 6.4 Hz and block L's 13-ustep
 * stroke at 3.25 Hz, so this path carries real measurement segments, not just
 * corner cases.
 *
 * MRES is 1/256 for the whole session and is never changed: it fixes the
 * step-to-position mapping and the electrical phase reference (procedure
 * section 2).
 *
 * During the run all timing is retained in RAM; UART output resumes only
 * after SEQ_DONE.  Before flashing, perform the E0 warm-up.  Then 'g' to
 * home, centre and run, or 'p'/'r' when already homed and parked at 60.000 mm.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/rmt_tx.h"
#include "f0b_config.h"

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

constexpr double USTEP_PER_MM =
    (double)CFG_STEPS_PER_REV * (double)CFG_MRES / CFG_LEAD_MM;
constexpr double ACCEL_U_S2 = CFG_A_MM_S2 * USTEP_PER_MM;
constexpr double TRANSIT_RATE_HZ = CFG_V_MM_S * USTEP_PER_MM;
// Below this rate one period will not fit in an RMT symbol.
constexpr double RMT_MIN_RATE_HZ = 1000000.0 / (double)CFG_RMT_MAX_PERIOD_US;

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
  uint16_t mscntStart;
  uint16_t mscntEnd;
  uint32_t drvStatus;
};

SegRecord records[CFG_SEG_COUNT];
uint16_t recordCount = 0;
int32_t positionU = CFG_PARK_U;
uint16_t currentBlock = 0;
bool driverReady = false;
bool rmtReady = false;
bool positionArmed = false;
bool abortRequested = false;
bool unsafeFault = false;
bool running = false;
bool timeLimited = false;
bool lastDirPositive = true;
bool directionKnown = false;
uint32_t sequenceStartedMs = 0;
uint16_t homeNegMscnt = 0;
uint16_t homePosMscnt = 0;
uint16_t centreMscnt = 0;
bool homePhaseValid = false;

inline double mmOf(int32_t u) { return (double)u / USTEP_PER_MM; }
inline int32_t ustepsOf(double mm) { return (int32_t)llround(mm * USTEP_PER_MM); }

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
  while (Serial.available()) {
    if ((char)Serial.read() == 3) abortRequested = true;
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

// Wait for an absolute micros() deadline, yielding to the RTOS while there is
// time to spare so the task watchdog is fed during the long slow-ramp gaps.
bool waitUntilUs(uint32_t deadlineUs) {
  for (;;) {
    const int32_t remaining = (int32_t)(deadlineUs - micros());
    if (remaining <= 0) return true;
    if (pollAbort()) return false;
    if (remaining > 3000) delay(1);
    else delayMicroseconds((uint32_t)remaining);
  }
}

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

// One safety gate for every pulse-emitting path.
bool stepAllowed(bool positive, uint32_t count, bool seeking) {
  if (switchAhead(positive)) {
    if (!seeking) unsafeFault = true;
    return false;
  }
  const int32_t next = positionU + (positive ? (int32_t)count : -(int32_t)count);
  if (!seeking && (next > CFG_SOFT_MAX_U || next < CFG_SOFT_MIN_U)) {
    unsafeFault = true;
    return false;
  }
  return true;
}

// Emit one pulse and return without waiting out any gap.
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
  positionU += positive ? 1 : -1;
  return true;
}

// Slow path: one pulse per micros() deadline.  Used below RMT_MIN_RATE_HZ,
// where the acceleration phase is shorter than a single microstep and the
// profile is therefore constant-rate by construction.
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

// Buffered RMT path: trapezoidal profile peaking at rateHz.
uint32_t emitProfile(uint32_t pulses, bool positive, double rateHz,
                     uint32_t chunk = CFG_RMT_CHUNK, bool seeking = false) {
  if (!pulses) return 0;
  setDirection(positive);
  uint32_t nAccel = (uint32_t)(rateHz * rateHz / (2.0 * ACCEL_U_S2));
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
      if (i < nAccel) rate = sqrt(2.0 * ACCEL_U_S2 * (double)(i + 1));
      else if (i >= pulses - nAccel)
        rate = sqrt(2.0 * ACCEL_U_S2 * (double)(pulses - i));
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
    positionU += positive ? (int32_t)count : -(int32_t)count;
  }
  if (rmt_tx_wait_all_done(rmtChannel, -1) != ESP_OK) unsafeFault = true;
  return queued;
}

// Dispatch on rate: the only thing that decides the path.
uint32_t emitAtRate(uint32_t pulses, bool positive, double rateHz) {
  if (!pulses) return 0;
  if (rateHz >= RMT_MIN_RATE_HZ) return emitProfile(pulses, positive, rateHz);
  return emitSlowRate(pulses, positive, rateHz);
}

bool travelTo(int32_t targetU, double rateHz) {
  if (!positionArmed || targetU < CFG_SOFT_MIN_U || targetU > CFG_SOFT_MAX_U) {
    unsafeFault = true;
    return false;
  }
  const int32_t delta = targetU - positionU;
  if (!delta) return true;
  const bool positive = delta > 0;
  const uint32_t wanted = (uint32_t)(positive ? delta : -delta);
  const uint32_t emitted = emitAtRate(wanted, positive, rateHz);
  return emitted == wanted && !unsafeFault && !abortRequested;
}

// n single microsteps towards the target, dwelling after each.
bool stepTo(int32_t targetU, uint32_t dwellMs) {
  const int32_t delta = targetU - positionU;
  if (!delta) return true;
  const bool positive = delta > 0;
  uint32_t remaining = (uint32_t)(positive ? delta : -delta);
  setDirection(positive);
  while (remaining--) {
    if (pollAbort()) return false;
    if (!emitSinglePulse(positive, false)) return false;
    if (!waitMs(dwellMs)) return false;
  }
  return !unsafeFault;
}

// Homing is identical to the proven E0/E1-E4 implementation.
constexpr double HOME_SEEK_MM_S = 2.0;
constexpr double HOME_SLOW_MM_S = 0.15;
constexpr double HOME_BACKOFF_MM = 2.0;
constexpr uint8_t HOME_TOUCHES = 3;
constexpr double HOME_TOL_MM = 0.02;
constexpr double HOME_MAX_TRAVEL_MM = 200.0;
constexpr double END_MARGIN_MM = 2.0;
constexpr uint32_t HOME_CHUNK = 16;

void jog(double mm, double vMmS) {
  const int32_t delta = ustepsOf(mm);
  emitProfile((uint32_t)(delta >= 0 ? delta : -delta), delta >= 0,
              vMmS * USTEP_PER_MM, HOME_CHUNK, true);
}

bool homeEnd(bool positive, int32_t *trigger) {
  const double sign = positive ? 1.0 : -1.0;
  if (switchAhead(positive)) jog(-sign * HOME_BACKOFF_MM, HOME_SLOW_MM_S);
  jog(sign * HOME_MAX_TRAVEL_MM, HOME_SEEK_MM_S);
  if (abortRequested || !switchAhead(positive)) return false;
  int32_t lo = 0, hi = 0;
  for (uint8_t touch = 0; touch < HOME_TOUCHES; ++touch) {
    jog(-sign * HOME_BACKOFF_MM, HOME_SLOW_MM_S);
    jog(sign * 3.0 * HOME_BACKOFF_MM, HOME_SLOW_MM_S);
    if (abortRequested || !switchAhead(positive)) return false;
    if (touch == 0 || positionU < lo) lo = positionU;
    if (touch == 0 || positionU > hi) hi = positionU;
  }
  Serial.printf("# HOME_END,%s,spread_mm=%.4f\n",
                positive ? "pos" : "neg", mmOf(hi - lo));
  if (mmOf(hi - lo) > HOME_TOL_MM) return false;
  *trigger = positive ? lo : hi;
  return true;
}

bool homeAndCentre() {
  positionArmed = false;
  abortRequested = unsafeFault = false;
  digitalWrite(EN_PIN, LOW);
  int32_t negTrigger = 0, posTrigger = 0;
  if (!homeEnd(false, &negTrigger)) {
    Serial.println("# ERR,reason=home_neg");
    return false;
  }
  homeNegMscnt = (uint16_t)driver.MSCNT();
  jog(HOME_BACKOFF_MM, HOME_SLOW_MM_S);
  if (!homeEnd(true, &posTrigger)) {
    Serial.println("# ERR,reason=home_pos");
    return false;
  }
  homePosMscnt = (uint16_t)driver.MSCNT();
  positionU -= negTrigger;
  const int32_t travel = posTrigger - negTrigger;
  Serial.printf("# HOMED,travel_mm=%.4f,soft_mm=%.3f..%.3f\n",
                mmOf(travel), mmOf(CFG_SOFT_MIN_U), mmOf(CFG_SOFT_MAX_U));
  if (CFG_SOFT_MIN_U < ustepsOf(END_MARGIN_MM) ||
      CFG_SOFT_MAX_U > travel - ustepsOf(END_MARGIN_MM)) {
    Serial.println("# ERR,reason=travel_short_for_soft_limits");
    return false;
  }
  positionArmed = true;
  // Approach the park point from below, so the session starts with a wipe in
  // the same direction every block assumes.
  const int32_t pre = CFG_PARK_U - ustepsOf(3.0);
  if (!travelTo(pre, TRANSIT_RATE_HZ) || !travelTo(CFG_PARK_U, TRANSIT_RATE_HZ))
    return false;
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
  // MRES field: 0 = 1/256, 1 = 1/128, ... 8 = full step.
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
                        !(chop & (1UL << 28)) &&
                        (((chop & mask) >> 24) == code) &&
                        driver.ihold() == driver.irun();
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,status=0x%08lX,"
                "mres=%u,irun=%u,ihold=%u,tpowerdown=%u,iholddelay=%u\n",
                (unsigned long)gconf, (unsigned long)chop,
                (unsigned long)status, CFG_MRES, driver.irun(), driver.ihold(),
                CFG_TPOWERDOWN, CFG_IHOLDDELAY);
  if (!readback) return false;
  digitalWrite(EN_PIN, LOW);
  return true;
}

bool driverStatusFault(uint32_t status) {
  return status & ((1UL << 1) | (1UL << 2) | (1UL << 3) |
                   (1UL << 4) | (1UL << 5));
}

bool executeSegment(const F0bSeg &seg) {
  if (recordCount >= CFG_SEG_COUNT) {
    unsafeFault = true;
    return false;
  }
  SegRecord &record = records[recordCount];
  record.row = seg.row;
  record.mscntStart = (uint16_t)driver.MSCNT();
  record.startMs = millis();
  bool ok = true;
  switch (seg.type) {
    case F0B_MOVE:
      ok = travelTo(seg.u_target, TRANSIT_RATE_HZ);
      break;
    case F0B_RAMP:
      ok = travelTo(seg.u_target, (double)seg.rate_mhz / 1000.0);
      break;
    case F0B_STEP:
      ok = stepTo(seg.u_target, seg.dwell_ms);
      break;
    case F0B_DWELL:
      ok = waitMs(seg.dwell_ms);
      break;
    case F0B_MARK:
      syncWrite(true);
      ok = waitMs(seg.dwell_ms);
      syncWrite(false);
      break;
    default:
      unsafeFault = true;
      ok = false;
      break;
  }
  record.endMs = millis();
  record.mscntEnd = (uint16_t)driver.MSCNT();
  record.drvStatus = driver.DRV_STATUS();
  ++recordCount;
  if (driverStatusFault(record.drvStatus)) {
    unsafeFault = true;
    return false;
  }
  return ok;
}

void printHeader() {
  Serial.printf("# SESSION,experiment=%s,procedure=v2,config_sha256=%s,"
                "table_sha256=%s,segments=%u,marks=%u,motion=RMT+timed\n",
                CFG_EXPERIMENT, CFG_HASH, CFG_TABLE_HASH, CFG_SEG_COUNT,
                CFG_MARK_COUNT);
  Serial.printf("# SETTINGS,lead_mm=%.3f,steps_rev=%u,mres=%u,chopper=SpreadCycle,"
                "run_ma=%u,hold_ma=%u,v_mm_s=%.3f,a_mm_s2=%.3f,"
                "settle_ms=%lu,record_ms=%lu,ustep_per_mm=%.1f,"
                "rmt_min_rate_hz=%.2f,aquadb_clock_ns=%u\n",
                CFG_LEAD_MM, CFG_STEPS_PER_REV, CFG_MRES, CFG_RUN_CURRENT_MA,
                CFG_HOLD_CURRENT_MA, CFG_V_MM_S, CFG_A_MM_S2,
                (unsigned long)CFG_SETTLE_MS, (unsigned long)CFG_RECORD_MS,
                USTEP_PER_MM, RMT_MIN_RATE_HZ, CFG_AQUADB_CLOCK_NS);
  Serial.printf("# HOME_PHASE,valid=%u,home_neg_MSCNT=%u,home_pos_MSCNT=%u,"
                "centre_MSCNT=%u,x_usteps=%ld,x_mm=%.6f\n",
                homePhaseValid, homeNegMscnt, homePosMscnt, centreMscnt,
                (long)positionU, mmOf(positionU));
}

const char *typeName(uint8_t type) {
  switch (type) {
    case F0B_MOVE: return "move";
    case F0B_RAMP: return "ramp";
    case F0B_STEP: return "step";
    case F0B_DWELL: return "dwell";
    case F0B_MARK: return "mark";
    default: return "?";
  }
}

void dumpTiming() {
  Serial.println("# SEG_COLUMNS,row,block,seg,u_target,rate_mhz,dwell_ms,"
                 "el_off,MSCNT_start,MSCNT_end,start_ms,end_ms,drv_status");
  for (uint16_t i = 0; i < recordCount; ++i) {
    const SegRecord &record = records[i];
    const F0bSeg &seg = CFG_SEGS[record.row];
    const int16_t elOff =
        (int16_t)((seg.u_target % CFG_MOTOR_CYCLE_U + CFG_MOTOR_CYCLE_U +
                   CFG_MOTOR_CYCLE_U / 2) % CFG_MOTOR_CYCLE_U) -
        (int16_t)(CFG_MOTOR_CYCLE_U / 2);
    Serial.printf("# SEG,%u,%s,%s,%ld,%lu,%lu,%d,%u,%u,%lu,%lu,0x%08lX\n",
                  seg.row, CFG_BLOCK_NAMES[seg.block], typeName(seg.type),
                  (long)seg.u_target, (unsigned long)seg.rate_mhz,
                  (unsigned long)seg.dwell_ms, elOff, record.mscntStart,
                  record.mscntEnd, (unsigned long)record.startMs,
                  (unsigned long)record.endMs,
                  (unsigned long)record.drvStatus);
  }
}

void runF0b() {
  if (!driverReady || !rmtReady || !positionArmed || running) {
    Serial.printf("# ERR,reason=not_ready,driver=%u,rmt=%u,armed=%u,running=%u\n",
                  driverReady, rmtReady, positionArmed, running);
    return;
  }
  running = true;
  abortRequested = unsafeFault = timeLimited = false;
  recordCount = 0;
  sequenceStartedMs = millis();
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
                "complete=%u,aborted=%u,unsafe=%u,time_limited=%u,block=%s\n",
                (unsigned long)elapsed, (unsigned long)CFG_PLAN_DURATION_MS,
                recordCount, complete, abortRequested, unsafeFault, timeLimited,
                CFG_BLOCK_NAMES[currentBlock]);
  dumpTiming();
  Serial.printf("# SUMMARY,config_sha256=%s,table_sha256=%s,home_neg_MSCNT=%u,"
                "final_x_usteps=%ld,final_x_mm=%.6f,records_dumped=%u\n",
                CFG_HASH, CFG_TABLE_HASH, homeNegMscnt, (long)positionU,
                mmOf(positionU), recordCount);
  if (!complete || unsafeFault) {
    digitalWrite(EN_PIN, HIGH);
    positionArmed = false;
  }
  running = false;
}

void printMenu() {
  Serial.println("\n--- F0b friction scoping ---");
  Serial.printf("driver=%u rmt=%u armed=%u position_mm=%.6f\n",
                driverReady, rmtReady, positionArmed, mmOf(positionU));
  Serial.println("g  home, centre, run        h  home and centre only");
  Serial.println("p  arm manually at x0        r  run generated table");
  Serial.println("t  table/hash summary        s  status");
  Serial.println("0  release motor             Ctrl+C aborts");
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
  Serial.printf("# BOOT,fw=f0b_friction_scoping,config_sha256=%s\n", CFG_HASH);
  rmtReady = beginRmt();
  driverReady = configureDriver();
  Serial.printf("# READY,driver=%u,rmt=%u,precondition_armed=0\n",
                driverReady, rmtReady);
  Serial.println("# PRECONDITION_REQUIRED,E0_warmup_complete=1,IDS_recording=1,"
                 "scope_record_min=60");
  printMenu();
}

void loop() {
  if (!Serial.available()) {
    delay(5);
    return;
  }
  const char key = (char)Serial.read();
  if (key == '\r' || key == '\n') return;
  if (key == 3) {
    abortRequested = true;
    Serial.println("# ABORT_REQUEST");
    return;
  }
  switch (key) {
    case 'p':
      if (!driverReady || !rmtReady || switchPressed(SW_NEG_PIN) ||
          switchPressed(SW_POS_PIN)) {
        Serial.println("# ERR,reason=cannot_arm");
      } else {
        digitalWrite(EN_PIN, LOW);
        positionU = CFG_PARK_U;
        positionArmed = true;
        abortRequested = unsafeFault = false;
        centreMscnt = (uint16_t)driver.MSCNT();
        homePhaseValid = false;
        Serial.printf("# PRECONDITION_ARMED,x_usteps=%ld,x_mm=%.6f,MSCNT=%u,"
                      "home_phase_valid=0\n",
                      (long)positionU, mmOf(positionU), centreMscnt);
      }
      break;
    case 'r': runF0b(); break;
    case 'h': homeAndCentre(); break;
    case 'g': if (homeAndCentre()) runF0b(); break;
    case 's':
      Serial.printf("# STATUS,block=%s,armed=%u,x_usteps=%ld,x_mm=%.6f,"
                    "neg=%u,pos=%u,MSCNT=%u,status=0x%08lX\n",
                    CFG_BLOCK_NAMES[currentBlock], positionArmed,
                    (long)positionU, mmOf(positionU), switchPressed(SW_NEG_PIN),
                    switchPressed(SW_POS_PIN), (uint16_t)driver.MSCNT(),
                    (unsigned long)driver.DRV_STATUS());
      break;
    case 't':
      Serial.printf("# PLAN,segments=%u,marks=%u,planned_ms=%lu,watchdog_ms=%lu,"
                    "x0_usteps=%ld,config_sha256=%s,table_sha256=%s\n",
                    CFG_SEG_COUNT, CFG_MARK_COUNT,
                    (unsigned long)CFG_PLAN_DURATION_MS,
                    (unsigned long)CFG_WATCHDOG_MS, (long)CFG_X0_U,
                    CFG_HASH, CFG_TABLE_HASH);
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

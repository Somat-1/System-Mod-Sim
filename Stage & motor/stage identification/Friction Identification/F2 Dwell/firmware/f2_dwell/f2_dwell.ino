/*
 * F2 dwell (frictional ageing and creep) - ESP32-S3 + TMC2209.
 *
 * Same engine as F0b and F1, same five segment types, two test positions:
 *
 *   MOVE   trapezoidal transit at CFG_V_MM_S (1 mm/s)
 *   RAMP   trapezoidal profile peaking at the segment's own rate
 *   STEP   n single microsteps with a dwell after each
 *   DWELL  hold in place, motor energised at hold current
 *   MARK   the recorded window of a reference stop
 *
 * F2 emits no STEP segments; the path is retained so that the motion engine
 * stays identical to the one F0b and F1 use.
 *
 * WHAT IS DIFFERENT IN F2.  Here the DWELL *is* the measurement, not an
 * interval between measurements, so:
 *
 *   - every DWELL and MARK is timed against a micros() deadline rather than
 *     a 2 ms millis() poll, and the achieved duration is recorded;
 *   - each segment carries its 'kind', which separates the trial's timed
 *     dwell from the 0.5 s settle before it, the 2 s hold after the ramp and
 *     the reference stops - all four are the same segment type;
 *   - the run ends with a DWELL_CHECK line: the largest timing error over
 *     the trial dwells, against CFG_DWELL_TOL_MS (procedure section 5);
 *   - 'd' runs that timer check on its own, without moving the stage, so the
 *     pre-run checklist can be completed before the session starts.
 *
 * IRUN = IHOLD, TPOWERDOWN = 0 and IHOLDDELAY = 0, so the current during a
 * 300 s dwell is the current during motion; configureDriver() reads this
 * back and refuses to run if the driver disagrees.  The motor is never
 * released between trials: releasing it would change the hold state that
 * the dwell is supposed to age.
 *
 * Rates at or above ~33 Hz use the buffered RMT path; below that one period
 * will not fit in an RMT symbol, so pulses go out one at a time against a
 * micros() deadline.  F2's post-dwell ramp is 12.8 Hz, so every measurement
 * ramp in this session runs on the timed path.
 *
 * PROVISIONAL TABLE.  The post-dwell ramp is 2*s90 long and s90 comes from
 * F0b.  Until F0b is processed, f2_config.h is generated from the
 * procedure's placeholder defaults and CFG_PROVISIONAL is 1.  F2 is also
 * conditional: the full grid is justified only if F0b's dwell probe (1 s vs
 * 60 s) exceeded 3*d_min, otherwise the reduced variant is the one to run.
 * Re-run f2_make_table.py after F0b and reflash.
 *
 * MRES is 1/256 for the whole session and is never changed.  During the run
 * all timing is retained in RAM; UART output resumes only after SEQ_DONE.
 * 'g' homes, centres and runs; 'p'/'r' when already homed and parked at x_B.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/rmt_tx.h"
#include "f2_config.h"

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

/* kind indices, in the order f2_make_table.py writes CFG_KIND_NAMES */
constexpr uint8_t KIND_TRANSIT = 0;
constexpr uint8_t KIND_SETUP = 1;
constexpr uint8_t KIND_SETTLE = 2;
constexpr uint8_t KIND_DWELL = 3;
constexpr uint8_t KIND_MEASURE = 4;
constexpr uint8_t KIND_HOLD = 5;
constexpr uint8_t KIND_REF = 6;

constexpr double USTEP_PER_MM =
    (double)CFG_STEPS_PER_REV * (double)CFG_MRES / CFG_LEAD_MM;
constexpr double ACCEL_U_S2 = CFG_A_MM_S2 * USTEP_PER_MM;
constexpr double TRANSIT_RATE_HZ = CFG_V_MM_S * USTEP_PER_MM;
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
  uint32_t heldUs;        /* achieved duration of a DWELL or MARK */
  uint16_t mscntStart;
  uint16_t mscntEnd;
  uint32_t drvStatus;
};

SegRecord records[CFG_SEG_COUNT];
uint16_t recordCount = 0;
int32_t positionU = CFG_PARK_U;
uint16_t currentBlock = 0;
uint16_t currentSub = 0;
uint16_t currentTrial = 0;
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
uint16_t runIrun = 0;
uint16_t runIhold = 0;

/* dwell timing, accumulated over the trial dwells only */
int32_t worstDwellErrorUs = 0;
uint16_t worstDwellRow = 0;
uint16_t dwellsTimed = 0;
uint16_t dwellsOutOfTolerance = 0;

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

bool waitUntilUs(uint32_t deadlineUs) {
  for (;;) {
    const int32_t remaining = (int32_t)(deadlineUs - micros());
    if (remaining <= 0) return true;
    if (pollAbort()) return false;
    if (remaining > 3000) delay(1);
    else delayMicroseconds((uint32_t)remaining);
  }
}

/*
 * Hold position for exactly durationMs and report what was achieved.
 * The deadline is taken once, in microseconds, so the error does not
 * accumulate with the length of the dwell; the longest dwell in the table is
 * CFG_LONGEST_DWELL_MS, well inside the 71 min micros() period and inside
 * the int32 range used for the remaining time.  No step pulse is emitted,
 * and the driver is left energised at hold current.
 */
bool holdExactly(uint32_t durationMs, uint32_t *achievedUs) {
  const uint32_t startUs = micros();
  const uint32_t deadline = startUs + durationMs * 1000UL;
  bool ok = true;
  for (;;) {
    const int32_t remaining = (int32_t)(deadline - micros());
    if (remaining <= 0) break;
    if (pollAbort()) {
      ok = false;
      break;
    }
    if (remaining > 4000) delay(1);
    else delayMicroseconds((uint32_t)remaining);
  }
  *achievedUs = micros() - startUs;
  return ok;
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
  driver.ihold(driver.irun());          /* IRUN = IHOLD: section 2 */
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
  runIrun = driver.irun();
  runIhold = driver.ihold();
  const bool readback = (gconf & (1UL << 2)) && !(status & (1UL << 30)) &&
                        !(chop & (1UL << 28)) &&
                        (((chop & mask) >> 24) == code) &&
                        runIhold == runIrun;
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,status=0x%08lX,"
                "mres=%u,irun=%u,ihold=%u,irun_eq_ihold=%u,tpowerdown=%u,"
                "iholddelay=%u\n",
                (unsigned long)gconf, (unsigned long)chop,
                (unsigned long)status, CFG_MRES, runIrun, runIhold,
                runIrun == runIhold, CFG_TPOWERDOWN, CFG_IHOLDDELAY);
  if (!readback) {
    Serial.println("# ERR,reason=driver_readback,hold_current_must_equal_run");
    return false;
  }
  digitalWrite(EN_PIN, LOW);
  return true;
}

bool driverStatusFault(uint32_t status) {
  return status & ((1UL << 1) | (1UL << 2) | (1UL << 3) |
                   (1UL << 4) | (1UL << 5));
}

bool executeSegment(const F2Seg &seg) {
  if (recordCount >= CFG_SEG_COUNT) {
    unsafeFault = true;
    return false;
  }
  SegRecord &record = records[recordCount];
  record.row = seg.row;
  record.heldUs = 0;
  record.mscntStart = (uint16_t)driver.MSCNT();
  record.startMs = millis();
  bool ok = true;
  switch (seg.type) {
    case F2_MOVE:
      ok = travelTo(seg.u_target, TRANSIT_RATE_HZ);
      break;
    case F2_RAMP:
      ok = travelTo(seg.u_target, (double)seg.rate_mhz / 1000.0);
      break;
    case F2_STEP:                       /* unused in F2; kept for engine parity */
      ok = stepTo(seg.u_target, seg.dwell_ms);
      break;
    case F2_DWELL:
      ok = holdExactly(seg.dwell_ms, &record.heldUs);
      break;
    case F2_MARK:
      syncWrite(true);
      ok = holdExactly(seg.dwell_ms, &record.heldUs);
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
  /* The trial dwell is the measurement: hold its timing to account. */
  if (seg.kind == KIND_DWELL && ok) {
    const int32_t errorUs =
        (int32_t)(record.heldUs - seg.dwell_ms * 1000UL);
    ++dwellsTimed;
    if ((uint32_t)abs(errorUs) > CFG_DWELL_TOL_MS * 1000UL)
      ++dwellsOutOfTolerance;
    if (abs(errorUs) > abs(worstDwellErrorUs)) {
      worstDwellErrorUs = errorUs;
      worstDwellRow = seg.row;
    }
  }
  /* The carriage must not creep away under a loss of hold current. */
  if (driverStatusFault(record.drvStatus)) {
    unsafeFault = true;
    return false;
  }
  return ok;
}

void printHeader() {
  Serial.printf("# SESSION,experiment=%s,variant=%s,config_sha256=%s,"
                "table_sha256=%s,segments=%u,trials=%u,dwell_windows=%u,"
                "marks=%u,provisional=%u,s90_um=%.4f,d_min_um=%.4f,"
                "motion=RMT+timed\n",
                CFG_EXPERIMENT, CFG_VARIANT, CFG_HASH, CFG_TABLE_HASH,
                CFG_SEG_COUNT, CFG_TRIAL_COUNT, CFG_DWELL_WINDOWS,
                CFG_MARK_COUNT, CFG_PROVISIONAL, CFG_S90_UM, CFG_D_MIN_UM);
  if (CFG_PROVISIONAL)
    Serial.println("# WARNING,provisional_post_dwell_ramp=1,"
                   "regenerate_after_F0b=1,f2_conditional_on_dwell_probe=1");
  Serial.printf("# SETTINGS,lead_mm=%.3f,steps_rev=%u,mres=%u,chopper=SpreadCycle,"
                "intpol=0,run_ma=%u,hold_ma=%u,irun=%u,ihold=%u,v_mm_s=%.3f,"
                "a_mm_s2=%.3f,ustep_per_mm=%.1f,rmt_min_rate_hz=%.2f,"
                "aquadb_clock_ns=%u\n",
                CFG_LEAD_MM, CFG_STEPS_PER_REV, CFG_MRES, CFG_RUN_CURRENT_MA,
                CFG_HOLD_CURRENT_MA, runIrun, runIhold, CFG_V_MM_S,
                CFG_A_MM_S2, USTEP_PER_MM, RMT_MIN_RATE_HZ,
                CFG_AQUADB_CLOCK_NS);
  Serial.printf("# DWELL_PLAN,seed=%llu,settle_ms=%lu,hold_ms=%lu,"
                "longest_dwell_ms=%lu,tolerance_ms=%lu,ramp_usteps=%lu,"
                "ramp_rate_hz=%.3f\n",
                (unsigned long long)CFG_SHUFFLE_SEED,
                (unsigned long)CFG_PRE_DWELL_SETTLE_MS,
                (unsigned long)CFG_POST_RAMP_HOLD_MS,
                (unsigned long)CFG_LONGEST_DWELL_MS,
                (unsigned long)CFG_DWELL_TOL_MS,
                (unsigned long)CFG_RAMP_USTEPS, CFG_RAMP_RATE_HZ);
  Serial.printf("# POSITIONS,x_A_usteps=%ld,x_B_usteps=%ld\n",
                (long)CFG_X_A_U, (long)CFG_X_B_U);
  Serial.printf("# HOME_PHASE,valid=%u,home_neg_MSCNT=%u,home_pos_MSCNT=%u,"
                "centre_MSCNT=%u,x_usteps=%ld,x_mm=%.6f\n",
                homePhaseValid, homeNegMscnt, homePosMscnt, centreMscnt,
                (long)positionU, mmOf(positionU));
}

const char *typeName(uint8_t type) {
  switch (type) {
    case F2_MOVE: return "move";
    case F2_RAMP: return "ramp";
    case F2_STEP: return "step";
    case F2_DWELL: return "dwell";
    case F2_MARK: return "mark";
    default: return "?";
  }
}

void dumpTiming() {
  Serial.println("# SEG_COLUMNS,row,block,sub,kind,seg,trial,u_target,rate_mhz,"
                 "dwell_ms,held_us,dwell_err_us,el_off,MSCNT_start,MSCNT_end,"
                 "start_ms,end_ms,drv_status");
  for (uint16_t i = 0; i < recordCount; ++i) {
    const SegRecord &record = records[i];
    const F2Seg &seg = CFG_SEGS[record.row];
    const int16_t elOff =
        (int16_t)((seg.u_target % CFG_MOTOR_CYCLE_U + CFG_MOTOR_CYCLE_U +
                   CFG_MOTOR_CYCLE_U / 2) % CFG_MOTOR_CYCLE_U) -
        (int16_t)(CFG_MOTOR_CYCLE_U / 2);
    const int32_t errorUs =
        (seg.type == F2_DWELL || seg.type == F2_MARK)
            ? (int32_t)(record.heldUs - seg.dwell_ms * 1000UL)
            : 0;
    Serial.printf("# SEG,%u,%s,%s,%s,%s,%u,%ld,%lu,%lu,%lu,%ld,%d,%u,%u,%lu,"
                  "%lu,0x%08lX\n",
                  seg.row, CFG_BLOCK_NAMES[seg.block], CFG_SUB_NAMES[seg.sub],
                  CFG_KIND_NAMES[seg.kind], typeName(seg.type), seg.trial,
                  (long)seg.u_target, (unsigned long)seg.rate_mhz,
                  (unsigned long)seg.dwell_ms, (unsigned long)record.heldUs,
                  (long)errorUs, elOff, record.mscntStart, record.mscntEnd,
                  (unsigned long)record.startMs, (unsigned long)record.endMs,
                  (unsigned long)record.drvStatus);
  }
}

void printDwellCheck() {
  Serial.printf("# DWELL_CHECK,dwells_timed=%u,planned=%u,tolerance_ms=%lu,"
                "worst_error_us=%ld,worst_row=%u,out_of_tolerance=%u,pass=%u\n",
                dwellsTimed, CFG_DWELL_WINDOWS,
                (unsigned long)CFG_DWELL_TOL_MS, (long)worstDwellErrorUs,
                worstDwellRow, dwellsOutOfTolerance,
                (dwellsOutOfTolerance == 0 &&
                 dwellsTimed == CFG_DWELL_WINDOWS) ? 1 : 0);
}

/*
 * Pre-run check, section 5: verify the dwell timer without moving the stage.
 * The stage stays where it is, the driver stays energised, and the achieved
 * durations of a short, a medium and the session's longest dwell are
 * reported against the tolerance.
 */
void dwellTimerSelfTest() {
  const uint32_t trials[] = {100UL, 2000UL, 30000UL, CFG_LONGEST_DWELL_MS};
  Serial.printf("# DWELL_SELFTEST_START,tolerance_ms=%lu,no_motion=1\n",
                (unsigned long)CFG_DWELL_TOL_MS);
  int32_t worst = 0;
  bool ok = true;
  for (uint8_t i = 0; i < sizeof(trials) / sizeof(trials[0]); ++i) {
    uint32_t heldUs = 0;
    const int32_t before = positionU;
    const bool done = holdExactly(trials[i], &heldUs);
    const int32_t errorUs = (int32_t)(heldUs - trials[i] * 1000UL);
    if (abs(errorUs) > abs(worst)) worst = errorUs;
    if (!done || positionU != before ||
        (uint32_t)abs(errorUs) > CFG_DWELL_TOL_MS * 1000UL)
      ok = false;
    Serial.printf("# DWELL_SELFTEST,planned_ms=%lu,held_us=%lu,err_us=%ld,"
                  "moved=%u,aborted=%u\n",
                  (unsigned long)trials[i], (unsigned long)heldUs,
                  (long)errorUs, positionU != before, !done);
    if (!done) break;
  }
  Serial.printf("# DWELL_SELFTEST_DONE,worst_error_us=%ld,pass=%u\n",
                (long)worst, ok);
  abortRequested = false;
}

void runF2() {
  if (!driverReady || !rmtReady || !positionArmed || running) {
    Serial.printf("# ERR,reason=not_ready,driver=%u,rmt=%u,armed=%u,running=%u\n",
                  driverReady, rmtReady, positionArmed, running);
    return;
  }
  running = true;
  abortRequested = unsafeFault = timeLimited = false;
  recordCount = 0;
  dwellsTimed = dwellsOutOfTolerance = worstDwellRow = 0;
  worstDwellErrorUs = 0;
  sequenceStartedMs = millis();
  centreMscnt = (uint16_t)driver.MSCNT();
  printHeader();
  Serial.println("# SEQ_START,uart_quiet=1,timing_buffer=RAM,motor_energised=1");
  bool ok = true;
  for (uint16_t i = 0; ok && i < CFG_SEG_COUNT; ++i) {
    /* The watchdog is tested between segments, so the longest dwell can
       overrun the gate by at most CFG_LONGEST_DWELL_MS; f2_make_table.py
       checks that gate + longest dwell still fits the scope record. */
    if (millis() - sequenceStartedMs >= CFG_WATCHDOG_MS) {
      timeLimited = true;
      ok = false;
      break;
    }
    currentBlock = CFG_SEGS[i].block;
    currentSub = CFG_SEGS[i].sub;
    currentTrial = CFG_SEGS[i].trial;
    ok = executeSegment(CFG_SEGS[i]);
    if (abortRequested || unsafeFault) ok = false;
  }
  const uint32_t elapsed = millis() - sequenceStartedMs;
  const bool complete = ok && recordCount == CFG_SEG_COUNT;
  Serial.printf("# SEQ_DONE,elapsed_ms=%lu,planned_ms=%lu,segments=%u,"
                "complete=%u,aborted=%u,unsafe=%u,time_limited=%u,"
                "block=%s,sub=%s,trial=%u\n",
                (unsigned long)elapsed, (unsigned long)CFG_PLAN_DURATION_MS,
                recordCount, complete, abortRequested, unsafeFault, timeLimited,
                CFG_BLOCK_NAMES[currentBlock], CFG_SUB_NAMES[currentSub],
                currentTrial);
  dumpTiming();
  printDwellCheck();
  Serial.printf("# SUMMARY,config_sha256=%s,table_sha256=%s,variant=%s,"
                "provisional=%u,seed=%llu,irun=%u,ihold=%u,home_neg_MSCNT=%u,"
                "final_x_usteps=%ld,final_x_mm=%.6f,records_dumped=%u\n",
                CFG_HASH, CFG_TABLE_HASH, CFG_VARIANT, CFG_PROVISIONAL,
                (unsigned long long)CFG_SHUFFLE_SEED, runIrun, runIhold,
                homeNegMscnt, (long)positionU, mmOf(positionU), recordCount);
  /* The motor is left energised on a clean run: the hold state is part of
     what F2 measures, and releasing it would wipe it. */
  if (!complete || unsafeFault) {
    digitalWrite(EN_PIN, HIGH);
    positionArmed = false;
  }
  running = false;
}

void printMenu() {
  Serial.println("\n--- F2 dwell (frictional ageing and creep) ---");
  Serial.printf("variant=%s driver=%u rmt=%u armed=%u position_mm=%.6f "
                "provisional=%u\n",
                CFG_VARIANT, driverReady, rmtReady, positionArmed,
                mmOf(positionU), CFG_PROVISIONAL);
  Serial.println("g  home, centre, run        h  home and centre only");
  Serial.println("p  arm manually at x_B       r  run generated table");
  Serial.println("d  dwell timer check         t  table/hash summary");
  Serial.println("s  status                    0  release motor");
  Serial.println("Ctrl+C aborts");
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
  Serial.printf("# BOOT,fw=f2_dwell,variant=%s,config_sha256=%s,provisional=%u\n",
                CFG_VARIANT, CFG_HASH, CFG_PROVISIONAL);
  if (CFG_PROVISIONAL)
    Serial.printf("# WARNING,post_dwell_ramp_from_placeholder_s90_um=%.4f,"
                  "d_min_um=%.4f,rerun_f2_make_table_after_F0b=1\n",
                  CFG_S90_UM, CFG_D_MIN_UM);
  rmtReady = beginRmt();
  driverReady = configureDriver();
  Serial.printf("# READY,driver=%u,rmt=%u,precondition_armed=0\n",
                driverReady, rmtReady);
  Serial.println("# PRECONDITION_REQUIRED,E0_warmup_complete=1,IDS_recording=1,"
                 "scope_record_min=60,F0b_processed=1,dwell_probe_checked=1,"
                 "temperature_logging=1");
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
    case 'r': runF2(); break;
    case 'h': homeAndCentre(); break;
    case 'g': if (homeAndCentre()) runF2(); break;
    case 'd': dwellTimerSelfTest(); break;
    case 's':
      Serial.printf("# STATUS,block=%s,sub=%s,trial=%u,armed=%u,x_usteps=%ld,"
                    "x_mm=%.6f,neg=%u,pos=%u,MSCNT=%u,irun=%u,ihold=%u,"
                    "status=0x%08lX\n",
                    CFG_BLOCK_NAMES[currentBlock], CFG_SUB_NAMES[currentSub],
                    currentTrial, positionArmed, (long)positionU,
                    mmOf(positionU), switchPressed(SW_NEG_PIN),
                    switchPressed(SW_POS_PIN), (uint16_t)driver.MSCNT(),
                    runIrun, runIhold, (unsigned long)driver.DRV_STATUS());
      break;
    case 't':
      Serial.printf("# PLAN,variant=%s,segments=%u,trials=%u,dwell_windows=%u,"
                    "marks=%u,planned_ms=%lu,watchdog_ms=%lu,longest_dwell_ms=%lu,"
                    "seed=%llu,provisional=%u,s90_um=%.4f,config_sha256=%s,"
                    "table_sha256=%s\n",
                    CFG_VARIANT, CFG_SEG_COUNT, CFG_TRIAL_COUNT,
                    CFG_DWELL_WINDOWS, CFG_MARK_COUNT,
                    (unsigned long)CFG_PLAN_DURATION_MS,
                    (unsigned long)CFG_WATCHDOG_MS,
                    (unsigned long)CFG_LONGEST_DWELL_MS,
                    (unsigned long long)CFG_SHUFFLE_SEED, CFG_PROVISIONAL,
                    CFG_S90_UM, CFG_HASH, CFG_TABLE_HASH);
      break;
    case '0':
      digitalWrite(EN_PIN, HIGH);
      positionArmed = false;
      Serial.println("# MOTOR_RELEASED,coordinate_invalid=1,hold_state_wiped=1");
      break;
    case '?': break;
    default: Serial.println("# ERR,reason=unknown_command"); break;
  }
  printMenu();
}

/*
 * Settling matrix across microstep sizes - ESP32-S3 + TMC2209 + end switches.
 *
 * NOT YET COMPILED OR RUN ON HARDWARE. No toolchain was available where this
 * was written. Compile it, set the switch pins, and bring it up with the
 * motor uncoupled.
 *
 * The question this rig answers: when the stage is told to move a given
 * distance, how far does it actually end up moving, and does the microstep
 * size change that? So the measurement is a DISPLACEMENT between two settled
 * stops, not a position against a datum. Nothing here needs an absolute
 * reference or a drift correction.
 *
 * The board cannot see the stage - the interferometer and the EL5101 counter
 * record separately - so its job is to command the moves and to mark the
 * recording precisely enough that each settled stop can be found afterwards.
 * Every stop emits
 *
 *     # STOP,<row>,<tag>,d=<u256>,dir=<+-1>,measured=<0|1>,pos=<u256>,mscnt=<n>
 *
 * at the moment the dwell ends, and the block markers put a distinctive
 * excursion in the recording so the two clocks can be aligned.
 *
 * POSITIONS are int64 u256 = 1/256 full step. 100 full steps = 1 mm, so
 * 1 u256 = 39.06 nm. After homing, 0 is the centre of the measured travel,
 * snapped to a full step.
 *
 * WIRING - the Stepper Isolated bench (see ../../../../CONNECTIONS.md):
 *
 *   STEP GPIO7   DIR GPIO8   EN GPIO9 (ENN, active low)
 *   TMC UART TX GPIO17 (via 1k) / RX GPIO18
 *   END SWITCH neg GPIO10   pos GPIO11
 *
 * Switches are 10T85 microswitches, COM to common ground and NO to the GPIO,
 * so the pins use the internal pull-up and read LOW when pressed. Not
 * fail-safe: a cut wire reads as clear. Moving each to the NC contact and
 * setting SWITCH_ACTIVE_LEVEL HIGH would stop the stage on a broken cable.
 *
 * Nothing moves until a mode is chosen from the menu.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "settling_plan.h"

// ---- pins -----------------------------------------------------------------
constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr uint8_t SW_NEG_PIN = 10;
constexpr uint8_t SW_POS_PIN = 11;
constexpr int SWITCH_ACTIVE_LEVEL = LOW;

// ---- driver ---------------------------------------------------------------
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t RUN_CURRENT_MA = 400;

// ---- geometry and rates ---------------------------------------------------
constexpr int64_t U256 = 256;
constexpr int32_t FULL_STEPS_PER_MM = 100;
constexpr int64_t U256_PER_MM = U256 * FULL_STEPS_PER_MM;
constexpr uint32_t MAX_PULSE_RATE_HZ = 20000;
constexpr float MIN_PULSE_RATE_HZ = 2.0F;
constexpr uint32_t ACCEL_FULLSTEPS_S2 = 4000;   // 40 mm/s^2
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;
constexpr uint32_t SAMPLE_INTERVAL_MS = 20;

// ---- homing ---------------------------------------------------------------
constexpr uint16_t HOME_FAST_FS = 200;
constexpr uint16_t HOME_SLOW_FS = 15;
constexpr uint8_t  HOME_TOUCHES = 3;
constexpr int64_t  HOME_BACKOFF_U256 = 200 * U256;      // 2 mm
constexpr int64_t  HOME_TOLERANCE_U256 = 2 * U256;      // 20 um
constexpr int64_t  HOME_MARGIN_U256 = 200 * U256;       // 2 mm
constexpr int64_t  HOME_MAX_TRAVEL_U256 = 20000 * U256; // 200 mm
constexpr uint8_t  HOME_RETRIES = 1;
constexpr uint32_t SWITCH_DEBOUNCE_US = 3000;
constexpr uint16_t HOME_MRES = 16;

// A full step sits at MSCNT 128, 384, 640, 896 - both phases equal, 45 deg
// electrical. MRES may only change there; confirm against the datasheet.
constexpr uint16_t MSCNT_FULL_STEP_OFFSET = 128;
constexpr uint16_t MSCNT_FULL_STEP_PERIOD = 256;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

int64_t posU256 = 0;
int64_t stepU256 = U256 / HOME_MRES;
uint16_t mres = HOME_MRES;
uint16_t currentMa = RUN_CURRENT_MA;
bool driverReady = false, homed = false, running = false, abortRequested = false;
int64_t limNegU256 = 0, limPosU256 = 0, swNegU256 = 0, swPosU256 = 0;
bool lastDirPositive = true, dirKnown = false;
uint32_t nextSampleMs = 0;
int32_t currentRow = -1;
String inBuf;
String pendingLine;
char pendingKey = 0;

enum StopReason : uint8_t { STOP_NONE = 0, STOP_ABORT, STOP_SWITCH,
                            STOP_LIMIT, STOP_ERROR };
StopReason stopReason = STOP_NONE;
enum GuardMode : uint8_t { GUARD_NORMAL, GUARD_SEEK, GUARD_RETRACT };

// ----------------------------------------------------------- end switches

bool switchPressed(uint8_t pin) {
  if (digitalRead(pin) != SWITCH_ACTIVE_LEVEL) return false;
  const uint32_t until = micros() + SWITCH_DEBOUNCE_US;
  while ((int32_t)(micros() - until) < 0) {
    if (digitalRead(pin) != SWITCH_ACTIVE_LEVEL) return false;
  }
  return true;
}

bool switchNeg() { return switchPressed(SW_NEG_PIN); }
bool switchPos() { return switchPressed(SW_POS_PIN); }
bool switchAhead(bool positive) { return positive ? switchPos() : switchNeg(); }

// ------------------------------------------------------------------ input

/** A menu key is one the menu acts on the instant it arrives, with no Enter.
 *  Anything else starts a line, so PING / DRIVER / HOME / ABORT still work. */
bool isMenuKey(char ch) { return (ch >= '0' && ch <= '9') || ch == '?'; }

void pumpSerial() {
  while (Serial.available()) {
    const char ch = (char)Serial.read();
    if (ch == 3) { abortRequested = true; inBuf = ""; continue; }  // Ctrl+C
    if (ch == '\r' || ch == '\n') {
      String line = inBuf;
      inBuf = "";
      line.trim();
      if (line.length()) {
        if (line.equalsIgnoreCase("ABORT")) abortRequested = true;
        else pendingLine = line;
      }
      continue;
    }
    if (inBuf.length() == 0 && isMenuKey(ch)) {
      pendingKey = ch;
      continue;
    }
    if (inBuf.length() < 180) inBuf += ch;
  }
}

String readLine(const char *prompt, uint32_t timeoutMs = 30000) {
  Serial.print(prompt);
  inBuf = "";
  pendingKey = 0;
  const uint32_t until = millis() + timeoutMs;
  while ((int32_t)(millis() - until) < 0) {
    while (Serial.available()) {
      const char ch = (char)Serial.read();
      if (ch == '\r') continue;
      if (ch == '\n') {
        Serial.println();
        String s = inBuf;
        inBuf = "";
        s.trim();
        return s;
      }
      if (ch == 27) { inBuf = ""; Serial.println(); return ""; }
      if ((ch == 8 || ch == 127) && inBuf.length()) {
        inBuf.remove(inBuf.length() - 1);
        Serial.print("\b \b");
        continue;
      }
      if (inBuf.length() < 60) { inBuf += ch; Serial.print(ch); }
    }
    delay(2);
  }
  Serial.println(" (timed out)");
  inBuf = "";
  return "";
}

// ----------------------------------------------------------------- driver

uint32_t readMscnt() { return driver.MSCNT(); }

bool mscntOnFullStep(uint32_t m) {
  return ((m + MSCNT_FULL_STEP_PERIOD - MSCNT_FULL_STEP_OFFSET)
          % MSCNT_FULL_STEP_PERIOD) == 0;
}

bool applyMres(uint16_t value, bool requireFullStep) {
  if (value < 1 || value > 256 || (value & (value - 1))) {
    Serial.println("# ERR,reason=bad_mres");
    return false;
  }
  if (requireFullStep) {
    if (posU256 % U256 != 0) {
      Serial.printf("# ERR,reason=mres_not_at_full_step,pos=%lld\n",
                    (long long)posU256);
      return false;
    }
    const uint32_t m = readMscnt();
    if (!mscntOnFullStep(m)) {
      Serial.printf("# ERR,reason=mscnt_off_full_step,mscnt=%lu\n",
                    (unsigned long)m);
      return false;
    }
  }
  uint8_t code = 8;
  for (uint16_t v = value; v > 1; v >>= 1) --code;
  const uint32_t MASK = 0x0F000000UL;
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    const uint32_t chop = driver.CHOPCONF();
    driver.CHOPCONF((chop & ~MASK) | ((uint32_t)code << 24));
    if (((driver.CHOPCONF() & MASK) >> 24) == code) {
      mres = value;
      stepU256 = U256 / value;
      return true;
    }
    delay(5);
  }
  Serial.println("# ERR,reason=mres_write_failed");
  return false;
}

void reportDriver() {
  const uint32_t g = driver.GCONF(), c = driver.CHOPCONF(), d = driver.DRV_STATUS();
  Serial.printf("# DRIVER,spreadcycle=%u,intpol=%u,stealth=%u,mres=%u,"
                "mscnt=%lu,irun=%u,ihold=%u,current_ma=%u\n",
                (g & (1UL << 2)) ? 1u : 0u, (c & (1UL << 28)) ? 1u : 0u,
                (d & (1UL << 30)) ? 1u : 0u, mres,
                (unsigned long)readMscnt(), driver.irun(), driver.ihold(),
                currentMa);
}

bool configureDriver() {
  digitalWrite(EN_PIN, HIGH);
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(100);
  const uint32_t ioin = driver.IOIN();
  Serial.printf("# DRIVER_ID,ioin=0x%08lX,version=0x%02X,expected=0x21\n",
                (unsigned long)ioin, (unsigned)((ioin >> 24) & 0xFF));
  driver.begin();
  driver.toff(5);
  driver.blank_time(24);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.rms_current(currentMa);
  driver.ihold(driver.irun());        // no standstill reduction
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);
  driver.semin(0);
  driver.TCOOLTHRS(0);
  driver.pwm_autoscale(false);

  uint8_t connection = 2;
  for (uint8_t i = 0; i < 5 && connection; ++i) {
    connection = driver.test_connection();
    if (connection) delay(50);
  }
  if (connection) { Serial.println("# ERR,reason=driver_no_reply"); return false; }
  if (!applyMres(HOME_MRES, false)) return false;

  const uint32_t g = driver.GCONF(), c = driver.CHOPCONF(), d = driver.DRV_STATUS();
  if (!(g & (1UL << 2)) || (c & (1UL << 28)) || (d & (1UL << 30))) {
    Serial.println("# ERR,reason=driver_mode_not_confirmed");
    return false;
  }
  digitalWrite(EN_PIN, LOW);
  return true;
}

// -------------------------------------------------------------- telemetry

void serviceTelemetry() {
  const uint32_t now = millis();
  if ((int32_t)(now - nextSampleMs) >= 0) {
    Serial.printf("# S,%lu,%lld,%lu,%ld\n", (unsigned long)now,
                  (long long)posU256, (unsigned long)readMscnt(),
                  (long)currentRow);
    nextSampleMs = now + SAMPLE_INTERVAL_MS;
  }
}

// ----------------------------------------------------------------- motion

bool motionGuard(bool positive, GuardMode mode) {
  if (abortRequested) { stopReason = STOP_ABORT; return false; }
  if (mode == GUARD_NORMAL) {
    if (switchNeg() || switchPos()) { stopReason = STOP_SWITCH; return false; }
    if (homed) {
      const int64_t next = posU256 + (positive ? stepU256 : -stepU256);
      if (next > limPosU256 || next < limNegU256) {
        stopReason = STOP_LIMIT;
        return false;
      }
    }
  } else if (mode == GUARD_SEEK) {
    if (switchAhead(positive)) { stopReason = STOP_SWITCH; return false; }
  }
  return true;
}

void pulseOnce(bool positive) {
  if (!dirKnown || positive != lastDirPositive) {
    digitalWrite(DIR_PIN, positive ? HIGH : LOW);
    delayMicroseconds(DIR_SETUP_US);
    lastDirPositive = positive;
    dirKnown = true;
  }
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(STEP_HIGH_US);
  digitalWrite(STEP_PIN, LOW);
  posU256 += positive ? stepU256 : -stepU256;
}

bool dwellMs(uint32_t ms) {
  const uint32_t until = millis() + ms;
  while ((int32_t)(millis() - until) < 0) {
    serviceTelemetry();
    pumpSerial();
    if (abortRequested) { stopReason = STOP_ABORT; return false; }
    delay(1);
  }
  return true;
}

bool moveToU256(int64_t target, uint16_t speedFs,
                GuardMode mode = GUARD_NORMAL) {
  const int64_t delta = target - posU256;
  if (delta == 0) return true;
  const bool positive = delta > 0;
  const int64_t magnitude = positive ? delta : -delta;
  if (magnitude % stepU256 != 0) {
    Serial.printf("# ERR,reason=target_off_grid,target=%lld,mres=%u\n",
                  (long long)target, mres);
    stopReason = STOP_ERROR;
    return false;
  }
  if (speedFs == 0 || (uint32_t)speedFs * mres > MAX_PULSE_RATE_HZ) {
    Serial.printf("# ERR,reason=pulse_rate,speed_fs=%u,mres=%u\n",
                  speedFs, mres);
    stopReason = STOP_ERROR;
    return false;
  }
  const uint32_t pulses = (uint32_t)(magnitude / stepU256);
  const float vMax = (float)speedFs * (float)mres;
  const float accel = (float)ACCEL_FULLSTEPS_S2 * (float)mres;
  uint32_t nAcc = (uint32_t)(vMax * vMax / (2.0F * accel));
  if (2 * nAcc > pulses) nAcc = pulses / 2;

  for (uint32_t i = 0; i < pulses; ++i) {
    if (!motionGuard(positive, mode)) return false;
    pulseOnce(positive);
    float v;
    if (i < nAcc)                v = sqrtf(2.0F * accel * (float)(i + 1));
    else if (i >= pulses - nAcc) v = sqrtf(2.0F * accel * (float)(pulses - i));
    else                         v = vMax;
    if (v < MIN_PULSE_RATE_HZ) v = MIN_PULSE_RATE_HZ;
    if (v > (float)MAX_PULSE_RATE_HZ) v = (float)MAX_PULSE_RATE_HZ;
    const uint32_t due = micros() + (uint32_t)(1000000.0F / v);
    while ((int32_t)(micros() - due) < 0) {
      serviceTelemetry();
      pumpSerial();
      if (abortRequested) { stopReason = STOP_ABORT; return false; }
    }
  }
  return true;
}

// ----------------------------------------------------------------- homing

bool retractFrom(bool positive, int64_t distance) {
  stopReason = STOP_NONE;
  return moveToU256(posU256 + (positive ? -distance : distance),
                    HOME_SLOW_FS, GUARD_RETRACT);
}

bool homeOneEnd(bool positive, int64_t *triggerOut) {
  const char *name = positive ? "pos" : "neg";
  for (uint8_t attempt = 0; attempt <= HOME_RETRIES; ++attempt) {
    if (switchAhead(positive) && !retractFrom(positive, HOME_BACKOFF_U256))
      return false;
    Serial.printf("# HOMING,%s,seek\n", name);
    stopReason = STOP_NONE;
    moveToU256(posU256 + (positive ? HOME_MAX_TRAVEL_U256
                                   : -HOME_MAX_TRAVEL_U256),
               HOME_FAST_FS, GUARD_SEEK);
    if (abortRequested) return false;
    if (!switchAhead(positive)) {
      Serial.printf("# ERR,reason=no_switch_within_travel,end=%s\n", name);
      return false;
    }
    int64_t touch[HOME_TOUCHES];
    bool good = true;
    for (uint8_t t = 0; t < HOME_TOUCHES; ++t) {
      if (!retractFrom(positive, HOME_BACKOFF_U256)) return false;
      if (switchAhead(positive)) {
        Serial.printf("# ERR,reason=switch_still_closed,end=%s\n", name);
        return false;
      }
      stopReason = STOP_NONE;
      moveToU256(posU256 + (positive ? HOME_BACKOFF_U256 * 3
                                     : -HOME_BACKOFF_U256 * 3),
                 HOME_SLOW_FS, GUARD_SEEK);
      if (abortRequested) return false;
      if (!switchAhead(positive)) { good = false; break; }
      touch[t] = posU256;
      Serial.printf("# HOMING,%s,touch=%u,pos=%lld\n", name,
                    (unsigned)(t + 1), (long long)touch[t]);
    }
    if (!good) continue;
    int64_t lo = touch[0], hi = touch[0];
    for (uint8_t t = 1; t < HOME_TOUCHES; ++t) {
      if (touch[t] < lo) lo = touch[t];
      if (touch[t] > hi) hi = touch[t];
    }
    Serial.printf("# HOMING,%s,spread=%lld,tolerance=%lld\n", name,
                  (long long)(hi - lo), (long long)HOME_TOLERANCE_U256);
    if (hi - lo <= HOME_TOLERANCE_U256) {
      *triggerOut = positive ? lo : hi;
      return true;
    }
  }
  Serial.printf("# ERR,reason=home_spread_out_of_tolerance,end=%s\n", name);
  return false;
}

bool homeBothEnds() {
  homed = false;
  abortRequested = false;
  if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return false; }
  if (!applyMres(HOME_MRES, false)) return false;
  Serial.println("# HOMING,start");
  int64_t negTrig = 0, posTrig = 0;
  if (!homeOneEnd(false, &negTrig)) return false;
  if (!retractFrom(false, HOME_BACKOFF_U256)) return false;
  if (!homeOneEnd(true, &posTrig)) return false;
  if (!retractFrom(true, HOME_BACKOFF_U256)) return false;

  const int64_t travel = posTrig - negTrig;
  if (travel <= 2 * HOME_MARGIN_U256) {
    Serial.printf("# ERR,reason=travel_too_short,travel=%lld\n",
                  (long long)travel);
    return false;
  }
  int64_t centre = ((negTrig + travel / 2) / U256) * U256;
  posU256 -= centre;
  swNegU256 = negTrig - centre;
  swPosU256 = posTrig - centre;
  limNegU256 = ((swNegU256 + HOME_MARGIN_U256 + U256 - 1) / U256) * U256;
  limPosU256 = ((swPosU256 - HOME_MARGIN_U256) / U256) * U256;
  homed = true;
  Serial.printf("# HOMED,travel_u256=%lld,travel_mm=%.4f,lim_neg=%lld,"
                "lim_pos=%lld\n", (long long)travel,
                (double)travel / (double)U256_PER_MM,
                (long long)limNegU256, (long long)limPosU256);
  return true;
}

/** Nudge to the nearest full step so an MRES change is legal. */
bool alignToFullStep() {
  // MSCNT decides where an electrical full step is, not our pulse count. At
  // reset MSCNT is 0, which is a phase-zero position, NOT one of the
  // both-phases-equal points at 128/384/640/896 - so stepping has to walk it
  // onto the grid before any MRES change is legal. Step bookkeeping and MSCNT
  // then advance together (both by 256 per full step), so aligning once here
  // keeps them locked for the whole run.
  const int64_t step = stepU256;
  for (uint16_t i = 0; i < 2 * MSCNT_FULL_STEP_PERIOD; ++i) {
    if (mscntOnFullStep(readMscnt())) break;
    if (!motionGuard(true, GUARD_NORMAL)) return false;
    pulseOnce(true);
    delay(3);
  }
  const uint32_t m = readMscnt();
  if (!mscntOnFullStep(m)) {
    Serial.printf("# ERR,reason=cannot_align_mscnt,mscnt=%lu,mres=%u\n",
                  (unsigned long)m, mres);
    return false;
  }
  // Re-datum so this electrical full step is a whole step in the bookkeeping.
  // The shift is under one full step (10 um) against 2 mm end margins, and
  // every measurement is a difference, so it costs nothing.
  const int64_t before = posU256;
  posU256 = ((posU256 + U256 / 2) / U256) * U256;
  Serial.printf("# ALIGNED,pos=%lld,mscnt=%lu,on_grid=1,datum_shift=%lld,"
                "step_u256=%lld\n",
                (long long)posU256, (unsigned long)m,
                (long long)(posU256 - before), (long long)step);
  return true;
}

/* Prove the microstep size really changed, not just the register.
 *
 * MSCNT is the driver's own index into the sine table - 1024 entries over one
 * electrical period, 256 per full step - so one STEP pulse must advance it by
 * exactly 256/MRES. A CHOPCONF readback only shows the register holds the
 * code; this shows the driver is stepping on the grid that code asks for. If
 * MicroPlyer (intpol) were on, every size would interpolate to 1/256 and this
 * test would fail for the coarse sizes - which is the point.
 *
 * Five pulses out and five back, so position and phase both end where they
 * started.
 */
bool mscntDeltaMatches(uint32_t before, uint32_t after, uint32_t expect) {
  const uint32_t d = (after + 1024u - before) % 1024u;
  const uint32_t e = expect % 1024u;
  return d == e || d == (1024u - e) % 1024u;   // either direction of travel
}

bool verifyMresSwitching() {
  if (!driverReady) {
    Serial.println("# ERR,reason=driver_not_ready");
    return false;
  }
  if (!alignToFullStep()) return false;

  uint16_t sizes[8];
  uint8_t n = 0;
  for (uint16_t i = 0; i < PLAN_ROWS && n < 8; ++i) {
    bool seen = false;
    for (uint8_t j = 0; j < n; ++j) seen |= (sizes[j] == PLAN[i].mres);
    if (!seen) sizes[n++] = PLAN[i].mres;
  }

  Serial.printf("\n# MRES_CHECK_START,sizes=%u\n", n);
  Serial.println(F("  MRES   per pulse   1 pulse   5 pulses   closes   verdict"));
  bool allOk = true;
  for (uint8_t i = 0; i < n; ++i) {
    const uint16_t m = sizes[i];
    if (!applyMres(m, true)) { allOk = false; continue; }
    const uint32_t expect = MSCNT_FULL_STEP_PERIOD / m;
    const uint32_t before = readMscnt();
    pulseOnce(true);
    delay(5);
    const uint32_t after1 = readMscnt();
    for (uint8_t k = 0; k < 4; ++k) { pulseOnce(true); delay(5); }
    const uint32_t after5 = readMscnt();
    for (uint8_t k = 0; k < 5; ++k) { pulseOnce(false); delay(5); }
    const uint32_t back = readMscnt();

    const bool ok1 = mscntDeltaMatches(before, after1, expect);
    const bool ok5 = mscntDeltaMatches(before, after5, 5 * expect);
    const bool okBack = (back == before);
    const bool ok = ok1 && ok5 && okBack;
    allOk &= ok;
    Serial.printf("  1/%-4u %9lu %9lu %10lu %8s   %s\n", m,
                  (unsigned long)expect,
                  (unsigned long)((after1 + 1024u - before) % 1024u),
                  (unsigned long)((after5 + 1024u - before) % 1024u),
                  okBack ? "yes" : "NO", ok ? "OK" : "MISMATCH");
    if (!ok) {
      Serial.printf("# MRES_CHECK_FAIL,mres=%u,expect=%lu,before=%lu,"
                    "after1=%lu,after5=%lu,back=%lu\n", m,
                    (unsigned long)expect, (unsigned long)before,
                    (unsigned long)after1, (unsigned long)after5,
                    (unsigned long)back);
    }
  }
  Serial.printf("# MRES_CHECK,%s\n", allOk ? "pass" : "FAIL");
  if (!allOk) {
    Serial.println("Microstep size is not following the register. Check "
                   "GCONF.mstep_reg_select and that intpol is off before "
                   "running the matrix - the comparison is meaningless "
                   "otherwise.");
  }
  return allOk;
}

// ------------------------------------------------------------- the matrix

bool executeRow(const PlanRow &r, uint16_t index) {
  currentRow = index;
  if (r.kind == PLAN_SET_MRES) {
    if (!applyMres(r.mres, true)) return false;
    Serial.printf("# MRES,%u,mscnt=%lu\n", mres, (unsigned long)readMscnt());
  } else {
    if (r.mres != mres && !applyMres(r.mres, true)) return false;
    Serial.printf("# MOVE,%u,%s,target=%ld,d=%ld,dir=%d,mres=%u\n", index,
                  r.tag, (long)r.target_u256, (long)r.d_u256,
                  (int)r.direction, r.mres);
    if (!moveToU256(r.target_u256, r.speed_fs)) return false;
  }
  if (!dwellMs(r.dwell_ms)) return false;
  // The dwell has ended: this is the settled stop the analysis reads.
  Serial.printf("# STOP,%u,%s,d=%ld,dir=%d,measured=%u,pos=%lld,mscnt=%lu\n",
                index, r.tag, (long)r.d_u256, (int)r.direction,
                (unsigned)r.measured, (long long)posU256,
                (unsigned long)readMscnt());
  return true;
}

/* The plan was written against a nominal travel. Check it fits the travel
 * this stage actually has before moving, rather than failing on STOP_LIMIT
 * two thirds of the way through a 25 minute run. */
bool planFitsTravel() {
  int64_t lo = 0, hi = 0;
  for (uint16_t i = 0; i < PLAN_ROWS; ++i) {
    if (PLAN[i].kind != PLAN_MOVE) continue;
    const int64_t t = PLAN[i].target_u256;
    if (t < lo) lo = t;
    if (t > hi) hi = t;
  }
  Serial.printf("# PLAN_EXTENT,min_mm=%.3f,max_mm=%.3f,limits_mm=%.3f..%.3f\n",
                (double)lo / (double)U256_PER_MM,
                (double)hi / (double)U256_PER_MM,
                (double)limNegU256 / (double)U256_PER_MM,
                (double)limPosU256 / (double)U256_PER_MM);
  if (lo < limNegU256 || hi > limPosU256) {
    Serial.printf("# ERR,reason=plan_exceeds_travel,short_by_neg_mm=%.3f,"
                  "short_by_pos_mm=%.3f\n",
                  lo < limNegU256
                      ? (double)(limNegU256 - lo) / (double)U256_PER_MM : 0.0,
                  hi > limPosU256
                      ? (double)(hi - limPosU256) / (double)U256_PER_MM : 0.0);
    return false;
  }
  return true;
}

void runMatrix() {
  if (!homed) {
    Serial.println("\nThe matrix needs the travel measured first. Run 1 or 2.");
    return;
  }
  if (!planFitsTravel()) {
    Serial.println("The plan does not fit this travel. Shorten the long "
                   "distances in scripts/settling_design.py and regenerate "
                   "the header.");
    return;
  }
  Serial.printf("\nSettling matrix: %u rows, %u measured displacements, "
                "about %lu min.\nThe extent printed above is where the stage "
                "will travel. Type RUN to start.\n",
                PLAN_ROWS, PLAN_MEASURED,
                (unsigned long)((PLAN_SECONDS + 30) / 60));
  if (!readLine("> ").equalsIgnoreCase("RUN")) {
    Serial.println("Cancelled.");
    return;
  }
  if (!moveToU256(0, HOME_SLOW_FS)) {
    Serial.println("# ERR,reason=x_start_failed");
    return;
  }
  // Prove the microstep size follows the register before spending 32 min
  // comparing microstep sizes. verifyMresSwitching() aligns MSCNT itself.
  if (!verifyMresSwitching()) {
    Serial.println("Not running the matrix: microstep switching is not "
                   "verified.");
    return;
  }
  running = true;
  abortRequested = false;
  stopReason = STOP_NONE;
  // Anything typed during the run is noise, not a queued menu choice.
  pendingLine = "";
  pendingKey = 0;
  digitalWrite(LED_PIN, HIGH);
  Serial.printf("# RUN_START,rows=%u,measured=%u,pos=%lld\n",
                PLAN_ROWS, PLAN_MEASURED, (long long)posU256);
  bool ok = true;
  for (uint16_t i = 0; i < PLAN_ROWS && !abortRequested; ++i) {
    if (!executeRow(PLAN[i], i)) { ok = false; break; }
  }
  digitalWrite(LED_PIN, LOW);
  currentRow = -1;
  running = false;
  pendingLine = "";
  pendingKey = 0;
  if (abortRequested) Serial.printf("# RUN_ABORTED,pos=%lld\n", (long long)posU256);
  else if (!ok)       Serial.printf("# RUN_FAILED,pos=%lld,reason=%u\n",
                                    (long long)posU256, (unsigned)stopReason);
  else                Serial.printf("# RUN_COMPLETE,pos=%lld\n", (long long)posU256);
}

// -------------------------------------------------------------------- menu

void printMenu() {
  Serial.println(F("\n================================================================"));
  Serial.println(F(" Settling matrix across microstep sizes"));
  Serial.printf(  " driver=%u homed=%u   1/%u   %u mA\n",
                driverReady ? 1u : 0u, homed ? 1u : 0u, mres, currentMa);
  if (homed) {
    Serial.printf(" travel %.4f mm   limits %.3f .. %.3f mm   pos %.4f mm\n",
                  (double)(swPosU256 - swNegU256) / (double)U256_PER_MM,
                  (double)limNegU256 / (double)U256_PER_MM,
                  (double)limPosU256 / (double)U256_PER_MM,
                  (double)posU256 / (double)U256_PER_MM);
  } else {
    Serial.println(F(" travel not measured yet - home before running"));
  }
  Serial.println(F("================================================================"));
  Serial.printf (  "  1  Home + run matrix   home, then %u rows (~%lu min)\n",
                PLAN_ROWS, (unsigned long)((PLAN_SECONDS + 30) / 60));
  Serial.println(F("  2  Home only"));
  Serial.println(F("  3  Run matrix          uses the travel already measured"));
  Serial.println(F("  4  Go to centre        x_start, on a full step"));
  Serial.println(F("  5  Driver readback"));
  Serial.println(F("  7  Verify MRES switching   MSCNT per pulse, no travel"));
  Serial.println(F("  6  Switch state"));
  Serial.println(F("  0  Release motor       EN high, outputs off"));
  Serial.println(F("  ?  Reprint this menu"));
  Serial.print  (F("\n> "));
}

void dispatchMenu(char key) {
  switch (key) {
    case '?': printMenu(); return;
    case '1': if (homeBothEnds()) runMatrix(); break;
    case '2': if (homeBothEnds()) Serial.println("Homed."); break;
    case '3': runMatrix(); break;
    case '4':
      if (homed && moveToU256(0, HOME_SLOW_FS)) alignToFullStep();
      else Serial.println("# ERR,reason=not_homed");
      break;
    case '5': reportDriver(); break;
    case '7': verifyMresSwitching(); break;
    case '6':
      Serial.printf("# SWITCHES,neg=%u,pos=%u\n", switchNeg() ? 1u : 0u,
                    switchPos() ? 1u : 0u);
      break;
    case '0':
      digitalWrite(EN_PIN, HIGH);
      homed = false;
      Serial.println("Outputs off. Travel invalidated - home again.");
      break;
    default: Serial.println("Unknown key."); break;
  }
  printMenu();
}

// ------------------------------------------------------------------- setup

void setup() {
  pinMode(EN_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  pinMode(LED_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, HIGH);
  digitalWrite(LED_PIN, LOW);
  pinMode(SW_NEG_PIN, INPUT_PULLUP);
  pinMode(SW_POS_PIN, INPUT_PULLUP);

  Serial.begin(115200);
  const uint32_t start = millis();
  while (!Serial && millis() - start < 3000) delay(10);
  delay(200);

  Serial.println("# BOOT,fw=settling_test");
  Serial.printf("# PINS,step=%u,dir=%u,en=%u,uart_rx=%u,uart_tx=%u,"
                "sw_neg=%u,sw_pos=%u\n",
                STEP_PIN, DIR_PIN, EN_PIN, UART_RX_PIN, UART_TX_PIN,
                SW_NEG_PIN, SW_POS_PIN);
  Serial.printf("# SWITCHES,neg=%u,pos=%u,active_level=%s\n",
                switchNeg() ? 1u : 0u, switchPos() ? 1u : 0u,
                SWITCH_ACTIVE_LEVEL == LOW ? "low" : "high");
  driverReady = configureDriver();
  if (!driverReady) Serial.println("# ERR,reason=driver_config_failed");
  reportDriver();
  Serial.printf("# READY,driver=%u,rows=%u,measured=%u\n",
                driverReady ? 1u : 0u, PLAN_ROWS, PLAN_MEASURED);
  nextSampleMs = millis();
  printMenu();
}

void loop() {
  pumpSerial();
  if (abortRequested && !running) {
    abortRequested = false;
    Serial.println("# ABORT");
  }
  if (pendingKey) {
    const char key = pendingKey;
    pendingKey = 0;
    Serial.println(key);
    dispatchMenu(key);
  }
  if (pendingLine.length()) {
    const String line = pendingLine;
    pendingLine = "";
    if (line.equalsIgnoreCase("PING")) {
      Serial.printf("# PONG,fw=settling_test,driver=%u,homed=%u,rows=%u\n",
                    driverReady ? 1u : 0u, homed ? 1u : 0u, PLAN_ROWS);
    } else if (line.equalsIgnoreCase("DRIVER")) {
      reportDriver();
    } else if (line.equalsIgnoreCase("HOME")) {
      homeBothEnds();
    } else {
      Serial.println("# ERR,reason=unknown_command");
      Serial.print("> ");
    }
  }
  delay(2);
}

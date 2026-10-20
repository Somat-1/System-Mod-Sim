/*
 * MRES test runner - ESP32-S3 + TMC2209 + quadrature encoder + two end switches.
 *
 * NOT YET COMPILED OR RUN ON HARDWARE. No toolchain was available where this
 * was written. Compile, then bring it up with the motor uncoupled.
 *
 * Two parts, picked from an interactive menu on the USB CDC terminal:
 *
 *   Part 1  Home + warm-up. Learns both ends of travel by touching each switch
 *           several times, puts x_start at the centre of what it measured, then
 *           shuttles continuously until stopped. Current, speed, MRES and
 *           amplitude are settable live while it runs.
 *
 *   Part 2  The microstep-resolution sequence. Homes to the centre, then runs
 *           the plan compiled into mres_plan.h standalone - no host needed.
 *
 *   Part 3  Host plan mode. The PLAN ADD / PLAN RUN protocol of ../README.md,
 *           so scripts/run_mres.py can drive a plan without reflashing. Same
 *           executor as part 2, different source of rows.
 *
 * POSITIONS are int64 in u256 = 1/256 full step. 100 full steps = 1 mm, so
 * 1 u256 = 39.06 nm. After homing, 0 is the centre of the measured travel,
 * snapped to a full step, and both soft limits are relative to it.
 *
 * WIRING - the Stepper Isolated bench (see ../../../../CONNECTIONS.md):
 *
 *   STEP GPIO7   DIR GPIO8   EN GPIO9 (ENN, active low)
 *   TMC UART TX GPIO17 (via 1k) / RX GPIO18
 *   ENCODER  A GPIO4   B GPIO5   Z GPIO6
 *   END SWITCH neg GPIO15   pos GPIO16     <-- PLACEHOLDER, set before use
 *
 * The end switches are 10T85 microswitches with COM on common ground and NO to
 * the GPIO, so the pins use the internal pull-up and read LOW when pressed.
 * This is not fail-safe: a broken switch wire reads as "clear" and the stage
 * would drive into the hard stop. Moving each wire to the NC contact and
 * setting SWITCH_ACTIVE_LEVEL to HIGH would make a cut cable stop the stage.
 *
 * Nothing moves until a mode is chosen. Both switches are checked before every
 * pulse in every mode; outside homing a trigger stops motion immediately.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/pulse_cnt.h"
#include "mres_plan.h"

// ---- pins -----------------------------------------------------------------
constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;            // TMC ENN: HIGH = outputs off
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr gpio_num_t ENC_A = GPIO_NUM_4;
constexpr gpio_num_t ENC_B = GPIO_NUM_5;
constexpr gpio_num_t ENC_Z = GPIO_NUM_6;

// PLACEHOLDER pin numbers - confirm against the bench before the first run.
constexpr uint8_t SW_NEG_PIN = 15;       // switch at the negative end
constexpr uint8_t SW_POS_PIN = 16;       // switch at the positive end
constexpr int SWITCH_ACTIVE_LEVEL = LOW; // NO to GND + internal pull-up

// ---- driver ---------------------------------------------------------------
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t DEFAULT_CURRENT_MA = 360;

// ---- geometry and rates ---------------------------------------------------
constexpr int64_t U256 = 256;                  // u256 per full step
constexpr int32_t FULL_STEPS_PER_MM = 100;
constexpr int64_t U256_PER_MM = U256 * FULL_STEPS_PER_MM;
constexpr uint32_t MAX_PULSE_RATE_HZ = 20000;
constexpr float MIN_PULSE_RATE_HZ = 2.0F;
constexpr uint32_t ACCEL_FULLSTEPS_S2 = 4000;  // 40 mm/s^2
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;
constexpr uint32_t STAIR_PULSE_GAP_US = 200;
constexpr uint32_t SAMPLE_INTERVAL_MS = 20;    // 50 Hz telemetry
constexpr size_t MAX_UPLOADED_ROWS = 512;

// ---- homing ---------------------------------------------------------------
// One fast seek finds each switch, then HOME_TOUCHES slow touches measure it.
// The slow touches must agree within HOME_TOLERANCE_U256 or homing fails: a
// bouncing or soft contact shows up as a spread, and accepting it silently
// would put the soft limit in the wrong place.
constexpr uint16_t HOME_FAST_FS = 200;                 // 2 mm/s seek
constexpr uint16_t HOME_SLOW_FS = 15;                  // 0.15 mm/s touch
constexpr uint8_t  HOME_TOUCHES = 3;
constexpr int64_t  HOME_BACKOFF_U256 = 50 * U256;      // 0.5 mm
constexpr int64_t  HOME_TOLERANCE_U256 = 2 * U256;     // 20 um
constexpr int64_t  HOME_MARGIN_U256 = 50 * U256;       // 0.5 mm inside the switch
constexpr int64_t  HOME_MAX_TRAVEL_U256 = 6000 * U256; // 60 mm before giving up
constexpr uint8_t  HOME_RETRIES = 1;
constexpr uint32_t SWITCH_DEBOUNCE_US = 3000;
constexpr uint16_t HOME_MRES = 16;                     // homing runs at 1/16

// ---- warm-up --------------------------------------------------------------
constexpr uint16_t WARMUP_SPEED_FS = 70;       // the rig's MODERATE, 0.7 mm/s
constexpr uint32_t WARMUP_END_DWELL_MS = 200;

// The TMC sine table puts a full step at MSCNT 128, 384, 640 and 896 - both
// phases equal, 45 deg electrical. MRES may only change there. Confirm against
// the datasheet for your table before trusting a block that changes MRES.
constexpr uint16_t MSCNT_FULL_STEP_OFFSET = 128;
constexpr uint16_t MSCNT_FULL_STEP_PERIOD = 256;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);
pcnt_unit_handle_t pcnt = nullptr;

// ---- state ----------------------------------------------------------------
struct Settings {
  uint16_t currentMa = DEFAULT_CURRENT_MA;
  uint16_t mres = HOME_MRES;
  uint16_t warmupSpeedFs = WARMUP_SPEED_FS;
  int64_t warmupAmplitudeU256 = 0;   // peak-to-peak; 0 = full travel
  float iholdMultiplier = 1.0F;
  uint8_t tpowerdown = 0;
  uint8_t iholddelay = 0;
  bool intpol = false;
} cfg;

int64_t posU256 = 0;
int64_t stepU256 = U256 / HOME_MRES;
bool driverReady = false, encoderReady = false;
bool homed = false, running = false, abortRequested = false;
int64_t limNegU256 = 0, limPosU256 = 0;        // soft limits, centre = 0
int64_t swNegU256 = 0, swPosU256 = 0;          // where each switch triggered
bool lastDirPositive = true, dirKnown = false;
uint32_t nextSampleMs = 0;
int32_t currentRow = -1;
volatile int32_t indexCount = 0;

PlanRow uploaded[MAX_UPLOADED_ROWS];
size_t uploadedCount = 0;

String inBuf;
String pendingLine;
char pendingKey = 0;
bool hotkeys = false;          // true while the warm-up is running

enum StopReason : uint8_t {
  STOP_NONE = 0, STOP_ABORT, STOP_SWITCH, STOP_LIMIT, STOP_KEY, STOP_ERROR
};
StopReason stopReason = STOP_NONE;

/** How the switches are treated for one move. */
enum GuardMode : uint8_t {
  GUARD_NORMAL,    // either switch is a fault: stop
  GUARD_SEEK,      // the switch ahead is the goal: stop, no fault
  GUARD_RETRACT    // moving off a pressed switch: ignore both
};

// ---------------------------------------------------------------- encoder

void IRAM_ATTR onIndex() { indexCount++; }

int32_t encoderCount() {
  int value = 0;
  if (pcnt) pcnt_unit_get_count(pcnt, &value);
  return value;
}

bool setupEncoder() {
  pcnt_unit_config_t uc = {};
  uc.high_limit = 30000;
  uc.low_limit = -30000;
  uc.flags.accum_count = 1;
  if (pcnt_new_unit(&uc, &pcnt) != ESP_OK) return false;
  pcnt_glitch_filter_config_t filt = {};
  filt.max_glitch_ns = 1000;
  pcnt_unit_set_glitch_filter(pcnt, &filt);

  pcnt_chan_config_t a = {}, b = {};
  a.edge_gpio_num = ENC_A; a.level_gpio_num = ENC_B;
  b.edge_gpio_num = ENC_B; b.level_gpio_num = ENC_A;
  pcnt_channel_handle_t cha = nullptr, chb = nullptr;
  if (pcnt_new_channel(pcnt, &a, &cha) != ESP_OK) return false;
  if (pcnt_new_channel(pcnt, &b, &chb) != ESP_OK) return false;
  pcnt_channel_set_edge_action(cha, PCNT_CHANNEL_EDGE_ACTION_DECREASE,
                               PCNT_CHANNEL_EDGE_ACTION_INCREASE);
  pcnt_channel_set_level_action(cha, PCNT_CHANNEL_LEVEL_ACTION_KEEP,
                                PCNT_CHANNEL_LEVEL_ACTION_INVERSE);
  pcnt_channel_set_edge_action(chb, PCNT_CHANNEL_EDGE_ACTION_INCREASE,
                               PCNT_CHANNEL_EDGE_ACTION_DECREASE);
  pcnt_channel_set_level_action(chb, PCNT_CHANNEL_LEVEL_ACTION_KEEP,
                                PCNT_CHANNEL_LEVEL_ACTION_INVERSE);
  gpio_pullup_en(ENC_A);
  gpio_pullup_en(ENC_B);
  // accum_count only carries the count past the limits if they are watch
  // points; without these the count silently restarts every 30000 counts.
  pcnt_unit_add_watch_point(pcnt, 30000);
  pcnt_unit_add_watch_point(pcnt, -30000);
  pcnt_unit_enable(pcnt);
  pcnt_unit_clear_count(pcnt);
  pcnt_unit_start(pcnt);
  pinMode(ENC_Z, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(ENC_Z), onIndex, RISING);
  return true;
}

// ----------------------------------------------------------- end switches

/** Debounced: the pin must hold the active level for the whole window. */
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

void pumpSerial() {
  while (Serial.available()) {
    const char ch = (char)Serial.read();
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
    if (inBuf.length() < 180) inBuf += ch;
    // While the warm-up runs, a lone printable character is a hotkey, so a
    // plain serial monitor that needs Enter and a raw terminal both work.
    if (hotkeys && inBuf.length() == 1) {
      pendingKey = inBuf[0];
      inBuf = "";
    }
  }
}

char takeKey() {
  pumpSerial();
  const char k = pendingKey;
  pendingKey = 0;
  return k;
}

/** Blocking prompt. ESC cancels and returns "". Hotkeys are off meanwhile. */
String readLine(const char *prompt, uint32_t timeoutMs = 30000) {
  const bool saved = hotkeys;
  hotkeys = false;
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
        hotkeys = saved;
        s.trim();
        return s;
      }
      if (ch == 27) {
        inBuf = "";
        Serial.println(" (cancelled)");
        hotkeys = saved;
        return "";
      }
      if (ch == 8 || ch == 127) {
        if (inBuf.length()) {
          inBuf.remove(inBuf.length() - 1);
          Serial.print("\b \b");
        }
        continue;
      }
      if (inBuf.length() < 60) {
        inBuf += ch;
        Serial.print(ch);
      }
    }
    delay(2);
  }
  Serial.println(" (timed out)");
  inBuf = "";
  hotkeys = saved;
  return "";
}

// ----------------------------------------------------------------- driver

uint32_t readMscnt() { return driver.MSCNT(); }

bool mscntOnFullStep(uint32_t mscnt) {
  return ((mscnt + MSCNT_FULL_STEP_PERIOD - MSCNT_FULL_STEP_OFFSET)
          % MSCNT_FULL_STEP_PERIOD) == 0;
}

void applyCurrent(uint16_t ma) {
  driver.rms_current(ma);
  const uint8_t irun = driver.irun();
  long ihold = lroundf(irun * cfg.iholdMultiplier);
  if (ihold < 0) ihold = 0;
  if (ihold > 31) ihold = 31;
  driver.ihold((uint8_t)ihold);
  cfg.currentMa = ma;
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
    const uint32_t back = driver.CHOPCONF();
    if (((back & MASK) >> 24) == code) {
      cfg.mres = value;
      stepU256 = U256 / value;
      return true;
    }
    delay(5);
  }
  Serial.println("# ERR,reason=mres_write_failed");
  return false;
}

void reportDriver() {
  const uint32_t g = driver.GCONF();
  const uint32_t c = driver.CHOPCONF();
  const uint32_t d = driver.DRV_STATUS();
  Serial.printf("# DRIVER,spreadcycle=%u,intpol=%u,stealth=%u,ihold_mult=%.2f,"
                "tpowerdown=%u,iholddelay=%u,mres=%u,mscnt=%lu,irun=%u,ihold=%u,"
                "current_ma=%u\n",
                (g & (1UL << 2)) ? 1u : 0u, (c & (1UL << 28)) ? 1u : 0u,
                (d & (1UL << 30)) ? 1u : 0u, cfg.iholdMultiplier,
                cfg.tpowerdown, cfg.iholddelay, cfg.mres,
                (unsigned long)readMscnt(), driver.irun(), driver.ihold(),
                cfg.currentMa);
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
  applyCurrent(cfg.currentMa);
  driver.iholddelay(cfg.iholddelay);
  driver.TPOWERDOWN(cfg.tpowerdown);
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(cfg.intpol);
  driver.semin(0);
  driver.TCOOLTHRS(0);
  driver.pwm_autoscale(false);

  uint8_t connection = 2;
  for (uint8_t i = 0; i < 5 && connection; ++i) {
    connection = driver.test_connection();
    if (connection) delay(50);
  }
  if (connection) {
    Serial.println("# ERR,reason=driver_no_reply");
    return false;
  }
  if (!applyMres(cfg.mres, false)) return false;

  // SpreadCycle on, StealthChop inactive, and intpol as asked for.
  const uint32_t g = driver.GCONF();
  const uint32_t c = driver.CHOPCONF();
  const uint32_t d = driver.DRV_STATUS();
  const bool spread = (g & (1UL << 2)) != 0;
  const bool intpol = (c & (1UL << 28)) != 0;
  const bool stealth = (d & (1UL << 30)) != 0;
  if (!spread || stealth || intpol != cfg.intpol) {
    Serial.printf("# ERR,reason=driver_mode_not_confirmed,spreadcycle=%u,"
                  "stealth=%u,intpol=%u\n", spread, stealth, intpol);
    return false;
  }
  digitalWrite(EN_PIN, LOW);
  return true;
}

// -------------------------------------------------------------- telemetry

void emitSample() {
  Serial.printf("# S,%lu,%lld,%ld,%lu,%ld\n", (unsigned long)millis(),
                (long long)posU256, (long)encoderCount(),
                (unsigned long)readMscnt(), (long)currentRow);
}

void serviceTelemetry() {
  const uint32_t now = millis();
  if ((int32_t)(now - nextSampleMs) >= 0) {
    emitSample();
    nextSampleMs = now + SAMPLE_INTERVAL_MS;
  }
}

// ----------------------------------------------------------------- motion

/** Checked before every pulse. */
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

/** Wait, streaming telemetry and watching for ABORT and (optionally) a key. */
bool dwellMs(uint32_t ms, bool watchKeys = false) {
  const uint32_t until = millis() + ms;
  while ((int32_t)(millis() - until) < 0) {
    serviceTelemetry();
    pumpSerial();
    if (abortRequested) { stopReason = STOP_ABORT; return false; }
    if (watchKeys && pendingKey) { stopReason = STOP_KEY; return false; }
    delay(1);
  }
  return true;
}

/**
 * Trapezoidal move to an absolute u256 target at the current MRES.
 *
 * The target must sit on the current microstep grid. The plan guarantees it,
 * and a mismatch is a design error rather than something to round away.
 */
bool moveToU256(int64_t target, uint16_t speedFs,
                GuardMode mode = GUARD_NORMAL, bool watchKeys = false) {
  const int64_t delta = target - posU256;
  if (delta == 0) return true;
  const bool positive = delta > 0;
  const int64_t magnitude = positive ? delta : -delta;
  if (magnitude % stepU256 != 0) {
    Serial.printf("# ERR,reason=target_off_grid,target=%lld,mres=%u\n",
                  (long long)target, cfg.mres);
    stopReason = STOP_ERROR;
    return false;
  }
  if (speedFs == 0 || (uint32_t)speedFs * cfg.mres > MAX_PULSE_RATE_HZ) {
    Serial.printf("# ERR,reason=pulse_rate,speed_fs=%u,mres=%u\n",
                  speedFs, cfg.mres);
    stopReason = STOP_ERROR;
    return false;
  }
  const uint32_t pulses = (uint32_t)(magnitude / stepU256);

  const float vMax = (float)speedFs * (float)cfg.mres;
  const float accel = (float)ACCEL_FULLSTEPS_S2 * (float)cfg.mres;
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
    const uint32_t interval = (uint32_t)(1000000.0F / v);

    const uint32_t due = micros() + interval;
    while ((int32_t)(micros() - due) < 0) {
      serviceTelemetry();
      pumpSerial();
      if (abortRequested) { stopReason = STOP_ABORT; return false; }
      if (watchKeys && pendingKey) { stopReason = STOP_KEY; return false; }
    }
  }
  return true;
}

/** Staircase: n_inc increments of inc_usteps microsteps, dwelling after each. */
bool stairToU256(int64_t target, uint16_t nInc, uint16_t incUsteps,
                 uint32_t incDwellMs) {
  const int64_t delta = target - posU256;
  const int64_t perInc = (int64_t)incUsteps * stepU256;
  const bool positive = delta > 0;
  const int64_t magnitude = positive ? delta : -delta;
  if (nInc == 0 || incUsteps == 0 || magnitude != (int64_t)nInc * perInc) {
    Serial.printf("# ERR,reason=stair_mismatch,delta=%lld,n_inc=%u,"
                  "inc_usteps=%u,mres=%u\n",
                  (long long)delta, nInc, incUsteps, cfg.mres);
    stopReason = STOP_ERROR;
    return false;
  }
  for (uint16_t k = 0; k < nInc; ++k) {
    for (uint16_t p = 0; p < incUsteps; ++p) {
      if (!motionGuard(positive, GUARD_NORMAL)) return false;
      pulseOnce(positive);
      delayMicroseconds(STAIR_PULSE_GAP_US);
    }
    if (!dwellMs(incDwellMs)) return false;
  }
  return true;
}

// ----------------------------------------------------------------- homing

/** Move clear of a switch we are standing on. */
bool retractFrom(bool positive, int64_t distance) {
  const int64_t away = posU256 + (positive ? -distance : distance);
  stopReason = STOP_NONE;
  return moveToU256(away, HOME_SLOW_FS, GUARD_RETRACT);
}

/**
 * Touch one end repeatedly and return the position its switch triggers at.
 *
 * Returns the most conservative of the slow touches - the one that stopped
 * soonest - so the soft limit never ends up outside any measured trigger.
 */
bool homeOneEnd(bool positive, int64_t *triggerOut) {
  const char *name = positive ? "pos" : "neg";
  for (uint8_t attempt = 0; attempt <= HOME_RETRIES; ++attempt) {
    if (switchAhead(positive) && !retractFrom(positive, HOME_BACKOFF_U256))
      return false;

    Serial.printf("# HOMING,%s,seek\n", name);
    stopReason = STOP_NONE;
    const int64_t far = posU256 + (positive ? HOME_MAX_TRAVEL_U256
                                            : -HOME_MAX_TRAVEL_U256);
    moveToU256(far, HOME_FAST_FS, GUARD_SEEK);
    if (abortRequested) { stopReason = STOP_ABORT; return false; }
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
      const int64_t creep = posU256 + (positive ? HOME_BACKOFF_U256 * 3
                                                : -HOME_BACKOFF_U256 * 3);
      moveToU256(creep, HOME_SLOW_FS, GUARD_SEEK);
      if (abortRequested) { stopReason = STOP_ABORT; return false; }
      if (!switchAhead(positive)) { good = false; break; }
      touch[t] = posU256;
      Serial.printf("# HOMING,%s,touch=%u,pos=%lld\n", name,
                    (unsigned)(t + 1), (long long)touch[t]);
    }
    if (!good) {
      Serial.printf("# HOMING,%s,missed_touch,retry\n", name);
      continue;
    }

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
    Serial.printf("# HOMING,%s,spread_too_wide,retry\n", name);
  }
  Serial.printf("# ERR,reason=home_spread_out_of_tolerance,end=%s\n", name);
  return false;
}

/** Home both ends, put 0 at the centre and set the soft limits inside. */
bool homeBothEnds() {
  homed = false;
  abortRequested = false;
  stopReason = STOP_NONE;
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

  // Put 0 at the centre, snapped to a full step so part 2 can start there.
  int64_t centre = negTrig + travel / 2;
  centre = (centre / U256) * U256;
  posU256 -= centre;
  swNegU256 = negTrig - centre;
  swPosU256 = posTrig - centre;
  // Soft limits a margin inside each switch, themselves on a full step.
  limNegU256 = ((swNegU256 + HOME_MARGIN_U256 + U256 - 1) / U256) * U256;
  limPosU256 = ((swPosU256 - HOME_MARGIN_U256) / U256) * U256;
  homed = true;
  if (cfg.warmupAmplitudeU256 == 0)
    cfg.warmupAmplitudeU256 = limPosU256 - limNegU256;

  Serial.printf("# HOMED,travel_u256=%lld,travel_mm=%.4f,sw_neg=%lld,sw_pos=%lld,"
                "lim_neg=%lld,lim_pos=%lld\n",
                (long long)travel, (double)travel / (double)U256_PER_MM,
                (long long)swNegU256, (long long)swPosU256,
                (long long)limNegU256, (long long)limPosU256);
  Serial.println("# HOMING,done");
  return true;
}

bool goToCentre() {
  if (!homed) { Serial.println("# ERR,reason=not_homed"); return false; }
  if (!applyMres(HOME_MRES, false)) return false;
  if (!moveToU256(0, cfg.warmupSpeedFs)) return false;
  return dwellMs(500);
}

// ---------------------------------------------------------------- warm-up

void printWarmupHelp() {
  Serial.println(F("\n  keys:  p pause/resume   s stop   c current   v speed   "
                   "m MRES   a amplitude   ? status"));
}

/** Continuous shuttle between the soft limits. Runs until stopped. */
void runWarmup() {
  if (!homed) { Serial.println("Not homed - run homing first."); return; }
  Serial.printf("\nWarm-up: %.3f mm peak-to-peak at %u FS/s, 1/%u, %u mA\n",
                (double)cfg.warmupAmplitudeU256 / (double)U256_PER_MM,
                cfg.warmupSpeedFs, cfg.mres, cfg.currentMa);
  printWarmupHelp();

  bool paused = false, positive = true;
  uint32_t cycles = 0;
  const uint32_t started = millis();
  running = true;
  hotkeys = true;
  abortRequested = false;

  while (running) {
    if (!paused) {
      // Endpoints snap to a full step so an MRES change is legal there.
      int64_t half = (cfg.warmupAmplitudeU256 / 2 / U256) * U256;
      int64_t target = positive ? half : -half;
      if (target > limPosU256) target = (limPosU256 / U256) * U256;
      if (target < limNegU256) target = ((limNegU256 + U256 - 1) / U256) * U256;

      stopReason = STOP_NONE;
      const bool ok = moveToU256(target, cfg.warmupSpeedFs, GUARD_NORMAL, true);
      if (!ok && stopReason != STOP_KEY) {
        Serial.printf("\n# WARMUP_STOPPED,reason=%u,pos=%lld\n",
                      (unsigned)stopReason, (long long)posU256);
        break;
      }
      if (ok) {
        positive = !positive;
        if (positive) ++cycles;
        dwellMs(WARMUP_END_DWELL_MS, true);
      }
    } else {
      dwellMs(50, true);
    }

    const char key = takeKey();
    if (!key) continue;
    switch (key) {
      case 'p':
      case ' ':
        paused = !paused;
        Serial.printf("\n  %s at %.4f mm\n", paused ? "paused" : "resumed",
                      (double)posU256 / (double)U256_PER_MM);
        break;
      case 's':
      case 'q':
        running = false;
        break;
      case '?':
        Serial.printf("\n  cycles=%lu elapsed=%lus pos=%.4f mm enc=%ld "
                      "1/%u %u mA %u FS/s amp=%.3f mm\n",
                      (unsigned long)cycles,
                      (unsigned long)((millis() - started) / 1000),
                      (double)posU256 / (double)U256_PER_MM,
                      (long)encoderCount(), cfg.mres, cfg.currentMa,
                      cfg.warmupSpeedFs,
                      (double)cfg.warmupAmplitudeU256 / (double)U256_PER_MM);
        break;
      case 'c': {
        const long v = readLine("\n  run current mA (50-1200): ").toInt();
        if (v >= 50 && v <= 1200) {
          applyCurrent((uint16_t)v);
          Serial.printf("  current = %u mA\n", cfg.currentMa);
        } else {
          Serial.println("  unchanged");
        }
        break;
      }
      case 'v': {
        const long v = readLine("\n  speed full-step/s (1-800): ").toInt();
        if (v >= 1 && v <= 800 && (uint32_t)v * cfg.mres <= MAX_PULSE_RATE_HZ) {
          cfg.warmupSpeedFs = (uint16_t)v;
          Serial.printf("  speed = %u FS/s\n", cfg.warmupSpeedFs);
        } else {
          Serial.printf("  rejected: %ld x 1/%u would exceed %lu Hz\n",
                        v, cfg.mres, (unsigned long)MAX_PULSE_RATE_HZ);
        }
        break;
      }
      case 'm': {
        const long v = readLine("\n  MRES (1,2,4,...,256): ").toInt();
        if (v < 1 || v > 256 || (v & (v - 1))) {
          Serial.println("  rejected: not a power of two");
        } else if ((uint32_t)cfg.warmupSpeedFs * v > MAX_PULSE_RATE_HZ) {
          Serial.printf("  rejected: %u FS/s x 1/%ld would exceed %lu Hz\n",
                        cfg.warmupSpeedFs, v, (unsigned long)MAX_PULSE_RATE_HZ);
        } else if (posU256 % U256 != 0) {
          Serial.println("  refused: not on a full step (pause at an endpoint)");
        } else if (applyMres((uint16_t)v, true)) {
          Serial.printf("  MRES = 1/%u\n", cfg.mres);
        }
        break;
      }
      case 'a': {
        const double mm = readLine("\n  peak-to-peak mm (0 = full travel): ").toDouble();
        const int64_t full = limPosU256 - limNegU256;
        const int64_t want = (int64_t)llround(mm * (double)U256_PER_MM);
        if (mm <= 0 || want > full) {
          cfg.warmupAmplitudeU256 = full;
          if (mm > 0) Serial.println("  clamped to full travel");
        } else {
          cfg.warmupAmplitudeU256 = want;
        }
        Serial.printf("  amplitude = %.3f mm peak-to-peak\n",
                      (double)cfg.warmupAmplitudeU256 / (double)U256_PER_MM);
        break;
      }
      default:
        printWarmupHelp();
        break;
    }
  }
  hotkeys = false;
  running = false;
  Serial.printf("\n# WARMUP_END,cycles=%lu,pos=%lld\n", (unsigned long)cycles,
                (long long)posU256);
}

// ----------------------------------------------------------- plan executor

bool executeRow(const PlanRow &r, uint16_t index) {
  currentRow = index;
  const char *kindName = r.kind == PLAN_MOVE ? "move"
                       : (r.kind == PLAN_STAIR ? "stair" : "set_mres");
  Serial.printf("# MOVE,%u,%s,%s,target=%ld,mres=%u\n", index, kindName,
                r.tag, (long)r.target_u256, r.mres);

  if (r.kind == PLAN_SET_MRES) {
    if (!applyMres(r.mres, true)) return false;
    Serial.printf("# MRES,%u,mscnt=%lu\n", cfg.mres,
                  (unsigned long)readMscnt());
    if (!dwellMs(r.dwell_ms)) return false;
  } else {
    if (r.mres != cfg.mres && !applyMres(r.mres, true)) return false;
    if (r.kind == PLAN_STAIR) {
      if (!stairToU256(r.target_u256, r.n_inc, r.inc_usteps, r.inc_dwell_ms))
        return false;
    } else if (!moveToU256(r.target_u256, r.speed_fs)) {
      return false;
    }
    if (!dwellMs(r.dwell_ms)) return false;
  }
  Serial.printf("# ARRIVED,%u,pos=%lld,enc=%ld,mscnt=%lu\n", index,
                (long long)posU256, (long)encoderCount(),
                (unsigned long)readMscnt());
  return true;
}

void runPlan(const PlanRow *rows, uint16_t count, const char *what) {
  if (!driverReady || !encoderReady) {
    Serial.println("# ERR,reason=not_ready");
    return;
  }
  if (count == 0) {
    Serial.println("# ERR,reason=plan_empty");
    return;
  }
  running = true;
  abortRequested = false;
  stopReason = STOP_NONE;
  digitalWrite(LED_PIN, HIGH);
  Serial.printf("# RUN_START,moves=%u,pos=%lld,what=%s\n", count,
                (long long)posU256, what);
  bool ok = true;
  for (uint16_t i = 0; i < count && !abortRequested; ++i) {
    if (!executeRow(rows[i], i)) { ok = false; break; }
  }
  digitalWrite(LED_PIN, LOW);
  currentRow = -1;
  running = false;
  if (abortRequested)
    Serial.printf("# RUN_ABORTED,pos=%lld\n", (long long)posU256);
  else if (!ok)
    Serial.printf("# RUN_FAILED,pos=%lld,reason=%u\n", (long long)posU256,
                  (unsigned)stopReason);
  else
    Serial.printf("# RUN_COMPLETE,pos=%lld\n", (long long)posU256);
}

/** Part 2: centre the stage, then run the compiled-in sequence. */
void runSequence() {
  if (!homed) {
    Serial.println("\nPart 2 needs the travel measured first. Run 1 or 4.");
    return;
  }
  Serial.printf("\nPart 2: %u rows, about %lu min. The stage moves +/-2 mm\n"
                "about the centre. Type RUN to start, anything else cancels.\n",
                MRES_PLAN_ROWS, (unsigned long)(MRES_PLAN_SECONDS / 60));
  if (!readLine("> ").equalsIgnoreCase("RUN")) {
    Serial.println("Cancelled.");
    return;
  }
  if (!goToCentre()) { Serial.println("# ERR,reason=centre_failed"); return; }
  if (posU256 % U256 != 0) {
    Serial.println("# ERR,reason=x_start_not_full_step");
    return;
  }
  runPlan(MRES_PLAN, MRES_PLAN_ROWS, "embedded");
}

// ----------------------------------------------------------- host protocol

bool handleHostCommand(const String &line) {
  if (line.equalsIgnoreCase("PING")) {
    Serial.printf("# PONG,fw=mres_test,ready=%u,driver=%u,encoder=%u,homed=%u,"
                  "moves=%u\n",
                  (driverReady && encoderReady) ? 1u : 0u,
                  driverReady ? 1u : 0u, encoderReady ? 1u : 0u,
                  homed ? 1u : 0u, (unsigned)uploadedCount);
    return true;
  }
  if (line.equalsIgnoreCase("DRIVER")) { reportDriver(); return true; }
  if (line.startsWith("HOLD ")) {
    float mult = 1.0F;
    int tpd = 0, hd = 0;
    if (sscanf(line.c_str() + 5, "%f %d %d", &mult, &tpd, &hd) != 3) {
      Serial.println("# ERR,reason=bad_hold_args");
      return true;
    }
    cfg.iholdMultiplier = mult;
    cfg.tpowerdown = (uint8_t)constrain(tpd, 0, 255);
    cfg.iholddelay = (uint8_t)constrain(hd, 0, 15);
    driver.TPOWERDOWN(cfg.tpowerdown);
    driver.iholddelay(cfg.iholddelay);
    applyCurrent(cfg.currentMa);
    Serial.printf("# HOLD_OK,ihold_mult=%.2f,tpowerdown=%u,iholddelay=%u,"
                  "ihold=%u\n", cfg.iholdMultiplier, cfg.tpowerdown,
                  cfg.iholddelay, driver.ihold());
    return true;
  }
  if (line.startsWith("INTPOL ")) {
    cfg.intpol = line.substring(7).toInt() != 0;
    driver.intpol(cfg.intpol);
    Serial.printf("# INTPOL_OK,intpol=%u\n", cfg.intpol ? 1u : 0u);
    return true;
  }
  if (line.equalsIgnoreCase("ZERO")) {
    if (running) { Serial.println("# ERR,reason=running"); return true; }
    posU256 = 0;
    indexCount = 0;
    if (pcnt) pcnt_unit_clear_count(pcnt);
    Serial.printf("# ZERO_OK,pos=0,enc=%ld\n", (long)encoderCount());
    return true;
  }
  if (line.equalsIgnoreCase("STATUS")) {
    Serial.printf("# STATUS,pos=%lld,enc=%ld,mres=%u,current_ma=%u,homed=%u,"
                  "lim_neg=%lld,lim_pos=%lld,moves=%u,running=%u\n",
                  (long long)posU256, (long)encoderCount(), cfg.mres,
                  cfg.currentMa, homed ? 1u : 0u, (long long)limNegU256,
                  (long long)limPosU256, (unsigned)uploadedCount,
                  running ? 1u : 0u);
    return true;
  }
  if (line.equalsIgnoreCase("PLAN CLEAR")) {
    uploadedCount = 0;
    Serial.println("# PLAN_CLEARED");
    return true;
  }
  if (line.startsWith("PLAN ADD ")) {
    if (uploadedCount >= MAX_UPLOADED_ROWS) {
      Serial.println("# ERR,reason=plan_full");
      return true;
    }
    char kind[12] = "", tag[24] = "";
    long target = 0;
    unsigned speed = 0, dwell = 0, m = 0, nInc = 0, incU = 0, incDwell = 0;
    const int n = sscanf(line.c_str() + 9, "%11s %ld %u %u %u %u %u %u %23s",
                         kind, &target, &speed, &dwell, &m, &nInc, &incU,
                         &incDwell, tag);
    if (n != 9) { Serial.println("# ERR,reason=bad_plan_fields"); return true; }
    uint8_t k = 255;
    if (!strcmp(kind, "move")) k = PLAN_MOVE;
    else if (!strcmp(kind, "stair")) k = PLAN_STAIR;
    else if (!strcmp(kind, "set_mres")) k = PLAN_SET_MRES;
    if (k == 255) { Serial.println("# ERR,reason=bad_kind"); return true; }
    if (m < 1 || m > 256 || (m & (m - 1))) {
      Serial.println("# ERR,reason=bad_mres");
      return true;
    }
    PlanRow &r = uploaded[uploadedCount];
    r.kind = k;
    r.target_u256 = (int32_t)target;
    r.speed_fs = (uint16_t)speed;
    r.dwell_ms = dwell;
    r.mres = (uint16_t)m;
    r.n_inc = (uint16_t)nInc;
    r.inc_usteps = (uint16_t)incU;
    r.inc_dwell_ms = incDwell;
    strncpy(r.tag, tag, sizeof(r.tag) - 1);
    r.tag[sizeof(r.tag) - 1] = '\0';
    ++uploadedCount;
    Serial.printf("# PLAN_ADD,index=%u,count=%u,tag=%s\n",
                  (unsigned)(uploadedCount - 1), (unsigned)uploadedCount, r.tag);
    return true;
  }
  if (line.equalsIgnoreCase("PLAN RUN")) {
    runPlan(uploaded, (uint16_t)uploadedCount, "uploaded");
    return true;
  }
  if (line.equalsIgnoreCase("HOME")) {
    if (homeBothEnds())
      Serial.printf("# HOME_OK,lim_neg=%lld,lim_pos=%lld\n",
                    (long long)limNegU256, (long long)limPosU256);
    else
      Serial.println("# ERR,reason=home_failed");
    return true;
  }
  return false;
}

// -------------------------------------------------------------------- menu

void printMenu() {
  Serial.println(F("\n================================================================"));
  Serial.println(F(" MRES test firmware - microstep-resolution selection"));
  Serial.printf(  " driver=%u encoder=%u homed=%u   1/%u   %u mA\n",
                driverReady ? 1u : 0u, encoderReady ? 1u : 0u, homed ? 1u : 0u,
                cfg.mres, cfg.currentMa);
  if (homed) {
    Serial.printf(" travel %.4f mm   soft limits %.3f .. %.3f mm   pos %.4f mm\n",
                  (double)(swPosU256 - swNegU256) / (double)U256_PER_MM,
                  (double)limNegU256 / (double)U256_PER_MM,
                  (double)limPosU256 / (double)U256_PER_MM,
                  (double)posU256 / (double)U256_PER_MM);
  } else {
    Serial.println(F(" travel not measured yet - run 1 or 4 before 2"));
  }
  Serial.println(F("================================================================"));
  Serial.println(F("  1  Home + warm-up    measure both ends, then shuttle"));
  Serial.printf (  "  2  MRES sequence     home to centre, run %u rows (~%lu min)\n",
                MRES_PLAN_ROWS, (unsigned long)(MRES_PLAN_SECONDS / 60));
  Serial.println(F("  3  Host plan mode    PLAN ADD / PLAN RUN for run_mres.py"));
  Serial.println(F("  4  Home only         measure both ends, then stop"));
  Serial.println(F("  5  Go to centre      move to x_start on a full step"));
  Serial.println(F("  6  Settings          current, speed, MRES, amplitude"));
  Serial.println(F("  7  Driver readback   the # DRIVER line"));
  Serial.println(F("  8  Switch state      read both end switches now"));
  Serial.println(F("  0  Release motor     EN high, outputs off"));
  Serial.println(F("  ?  This menu"));
  Serial.print  (F("\n> "));
}

void settingsMenu() {
  for (;;) {
    Serial.printf("\n  current   %u mA\n  speed     %u full-step/s\n"
                  "  MRES      1/%u\n  amplitude %.3f mm peak-to-peak\n"
                  "  ihold x%.2f   tpowerdown %u   iholddelay %u   intpol %u\n",
                  cfg.currentMa, cfg.warmupSpeedFs, cfg.mres,
                  (double)cfg.warmupAmplitudeU256 / (double)U256_PER_MM,
                  cfg.iholdMultiplier, cfg.tpowerdown, cfg.iholddelay,
                  cfg.intpol ? 1u : 0u);
    const String s = readLine("  c current  v speed  m MRES  a amplitude  "
                              "Enter to go back > ");
    if (s.length() == 0) return;
    switch (s[0]) {
      case 'c': {
        const long v = readLine("  run current mA (50-1200): ").toInt();
        if (v >= 50 && v <= 1200) applyCurrent((uint16_t)v);
        break;
      }
      case 'v': {
        const long v = readLine("  speed full-step/s (1-800): ").toInt();
        if (v >= 1 && v <= 800 && (uint32_t)v * cfg.mres <= MAX_PULSE_RATE_HZ)
          cfg.warmupSpeedFs = (uint16_t)v;
        else Serial.println("  rejected");
        break;
      }
      case 'm': {
        const long v = readLine("  MRES (1,2,4,...,256): ").toInt();
        if (v >= 1 && v <= 256 && (v & (v - 1)) == 0)
          applyMres((uint16_t)v, false);
        else Serial.println("  rejected");
        break;
      }
      case 'a': {
        if (!homed) { Serial.println("  home first"); break; }
        const double mm = readLine("  peak-to-peak mm (0 = full travel): ").toDouble();
        const int64_t full = limPosU256 - limNegU256;
        const int64_t want = (int64_t)llround(mm * (double)U256_PER_MM);
        cfg.warmupAmplitudeU256 = (mm <= 0 || want > full) ? full : want;
        break;
      }
      default:
        return;
    }
  }
}

void dispatchMenu(char key) {
  switch (key) {
    case '1':
      if (homeBothEnds()) runWarmup();
      else Serial.println("Homing failed - warm-up not started.");
      break;
    case '2': runSequence(); break;
    case '3':
      Serial.println(F("\nHost plan mode: PING DRIVER HOLD INTPOL ZERO STATUS\n"
                       "PLAN CLEAR / PLAN ADD ... / PLAN RUN, ABORT.\n"
                       "run_mres.py can take it from here. ? for the menu."));
      return;                       // no menu reprint, the host is talking
    case '4': if (homeBothEnds()) Serial.println("Homed."); break;
    case '5':
      if (goToCentre())
        Serial.printf("At centre, pos=%lld\n", (long long)posU256);
      break;
    case '6': settingsMenu(); break;
    case '7': reportDriver(); break;
    case '8':
      Serial.printf("# SWITCHES,neg=%u,pos=%u\n", switchNeg() ? 1u : 0u,
                    switchPos() ? 1u : 0u);
      break;
    case '0':
      digitalWrite(EN_PIN, HIGH);
      homed = false;
      Serial.println("Outputs off. Travel invalidated - home again before moving.");
      break;
    case '?': break;
    default:
      Serial.println("Unknown key.");
      break;
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

  Serial.println("# BOOT,fw=mres_test");
  Serial.printf("# PINS,step=%u,dir=%u,en=%u,uart_rx=%u,uart_tx=%u,"
                "enc_a=%d,enc_b=%d,enc_z=%d,sw_neg=%u,sw_pos=%u\n",
                STEP_PIN, DIR_PIN, EN_PIN, UART_RX_PIN, UART_TX_PIN,
                (int)ENC_A, (int)ENC_B, (int)ENC_Z, SW_NEG_PIN, SW_POS_PIN);
  Serial.printf("# SWITCHES,neg=%u,pos=%u,active_level=%s\n",
                switchNeg() ? 1u : 0u, switchPos() ? 1u : 0u,
                SWITCH_ACTIVE_LEVEL == LOW ? "low" : "high");

  encoderReady = setupEncoder();
  if (!encoderReady) Serial.println("# ERR,reason=pcnt_init_failed");
  driverReady = configureDriver();
  if (!driverReady) Serial.println("# ERR,reason=driver_config_failed");
  reportDriver();
  Serial.printf("# READY,driver=%u,encoder=%u,plan_rows=%u\n",
                driverReady ? 1u : 0u, encoderReady ? 1u : 0u, MRES_PLAN_ROWS);
  nextSampleMs = millis();
  printMenu();
}

void loop() {
  pumpSerial();
  if (abortRequested && !running) {
    abortRequested = false;
    Serial.println("# ABORT");
  }
  if (pendingLine.length()) {
    const String line = pendingLine;
    pendingLine = "";
    if (handleHostCommand(line)) return;
    if (line.length() == 1) { dispatchMenu(line[0]); return; }
    Serial.println("# ERR,reason=unknown_command");
    Serial.print("> ");
  }
  delay(2);
}

/*
 * ESP32-S3/TMC2209 v5 centered, bidirectional MRES campaign.
 *
 * STEP GPIO5, DIR GPIO6, UART ESP_RX GPIO18 / ESP_TX GPIO17, RGB LED GPIO48,
 * EN/ENN externally grounded -- same verified pin map as v4.
 *
 * This differs from the v4 campaign in three ways, all driven by what v4 and
 * the manual range-finding tool established:
 *
 *   1. R_SENSE_OHM is 0.11, the resistor actually fitted on this board. Every
 *      earlier sketch in this project used 0.03 (a different board's value),
 *      which under-delivered current by roughly 3x (~122 mA actual against a
 *      360 mA request). This sketch targets 360 mA RMS and actually delivers
 *      close to it.
 *
 *   2. There are no end-of-travel sensors. The range used here
 *      (REUSED_RANGE_MIN/MAX_FULL_STEPS = 2050 / 6450 full steps) is reused
 *      from the tmp_v4_manual_range_sweep session that immediately preceded
 *      this firmware, where the stage was confirmed sitting AT max when that
 *      was recorded. Flashing does not move the motor, so this firmware's own
 *      boot position (0, in its own frame -- the position counter always
 *      resets at boot, the physical stage does not) is taken to be that same
 *      physical spot, letting the range be reapplied without rejogging.
 *      THIS IS AN ASSUMPTION, NOT A MEASUREMENT: it silently produces a wrong
 *      range if the stage was jogged, by any firmware, between that session
 *      and this boot. See applyReusedRangeAndCenter().
 *
 *      Given that assumption, setup() applies the range immediately, prints
 *      it, gives a 5-second ABORT-able countdown, then drives to its centre
 *      automatically -- before RUN, not as part of it -- so the operator has
 *      a real window to catch a wrong assumption before the first automatic
 *      move happens. Manual J+/J-/NUDGE+/NUDGE-/SETMIN/SETMAX/CLEARRANGE
 *      (same commands as tmp_v4_manual_range_sweep) remain available at any
 *      time afterward to correct the range by hand if it turns out wrong.
 *
 *   3. Once centred, RUN re-centres (a no-op unless the stage was jogged
 *      since boot; position 0 stays the origin every block below returns to,
 *      exactly like v4's CAMPAIGN_START origin), then repeats the v4 MRES
 *      sequence (1, 4, 16, 32) with the oscillation block unchanged, but with
 *      the 25 mm one-way trajectories replaced by TRAJECTORY_FULL_STEPS
 *      (1200 full steps = 12 mm) out-and-back excursions run in BOTH
 *      directions from centre, not v4's single positive-only excursion. Both
 *      the DIRECT and INDIVIDUAL command-generation modes are kept, at the
 *      same three speeds as v4. Centring itself, and the pre-run throughput
 *      preflight, run at MRES 32 (the session's finest resolution) so the
 *      move lands on an exact microstep regardless of where SETMIN/SETMAX
 *      were captured -- those may not fall on a coarser MRES's step grid,
 *      since NUDGE is deliberately fine enough to creep up to a hard stop.
 *
 *      12 mm was chosen, not v4's 25 mm, because the actual measured range
 *      here is only 44 mm wide (22 mm on each side of centre): a 25 mm
 *      one-way excursion does not fit at all. 12 mm each way leaves 10 mm of
 *      physical clearance to either hard stop (about 45% margin) and keeps
 *      the full four-MRES, two-direction, two-mode, three-speed sequence to
 *      roughly 42 minutes, under the 55-minute recording budget v4 used.
 *      TRAJECTORY_FULL_STEPS is a constant below if a different distance is
 *      wanted -- confirm it is safely inside SETMIN/SETMAX before using it.
 *
 * Chopper mode is forced to SpreadCycle with interpolation off and verified
 * by readback, identically to v4 and to the manual tool, for the same
 * reason: the TMC2209's power-on defaults corrupt commanded microsteps badly
 * enough to invalidate this measurement.
 *
 * Commands over USB CDC:
 *   J+ [n] / J- [n]    Jog n full steps (default 5); clamped to the reused
 *                       range unless it's been cleared.
 *   NUDGE+ / NUDGE-     Single microstep, for the final approach to a stop.
 *   SETMIN / SETMAX     Overwrite the reused limit with the current position.
 *   CLEARRANGE          Forget the range entirely (not just override a side).
 *   RUN                 Preflight, re-centre (a no-op unless jogged since
 *                       boot), then run the full MRES campaign. Refused
 *                       until a range is set.
 *   PAUSE / RESUME      Pause any in-progress motion or run; resume it.
 *   ABORT               Stop now, not resumable.
 *   STATUS / HELP
 *
 * Only ABORT, PAUSE, RESUME, STATUS, and HELP are accepted while busy;
 * anything else is rejected with "# BUSY".
 */

#include <Arduino.h>
#include <TMCStepper.h>

constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.11F;  // Actual resistor on this board.
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t CURRENT_RMS_MA = 360;

constexpr uint16_t FULL_STEPS_PER_MM = 100;  // 2 mm/rev lead, 200 full steps/rev.
constexpr uint16_t JOG_MRES = 16;            // Fixed resolution during manual jogging.

constexpr uint16_t MRES_VALUES[] = {1, 4, 16, 32};
constexpr size_t MRES_COUNT = sizeof(MRES_VALUES) / sizeof(MRES_VALUES[0]);
constexpr uint32_t RATE_MILLIHZ[] = {27500, 70000, 200000};
constexpr const char *RATE_NAMES[] = {"SLOW", "MODERATE", "FAST"};
constexpr size_t RATE_COUNT = sizeof(RATE_MILLIHZ) / sizeof(RATE_MILLIHZ[0]);

// 12 mm each way from centre; see header for the margin/duration math.
constexpr uint32_t TRAJECTORY_FULL_STEPS = 1200UL;
constexpr uint32_t CENTERING_RATE_FULL_STEPS_S = 70;

// Reused from the tmp_v4_manual_range_sweep session immediately preceding
// this firmware: SETMIN/SETMAX reported min_full_steps=2050,
// max_full_steps=6450 with the stage confirmed sitting AT max. Flashing does
// not move the motor, so boot position (0, in this firmware's own frame) is
// taken to be that same physical spot -- see applyReusedRangeAndCenter().
// This is an assumption, not a measurement: it silently breaks if the stage
// was jogged between that session and this boot.
constexpr int32_t REUSED_RANGE_MIN_FULL_STEPS = 2050;
constexpr int32_t REUSED_RANGE_MAX_FULL_STEPS = 6450;
constexpr uint32_t AUTO_CENTER_COUNTDOWN_MS = 5000;

constexpr uint8_t OSCILLATION_CYCLES = 15;
constexpr uint32_t OSCILLATION_HALF_DWELL_MS = 1000;
constexpr uint32_t OSCILLATION_RATE_MILLIHZ = 250000;
constexpr uint32_t ENDPOINT_DWELL_MS = 1000;
constexpr uint32_t MARKER_RATE_MILLIHZ = 150000;
constexpr uint32_t MARKER_REVERSE_DWELL_MS = 1000;
constexpr uint32_t MARKER_SETTLE_MS = 500;
constexpr uint32_t CAMPAIGN_LEAD_MS = 2000;
constexpr uint32_t CAMPAIGN_TAIL_MS = 2000;
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;

constexpr int32_t JOG_DEFAULT_FULL_STEPS = 5;
constexpr uint32_t JOG_RATE_FULL_STEPS_S = 40;
constexpr uint32_t NUDGE_RATE_FULL_STEPS_S = 8;

constexpr uint32_t PREFLIGHT_RATE_MILLIHZ[] = {
    27500, 70000, 100000, 150000, 200000, 250000,
    300000, 400000, 500000, 750000, 1000000};
constexpr size_t PREFLIGHT_RATE_COUNT =
    sizeof(PREFLIGHT_RATE_MILLIHZ) / sizeof(PREFLIGHT_RATE_MILLIHZ[0]);
constexpr uint8_t PREFLIGHT_COMMANDS = 20;
constexpr uint32_t REQUIRED_COMMAND_RATE_MILLIHZ = 200000;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

bool driverReady = false;
bool busy = false;
volatile bool abortRequested = false;
volatile bool pauseRequested = false;

uint16_t currentMres = 0;
int64_t positionU32 = 0;  // 1 unit = 1/32 full step, absolute since boot.
bool haveMin = false, haveMax = false, rangeSet = false;
int64_t rangeMinCandidateU32 = 0, rangeMaxCandidateU32 = 0;
int64_t rangeMinU32 = 0, rangeMaxU32 = 0, rangeCenterU32 = 0;

uint32_t activeRateMilliHz = 0;
uint32_t maximumSupportedMilliHz = 0;
uint8_t runIndex = 0;
uint8_t markerIndex = 0;
String currentBlock = "SESSION";

void colour(uint8_t red, uint8_t green, uint8_t blue) {
  rgbLedWrite(LED_PIN, red, green, blue);
}

void logEvent(const char *event, const char *mode, uint32_t rateMilliHz,
              uint8_t iteration, int32_t markerAmplitude, const char *detail) {
  Serial.printf("%llu,%s,%u,%u,%s,%s,%lu,%u,%ld,%lld,%s\n",
                static_cast<unsigned long long>(esp_timer_get_time()), event,
                runIndex, currentMres, currentBlock.c_str(), mode,
                static_cast<unsigned long>(rateMilliHz), iteration,
                static_cast<long>(markerAmplitude),
                static_cast<long long>(positionU32), detail);
}

void printHeader() {
  Serial.println(
      "timestamp_us,event,run_index,mres,block,mode,rate_millihz,iteration,"
      "marker_amplitude_full_steps,position_u32,detail");
}

void printHelp() {
  Serial.println("# HELP,J+ [n] / J- [n]    jog n full steps (default 5)");
  Serial.println("# HELP,NUDGE+ / NUDGE-    jog one microstep");
  Serial.println("# HELP,SETMIN / SETMAX    overwrite the reused limit with current position");
  Serial.println("# HELP,CLEARRANGE         forget the range entirely");
  Serial.println("# HELP,RUN                preflight, re-center, run the MRES campaign");
  Serial.println("# HELP,range is reused from the last session and auto-centered at boot --");
  Serial.println("# HELP,see the REUSED_RANGE_ASSUMED/WARNING lines printed above at boot");
  Serial.println("# HELP,PAUSE / RESUME     pause/resume in-progress motion");
  Serial.println("# HELP,ABORT              stop now, not resumable");
  Serial.println("# HELP,STATUS / HELP");
}

void printStatus() {
  Serial.printf(
      "# STATUS,position_full_steps=%lld,position_mm=%.4f,range_set=%u,"
      "min_full_steps=%lld,max_full_steps=%lld,center_full_steps=%lld,"
      "trajectory_half_range_mm=%.2f,mres=%u,busy=%u,paused=%u,driver=%u,"
      "run_ready=%u\n",
      static_cast<long long>(positionU32 / 32),
      static_cast<double>(positionU32) / 32.0 / FULL_STEPS_PER_MM,
      static_cast<unsigned>(rangeSet),
      rangeSet ? static_cast<long long>(rangeMinU32 / 32) : 0,
      rangeSet ? static_cast<long long>(rangeMaxU32 / 32) : 0,
      rangeSet ? static_cast<long long>(rangeCenterU32 / 32) : 0,
      static_cast<double>(TRAJECTORY_FULL_STEPS) / FULL_STEPS_PER_MM,
      currentMres, static_cast<unsigned>(busy), static_cast<unsigned>(pauseRequested),
      static_cast<unsigned>(driverReady),
      static_cast<unsigned>(rangeSet && !busy));
}

void doJog(int direction, const String &args);
void doNudge(int direction);
void setLimit(bool isMin);
bool runCampaign();

// Called frequently (idle loop, and every ~32 pulses / each dwell tick during
// motion) so ABORT/PAUSE/RESUME/STATUS/HELP always land promptly. Anything
// else is rejected while busy rather than queued.
void pollSerial() {
  if (!Serial.available()) return;
  String raw = Serial.readStringUntil('\n');
  raw.trim();
  if (raw.length() == 0) return;
  String upper = raw;
  upper.toUpperCase();
  int sp = upper.indexOf(' ');
  String cmd = sp < 0 ? upper : upper.substring(0, sp);
  String args = sp < 0 ? String("") : raw.substring(sp + 1);
  args.trim();

  if (cmd == "ABORT") {
    abortRequested = true;
    pauseRequested = false;
    Serial.println("# ABORT_REQUESTED");
    return;
  }
  if (cmd == "PAUSE") {
    pauseRequested = true;
    Serial.println("# PAUSE_REQUESTED");
    return;
  }
  if (cmd == "RESUME") {
    pauseRequested = false;
    Serial.println("# RESUME_REQUESTED");
    return;
  }
  if (cmd == "STATUS") {
    printStatus();
    return;
  }
  if (cmd == "HELP") {
    printHelp();
    return;
  }

  if (busy) {
    Serial.printf("# BUSY,command_ignored=%s\n", cmd.c_str());
    return;
  }

  if (cmd == "J+" || cmd == "J-") {
    doJog(cmd == "J+" ? 1 : -1, args);
  } else if (cmd == "NUDGE+" || cmd == "NUDGE-") {
    doNudge(cmd == "NUDGE+" ? 1 : -1);
  } else if (cmd == "SETMIN") {
    setLimit(true);
  } else if (cmd == "SETMAX") {
    setLimit(false);
  } else if (cmd == "CLEARRANGE") {
    rangeSet = false;
    haveMin = false;
    haveMax = false;
    Serial.println("# RANGE_CLEARED");
  } else if (cmd == "RUN") {
    if (!rangeSet) {
      Serial.println("# RUN_REJECTED,reason=range_not_set");
    } else {
      runCampaign();
    }
  } else if (cmd.length() > 0) {
    Serial.printf("# UNKNOWN_COMMAND=%s\n", cmd.c_str());
  }
}

bool waitWhilePaused() {
  while (pauseRequested && !abortRequested) {
    pollSerial();
    delay(5);
  }
  return !abortRequested;
}

bool cancellableDwell(uint32_t durationMs) {
  const uint32_t start = millis();
  while (millis() - start < durationMs) {
    pollSerial();
    if (abortRequested) return false;
    if (pauseRequested && !waitWhilePaused()) return false;
    delay(2);
  }
  return true;
}

bool configureMotionRate(uint32_t fullStepRateMilliHz) {
  const uint64_t pulseRateMilliHz =
      static_cast<uint64_t>(fullStepRateMilliHz) * currentMres;
  if (pulseRateMilliHz == 0 || pulseRateMilliHz > UINT32_MAX) return false;
  const uint64_t periodUs = 1000000000ULL / pulseRateMilliHz;
  if (periodUs <= STEP_HIGH_US + 2) return false;
  activeRateMilliHz = fullStepRateMilliHz;
  return true;
}

// Updates positionU32 per pulse actually sent, not on completion, so
// position stays correct even when cut short by ABORT.
bool moveMicrosteps(int32_t signedMicrosteps) {
  if (signedMicrosteps == 0) return true;
  if (activeRateMilliHz == 0 || currentMres == 0) return false;
  const bool positive = signedMicrosteps > 0;
  const uint32_t count = positive
      ? static_cast<uint32_t>(signedMicrosteps)
      : static_cast<uint32_t>(-static_cast<int64_t>(signedMicrosteps));
  const uint64_t pulseRateMilliHz =
      static_cast<uint64_t>(activeRateMilliHz) * currentMres;
  const uint32_t periodUs = static_cast<uint32_t>(
      (1000000000ULL + pulseRateMilliHz / 2) / pulseRateMilliHz);
  const int32_t u32Step = 32 / currentMres;
  digitalWrite(DIR_PIN, positive ? HIGH : LOW);
  delayMicroseconds(DIR_SETUP_US);
  for (uint32_t pulse = 0; pulse < count; ++pulse) {
    const uint64_t startedUs = esp_timer_get_time();
    digitalWrite(STEP_PIN, HIGH);
    delayMicroseconds(STEP_HIGH_US);
    digitalWrite(STEP_PIN, LOW);
    positionU32 += positive ? u32Step : -u32Step;
    while (esp_timer_get_time() - startedUs < periodUs) {
      const uint64_t remaining = periodUs - (esp_timer_get_time() - startedUs);
      if (remaining > 1500) delayMicroseconds(500);
    }
    if ((pulse & 31U) == 31U) {
      pollSerial();
      if (pauseRequested && !waitWhilePaused()) return false;
      if (abortRequested) return false;
    }
  }
  return !abortRequested;
}

bool moveFullStepsDirect(int32_t signedFullSteps, uint32_t rateMilliHz) {
  if (!configureMotionRate(rateMilliHz)) return false;
  const int64_t pulses = static_cast<int64_t>(signedFullSteps) * currentMres;
  if (pulses > INT32_MAX || pulses < INT32_MIN) return false;
  return moveMicrosteps(static_cast<int32_t>(pulses));
}

bool commandOneFullStep(int direction) {
  const int32_t pulses = direction > 0 ? currentMres : -currentMres;
  if (!moveMicrosteps(pulses)) return false;
  return !abortRequested;
}

void doJog(int direction, const String &args) {
  long n = args.length() ? args.toInt() : JOG_DEFAULT_FULL_STEPS;
  if (n <= 0) n = JOG_DEFAULT_FULL_STEPS;
  int64_t target = positionU32 + static_cast<int64_t>(direction) * n * 32;
  if (rangeSet) {
    if (target < rangeMinU32) target = rangeMinU32;
    if (target > rangeMaxU32) target = rangeMaxU32;
  }
  const int32_t deltaFullSteps = static_cast<int32_t>((target - positionU32) / 32);
  if (deltaFullSteps == 0) {
    Serial.println("# JOG_NOOP,at_limit_or_zero");
    return;
  }
  busy = true;
  abortRequested = false;
  colour(0, 20, 0);
  const bool ok = moveFullStepsDirect(deltaFullSteps, JOG_RATE_FULL_STEPS_S * 1000);
  busy = false;
  colour(0, 18, 18);
  Serial.printf("# JOG_%s,requested_full_steps=%ld,applied_full_steps=%ld,ok=%u\n",
                direction > 0 ? "FWD" : "REV", n, static_cast<long>(deltaFullSteps),
                static_cast<unsigned>(ok));
  printStatus();
}

void doNudge(int direction) {
  const int32_t step = 32 / currentMres;
  const int64_t target = positionU32 + direction * step;
  if (rangeSet) {
    if (target < rangeMinU32) {
      Serial.println("# NUDGE_BLOCKED,at_min");
      return;
    }
    if (target > rangeMaxU32) {
      Serial.println("# NUDGE_BLOCKED,at_max");
      return;
    }
  }
  busy = true;
  abortRequested = false;
  colour(0, 20, 0);
  bool ok = false;
  // configureMotionRate expects a milliHz-scaled rate (Hz x 1000); the
  // constants below are plain Hz, so they must be scaled at every call site.
  if (configureMotionRate(NUDGE_RATE_FULL_STEPS_S * 1000)) {
    ok = moveMicrosteps(direction);
  }
  busy = false;
  colour(0, 18, 18);
  Serial.printf("# NUDGE_%s,ok=%u\n", direction > 0 ? "FWD" : "REV", static_cast<unsigned>(ok));
  printStatus();
}

void setLimit(bool isMin) {
  if (isMin) {
    rangeMinCandidateU32 = positionU32;
    haveMin = true;
  } else {
    rangeMaxCandidateU32 = positionU32;
    haveMax = true;
  }
  if (haveMin && haveMax) {
    rangeMinU32 = min(rangeMinCandidateU32, rangeMaxCandidateU32);
    rangeMaxU32 = max(rangeMinCandidateU32, rangeMaxCandidateU32);
    rangeCenterU32 = (rangeMinU32 + rangeMaxU32) / 2;
    rangeSet = true;
    const double travelMm =
        static_cast<double>(rangeMaxU32 - rangeMinU32) / 32.0 / FULL_STEPS_PER_MM;
    const double halfMm = travelMm / 2.0;
    const double marginMm = halfMm - static_cast<double>(TRAJECTORY_FULL_STEPS) / FULL_STEPS_PER_MM;
    Serial.printf(
        "# RANGE_SET,min_full_steps=%lld,max_full_steps=%lld,center_full_steps=%lld,"
        "travel_mm=%.3f,trajectory_margin_mm=%.3f\n",
        static_cast<long long>(rangeMinU32 / 32), static_cast<long long>(rangeMaxU32 / 32),
        static_cast<long long>(rangeCenterU32 / 32), travelMm, marginMm);
    if (marginMm < 0) {
      Serial.println("# WARNING,TRAJECTORY_FULL_STEPS_DOES_NOT_FIT_IN_RANGE");
    }
  } else {
    Serial.printf("# %s_CANDIDATE_SET,position_full_steps=%lld (send %s to complete range)\n",
                  isMin ? "MIN" : "MAX", static_cast<long long>(positionU32 / 32),
                  isMin ? "SETMAX" : "SETMIN");
  }
}

// See esp32_v4_mres_trajectory_campaign.ino: a bare CHOPCONF read is
// ambiguous because TMCStepper's UART read gives up after 2 retries and
// returns 0 with CRCerror set, which looks like a genuine MRES-code-0
// readback.
bool readChopconfChecked(uint32_t &value) {
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    const uint32_t candidate = driver.CHOPCONF();
    if (!driver.CRCerror && (candidate & 0x0FUL) != 0UL) {
      value = candidate;
      return true;
    }
    delay(5);
  }
  return false;
}

bool setMres(uint16_t mres) {
  uint8_t code = 8;
  for (uint16_t value = mres; value > 1; value >>= 1) --code;
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  for (uint8_t attempt = 1; attempt <= 5; ++attempt) {
    uint32_t chopconf = 0;
    if (!readChopconfChecked(chopconf)) {
      Serial.printf("# MRES_READ_RETRY,attempt=%u,stage=pre_write\n", attempt);
      delay(10);
      continue;
    }
    driver.CHOPCONF((chopconf & ~MRES_MASK) | (static_cast<uint32_t>(code) << 24));
    uint32_t verify = 0;
    if (!readChopconfChecked(verify)) {
      Serial.printf("# MRES_READ_RETRY,attempt=%u,stage=readback\n", attempt);
      delay(10);
      continue;
    }
    const uint8_t readback = (verify & MRES_MASK) >> 24;
    if (readback == code) {
      currentMres = mres;
      Serial.printf("# MRES_OK,mres=%u,code=%u,attempt=%u,chopconf=0x%08lX\n", mres, code,
                    attempt, static_cast<unsigned long>(verify));
      return true;
    }
    Serial.printf("# MRES_MISMATCH_RETRY,attempt=%u,wrote=%u,read=%u,chopconf=0x%08lX\n",
                  attempt, code, readback, static_cast<unsigned long>(verify));
    delay(10);
  }
  Serial.printf("# MRES_READBACK_FAILED,wrote=%u,attempts=5\n", code);
  return false;
}

bool configureDriver() {
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();
  driver.toff(5);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.rms_current(CURRENT_RMS_MA);
  driver.ihold(driver.irun());
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);

  const uint32_t gconf = driver.GCONF();
  const uint32_t chopconf = driver.CHOPCONF();
  const bool spreadCycleOn = ((gconf >> 2) & 1U) != 0U;
  const bool interpolationOff = ((chopconf >> 28) & 1U) == 0U;
  Serial.printf("# CHOPPER_CONFIG,gconf=0x%08lX,chopconf=0x%08lX,en_spreadcycle=%u,intpol=%u\n",
                static_cast<unsigned long>(gconf), static_cast<unsigned long>(chopconf),
                static_cast<unsigned>(spreadCycleOn), static_cast<unsigned>(!interpolationOff));
  if (!spreadCycleOn || !interpolationOff) {
    Serial.println("# CHOPPER_CONFIG_FAILED,refusing_to_run");
    return false;
  }
  const uint32_t drvStatus = driver.DRV_STATUS();
  const bool stealthActive = ((drvStatus >> 30) & 1U) != 0U;
  if (stealthActive) {
    Serial.println("# STEALTHCHOP_STILL_ACTIVE,refusing_to_run");
    return false;
  }
  bool uartOk = false;
  for (uint8_t attempt = 0; attempt < 10; ++attempt) {
    if (driver.test_connection() == 0 && driver.version() == 0x21) {
      uartOk = true;
      break;
    }
    delay(150);
  }
  if (!uartOk) {
    Serial.println("# TMC_UART_FAILED");
    return false;
  }
  Serial.println("# TMC_UART_OK,version=0x21,rx=18,tx=17");
  Serial.printf("# CURRENT,requested_rms_ma=%u,irun_code=%u,vsense=%u,actual_rms_ma=%u\n",
                CURRENT_RMS_MA, driver.irun(), driver.vsense(), driver.rms_current());
  return setMres(JOG_MRES);
}

bool checkOrigin(const char *where) {
  if (positionU32 == 0) return true;
  Serial.printf("# ORIGIN_FAILED,where=%s,position_u32=%lld\n", where,
                static_cast<long long>(positionU32));
  return false;
}

uint16_t markerAmplitude() {
  ++markerIndex;
  return 8 + 4 * markerIndex;
}

bool runMarker(const String &label) {
  currentBlock = "MARKER_" + label;
  const uint16_t amplitude = markerAmplitude();
  colour(18, 0, 24);
  logEvent("BLOCK_START", "MARKER", MARKER_RATE_MILLIHZ, 1, amplitude, "negative_then_return");
  if (!moveFullStepsDirect(-static_cast<int32_t>(amplitude), MARKER_RATE_MILLIHZ)) return false;
  if (!cancellableDwell(MARKER_REVERSE_DWELL_MS)) return false;
  if (!moveFullStepsDirect(amplitude, MARKER_RATE_MILLIHZ)) return false;
  if (!cancellableDwell(MARKER_SETTLE_MS)) return false;
  logEvent("BLOCK_END", "MARKER", MARKER_RATE_MILLIHZ, 1, amplitude, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runOscillation() {
  currentBlock = "OSCILLATION_MRES_" + String(currentMres) + "_15CYCLES_30S";
  colour(0, 24, 8);
  logEvent("BLOCK_START", "OSCILLATION", OSCILLATION_RATE_MILLIHZ, 1, 0,
           "15_cycles;one_microstep;1s_each_endpoint");
  if (!configureMotionRate(OSCILLATION_RATE_MILLIHZ)) return false;
  for (uint8_t cycle = 1; cycle <= OSCILLATION_CYCLES; ++cycle) {
    if (!moveMicrosteps(1)) return false;
    logEvent("OSC_FORWARD", "OSCILLATION", OSCILLATION_RATE_MILLIHZ, cycle, 0, "one_microstep");
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
    if (!moveMicrosteps(-1)) return false;
    logEvent("OSC_RETURN", "OSCILLATION", OSCILLATION_RATE_MILLIHZ, cycle, 0, "one_microstep");
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
  }
  logEvent("BLOCK_END", "OSCILLATION", OSCILLATION_RATE_MILLIHZ, 1, 0, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runDirectTrajectory(int sign, const char *rateName, uint32_t rateMilliHz) {
  const char *dirName = sign > 0 ? "POS" : "NEG";
  currentBlock =
      "TRAJECTORY_MRES_" + String(currentMres) + "_DIRECT_" + dirName + "_" + rateName;
  colour(0, 8, 28);
  logEvent("BLOCK_START", "DIRECT", rateMilliHz, 1, 0,
           "12mm_from_center;one_command_per_leg");
  const int32_t signedSteps = sign * static_cast<int32_t>(TRAJECTORY_FULL_STEPS);
  if (!moveFullStepsDirect(signedSteps, rateMilliHz)) return false;
  logEvent("ENDPOINT", "DIRECT", rateMilliHz, 1, 0, "outbound");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  if (!moveFullStepsDirect(-signedSteps, rateMilliHz)) return false;
  logEvent("ENDPOINT", "DIRECT", rateMilliHz, 1, 0, "origin");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  logEvent("BLOCK_END", "DIRECT", rateMilliHz, 1, 0, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runIndividualLeg(int direction, uint32_t rateMilliHz, const char *legName) {
  if (!configureMotionRate(rateMilliHz)) return false;
  const uint64_t startUs = esp_timer_get_time();
  for (uint32_t index = 0; index < TRAJECTORY_FULL_STEPS; ++index) {
    if (!commandOneFullStep(direction)) return false;
    if ((index & 31U) == 31U) {
      pollSerial();
      if (pauseRequested && !waitWhilePaused()) return false;
      if (abortRequested) return false;
    }
  }
  const uint64_t elapsedUs = esp_timer_get_time() - startUs;
  const double achieved = 1.0e6 * static_cast<double>(TRAJECTORY_FULL_STEPS) / elapsedUs;
  char detail[96];
  snprintf(detail, sizeof(detail), "%s;achieved_full_step_commands_s=%.3f", legName, achieved);
  logEvent("INDIVIDUAL_LEG_END", "INDIVIDUAL", rateMilliHz, 1, 0, detail);
  return true;
}

bool runIndividualTrajectory(int sign, const char *rateName, uint32_t rateMilliHz) {
  const char *dirName = sign > 0 ? "POS" : "NEG";
  currentBlock =
      "TRAJECTORY_MRES_" + String(currentMres) + "_INDIVIDUAL_" + dirName + "_" + rateName;
  colour(28, 7, 0);
  logEvent("BLOCK_START", "INDIVIDUAL", rateMilliHz, 1, 0,
           "12mm_from_center;one_full_step_command_at_a_time");
  if (!runIndividualLeg(sign, rateMilliHz, "outbound")) return false;
  logEvent("ENDPOINT", "INDIVIDUAL", rateMilliHz, 1, 0, "outbound");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  if (!runIndividualLeg(-sign, rateMilliHz, "return")) return false;
  logEvent("ENDPOINT", "INDIVIDUAL", rateMilliHz, 1, 0, "origin");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  logEvent("BLOCK_END", "INDIVIDUAL", rateMilliHz, 1, 0, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runThroughputPreflight() {
  currentBlock = "THROUGHPUT_PREFLIGHT_UNRECORDED";
  runIndex = 0;
  if (!setMres(32)) return false;
  const int64_t preflightStart = positionU32;
  colour(24, 14, 0);
  logEvent("PREFLIGHT_START", "INDIVIDUAL", 0, 0, 0, "alternating_one_full_step_commands");
  maximumSupportedMilliHz = 0;
  for (size_t rateIndex = 0; rateIndex < PREFLIGHT_RATE_COUNT; ++rateIndex) {
    const uint32_t requested = PREFLIGHT_RATE_MILLIHZ[rateIndex];
    if (!configureMotionRate(requested)) {
      Serial.printf("# PREFLIGHT_RATE_CONFIG_FAILED,rate_millihz=%lu\n",
                    static_cast<unsigned long>(requested));
      return false;
    }
    const uint64_t startUs = esp_timer_get_time();
    for (uint8_t command = 0; command < PREFLIGHT_COMMANDS; ++command) {
      const int direction = (command & 1U) == 0 ? 1 : -1;
      if (!commandOneFullStep(direction)) {
        Serial.printf("# PREFLIGHT_MOVE_FAILED,rate_millihz=%lu,command=%u\n",
                      static_cast<unsigned long>(requested), command);
        return false;
      }
    }
    const uint64_t elapsedUs = esp_timer_get_time() - startUs;
    const uint32_t achievedMilliHz =
        static_cast<uint32_t>(1.0e9 * static_cast<double>(PREFLIGHT_COMMANDS) / elapsedUs);
    const bool supported = achievedMilliHz >= requested * 95UL / 100UL;
    char detail[112];
    snprintf(detail, sizeof(detail), "requested_hz=%.3f;achieved_hz=%.3f;supported=%u",
             requested / 1000.0, achievedMilliHz / 1000.0, supported);
    logEvent("PREFLIGHT_RESULT", "INDIVIDUAL", requested, rateIndex + 1, 0, detail);
    if (!supported) break;
    maximumSupportedMilliHz = requested;
  }
  if (positionU32 != preflightStart) {
    Serial.printf("# PREFLIGHT_ORIGIN_FAILED,position_u32=%lld,expected=%lld\n",
                  static_cast<long long>(positionU32), static_cast<long long>(preflightStart));
    return false;
  }
  bool throughputReady = maximumSupportedMilliHz >= REQUIRED_COMMAND_RATE_MILLIHZ;
  char detail[96];
  snprintf(detail, sizeof(detail), "maximum_tested_supported_hz=%.3f;ready=%u",
           maximumSupportedMilliHz / 1000.0, throughputReady);
  logEvent("PREFLIGHT_END", "INDIVIDUAL", maximumSupportedMilliHz, 0, 0, detail);
  return throughputReady;
}

bool moveToCenter() {
  currentBlock = "MOVE_TO_CENTER";
  // Still at MRES 32 from the preflight: 1 u32 unit == 1 microstep there, so
  // the centre (which may not land on a coarser MRES's step grid, since
  // SETMIN/SETMAX were captured mid-NUDGE) is always exactly reachable.
  const int64_t deltaU32 = rangeCenterU32 - positionU32;
  if (deltaU32 < INT32_MIN || deltaU32 > INT32_MAX) return false;
  colour(0, 8, 28);
  // configureMotionRate expects a milliHz-scaled rate (Hz x 1000), matching
  // RATE_MILLIHZ elsewhere in this file -- CENTERING_RATE_FULL_STEPS_S is a
  // plain Hz constant, so it must be scaled here too. Forgetting this made
  // the centering move run 1000x slower than intended on first flash
  // (~2.24 pulses/s instead of ~2240 pulses/s at MRES 32).
  logEvent("CENTERING_START", "BASELINE", CENTERING_RATE_FULL_STEPS_S * 1000, 1, 0,
           "moving_to_range_center");
  if (!configureMotionRate(CENTERING_RATE_FULL_STEPS_S * 1000)) return false;
  if (!moveMicrosteps(static_cast<int32_t>(deltaU32))) return false;
  if (positionU32 != rangeCenterU32) {
    Serial.printf("# CENTERING_FAILED,position_u32=%lld,expected=%lld\n",
                  static_cast<long long>(positionU32),
                  static_cast<long long>(rangeCenterU32));
    return false;
  }
  positionU32 = 0;  // Redefine the origin as the range center from here on.
  logEvent("CENTERING_COMPLETE", "BASELINE", CENTERING_RATE_FULL_STEPS_S * 1000, 1, 0, "origin");
  return true;
}

// Like cancellableDwell, but prints a once-per-second countdown so an
// operator watching the terminal has a clear window to ABORT before an
// automatic move that was never physically confirmed (see
// applyReusedRangeAndCenter) actually starts.
bool countdownWithAbort(uint32_t totalMs, const char *label) {
  const uint32_t start = millis();
  uint32_t lastPrintedS = UINT32_MAX;
  while (millis() - start < totalMs) {
    pollSerial();
    if (abortRequested) return false;
    if (pauseRequested && !waitWhilePaused()) return false;
    const uint32_t remainingS = (totalMs - (millis() - start) + 999) / 1000;
    if (remainingS != lastPrintedS) {
      Serial.printf("# %s,seconds_remaining=%lu,send_ABORT_to_cancel\n", label,
                    static_cast<unsigned long>(remainingS));
      lastPrintedS = remainingS;
    }
    delay(20);
  }
  return true;
}

// Applies the range reused from the last manual-tool session (see
// REUSED_RANGE_MIN/MAX_FULL_STEPS) and, after a countdown the operator can
// ABORT, moves to its center automatically. Failure here (readback failure,
// or an operator ABORT) is not treated as fatal: it falls back to the idle
// jog/SETMIN/SETMAX/RUN command set so the range can be corrected by hand.
void applyReusedRangeAndCenter() {
  rangeMaxU32 = 0;  // Assumed boot position == the session's confirmed max.
  rangeMinU32 = static_cast<int64_t>(REUSED_RANGE_MIN_FULL_STEPS - REUSED_RANGE_MAX_FULL_STEPS) * 32;
  rangeCenterU32 = (rangeMinU32 + rangeMaxU32) / 2;
  rangeSet = true;
  haveMin = true;
  haveMax = true;
  const double travelMm =
      static_cast<double>(rangeMaxU32 - rangeMinU32) / 32.0 / FULL_STEPS_PER_MM;
  Serial.printf(
      "# REUSED_RANGE_ASSUMED,assumes_no_motion_since_last_session,"
      "min_full_steps=%lld,max_full_steps=%lld,center_full_steps=%lld,travel_mm=%.3f\n",
      static_cast<long long>(rangeMinU32 / 32), static_cast<long long>(rangeMaxU32 / 32),
      static_cast<long long>(rangeCenterU32 / 32), travelMm);
  Serial.println("# WARNING,VERIFY_STAGE_WAS_ACTUALLY_AT_MAX_BEFORE_TRUSTING_THIS");

  busy = true;
  abortRequested = false;
  colour(20, 16, 0);
  if (!countdownWithAbort(AUTO_CENTER_COUNTDOWN_MS, "AUTO_CENTER_STARTING")) {
    busy = false;
    colour(0, 18, 18);
    Serial.println("# AUTO_CENTER_ABORTED,manual_control_available");
    return;
  }
  if (!setMres(32) || !moveToCenter()) {
    busy = false;
    colour(32, 0, 0);
    Serial.println("# AUTO_CENTER_FAILED,manual_control_available");
    colour(0, 18, 18);
    return;
  }
  // moveToCenter() just redefined position 0 as the new origin; rebase the
  // recorded range onto that same origin so it stays valid for jogging and
  // for RUN's own (now a no-op) re-centering step.
  rangeMinU32 -= rangeCenterU32;
  rangeMaxU32 -= rangeCenterU32;
  rangeCenterU32 = 0;
  setMres(JOG_MRES);
  busy = false;
  colour(0, 18, 18);
  Serial.println("# AUTO_CENTER_COMPLETE");
}

bool runCampaign() {
  busy = true;
  abortRequested = false;
  pauseRequested = false;
  colour(0, 20, 0);

  if (!runThroughputPreflight()) {
    busy = false;
    colour(32, 0, 0);
    Serial.println("# CAMPAIGN_ABORTED_OR_FAILED,stage=preflight");
    return false;
  }
  if (!moveToCenter()) {
    busy = false;
    colour(32, 0, 0);
    Serial.println("# CAMPAIGN_ABORTED_OR_FAILED,stage=centering");
    return false;
  }

  markerIndex = 0;
  currentBlock = "CAMPAIGN";
  logEvent("CAMPAIGN_START", "BASELINE", 0, 1, 0, "centered;bidirectional");
  if (!cancellableDwell(CAMPAIGN_LEAD_MS)) {
    busy = false;
    return false;
  }

  for (size_t mresIndex = 0; mresIndex < MRES_COUNT; ++mresIndex) {
    runIndex = mresIndex + 1;
    if (!setMres(MRES_VALUES[mresIndex])) {
      busy = false;
      return false;
    }
    logEvent("RUN_CONFIG", "BASELINE", 0, 1, 0, "mres_configured");
    if (!runMarker("CONFIG_" + String(runIndex) + "_MRES_" + String(currentMres))) {
      busy = false;
      return false;
    }
    const String oscillation = "OSCILLATION_MRES_" + String(currentMres) + "_15CYCLES_30S";
    if (!runMarker(oscillation)) {
      busy = false;
      return false;
    }
    if (!runOscillation()) {
      busy = false;
      return false;
    }

    for (int sign = 1; sign >= -1; sign -= 2) {
      for (size_t mode = 0; mode < 2; ++mode) {
        for (size_t rateIndex = 0; rateIndex < RATE_COUNT; ++rateIndex) {
          const String modeName = mode == 0 ? "DIRECT" : "INDIVIDUAL";
          const String dirName = sign > 0 ? "POS" : "NEG";
          const String label = "TRAJECTORY_MRES_" + String(currentMres) + "_" + modeName + "_" +
                               dirName + "_" + RATE_NAMES[rateIndex];
          if (!runMarker(label)) {
            busy = false;
            return false;
          }
          const bool ok = mode == 0
              ? runDirectTrajectory(sign, RATE_NAMES[rateIndex], RATE_MILLIHZ[rateIndex])
              : runIndividualTrajectory(sign, RATE_NAMES[rateIndex], RATE_MILLIHZ[rateIndex]);
          if (!ok) {
            busy = false;
            return false;
          }
        }
      }
    }
    logEvent("RUN_COMPLETE", "BASELINE", 0, 1, 0, "origin");
    if (!checkOrigin("run_complete")) {
      busy = false;
      return false;
    }
  }

  currentBlock = "CAMPAIGN";
  if (!cancellableDwell(CAMPAIGN_TAIL_MS)) {
    busy = false;
    return false;
  }
  logEvent("CAMPAIGN_COMPLETE", "BASELINE", 0, 1, 0, "origin");
  busy = false;
  colour(0, 22, 22);
  Serial.println("# CAMPAIGN_FINISHED");
  return true;
}

void failWithCode(uint8_t code) {
  Serial.printf("# FAILURE_CODE=%u\n", code);
  while (true) {
    for (uint8_t blink = 0; blink < code; ++blink) {
      colour(32, 0, 0);
      delay(250);
      colour(0, 0, 0);
      delay(250);
    }
    delay(1500);
  }
}

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);

  colour(24, 24, 24);
  delay(1000);
  colour(0, 0, 0);

  Serial.begin(115200);
  Serial.setTimeout(20);
  const uint32_t serialStart = millis();
  while (!Serial && millis() - serialStart < 3000) delay(10);
  delay(4000);  // Let the USB-CDC capture attach before anything is printed.

  Serial.println("# ESP32_V5_CENTERED_MRES_CAMPAIGN (reuses last session's range; auto-centers; waits for RUN)");
  Serial.println("# STEP=GPIO5 DIR=GPIO6 UART_RX=GPIO18 UART_TX=GPIO17 EN_EXTERNALLY_GROUNDED");
  printHelp();
  printHeader();

  driverReady = configureDriver();
  if (!driverReady) failWithCode(1);

  applyReusedRangeAndCenter();

  colour(0, 18, 18);
  Serial.println("# READY,send_RUN_or_use_J+/J-/SETMIN/SETMAX_to_override_the_reused_range");
  printStatus();
}

void loop() {
  pollSerial();
  delay(5);
}

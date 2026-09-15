/*
 * ESP32-S3/TMC2209 v5 microstepping-only sequence, extended MRES ladder.
 *
 * STEP GPIO5, DIR GPIO6, UART ESP_RX GPIO18 / ESP_TX GPIO17, RGB LED GPIO48,
 * EN/ENN externally grounded -- same verified pin map as v4/v5.
 *
 * Identical to esp32_v5_microstepping_only except for the MRES ladder: this
 * one covers every doubling from full step to 1/64, 1, 2, 4, 8, 16, 32, 64,
 * instead of just 1, 4, 16, 32. Same rationale as before -- for each MRES, a
 * small CONFIG marker, a slightly larger OSCILLATION marker (amplitude grows
 * by 4 full steps every marker across the whole campaign), then 15 cycles of
 * one-microstep forward/return, 1 s dwell each side. No 25 mm/12 mm
 * trajectory legs, no DIRECT vs INDIVIDUAL command-mode comparison, no
 * throughput preflight (that only mattered for the INDIVIDUAL trajectory
 * mode, which does not exist here).
 *
 * There is deliberately no SETMIN/SETMAX/range/auto-centering step. The
 * largest single move this sketch ever makes is a marker of a few tens of
 * full steps (well under 1 mm, growing to 8+4*14=64 full steps = 0.64 mm at
 * the very last marker); the oscillation itself is +/-1 microstep. That is
 * safe to run from wherever the stage happens to be sitting when this boots,
 * which is why it starts immediately rather than waiting for a range to be
 * established the way esp32_v5_centered_mres_campaign does. Position 0 is
 * simply wherever the stage was at boot.
 *
 * Seven MRES values instead of four means 14 markers instead of 8 and a
 * correspondingly longer run: ~34 s/block x 7 blocks =~ 4 minutes total.
 *
 * R_SENSE_OHM is 0.11, the resistor actually fitted on this board. Chopper
 * mode is forced to SpreadCycle with interpolation off and verified by
 * readback, for the same reason as v4/v5: the TMC2209's power-on defaults
 * corrupt commanded microsteps badly enough to invalidate this measurement.
 *
 * Runs once automatically after boot. Commands over USB CDC: RUN (repeat),
 * PAUSE/RESUME, ABORT, STATUS, HELP.
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

constexpr uint16_t MRES_VALUES[] = {1, 2, 4, 8, 16, 32, 64};
constexpr size_t MRES_COUNT = sizeof(MRES_VALUES) / sizeof(MRES_VALUES[0]);

constexpr uint8_t OSCILLATION_CYCLES = 15;
constexpr uint32_t OSCILLATION_HALF_DWELL_MS = 1000;
constexpr uint32_t OSCILLATION_RATE_MILLIHZ = 250000;
constexpr uint32_t MARKER_RATE_MILLIHZ = 150000;
constexpr uint32_t MARKER_REVERSE_DWELL_MS = 1000;
constexpr uint32_t MARKER_SETTLE_MS = 500;
constexpr uint32_t CAMPAIGN_LEAD_MS = 2000;
constexpr uint32_t CAMPAIGN_TAIL_MS = 2000;
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

bool driverReady = false;
bool busy = false;
volatile bool abortRequested = false;
volatile bool pauseRequested = false;

uint16_t currentMres = 0;
int64_t positionU32 = 0;  // 1 unit = 1/32 full step, absolute since boot.
uint32_t activeRateMilliHz = 0;
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

void printStatus() {
  Serial.printf(
      "# STATUS,position_full_steps=%lld,mres=%u,busy=%u,paused=%u,driver=%u\n",
      static_cast<long long>(positionU32 / 32), currentMres,
      static_cast<unsigned>(busy), static_cast<unsigned>(pauseRequested),
      static_cast<unsigned>(driverReady));
}

bool runCampaign();

void pollSerial() {
  if (!Serial.available()) return;
  String raw = Serial.readStringUntil('\n');
  raw.trim();
  if (raw.length() == 0) return;
  String cmd = raw;
  cmd.toUpperCase();

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
    Serial.println("# HELP,RUN (repeat) / PAUSE / RESUME / ABORT / STATUS / HELP");
    return;
  }
  if (cmd == "RUN") {
    if (busy) {
      Serial.println("# RUN_REJECTED,reason=busy");
    } else {
      runCampaign();
    }
    return;
  }
  if (cmd.length() > 0) {
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
  return true;
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

bool runCampaign() {
  busy = true;
  abortRequested = false;
  pauseRequested = false;
  markerIndex = 0;
  positionU32 = 0;
  colour(0, 20, 0);
  currentBlock = "CAMPAIGN";
  logEvent("CAMPAIGN_START", "BASELINE", 0, 1, 0, "microstepping_only;no_trajectories;ladder_1to64");
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
  Serial.println("# CAMPAIGN_FINISHED; SEND_RUN_TO_REPEAT");
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

  Serial.println("# ESP32_V5_MICROSTEPPING_1TO64 (no trajectories; auto-runs at boot)");
  Serial.println("# STEP=GPIO5 DIR=GPIO6 UART_RX=GPIO18 UART_TX=GPIO17 EN_EXTERNALLY_GROUNDED");
  Serial.println("# EN_EXTERNALLY_GROUNDED; AUTO_RUN; COMMANDS=RUN,PAUSE,RESUME,ABORT,STATUS,HELP");
  printHeader();

  driverReady = configureDriver();
  if (!driverReady) failWithCode(1);

  if (!runCampaign()) {
    colour(32, 0, 0);
    Serial.println("# CAMPAIGN_ABORTED_OR_FAILED");
  }
}

void loop() {
  pollSerial();
  delay(5);
}

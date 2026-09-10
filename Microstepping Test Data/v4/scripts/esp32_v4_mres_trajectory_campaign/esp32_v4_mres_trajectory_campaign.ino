/*
 * ESP32-S3/TMC2209 v4 MRES and 25 mm trajectory campaign.
 *
 * STEP GPIO5, DIR GPIO6, UART ESP_RX GPIO18 / ESP_TX GPIO17,
 * RGB LED GPIO48, EN/ENN externally grounded.
 *
 * Baseline motion only: this sketch does not configure StealthChop,
 * SpreadCycle, StallGuard, or CoolStep. It configures the UART interface,
 * current, and MRES required by the experiment, and explicitly disables
 * MicroPlyer interpolation.
 *
 * Interpolation is disabled deliberately rather than left alone. CHOPCONF
 * bit 28 (intpol) powers up SET on the TMC2209, so leaving it unconfigured
 * means interpolation is ON: every commanded microstep is expanded into 256
 * sub-steps and smeared across the interval to the next step. Measured on
 * 2026-09-10, that turned a commanded one-microstep move at MRES=1 into a
 * 13-count creep followed by the remaining 243 counts over the next second,
 * instead of a single clean 256-count move. That defeats the point of the
 * oscillation block, so intpol is now forced off.
 *
 * The first campaign starts automatically after preflight. Commands over USB
 * CDC remain available afterward: RUN (repeat), STATUS, ABORT.
 */

#include <Arduino.h>
#include <TMCStepper.h>

constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.03F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t CURRENT_RMS_MA = 360;

constexpr uint16_t FULL_STEPS_PER_REV = 200;
constexpr uint16_t FULL_STEPS_PER_MM = 100;
constexpr uint32_t TRAJECTORY_FULL_STEPS = 25UL * FULL_STEPS_PER_MM;
constexpr uint16_t MRES_VALUES[] = {1, 4, 16, 32};
constexpr size_t MRES_COUNT = sizeof(MRES_VALUES) / sizeof(MRES_VALUES[0]);
constexpr uint32_t RATE_MILLIHZ[] = {27500, 70000, 200000};
constexpr const char *RATE_NAMES[] = {"SLOW", "MODERATE", "FAST"};
constexpr size_t RATE_COUNT = sizeof(RATE_MILLIHZ) / sizeof(RATE_MILLIHZ[0]);

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
bool throughputReady = false;
bool running = false;
bool abortRequested = false;
uint16_t currentMres = 0;
uint8_t runIndex = 0;
uint8_t markerIndex = 0;
int64_t positionU32 = 0;  // 1 unit = 1/32 full step.
uint32_t maximumSupportedMilliHz = 0;
uint32_t activeRateMilliHz = 0;
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

bool pollAbort() {
  if (!Serial.available()) return abortRequested;
  String command = Serial.readStringUntil('\n');
  command.trim();
  command.toUpperCase();
  if (command == "ABORT") {
    abortRequested = true;
    Serial.println("# ABORT_REQUESTED");
  }
  return abortRequested;
}

bool cancellableDwell(uint32_t durationMs) {
  const uint32_t start = millis();
  while (millis() - start < durationMs) {
    if (pollAbort()) return false;
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
  digitalWrite(DIR_PIN, positive ? HIGH : LOW);
  delayMicroseconds(DIR_SETUP_US);
  for (uint32_t pulse = 0; pulse < count; ++pulse) {
    const uint64_t startedUs = esp_timer_get_time();
    digitalWrite(STEP_PIN, HIGH);
    delayMicroseconds(STEP_HIGH_US);
    digitalWrite(STEP_PIN, LOW);
    while (esp_timer_get_time() - startedUs < periodUs) {
      const uint64_t remaining =
          periodUs - (esp_timer_get_time() - startedUs);
      if (remaining > 1500) delayMicroseconds(500);
    }
    if ((pulse & 31U) == 31U && pollAbort()) return false;
  }
  return !abortRequested;
}

bool moveFullStepsDirect(int32_t signedFullSteps, uint32_t rateMilliHz) {
  if (!configureMotionRate(rateMilliHz)) return false;
  const int64_t pulses = static_cast<int64_t>(signedFullSteps) * currentMres;
  if (pulses > INT32_MAX || pulses < INT32_MIN) return false;
  if (!moveMicrosteps(static_cast<int32_t>(pulses))) return false;
  positionU32 += static_cast<int64_t>(signedFullSteps) * 32;
  return true;
}

bool moveOneMicrostep(int direction) {
  if (!configureMotionRate(OSCILLATION_RATE_MILLIHZ)) return false;
  if (!moveMicrosteps(direction > 0 ? 1 : -1)) return false;
  positionU32 += direction * (32 / currentMres);
  return true;
}

bool commandOneFullStep(int direction) {
  const int32_t pulses = direction > 0 ? currentMres : -currentMres;
  // Each call emits exactly one physical full step, then returns before the
  // next command is submitted. The rate is configured once per leg.
  if (!moveMicrosteps(pulses)) return false;
  positionU32 += direction * 32;
  return !pollAbort();
}

// TMCStepper's UART read gives up after max_retries (2) and then returns 0
// with CRCerror set. A bare CHOPCONF read is therefore ambiguous: 0 could be
// the register contents or a dead transaction. Two things go wrong if that is
// not checked. A failed *readback* looks like MRES code 0 and aborts the
// campaign, which is what ended the 2026-09-10 17:09 run one block into
// MRES=4. Worse, a failed *initial* read makes the read-modify-write below
// store toff=0, silently disabling the driver while the readback still
// reports the MRES that was asked for. So validate every read: CRC clean, and
// toff non-zero, which is always true here because configureDriver sets it.
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
    driver.CHOPCONF((chopconf & ~MRES_MASK) |
                    (static_cast<uint32_t>(code) << 24));
    uint32_t verify = 0;
    if (!readChopconfChecked(verify)) {
      Serial.printf("# MRES_READ_RETRY,attempt=%u,stage=readback\n", attempt);
      delay(10);
      continue;
    }
    const uint8_t readback = (verify & MRES_MASK) >> 24;
    if (readback == code) {
      currentMres = mres;
      Serial.printf("# MRES_OK,mres=%u,code=%u,attempt=%u,chopconf=0x%08lX\n",
                    mres, code, attempt,
                    static_cast<unsigned long>(verify));
      return true;
    }
    Serial.printf("# MRES_MISMATCH_RETRY,attempt=%u,wrote=%u,read=%u,"
                  "chopconf=0x%08lX\n",
                  attempt, code, readback,
                  static_cast<unsigned long>(verify));
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
  driver.intpol(false);  // Powers up enabled; see header note.
  const uint32_t chopconf = driver.CHOPCONF();
  if (((chopconf >> 28) & 1U) != 0U) {
    Serial.printf("# INTPOL_DISABLE_FAILED,chopconf=0x%08lX\n",
                  static_cast<unsigned long>(chopconf));
    return false;
  }
  Serial.printf("# INTPOL_DISABLED,chopconf=0x%08lX\n",
                static_cast<unsigned long>(chopconf));
  for (uint8_t attempt = 0; attempt < 10; ++attempt) {
    if (driver.test_connection() == 0 && driver.version() == 0x21) {
      Serial.println("# TMC_UART_OK,version=0x21,rx=18,tx=17");
      return true;
    }
    delay(150);
  }
  Serial.println("# TMC_UART_FAILED");
  return false;
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
  logEvent("BLOCK_START", "MARKER", MARKER_RATE_MILLIHZ, 1, amplitude,
           "negative_then_return");
  if (!moveFullStepsDirect(-static_cast<int32_t>(amplitude),
                           MARKER_RATE_MILLIHZ)) return false;
  if (!cancellableDwell(MARKER_REVERSE_DWELL_MS)) return false;
  if (!moveFullStepsDirect(amplitude, MARKER_RATE_MILLIHZ)) return false;
  if (!cancellableDwell(MARKER_SETTLE_MS)) return false;
  logEvent("BLOCK_END", "MARKER", MARKER_RATE_MILLIHZ, 1, amplitude,
           "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runOscillation() {
  currentBlock = "OSCILLATION_MRES_" + String(currentMres) +
                 "_15CYCLES_30S";
  colour(0, 24, 8);
  logEvent("BLOCK_START", "OSCILLATION", OSCILLATION_RATE_MILLIHZ, 1, 0,
           "15_cycles;one_microstep;1s_each_endpoint");
  for (uint8_t cycle = 1; cycle <= OSCILLATION_CYCLES; ++cycle) {
    if (!moveOneMicrostep(1)) return false;
    logEvent("OSC_FORWARD", "OSCILLATION", OSCILLATION_RATE_MILLIHZ,
             cycle, 0, "one_microstep");
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
    if (!moveOneMicrostep(-1)) return false;
    logEvent("OSC_RETURN", "OSCILLATION", OSCILLATION_RATE_MILLIHZ,
             cycle, 0, "one_microstep");
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
  }
  logEvent("BLOCK_END", "OSCILLATION", OSCILLATION_RATE_MILLIHZ, 1, 0,
           "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runDirectTrajectory(const char *rateName, uint32_t rateMilliHz) {
  currentBlock = "TRAJECTORY_MRES_" + String(currentMres) + "_DIRECT_" +
                 rateName;
  colour(0, 8, 28);
  logEvent("BLOCK_START", "DIRECT", rateMilliHz, 1, 0,
           "25mm;one_command_per_leg");
  if (!moveFullStepsDirect(TRAJECTORY_FULL_STEPS, rateMilliHz)) return false;
  logEvent("ENDPOINT", "DIRECT", rateMilliHz, 1, 0, "plus_25mm");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  if (!moveFullStepsDirect(-static_cast<int32_t>(TRAJECTORY_FULL_STEPS),
                           rateMilliHz)) return false;
  logEvent("ENDPOINT", "DIRECT", rateMilliHz, 1, 0, "origin");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  logEvent("BLOCK_END", "DIRECT", rateMilliHz, 1, 0, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runIndividualLeg(int direction, uint32_t rateMilliHz,
                      const char *legName) {
  if (!configureMotionRate(rateMilliHz)) return false;
  const uint64_t startUs = esp_timer_get_time();
  for (uint32_t index = 0; index < TRAJECTORY_FULL_STEPS; ++index) {
    if (!commandOneFullStep(direction)) return false;
  }
  const uint64_t elapsedUs = esp_timer_get_time() - startUs;
  const double achieved =
      1.0e6 * static_cast<double>(TRAJECTORY_FULL_STEPS) / elapsedUs;
  char detail[96];
  snprintf(detail, sizeof(detail), "%s;achieved_full_step_commands_s=%.3f",
           legName, achieved);
  logEvent("INDIVIDUAL_LEG_END", "INDIVIDUAL", rateMilliHz, 1, 0, detail);
  return true;
}

bool runIndividualTrajectory(const char *rateName, uint32_t rateMilliHz) {
  currentBlock = "TRAJECTORY_MRES_" + String(currentMres) + "_INDIVIDUAL_" +
                 rateName;
  colour(28, 7, 0);
  logEvent("BLOCK_START", "INDIVIDUAL", rateMilliHz, 1, 0,
           "25mm;2500_one_full_step_commands_per_leg");
  if (!runIndividualLeg(1, rateMilliHz, "outbound")) return false;
  logEvent("ENDPOINT", "INDIVIDUAL", rateMilliHz, 1, 0, "plus_25mm");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  if (!runIndividualLeg(-1, rateMilliHz, "return")) return false;
  logEvent("ENDPOINT", "INDIVIDUAL", rateMilliHz, 1, 0, "origin");
  if (!cancellableDwell(ENDPOINT_DWELL_MS)) return false;
  logEvent("BLOCK_END", "INDIVIDUAL", rateMilliHz, 1, 0, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runThroughputPreflight() {
  currentBlock = "THROUGHPUT_PREFLIGHT_UNRECORDED";
  runIndex = 0;
  positionU32 = 0;
  if (!setMres(32)) return false;
  colour(24, 14, 0);
  logEvent("PREFLIGHT_START", "INDIVIDUAL", 0, 0, 0,
           "alternating_one_full_step_commands");
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
    const uint32_t achievedMilliHz = static_cast<uint32_t>(
        1.0e9 * static_cast<double>(PREFLIGHT_COMMANDS) / elapsedUs);
    const bool supported = achievedMilliHz >= requested * 95UL / 100UL;
    char detail[112];
    snprintf(detail, sizeof(detail),
             "requested_hz=%.3f;achieved_hz=%.3f;supported=%u",
             requested / 1000.0, achievedMilliHz / 1000.0, supported);
    logEvent("PREFLIGHT_RESULT", "INDIVIDUAL", requested, rateIndex + 1,
             0, detail);
    if (!supported) break;
    maximumSupportedMilliHz = requested;
  }
  if (!checkOrigin("throughput_preflight")) return false;
  throughputReady = maximumSupportedMilliHz >= REQUIRED_COMMAND_RATE_MILLIHZ;
  char detail[96];
  snprintf(detail, sizeof(detail), "maximum_tested_supported_hz=%.3f;ready=%u",
           maximumSupportedMilliHz / 1000.0, throughputReady);
  logEvent("PREFLIGHT_END", "INDIVIDUAL", maximumSupportedMilliHz, 0, 0,
           detail);
  return throughputReady;
}

bool runCampaign() {
  markerIndex = 0;
  positionU32 = 0;
  abortRequested = false;
  running = true;
  colour(0, 20, 0);
  currentBlock = "CAMPAIGN";
  logEvent("CAMPAIGN_START", "BASELINE", 0, 1, 0,
           "no_chopper_or_stallguard_configuration");
  if (!cancellableDwell(CAMPAIGN_LEAD_MS)) return false;

  for (size_t mresIndex = 0; mresIndex < MRES_COUNT; ++mresIndex) {
    runIndex = mresIndex + 1;
    if (!setMres(MRES_VALUES[mresIndex])) return false;
    logEvent("RUN_CONFIG", "BASELINE", 0, 1, 0, "mres_configured");
    if (!runMarker("CONFIG_" + String(runIndex) + "_MRES_" +
                   String(currentMres))) return false;
    const String oscillation = "OSCILLATION_MRES_" + String(currentMres) +
                               "_15CYCLES_30S";
    if (!runMarker(oscillation)) return false;
    if (!runOscillation()) return false;

    for (size_t mode = 0; mode < 2; ++mode) {
      for (size_t rateIndex = 0; rateIndex < RATE_COUNT; ++rateIndex) {
        const String modeName = mode == 0 ? "DIRECT" : "INDIVIDUAL";
        const String label = "TRAJECTORY_MRES_" + String(currentMres) + "_" +
                             modeName + "_" + RATE_NAMES[rateIndex];
        if (!runMarker(label)) return false;
        const bool ok = mode == 0
            ? runDirectTrajectory(RATE_NAMES[rateIndex], RATE_MILLIHZ[rateIndex])
            : runIndividualTrajectory(RATE_NAMES[rateIndex],
                                      RATE_MILLIHZ[rateIndex]);
        if (!ok) return false;
      }
    }
    logEvent("RUN_COMPLETE", "BASELINE", 0, 1, 0, "origin");
    if (!checkOrigin("run_complete")) return false;
  }

  currentBlock = "CAMPAIGN";
  if (!cancellableDwell(CAMPAIGN_TAIL_MS)) return false;
  logEvent("CAMPAIGN_COMPLETE", "BASELINE", 0, 1, 0, "origin");
  running = false;
  colour(0, 22, 22);
  return true;
}

void failWithCode(uint8_t code) {
  running = false;
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

  colour(24, 24, 24);  // One white boot blink.
  delay(1000);
  colour(0, 0, 0);

  Serial.begin(115200);
  Serial.setTimeout(20);
  const uint32_t serialStart = millis();
  while (!Serial && millis() - serialStart < 3000) delay(10);
  Serial.println("# ESP32_V4_MRES_TRAJECTORY_CAMPAIGN");
  Serial.println("# EN_EXTERNALLY_GROUNDED; AUTO_RUN; COMMANDS=RUN,STATUS,ABORT");
  printHeader();

  driverReady = configureDriver();
  if (!driverReady) failWithCode(1);
  if (!runThroughputPreflight()) failWithCode(2);

  colour(0, 20, 20);
  Serial.printf("# PREFLIGHT_PASSED,maximum_tested_supported_hz=%.3f,AUTO_RUN\n",
                maximumSupportedMilliHz / 1000.0);
  if (!runCampaign()) {
    running = false;
    colour(32, 0, 0);
    Serial.println("# CAMPAIGN_ABORTED_OR_FAILED");
  } else {
    Serial.println("# CAMPAIGN_FINISHED; SEND_RUN_TO_REPEAT");
  }
}

void loop() {
  if (!Serial.available()) {
    delay(5);
    return;
  }
  String command = Serial.readStringUntil('\n');
  command.trim();
  command.toUpperCase();
  if (command == "RUN") {
    if (!driverReady || !throughputReady || running) {
      Serial.println("# RUN_REJECTED");
      return;
    }
    const bool ok = runCampaign();
    running = false;
    if (!ok) {
      colour(32, 0, 0);
      Serial.println("# CAMPAIGN_ABORTED_OR_FAILED");
    } else {
      Serial.println("# CAMPAIGN_FINISHED; SEND_RUN_TO_REPEAT");
    }
  } else if (command == "STATUS") {
    Serial.printf(
        "# STATUS,driver=%u,throughput=%u,running=%u,aborted=%u,mres=%u,"
        "position_u32=%lld,max_supported_hz=%.3f\n",
        driverReady, throughputReady, running, abortRequested, currentMres,
        static_cast<long long>(positionU32),
        maximumSupportedMilliHz / 1000.0);
  } else if (command == "ABORT") {
    abortRequested = true;
    Serial.println("# ABORT_REQUESTED");
  } else if (command.length() > 0) {
    Serial.printf("# UNKNOWN_COMMAND=%s\n", command.c_str());
  }
}

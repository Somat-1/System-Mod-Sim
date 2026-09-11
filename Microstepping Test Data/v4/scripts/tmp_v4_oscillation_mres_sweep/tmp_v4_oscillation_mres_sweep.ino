/*
 * TEMPORARY diagnostic sketch -- not part of the v4 campaign.
 *
 * One-microstep back-and-forth oscillation swept across the full binary MRES
 * ladder: 1, 2, 4, 8, 16, 32. Each resolution runs for 20 s with a 2 s dwell
 * at each end of the cycle, so one cycle is 4 s and each block is 5 cycles.
 *
 * Differences from tmp_v4_oscillation_check:
 *   - six MRES values instead of four (adds 2 and 8);
 *   - 2 s half-dwell instead of 1 s;
 *   - fixed 20 s per resolution instead of a fixed cycle count;
 *   - oscillation only -- no step-integrity, pin-hold, or reference blocks.
 *
 * Driver configuration matches the corrected campaign firmware: SpreadCycle
 * forced, MicroPlyer interpolation disabled, both verified by readback before
 * the sweep will start. Neither may be left at its power-on default, because
 * on the TMC2209 en_spreadCycle resets to 0 (StealthChop) and intpol resets
 * to 1 (interpolation on) -- the two settings this measurement must not use.
 * StallGuard and CoolStep are left untouched.
 *
 * Every commanded microstep is checked against the driver's own MSCNT
 * register, which should advance by exactly 256/MRES: 256, 128, 64, 32, 16, 8
 * for the six resolutions in order. A mismatch means the step did not land as
 * commanded, independently of what the firmware believes it emitted.
 *
 * STEP GPIO5, DIR GPIO6, UART ESP_RX GPIO18 / ESP_TX GPIO17, RGB LED GPIO48,
 * EN/ENN externally grounded.
 *
 * Runs once automatically after boot. Commands: RUN (repeat), STATUS, REGS,
 * ABORT.
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

// Full binary ladder. positionU32 counts 1/32 full step, so 32/MRES stays an
// integer for every entry, and 256/MRES likewise for the expected MSCNT step.
constexpr uint16_t MRES_VALUES[] = {1, 2, 4, 8, 16, 32};
constexpr size_t MRES_COUNT = sizeof(MRES_VALUES) / sizeof(MRES_VALUES[0]);

constexpr uint32_t OSCILLATION_BLOCK_MS = 20000;
constexpr uint32_t OSCILLATION_HALF_DWELL_MS = 2000;
constexpr uint32_t OSCILLATION_CYCLE_MS = 2UL * OSCILLATION_HALF_DWELL_MS;
constexpr uint8_t OSCILLATION_CYCLES =
    OSCILLATION_BLOCK_MS / OSCILLATION_CYCLE_MS;  // 20000 / 4000 = 5
constexpr uint32_t OSCILLATION_RATE_MILLIHZ = 250000;

static_assert(OSCILLATION_CYCLES > 0, "block shorter than one cycle");
static_assert(OSCILLATION_BLOCK_MS % OSCILLATION_CYCLE_MS == 0,
              "block duration must be a whole number of cycles");

constexpr uint16_t MARKER_FULL_STEPS = 20;
constexpr uint32_t MARKER_RATE_MILLIHZ = 150000;
constexpr uint32_t MARKER_REVERSE_DWELL_MS = 1000;
constexpr uint32_t MARKER_SETTLE_MS = 500;
constexpr uint32_t RUN_LEAD_MS = 2000;
constexpr uint32_t RUN_TAIL_MS = 2000;
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

bool driverReady = false;
bool running = false;
bool abortRequested = false;
uint16_t currentMres = 0;
uint8_t runIndex = 0;
int64_t positionU32 = 0;  // 1 unit = 1/32 full step.
uint32_t activeRateMilliHz = 0;
String currentBlock = "SESSION";

void colour(uint8_t red, uint8_t green, uint8_t blue) {
  rgbLedWrite(LED_PIN, red, green, blue);
}

void printHeader() {
  Serial.println(
      "timestamp_us,event,run_index,mres,block,iteration,position_u32,"
      "mscnt_before,mscnt_after,mscnt_delta,mscnt_delta_expected,"
      "sg_result,drv_status_hex,stst,cs_actual,detail");
}

// MSCNT indexes a 1024-entry sine table spanning four full steps, so one
// commanded microstep advances it by 256 / MRES.
uint16_t expectedMscntDelta() {
  return currentMres == 0 ? 0 : static_cast<uint16_t>(256U / currentMres);
}

uint16_t readMscntSettled() {
  delay(3);
  return driver.MSCNT();
}

void logStep(const char *event, uint8_t iteration, uint16_t mscntBefore,
             uint16_t mscntAfter, int32_t microsteps, const char *detail) {
  int32_t delta = static_cast<int32_t>(mscntAfter) -
                  static_cast<int32_t>(mscntBefore);
  if (delta > 512) delta -= 1024;   // Wrapped down through 0.
  if (delta < -512) delta += 1024;  // Wrapped up through 1023.
  const int32_t expected = microsteps * expectedMscntDelta();
  const uint32_t drvStatus = driver.DRV_STATUS();
  Serial.printf("%llu,%s,%u,%u,%s,%u,%lld,%u,%u,%ld,%ld,%u,0x%08lX,%u,%u,%s\n",
                static_cast<unsigned long long>(esp_timer_get_time()), event,
                runIndex, currentMres, currentBlock.c_str(), iteration,
                static_cast<long long>(positionU32), mscntBefore, mscntAfter,
                static_cast<long>(delta), static_cast<long>(expected),
                driver.SG_RESULT(), static_cast<unsigned long>(drvStatus),
                driver.stst(), driver.cs_actual(), detail);
}

void logEvent(const char *event, uint8_t iteration, const char *detail) {
  Serial.printf("%llu,%s,%u,%u,%s,%u,%lld,,,,,,,,,%s\n",
                static_cast<unsigned long long>(esp_timer_get_time()), event,
                runIndex, currentMres, currentBlock.c_str(), iteration,
                static_cast<long long>(positionU32), detail);
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
      const uint64_t remaining = periodUs - (esp_timer_get_time() - startedUs);
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

bool stepOneMicrostep(int direction, const char *event, uint8_t iteration,
                      const char *detail) {
  if (!configureMotionRate(OSCILLATION_RATE_MILLIHZ)) return false;
  const uint16_t mscntBefore = driver.MSCNT();
  if (!moveMicrosteps(direction > 0 ? 1 : -1)) return false;
  positionU32 += direction * (32 / currentMres);
  const uint16_t mscntAfter = readMscntSettled();
  logStep(event, iteration, mscntBefore, mscntAfter, direction, detail);
  return true;
}

void dumpRegisters(const char *where) {
  const uint32_t gconf = driver.GCONF();
  const uint32_t chopconf = driver.CHOPCONF();
  const uint32_t drvStatus = driver.DRV_STATUS();
  Serial.printf("# REGS,%s,gconf=0x%08lX,chopconf=0x%08lX,pwmconf=0x%08lX,"
                "drv_status=0x%08lX,en_spreadcycle=%u,intpol=%u,stealth=%u,"
                "mscnt=%u,mres=%u,expected_mscnt_delta=%u\n",
                where, static_cast<unsigned long>(gconf),
                static_cast<unsigned long>(chopconf),
                static_cast<unsigned long>(driver.PWMCONF()),
                static_cast<unsigned long>(drvStatus),
                static_cast<unsigned>((gconf >> 2) & 1U),
                static_cast<unsigned>((chopconf >> 28) & 1U),
                static_cast<unsigned>((drvStatus >> 30) & 1U),
                driver.MSCNT(), currentMres, expectedMscntDelta());
}

// See the campaign firmware: TMCStepper returns 0 with CRCerror set once its
// retries are exhausted, so a bare CHOPCONF read cannot be trusted. toff is
// always non-zero here because configureDriver sets it.
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
      // The MRES write is a read-modify-write of CHOPCONF, so confirm it did
      // not disturb interpolation on its way through.
      if (((verify >> 28) & 1U) != 0U) {
        Serial.printf("# INTPOL_REENABLED_BY_MRES_WRITE,chopconf=0x%08lX\n",
                      static_cast<unsigned long>(verify));
        return false;
      }
      Serial.printf("# MRES_OK,mres=%u,code=%u,attempt=%u,chopconf=0x%08lX,"
                    "expected_mscnt_delta=%u\n",
                    mres, code, attempt,
                    static_cast<unsigned long>(verify), 256U / mres);
      return true;
    }
    Serial.printf("# MRES_MISMATCH_RETRY,attempt=%u,wrote=%u,read=%u\n",
                  attempt, code, readback);
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

  // Forced and verified, never left at default: en_spreadCycle resets to 0
  // (StealthChop) and intpol resets to 1 (MicroPlyer on). TPWMTHRS=0 blocks
  // any velocity-dependent switchover back to StealthChop mid-sweep.
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);

  for (uint8_t attempt = 0; attempt < 10; ++attempt) {
    if (driver.test_connection() == 0 && driver.version() == 0x21) {
      Serial.println("# TMC_UART_OK,version=0x21,rx=18,tx=17");
      break;
    }
    if (attempt == 9) {
      Serial.println("# TMC_UART_FAILED");
      return false;
    }
    delay(150);
  }

  const uint32_t gconf = driver.GCONF();
  const uint32_t chopconf = driver.CHOPCONF();
  const uint32_t drvStatus = driver.DRV_STATUS();
  const bool spreadCycleOn = ((gconf >> 2) & 1U) != 0U;
  const bool interpolationOn = ((chopconf >> 28) & 1U) != 0U;
  const bool stealthActive = ((drvStatus >> 30) & 1U) != 0U;
  dumpRegisters("startup");
  if (!spreadCycleOn || interpolationOn || stealthActive) {
    Serial.printf("# CHOPPER_CONFIG_FAILED,spreadcycle=%u,intpol=%u,"
                  "stealth=%u,refusing_to_run\n",
                  static_cast<unsigned>(spreadCycleOn),
                  static_cast<unsigned>(interpolationOn),
                  static_cast<unsigned>(stealthActive));
    return false;
  }
  Serial.println("# CHOPPER_CONFIG_OK,spreadcycle=1,intpol=0,stealth=0");
  return true;
}

bool checkOrigin(const char *where) {
  if (positionU32 == 0) return true;
  Serial.printf("# ORIGIN_FAILED,where=%s,position_u32=%lld\n", where,
                static_cast<long long>(positionU32));
  return false;
}

bool runMarker(const String &label) {
  const String previousBlock = currentBlock;
  currentBlock = "MARKER_" + label;
  colour(18, 0, 24);
  logEvent("BLOCK_START", 1, "negative_then_return;20_full_steps");
  if (!moveFullStepsDirect(-static_cast<int32_t>(MARKER_FULL_STEPS),
                           MARKER_RATE_MILLIHZ)) return false;
  if (!cancellableDwell(MARKER_REVERSE_DWELL_MS)) return false;
  if (!moveFullStepsDirect(MARKER_FULL_STEPS, MARKER_RATE_MILLIHZ))
    return false;
  if (!cancellableDwell(MARKER_SETTLE_MS)) return false;
  logEvent("BLOCK_END", 1, "origin");
  const bool ok = checkOrigin(currentBlock.c_str());
  currentBlock = previousBlock;
  return ok;
}

bool runOscillation() {
  currentBlock = "OSCILLATION_MRES_" + String(currentMres) + "_" +
                 String(OSCILLATION_CYCLES) + "CYCLES_20S";
  colour(0, 24, 8);
  char detail[128];
  snprintf(detail, sizeof(detail),
           "cycles=%u;half_dwell_ms=%lu;block_ms=%lu;expected_mscnt_delta=%u",
           OSCILLATION_CYCLES,
           static_cast<unsigned long>(OSCILLATION_HALF_DWELL_MS),
           static_cast<unsigned long>(OSCILLATION_BLOCK_MS),
           expectedMscntDelta());
  logEvent("BLOCK_START", 1, detail);
  const uint64_t startedUs = esp_timer_get_time();
  for (uint8_t cycle = 1; cycle <= OSCILLATION_CYCLES; ++cycle) {
    if (!stepOneMicrostep(1, "OSC_FORWARD", cycle, "one_microstep"))
      return false;
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
    if (!stepOneMicrostep(-1, "OSC_RETURN", cycle, "one_microstep"))
      return false;
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
  }
  char elapsed[64];
  snprintf(elapsed, sizeof(elapsed), "origin;elapsed_ms=%llu",
           static_cast<unsigned long long>(
               (esp_timer_get_time() - startedUs) / 1000ULL));
  logEvent("BLOCK_END", 1, elapsed);
  return checkOrigin(currentBlock.c_str());
}

bool runSweep() {
  positionU32 = 0;
  abortRequested = false;
  running = true;
  colour(0, 20, 0);
  currentBlock = "SWEEP";
  runIndex = 0;
  logEvent("SWEEP_START", 1,
           "oscillation_only;spreadcycle;intpol_off;mres_1_2_4_8_16_32");
  if (!cancellableDwell(RUN_LEAD_MS)) return false;

  for (size_t mresIndex = 0; mresIndex < MRES_COUNT; ++mresIndex) {
    runIndex = mresIndex + 1;
    if (!setMres(MRES_VALUES[mresIndex])) return false;
    currentBlock = "CONFIG_MRES_" + String(currentMres);
    logEvent("RUN_CONFIG", 1, "mres_configured");
    if (!runMarker(currentBlock)) return false;
    if (!runOscillation()) return false;
    currentBlock = "SWEEP";
    logEvent("RUN_COMPLETE", 1, "origin");
    if (!checkOrigin("run_complete")) return false;
  }

  currentBlock = "SWEEP";
  if (!cancellableDwell(RUN_TAIL_MS)) return false;
  dumpRegisters("final");
  logEvent("SWEEP_COMPLETE", 1, "origin");
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
  // The S3's native USB-CDC only flushes once a host attaches, so anything
  // printed during re-enumeration is lost. Give a capture time to open.
  delay(4000);
  Serial.println("# TMP_V4_OSCILLATION_MRES_SWEEP (temporary diagnostic)");
  Serial.println("# MRES=1,2,4,8,16,32; 20s per size; 2s half-dwell; "
                 "SpreadCycle; MicroPlyer OFF");
  Serial.println("# EN_EXTERNALLY_GROUNDED; AUTO_RUN; "
                 "COMMANDS=RUN,STATUS,REGS,ABORT");
  const esp_reset_reason_t resetReason = esp_reset_reason();
  Serial.printf("# RESET_REASON_CODE=%d\n", static_cast<int>(resetReason));
  printHeader();

  driverReady = configureDriver();
  if (!driverReady) failWithCode(1);

  if (!runSweep()) {
    running = false;
    colour(32, 0, 0);
    Serial.println("# SWEEP_ABORTED_OR_FAILED");
  } else {
    Serial.println("# SWEEP_FINISHED; SEND_RUN_TO_REPEAT");
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
    if (!driverReady || running) {
      Serial.println("# RUN_REJECTED");
      return;
    }
    const bool ok = runSweep();
    running = false;
    if (!ok) {
      colour(32, 0, 0);
      Serial.println("# SWEEP_ABORTED_OR_FAILED");
    } else {
      Serial.println("# SWEEP_FINISHED; SEND_RUN_TO_REPEAT");
    }
  } else if (command == "REGS") {
    dumpRegisters("on_demand");
  } else if (command == "STATUS") {
    Serial.printf("# STATUS,driver=%u,running=%u,aborted=%u,mres=%u,"
                  "position_u32=%lld\n",
                  driverReady, running, abortRequested, currentMres,
                  static_cast<long long>(positionU32));
  } else if (command == "ABORT") {
    abortRequested = true;
    Serial.println("# ABORT_REQUESTED");
  } else if (command.length() > 0) {
    Serial.printf("# UNKNOWN_COMMAND=%s\n", command.c_str());
  }
}

/*
 * TEMPORARY diagnostic sketch -- not part of the v4 campaign.
 *
 * Isolates the one-microstep back-and-forth oscillation block of
 * scripts/esp32_v4_mres_trajectory_campaign, at the same MRES values
 * (1, 4, 16, 32), the same 250 full-step/s move rate, and the same 1 s
 * dwells, but with 5 cycles instead of 15 and no 25 mm trajectories. The
 * whole run is about one minute.
 *
 * Pin map, driver configuration, and MRES write/readback are copied verbatim
 * from the campaign firmware so this reproduces its conditions exactly. Like
 * the campaign it deliberately leaves StealthChop, SpreadCycle, StallGuard,
 * CoolStep, and interpolation at their power-on defaults.
 *
 * Added diagnostics, which the campaign does not log:
 *   - MSCNT (driver microstep counter) before and after every single step,
 *     with the delta the step should have produced. MSCNT moving proves the
 *     STEP pulse reached the driver and the driver advanced its sine table;
 *     a rotor that does not follow a moving MSCNT is a torque/friction/
 *     coupling problem, not a pulse-generation problem.
 *   - SG_RESULT, DRV_STATUS, the stst (standstill) flag, and cs_actual.
 *   - A one-full-step out-and-back reference move after each MRES block, so
 *     a visible move of known size sits next to the microstep zig-zags.
 *
 * STEP GPIO5, DIR GPIO6, UART ESP_RX GPIO18 / ESP_TX GPIO17, RGB LED GPIO48,
 * EN/ENN externally grounded.
 *
 * Runs once automatically after boot. Commands: RUN (repeat), STATUS, ABORT.
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

constexpr uint16_t MRES_VALUES[] = {1, 4, 16, 32};
constexpr size_t MRES_COUNT = sizeof(MRES_VALUES) / sizeof(MRES_VALUES[0]);

// Identical to the campaign except for the cycle count.
constexpr uint8_t OSCILLATION_CYCLES = 5;
constexpr uint32_t OSCILLATION_HALF_DWELL_MS = 1000;
constexpr uint32_t OSCILLATION_RATE_MILLIHZ = 250000;

constexpr uint16_t MARKER_FULL_STEPS = 20;
constexpr uint32_t MARKER_RATE_MILLIHZ = 150000;
constexpr uint32_t MARKER_REVERSE_DWELL_MS = 1000;
constexpr uint32_t MARKER_SETTLE_MS = 500;
constexpr uint32_t REFERENCE_DWELL_MS = 1000;
constexpr uint32_t RUN_LEAD_MS = 2000;
constexpr uint32_t RUN_TAIL_MS = 2000;
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;

// The TMC2209 powers up with MicroPlyer interpolation ENABLED (CHOPCONF bit
// 28). Leaving it unconfigured, as the campaign firmware does, therefore does
// not mean "off": each commanded microstep is interpolated into 256 sub-steps
// and spread over the interval to the next step. Set this true to turn it off
// and get one crisp move per commanded microstep.
constexpr bool DISABLE_INTERPOLATION = true;

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

// A second MSCNT read issued straight after the first returns a partially
// settled value, so let the UART line and the sequencer settle before
// sampling. Without this the post-step reading understates the move badly
// (13 counts instead of 256 at MRES=1).
uint16_t readMscntSettled() {
  delay(3);
  return driver.MSCNT();
}

// MSCNT indexes a 1024-entry sine table spanning four full steps, so one
// commanded microstep advances it by 256 / MRES.
uint16_t expectedMscntDelta() {
  return currentMres == 0 ? 0 : static_cast<uint16_t>(256U / currentMres);
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

// One microstep, with MSCNT sampled either side of the pulse.
bool stepOneMicrostep(int direction, const char *event, uint8_t iteration,
                      const char *detail) {
  if (!configureMotionRate(OSCILLATION_RATE_MILLIHZ)) return false;
  const uint16_t mscntBefore = driver.MSCNT();
  if (!moveMicrosteps(direction > 0 ? 1 : -1)) return false;
  positionU32 += direction * (32 / currentMres);
  const uint16_t mscntAfter = driver.MSCNT();
  logStep(event, iteration, mscntBefore, mscntAfter, direction, detail);
  return true;
}

bool setMres(uint16_t mres) {
  uint8_t code = 8;
  for (uint16_t value = mres; value > 1; value >>= 1) --code;
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  uint32_t chopconf = driver.CHOPCONF();
  driver.CHOPCONF((chopconf & ~MRES_MASK) |
                  (static_cast<uint32_t>(code) << 24));
  const uint8_t readback = (driver.CHOPCONF() & MRES_MASK) >> 24;
  if (readback != code) {
    Serial.printf("# MRES_READBACK_FAILED,wrote=%u,read=%u\n", code, readback);
    return false;
  }
  currentMres = mres;
  Serial.printf("# MRES_OK,mres=%u,code=%u,expected_mscnt_delta=%u\n", mres,
                code, 256U / mres);
  return true;
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
  if (DISABLE_INTERPOLATION) {
    driver.intpol(false);
    Serial.printf("# INTPOL_DISABLED,chopconf=0x%08lX,intpol_bit=%u\n",
                  static_cast<unsigned long>(driver.CHOPCONF()),
                  static_cast<unsigned>((driver.CHOPCONF() >> 28) & 1U));
  } else {
    Serial.printf("# INTPOL_LEFT_AT_DEFAULT,chopconf=0x%08lX,intpol_bit=%u\n",
                  static_cast<unsigned long>(driver.CHOPCONF()),
                  static_cast<unsigned>((driver.CHOPCONF() >> 28) & 1U));
  }
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

// IOIN reports the live logic level on the driver's own input pins. Reading it
// with STEP driven high and then low says whether the pulse is arriving at the
// chip at all, which separates a wiring fault from a firmware fault.
void dumpRegisters(const char *where) {
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);
  delayMicroseconds(200);
  const uint32_t ioinLow = driver.IOIN();
  // STEP is held high for the whole UART transaction, not just a pulse, so
  // the driver cannot miss it if the line is connected.
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(200);
  const uint32_t ioinStepHigh = driver.IOIN();
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, HIGH);
  delayMicroseconds(200);
  const uint32_t ioinDirHigh = driver.IOIN();
  digitalWrite(DIR_PIN, LOW);
  delayMicroseconds(200);
  // digitalRead() on an output pin returns the actual pad level on the ESP32,
  // so a pin that cannot reach a high level (shorted, or driven by something
  // else) is distinguishable here from an intact pad behind a broken wire.
  digitalWrite(STEP_PIN, HIGH);
  digitalWrite(DIR_PIN, HIGH);
  delayMicroseconds(200);
  const int stepPadHigh = digitalRead(STEP_PIN);
  const int dirPadHigh = digitalRead(DIR_PIN);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);
  delayMicroseconds(200);
  const int stepPadLow = digitalRead(STEP_PIN);
  const int dirPadLow = digitalRead(DIR_PIN);
  Serial.printf("# PADTEST,%s,step_pad_low=%d,step_pad_high=%d,"
                "dir_pad_low=%d,dir_pad_high=%d,esp_step_pad_ok=%u,"
                "esp_dir_pad_ok=%u\n",
                where, stepPadLow, stepPadHigh, dirPadLow, dirPadHigh,
                static_cast<unsigned>(stepPadHigh == 1 && stepPadLow == 0),
                static_cast<unsigned>(dirPadHigh == 1 && dirPadLow == 0));
  Serial.printf("# PINTEST,%s,ioin_low=0x%08lX,ioin_step_high=0x%08lX,"
                "ioin_dir_high=0x%08lX,step_bit_low=%u,step_bit_high=%u,"
                "dir_bit_low=%u,dir_bit_high=%u,step_line_ok=%u,"
                "dir_line_ok=%u\n",
                where, static_cast<unsigned long>(ioinLow),
                static_cast<unsigned long>(ioinStepHigh),
                static_cast<unsigned long>(ioinDirHigh),
                static_cast<unsigned>((ioinLow >> 7) & 1U),
                static_cast<unsigned>((ioinStepHigh >> 7) & 1U),
                static_cast<unsigned>((ioinLow >> 9) & 1U),
                static_cast<unsigned>((ioinDirHigh >> 9) & 1U),
                static_cast<unsigned>(((ioinStepHigh >> 7) & 1U) &&
                                      !((ioinLow >> 7) & 1U)),
                static_cast<unsigned>(((ioinDirHigh >> 9) & 1U) &&
                                      !((ioinLow >> 9) & 1U)));
  Serial.printf("# REGS,%s,mscnt=%u,tstep=%lu,gconf=0x%08lX,chopconf=0x%08lX,"
                "drv_status=0x%08lX,ioin_step_low=0x%08lX,"
                "ioin_step_high=0x%08lX,ioin_enn=%u,ioin_step_bit_low=%u,"
                "ioin_step_bit_high=%u,ioin_dir_bit=%u,ioin_ms1=%u,ioin_ms2=%u,"
                "ioin_pdn=%u,ioin_sel_a=%u\n",
                where, driver.MSCNT(),
                static_cast<unsigned long>(driver.TSTEP()),
                static_cast<unsigned long>(driver.GCONF()),
                static_cast<unsigned long>(driver.CHOPCONF()),
                static_cast<unsigned long>(driver.DRV_STATUS()),
                static_cast<unsigned long>(ioinLow),
                static_cast<unsigned long>(ioinStepHigh),
                static_cast<unsigned>((ioinLow >> 0) & 1U),
                static_cast<unsigned>((ioinLow >> 7) & 1U),
                static_cast<unsigned>((ioinStepHigh >> 7) & 1U),
                static_cast<unsigned>((ioinLow >> 9) & 1U),
                static_cast<unsigned>((ioinLow >> 2) & 1U),
                static_cast<unsigned>((ioinLow >> 3) & 1U),
                static_cast<unsigned>((ioinLow >> 6) & 1U),
                static_cast<unsigned>((ioinLow >> 8) & 1U));
}

// Holds STEP and DIR in static states long enough to read with a multimeter,
// announcing each phase over serial and on the LED so the probe can be placed
// without watching the terminal. Each phase also reports what the driver sees
// on its own input pins, so an intermittent wire shows up as a phase where the
// ESP pad is high but ioin disagrees.
constexpr uint32_t HOLD_PHASE_MS = 8000;

struct HoldPhase {
  const char *name;
  uint8_t stepLevel;
  uint8_t dirLevel;
  uint8_t red, green, blue;
};

constexpr HoldPhase HOLD_PHASES[] = {
    {"STEP_LOW_DIR_LOW", LOW, LOW, 0, 0, 30},    // blue
    {"STEP_HIGH_DIR_LOW", HIGH, LOW, 30, 0, 0},  // red
    {"STEP_LOW_DIR_HIGH", LOW, HIGH, 0, 30, 0},  // green
    {"STEP_HIGH_DIR_HIGH", HIGH, HIGH, 26, 22, 0},  // yellow
};
constexpr size_t HOLD_PHASE_COUNT =
    sizeof(HOLD_PHASES) / sizeof(HOLD_PHASES[0]);

bool runPinHold(uint8_t cycles) {
  Serial.printf("# PINHOLD_START,cycles=%u,phase_ms=%lu,step_pin=%u,"
                "dir_pin=%u\n",
                cycles, static_cast<unsigned long>(HOLD_PHASE_MS), STEP_PIN,
                DIR_PIN);
  Serial.println("# PINHOLD_LED,blue=STEP_LOW/DIR_LOW,red=STEP_HIGH/DIR_LOW,"
                 "green=STEP_LOW/DIR_HIGH,yellow=STEP_HIGH/DIR_HIGH");
  for (uint8_t cycle = 1; cycle <= cycles; ++cycle) {
    for (size_t index = 0; index < HOLD_PHASE_COUNT; ++index) {
      const HoldPhase &phase = HOLD_PHASES[index];
      digitalWrite(STEP_PIN, phase.stepLevel);
      digitalWrite(DIR_PIN, phase.dirLevel);
      colour(phase.red, phase.green, phase.blue);
      delayMicroseconds(500);
      const uint32_t ioin = driver.IOIN();
      Serial.printf("# PINHOLD,cycle=%u,phase=%s,commanded_step=%u,"
                    "commanded_dir=%u,esp_step_pad=%d,esp_dir_pad=%d,"
                    "ioin=0x%08lX,ioin_step_bit=%u,ioin_dir_bit=%u,"
                    "ioin_enn=%u,measure_now_for_ms=%lu\n",
                    cycle, phase.name, phase.stepLevel, phase.dirLevel,
                    digitalRead(STEP_PIN), digitalRead(DIR_PIN),
                    static_cast<unsigned long>(ioin),
                    static_cast<unsigned>((ioin >> 7) & 1U),
                    static_cast<unsigned>((ioin >> 9) & 1U),
                    static_cast<unsigned>((ioin >> 0) & 1U),
                    static_cast<unsigned long>(HOLD_PHASE_MS));
      if (!cancellableDwell(HOLD_PHASE_MS)) return false;
    }
  }
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);
  Serial.println("# PINHOLD_END");
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

// One full step out and back at the same rate as the microsteps: a
// known-size reference next to the microstep zig-zags.
bool runReferenceFullStep() {
  currentBlock = "REFERENCE_FULLSTEP_MRES_" + String(currentMres);
  colour(0, 8, 28);
  logEvent("BLOCK_START", 1, "one_full_step_out_and_back");
  if (!configureMotionRate(OSCILLATION_RATE_MILLIHZ)) return false;
  const uint16_t beforeOut = driver.MSCNT();
  if (!moveMicrosteps(static_cast<int32_t>(currentMres))) return false;
  positionU32 += 32;
  logStep("REFERENCE_FORWARD", 1, beforeOut, readMscntSettled(),
          static_cast<int32_t>(currentMres), "one_full_step");
  if (!cancellableDwell(REFERENCE_DWELL_MS)) return false;
  const uint16_t beforeBack = driver.MSCNT();
  if (!moveMicrosteps(-static_cast<int32_t>(currentMres))) return false;
  positionU32 -= 32;
  logStep("REFERENCE_RETURN", 1, beforeBack, readMscntSettled(),
          -static_cast<int32_t>(currentMres), "one_full_step");
  if (!cancellableDwell(REFERENCE_DWELL_MS)) return false;
  logEvent("BLOCK_END", 1, "origin");
  return checkOrigin(currentBlock.c_str());
}

bool runOscillation() {
  currentBlock = "OSCILLATION_MRES_" + String(currentMres) + "_" +
                 String(OSCILLATION_CYCLES) + "CYCLES";
  colour(0, 24, 8);
  logEvent("BLOCK_START", 1,
           "one_microstep_forward;1s;one_microstep_back;1s");
  for (uint8_t cycle = 1; cycle <= OSCILLATION_CYCLES; ++cycle) {
    if (!stepOneMicrostep(1, "OSC_FORWARD", cycle, "one_microstep"))
      return false;
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
    if (!stepOneMicrostep(-1, "OSC_RETURN", cycle, "one_microstep"))
      return false;
    if (!cancellableDwell(OSCILLATION_HALF_DWELL_MS)) return false;
  }
  logEvent("BLOCK_END", 1, "origin");
  return checkOrigin(currentBlock.c_str());
}

// Does a STEP edge reach the driver at all? Ten microsteps in one direction
// cannot be hidden by MSCNT quantisation the way a single step might be, and
// TSTEP (the driver's own measured interval between step pulses) is an
// independent witness: it saturates at 0xFFFFF when no edges arrive.
bool runStepIntegrityTest() {
  currentBlock = "STEP_INTEGRITY";
  colour(28, 7, 0);
  if (!setMres(32)) return false;
  if (!configureMotionRate(OSCILLATION_RATE_MILLIHZ)) return false;
  Serial.printf("# STEP_INTEGRITY_START,step_pin=%u,dir_pin=%u,mres=32,"
                "expected_delta_per_microstep=8\n",
                STEP_PIN, DIR_PIN);
  dumpRegisters("baseline");

  for (int direction = 1; direction >= -1; direction -= 2) {
    for (uint8_t index = 1; index <= 10; ++index) {
      const uint16_t before = driver.MSCNT();
      if (!moveMicrosteps(direction)) return false;
      positionU32 += direction * (32 / currentMres);
      const uint16_t after = readMscntSettled();
      Serial.printf("# ONEWAY,dir=%d,n=%u,mscnt_before=%u,mscnt_after=%u,"
                    "tstep=%lu\n",
                    direction, index, before, after,
                    static_cast<unsigned long>(driver.TSTEP()));
      if (!cancellableDwell(150)) return false;
    }
  }

  // A continuous burst: TSTEP read straight afterwards should be close to the
  // 125 us commanded pulse period if the edges are landing.
  const uint16_t burstBefore = driver.MSCNT();
  if (!moveMicrosteps(64)) return false;
  positionU32 += 64 * (32 / currentMres);
  const uint16_t burstAfter = readMscntSettled();
  Serial.printf("# BURST,microsteps=64,mscnt_before=%u,mscnt_after=%u,"
                "expected_delta=512,tstep=%lu\n",
                burstBefore, burstAfter,
                static_cast<unsigned long>(driver.TSTEP()));
  if (!cancellableDwell(500)) return false;
  if (!moveMicrosteps(-64)) return false;
  positionU32 -= 64 * (32 / currentMres);
  Serial.printf("# STEP_INTEGRITY_END,mscnt=%u\n", driver.MSCNT());
  return checkOrigin("step_integrity");
}

bool runCheck() {
  positionU32 = 0;
  abortRequested = false;
  running = true;
  colour(0, 20, 0);
  currentBlock = "CHECK";
  runIndex = 0;
  logEvent("CHECK_START", 1, "oscillation_only;baseline_motion");
  if (!runStepIntegrityTest()) return false;
  if (!cancellableDwell(RUN_LEAD_MS)) return false;

  for (size_t mresIndex = 0; mresIndex < MRES_COUNT; ++mresIndex) {
    runIndex = mresIndex + 1;
    if (!setMres(MRES_VALUES[mresIndex])) return false;
    currentBlock = "CONFIG_MRES_" + String(currentMres);
    logEvent("RUN_CONFIG", 1, "mres_configured");
    if (!runMarker(currentBlock)) return false;
    if (!runOscillation()) return false;
    if (!runReferenceFullStep()) return false;
    currentBlock = "CHECK";
    logEvent("RUN_COMPLETE", 1, "origin");
    if (!checkOrigin("run_complete")) return false;
  }

  currentBlock = "CHECK";
  if (!cancellableDwell(RUN_TAIL_MS)) return false;
  dumpRegisters("final");
  logEvent("CHECK_COMPLETE", 1, "origin");
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
  // printed during re-enumeration is lost. Give the capture time to open.
  delay(4000);
  Serial.println("# TMP_V4_OSCILLATION_CHECK (temporary diagnostic)");
  Serial.println("# EN_EXTERNALLY_GROUNDED; AUTO_RUN; "
                 "COMMANDS=RUN,STATUS,REGS,PULSE,HOLD,ABORT");
  // A reset part-way through the sequence looks exactly like the firmware
  // "stopping": the run restarts and the LED blinks white again. The reset
  // reason distinguishes a brownout or panic from an ordinary power-up.
  const esp_reset_reason_t resetReason = esp_reset_reason();
  const char *resetName = "OTHER";
  switch (resetReason) {
    case ESP_RST_POWERON: resetName = "POWERON"; break;
    case ESP_RST_EXT: resetName = "EXTERNAL_PIN"; break;
    case ESP_RST_SW: resetName = "SOFTWARE"; break;
    case ESP_RST_PANIC: resetName = "PANIC_EXCEPTION"; break;
    case ESP_RST_INT_WDT: resetName = "INTERRUPT_WATCHDOG"; break;
    case ESP_RST_TASK_WDT: resetName = "TASK_WATCHDOG"; break;
    case ESP_RST_WDT: resetName = "OTHER_WATCHDOG"; break;
    case ESP_RST_BROWNOUT: resetName = "BROWNOUT"; break;
    case ESP_RST_USB: resetName = "USB_PERIPHERAL"; break;
    default: break;
  }
  Serial.printf("# RESET_REASON=%s,code=%d\n", resetName,
                static_cast<int>(resetReason));
  printHeader();

  driverReady = configureDriver();
  if (!driverReady) failWithCode(1);

  // Static-level phases first, so the pads can be measured with a multimeter
  // before any pulsing starts.
  if (!runPinHold(1)) {
    Serial.println("# PINHOLD_ABORTED");
  }

  if (!runCheck()) {
    running = false;
    colour(32, 0, 0);
    Serial.println("# CHECK_ABORTED_OR_FAILED");
  } else {
    Serial.println("# CHECK_FINISHED; SEND_RUN_TO_REPEAT");
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
    const bool ok = runCheck();
    running = false;
    if (!ok) {
      colour(32, 0, 0);
      Serial.println("# CHECK_ABORTED_OR_FAILED");
    } else {
      Serial.println("# CHECK_FINISHED; SEND_RUN_TO_REPEAT");
    }
  } else if (command == "STATUS") {
    Serial.printf("# STATUS,driver=%u,running=%u,aborted=%u,mres=%u,"
                  "position_u32=%lld\n",
                  driverReady, running, abortRequested, currentMres,
                  static_cast<long long>(positionU32));
  } else if (command == "HOLD") {
    abortRequested = false;
    runPinHold(4);
  } else if (command == "REGS") {
    dumpRegisters("on_demand");
  } else if (command == "PULSE") {
    if (currentMres == 0) setMres(32);
    if (configureMotionRate(OSCILLATION_RATE_MILLIHZ)) {
      const uint16_t before = driver.MSCNT();
      moveMicrosteps(64);
      const uint16_t after = readMscntSettled();
      moveMicrosteps(-64);
      Serial.printf("# PULSE,microsteps=64,mscnt_before=%u,mscnt_after=%u,"
                    "expected_delta=%u\n",
                    before, after, 64U * (256U / currentMres));
      dumpRegisters("after_pulse");
    }
  } else if (command == "ABORT") {
    abortRequested = true;
    Serial.println("# ABORT_REQUESTED");
  } else if (command.length() > 0) {
    Serial.printf("# UNKNOWN_COMMAND=%s\n", command.c_str());
  }
}

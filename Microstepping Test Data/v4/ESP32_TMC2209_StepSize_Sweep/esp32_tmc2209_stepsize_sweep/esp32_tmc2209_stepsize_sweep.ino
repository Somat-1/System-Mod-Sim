// ESP32-S3 + TMC2209: step-size sweep with StallGuard/DRV_STATUS logging.
//
// Purpose: for each microstep resolution MRES in {1,2,4,8,16,32,64}, step
// the motor back and forth by exactly one microstep every 2 s for
// STEP_BLOCK_DURATION_S, bracketed by a large (MARKER_FULL_STEPS-sized)
// out-and-back marker move before/after each block, so the resulting
// position trace segments trivially by eye: a big signature jump,
// followed by ~2 mm-scale... no, by a run of tiny single-microstep
// zig-zags whose PHYSICAL size shrinks as MRES rises (that shrinkage is
// exactly what this test is checking), followed by the next big marker.
//
// Driver setup, pin map, UART instantiation, and MRES-via-raw-CHOPCONF
// convention are carried over unchanged from
// ../../../ESP32S3_TMC2209_Chirp_Test/ESP32S3_TMC2209_Chirp_Test.ino and
// ../../../Microstepping Test Data/v2/scripts/run_identification_esp32_tmc2209/
// run_identification_esp32_tmc2209.ino -- both already run the driver in
// SpreadCycle-only mode (en_spreadCycle(true), intpol(false), constant
// IHOLD=IRUN, no pwm_autoscale/TCOOLTHRS/PWMCONF/CoolStep calls anywhere).
// That is reproduced verbatim here: StealthChop and CoolStep are never
// enabled, so neither can perturb microstep positioning during the sweep.
//
// StallGuard/DRV_STATUS logging is NEW -- neither existing sketch reads
// these registers. Logged once per single-microstep step and once per
// marker move via the TMCStepper library's SG_RESULT()/DRV_STATUS()
// (and per-bit accessors) methods.

#include <TMCStepper.h>
#include <esp_timer.h>

// --- Pin map (identical wiring to the existing two sketches) ----------
constexpr int STEP_PIN = 6;
constexpr int DIR_PIN = 7;
constexpr int EN_PIN = 5;
constexpr int UART_TX_PIN = 17;
constexpr int UART_RX_PIN = 18;
// MS1/MS2 tied to GND -> driver UART address 0 (see wiring note in the
// existing run_identification_esp32_tmc2209.ino). EN needs an external
// 10 kOhm pull-up to 3V3 so the driver stays disabled before setup() runs.

constexpr float R_SENSE_OHM = 0.03F;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;
constexpr uint16_t CURRENT_RMS_MA = 360;  // matches ESP32S3_TMC2209_Chirp_Test.ino

// --- Step timing (bit-banged, matching the chirp test's stepAt() style;
// RMT is unnecessary here -- steps are seconds apart, not a fast burst) --
constexpr uint32_t STEP_HIGH_US = 2;
constexpr uint32_t DIR_SETUP_US = 5;
constexpr uint32_t MARKER_PULSE_PERIOD_US = 1000;  // 1000 pulses/s during a marker move

// --- Sequence parameters ------------------------------------------------
const uint16_t MRES_SEQUENCE[] = {1, 2, 4, 8, 16, 32, 64};
constexpr size_t MRES_COUNT = sizeof(MRES_SEQUENCE) / sizeof(MRES_SEQUENCE[0]);
constexpr uint32_t STEP_INTERVAL_MS = 2000;      // gap between consecutive single steps
constexpr uint32_t STEP_BLOCK_DURATION_S = 35;    // per MRES, within the requested 30-40 s
// Marker size is fixed in FULL-STEP-equivalents (not microsteps), so its
// PHYSICAL travel is identical across every block regardless of the
// current MRES -- always a large, easily recognizable jump next to the
// (increasingly small) single-microstep zig-zags under test.
constexpr uint32_t MARKER_FULL_STEPS = 20;
constexpr uint32_t MARKER_SETTLE_MS = 300;
// Hard travel guard, in the FINEST microstep unit (1/64 full step) so it
// stays valid across every MRES in the sweep. 20 full-step marker moves
// out and back, worst case at MRES=64, is 20*64=1280 of these units; give
// ample margin beyond that.
constexpr int64_t MAX_ABS_POSITION_U64 = 20000;

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);

bool driverConfigured = false;
volatile bool aborted = false;
int64_t positionInFinestUnits = 0;  // signed accumulator, units of 1/64 full step
uint16_t currentMres = 0;

void logHeader() {
  Serial.println(
      "timestamp_us,event,mres,direction,step_kind,position_microsteps,"
      "position_full_steps_equiv,sg_result,drv_status_hex,stst,otpw,ot,"
      "s2ga,s2gb,ola,olb,stallguard_flag,cs_actual");
}

void logEvent(const char *event, int direction, const char *stepKind) {
  const int64_t positionMicrosteps =
      currentMres > 0 ? positionInFinestUnits * currentMres / 64 : 0;
  const double positionFullSteps =
      static_cast<double>(positionInFinestUnits) / 64.0;
  uint16_t sgResult = 0;
  uint32_t drvStatus = 0;
  bool stst = false, otpw = false, ot = false;
  bool s2ga = false, s2gb = false, ola = false, olb = false;
  bool stallguardFlag = false;
  uint8_t csActual = 0;
  if (driverConfigured) {
    sgResult = driver.SG_RESULT();
    drvStatus = driver.DRV_STATUS();
    stst = driver.stst();
    otpw = driver.otpw();
    ot = driver.ot();
    s2ga = driver.s2ga();
    s2gb = driver.s2gb();
    ola = driver.ola();
    olb = driver.olb();
    stallguardFlag = driver.stallguard();
    csActual = driver.cs_actual();
  }
  Serial.printf(
      "%lld,%s,%u,%d,%s,%lld,%.4f,%u,0x%08lX,%d,%d,%d,%d,%d,%d,%d,%d,%u\n",
      static_cast<long long>(esp_timer_get_time()), event, currentMres,
      direction, stepKind, static_cast<long long>(positionMicrosteps),
      positionFullSteps, sgResult, static_cast<unsigned long>(drvStatus),
      stst, otpw, ot, s2ga, s2gb, ola, olb, stallguardFlag, csActual);
}

bool pollAbort() {
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    line.trim();
    line.toUpperCase();
    if (line == "ABORT") {
      aborted = true;
      Serial.println("# ABORT received.");
    }
  }
  return aborted;
}

void dwellMs(uint32_t duration_ms) {
  const uint32_t start = millis();
  while (millis() - start < duration_ms) {
    if (pollAbort()) return;
    delay(5);
  }
}

// One microstep pulse at whatever MRES is currently configured.
// direction: +1 or -1.
bool singleStep(int direction) {
  digitalWrite(DIR_PIN, direction > 0 ? HIGH : LOW);
  delayMicroseconds(DIR_SETUP_US);
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(STEP_HIGH_US);
  digitalWrite(STEP_PIN, LOW);
  const int64_t finestUnitsPerMicrostep = 64 / currentMres;
  const int64_t candidate =
      positionInFinestUnits + direction * finestUnitsPerMicrostep;
  if (llabs(candidate) > MAX_ABS_POSITION_U64) {
    Serial.println("# Travel guard rejected step; aborting.");
    aborted = true;
    return false;
  }
  positionInFinestUnits = candidate;
  return true;
}

// A fast run of `count` microsteps in one direction, for marker moves.
bool pulseTrain(int direction, uint32_t count, const char *stepKind) {
  for (uint32_t i = 0; i < count; ++i) {
    if (!singleStep(direction)) return false;
    delayMicroseconds(MARKER_PULSE_PERIOD_US);
    if (pollAbort()) return false;
  }
  logEvent("MOVE", direction, stepKind);
  return true;
}

bool runMarker(const char *label) {
  Serial.printf("# Marker: %s\n", label);
  const uint32_t microsteps = MARKER_FULL_STEPS * currentMres;
  logEvent("MARKER_START", 0, "marker");
  if (!pulseTrain(-1, microsteps, "marker_negative")) return false;
  dwellMs(MARKER_SETTLE_MS);
  if (aborted) return false;
  if (!pulseTrain(+1, microsteps, "marker_return")) return false;
  dwellMs(MARKER_SETTLE_MS);
  logEvent("MARKER_END", 0, "marker");
  return !aborted;
}

bool checkAtOrigin(const char *where) {
  if (positionInFinestUnits != 0) {
    Serial.printf(
        "# Origin check FAILED at %s: position=%lld (finest units)\n",
        where, static_cast<long long>(positionInFinestUnits));
    aborted = true;
    return false;
  }
  return true;
}

// Raw CHOPCONF.MRES write + readback, same convention as both existing
// sketches (TMCStepper 0.7.3's microsteps() call is unreliable on this
// toolchain -- see their comments). mresCode = 8 - log2(mres).
bool setMres(uint16_t mres) {
  uint8_t mresCode = 8;
  uint16_t v = mres;
  while (v > 1) {
    v >>= 1;
    mresCode--;
  }
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  uint32_t chopconf = driver.CHOPCONF();
  chopconf = (chopconf & ~MRES_MASK) | (static_cast<uint32_t>(mresCode) << 24);
  driver.CHOPCONF(chopconf);
  const uint8_t readbackCode =
      static_cast<uint8_t>((driver.CHOPCONF() & MRES_MASK) >> 24);
  if (readbackCode != mresCode) {
    Serial.printf(
        "# CHOPCONF MRES readback mismatch: wrote %u, read %u\n",
        mresCode, readbackCode);
    return false;
  }
  currentMres = mres;
  Serial.printf("# MRES set to %u (CHOPCONF code %u)\n", mres, mresCode);
  logEvent("MRES_SET", 0, "config");
  return true;
}

bool configureDriver() {
  DriverSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();
  driver.toff(5);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  // StealthChop OFF (SpreadCycle only) -- required so microstep
  // positioning isn't perturbed by StealthChop's adaptive voltage-mode
  // chopper. pwm_autoscale()/TCOOLTHRS()/PWMCONF() (StealthChop/CoolStep)
  // are intentionally never called, matching both existing sketches.
  driver.en_spreadCycle(true);
  driver.intpol(false);
  driver.rms_current(CURRENT_RMS_MA);
  driver.ihold(driver.irun());  // constant current, no standstill reduction
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  uint8_t connection = 2;
  for (uint8_t attempt = 1; attempt <= 10; ++attempt) {
    connection = driver.test_connection();
    if (connection == 0) break;
    delay(200);
  }
  if (connection != 0) {
    Serial.println("# TMC2209 UART connection failed.");
    return false;
  }
  Serial.println("# TMC2209 UART connection OK.");
  return true;
}

bool runStepSizeBlock(uint16_t mres) {
  if (!setMres(mres)) return false;
  if (!runMarker("block_start")) return false;
  if (!checkAtOrigin("block_start")) return false;

  Serial.printf(
      "# Oscillating 1 microstep every %lu ms for %lu s at MRES=%u\n",
      static_cast<unsigned long>(STEP_INTERVAL_MS),
      static_cast<unsigned long>(STEP_BLOCK_DURATION_S), mres);
  const uint32_t blockStart = millis();
  int direction = +1;
  while ((millis() - blockStart) < STEP_BLOCK_DURATION_S * 1000UL) {
    if (!singleStep(direction)) return false;
    logEvent("STEP", direction, "oscillation");
    direction = -direction;
    dwellMs(STEP_INTERVAL_MS);
    if (aborted) return false;
  }
  // The alternating pattern returns to the origin after an even number of
  // steps; if the block duration produced an odd count, take one more
  // corrective step before the closing marker.
  if (positionInFinestUnits != 0) {
    const int correctiveDirection = positionInFinestUnits > 0 ? -1 : +1;
    if (!singleStep(correctiveDirection)) return false;
    logEvent("STEP", correctiveDirection, "corrective");
  }
  if (!checkAtOrigin("block_before_end_marker")) return false;

  if (!runMarker("block_end")) return false;
  return checkAtOrigin("block_end");
}

bool runSequence() {
  Serial.println("# Sequence start.");
  for (size_t i = 0; i < MRES_COUNT; ++i) {
    if (!runStepSizeBlock(MRES_SEQUENCE[i])) return false;
  }
  Serial.println("# Sequence complete.");
  return true;
}

void safeStop() {
  digitalWrite(EN_PIN, HIGH);  // active-LOW enable -> HIGH disables the driver
  Serial.println("# Driver disabled (safe stop).");
}

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  pinMode(EN_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);  // disabled until configuration succeeds
  digitalWrite(STEP_PIN, LOW);

  Serial.begin(CONSOLE_BAUD);
  Serial.setTimeout(20);
  const uint32_t waitStart = millis();
  while (!Serial && millis() - waitStart < 3000) {
    delay(10);
  }

  Serial.println("# ESP32-S3 + TMC2209 step-size sweep with StallGuard logging.");
  Serial.println("# Commands: CHECK, RUN, ABORT");
  logHeader();

  driverConfigured = configureDriver();
  if (driverConfigured) {
    digitalWrite(EN_PIN, LOW);  // enable the driver
    Serial.println("# Driver configured and enabled. Send RUN to start.");
  } else {
    Serial.println("# Driver configuration failed; RUN is disabled.");
  }
}

void loop() {
  if (!Serial.available()) {
    delay(5);
    return;
  }
  String line = Serial.readStringUntil('\n');
  line.trim();
  line.toUpperCase();

  if (line == "CHECK") {
    Serial.printf("# driverConfigured=%d aborted=%d position_finest_units=%lld\n",
                  driverConfigured, aborted,
                  static_cast<long long>(positionInFinestUnits));
  } else if (line == "RUN") {
    if (!driverConfigured) {
      Serial.println("# Cannot RUN: driver not configured.");
      return;
    }
    aborted = false;
    positionInFinestUnits = 0;
    const bool ok = runSequence();
    if (!ok || aborted) {
      Serial.println("# Sequence aborted or failed.");
    }
    safeStop();
  } else if (line == "ABORT") {
    aborted = true;
    Serial.println("# ABORT received.");
  } else if (line.length() > 0) {
    Serial.printf("# Unrecognized command: %s\n", line.c_str());
  }
}

/*
 * ESP32-S3 + TMC2209 v2 resonance-band chirp: cruise amplitude + smooth notch.
 *
 * UNFLASHED / UNTESTED. This is a fresh implementation, not a port of the
 * v1 sketch, and it has not been bench-validated on real hardware. Two
 * things must happen before a real run:
 *
 *   1. Confirm every constant in "Pin assignment (PLACEHOLDER)" and
 *      "Mirrors chirp_v2_schedule.py" below against the actual wiring and
 *      the current chirp_v2_schedule.py / SWEEP_CONFIG_LOG.md -- this file
 *      hardcodes copies of those values because firmware cannot import the
 *      Python module, so nothing here updates itself when the model
 *      changes.
 *   2. Bench-test the achievable STEP pulse rate (e.g. drive STEP_PIN
 *      alone, alternating direction, and count edges over a timed window)
 *      before trusting STEP_RATE_CLAMP_HZ. v1's original digitalWrite()
 *      busy-wait approach implied only ~111 kHz was safe; this file
 *      switches to direct GPIO register writes for the hot path
 *      specifically to buy headroom above that, but that fix is reasoned
 *      about, not measured. See chirp_v2_schedule.py's STEP_RATE_CLAMP_HZ
 *      comment and README.md Section 3 for the analysis this responds to.
 *
 * Command generation is a position-tracking control loop, not a fixed-step
 * alternation like v1: commandedStateAt(elapsed_s) analytically evaluates
 * the exact ideal microstep position for ANY elapsed time (closed-form
 * phase integration, no accumulated numerical drift), and the main loop
 * continuously steps the real position toward that target, one microstep
 * at a time, rate-limited by MIN_STEP_PERIOD_US. This is necessary because
 * v2's amplitude varies with frequency (cruise * notch, clamped), unlike
 * v1's fixed +-1 microstep alternation.
 *
 * Framework: Arduino-ESP32 3.x; library: TMCStepper.
 * Serial commands at 115200 baud: CHECK, RUN, ABORT.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <esp_timer.h>
#include <cmath>
#include <cstring>
#include "soc/gpio_struct.h"

namespace chirp_v2 {

// --- Pin assignment (PLACEHOLDER) -----------------------------------------
// Reused from v1's documented wiring (ESP32S3_TMC2209_Chirp_Test/README.md
// "Wiring" table) as a starting point ONLY. CONFIRM against the actual rig
// before flashing -- these are not re-verified for this file.
constexpr uint8_t STEP_PIN = 6;
constexpr uint8_t DIR_PIN = 7;
constexpr uint8_t EN_PIN = 5;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t TRIG_OUT_PIN = 1;
constexpr uint8_t TRIG_ECHO_PIN = 2;
// Direct GPIO_OUT register access (fast path for STEP) only covers pins
// 0-31 on the ESP32-S3. If STEP_PIN ever moves to 32+, stepPulse() below
// must switch to GPIO.out1_w1ts/out1_w1tc instead.
static_assert(STEP_PIN < 32, "STEP_PIN >= 32 needs the GPIO.out1_* registers");

constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.03F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;

// --- Mirrors chirp_v2_schedule.py -- keep these in sync by hand ----------

constexpr uint16_t MRES = 16;               // 1/16 step. See README Section 2.
constexpr uint8_t MRES_REGISTER_CODE = 4;   // TMC2209 MRES field code for 1/16.
constexpr uint16_t CURRENT_RMS_MA = 400;

constexpr float F_N_HZ_PLACEHOLDER = 176.7F;      // REPLACE with measured f_n.
constexpr float NOTCH_HALFWIDTH_HZ = 56.7F;
constexpr float NOTCH_DEPTH_RATIO = 0.10F;

constexpr float CRUISE_FULL_STEPS = 3.0F;
constexpr float CRUISE_MICROSTEPS = CRUISE_FULL_STEPS * MRES;  // 48

constexpr float STROKE_CLAMP_MM = 0.5F;
constexpr float FULL_STEP_MM = 0.010F;  // 2 mm lead / 200 full steps/rev.
constexpr float STROKE_CLAMP_MICROSTEPS =
    (STROKE_CLAMP_MM / FULL_STEP_MM) * MRES;  // 800; inactive at cruise=3.

// UNVERIFIED -- see the file header and chirp_v2_schedule.py's comment on
// this same constant. This is a peak MICROSTEP PULSE rate, not a full-step
// rate.
constexpr float STEP_RATE_CLAMP_HZ = 200000.0F;

constexpr float LOG_START_HZ = 1.0F;
constexpr float LOG_END_HZ = 60.0F;
constexpr double LOG_DURATION_S = 120.0;

constexpr float LINEAR_START_HZ = LOG_END_HZ;
constexpr float LINEAR_END_HZ = 1000.0F;
constexpr double LINEAR_DURATION_S = 120.0;

constexpr double SWEEP_DURATION_S = LOG_DURATION_S + LINEAR_DURATION_S;

constexpr double PRE_ROLL_S = 30.0;
constexpr double SYNC_BURST_S = 0.1;
constexpr double SYNC_GAP_S = 0.1;
constexpr int SYNC_BURST_COUNT = 3;
constexpr double SYNC_MARKER_DURATION_S =
    SYNC_BURST_COUNT * SYNC_BURST_S + (SYNC_BURST_COUNT - 1) * SYNC_GAP_S;
constexpr float SYNC_MARKER_FREQ_HZ = 750.0F;
constexpr float SYNC_MARKER_AMPLITUDE_FULL_STEPS = 2.0F;
constexpr float SYNC_MARKER_AMPLITUDE_MICROSTEPS =
    SYNC_MARKER_AMPLITUDE_FULL_STEPS * MRES;  // 32

constexpr double MID_DWELL_S = 5.0;
constexpr double TAIL_S = 30.0;

constexpr double LEVEL_DURATION_S =
    PRE_ROLL_S + SYNC_MARKER_DURATION_S + SWEEP_DURATION_S + MID_DWELL_S
    + SYNC_MARKER_DURATION_S + SWEEP_DURATION_S + TAIL_S;  // 546 s

// Defensive sanity ceiling: if a bug ever made the target exceed this,
// abort rather than drive the stage further than the design intends.
constexpr int32_t MAX_ABS_TARGET_MICROSTEPS =
    static_cast<int32_t>(CRUISE_MICROSTEPS * 1.2F);

// --- Timing floors (real, not yet bench-confirmed against hardware) ------

constexpr uint32_t STEP_HIGH_US = 2;
constexpr uint32_t DIR_SETUP_US = 2;
// Minimum spacing between STEP edges, regardless of how fast the control
// loop iterates. 3 us -> ~333 kHz physical ceiling, meant to sit above
// STEP_RATE_CLAMP_HZ (200 kHz) with margin. This does NOT by itself prove
// the loop can keep up -- only that it will never exceed this floor.
constexpr int64_t MIN_STEP_PERIOD_US = 3;

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);

volatile bool aborted = false;
bool driverConfigured = false;
int32_t currentMicrosteps = 0;
int8_t lastDirection = 0;  // -1, 0 (unknown/first step), or +1.
int64_t lastEdgeUs = 0;

void logEvent(const char *event, const char *segment, float frequencyHz) {
  Serial.printf("%lld,%s,%s,%.6f,%ld\n",
                static_cast<long long>(esp_timer_get_time()), event, segment,
                static_cast<double>(frequencyHz),
                static_cast<long>(currentMicrosteps));
}

void pollAbort() {
  while (Serial.available() > 0) {
    String command = Serial.readStringUntil('\n');
    command.trim();
    command.toUpperCase();
    if (command == "ABORT") {
      aborted = true;
      logEvent("ABORT_REQUEST", "SERIAL", 0.0F);
    }
  }
}

// --- Closed-form phase (no accumulated integration drift) ---------------

double phaseLogUpRad(double t) {
  const double k = static_cast<double>(LOG_END_HZ) / LOG_START_HZ;
  return 2.0 * PI * LOG_START_HZ * LOG_DURATION_S / log(k)
         * (pow(k, t / LOG_DURATION_S) - 1.0);
}

double phaseLinUpRad(double tPrime) {
  const double base = phaseLogUpRad(LOG_DURATION_S);
  return base + 2.0 * PI * (
      static_cast<double>(LINEAR_START_HZ) * tPrime
      + (static_cast<double>(LINEAR_END_HZ) - LINEAR_START_HZ)
        / (2.0 * LINEAR_DURATION_S) * tPrime * tPrime
  );
}

double phaseLinDownRad(double tau) {
  return 2.0 * PI * (
      static_cast<double>(LINEAR_END_HZ) * tau
      - (static_cast<double>(LINEAR_END_HZ) - LINEAR_START_HZ)
        / (2.0 * LINEAR_DURATION_S) * tau * tau
  );
}

double phaseLogDownRad(double tauPrime) {
  const double base = phaseLinDownRad(LINEAR_DURATION_S);
  const double kDown = static_cast<double>(LOG_START_HZ) / LOG_END_HZ;  // < 1
  return base + 2.0 * PI * LOG_END_HZ * LOG_DURATION_S / log(kDown)
         * (pow(kDown, tauPrime / LOG_DURATION_S) - 1.0);
}

// --- Amplitude: cruise * smooth notch, clamped twice ----------------------

float resonanceNotchTaper(float frequencyHz) {
  float x = (frequencyHz - F_N_HZ_PLACEHOLDER) / NOTCH_HALFWIDTH_HZ;
  x = constrain(x, -1.0F, 1.0F);
  const float window = 0.5F * (1.0F + cosf(PI * x));
  return 1.0F - (1.0F - NOTCH_DEPTH_RATIO) * window;
}

float amplitudeMicrosteps(float frequencyHz) {
  const float tapered = CRUISE_MICROSTEPS * resonanceNotchTaper(frequencyHz);
  const float omega = 2.0F * PI * frequencyHz;
  const float aStepRate = STEP_RATE_CLAMP_HZ / omega;
  float amplitude = tapered;
  if (STROKE_CLAMP_MICROSTEPS < amplitude) amplitude = STROKE_CLAMP_MICROSTEPS;
  if (aStepRate < amplitude) amplitude = aStepRate;
  return amplitude;
}

// --- Sync marker -----------------------------------------------------------

int32_t markerMicrostepsAt(double localS) {
  const double period = SYNC_BURST_S + SYNC_GAP_S;
  for (int burst = 0; burst < SYNC_BURST_COUNT; ++burst) {
    const double burstLo = burst * period;
    const double burstHi = burstLo + SYNC_BURST_S;
    if (localS >= burstLo && localS < burstHi) {
      const double t = localS - burstLo;
      return lround(SYNC_MARKER_AMPLITUDE_MICROSTEPS
                     * sin(2.0 * PI * SYNC_MARKER_FREQ_HZ * t));
    }
  }
  return 0;
}

// --- Full-level commanded state --------------------------------------------

struct CommandState {
  int32_t targetMicrosteps;
  float frequencyHz;  // NAN where there is no defined instantaneous tone.
  const char *segment;
};

CommandState commandedStateAt(double elapsedS) {
  double t = elapsedS;

  if (t < PRE_ROLL_S) return {0, NAN, "PRE_ROLL"};
  t -= PRE_ROLL_S;

  if (t < SYNC_MARKER_DURATION_S) {
    return {markerMicrostepsAt(t), NAN, "MARKER_1"};
  }
  t -= SYNC_MARKER_DURATION_S;

  if (t < SWEEP_DURATION_S) {
    float f;
    double phase;
    if (t < LOG_DURATION_S) {
      f = LOG_START_HZ * powf(LOG_END_HZ / LOG_START_HZ,
                               static_cast<float>(t / LOG_DURATION_S));
      phase = phaseLogUpRad(t);
    } else {
      const double t2 = t - LOG_DURATION_S;
      f = LINEAR_START_HZ + (LINEAR_END_HZ - LINEAR_START_HZ)
                             * static_cast<float>(t2 / LINEAR_DURATION_S);
      phase = phaseLinUpRad(t2);
    }
    const float amplitude = amplitudeMicrosteps(f);
    return {lround(amplitude * sin(phase)), f, "UP_SWEEP"};
  }
  t -= SWEEP_DURATION_S;

  if (t < MID_DWELL_S) return {0, NAN, "MID_DWELL"};
  t -= MID_DWELL_S;

  if (t < SYNC_MARKER_DURATION_S) {
    return {markerMicrostepsAt(t), NAN, "MARKER_2"};
  }
  t -= SYNC_MARKER_DURATION_S;

  if (t < SWEEP_DURATION_S) {
    float f;
    double phase;
    if (t < LINEAR_DURATION_S) {
      f = LINEAR_END_HZ - (LINEAR_END_HZ - LINEAR_START_HZ)
                           * static_cast<float>(t / LINEAR_DURATION_S);
      phase = phaseLinDownRad(t);
    } else {
      const double t2 = t - LINEAR_DURATION_S;
      f = LOG_END_HZ * powf(LOG_START_HZ / LOG_END_HZ,
                             static_cast<float>(t2 / LOG_DURATION_S));
      phase = phaseLogDownRad(t2);
    }
    const float amplitude = amplitudeMicrosteps(f);
    return {lround(amplitude * sin(phase)), f, "DOWN_SWEEP"};
  }
  t -= SWEEP_DURATION_S;

  if (t < TAIL_S) return {0, NAN, "TAIL"};
  return {0, NAN, "DONE"};
}

// --- Driver setup ------------------------------------------------------

void dumpRegisters(const char *when) {
  Serial.printf("# REGDUMP %s: GCONF=0x%08lX CHOPCONF=0x%08lX "
                "IHOLD_IRUN=0x%08lX TPWMTHRS=0x%08lX PWMCONF=0x%08lX\n",
                when,
                static_cast<unsigned long>(driver.GCONF()),
                static_cast<unsigned long>(driver.CHOPCONF()),
                static_cast<unsigned long>(driver.IHOLD_IRUN()),
                static_cast<unsigned long>(driver.TPWMTHRS()),
                static_cast<unsigned long>(driver.PWMCONF()));
}

bool configureDriver() {
  DriverSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();

  // Fixed chopper timing -- tune once for this motor, then do not touch.
  // NOT COMPILE-VERIFIED: hstrt()/hend()/tbl() are TMCStepper's field-name
  // accessors as remembered, not confirmed against the installed library
  // version. If these fail to build, check TMCStepper.h for the actual
  // CHOPCONF field setter names before renaming blindly.
  driver.toff(5);
  driver.tbl(2);     // TBL index (blank time); v1 left this at reset default.
  driver.hstrt(5);
  driver.hend(2);

  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.rms_current(CURRENT_RMS_MA);

  // Force SpreadCycle for the whole run: both belt (threshold) and
  // suspenders (explicit enable), so no mode switch can occur mid-sweep.
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);

  // Standstill current must not drop -- the excitation is centered there.
  driver.ihold(driver.irun());
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  driver.intpol(false);

  uint8_t connection = 2;
  for (uint8_t attempt = 1; attempt <= 10; ++attempt) {
    connection = driver.test_connection();
    Serial.printf("# TMC connection attempt %u: %u (0 means OK)\n",
                  attempt, connection);
    if (connection == 0) break;
    delay(200);
  }
  if (connection != 0) return false;

  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  uint32_t chopconf = driver.CHOPCONF();
  chopconf = (chopconf & ~MRES_MASK)
             | (static_cast<uint32_t>(MRES_REGISTER_CODE) << 24);
  driver.CHOPCONF(chopconf);
  const uint8_t readbackCode = static_cast<uint8_t>(
      (driver.CHOPCONF() & MRES_MASK) >> 24);
  if (readbackCode != MRES_REGISTER_CODE) {
    Serial.printf("# MRES readback failed: requested code=%u readback=%u\n",
                  MRES_REGISTER_CODE, readbackCode);
    return false;
  }

  driverConfigured = true;
  return true;
}

// --- Fast STEP pulse (direct GPIO register write) -------------------------

// NOT COMPILE-VERIFIED: GPIO.out_w1ts/out_w1tc direct assignment is the
// common fast-GPIO idiom for ESP32 Arduino cores, but the exact field type
// in soc/gpio_struct.h has varied release to release. If this does not
// compile as a direct assignment, try GPIO.out_w1ts.val = ... instead.
inline void stepPulseFast() {
  GPIO.out_w1ts = (1UL << STEP_PIN);
  delayMicroseconds(STEP_HIGH_US);
  GPIO.out_w1tc = (1UL << STEP_PIN);
}

bool stepToward(int32_t target, const char *segment, float frequencyHz) {
  if (aborted) return false;
  if (target == currentMicrosteps) return true;

  if (abs(target) > MAX_ABS_TARGET_MICROSTEPS) {
    Serial.printf("# Safety abort: target %ld microsteps exceeds sanity "
                  "ceiling %ld (segment=%s, f=%.3f Hz)\n",
                  static_cast<long>(target),
                  static_cast<long>(MAX_ABS_TARGET_MICROSTEPS), segment,
                  static_cast<double>(frequencyHz));
    aborted = true;
    return false;
  }

  const int8_t direction = (target > currentMicrosteps) ? 1 : -1;
  const int64_t nowUs = esp_timer_get_time();
  if (nowUs - lastEdgeUs < MIN_STEP_PERIOD_US) return true;  // rate-limited

  if (direction != lastDirection) {
    digitalWrite(DIR_PIN, direction > 0 ? HIGH : LOW);
    delayMicroseconds(DIR_SETUP_US);
    lastDirection = direction;
  }
  stepPulseFast();
  currentMicrosteps += direction;
  lastEdgeUs = esp_timer_get_time();
  return true;
}

void returnToOrigin() {
  const bool previousAbort = aborted;
  aborted = false;
  while (currentMicrosteps != 0) {
    const int8_t direction = (currentMicrosteps > 0) ? -1 : 1;
    if (direction != lastDirection) {
      digitalWrite(DIR_PIN, direction > 0 ? HIGH : LOW);
      delayMicroseconds(DIR_SETUP_US);
      lastDirection = direction;
    }
    while (esp_timer_get_time() - lastEdgeUs < MIN_STEP_PERIOD_US) {}
    stepPulseFast();
    currentMicrosteps += direction;
    lastEdgeUs = esp_timer_get_time();
  }
  aborted = previousAbort;
}

// --- Top-level run -------------------------------------------------------

void printPlan() {
  Serial.println("# v2 chirp plan: cruise amplitude + smooth notch (no e_max)");
  Serial.printf("# MRES=%u (code %u) current_rms_mA=%u\n", MRES,
                MRES_REGISTER_CODE, CURRENT_RMS_MA);
  Serial.printf("# cruise=%.2f full steps (%.1f microsteps)\n",
                static_cast<double>(CRUISE_FULL_STEPS),
                static_cast<double>(CRUISE_MICROSTEPS));
  Serial.printf("# notch: f_n_placeholder=%.1f Hz halfwidth=%.1f Hz "
                "depth_ratio=%.2f -> band %.1f-%.1f Hz\n",
                static_cast<double>(F_N_HZ_PLACEHOLDER),
                static_cast<double>(NOTCH_HALFWIDTH_HZ),
                static_cast<double>(NOTCH_DEPTH_RATIO),
                static_cast<double>(F_N_HZ_PLACEHOLDER - NOTCH_HALFWIDTH_HZ),
                static_cast<double>(F_N_HZ_PLACEHOLDER + NOTCH_HALFWIDTH_HZ));
  Serial.printf("# sweep: log %.0f->%.0f Hz in %.0f s, linear %.0f->%.0f Hz "
                "in %.0f s, each direction\n",
                static_cast<double>(LOG_START_HZ),
                static_cast<double>(LOG_END_HZ), LOG_DURATION_S,
                static_cast<double>(LINEAR_START_HZ),
                static_cast<double>(LINEAR_END_HZ), LINEAR_DURATION_S);
  Serial.printf("# level duration: %.1f s (%.2f min)\n", LEVEL_DURATION_S,
                LEVEL_DURATION_S / 60.0);
  Serial.printf("# step_rate_clamp_hz=%.0f (UNVERIFIED on this firmware -- "
                "bench-test before trusting)\n",
                static_cast<double>(STEP_RATE_CLAMP_HZ));
  Serial.printf("# TMC connection/configuration: %s; motor remains disabled\n",
                driverConfigured ? "OK" : "FAILED");
}

void runChirp() {
  if (!driverConfigured) {
    Serial.println("# RUN rejected: TMC2209 configuration did not pass.");
    return;
  }
  if (currentMicrosteps != 0) {
    Serial.println("# RUN rejected: logical position is not at origin.");
    return;
  }

  dumpRegisters("PRE_RUN");
  aborted = false;
  lastDirection = 0;
  lastEdgeUs = esp_timer_get_time();
  digitalWrite(EN_PIN, LOW);
  delay(50);
  digitalWrite(TRIG_OUT_PIN, HIGH);
  delayMicroseconds(5);
  logEvent("CHIRP_START", "FULL", 0.0F);

  const int64_t startUs = esp_timer_get_time();
  const char *lastSegment = "";
  while (!aborted) {
    const double elapsedS = (esp_timer_get_time() - startUs) * 1.0e-6;
    if (elapsedS >= LEVEL_DURATION_S) break;

    const CommandState state = commandedStateAt(elapsedS);
    if (strcmp(state.segment, lastSegment) != 0) {
      logEvent("SEGMENT_START", state.segment, state.frequencyHz);
      lastSegment = state.segment;
    }
    stepToward(state.targetMicrosteps, state.segment, state.frequencyHz);
    pollAbort();
  }

  returnToOrigin();
  digitalWrite(TRIG_OUT_PIN, LOW);
  delayMicroseconds(5);
  logEvent(aborted ? "CHIRP_ABORTED" : "CHIRP_COMPLETE", "FULL", 0.0F);
  dumpRegisters("POST_RUN");
  digitalWrite(EN_PIN, HIGH);
  Serial.println(aborted ? "# Chirp aborted; motor disabled."
                          : "# Chirp complete; motor disabled.");
}

}  // namespace chirp_v2

using namespace chirp_v2;

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  pinMode(EN_PIN, OUTPUT);
  pinMode(TRIG_OUT_PIN, OUTPUT);
  pinMode(TRIG_ECHO_PIN, INPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);
  digitalWrite(EN_PIN, HIGH);
  digitalWrite(TRIG_OUT_PIN, LOW);

  Serial.begin(CONSOLE_BAUD);
  Serial.setTimeout(20);
  const uint32_t waitStart = millis();
  while (!Serial && millis() - waitStart < 3000) delay(10);
  Serial.println("# ESP32-S3/TMC2209 v2 chirp (cruise + smooth notch) -- "
                  "UNFLASHED DESIGN, confirm pins and bench-test timing "
                  "before a real run");
  driverConfigured = configureDriver();
  digitalWrite(EN_PIN, HIGH);
  Serial.println("timestamp_us,event,segment,frequency_hz,position_microsteps");
  Serial.println("# Enter CHECK or RUN. Send ABORT during motion.");
}

void loop() {
  if (!Serial.available()) {
    delay(2);
    return;
  }
  String command = Serial.readStringUntil('\n');
  command.trim();
  command.toUpperCase();
  if (command == "CHECK") {
    printPlan();
  } else if (command == "RUN") {
    runChirp();
  } else if (command == "ABORT") {
    aborted = true;
  } else {
    Serial.println("# Unknown command. Use CHECK, RUN, or ABORT.");
  }
}

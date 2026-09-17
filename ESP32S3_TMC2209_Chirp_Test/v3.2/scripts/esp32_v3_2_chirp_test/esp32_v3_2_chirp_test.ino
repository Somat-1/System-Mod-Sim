/*
 * ESP32-S3 + TMC2209 v3.2 resonance-band chirp: full sweep, aggressive
 * notch. Dual-core step generation.
 *
 * NOT the recommended default -- see
 * ESP32S3_TMC2209_Chirp_Test/v3/scripts/esp32_v3_chirp_test instead, which
 * replaces this file's notch with a genuine 150-250 Hz frequency gap (zero
 * excitation, not reduced amplitude). This file exists for reference/
 * comparison: it is close to the configuration that was actually running
 * during the second 2026-09-15 hardware incident (uncontrolled motor
 * motion, reported live around 200 Hz) -- full 1-1000 Hz sweep, notch
 * widened to a 65-285 Hz band at 97% suppression at its center. That
 * amplitude reduction, however aggressive, was not sufficient to prevent
 * the incident. Treat any run of this file with the same caution as that
 * event: stay at the rig, ready to cut power, and expect that 150-250 Hz
 * may still cause uncontrolled motion despite the notch.
 *
 * Forked from
 * ESP32S3_TMC2209_Chirp_Test/v2/scripts/esp32_v2_chirp_test_rsense_fixed,
 * which is v2's originally-validated sequence with R_SENSE corrected to
 * 0.11 -- that correction, the indefinite TMC-connection retry, the 7 s
 * auto-start countdown, and the LED phase indicator all carry over
 * unchanged. Only the notch constants differ from that baseline:
 *
 *   F_N_HZ_PLACEHOLDER   176.7 -> 175.0  (measured, not the model value)
 *   NOTCH_HALFWIDTH_HZ    56.7 -> 110.0  (band 120-233 Hz -> 65-285 Hz)
 *   NOTCH_DEPTH_RATIO      0.10 -> 0.03  (90% -> 97% suppression at center)
 *
 * Everything else -- the full log(1-60 Hz)/linear(60-1000 Hz) sweep law,
 * cruise amplitude, current, MRES -- is unchanged from v2.
 *
 * Serial commands at 115200 baud: CHECK, RUN, ABORT.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <esp_timer.h>
#include <esp_rom_sys.h>
#include <esp_task_wdt.h>
#include <cmath>
#include <cstring>
#include "soc/gpio_struct.h"

namespace chirp_v3_2 {

// --- Pin assignment: confirmed against the working v4 MRES trajectory
// campaign sketch's wiring for this rig. ---
constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;  // Onboard RGB LED; same pin as v4/v5/v6.
// No EN_PIN: EN/ENN is physically grounded on this rig (always enabled).
static_assert(STEP_PIN < 32, "STEP_PIN >= 32 needs GPIO.out1_*");
static_assert(DIR_PIN < 32, "DIR_PIN >= 32 needs GPIO.out1_*");

constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.11F;  // Actual resistor on this board; was 0.03.
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;

// --- Mirrors chirp_v2_schedule.py -- keep these in sync by hand ----------

constexpr uint16_t MRES = 16;               // 1/16 step. See README Section 2.
constexpr uint8_t MRES_REGISTER_CODE = 4;   // TMC2209 MRES field code for 1/16.
constexpr uint16_t CURRENT_RMS_MA = 400;

// Widened/deepened from v2's 176.7/56.7/0.10 -- see file header. Verified
// on real hardware NOT sufficient to prevent uncontrolled motion around
// 200 Hz; kept here for reference/comparison only.
constexpr float F_N_HZ_PLACEHOLDER = 175.0F;
constexpr float NOTCH_HALFWIDTH_HZ = 110.0F;
constexpr float NOTCH_DEPTH_RATIO = 0.03F;

constexpr float CRUISE_FULL_STEPS = 3.0F;
constexpr float CRUISE_MICROSTEPS = CRUISE_FULL_STEPS * MRES;  // 48

constexpr float STROKE_CLAMP_MM = 0.5F;
constexpr float FULL_STEP_MM = 0.010F;  // 2 mm lead / 200 full steps/rev.
constexpr float STROKE_CLAMP_MICROSTEPS =
    (STROKE_CLAMP_MM / FULL_STEP_MM) * MRES;  // 800; inactive at cruise=3.

// UNVERIFIED -- this is a peak MICROSTEP PULSE rate, not a full-step rate.
// Benchmarking showed the dual-core architecture faithfully tracks the
// intended waveform well below this ceiling at MRES=16, so this clamp is
// conservative rather than binding in practice.
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

constexpr uint32_t AUTO_RUN_COUNTDOWN_MS = 7000;  // Margin once TMC is detected.
constexpr uint32_t CONNECTION_RETRY_DELAY_MS = 1000;

// --- Timing floors -----------------------------------------------------

constexpr uint32_t STEP_HIGH_US = 1;
constexpr uint32_t DIR_SETUP_US = 2;
constexpr int64_t MIN_STEP_PERIOD_US = 3;  // benchmarked: not the binding
                                            // constraint once split across
                                            // both cores at MRES=16.

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);

// --- Shared state between core 1 (setup/loop, publishes targetMicrosteps)
// and core 0 (stepTaskFn, owns everything else). Single 32-bit-or-smaller
// volatile reads/writes are atomic on Xtensa, so no lock is needed for
// this producer/consumer pattern. ---
volatile int32_t targetMicrosteps = 0;
volatile int32_t currentMicrosteps = 0;
volatile int8_t lastDirection = 0;
volatile bool stepTaskRunning = false;
volatile bool aborted = false;

bool driverConfigured = false;
TaskHandle_t stepTaskHandle = nullptr;

// Phase indicator, since serial has proven unreliable to watch live on this
// rig: yellow=waiting for/configuring the TMC2209, blue=detected, counting
// down to auto-start, green=running, cyan=done, red=failed/aborted.
void colour(uint8_t red, uint8_t green, uint8_t blue) {
  rgbLedWrite(LED_PIN, red, green, blue);
}

void logEvent(const char *event, const char *segment, float frequencyHz) {
  Serial.printf("%lld,%s,%s,%.6f,%ld\n",
                static_cast<long long>(esp_timer_get_time()), event, segment,
                static_cast<double>(frequencyHz),
                static_cast<long>(currentMicrosteps));
}

// Accumulates characters across as many calls as it takes -- no timeout,
// so it works regardless of typing speed (readStringUntil's short timeout
// was truncating slowly-typed commands into unrecognized single
// characters before a full line ever arrived).
String serialLineBuffer;

bool readSerialLine(String &outLine) {
  while (Serial.available() > 0) {
    const char c = static_cast<char>(Serial.read());
    if (c == '\r') continue;
    if (c == '\n') {
      outLine = serialLineBuffer;
      serialLineBuffer = "";
      outLine.trim();
      outLine.toUpperCase();
      return true;
    }
    serialLineBuffer += c;
  }
  return false;
}

void pollAbort() {
  String command;
  while (readSerialLine(command)) {
    if (command == "ABORT") {
      aborted = true;
      logEvent("ABORT_REQUEST", "SERIAL", 0.0F);
    }
  }
}

// --- Core 0: dedicated step-issuing task. No floating point, ever. -------
void stepTaskFn(void *pvParameters) {
  (void)pvParameters;
  int64_t lastEdgeUs = esp_timer_get_time();
  for (;;) {
    if (!stepTaskRunning) {
      vTaskDelay(1);
      lastEdgeUs = esp_timer_get_time();
      continue;
    }
    const int32_t target = targetMicrosteps;
    const int32_t current = currentMicrosteps;
    if (target != current) {
      const int64_t now = esp_timer_get_time();
      if (now - lastEdgeUs >= MIN_STEP_PERIOD_US) {
        const int8_t direction = (target > current) ? 1 : -1;
        if (direction != lastDirection) {
          digitalWrite(DIR_PIN, direction > 0 ? HIGH : LOW);
          esp_rom_delay_us(DIR_SETUP_US);
          lastDirection = direction;
        }
        GPIO.out_w1ts = (1UL << STEP_PIN);
        esp_rom_delay_us(STEP_HIGH_US);
        GPIO.out_w1tc = (1UL << STEP_PIN);
        currentMicrosteps = current + direction;
        lastEdgeUs = now;
      }
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
  driver.toff(5);
  driver.tbl(2);
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

  // GCONF says what was asked for; DRV_STATUS bit 30 says what the chopper
  // is actually doing. This project has hit the gap between those twice
  // before (see the v4 MRES campaign's history: "Force SpreadCycle...it was
  // running StealthChop by default") -- the chirp codebase never carried
  // that verification over until now. Check both, refuse to run on a
  // mismatch, rather than trust the write.
  const uint32_t gconf = driver.GCONF();
  const bool spreadCycleRequested = ((gconf >> 2) & 1U) != 0U;
  Serial.printf("# CHOPPER_CONFIG,gconf=0x%08lX,en_spreadcycle_requested=%u\n",
                static_cast<unsigned long>(gconf),
                static_cast<unsigned>(spreadCycleRequested));
  if (!spreadCycleRequested) {
    Serial.println("# CHOPPER_CONFIG_FAILED,refusing_to_run");
    return false;
  }
  const uint32_t drvStatus = driver.DRV_STATUS();
  const bool stealthActive = ((drvStatus >> 30) & 1U) != 0U;
  Serial.printf("# CHOPPER_ACTIVE,drv_status=0x%08lX,stealth_bit=%u\n",
                static_cast<unsigned long>(drvStatus),
                static_cast<unsigned>(stealthActive));
  if (stealthActive) {
    Serial.println("# STEALTHCHOP_STILL_ACTIVE,refusing_to_run");
    return false;
  }

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

// --- Top-level run -------------------------------------------------------

void printPlan() {
  Serial.println("# v3.2 chirp plan: full sweep, AGGRESSIVE NOTCH (not the "
                  "recommended default -- see v3's hard gap)");
  Serial.println("# dual-core step generation: core 0 issues edges, core 1 "
                  "computes the target");
  Serial.printf("# MRES=%u (code %u) current_rms_mA=%u (R_SENSE_OHM=%.2f)\n",
                MRES, MRES_REGISTER_CODE, CURRENT_RMS_MA,
                static_cast<double>(R_SENSE_OHM));
  Serial.printf("# cruise=%.2f full steps (%.1f microsteps)\n",
                static_cast<double>(CRUISE_FULL_STEPS),
                static_cast<double>(CRUISE_MICROSTEPS));
  Serial.printf("# notch: f_n=%.1f Hz halfwidth=%.1f Hz depth_ratio=%.2f -> "
                "band %.1f-%.1f Hz (97%% suppressed at center; NOT verified "
                "safe -- see file header)\n",
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
  Serial.printf("# step_rate_clamp_hz=%.0f (conservative; dual-core "
                "benchmarking tracked the intended waveform well below "
                "this at MRES=%u)\n",
                static_cast<double>(STEP_RATE_CLAMP_HZ), MRES);
  Serial.printf("# TMC connection/configuration: %s; EN is hardwired -- "
                "driver stays enabled as long as it is powered\n",
                driverConfigured ? "OK" : "FAILED");
}

void returnToOrigin() {
  const bool previousAbort = aborted;
  aborted = false;
  targetMicrosteps = 0;
  stepTaskRunning = true;
  while (currentMicrosteps != 0) {
    delay(1);
  }
  aborted = previousAbort;
}

void runChirp() {
  if (!driverConfigured) {
    colour(32, 0, 0);
    Serial.println("# RUN rejected: TMC2209 configuration did not pass.");
    return;
  }
  if (currentMicrosteps != 0) {
    colour(32, 0, 0);
    Serial.println("# RUN rejected: logical position is not at origin.");
    return;
  }

  dumpRegisters("PRE_RUN");
  aborted = false;
  lastDirection = 0;
  stepTaskRunning = true;
  logEvent("CHIRP_START", "FULL", 0.0F);
  colour(0, 32, 0);  // green: running.

  const int64_t startUs = esp_timer_get_time();
  const char *lastSegment = "";
  while (!aborted) {
    const double elapsedS = (esp_timer_get_time() - startUs) * 1.0e-6;
    if (elapsedS >= LEVEL_DURATION_S) break;

    const CommandState state = commandedStateAt(elapsedS);
    if (abs(state.targetMicrosteps) > MAX_ABS_TARGET_MICROSTEPS) {
      Serial.printf("# Safety abort: target %ld microsteps exceeds sanity "
                    "ceiling %ld (segment=%s, f=%.3f Hz)\n",
                    static_cast<long>(state.targetMicrosteps),
                    static_cast<long>(MAX_ABS_TARGET_MICROSTEPS),
                    state.segment, static_cast<double>(state.frequencyHz));
      aborted = true;
      colour(32, 0, 0);
      break;
    }
    if (strcmp(state.segment, lastSegment) != 0) {
      logEvent("SEGMENT_START", state.segment, state.frequencyHz);
      lastSegment = state.segment;
    }
    targetMicrosteps = state.targetMicrosteps;
    pollAbort();
  }
  if (aborted) colour(32, 0, 0);

  returnToOrigin();
  logEvent(aborted ? "CHIRP_ABORTED" : "CHIRP_COMPLETE", "FULL", 0.0F);
  dumpRegisters("POST_RUN");
  stepTaskRunning = false;
  colour(aborted ? 32 : 0, aborted ? 0 : 32, aborted ? 0 : 32);  // red / cyan.
  Serial.println(aborted
                      ? "# Chirp aborted; motor returned to origin. EN is "
                        "hardwired, driver remains enabled."
                      : "# Chirp complete; motor returned to origin. EN is "
                        "hardwired, driver remains enabled.");
}

}  // namespace chirp_v3_2

using namespace chirp_v3_2;

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);

  colour(24, 24, 24);  // One white boot blink.
  delay(500);
  colour(0, 0, 0);

  Serial.begin(CONSOLE_BAUD);
  Serial.setTimeout(20);
  const uint32_t waitStart = millis();
  while (!Serial && millis() - waitStart < 3000) delay(10);
  Serial.println("# ESP32-S3/TMC2209 v3.2 chirp (full sweep, aggressive "
                  "notch -- NOT the recommended default, see v3)");

  esp_task_wdt_deinit();
  xTaskCreatePinnedToCore(stepTaskFn, "stepTask", 4096, nullptr, 1,
                          &stepTaskHandle, 0);  // pin to core 0.
  delay(50);

  Serial.println("# Waiting for TMC2209 (retrying every ~1 s; connect/power "
                  "it now if not already). LED is yellow while waiting.");
  colour(24, 24, 0);  // yellow: waiting for / configuring the driver.
  uint32_t outerAttempt = 0;
  while (!(driverConfigured = configureDriver())) {
    ++outerAttempt;
    Serial.printf("# Still waiting for TMC2209 (attempt %lu)...\n",
                  static_cast<unsigned long>(outerAttempt));
    delay(CONNECTION_RETRY_DELAY_MS);
  }
  Serial.println("# TMC2209 detected and configured.");
  printPlan();
  Serial.println("timestamp_us,event,segment,frequency_hz,position_microsteps");

  colour(0, 0, 32);  // blue: detected, counting down to auto-start.
  Serial.printf("# Starting automatically in %lu ms. Send ABORT now to cancel "
                "(RUN will still work manually afterward).\n",
                static_cast<unsigned long>(AUTO_RUN_COUNTDOWN_MS));
  const uint32_t countdownStart = millis();
  while (millis() - countdownStart < AUTO_RUN_COUNTDOWN_MS) {
    pollAbort();
    if (aborted) {
      Serial.println("# Auto-start cancelled. Send RUN manually when ready.");
      colour(0, 0, 0);
      return;
    }
    delay(20);
  }
  runChirp();
}

void loop() {
  String command;
  if (!readSerialLine(command)) {
    delay(2);
    return;
  }
  if (command == "CHECK") {
    printPlan();
  } else if (command == "RUN") {
    runChirp();
  } else if (command == "ABORT") {
    aborted = true;
  } else if (command.length() > 0) {
    Serial.println("# Unknown command. Use CHECK, RUN, or ABORT.");
  }
}

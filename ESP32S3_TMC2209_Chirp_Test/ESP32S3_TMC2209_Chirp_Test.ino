/*
 * ESP32-S3 + TMC2209 bounded alternating-step chirp.
 *
 * The commanded position toggles between the starting position and one
 * microstep, then returns after every cycle. This excites the stage without
 * accumulating one-direction travel. The sweep is 0/idle -> 1..1000 Hz up,
 * then 1000..1 Hz down -> 0/idle. A configurable resonance exclusion band
 * suppresses STEP commands around the modeled motor-dominated resonance.
 *
 * Framework: Arduino-ESP32 3.x; library: TMCStepper.
 * Serial commands at 115200 baud: CHECK, RUN, ABORT.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <esp_timer.h>

namespace chirp_test {

// Copied from the validated v2 ESP32-S3/TMC2209 runner.
constexpr uint8_t STEP_PIN = 6;
constexpr uint8_t DIR_PIN = 7;
constexpr uint8_t EN_PIN = 5;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t TRIG_OUT_PIN = 1;
constexpr uint8_t TRIG_ECHO_PIN = 2;

constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.03F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;

// Conservative excitation settings. With a 2 mm pitch and 200 full
// steps/revolution, MRES=4 corresponds to a theoretical 2.5 um increment.
constexpr uint16_t MRES = 4;
constexpr uint16_t CURRENT_RMS_MA = 360;
constexpr uint32_t STEP_HIGH_US = 2;
constexpr uint32_t DIR_SETUP_US = 5;

// True zero frequency is represented by the stationary lead/tail dwells.
constexpr float SWEEP_START_HZ = 1.0F;
constexpr float SWEEP_END_HZ = 1000.0F;
constexpr float SWEEP_DURATION_S = 30.0F;
constexpr uint32_t LEAD_IN_MS = 2000;
constexpr uint32_t TURNAROUND_MS = 500;
constexpr uint32_t TAIL_MS = 2000;

// Default safety notch around the modeled motor-dominated first resonance.
// No excitation data are produced inside this band. Adjust only after a
// low-amplitude pilot establishes the physical motor resonance.
constexpr bool RESONANCE_NOTCH_ENABLED = true;
constexpr float NOTCH_LOW_HZ = 120.0F;
constexpr float NOTCH_HIGH_HZ = 230.0F;

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);

volatile bool aborted = false;
bool driverConfigured = false;
int8_t logicalPosition = 0;  // Bounded to 0 or +1 microstep.

void logEvent(const char *event, const char *segment, float frequencyHz) {
  Serial.printf("%lld,%s,%s,%.6f,%d,%d\n",
                static_cast<long long>(esp_timer_get_time()), event, segment,
                static_cast<double>(frequencyHz), logicalPosition,
                digitalRead(TRIG_ECHO_PIN));
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

bool waitUntilUs(int64_t deadlineUs) {
  while (!aborted) {
    const int64_t remainingUs = deadlineUs - esp_timer_get_time();
    if (remainingUs <= 0) return true;
    pollAbort();
    if (remainingUs > 2500) {
      delay(1);
    } else if (remainingUs > 200) {
      delayMicroseconds(static_cast<uint32_t>(remainingUs - 100));
    }
  }
  return false;
}

bool dwellAbortable(uint32_t durationMs, const char *label) {
  const int64_t deadlineUs =
      esp_timer_get_time() + static_cast<int64_t>(durationMs) * 1000;
  logEvent("DWELL_START", label, 0.0F);
  const bool completed = waitUntilUs(deadlineUs);
  logEvent("DWELL_END", label, 0.0F);
  return completed;
}

bool configureDriver() {
  DriverSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();
  driver.toff(5);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.en_spreadCycle(true);
  driver.intpol(false);
  driver.rms_current(CURRENT_RMS_MA);
  driver.ihold(driver.irun());
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  uint8_t connection = 2;
  for (uint8_t attempt = 1; attempt <= 10; ++attempt) {
    connection = driver.test_connection();
    Serial.printf("# TMC connection attempt %u: %u (0 means OK)\n",
                  attempt, connection);
    if (connection == 0) break;
    delay(200);
  }
  if (connection != 0) return false;

  const uint8_t mresCode = MRES == 16 ? 4 : MRES == 4 ? 6 :
                           MRES == 2 ? 7 : 8;
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  uint32_t chopconf = driver.CHOPCONF();
  chopconf = (chopconf & ~MRES_MASK) |
             (static_cast<uint32_t>(mresCode) << 24);
  driver.CHOPCONF(chopconf);
  const uint8_t readbackCode = static_cast<uint8_t>(
      (driver.CHOPCONF() & MRES_MASK) >> 24);
  if (readbackCode != mresCode) {
    Serial.printf("# MRES readback failed: requested code=%u readback=%u\n",
                  mresCode, readbackCode);
    return false;
  }
  driverConfigured = true;
  return true;
}

bool stepAt(int direction, int64_t edgeDeadlineUs) {
  if (aborted) return false;
  if (!((direction > 0 && logicalPosition == 0) ||
        (direction < 0 && logicalPosition == 1))) {
    Serial.println("# Position-bound violation; aborting.");
    aborted = true;
    return false;
  }
  digitalWrite(DIR_PIN, direction > 0 ? HIGH : LOW);
  const int64_t setupDeadline = edgeDeadlineUs - DIR_SETUP_US;
  if (!waitUntilUs(setupDeadline)) return false;
  while (esp_timer_get_time() < edgeDeadlineUs) {}
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(STEP_HIGH_US);
  digitalWrite(STEP_PIN, LOW);
  logicalPosition += direction > 0 ? 1 : -1;
  return true;
}

void returnToOrigin() {
  if (logicalPosition == 1) {
    const bool previousAbort = aborted;
    aborted = false;
    stepAt(-1, esp_timer_get_time() + DIR_SETUP_US + 20);
    aborted = previousAbort;
  }
}

float sweepFrequency(float elapsedS, bool ascending) {
  const float fraction = constrain(elapsedS / SWEEP_DURATION_S, 0.0F, 1.0F);
  const float span = SWEEP_END_HZ - SWEEP_START_HZ;
  return ascending ? SWEEP_START_HZ + span * fraction
                   : SWEEP_END_HZ - span * fraction;
}

bool frequencyIsNotched(float frequencyHz) {
  return RESONANCE_NOTCH_ENABLED && frequencyHz >= NOTCH_LOW_HZ &&
         frequencyHz <= NOTCH_HIGH_HZ;
}

bool runSweep(const char *segment, bool ascending) {
  const int64_t startUs = esp_timer_get_time();
  const int64_t endUs =
      startUs + static_cast<int64_t>(SWEEP_DURATION_S * 1000000.0F);
  int64_t nextEdgeUs = startUs;
  bool notchActive = false;
  logEvent("SWEEP_START", segment,
           ascending ? SWEEP_START_HZ : SWEEP_END_HZ);

  while (!aborted) {
    const int64_t nowUs = esp_timer_get_time();
    if (nowUs >= endUs) break;
    const float elapsedS = static_cast<float>(nowUs - startUs) * 1.0e-6F;
    const float frequencyHz = sweepFrequency(elapsedS, ascending);

    if (frequencyIsNotched(frequencyHz)) {
      if (!notchActive) {
        returnToOrigin();
        notchActive = true;
        logEvent("NOTCH_ENTER", segment, frequencyHz);
      }
      pollAbort();
      delay(1);
      nextEdgeUs = esp_timer_get_time();
      continue;
    }

    if (notchActive) {
      notchActive = false;
      logEvent("NOTCH_EXIT", segment, frequencyHz);
      nextEdgeUs = esp_timer_get_time();
    }

    const int64_t halfPeriodUs = static_cast<int64_t>(
        llround(500000.0 / static_cast<double>(frequencyHz)));
    const int64_t minimumIntervalUs = DIR_SETUP_US + STEP_HIGH_US + 2;
    nextEdgeUs += (
        halfPeriodUs > minimumIntervalUs ? halfPeriodUs : minimumIntervalUs);
    if (nextEdgeUs >= endUs) break;
    const int direction = logicalPosition == 0 ? 1 : -1;
    if (!stepAt(direction, nextEdgeUs)) break;
  }

  returnToOrigin();
  const float finalFrequency = ascending ? SWEEP_END_HZ : SWEEP_START_HZ;
  logEvent(aborted ? "SWEEP_ABORTED" : "SWEEP_END", segment,
           finalFrequency);
  return !aborted;
}

void printPlan() {
  Serial.println("# Bounded alternating-step chirp plan");
  Serial.printf("# MRES=%u current_rms_mA=%u theoretical_increment_um=2.5\n",
                MRES, CURRENT_RMS_MA);
  Serial.printf("# stationary %lu ms; up %.1f -> %.1f Hz in %.1f s\n",
                static_cast<unsigned long>(LEAD_IN_MS),
                static_cast<double>(SWEEP_START_HZ),
                static_cast<double>(SWEEP_END_HZ),
                static_cast<double>(SWEEP_DURATION_S));
  Serial.printf("# down %.1f -> %.1f Hz in %.1f s; stationary %lu ms\n",
                static_cast<double>(SWEEP_END_HZ),
                static_cast<double>(SWEEP_START_HZ),
                static_cast<double>(SWEEP_DURATION_S),
                static_cast<unsigned long>(TAIL_MS));
  Serial.printf("# resonance_notch=%s band=%.1f..%.1f Hz\n",
                RESONANCE_NOTCH_ENABLED ? "ON" : "OFF",
                static_cast<double>(NOTCH_LOW_HZ),
                static_cast<double>(NOTCH_HIGH_HZ));
  Serial.println("# position bound: 0 to +1 configured microstep; net zero");
  Serial.printf("# TMC connection/configuration: %s; motor remains disabled\n",
                driverConfigured ? "OK" : "FAILED");
}

void runChirp() {
  if (!driverConfigured) {
    Serial.println("# RUN rejected: TMC2209 configuration did not pass.");
    return;
  }
  if (logicalPosition != 0) {
    Serial.println("# RUN rejected: logical position is not at origin.");
    return;
  }

  aborted = false;
  digitalWrite(EN_PIN, LOW);
  delay(50);
  digitalWrite(TRIG_OUT_PIN, HIGH);
  delayMicroseconds(5);
  logEvent("CHIRP_START", "FULL", 0.0F);

  bool ok = dwellAbortable(LEAD_IN_MS, "LEAD_IN");
  if (ok) ok = runSweep("UP", true);
  if (ok) ok = dwellAbortable(TURNAROUND_MS, "TURNAROUND");
  if (ok) ok = runSweep("DOWN", false);
  if (ok) ok = dwellAbortable(TAIL_MS, "TAIL");

  returnToOrigin();
  digitalWrite(TRIG_OUT_PIN, LOW);
  delayMicroseconds(5);
  logEvent(ok ? "CHIRP_COMPLETE" : "CHIRP_ABORTED", "FULL", 0.0F);
  digitalWrite(EN_PIN, HIGH);
  Serial.println(ok ? "# Chirp complete; motor disabled."
                    : "# Chirp aborted; motor disabled.");
}

}  // namespace chirp_test

using namespace chirp_test;

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
  Serial.println("# ESP32-S3/TMC2209 bounded chirp test");
  driverConfigured = configureDriver();
  digitalWrite(EN_PIN, HIGH);
  Serial.println("timestamp_us,event,segment,frequency_hz,position_microsteps,trigger_echo");
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

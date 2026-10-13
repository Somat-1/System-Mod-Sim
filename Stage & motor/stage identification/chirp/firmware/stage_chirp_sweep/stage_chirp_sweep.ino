/*
 * Mechanical chirp sweep — ESP32-S3 + TMC2209, two selectable profiles.
 *
 *   RUN FULL   1 Hz -> 1 kHz and back, with a hard 150-250 Hz gap
 *   RUN HIGH   250 Hz -> 1 kHz and back, continuous, started from a dwell
 *   CHECK      configure and report, move nothing
 *   ABORT      stop now (also cancels the start countdown)
 *
 * Structure and constants are carried over from the hardware-validated
 * ESP32S3_TMC2209_Chirp_Test/v3; see ../scripts/chirp_design.py for the
 * profiles as data and ../README.md for why there are two.
 *
 * ---------------------------------------------------------------------------
 * SAFETY — 150 to 250 Hz
 *
 * Two incidents on 2026-09-15 produced uncontrolled motor motion, reported live
 * at around 200 Hz. A smooth amplitude notch did not prevent them: not at 90%
 * suppression, and not widened to 65-285 Hz at 97%. Attenuating that band was
 * not sufficient however aggressively it was done.
 *
 * So FULL does not attenuate it. It does not visit it. The up sweep runs
 * 1->150 Hz, changes target to 250 Hz in one step, and continues 250->1000 Hz;
 * the down sweep mirrors that. Do not "close the gap to get the missing data".
 *
 * Stay at the rig with the power within reach for any run of FULL.
 * ---------------------------------------------------------------------------
 *
 * Step generation is split across both cores, which is not decoration. Three
 * approaches were benchmarked on this board for v3: a polled single-core loop
 * flatlined near 55 kHz because the waveform maths ran inline with GPIO
 * issuance. Here core 1 computes the target position from the chirp law and
 * core 0 does nothing but drive STEP toward whatever target it last saw,
 * rate-limited by MIN_STEP_PERIOD_US.
 *
 * WIRING — ../../../CONNECTIONS.md. STEP GPIO5, DIR GPIO6, TMC UART TX 17 /
 * RX 18, EN grounded, LED GPIO48. The encoder is not read here: at 1 kHz the
 * useful measurement is an accelerometer, and the screw encoder cannot resolve
 * a 3 full-step oscillation meaningfully at that rate.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <esp_timer.h>
#include <esp_rom_sys.h>
#include <atomic>
#include <cmath>

constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr uint8_t DRIVER_ADDRESS = 0;

/* R_SENSE DIFFERS BETWEEN BOARDS ON THIS RIG — CHECK BEFORE RUNNING.
 * 0.11 is what the v6 campaign and everything else under "Stage & motor" uses,
 * and is the default here for consistency. The v3 chirp used 0.10 for a board
 * marked R100, noting its 0.11 predecessor had died. Wrong value = wrong coil
 * current by roughly 10%, which at 400 mA is enough to change what the motor
 * does at resonance. Read the markings on the board actually fitted. */
constexpr float R_SENSE_OHM = 0.11F;

constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t MRES = 16;
constexpr uint16_t CURRENT_RMS_MA = 400;

constexpr float CRUISE_FULL_STEPS = 3.0F;
constexpr float CRUISE_MICROSTEPS = CRUISE_FULL_STEPS * MRES;     // 48
constexpr float STROKE_CLAMP_MICROSTEPS = 800.0F;                 // 0.5 mm
constexpr uint32_t MIN_STEP_PERIOD_US = 5;                        // 200 kHz ceiling

constexpr float SWEEP_LOW_START_HZ = 1.0F;
constexpr float SWEEP_GAP_LOW_HZ = 150.0F;
constexpr float SWEEP_GAP_HIGH_HZ = 250.0F;
constexpr float SWEEP_HIGH_END_HZ = 1000.0F;
constexpr double SWEEP_RATE_HZ_PER_S = 7.8125;

constexpr double PRE_ROLL_S = 30.0;
constexpr double TAIL_S = 5.0;
constexpr double MID_DWELL_S = 2.0;
constexpr double SETTLE_S = 3.0;
constexpr double SYNC_BURST_S = 0.1;
constexpr double SYNC_GAP_S = 0.1;
constexpr int SYNC_BURSTS = 3;
constexpr float SYNC_FREQ_HZ = 750.0F;
constexpr float SYNC_AMPLITUDE_MICROSTEPS = 2.0F * MRES;
constexpr uint32_t COUNTDOWN_S = 7;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

// Shared between the cores: core 1 writes, core 0 reads. Only ever one writer
// and one reader, and the value is a plain integer, so relaxed ordering on an
// atomic is enough - no lock, which must not exist on the step path.
std::atomic<int32_t> targetMicrosteps{0};
std::atomic<bool> stepTaskRun{false};
volatile bool abortRequested = false;
int32_t actualMicrosteps = 0;
bool driverReady = false;

enum Phase : uint8_t { IDLE, COUNTDOWN, RUNNING, DONE, FAILED };
Phase phase = IDLE;

// ------------------------------------------------------------------- LED
void led(uint8_t r, uint8_t g, uint8_t b) {
#ifdef RGB_BUILTIN
  neopixelWrite(LED_PIN, r, g, b);
#else
  digitalWrite(LED_PIN, (r || g || b) ? HIGH : LOW);
#endif
}

// ---------------------------------------------------------------- driver
template <typename Reader>
bool readTwice(Reader read, uint32_t &out) {
  for (uint8_t a = 0; a < 4; ++a) {
    const uint32_t x = read(), y = read();
    if (x == y && x != 0xFFFFFFFFUL) { out = x; return true; }
    delay(5);
  }
  return false;
}

bool configureDriver() {
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(20);
  driver.begin();
  driver.toff(5);
  driver.blank_time(24);
  driver.rms_current(CURRENT_RMS_MA);
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);        // resets to ON; would smear every microstep
  driver.pwm_autoscale(false);
  driver.microsteps(MRES);

  uint32_t g = 0, c = 0, d = 0;
  if (!readTwice([] { return driver.GCONF(); }, g) ||
      !readTwice([] { return driver.CHOPCONF(); }, c) ||
      !readTwice([] { return driver.DRV_STATUS(); }, d)) {
    Serial.println("# ERR,reason=uart_readback_failed");
    return false;
  }
  const bool spread = g & (1UL << 2), intpol = c & (1UL << 28),
             stealth = d & (1UL << 30);
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,drv_status=0x%08lX,"
                "spreadcycle=%u,intpol=%u,stealth=%u,r_sense=%.3f\n",
                (unsigned long)g, (unsigned long)c, (unsigned long)d,
                spread, intpol, stealth, R_SENSE_OHM);
  if (!spread || intpol || stealth) {
    Serial.println("# ERR,reason=driver_mode_not_confirmed");
    return false;
  }
  return true;
}

// ------------------------------------------------------- core 0 step task
void stepTask(void *) {
  uint64_t lastEdge = 0;
  for (;;) {
    if (!stepTaskRun.load(std::memory_order_relaxed)) {
      vTaskDelay(1);
      continue;
    }
    const int32_t want = targetMicrosteps.load(std::memory_order_relaxed);
    if (want != actualMicrosteps) {
      const uint64_t now = esp_timer_get_time();
      if (now - lastEdge >= MIN_STEP_PERIOD_US) {
        const bool positive = want > actualMicrosteps;
        digitalWrite(DIR_PIN, positive ? HIGH : LOW);
        // DIR setup is satisfied by MIN_STEP_PERIOD_US since the previous edge.
        digitalWrite(STEP_PIN, HIGH);
        esp_rom_delay_us(2);
        digitalWrite(STEP_PIN, LOW);
        actualMicrosteps += positive ? 1 : -1;
        lastEdge = now;
      }
    }
    // Yield rarely: this task owns core 0 while a sweep is running, and the
    // watchdog only needs to see it occasionally.
    static uint32_t spins = 0;
    if (((++spins) & 0xFFFFU) == 0) vTaskDelay(1);
  }
}

// -------------------------------------------------------- the sweep itself

/** Amplitude in microsteps, clamped by stroke and by achievable step rate. */
float amplitudeMicrosteps(float f_hz) {
  float amp = CRUISE_MICROSTEPS;
  if (amp > STROKE_CLAMP_MICROSTEPS) amp = STROKE_CLAMP_MICROSTEPS;
  // A sine of amplitude A at f needs a peak rate of 2*pi*f*A microsteps/s.
  const float maxRate = 1e6F / MIN_STEP_PERIOD_US;
  const float needed = 2.0F * PI * f_hz * amp;
  if (needed > maxRate) amp = maxRate / (2.0F * PI * f_hz);
  return amp;
}

bool pollAbort() {
  static String buf;
  while (Serial.available()) {
    const char ch = (char)Serial.read();
    if (ch == '\n' || ch == '\r') {
      buf.trim();
      if (buf.equalsIgnoreCase("ABORT")) {
        abortRequested = true;
        buf = "";
        return true;
      }
      buf = "";
    } else if (buf.length() < 60) {
      buf += ch;
    }
  }
  return abortRequested;
}

/** Drive a sine at a frequency that ramps linearly from f0 to f1. */
bool sweepSegment(const char *label, double f0, double f1, double duration_s) {
  Serial.printf("# SEG,%s,f0=%.2f,f1=%.2f,dur=%.3f\n", label, f0, f1, duration_s);
  const uint64_t t0 = esp_timer_get_time();
  double phaseRad = 0.0;
  uint64_t lastUs = t0;
  uint32_t report = 0;
  for (;;) {
    const uint64_t now = esp_timer_get_time();
    const double t = (now - t0) * 1e-6;
    if (t >= duration_s) break;
    const double dt = (now - lastUs) * 1e-6;
    lastUs = now;
    const double f = (duration_s > 0) ? (f0 + (f1 - f0) * (t / duration_s)) : f0;
    // Integrate phase rather than evaluating sin(2*pi*f*t): with f changing,
    // the latter is not the signal it looks like - its instantaneous frequency
    // is f + t*df/dt, which here would be double the intended sweep rate.
    phaseRad += 2.0 * PI * f * dt;
    if (phaseRad > 2.0 * PI) phaseRad -= 2.0 * PI;
    const float amp = amplitudeMicrosteps((float)f);
    targetMicrosteps.store((int32_t)lroundf(amp * sinf((float)phaseRad)),
                           std::memory_order_relaxed);
    if (((++report) & 0x3FFFU) == 0) {
      Serial.printf("# F,%.3f,%.2f,%ld\n", t, f, (long)actualMicrosteps);
    }
    if (pollAbort()) return false;
  }
  return true;
}

bool dwellSegment(const char *label, double f_hz, double duration_s, bool excite) {
  Serial.printf("# SEG,%s,f=%.2f,dur=%.3f,excite=%u\n", label, f_hz, duration_s,
                excite);
  if (!excite) {
    targetMicrosteps.store(0, std::memory_order_relaxed);
    const uint64_t t0 = esp_timer_get_time();
    while ((esp_timer_get_time() - t0) * 1e-6 < duration_s) {
      if (pollAbort()) return false;
      delay(2);
    }
    return true;
  }
  return sweepSegment(label, f_hz, f_hz, duration_s);
}

bool syncMarker(const char *label) {
  Serial.printf("# SEG,%s,sync\n", label);
  for (int i = 0; i < SYNC_BURSTS; ++i) {
    if (!sweepSegment("sync_burst", SYNC_FREQ_HZ, SYNC_FREQ_HZ, SYNC_BURST_S))
      return false;
    if (i < SYNC_BURSTS - 1) {
      if (!dwellSegment("sync_gap", 0, SYNC_GAP_S, false)) return false;
    }
  }
  return true;
}

double sweepTime(double f0, double f1) {
  return fabs(f1 - f0) / SWEEP_RATE_HZ_PER_S;
}

bool runProfile(bool full) {
  phase = RUNNING;
  abortRequested = false;
  led(0, 60, 0);
  stepTaskRun.store(true, std::memory_order_relaxed);
  Serial.printf("# RUN_START,profile=%s,cruise_microsteps=%.0f,mres=%u,"
                "current_ma=%u\n", full ? "FULL" : "HIGH", CRUISE_MICROSTEPS,
                MRES, CURRENT_RMS_MA);

  bool ok = dwellSegment("pre_roll", 0, PRE_ROLL_S, false) && syncMarker("marker_1");

  if (ok && full) {
    ok = sweepSegment("up_low", SWEEP_LOW_START_HZ, SWEEP_GAP_LOW_HZ,
                      sweepTime(SWEEP_LOW_START_HZ, SWEEP_GAP_LOW_HZ));
    if (ok) Serial.printf("# GAP_JUMP,%.0f->%.0f\n", SWEEP_GAP_LOW_HZ,
                          SWEEP_GAP_HIGH_HZ);
    if (ok) ok = sweepSegment("up_high", SWEEP_GAP_HIGH_HZ, SWEEP_HIGH_END_HZ,
                              sweepTime(SWEEP_GAP_HIGH_HZ, SWEEP_HIGH_END_HZ));
    if (ok) ok = dwellSegment("mid_dwell", SWEEP_HIGH_END_HZ, MID_DWELL_S, false);
    if (ok) ok = sweepSegment("down_high", SWEEP_HIGH_END_HZ, SWEEP_GAP_HIGH_HZ,
                              sweepTime(SWEEP_HIGH_END_HZ, SWEEP_GAP_HIGH_HZ));
    if (ok) Serial.printf("# GAP_JUMP,%.0f->%.0f\n", SWEEP_GAP_HIGH_HZ,
                          SWEEP_GAP_LOW_HZ);
    if (ok) ok = sweepSegment("down_low", SWEEP_GAP_LOW_HZ, SWEEP_LOW_START_HZ,
                              sweepTime(SWEEP_GAP_LOW_HZ, SWEEP_LOW_START_HZ));
  } else if (ok) {
    // HIGH: settle at the start frequency first, so 250 Hz is reached by a
    // steady dwell rather than by a jump. That is the whole point of this
    // profile - it is what separates a real 250 Hz feature from the gap edge.
    ok = dwellSegment("settle_in", SWEEP_GAP_HIGH_HZ, SETTLE_S, true);
    if (ok) ok = sweepSegment("up", SWEEP_GAP_HIGH_HZ, SWEEP_HIGH_END_HZ,
                              sweepTime(SWEEP_GAP_HIGH_HZ, SWEEP_HIGH_END_HZ));
    if (ok) ok = dwellSegment("mid_dwell", SWEEP_HIGH_END_HZ, MID_DWELL_S, true);
    if (ok) ok = sweepSegment("down", SWEEP_HIGH_END_HZ, SWEEP_GAP_HIGH_HZ,
                              sweepTime(SWEEP_HIGH_END_HZ, SWEEP_GAP_HIGH_HZ));
    if (ok) ok = dwellSegment("settle_out", SWEEP_GAP_HIGH_HZ, SETTLE_S, true);
  }

  if (ok) ok = syncMarker("marker_2");
  if (ok) ok = dwellSegment("tail", 0, TAIL_S, false);

  // Return to the origin before releasing the step task, so the carriage is
  // left where it started rather than wherever the waveform stopped.
  targetMicrosteps.store(0, std::memory_order_relaxed);
  delay(200);
  stepTaskRun.store(false, std::memory_order_relaxed);

  phase = ok ? DONE : IDLE;
  led(ok ? 0 : 60, 0, ok ? 60 : 0);
  Serial.printf("# RUN_%s,residual_microsteps=%ld\n",
                ok ? "COMPLETE" : "ABORTED", (long)actualMicrosteps);
  return ok;
}

// --------------------------------------------------------------- commands

void handleCommand(String line) {
  line.trim();
  if (line.length() == 0) return;
  String verb = line, rest = "";
  const int sp = line.indexOf(' ');
  if (sp > 0) { verb = line.substring(0, sp); rest = line.substring(sp + 1); rest.trim(); }
  verb.toUpperCase();
  rest.toUpperCase();

  if (verb == "PING") {
    Serial.printf("# PONG,fw=stage_chirp_sweep,ready=%u\n", driverReady);
  } else if (verb == "CHECK") {
    Serial.printf("# CHECK,ready=%u,mres=%u,current_ma=%u,r_sense=%.3f,"
                  "cruise_microsteps=%.0f,gap=%.0f-%.0f\n",
                  driverReady, MRES, CURRENT_RMS_MA, R_SENSE_OHM,
                  CRUISE_MICROSTEPS, SWEEP_GAP_LOW_HZ, SWEEP_GAP_HIGH_HZ);
  } else if (verb == "ABORT") {
    abortRequested = true;
    Serial.println("# ABORT");
  } else if (verb == "RUN") {
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    if (rest != "FULL" && rest != "HIGH") {
      Serial.println("# ERR,reason=specify_RUN_FULL_or_RUN_HIGH");
      return;
    }
    const bool full = (rest == "FULL");
    if (full) {
      Serial.println("# WARNING,FULL sweeps to 1 kHz; 150-250 Hz is skipped "
                     "because it caused uncontrolled motion on 2026-09-15. "
                     "Stay at the rig with the power within reach.");
    }
    phase = COUNTDOWN;
    led(0, 0, 60);
    for (uint32_t i = COUNTDOWN_S; i > 0; --i) {
      Serial.printf("# COUNTDOWN,%lu\n", (unsigned long)i);
      const uint32_t t0 = millis();
      while (millis() - t0 < 1000) {
        if (pollAbort()) {
          phase = IDLE;
          led(60, 60, 0);
          Serial.println("# ABORTED,before_start");
          return;
        }
        delay(5);
      }
    }
    runProfile(full);
  } else {
    Serial.printf("# ERR,reason=unknown_command,verb=%s\n", verb.c_str());
  }
}

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, HIGH);
  pinMode(LED_PIN, OUTPUT);
  led(60, 60, 0);

  Serial.begin(115200);
  const uint32_t t0 = millis();
  while (!Serial && millis() - t0 < 3000) delay(10);
  delay(200);

  Serial.println("# BOOT,fw=stage_chirp_sweep");
  Serial.printf("# PINS,step=%u,dir=%u,uart_rx=%u,uart_tx=%u\n",
                STEP_PIN, DIR_PIN, UART_RX_PIN, UART_TX_PIN);
  Serial.printf("# GAP,%.0f-%.0f Hz not driven (2026-09-15 incidents)\n",
                SWEEP_GAP_LOW_HZ, SWEEP_GAP_HIGH_HZ);

  driverReady = configureDriver();
  led(driverReady ? 0 : 60, driverReady ? 0 : 0, driverReady ? 60 : 0);

  xTaskCreatePinnedToCore(stepTask, "step", 4096, nullptr, 20, nullptr, 0);

  Serial.printf("# READY,ready=%u\n", driverReady);
  Serial.println("# IDLE,RUN FULL | RUN HIGH | CHECK | ABORT  — nothing moves "
                 "until commanded");
}

void loop() {
  static String buf;
  while (Serial.available()) {
    const char ch = (char)Serial.read();
    if (ch == '\n' || ch == '\r') {
      if (buf.length()) { String l = buf; buf = ""; handleCommand(l); }
    } else if (buf.length() < 120) {
      buf += ch;
    }
  }
  delay(2);
}

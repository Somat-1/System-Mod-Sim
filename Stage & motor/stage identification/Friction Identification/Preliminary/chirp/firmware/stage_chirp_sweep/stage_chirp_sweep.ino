/*
 * Mechanical chirp sweep — ESP32-S3 + TMC2209, two selectable profiles.
 *
 *   RUN FULL   1 Hz -> 1 kHz and back, with a hard 150-250 Hz gap
 *   RUN HIGH   250 Hz -> 1 kHz and back, continuous, started from a dwell
 *   CHECK      configure and report, move nothing
 *   REINIT     retry the driver configuration and report why it fails
 *   WIGGLE     3 x +-100 full steps and back: a visible check that the motor
 *              follows; also runs automatically at the start of every RUN
 *   TEST f a s one fixed frequency f [Hz] at amplitude a [microsteps] for s
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
 * WIRING — the Stepper Isolated bench. STEP GPIO7, DIR GPIO8, EN GPIO9 (ENN,
 * active low, held HIGH until the driver is verified), TMC UART TX 17 /
 * RX 18, LED GPIO48, shaft encoder A/B/Z on GPIO4/5/6 (20000 counts/rev).
 *
 * ---------------------------------------------------------------------------
 * AMPLITUDE — constant +-2.5 full steps (20 microsteps at MRES 8); current
 * set per run (run_chirp.py --current-ma, 1000 mA for run 8, 2026-10-13)
 *
 * History: +-1.5 full steps / 400 mA gave a small accelerometer signal;
 * +-1.875 / 500 mA and +-1.875 / 1000 mA ran clean. +-2.5 is the operator's
 * choice for a better-defined recording and is ABOVE the 2-full-step lag at
 * which a stepper slips: a rotor that falls behind can lose sync and run
 * away. With no encoder fitted the sync-loss abort is off. Stay at the power
 * switch. MRES 8 keeps 2*pi*1000*20 = 126k microsteps/s at 1 kHz, so the
 * amplitude is never clamped.
 *
 * Constant across the whole band. History, 2026-10-13: a constant-acceleration
 * law (+-1 full step to 100 Hz, then 1/f^2, MRES 256) commanded at most +-0.16
 * full step at 250 Hz and +-0.01 at 1 kHz. That is inside this motor's measured
 * friction dead band (0.15-0.3 full steps, microstep resolution test), and a
 * HIGH run with it showed nothing on the accelerometer. +-1.5 full steps:
 *
 *  - is well above the dead band, so the rotor is actually driven.
 *
 * MEASUREMENT — every sweep window (0.1 s, at least one cycle) is demodulated
 * on the board against the sweep's own phase, for both the steps actually
 * issued and the encoder, and reported as
 *
 *   # W,t,f,cycles,f_issued,amp_target,amp_issued,amp_enc,phase_deg,lag_max,dc_err
 *
 * (amplitudes in microsteps at MRES 256; the encoder converted at 2.56
 * microsteps/count). amp_issued vs amp_target and f_issued vs f prove the step
 * generator kept up; amp_enc and phase are the response. A sync check aborts
 * the run if the rotor is ever more than 1.75 full steps from the command, or
 * its window mean drifts by more than half a full step.
 * ---------------------------------------------------------------------------
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <esp_timer.h>
#include <esp_rom_sys.h>
#include <atomic>
#include <cmath>
#include "driver/pulse_cnt.h"
#include "esp_task_wdt.h"

constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;     // TMC ENN: HIGH = outputs off
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr gpio_num_t ENC_A = GPIO_NUM_4;
constexpr gpio_num_t ENC_B = GPIO_NUM_5;
constexpr gpio_num_t ENC_Z = GPIO_NUM_6;
constexpr uint8_t DRIVER_ADDRESS = 0;

/* R_SENSE DIFFERS BETWEEN BOARDS ON THIS RIG — CHECK BEFORE RUNNING.
 * 0.10: the module now fitted is marked R100 (as in the v3 chirp, and in
 * Stepper_Encoder_Slow_Test). The v6 campaign's board used 0.11. Wrong value =
 * wrong coil current by roughly 10%, which at 400 mA is enough to change what
 * the motor does at resonance. Read the markings on the board actually fitted. */
constexpr float R_SENSE_OHM = 0.10F;

constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t MRES = 8;
constexpr uint16_t CURRENT_RMS_MA = 500;      // default; CURRENT changes it
uint16_t currentMa = CURRENT_RMS_MA;

constexpr float CRUISE_FULL_STEPS = 2.5F;
constexpr float CRUISE_MICROSTEPS = CRUISE_FULL_STEPS * MRES;     // 20
// No acceleration roll-off: the amplitude is constant at every frequency.
// (Set a finite corner to bring back A0*(F_CORNER/f)^2 above it.)
constexpr float F_CORNER_HZ = 1.0e9F;
constexpr uint32_t MIN_STEP_PERIOD_US = 5;                        // 200 kHz ceiling
constexpr float ENC_COUNTS_PER_FULLSTEP = 100.0F;                 // 20000 / 200
constexpr float US_PER_COUNT = MRES / ENC_COUNTS_PER_FULLSTEP;    // 2.56
// Above the +-1.875 full-step amplitude, so a rotor that merely fails to
// follow is not taken for a slip; a real slip is a 4-full-step jump, which
// the window-mean (DC) check catches.
constexpr float SYNC_LIMIT_US = 2.5F * MRES;                      // instantaneous
constexpr float DC_LIMIT_US = 0.5F * MRES;                        // window mean
constexpr double WINDOW_S = 0.1;

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

constexpr uint32_t COUNTDOWN_S = 7;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

// Shared between the cores: core 1 writes, core 0 reads. Only ever one writer
// and one reader, and the value is a plain integer, so relaxed ordering on an
// atomic is enough - no lock, which must not exist on the step path.
std::atomic<int32_t> targetMicrosteps{0};
std::atomic<bool> stepTaskRun{false};
volatile bool abortRequested = false;
volatile int32_t actualMicrosteps = 0;     // written by core 0 only
pcnt_unit_handle_t pcnt = nullptr;
double runT0 = 0;                          // seconds, for the W lines
const char *syncFault = nullptr;
bool encoderPresent = false;   // decided from the WIGGLE pre-check each run
int32_t wiggleEncMax = 0;

/** Phase-locked demodulation of one window, for the issued steps and the
 *  encoder, against the sweep's own phase. */
struct Window {
  double T = 0, ic = 0, qc = 0, ie = 0, qe = 0, sc = 0, se = 0, fsum = 0;
  float lagMax = 0;
  uint32_t cycles = 0, crossings = 0;
  int lastSign = 0;
  // lastSign carries over, so a crossing at the window boundary still counts.
  void reset() { const int s = lastSign; *this = Window(); lastSign = s; }
};

// -------------------------------------------------------------- encoder
bool setupEncoder() {
  pcnt_unit_config_t uc = {};
  uc.high_limit = 30000; uc.low_limit = -30000; uc.flags.accum_count = 1;
  if (pcnt_new_unit(&uc, &pcnt) != ESP_OK) return false;
  pcnt_glitch_filter_config_t filter = {};
  filter.max_glitch_ns = 1000;
  pcnt_unit_set_glitch_filter(pcnt, &filter);
  pcnt_chan_config_t ac = {}, bc = {};
  ac.edge_gpio_num = ENC_A; ac.level_gpio_num = ENC_B;
  bc.edge_gpio_num = ENC_B; bc.level_gpio_num = ENC_A;
  pcnt_channel_handle_t cha = nullptr, chb = nullptr;
  if (pcnt_new_channel(pcnt, &ac, &cha) != ESP_OK ||
      pcnt_new_channel(pcnt, &bc, &chb) != ESP_OK) return false;
  pcnt_channel_set_edge_action(cha, PCNT_CHANNEL_EDGE_ACTION_DECREASE,
                               PCNT_CHANNEL_EDGE_ACTION_INCREASE);
  pcnt_channel_set_level_action(cha, PCNT_CHANNEL_LEVEL_ACTION_KEEP,
                                PCNT_CHANNEL_LEVEL_ACTION_INVERSE);
  pcnt_channel_set_edge_action(chb, PCNT_CHANNEL_EDGE_ACTION_INCREASE,
                               PCNT_CHANNEL_EDGE_ACTION_DECREASE);
  pcnt_channel_set_level_action(chb, PCNT_CHANNEL_LEVEL_ACTION_KEEP,
                                PCNT_CHANNEL_LEVEL_ACTION_INVERSE);
  gpio_pullup_en(ENC_A); gpio_pullup_en(ENC_B);
  // Without these watch points accum_count does nothing and the count
  // restarts at 0 every 30000 counts.
  pcnt_unit_add_watch_point(pcnt, 30000);
  pcnt_unit_add_watch_point(pcnt, -30000);
  pcnt_unit_enable(pcnt); pcnt_unit_clear_count(pcnt); pcnt_unit_start(pcnt);
  return true;
}

int32_t encoderCount() {
  int v = 0;
  if (pcnt) pcnt_unit_get_count(pcnt, &v);
  return v;
}
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
  digitalWrite(EN_PIN, HIGH);
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(100);
  const uint32_t ioin = driver.IOIN();
  Serial.printf("# DRIVER_ID,ioin=0x%08lX,version=0x%02X,expected=0x21\n",
                (unsigned long)ioin, (unsigned)((ioin >> 24) & 0xFF));
  driver.begin();
  driver.toff(5);
  driver.blank_time(24);
  // Current comes from IRUN/IHOLD alone, not scaled by the module's VREF
  // potentiometer (I_scale_analog resets to 1 on this chip). Hold current is
  // run current with no power-down: the dwells are the measurements, and
  // halving the current part-way through each one would change what is
  // measured.
  driver.I_scale_analog(false);
  driver.hold_multiplier(1.0F);
  driver.rms_current(currentMa);
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);
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
             stealth = d & (1UL << 30), analog = g & 1UL;
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,drv_status=0x%08lX,"
                "spreadcycle=%u,intpol=%u,stealth=%u,r_sense=%.3f\n",
                (unsigned long)g, (unsigned long)c, (unsigned long)d,
                spread, intpol, stealth, R_SENSE_OHM);
  if (!spread || intpol || stealth || analog) {
    Serial.println("# ERR,reason=driver_mode_not_confirmed");
    return false;
  }
  digitalWrite(EN_PIN, LOW);
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
    // No yield while a sweep runs. A vTaskDelay(1) here used to stop stepping
    // for a whole 1 ms tick every 65536 passes, which the window check showed
    // as a 3-17% amplitude shortfall and up to 11% lost cycles above 600 Hz
    // (2026-10-13). setup() stops the task watchdog watching idle instead;
    // nothing else runs on core 0, and the task still sleeps when idle.
  }
}

// -------------------------------------------------------- the sweep itself

bool rateClamped = false;   // set if the step-rate guard ever had to act
float ampOverride = -1.0F;  // TEST: a fixed amplitude instead of the law

/** Commanded amplitude in microsteps: constant to F_CORNER, then 1/f^2. */
float amplitudeMicrosteps(float f_hz) {
  float amp = CRUISE_MICROSTEPS;
  if (ampOverride >= 0.0F) {
    amp = ampOverride;
  } else if (f_hz > F_CORNER_HZ) {
    const float r = F_CORNER_HZ / f_hz;
    amp *= r * r;
  }
  // A sine of amplitude A at f needs a peak rate of 2*pi*f*A microsteps/s.
  // The law above keeps this under the ceiling; the guard only reports.
  const float maxRate = 1e6F / MIN_STEP_PERIOD_US;
  const float needed = 2.0F * PI * f_hz * amp;
  if (needed > maxRate) { amp = maxRate / (2.0F * PI * f_hz); rateClamped = true; }
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

void emitWindow(Window &w) {
  if (w.T <= 0) return;
  const double f = w.fsum / w.T;
  const double ampC = 2.0 * sqrt(w.ic * w.ic + w.qc * w.qc) / w.T;
  const double ampE = 2.0 * sqrt(w.ie * w.ie + w.qe * w.qe) / w.T;
  double ph = (atan2(w.qe, w.ie) - atan2(w.qc, w.ic)) * 180.0 / PI;
  while (ph > 180.0) ph -= 360.0;
  while (ph <= -180.0) ph += 360.0;
  const double dc = (w.se - w.sc) / w.T;
  const double tNow = esp_timer_get_time() * 1e-6 - runT0;
  Serial.printf("# W,%.3f,%.2f,%lu,%.2f,%.3f,%.3f,%.3f,%.1f,%.0f,%.1f\n", tNow, f,
                (unsigned long)w.cycles, w.crossings / w.T,
                amplitudeMicrosteps((float)f), ampC, ampE, ph, w.lagMax, dc);
  if (encoderPresent && fabs(dc) > DC_LIMIT_US && !syncFault) syncFault = "dc_drift";
}

/** Drive a sine at a frequency that ramps linearly from f0 to f1. */
// Visible pre-check: three slow moves of +-WIGGLE_FULL_STEPS and back to 0,
// so whoever is at the rig can confirm the motor follows the STEP/DIR pulses
// before the sweep, whose amplitudes are too small to see.
constexpr float WIGGLE_FULL_STEPS = 100.0F;    // +-1 mm stage, +-180 deg shaft
constexpr float WIGGLE_SPEED_FS = 100.0F;      // full steps/s
constexpr int WIGGLE_REPS = 3;

bool rampTo(int32_t goal) {
  const int32_t start = targetMicrosteps.load(std::memory_order_relaxed);
  const double rate = WIGGLE_SPEED_FS * MRES;           // microsteps/s
  const double dur = fabs((double)(goal - start)) / rate;
  const uint64_t t0 = esp_timer_get_time();
  for (;;) {
    const double t = (esp_timer_get_time() - t0) * 1e-6;
    const double x = t >= dur ? 1.0 : t / dur;
    targetMicrosteps.store(start + (int32_t)lround((goal - start) * x),
                           std::memory_order_relaxed);
    const int32_t e = abs(encoderCount());
    if (e > wiggleEncMax) wiggleEncMax = e;
    if (pollAbort()) return false;
    if (x >= 1.0 && actualMicrosteps == goal) return true;
    delayMicroseconds(200);
  }
}

bool wiggle() {
  const int32_t n = (int32_t)(WIGGLE_FULL_STEPS * MRES);
  Serial.printf("# WIGGLE_START,full_steps=%.0f,speed_fs=%.0f,reps=%d\n",
                WIGGLE_FULL_STEPS, WIGGLE_SPEED_FS, WIGGLE_REPS);
  const bool wasRunning = stepTaskRun.load(std::memory_order_relaxed);
  wiggleEncMax = 0;
  stepTaskRun.store(true, std::memory_order_relaxed);
  bool ok = true;
  for (int r = 0; r < WIGGLE_REPS && ok; ++r) {
    ok = rampTo(n) && (delay(300), true) && rampTo(-n) && (delay(300), true);
    Serial.printf("# WIGGLE,rep=%d,issued=%ld\n", r + 1, (long)actualMicrosteps);
  }
  if (ok) ok = rampTo(0);
  if (!wasRunning) stepTaskRun.store(false, std::memory_order_relaxed);
  // +-100 full steps is ~10000 counts on the 20000 CPR shaft encoder; a
  // tenth of that is unmistakably an encoder, nothing at all means none.
  encoderPresent = wiggleEncMax > (int32_t)(0.1F * WIGGLE_FULL_STEPS * ENC_COUNTS_PER_FULLSTEP);
  Serial.printf("# WIGGLE_DONE,ok=%u,issued=%ld,enc_max=%ld,encoder=%s\n", ok,
                (long)actualMicrosteps, (long)wiggleEncMax,
                encoderPresent ? "present" : "absent_sync_check_off");
  return ok;
}

bool sweepSegment(const char *label, double f0, double f1, double duration_s) {
  Serial.printf("# SEG,%s,f0=%.2f,f1=%.2f,dur=%.3f\n", label, f0, f1, duration_s);
  const uint64_t t0 = esp_timer_get_time();
  double phaseRad = 0.0;
  uint64_t lastUs = t0;
  Window w;
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
    if (phaseRad > 2.0 * PI) { phaseRad -= 2.0 * PI; w.cycles++; }
    const float amp = amplitudeMicrosteps((float)f);
    const int32_t target = (int32_t)lroundf(amp * sinf((float)phaseRad));
    targetMicrosteps.store(target, std::memory_order_relaxed);

    // Measurement, against the same phase that generated the command.
    const int32_t issued = actualMicrosteps;
    const float enc = encoderCount() * US_PER_COUNT;
    const double s = sin(phaseRad), c = cos(phaseRad);
    w.T += dt; w.fsum += f * dt;
    w.ic += issued * s * dt; w.qc += issued * c * dt;
    w.ie += enc * s * dt;    w.qe += enc * c * dt;
    w.sc += issued * dt;     w.se += enc * dt;
    const float lag = fabsf((float)(target - issued));
    if (lag > w.lagMax) w.lagMax = lag;
    const int sign = issued > 0 ? 1 : (issued < 0 ? -1 : 0);
    if (sign) {
      if (sign > 0 && w.lastSign < 0) w.crossings++;
      w.lastSign = sign;
    }
    if (encoderPresent && fabsf(enc - issued) > SYNC_LIMIT_US && !syncFault)
      syncFault = "rotor_lag";
    if (w.T >= WINDOW_S && w.cycles >= 1) { emitWindow(w); w.reset(); }
    if (syncFault) {
      Serial.printf("# ERR,reason=sync_lost,check=%s,issued=%ld,enc_us=%.0f\n",
                    syncFault, (long)issued, enc);
      return false;
    }
    if (pollAbort()) return false;
  }
  emitWindow(w);
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
    // Markers follow the same amplitude law as the sweep.
    if (!sweepSegment("sync_burst", SYNC_FREQ_HZ, SYNC_FREQ_HZ, SYNC_BURST_S))
      return false;
    if (i < SYNC_BURSTS - 1) {
      if (!dwellSegment("sync_gap", 0, SYNC_GAP_S, false)) return false;
    }
  }
  return true;
}

/** The turn-around at 1 kHz, made unmistakable in an external recording:
 *  1 s silence, the 3-burst 750 Hz sync marker, 1 s silence, then the down
 *  sweep. */
bool topMarker() {
  Serial.printf("# TOP_1KHZ,t=%.3f,up_sweep_done\n",
                esp_timer_get_time() * 1e-6 - runT0);
  return dwellSegment("top_gap_1", 0, 1.0, false) && syncMarker("top_marker") &&
         dwellSegment("top_gap_2", 0, 1.0, false);
}

double sweepTime(double f0, double f1) {
  return fabs(f1 - f0) / SWEEP_RATE_HZ_PER_S;
}

bool runProfile(bool full) {
  phase = RUNNING;
  abortRequested = false;
  syncFault = nullptr;
  rateClamped = false;
  if (pcnt) pcnt_unit_clear_count(pcnt);       // encoder and command share 0
  runT0 = esp_timer_get_time() * 1e-6;
  led(0, 60, 0);
  stepTaskRun.store(true, std::memory_order_relaxed);
  Serial.printf("# RUN_START,profile=%s,a0_microsteps=%.0f,f_corner=%.0f,mres=%u,"
                "current_ma=%u,us_per_count=%.2f\n", full ? "FULL" : "HIGH",
                CRUISE_MICROSTEPS, F_CORNER_HZ, MRES, currentMa, US_PER_COUNT);

  if (pcnt) pcnt_unit_clear_count(pcnt);
  bool ok = wiggle();
  if (pcnt) pcnt_unit_clear_count(pcnt);       // sweep starts from encoder 0
  ok = ok && dwellSegment("pre_roll", 0, PRE_ROLL_S, false) &&
            syncMarker("marker_1");

  if (ok && full) {
    ok = sweepSegment("up_low", SWEEP_LOW_START_HZ, SWEEP_GAP_LOW_HZ,
                      sweepTime(SWEEP_LOW_START_HZ, SWEEP_GAP_LOW_HZ));
    if (ok) Serial.printf("# GAP_JUMP,%.0f->%.0f\n", SWEEP_GAP_LOW_HZ,
                          SWEEP_GAP_HIGH_HZ);
    if (ok) ok = sweepSegment("up_high", SWEEP_GAP_HIGH_HZ, SWEEP_HIGH_END_HZ,
                              sweepTime(SWEEP_GAP_HIGH_HZ, SWEEP_HIGH_END_HZ));
    if (ok) ok = topMarker();
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
    if (ok) ok = topMarker();
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
  // Driver temperature flags after the run: OTPW (warning) bit 0, OT bit 1.
  uint32_t ds = 0;
  if (readTwice([] { return driver.DRV_STATUS(); }, ds)) {
    Serial.printf("# DRV_STATUS,0x%08lX,otpw=%u,ot=%u,current_ma=%u\n",
                  (unsigned long)ds, (unsigned)(ds & 1UL), (unsigned)((ds >> 1) & 1UL),
                  currentMa);
  }
  Serial.printf("# RUN_%s,residual_microsteps=%ld,enc=%ld,rate_clamped=%u\n",
                ok ? "COMPLETE" : (syncFault ? "FAILED" : "ABORTED"),
                (long)actualMicrosteps, (long)encoderCount(), rateClamped);
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
                  "a0_microsteps=%.0f,f_corner=%.0f,amp_250=%.1f,amp_1000=%.2f,"
                  "gap=%.0f-%.0f,enc=%ld\n",
                  driverReady, MRES, currentMa, R_SENSE_OHM,
                  CRUISE_MICROSTEPS, F_CORNER_HZ, amplitudeMicrosteps(250),
                  amplitudeMicrosteps(1000), SWEEP_GAP_LOW_HZ, SWEEP_GAP_HIGH_HZ,
                  (long)encoderCount());
  } else if (verb == "CURRENT") {
    // CURRENT <mA>: coil current, RMS, 100-1200 (TMC2209 limit with cooling).
    const long v = rest.toInt();
    if (v < 100 || v > 1200) { Serial.println("# ERR,reason=bad_current"); return; }
    currentMa = (uint16_t)v;
    driver.rms_current(currentMa);
    Serial.printf("# CURRENT_OK,current_ma=%u,irun=%u,ihold=%u\n", currentMa,
                  driver.irun(), driver.ihold());
  } else if (verb == "WIGGLE") {
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    abortRequested = false;
    led(0, 60, 0);
    wiggle();
    led(0, 0, 60);
  } else if (verb == "TEST") {
    // TEST <f_hz> <amplitude_microsteps> <seconds>: one fixed frequency at a
    // chosen amplitude, for bring-up. Bounded: outside the 150-250 Hz gap,
    // at most +-1.5 full steps (below the 2-full-step slip limit), under the
    // step-rate ceiling with margin, at most 10 s.
    float f = 0, a = 0, secs = 0;
    if (!driverReady || sscanf(rest.c_str(), "%f %f %f", &f, &a, &secs) != 3 ||
        f < 1 || f > 1000 || (f > SWEEP_GAP_LOW_HZ && f < SWEEP_GAP_HIGH_HZ) ||
        a <= 0 || a > 1.5F * MRES || 2.0F * PI * f * a > 0.8F * 1e6F / MIN_STEP_PERIOD_US ||
        secs <= 0 || secs > 10) {
      Serial.println("# ERR,reason=bad_test_args_or_not_ready");
      return;
    }
    abortRequested = false; syncFault = nullptr;
    if (pcnt) pcnt_unit_clear_count(pcnt);
    runT0 = esp_timer_get_time() * 1e-6;
    ampOverride = a;   // sync checks follow the last WIGGLE's encoder finding
    stepTaskRun.store(true, std::memory_order_relaxed);
    led(0, 60, 0);
    const bool ok = sweepSegment("test", f, f, secs);
    targetMicrosteps.store(0, std::memory_order_relaxed);
    delay(200);
    stepTaskRun.store(false, std::memory_order_relaxed);
    ampOverride = -1.0F;
    led(0, 0, 60);
    Serial.printf("# TEST_DONE,ok=%u,residual_microsteps=%ld\n", ok,
                  (long)actualMicrosteps);
  } else if (verb == "REINIT") {
    // Retry the driver setup without a reset (boot messages are easily missed).
    driverReady = configureDriver() && pcnt != nullptr;
    Serial.printf("# REINIT_DONE,ready=%u,encoder=%u\n", driverReady, pcnt != nullptr);
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
    // A new RUN starts clean: an ABORT from an earlier run would otherwise
    // still be set and cancel this countdown at its first check.
    abortRequested = false;
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
  pinMode(EN_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);
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
  Serial.printf("# PINS,step=%u,dir=%u,en=%u,uart_rx=%u,uart_tx=%u\n",
                STEP_PIN, DIR_PIN, EN_PIN, UART_RX_PIN, UART_TX_PIN);
  Serial.printf("# GAP,%.0f-%.0f Hz not driven (2026-09-15 incidents)\n",
                SWEEP_GAP_LOW_HZ, SWEEP_GAP_HIGH_HZ);

  const bool encOk = setupEncoder();
  if (!encOk) Serial.println("# ERR,reason=pcnt_init_failed");
  driverReady = configureDriver() && encOk;
  led(driverReady ? 0 : 60, driverReady ? 0 : 0, driverReady ? 60 : 0);

  // The step task never yields during a sweep, so core 0's idle task cannot
  // run then. Stop the task watchdog watching the idle tasks (nothing else
  // needs core 0: no WiFi). Arduino's disableCore0WDT() is not used: on this
  // core version it left something calling esp_task_wdt_reset() on the
  // removed entry, flooding USB with "task not found" errors.
  esp_task_wdt_config_t wdt = {};
  wdt.timeout_ms = 5000;
  wdt.idle_core_mask = 0;
  wdt.trigger_panic = false;
  const esp_err_t wdtErr = esp_task_wdt_reconfigure(&wdt);
  xTaskCreatePinnedToCore(stepTask, "step", 4096, nullptr, 20, nullptr, 0);

  Serial.printf("# READY,ready=%u,wdt_reconfigure=%d\n", driverReady, (int)wdtErr);
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

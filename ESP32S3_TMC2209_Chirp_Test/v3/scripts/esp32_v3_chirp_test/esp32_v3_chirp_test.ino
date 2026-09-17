/*
 * ESP32-S3 + TMC2209 v3 resonance-band chirp: hard frequency gap, no notch.
 * Dual-core step generation.
 *
 * Forked from
 * ESP32S3_TMC2209_Chirp_Test/v2/scripts/esp32_v2_chirp_test_rsense_fixed
 * (itself v2's originally-validated sequence with R_SENSE corrected to
 * 0.11). One structural change from that baseline:
 *
 *   The sweep now has a genuine GAP from 150-250 Hz -- zero commanded
 *   excitation in that band, not reduced amplitude. An earlier v3 attempt
 *   (see git history) tried a smooth notch there instead: cruise amplitude
 *   tapered down to as little as 3% at the center frequency, widened
 *   afterward to cover 65-285 Hz at that same 97% suppression. Both
 *   versions of that notch still caused uncontrolled motor motion,
 *   reported live around 200 Hz -- reducing the amplitude in that band was
 *   not sufficient, however aggressively. So this version does not reduce
 *   amplitude there; it does not visit that band's frequencies at all.
 *   UP_SWEEP runs 1->150 Hz, jumps straight to 250 Hz, then continues
 *   250->1000 Hz; DOWN_SWEEP mirrors that. Both legs share one sweep rate
 *   (7.8125 Hz/s) so the chirp's character is consistent across the jump.
 *   The jump itself is an instantaneous target change, same as any other
 *   segment boundary in this file (e.g. PRE_ROLL -> MARKER_1) -- core 0's
 *   step task just steps toward whatever target it's given, as fast as
 *   MIN_STEP_PERIOD_US allows, so a jump bounded by the cruise amplitude
 *   (tens of microsteps) is not a new class of event.
 *
 *   The old smooth notch (F_N_HZ_PLACEHOLDER, NOTCH_HALFWIDTH_HZ,
 *   NOTCH_DEPTH_RATIO, resonanceNotchTaper()) is removed rather than kept
 *   alongside the gap -- it was shown not to matter, and removing it keeps
 *   amplitudeMicrosteps() to just the two clamps that are still load-bearing
 *   (stroke, step-rate).
 *
 * Everything else carried over from the R_SENSE-fixed baseline is
 * unchanged: R_SENSE_OHM=0.11, waits indefinitely for the TMC2209 to be
 * detected (retrying once a second) rather than failing after one attempt,
 * a 7 s countdown once detected before auto-starting (ABORT cancels it),
 * and an onboard LED phase indicator (yellow=waiting/configuring,
 * blue=detected/counting down, green=running, cyan=done, red=failed or
 * aborted) since serial has proven unreliable to watch live on this rig.
 *
 * Everything below this point is v2's own original header, unchanged:
 *
 * ---
 *
 * Validated on real hardware: a full RUN (pre-roll, both sync markers, both
 * sweep directions, mid-dwell, tail) completed segment-by-segment on
 * schedule with no watchdog reset and the driver connected/configured
 * correctly over the real UART link. This is a corrected/reworked version
 * of the original single-core, wrong-pin
 * ESP32S3_TMC2209_Chirp_Test/v2/scripts/esp32_v2_chirp_test/esp32_v2_chirp_test.ino,
 * built after benchmarking three step-generation approaches on this exact
 * board (scratch-pin tests, no TMC2209 involved, in a separate benchmark
 * sketch -- not committed to the repo):
 *
 *   1. Polled single-core loop (target computed and edge issued in the
 *      same loop iteration): real edge rate flatlined around ~55 kHz
 *      regardless of demand, because sin()/notch/clamp math ran inline
 *      with GPIO issuance.
 *   2. Single-core hardware-timer ISR (edge issuance decoupled from the
 *      math, but sharing one core): at a firing rate high enough to matter,
 *      the ISR itself starved the main loop's target computation instead
 *      (main-loop rate collapsed 5x), trading spatial staircasing for
 *      temporal staircasing.
 *   3. Dual-core (this file): core 0 runs a dedicated step-issuing task
 *      (integer compare + GPIO register write only, no floating point
 *      anywhere near it); core 1 runs the normal Arduino setup()/loop(),
 *      doing all the sin()/clamp math and publishing a target position.
 *      Neither starves the other.
 *
 * Two corrections from the original file, from real-hardware checks on
 * this board:
 *   - STEP/DIR/UART pins corrected to this rig's actual wiring (confirmed
 *     against the working v4 MRES trajectory campaign sketch): STEP=5,
 *     DIR=6, UART_RX=18, UART_TX=17. The original file's STEP=6/DIR=7/EN=5
 *     were placeholders reused from v1 and did not match this rig.
 *   - EN is physically grounded on this rig (always enabled), not wired to
 *     a GPIO -- so there is no EN_PIN here at all. The driver is enabled
 *     for as long as it is powered.
 *
 * Command generation is a position-tracking control loop:
 * commandedStateAt(elapsed_s) analytically evaluates the exact ideal
 * microstep position for ANY elapsed time (closed-form phase integration,
 * no accumulated numerical drift). Core 1 evaluates this and publishes the
 * result; core 0 continuously steps the real position toward whatever the
 * shared target currently is, rate-limited by MIN_STEP_PERIOD_US.
 *
 * Framework: Arduino-ESP32 3.x; library: TMCStepper.
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

namespace chirp_v3 {

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
constexpr float R_SENSE_OHM = 0.10F;  // R100 marking on this (new) board's
                                       // sense resistors, confirmed ~0.2 ohm
                                       // measured (includes lead resistance).
                                       // Was 0.11 for the old, now-dead board.
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;

// --- Mirrors chirp_v2_schedule.py -- keep these in sync by hand ----------

constexpr uint16_t MRES = 16;               // 1/16 step. See README Section 2.
constexpr uint8_t MRES_REGISTER_CODE = 4;   // TMC2209 MRES field code for 1/16.
constexpr uint16_t CURRENT_RMS_MA = 400;

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

// --- Sweep: two linear legs per direction, a hard gap between them ------
// See file header for why the gap replaced a smooth notch.
constexpr float SWEEP_LOW_START_HZ = 1.0F;
constexpr float SWEEP_GAP_LOW_HZ = 150.0F;    // Top of the low leg.
constexpr float SWEEP_GAP_HIGH_HZ = 250.0F;   // Bottom of the high leg.
constexpr float SWEEP_HIGH_END_HZ = 1000.0F;
constexpr double SWEEP_RATE_HZ_PER_S = 7.8125;  // (1000-250)/96, both legs share it.
constexpr double LOW_LEG_DURATION_S =
    (SWEEP_GAP_LOW_HZ - SWEEP_LOW_START_HZ) / SWEEP_RATE_HZ_PER_S;    // ~19.07 s
constexpr double HIGH_LEG_DURATION_S =
    (SWEEP_HIGH_END_HZ - SWEEP_GAP_HIGH_HZ) / SWEEP_RATE_HZ_PER_S;    // 96 s
constexpr double SWEEP_DURATION_S = LOW_LEG_DURATION_S + HIGH_LEG_DURATION_S;

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
    + SYNC_MARKER_DURATION_S + SWEEP_DURATION_S + TAIL_S;  // ~296 s

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

// Real chirp segments include long silent stretches (30 s pre-roll/tail,
// 5 s mid-dwell) where stepTaskRunning is true but target never changes.
// A tight loop with zero yield for 30 s starves core 0's IDLE task long
// enough to trip the ESP-IDF task watchdog and reboot the chip (confirmed
// on real hardware). disableCore0WDT() alone does not fully fix this on
// this IDF version (see setup()) -- the actual fix is esp_task_wdt_deinit()
// there.
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
// One linear-chirp formula, reused for every leg in both directions --
// passing (endHz, startHz, duration) instead of (startHz, endHz, duration)
// sweeps downward, so the same closed form covers all four legs.
double linearChirpPhaseRad(double t, double startHz, double endHz,
                           double durationS) {
  return 2.0 * PI * (startHz * t + (endHz - startHz) / (2.0 * durationS) * t * t);
}

// --- Amplitude: cruise, clamped by stroke and step-rate only ------------
// No resonance notch here -- see file header for why a hard frequency gap
// replaced it instead.

float amplitudeMicrosteps(float frequencyHz) {
  float amplitude = CRUISE_MICROSTEPS;
  const float omega = 2.0F * PI * frequencyHz;
  const float aStepRate = STEP_RATE_CLAMP_HZ / omega;
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

// One sweep direction, two legs with a hard gap between them. ascending
// picks 1->150 then 250->1000 (UP_SWEEP); !ascending picks 1000->250 then
// 150->1 (DOWN_SWEEP), i.e. the same two legs run high-to-low instead.
CommandState sweepStateAt(double t, bool ascending) {
  const float lowFrom = ascending ? SWEEP_LOW_START_HZ : SWEEP_GAP_LOW_HZ;
  const float lowTo = ascending ? SWEEP_GAP_LOW_HZ : SWEEP_LOW_START_HZ;
  const float highFrom = ascending ? SWEEP_GAP_HIGH_HZ : SWEEP_HIGH_END_HZ;
  const float highTo = ascending ? SWEEP_HIGH_END_HZ : SWEEP_GAP_HIGH_HZ;
  const char *segLow = ascending ? "UP_SWEEP_LOW" : "DOWN_SWEEP_LOW";
  const char *segHigh = ascending ? "UP_SWEEP_HIGH" : "DOWN_SWEEP_HIGH";

  // Ascending sweeps the low leg (1-150) first; descending sweeps the high
  // leg (1000-250) first -- each direction covers its nearer leg first.
  const bool lowLegFirst = ascending;
  const double firstLegS = lowLegFirst ? LOW_LEG_DURATION_S : HIGH_LEG_DURATION_S;

  if (t < firstLegS) {
    const float from = lowLegFirst ? lowFrom : highFrom;
    const float to = lowLegFirst ? lowTo : highTo;
    const float f = from + (to - from) * static_cast<float>(t / firstLegS);
    const double phase = linearChirpPhaseRad(t, from, to, firstLegS);
    const float amplitude = amplitudeMicrosteps(f);
    return {lround(amplitude * sin(phase)), f, lowLegFirst ? segLow : segHigh};
  }
  const double t2 = t - firstLegS;
  const double secondLegS = lowLegFirst ? HIGH_LEG_DURATION_S : LOW_LEG_DURATION_S;
  const float from = lowLegFirst ? highFrom : lowFrom;
  const float to = lowLegFirst ? highTo : lowTo;
  const float f = from + (to - from) * static_cast<float>(t2 / secondLegS);
  const double phase = linearChirpPhaseRad(t2, from, to, secondLegS);
  const float amplitude = amplitudeMicrosteps(f);
  return {lround(amplitude * sin(phase)), f, lowLegFirst ? segHigh : segLow};
}

CommandState commandedStateAt(double elapsedS) {
  double t = elapsedS;

  if (t < PRE_ROLL_S) return {0, NAN, "PRE_ROLL"};
  t -= PRE_ROLL_S;

  if (t < SYNC_MARKER_DURATION_S) {
    return {markerMicrostepsAt(t), NAN, "MARKER_1"};
  }
  t -= SYNC_MARKER_DURATION_S;

  if (t < SWEEP_DURATION_S) return sweepStateAt(t, true);
  t -= SWEEP_DURATION_S;

  if (t < MID_DWELL_S) return {0, NAN, "MID_DWELL"};
  t -= MID_DWELL_S;

  if (t < SYNC_MARKER_DURATION_S) {
    return {markerMicrostepsAt(t), NAN, "MARKER_2"};
  }
  t -= SYNC_MARKER_DURATION_S;

  if (t < SWEEP_DURATION_S) return sweepStateAt(t, false);
  t -= SWEEP_DURATION_S;

  if (t < TAIL_S) return {0, NAN, "TAIL"};
  return {0, NAN, "DONE"};
}

// --- Driver setup ------------------------------------------------------

void dumpRegisters(const char *when) {
  Serial.printf("# REGDUMP %s: GCONF=0x%08lX CHOPCONF=0x%08lX "
                "IHOLD_IRUN=0x%08lX TPWMTHRS=0x%08lX PWMCONF=0x%08lX "
                "TPOWERDOWN=0x%08lX\n",
                when,
                static_cast<unsigned long>(driver.GCONF()),
                static_cast<unsigned long>(driver.CHOPCONF()),
                static_cast<unsigned long>(driver.IHOLD_IRUN()),
                static_cast<unsigned long>(driver.TPWMTHRS()),
                static_cast<unsigned long>(driver.PWMCONF()),
                static_cast<unsigned long>(driver.TPOWERDOWN()));
  // DRV_STATUS, GSTAT, and IOIN are genuine read registers -- the chip
  // computes/latches these itself, unlike IHOLD_IRUN() above (confirmed in
  // TMCStepper's own source to be a local write-shadow return, never a UART
  // read; it can only ever echo what this firmware last sent, not what the
  // chip actually did with it). CS_ACTUAL is what the driver claims to be
  // outputting right now. GSTAT's drv_err specifically: on the TMC2209, a
  // short-circuit event latches this and disables the output stage entirely
  // until GSTAT is cleared or the chip resets -- fully consistent with
  // "everything else checks out, zero current anyway" if it's set. uv_cp
  // (charge pump undervoltage) would mean the high-side gate drive never
  // came up correctly, which would also fully explain it.
  const uint32_t drvStatus = driver.DRV_STATUS();
  Serial.printf("# DRV_STATUS %s: raw=0x%08lX stst=%u stealth=%u "
                "cs_actual=%u olA=%u olB=%u s2gA=%u s2gB=%u s2vsA=%u "
                "s2vsB=%u ot=%u otpw=%u\n",
                when, static_cast<unsigned long>(drvStatus),
                static_cast<unsigned>((drvStatus >> 31) & 1U),
                static_cast<unsigned>((drvStatus >> 30) & 1U),
                static_cast<unsigned>((drvStatus >> 16) & 0x1FU),
                static_cast<unsigned>((drvStatus >> 6) & 1U),
                static_cast<unsigned>((drvStatus >> 7) & 1U),
                static_cast<unsigned>((drvStatus >> 2) & 1U),
                static_cast<unsigned>((drvStatus >> 3) & 1U),
                static_cast<unsigned>((drvStatus >> 4) & 1U),
                static_cast<unsigned>((drvStatus >> 5) & 1U),
                static_cast<unsigned>((drvStatus >> 1) & 1U),
                static_cast<unsigned>(drvStatus & 1U));

  const uint8_t gstat = driver.GSTAT();
  Serial.printf("# GSTAT %s: raw=0x%02X reset=%u drv_err=%u uv_cp=%u\n",
                when, gstat,
                static_cast<unsigned>(gstat & 1U),
                static_cast<unsigned>((gstat >> 1) & 1U),
                static_cast<unsigned>((gstat >> 2) & 1U));

  const uint32_t ioin = driver.IOIN();
  Serial.printf("# IOIN %s: raw=0x%08lX enn_pin=%u step_pin=%u dir_pin=%u "
                "version=0x%02X\n",
                when, static_cast<unsigned long>(ioin),
                static_cast<unsigned>(ioin & 1U),
                static_cast<unsigned>((ioin >> 7) & 1U),
                static_cast<unsigned>((ioin >> 9) & 1U),
                static_cast<unsigned>((ioin >> 24) & 0xFFU));
}

bool configureDriver() {
  DriverSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();

  // GSTAT's drv_err/uv_cp are latching: once set (by a short-circuit event,
  // or a charge-pump undervoltage transient at power-up), they disable the
  // output stage until explicitly cleared, and stay latched across any
  // number of reconfigurations that don't clear them -- which is every
  // configureDriver() call before this one. Clear before anything else, and
  // log what was actually latched (informational: TMCStepper's GSTAT(x)
  // always writes 0b111 regardless of x, clearing all three bits).
  const uint8_t gstatBefore = driver.GSTAT();
  Serial.printf("# GSTAT_BEFORE_CLEAR,raw=0x%02X,reset=%u,drv_err=%u,uv_cp=%u\n",
                gstatBefore, static_cast<unsigned>(gstatBefore & 1U),
                static_cast<unsigned>((gstatBefore >> 1) & 1U),
                static_cast<unsigned>((gstatBefore >> 2) & 1U));
  driver.GSTAT(0);
  const uint8_t gstatAfter = driver.GSTAT();
  Serial.printf("# GSTAT_AFTER_CLEAR,raw=0x%02X,reset=%u,drv_err=%u,uv_cp=%u\n",
                gstatAfter, static_cast<unsigned>(gstatAfter & 1U),
                static_cast<unsigned>((gstatAfter >> 1) & 1U),
                static_cast<unsigned>((gstatAfter >> 2) & 1U));

  // Is uv_cp permanently stuck at 1, or does it flicker? A hard stuck-at-1
  // points at the charge pump being completely non-functional; flickering
  // would suggest something marginal (borderline cap value, noisy VM) --
  // clear-then-immediately-repoll, 20 times over ~1 s, and count each.
  uint8_t uvCpSetCount = 0;
  for (uint8_t i = 0; i < 20; ++i) {
    driver.GSTAT(0);
    delay(10);
    const uint8_t g = driver.GSTAT();
    if ((g >> 2) & 1U) ++uvCpSetCount;
  }
  Serial.printf("# UV_CP_POLL,set_in=%u_of_20_samples_over_~200ms\n",
                uvCpSetCount);

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
  Serial.println("# v3 chirp plan: hard 150-250 Hz gap, no amplitude notch");
  Serial.println("# dual-core step generation: core 0 issues edges, core 1 "
                  "computes the target");
  Serial.printf("# MRES=%u (code %u) current_rms_mA=%u (R_SENSE_OHM=%.2f)\n",
                MRES, MRES_REGISTER_CODE, CURRENT_RMS_MA,
                static_cast<double>(R_SENSE_OHM));
  Serial.printf("# cruise=%.2f full steps (%.1f microsteps)\n",
                static_cast<double>(CRUISE_FULL_STEPS),
                static_cast<double>(CRUISE_MICROSTEPS));
  Serial.printf("# sweep: %.0f->%.0f Hz (%.2f s), GAP %.0f-%.0f Hz (zero "
                "excitation), %.0f->%.0f Hz (%.1f s); mirrored on the way "
                "down; rate %.4f Hz/s both legs\n",
                static_cast<double>(SWEEP_LOW_START_HZ),
                static_cast<double>(SWEEP_GAP_LOW_HZ), LOW_LEG_DURATION_S,
                static_cast<double>(SWEEP_GAP_LOW_HZ),
                static_cast<double>(SWEEP_GAP_HIGH_HZ),
                static_cast<double>(SWEEP_GAP_HIGH_HZ),
                static_cast<double>(SWEEP_HIGH_END_HZ), HIGH_LEG_DURATION_S,
                SWEEP_RATE_HZ_PER_S);
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

}  // namespace chirp_v3

using namespace chirp_v3;

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
  Serial.println("# ESP32-S3/TMC2209 v3 chirp (hard 150-250 Hz gap, R_SENSE "
                  "0.11), dual-core step generation");

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

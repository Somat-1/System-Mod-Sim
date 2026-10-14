/*
 * ESP32-S3 / TMC2209 interactive jog firmware.
 *
 * Takes position commands over USB CDC and moves the stage to them. Its reason
 * to exist is the host script jog_stepper.py: set two positions, then shuttle
 * between them until told to stop.
 *
 * Wiring is the Stepper Isolated bench's, as used by
 * Stepper_Encoder_Slow_Test and TMC2209_Communication_Check:
 *
 *   STEP  GPIO7
 *   DIR   GPIO8
 *   EN    GPIO9  (TMC ENN, active low; held HIGH until the driver is
 *                 configured and verified, then LOW so the motor holds)
 *   UART  ESP_RX GPIO18 / ESP_TX GPIO17  (to TMC2209 PDN_UART via 1k)
 *   ENC   A GPIO4, B GPIO5, Z GPIO6
 *   LED   GPIO48 (onboard RGB)
 *
 * Two deliberate differences from the campaign sketches:
 *
 * 1. Nothing moves at boot. The campaign sketches start their trajectory
 *    automatically, which is only safe when the stage is known to have
 *    clearance. This one waits for a command, always.
 * 2. Step generation is non-blocking. The campaign sketches busy-wait inside
 *    moveMicrosteps() and can only poll for an abort every 32 pulses. Here one
 *    loop() pass emits at most one pulse, so STOP is acted on within a pulse
 *    period rather than at the end of a move.
 *
 * Soft limits are on by default and deliberately tight. The stage's measured
 * travel is 44 mm (tmp_v4_manual_range_sweep) but this firmware has no homing
 * and no end stops, so position 0 means only "wherever the stage was when it
 * was last zeroed". LIMIT exists so a mistyped target cannot drive the stage
 * into its own end stop at full speed.
 *
 * SpreadCycle is forced and MicroPlyer interpolation disabled, both verified by
 * readback, for the reasons in TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md: on
 * this chip both defaults are the wrong choice, so leaving the registers alone
 * actively selects StealthChop with interpolation on.
 *
 * Protocol: newline-terminated ASCII in, '#'-prefixed lines out. Positions are
 * in whole full steps (100 full steps = 1 mm on this rig). Rates are a
 * full-step rate in milliHz, matching the campaign sketches' convention.
 *
 *   PING                -> # PONG,fw=...
 *   STATUS              -> # STATUS,...
 *   POS?                -> # POS,...
 *   ZERO                 call the present position 0
 *   LIMIT <fullsteps>    soft limit, +/- about zero; 0 disables
 *   RATE <milliHz>       full-step rate
 *   MRES <1..256>        microstep resolution, idle only
 *   CURRENT <mA>         RMS coil current
 *   MOVETO <fullsteps>   single move
 *   STEP <n>             relative move of n microsteps at the present MRES
 *   JOG <a> <b>          remember two positions, in full steps
 *   JOGU <a> <b>         the same in 1/256 full-step units, so that a point
 *                        marked between full steps is kept exactly; both must
 *                        lie on the present microstep grid
 *   START                shuttle between them until STOP
 *   STOP                 stop at once, keep position
 *   REINIT               retry the driver configuration without a reset
 *   HOLD <0..100>        standstill current as % of run current (default 100)
 *   SEQ [s] [dwell_ms] [jump_fullsteps] [m1,m2,...]
 *                        the MRES sweep: for each MRES in turn, a marker jump
 *                        out and back, then the CALIB +1/-1 microstep
 *                        oscillation for s seconds. Defaults 30 s, 2000 ms,
 *                        4 full steps, MRES 1,2,4,8,16,32,64.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/pulse_cnt.h"

constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;     // TMC ENN: HIGH = outputs off
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
// Encoder, as on the Stepper Isolated bench.
constexpr gpio_num_t ENC_A = GPIO_NUM_4;
constexpr gpio_num_t ENC_B = GPIO_NUM_5;
constexpr gpio_num_t ENC_Z = GPIO_NUM_6;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;   // R100 sense resistors on this module.
constexpr uint32_t DRIVER_BAUD = 115200;

constexpr uint16_t FULL_STEPS_PER_MM = 100;
constexpr uint16_t DEFAULT_CURRENT_MA = 360;
constexpr uint16_t DEFAULT_MRES = 16;
constexpr uint32_t DEFAULT_RATE_MILLIHZ = 70000;   // 70 full steps/s = 0.7 mm/s
constexpr int32_t DEFAULT_LIMIT_FULLSTEPS = 0;     // off; set with LIMIT

constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;
constexpr uint32_t TELEMETRY_INTERVAL_MS = 100;
// One pulse per loop() pass, so the achievable pulse rate is bounded by how
// fast loop() runs. 20 kHz leaves a wide margin and is far above anything this
// stage is driven at; refusing higher rates is better than silently emitting a
// slower, uneven pulse train than was asked for.
constexpr uint32_t MAX_PULSE_RATE_HZ = 20000;

// Position is held in 1/256 full-step units so that every supported MRES
// divides it exactly and changing MRES cannot introduce a rounding drift.
constexpr int64_t U256_PER_FULLSTEP = 256;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

enum Mode : uint8_t { IDLE, MOVING, JOGGING, CALIB, SEQ };

Mode mode = IDLE;
int64_t posU256 = 0;
int64_t targetU256 = 0;
int64_t jogAU256 = 0;
int64_t jogBU256 = 0;
bool towardB = true;
int32_t limitFullSteps = DEFAULT_LIMIT_FULLSTEPS;
uint32_t rateMilliHz = DEFAULT_RATE_MILLIHZ;
uint16_t mres = DEFAULT_MRES;
uint16_t currentMa = DEFAULT_CURRENT_MA;
uint8_t holdPct = 100;          // standstill current as % of run current
bool driverReady = false;

uint32_t periodUs = 0;
int64_t stepU256 = U256_PER_FULLSTEP / DEFAULT_MRES;
uint64_t nextPulseUs = 0;
bool dirPositive = true;
uint32_t lastTelemetryMs = 0;

pcnt_unit_handle_t pcnt = nullptr;
volatile int32_t indexCount = 0;

// --- microstep calibration -------------------------------------------------
// One microstep forward, a long dwell, one microstep back, another dwell,
// repeated. The dwell is the measurement: it is the only time the stage is
// actually at rest. A +1/-1 pair returns to the same commanded position, so
// whatever the encoder reports at the end of the pair is the stage's error and
// not the command's.
uint32_t calibEndMs = 0;
uint32_t calibDwellMs = 2000;
uint32_t calibNextMs = 0;
uint32_t calibCycles = 0;
bool calibOutbound = true;     // the next microstep is the +1 of the pair

// --- MRES sweep ----------------------------------------------------------
// The calibration oscillation run at each MRES in turn, every block opened by
// a marker jump: J full steps out, a dwell, back, a dwell. The jump is several
// times the largest oscillation (one full step at MRES 1), so block boundaries
// are unmistakable in an encoder trace even without the event log. A final
// jump closes the last block. Every MRES divides a full step, so starting on a
// full step keeps the base point on the grid of every block.
enum SeqPhase : uint8_t {
  SEQ_APPLY, SEQ_JUMP_OUT, SEQ_HOLD_OUT, SEQ_JUMP_BACK, SEQ_HOLD_BACK, SEQ_OSC
};
constexpr uint8_t SEQ_MAX_BLOCKS = 9;
uint16_t seqMres[SEQ_MAX_BLOCKS] = {1, 2, 4, 8, 16, 32, 64};
uint8_t seqCount = 7;
uint8_t seqIdx = 0;            // == seqCount while doing the closing jump
SeqPhase seqPhase = SEQ_APPLY;
int64_t seqBaseU256 = 0;
int64_t seqJumpU256 = 0;
uint32_t seqBlockMs = 30000;
uint32_t seqHoldUntilMs = 0;
uint32_t seqOscEndMs = 0;
uint16_t seqRestoreMres = DEFAULT_MRES;

// --------------------------------------------------------------- encoder

void IRAM_ATTR onIndex() { indexCount++; }

bool setupEncoder() {
  pcnt_unit_config_t unit_cfg = {};
  unit_cfg.high_limit = 30000;
  unit_cfg.low_limit = -30000;
  unit_cfg.flags.accum_count = 1;
  if (pcnt_new_unit(&unit_cfg, &pcnt) != ESP_OK) return false;
  pcnt_glitch_filter_config_t filt = {};
  filt.max_glitch_ns = 1000;
  pcnt_unit_set_glitch_filter(pcnt, &filt);

  // Two channels wired opposite ways round is what makes the count x4 and
  // direction aware: each takes one line as its edge source and the other as
  // the level deciding up or down.
  pcnt_chan_config_t a_cfg = {};
  a_cfg.edge_gpio_num = ENC_A;
  a_cfg.level_gpio_num = ENC_B;
  pcnt_channel_handle_t cha = nullptr;
  if (pcnt_new_channel(pcnt, &a_cfg, &cha) != ESP_OK) return false;
  pcnt_chan_config_t b_cfg = {};
  b_cfg.edge_gpio_num = ENC_B;
  b_cfg.level_gpio_num = ENC_A;
  pcnt_channel_handle_t chb = nullptr;
  if (pcnt_new_channel(pcnt, &b_cfg, &chb) != ESP_OK) return false;

  pcnt_channel_set_edge_action(cha, PCNT_CHANNEL_EDGE_ACTION_DECREASE,
                               PCNT_CHANNEL_EDGE_ACTION_INCREASE);
  pcnt_channel_set_level_action(cha, PCNT_CHANNEL_LEVEL_ACTION_KEEP,
                                PCNT_CHANNEL_LEVEL_ACTION_INVERSE);
  pcnt_channel_set_edge_action(chb, PCNT_CHANNEL_EDGE_ACTION_INCREASE,
                               PCNT_CHANNEL_EDGE_ACTION_DECREASE);
  pcnt_channel_set_level_action(chb, PCNT_CHANNEL_LEVEL_ACTION_KEEP,
                                PCNT_CHANNEL_LEVEL_ACTION_INVERSE);
  gpio_pullup_en(ENC_A);
  gpio_pullup_en(ENC_B);
  // accum_count only carries the count across the limits if the limits are
  // watch points; without these the count silently restarts at 0 every
  // 30000 counts (1.5 shaft revolutions on the 20000 CPR encoder).
  pcnt_unit_add_watch_point(pcnt, 30000);
  pcnt_unit_add_watch_point(pcnt, -30000);
  pcnt_unit_enable(pcnt);
  pcnt_unit_clear_count(pcnt);
  pcnt_unit_start(pcnt);
  pinMode(ENC_Z, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(ENC_Z), onIndex, RISING);
  return true;
}

int32_t encoderCount() {
  int c = 0;
  if (pcnt) pcnt_unit_get_count(pcnt, &c);
  return (int32_t)c;
}

// ---------------------------------------------------------------- driver

// Read a register twice and insist the two agree. A TMC2209 UART read that
// fails returns garbage rather than signalling an error, and v4 lost a 41
// minute run to a single bad read being taken at face value; requiring two
// matching reads is what the campaign sketches settled on.
//
// Takes the read as a callable rather than a pointer-to-member: these getters
// are overloaded (uint32_t GCONF() beside void GCONF(uint32_t)) and declared
// across the library's class hierarchy, so a pointer-to-member spelling would
// depend on overload resolution and base-to-derived conversion going the way
// one hopes. A no-capture lambda sidesteps all of that.
template <typename Reader>
bool readTwice(Reader read, uint32_t &out) {
  for (uint8_t attempt = 0; attempt < 4; ++attempt) {
    const uint32_t a = read();
    const uint32_t b = read();
    if (a == b && a != 0xFFFFFFFFUL) {
      out = a;
      return true;
    }
    delay(5);
  }
  return false;
}

bool applyMres(uint16_t value) {
  uint8_t code = 8;
  for (uint16_t v = value; v > 1; v >>= 1) --code;
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  for (uint8_t attempt = 1; attempt <= 5; ++attempt) {
    uint32_t chopconf = 0;
    if (!readTwice([] { return driver.CHOPCONF(); }, chopconf)) {
      delay(10);
      continue;
    }
    driver.CHOPCONF((chopconf & ~MRES_MASK) |
                    (static_cast<uint32_t>(code) << 24));
    uint32_t verify = 0;
    if (!readTwice([] { return driver.CHOPCONF(); }, verify)) {
      delay(10);
      continue;
    }
    if (((verify & MRES_MASK) >> 24) == code) {
      mres = value;
      stepU256 = U256_PER_FULLSTEP / value;
      return true;
    }
    delay(10);
  }
  return false;
}

bool configureDriver() {
  digitalWrite(EN_PIN, HIGH);    // outputs off until the mode is proven
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(100);
  // Identify the chip before writing anything: a wrong or silent UART shows
  // up here as a version other than 0x21, which says far more than a later
  // readback failure does.
  const uint32_t ioin = driver.IOIN();
  Serial.printf("# DRIVER_ID,ioin=0x%08lX,version=0x%02X,expected=0x21\n",
                (unsigned long)ioin, (unsigned)((ioin >> 24) & 0xFF));
  driver.begin();
  driver.toff(5);
  driver.blank_time(24);
  // Current from IRUN/IHOLD alone, not scaled by the module's VREF pot
  // (I_scale_analog resets to 1). Hold current defaults to run current with no
  // power-down, so the dwells in CALIB/SEQ are held at full current; HOLD
  // changes the ratio deliberately.
  driver.I_scale_analog(false);
  driver.hold_multiplier(holdPct / 100.0F);
  driver.rms_current(currentMa);
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);
  driver.en_spreadCycle(true);   // default is StealthChop; we want current mode
  driver.TPWMTHRS(0);            // never hand back to StealthChop
  driver.intpol(false);          // default is ON, which smears every microstep
  driver.pwm_autoscale(false);
  if (!applyMres(mres)) {
    Serial.println("# ERR,reason=mres_write_failed");
    return false;
  }

  // Prove the configuration rather than assume it: a run that gets past here
  // is evidence the driver is in the mode we asked for.
  uint32_t gconf = 0, chopconf = 0, drvStatus = 0;
  if (!readTwice([] { return driver.GCONF(); }, gconf) ||
      !readTwice([] { return driver.CHOPCONF(); }, chopconf) ||
      !readTwice([] { return driver.DRV_STATUS(); }, drvStatus)) {
    Serial.println("# ERR,reason=uart_readback_failed");
    return false;
  }
  const bool analog = (gconf & 1UL) != 0;
  const bool spread = (gconf & (1UL << 2)) != 0;
  const bool intpol = (chopconf & (1UL << 28)) != 0;
  const bool stealth = (drvStatus & (1UL << 30)) != 0;
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,drv_status=0x%08lX,"
                "spreadcycle=%u,intpol=%u,stealth=%u,i_scale_analog=%u,"
                "irun=%u,ihold=%u\n",
                (unsigned long)gconf, (unsigned long)chopconf,
                (unsigned long)drvStatus, spread, intpol, stealth, analog,
                driver.irun(), driver.ihold());
  if (!spread || intpol || stealth || analog) {
    Serial.println("# ERR,reason=driver_mode_not_confirmed");
    return false;
  }
  digitalWrite(EN_PIN, LOW);     // enabled: the motor now holds position
  return true;
}

// ---------------------------------------------------------------- motion

bool setRate(uint32_t milliHz) {
  if (milliHz == 0) return false;
  const uint64_t pulseMilliHz = static_cast<uint64_t>(milliHz) * mres;
  if (pulseMilliHz > static_cast<uint64_t>(MAX_PULSE_RATE_HZ) * 1000ULL) return false;
  const uint64_t p = (1000000000ULL + pulseMilliHz / 2) / pulseMilliHz;
  if (p <= STEP_HIGH_US + 2 || p > 0xFFFFFFFFULL) return false;
  rateMilliHz = milliHz;
  periodUs = static_cast<uint32_t>(p);
  return true;
}

bool withinLimit(int64_t u256) {
  if (limitFullSteps <= 0) return true;
  const int64_t lim = static_cast<int64_t>(limitFullSteps) * U256_PER_FULLSTEP;
  return u256 >= -lim && u256 <= lim;
}

void stopMotion(const char *why) {
  mode = IDLE;
  Serial.printf("# STOP,pos=%ld,pos_u256=%lld,reason=%s\n",
                (long)(posU256 / U256_PER_FULLSTEP), (long long)posU256, why);
}

void emitOneMicrostep(bool positive) {
  if (positive != dirPositive) {
    digitalWrite(DIR_PIN, positive ? HIGH : LOW);
    dirPositive = positive;
    delayMicroseconds(DIR_SETUP_US);
  }
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(STEP_HIGH_US);
  digitalWrite(STEP_PIN, LOW);
  posU256 += positive ? stepU256 : -stepU256;
}

void serviceCalibration() {
  const uint32_t now = millis();
  if (now >= calibEndMs) {
    Serial.printf("# CALIB_DONE,cycles=%lu,pos=%ld,enc=%ld\n",
                  (unsigned long)calibCycles,
                  (long)(posU256 / U256_PER_FULLSTEP), (long)encoderCount());
    mode = IDLE;
    return;
  }
  if (now < calibNextMs) return;
  emitOneMicrostep(calibOutbound);
  calibNextMs = now + calibDwellMs;
  if (!calibOutbound) calibCycles++;
  Serial.printf("# CALIB,%c,cycle=%lu,pos_u256=%lld,enc=%ld,t_ms=%lu\n",
                calibOutbound ? '+' : '-', (unsigned long)calibCycles,
                (long long)posU256, (long)encoderCount(), (unsigned long)now);
  calibOutbound = !calibOutbound;
}

// One pass of non-blocking step generation toward `target`: at most one
// pulse. True once the stage is there. Stops on a soft limit, so callers must
// check `mode` if it returns false.
bool stepToward(int64_t target) {
  if (posU256 == target) return true;
  const bool positive = target > posU256;
  if (positive != dirPositive) {
    // Give the driver its direction setup time before the next edge.
    digitalWrite(DIR_PIN, positive ? HIGH : LOW);
    dirPositive = positive;
    nextPulseUs = esp_timer_get_time() + DIR_SETUP_US;
    return false;
  }

  const uint64_t now = esp_timer_get_time();
  if (now < nextPulseUs) return false;

  const int64_t next = posU256 + (positive ? stepU256 : -stepU256);
  if (!withinLimit(next)) {
    stopMotion("soft_limit");
    return false;
  }
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(STEP_HIGH_US);
  digitalWrite(STEP_PIN, LOW);
  posU256 = next;

  nextPulseUs += periodUs;
  // If the loop fell behind, resynchronise instead of trying to catch up with
  // a burst of back-to-back pulses the motor cannot follow.
  if (nextPulseUs < now) nextPulseUs = now + periodUs;
  return false;
}

void seqEvent(const char *what, const char *extra) {
  Serial.printf("# %s,block=%u,mres=%u,%spos_u256=%lld,enc=%ld,t_ms=%lu\n",
                what, seqIdx, mres, extra, (long long)posU256,
                (long)encoderCount(), (unsigned long)millis());
}

void serviceSequence() {
  const uint32_t now = millis();
  switch (seqPhase) {
    case SEQ_APPLY:
      if (seqIdx < seqCount) {
        if (!applyMres(seqMres[seqIdx]) || !setRate(rateMilliHz)) {
          Serial.println("# ERR,reason=seq_mres_failed");
          stopMotion("seq_mres_failed");
          return;
        }
        char extra[32];
        snprintf(extra, sizeof extra, "microstep_um=%.4f,", 10.0 / mres);
        seqEvent("SEQ_BLOCK", extra);
      }
      seqPhase = SEQ_JUMP_OUT;
      nextPulseUs = esp_timer_get_time();
      return;
    case SEQ_JUMP_OUT:
      if (!stepToward(seqBaseU256 + seqJumpU256)) return;
      seqEvent("SEQ_JUMP", "dir=out,");
      seqHoldUntilMs = now + calibDwellMs;
      seqPhase = SEQ_HOLD_OUT;
      return;
    case SEQ_HOLD_OUT:
      if (now < seqHoldUntilMs) return;
      seqPhase = SEQ_JUMP_BACK;
      nextPulseUs = esp_timer_get_time();
      return;
    case SEQ_JUMP_BACK:
      if (!stepToward(seqBaseU256)) return;
      seqEvent("SEQ_JUMP", "dir=back,");
      seqHoldUntilMs = now + calibDwellMs;
      seqPhase = SEQ_HOLD_BACK;
      return;
    case SEQ_HOLD_BACK:
      if (now < seqHoldUntilMs) return;
      if (seqIdx >= seqCount) {          // that was the closing jump
        applyMres(seqRestoreMres);
        setRate(rateMilliHz);
        Serial.printf("# SEQ_DONE,blocks=%u,pos_u256=%lld,enc=%ld,t_ms=%lu,"
                      "mres_restored=%u\n", seqCount, (long long)posU256,
                      (long)encoderCount(), (unsigned long)now, mres);
        mode = IDLE;
        return;
      }
      seqOscEndMs = now + seqBlockMs;
      calibNextMs = now;
      calibCycles = 0;
      calibOutbound = true;
      seqPhase = SEQ_OSC;
      return;
    case SEQ_OSC:
      // End only on a completed pair, so every block leaves the stage back
      // at the base point for the next marker jump.
      if (calibOutbound && now >= seqOscEndMs) {
        seqEvent("SEQ_OSC_END", "");   // the last dwell's settled reading
        seqIdx++;
        seqPhase = SEQ_APPLY;
        return;
      }
      if (now < calibNextMs) return;
      // Logged before the pulse: pos_u256/enc are where the stage settled at
      // the end of the dwell, and `next` is the microstep about to be made.
      {
        char extra[40];
        snprintf(extra, sizeof extra, "next=%c,cycle=%lu,",
                 calibOutbound ? '+' : '-',
                 (unsigned long)(calibCycles + 1));   // pair number
        seqEvent("SEQ_OSC", extra);
      }
      emitOneMicrostep(calibOutbound);
      calibNextMs = now + calibDwellMs;
      if (!calibOutbound) calibCycles++;
      calibOutbound = !calibOutbound;
      return;
  }
}

void serviceMotion() {
  if (mode == CALIB) { serviceCalibration(); return; }
  if (mode == SEQ) { serviceSequence(); return; }
  if (mode == IDLE || periodUs == 0) return;
  const int64_t target = (mode == JOGGING) ? (towardB ? jogBU256 : jogAU256)
                                           : targetU256;
  if (posU256 == target) {
    if (mode == JOGGING) {
      Serial.printf("# REACHED,pos=%ld,pos_u256=%lld,next=%c\n",
                    (long)(posU256 / U256_PER_FULLSTEP), (long long)posU256,
                    towardB ? 'A' : 'B');
      towardB = !towardB;
      return;
    }
    Serial.printf("# ARRIVED,pos=%ld,pos_u256=%lld\n",
                  (long)(posU256 / U256_PER_FULLSTEP), (long long)posU256);
    mode = IDLE;
    return;
  }
  stepToward(target);
}

const char *modeName() {
  switch (mode) {
    case MOVING: return "MOVING";
    case JOGGING: return "JOGGING";
    case CALIB: return "CALIB";
    case SEQ: return "SEQ";
    default: return "IDLE";
  }
}

// --------------------------------------------------------------- commands

void sendStatus() {
  Serial.printf("# STATUS,pos=%ld,pos_u256=%lld,a_u256=%lld,b_u256=%lld,"
                "target=%ld,a=%ld,b=%ld,mode=%s,"
                "rate_millihz=%lu,mres=%u,current_ma=%u,hold_pct=%u,limit=%ld,"
                "steps_per_mm=%u,enc=%ld,index=%ld\n",
                (long)(posU256 / U256_PER_FULLSTEP), (long long)posU256,
                (long long)jogAU256, (long long)jogBU256,
                (long)(targetU256 / U256_PER_FULLSTEP),
                (long)(jogAU256 / U256_PER_FULLSTEP),
                (long)(jogBU256 / U256_PER_FULLSTEP),
                modeName(),
                (unsigned long)rateMilliHz, mres, currentMa, holdPct,
                (long)limitFullSteps, FULL_STEPS_PER_MM,
                (long)encoderCount(), (long)indexCount);
}

void handleCommand(String line) {
  line.trim();
  if (line.length() == 0) return;
  String verb = line;
  String rest = "";
  const int sp = line.indexOf(' ');
  if (sp > 0) {
    verb = line.substring(0, sp);
    rest = line.substring(sp + 1);
    rest.trim();
  }
  verb.toUpperCase();

  if (verb == "PING") {
    Serial.printf("# PONG,fw=esp32s3_tmc2209_jog,ready=%u\n", driverReady);
    return;
  }
  if (verb == "STATUS") { sendStatus(); return; }
  if (verb == "POS?") {
    Serial.printf("# POS,pos=%ld,pos_u256=%lld\n",
                  (long)(posU256 / U256_PER_FULLSTEP), (long long)posU256);
    return;
  }
  if (verb == "STOP") { stopMotion("command"); return; }
  if (verb == "REINIT") {
    if (mode != IDLE) { Serial.println("# ERR,reason=busy"); return; }
    driverReady = configureDriver();
    Serial.printf("# REINIT_DONE,ready=%u\n", driverReady);
    return;
  }
  if (verb == "ENC?") {
    Serial.printf("# ENC,enc=%ld,index=%ld\n", (long)encoderCount(),
                  (long)indexCount);
    return;
  }
  if (verb == "CALIB") {
    // CALIB <duration_s> [dwell_ms] - the microstep calibration oscillation.
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    long secs = 30, dwell = 2000;
    sscanf(rest.c_str(), "%ld %ld", &secs, &dwell);
    if (secs <= 0 || secs > 3600 || dwell < 50 || dwell > 60000) {
      Serial.println("# ERR,reason=bad_calib_args");
      return;
    }
    calibDwellMs = (uint32_t)dwell;
    calibEndMs = millis() + (uint32_t)secs * 1000UL;
    calibNextMs = millis();
    calibCycles = 0;
    calibOutbound = true;
    mode = CALIB;
    Serial.printf("# CALIB_START,seconds=%ld,dwell_ms=%ld,mres=%u,"
                  "microstep_um=%.4f\n", secs, dwell, mres, 10.0 / mres);
    return;
  }

  if (verb == "SEQ") {
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    if (mode != IDLE) { Serial.println("# ERR,reason=busy"); return; }
    long secs = 30, dwell = 2000, jump = 4;
    char list[64] = {0};
    sscanf(rest.c_str(), "%ld %ld %ld %63s", &secs, &dwell, &jump, list);
    if (secs <= 0 || secs > 3600 || dwell < 50 || dwell > 60000 ||
        jump == 0 || jump > 10000 || jump < -10000) {
      Serial.println("# ERR,reason=bad_seq_args");
      return;
    }
    uint16_t parsed[SEQ_MAX_BLOCKS];
    uint8_t count = 0;
    if (list[0]) {
      for (char *tok = strtok(list, ","); tok; tok = strtok(nullptr, ",")) {
        const long m = atol(tok);
        if (count >= SEQ_MAX_BLOCKS || m < 1 || m > 256 || (m & (m - 1))) {
          Serial.println("# ERR,reason=bad_seq_mres_list");
          return;
        }
        parsed[count++] = (uint16_t)m;
      }
    } else {
      const uint16_t def[] = {1, 2, 4, 8, 16, 32, 64};
      for (uint16_t m : def) parsed[count++] = m;
    }
    uint16_t maxMres = 1;
    for (uint8_t i = 0; i < count; ++i) if (parsed[i] > maxMres) maxMres = parsed[i];
    if ((uint64_t)rateMilliHz * maxMres > (uint64_t)MAX_PULSE_RATE_HZ * 1000ULL) {
      Serial.printf("# ERR,reason=rate_too_high_for_mres,max_mres=%u\n", maxMres);
      return;
    }
    // Every MRES divides a full step, so a full-step base is on every grid.
    if (posU256 % U256_PER_FULLSTEP != 0) {
      Serial.println("# ERR,reason=not_on_full_step");
      return;
    }
    if (!withinLimit(posU256 + (int64_t)jump * U256_PER_FULLSTEP)) {
      Serial.println("# ERR,reason=outside_limit");
      return;
    }
    for (uint8_t i = 0; i < count; ++i) seqMres[i] = parsed[i];
    seqCount = count;
    seqIdx = 0;
    seqPhase = SEQ_APPLY;
    seqBaseU256 = posU256;
    seqJumpU256 = (int64_t)jump * U256_PER_FULLSTEP;
    seqBlockMs = (uint32_t)secs * 1000UL;
    calibDwellMs = (uint32_t)dwell;
    seqRestoreMres = mres;
    mode = SEQ;
    Serial.printf("# SEQ_START,blocks=%u,block_s=%ld,dwell_ms=%ld,jump=%ld,"
                  "base_u256=%lld,rate_millihz=%lu,t_ms=%lu\n", count, secs,
                  dwell, jump, (long long)posU256, (unsigned long)rateMilliHz,
                  (unsigned long)millis());
    return;
  }
  if (verb == "ZERO") {
    if (mode != IDLE) { Serial.println("# ERR,reason=busy"); return; }
    posU256 = 0;
    targetU256 = 0;
    Serial.println("# ZERO_OK,pos=0");
    return;
  }
  if (verb == "LIMIT") {
    const long v = rest.toInt();
    if (v < 0) { Serial.println("# ERR,reason=bad_limit"); return; }
    limitFullSteps = static_cast<int32_t>(v);
    Serial.printf("# LIMIT_OK,limit=%ld\n", (long)limitFullSteps);
    return;
  }
  if (verb == "RATE") {
    const long v = rest.toInt();
    if (v <= 0 || !setRate(static_cast<uint32_t>(v))) {
      Serial.println("# ERR,reason=bad_rate");
      return;
    }
    Serial.printf("# RATE_OK,rate_millihz=%lu,pulse_hz=%.1f\n",
                  (unsigned long)rateMilliHz, (double)rateMilliHz * mres / 1000.0);
    return;
  }
  if (verb == "CURRENT") {
    const long v = rest.toInt();
    if (v < 50 || v > 1200) { Serial.println("# ERR,reason=bad_current"); return; }
    currentMa = static_cast<uint16_t>(v);
    driver.rms_current(currentMa);
    Serial.printf("# CURRENT_OK,current_ma=%u\n", currentMa);
    return;
  }
  if (verb == "HOLD") {
    const long v = rest.toInt();
    if (rest.length() == 0 || v < 0 || v > 100) {
      Serial.println("# ERR,reason=bad_hold");
      return;
    }
    holdPct = static_cast<uint8_t>(v);
    driver.hold_multiplier(holdPct / 100.0F);
    driver.rms_current(currentMa);
    Serial.printf("# HOLD_OK,hold_pct=%u,irun=%u,ihold=%u\n", holdPct,
                  driver.irun(), driver.ihold());
    return;
  }
  if (verb == "MRES") {
    if (mode != IDLE) { Serial.println("# ERR,reason=busy"); return; }
    const long v = rest.toInt();
    bool power_of_two = v >= 1 && v <= 256 && (v & (v - 1)) == 0;
    if (!power_of_two) { Serial.println("# ERR,reason=bad_mres"); return; }
    if (!applyMres(static_cast<uint16_t>(v))) {
      Serial.println("# ERR,reason=mres_write_failed");
      return;
    }
    // Snap to the new microstep grid and re-derive the pulse period, which
    // depends on MRES as well as the full-step rate.
    posU256 = (posU256 / stepU256) * stepU256;
    if (!setRate(rateMilliHz)) setRate(DEFAULT_RATE_MILLIHZ);
    Serial.printf("# MRES_OK,mres=%u,pulse_hz=%.1f\n",
                  mres, (double)rateMilliHz * mres / 1000.0);
    return;
  }
  if (verb == "MOVETO") {
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    const int64_t t = static_cast<int64_t>(rest.toInt()) * U256_PER_FULLSTEP;
    if (!withinLimit(t)) { Serial.println("# ERR,reason=outside_limit"); return; }
    targetU256 = t;
    mode = MOVING;
    nextPulseUs = esp_timer_get_time();
    Serial.printf("# MOVETO,target=%ld\n", (long)(t / U256_PER_FULLSTEP));
    return;
  }
  if (verb == "STEP") {
    // Relative move of n microsteps at whatever MRES is set, so 'STEP 1' is
    // the smallest motion the driver can make.
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    if (mode != IDLE) { Serial.println("# ERR,reason=busy"); return; }
    const long n = rest.toInt();
    if (n == 0) { Serial.println("# ERR,reason=bad_step_count"); return; }
    const int64_t t = posU256 + static_cast<int64_t>(n) * stepU256;
    if (!withinLimit(t)) { Serial.println("# ERR,reason=outside_limit"); return; }
    targetU256 = t;
    mode = MOVING;
    nextPulseUs = esp_timer_get_time();
    Serial.printf("# STEP,n=%ld,mres=%u,target_u256=%lld\n", n, mres,
                  (long long)t);
    return;
  }
  if (verb == "JOGU") {
    // Endpoints in 1/256 full-step units. They must sit on the present
    // microstep grid, or the stage could never land exactly on them and would
    // dither either side of the target instead of reversing.
    long long a = 0, b = 0;
    if (sscanf(rest.c_str(), "%lld %lld", &a, &b) != 2) {
      Serial.println("# ERR,reason=jog_needs_two");
      return;
    }
    if (a % stepU256 != 0 || b % stepU256 != 0) {
      Serial.println("# ERR,reason=off_microstep_grid");
      return;
    }
    if (!withinLimit(a) || !withinLimit(b)) {
      Serial.println("# ERR,reason=outside_limit");
      return;
    }
    if (a == b) { Serial.println("# ERR,reason=equal_positions"); return; }
    jogAU256 = a;
    jogBU256 = b;
    Serial.printf("# JOG_SET,a_u256=%lld,b_u256=%lld\n", a, b);
    return;
  }
  if (verb == "JOG") {
    const int sp2 = rest.indexOf(' ');
    if (sp2 <= 0) { Serial.println("# ERR,reason=jog_needs_two"); return; }
    const int64_t a = static_cast<int64_t>(rest.substring(0, sp2).toInt()) * U256_PER_FULLSTEP;
    const int64_t b = static_cast<int64_t>(rest.substring(sp2 + 1).toInt()) * U256_PER_FULLSTEP;
    if (!withinLimit(a) || !withinLimit(b)) {
      Serial.println("# ERR,reason=outside_limit");
      return;
    }
    if (a == b) { Serial.println("# ERR,reason=equal_positions"); return; }
    jogAU256 = a;
    jogBU256 = b;
    Serial.printf("# JOG_SET,a=%ld,b=%ld\n", (long)(a / U256_PER_FULLSTEP),
                  (long)(b / U256_PER_FULLSTEP));
    return;
  }
  if (verb == "START") {
    if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
    if (jogAU256 == jogBU256) { Serial.println("# ERR,reason=no_jog_set"); return; }
    // Head for whichever end is further away, so starting from between the two
    // positions gives a full first leg rather than a stub.
    const int64_t dB = jogBU256 - posU256;
    const int64_t dA = jogAU256 - posU256;
    towardB = (dB < 0 ? -dB : dB) >= (dA < 0 ? -dA : dA);
    mode = JOGGING;
    nextPulseUs = esp_timer_get_time();
    Serial.printf("# START,a=%ld,b=%ld,first=%c\n",
                  (long)(jogAU256 / U256_PER_FULLSTEP),
                  (long)(jogBU256 / U256_PER_FULLSTEP), towardB ? 'B' : 'A');
    return;
  }
  Serial.printf("# ERR,reason=unknown_command,verb=%s\n", verb.c_str());
}

void pollSerial() {
  static String buffer;
  while (Serial.available()) {
    const char c = static_cast<char>(Serial.read());
    if (c == '\n' || c == '\r') {
      if (buffer.length()) {
        handleCommand(buffer);
        buffer = "";
      }
    } else if (buffer.length() < 120) {
      buffer += c;
    }
  }
}

// ------------------------------------------------------------------ setup

void setup() {
  pinMode(EN_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, HIGH);
  pinMode(LED_PIN, OUTPUT);
  digitalWrite(LED_PIN, LOW);

  Serial.begin(115200);
  const uint32_t waitStart = millis();
  while (!Serial && millis() - waitStart < 3000) delay(10);
  delay(200);

  Serial.println("# BOOT,fw=esp32s3_tmc2209_jog");
  Serial.printf("# PINS,step=%u,dir=%u,en=%u,uart_rx=%u,uart_tx=%u,"
                "enc_a=%d,enc_b=%d,enc_z=%d\n",
                STEP_PIN, DIR_PIN, EN_PIN, UART_RX_PIN, UART_TX_PIN,
                (int)ENC_A, (int)ENC_B, (int)ENC_Z);
  setRate(rateMilliHz);
  const bool encOk = setupEncoder();
  if (!encOk) Serial.println("# ERR,reason=pcnt_init_failed");
  driverReady = configureDriver();
  Serial.printf("# READY,ready=%u,encoder=%u,mres=%u,current_ma=%u,"
                "rate_millihz=%lu,limit=%ld,steps_per_mm=%u\n",
                driverReady, encOk, mres, currentMa, (unsigned long)rateMilliHz,
                (long)limitFullSteps, FULL_STEPS_PER_MM);
  Serial.println("# IDLE,nothing moves until commanded");
}

void loop() {
  pollSerial();
  serviceMotion();
  const uint32_t now = millis();
  if (mode != IDLE && now - lastTelemetryMs >= TELEMETRY_INTERVAL_MS) {
    lastTelemetryMs = now;
    Serial.printf("# POS,pos=%ld,pos_u256=%lld,enc=%ld,mode=%s\n",
                  (long)(posU256 / U256_PER_FULLSTEP), (long long)posU256,
                  (long)encoderCount(),
                  modeName());
  }
}

/*
 * Stage identification campaign runner — ESP32-S3 + TMC2209 + quadrature encoder.
 *
 * Executes the sequence defined in scripts/sequence_design.py and streams the
 * commanded and measured position throughout. The plan is **not** compiled in:
 * the host uploads it move by move over USB CDC before starting. One definition
 * of the campaign therefore exists, in Python, and this sketch is the thing that
 * carries it out — editing the plan never means reflashing.
 *
 *   PLAN CLEAR
 *   PLAN ADD <target> <speed> <dwell_ms> <mres> <micro_n> <micro_dwell_ms> <tag>
 *   PLAN RUN
 *
 * Positions are whole full steps from wherever ZERO was last called; 100 full
 * steps = 1 mm. Speeds are full steps/s. A move with micro_n > 0 is a staircase
 * leg: micro_n single microsteps, each followed by micro_dwell_ms, instead of a
 * continuous ramp.
 *
 * WIRING — see ../../../CONNECTIONS.md. The encoder has moved off GPIO5/6
 * because the TMC2209 uses those for STEP/DIR; both cannot have them.
 *
 *   STEP GPIO5   DIR GPIO6   TMC UART TX GPIO17 / RX GPIO18   EN grounded
 *   ENCODER  A GPIO4   B GPIO7   Z GPIO15
 *
 * The encoder is read by the PCNT peripheral, so counts are never missed
 * regardless of what loop() is doing — which matters here, because loop() is
 * also generating the step pulses.
 *
 * Nothing moves until PLAN RUN. ABORT stops immediately and is checked between
 * every pulse.
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/pulse_cnt.h"

// ---- pins (see CONNECTIONS.md) --------------------------------------------
constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr gpio_num_t ENC_A = GPIO_NUM_4;
constexpr gpio_num_t ENC_B = GPIO_NUM_7;
constexpr gpio_num_t ENC_Z = GPIO_NUM_15;

constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.11F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t DEFAULT_CURRENT_MA = 360;
constexpr uint16_t FULL_STEPS_PER_MM = 100;

constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;
constexpr int64_t U256 = 256;            // position unit = 1/256 full step
constexpr size_t MAX_MOVES = 512;
constexpr uint32_t SAMPLE_INTERVAL_MS = 20;   // 50 Hz telemetry
constexpr uint32_t MAX_PULSE_RATE_HZ = 20000;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);
pcnt_unit_handle_t pcnt = nullptr;

struct MoveSpec {
  int32_t target;          // full steps
  uint16_t speed;          // full steps/s
  uint16_t dwell_ms;
  uint16_t mres;
  uint16_t micro_n;        // >0 => staircase leg
  uint16_t micro_dwell_ms;
  char tag[20];
};

MoveSpec plan[MAX_MOVES];
size_t planCount = 0;

bool driverReady = false;
bool running = false;
bool abortRequested = false;
int64_t posU256 = 0;
uint16_t mres = 16;
uint16_t currentMa = DEFAULT_CURRENT_MA;
int64_t stepU256 = U256 / 16;
volatile int32_t indexCount = 0;

// ---------------------------------------------------------------- encoder

void IRAM_ATTR onIndex() { indexCount++; }

bool setupEncoder() {
  pcnt_unit_config_t unit_cfg = {};
  unit_cfg.high_limit = 30000;
  unit_cfg.low_limit = -30000;
  unit_cfg.flags.accum_count = 1;          // keep counting past the limits
  if (pcnt_new_unit(&unit_cfg, &pcnt) != ESP_OK) return false;

  pcnt_glitch_filter_config_t filt = {};
  filt.max_glitch_ns = 1000;
  pcnt_unit_set_glitch_filter(pcnt, &filt);

  // Two channels wired opposite ways round is what makes the count x4 and
  // direction-aware: each channel takes one line as its edge source and the
  // other as the level that decides up or down.
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
  pcnt_unit_enable(pcnt);
  pcnt_unit_clear_count(pcnt);
  pcnt_unit_start(pcnt);

  pinMode(ENC_Z, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(ENC_Z), onIndex, RISING);
  return true;
}

int32_t encoderCount() {
  int count = 0;
  if (pcnt) pcnt_unit_get_count(pcnt, &count);
  return (int32_t)count;
}

// ----------------------------------------------------------------- driver

template <typename Reader>
bool readTwice(Reader read, uint32_t &out) {
  for (uint8_t a = 0; a < 4; ++a) {
    const uint32_t x = read(), y = read();
    if (x == y && x != 0xFFFFFFFFUL) { out = x; return true; }
    delay(5);
  }
  return false;
}

bool applyMres(uint16_t value) {
  uint8_t code = 8;
  for (uint16_t v = value; v > 1; v >>= 1) --code;
  constexpr uint32_t MASK = 0x0F000000UL;
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    uint32_t chop = 0;
    if (!readTwice([] { return driver.CHOPCONF(); }, chop)) { delay(10); continue; }
    driver.CHOPCONF((chop & ~MASK) | ((uint32_t)code << 24));
    uint32_t back = 0;
    if (!readTwice([] { return driver.CHOPCONF(); }, back)) { delay(10); continue; }
    if (((back & MASK) >> 24) == code) {
      mres = value;
      stepU256 = U256 / value;
      return true;
    }
    delay(10);
  }
  return false;
}

bool configureDriver() {
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(20);
  driver.begin();
  driver.toff(5);
  driver.blank_time(24);
  driver.rms_current(currentMa);
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);
  driver.pwm_autoscale(false);
  if (!applyMres(mres)) { Serial.println("# ERR,reason=mres_write_failed"); return false; }

  uint32_t g = 0, c = 0, d = 0;
  if (!readTwice([] { return driver.GCONF(); }, g) ||
      !readTwice([] { return driver.CHOPCONF(); }, c) ||
      !readTwice([] { return driver.DRV_STATUS(); }, d)) {
    Serial.println("# ERR,reason=uart_readback_failed");
    return false;
  }
  const bool spread = g & (1UL << 2), intpol = c & (1UL << 28), stealth = d & (1UL << 30);
  Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,drv_status=0x%08lX,"
                "spreadcycle=%u,intpol=%u,stealth=%u\n",
                (unsigned long)g, (unsigned long)c, (unsigned long)d,
                spread, intpol, stealth);
  if (!spread || intpol || stealth) {
    Serial.println("# ERR,reason=driver_mode_not_confirmed");
    return false;
  }
  return true;
}

// ------------------------------------------------------------------ motion

bool pollAbort();

void emitSample(const char *phase) {
  Serial.printf("# S,%lu,%ld,%ld,%ld\n",
                (unsigned long)(esp_timer_get_time() / 1000),
                (long)(posU256 / U256), (long)encoderCount(), (long)indexCount);
  (void)phase;
}

/** Emit `count` microsteps at `periodUs`, sampling as it goes. */
bool pulseTrain(uint32_t count, bool positive, uint32_t periodUs) {
  digitalWrite(DIR_PIN, positive ? HIGH : LOW);
  delayMicroseconds(DIR_SETUP_US);
  static uint32_t lastSample = 0;
  for (uint32_t i = 0; i < count; ++i) {
    const uint64_t t0 = esp_timer_get_time();
    digitalWrite(STEP_PIN, HIGH);
    delayMicroseconds(STEP_HIGH_US);
    digitalWrite(STEP_PIN, LOW);
    posU256 += positive ? stepU256 : -stepU256;
    while ((uint32_t)(esp_timer_get_time() - t0) < periodUs) {
      const uint32_t now = millis();
      if (now - lastSample >= SAMPLE_INTERVAL_MS) { lastSample = now; emitSample("move"); }
      if ((i & 31U) == 31U && pollAbort()) return false;
    }
    if (pollAbort()) return false;
  }
  return true;
}

bool dwell(uint32_t ms) {
  const uint32_t t0 = millis();
  uint32_t lastSample = 0;
  while (millis() - t0 < ms) {
    if (millis() - lastSample >= SAMPLE_INTERVAL_MS) { lastSample = millis(); emitSample("dwell"); }
    if (pollAbort()) return false;
    delay(1);
  }
  return true;
}

bool runMove(const MoveSpec &m, size_t index) {
  if (m.mres != mres && !applyMres(m.mres)) return false;
  const int64_t target = (int64_t)m.target * U256;
  const int64_t delta = target - posU256;
  if (delta == 0) return dwell(m.dwell_ms);

  const bool positive = delta > 0;
  const uint64_t absd = positive ? delta : -delta;
  const uint32_t pulses = (uint32_t)(absd / stepU256);

  Serial.printf("# MOVE,%u,%s,target=%ld,speed=%u,mres=%u,pulses=%lu,micro=%u\n",
                (unsigned)index, m.tag, (long)m.target, m.speed, m.mres,
                (unsigned long)pulses, m.micro_n);

  if (m.micro_n > 0) {
    // Staircase leg: one microstep, then hold, repeated. The hold is the
    // measurement - a continuous ramp would never let the stage come to rest.
    for (uint16_t k = 0; k < m.micro_n; ++k) {
      if (!pulseTrain(1, positive, 1000)) return false;
      if (!dwell(m.micro_dwell_ms)) return false;
    }
  } else {
    const uint64_t pulseHz = (uint64_t)m.speed * mres;
    if (pulseHz == 0 || pulseHz > MAX_PULSE_RATE_HZ) {
      Serial.println("# ERR,reason=rate_out_of_range");
      return false;
    }
    const uint32_t periodUs = (uint32_t)(1000000ULL / pulseHz);
    if (!pulseTrain(pulses, positive, periodUs)) return false;
    if (!dwell(m.dwell_ms)) return false;
  }
  Serial.printf("# ARRIVED,%u,pos=%ld,enc=%ld\n",
                (unsigned)index, (long)(posU256 / U256), (long)encoderCount());
  return true;
}

// ---------------------------------------------------------------- commands

String inBuf;

bool pollAbort() {
  while (Serial.available()) {
    const char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      inBuf.trim();
      if (inBuf.equalsIgnoreCase("ABORT")) {
        abortRequested = true;
        inBuf = "";
        Serial.println("# ABORT");
        return true;
      }
      inBuf = "";
    } else if (inBuf.length() < 120) {
      inBuf += c;
    }
  }
  return abortRequested;
}

void runPlan() {
  if (!driverReady) { Serial.println("# ERR,reason=driver_not_ready"); return; }
  if (planCount == 0) { Serial.println("# ERR,reason=plan_empty"); return; }
  running = true;
  abortRequested = false;
  digitalWrite(LED_PIN, HIGH);
  Serial.printf("# RUN_START,moves=%u,pos=%ld,enc=%ld\n",
                (unsigned)planCount, (long)(posU256 / U256), (long)encoderCount());
  for (size_t i = 0; i < planCount && !abortRequested; ++i) {
    if (!runMove(plan[i], i)) break;
  }
  digitalWrite(LED_PIN, LOW);
  running = false;
  Serial.printf("# RUN_%s,pos=%ld,enc=%ld\n",
                abortRequested ? "ABORTED" : "COMPLETE",
                (long)(posU256 / U256), (long)encoderCount());
}

void handleCommand(String line) {
  line.trim();
  if (line.length() == 0) return;
  String verb = line;
  String rest = "";
  const int sp = line.indexOf(' ');
  if (sp > 0) { verb = line.substring(0, sp); rest = line.substring(sp + 1); rest.trim(); }
  verb.toUpperCase();

  if (verb == "PING") {
    Serial.printf("# PONG,fw=stage_identification,ready=%u,moves=%u\n",
                  driverReady, (unsigned)planCount);
  } else if (verb == "STATUS") {
    Serial.printf("# STATUS,pos=%ld,enc=%ld,index=%ld,mres=%u,current_ma=%u,"
                  "moves=%u,running=%u,steps_per_mm=%u\n",
                  (long)(posU256 / U256), (long)encoderCount(), (long)indexCount,
                  mres, currentMa, (unsigned)planCount, running, FULL_STEPS_PER_MM);
  } else if (verb == "ZERO") {
    posU256 = 0;
    if (pcnt) pcnt_unit_clear_count(pcnt);
    indexCount = 0;
    Serial.println("# ZERO_OK,pos=0,enc=0");
  } else if (verb == "CURRENT") {
    currentMa = (uint16_t)rest.toInt();
    driver.rms_current(currentMa);
    Serial.printf("# CURRENT_OK,current_ma=%u\n", currentMa);
  } else if (verb == "PLAN") {
    String sub = rest;
    int s2 = rest.indexOf(' ');
    String args = "";
    if (s2 > 0) { sub = rest.substring(0, s2); args = rest.substring(s2 + 1); }
    sub.toUpperCase();
    if (sub == "CLEAR") {
      planCount = 0;
      Serial.println("# PLAN_CLEARED");
    } else if (sub == "ADD") {
      if (planCount >= MAX_MOVES) { Serial.println("# ERR,reason=plan_full"); return; }
      MoveSpec m = {};
      char tag[20] = "";
      const int n = sscanf(args.c_str(), "%ld %hu %hu %hu %hu %hu %19s",
                           (long *)&m.target, &m.speed, &m.dwell_ms, &m.mres,
                           &m.micro_n, &m.micro_dwell_ms, tag);
      if (n < 6) { Serial.println("# ERR,reason=plan_add_fields"); return; }
      strncpy(m.tag, tag, sizeof(m.tag) - 1);
      plan[planCount++] = m;
      if (planCount % 50 == 0)
        Serial.printf("# PLAN_ADD,count=%u\n", (unsigned)planCount);
    } else if (sub == "RUN") {
      runPlan();
    } else {
      Serial.println("# ERR,reason=plan_subcommand");
    }
  } else if (verb == "ABORT") {
    abortRequested = true;
    Serial.println("# ABORT");
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
  digitalWrite(LED_PIN, LOW);

  Serial.begin(115200);
  const uint32_t t0 = millis();
  while (!Serial && millis() - t0 < 3000) delay(10);
  delay(200);

  Serial.println("# BOOT,fw=stage_identification");
  Serial.printf("# PINS,step=%u,dir=%u,uart_rx=%u,uart_tx=%u,enc_a=%d,enc_b=%d,enc_z=%d\n",
                STEP_PIN, DIR_PIN, UART_RX_PIN, UART_TX_PIN,
                (int)ENC_A, (int)ENC_B, (int)ENC_Z);
  const bool enc = setupEncoder();
  if (!enc) Serial.println("# ERR,reason=pcnt_init_failed");
  driverReady = configureDriver();
  Serial.printf("# READY,driver=%u,encoder=%u,mres=%u,current_ma=%u,steps_per_mm=%u\n",
                driverReady, enc, mres, currentMa, FULL_STEPS_PER_MM);
  Serial.println("# IDLE,upload a plan then PLAN RUN; nothing moves before that");
}

void loop() {
  while (Serial.available()) {
    const char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      if (inBuf.length()) { String l = inBuf; inBuf = ""; handleCommand(l); }
    } else if (inBuf.length() < 160) {
      inBuf += c;
    }
  }
  static uint32_t last = 0;
  if (!running && millis() - last > 1000) { last = millis(); }
}

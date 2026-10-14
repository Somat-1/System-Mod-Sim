/* Rotor-side motor identification: ESP32-S3 + TMC2209 + shaft encoder.
 * Native microsteps only. SpreadCycle is forced, MicroPlyer is disabled, and
 * the plan cannot run unless register readback confirms both conditions and
 * DRV_STATUS confirms that StealthChop is inactive.
 */
#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/pulse_cnt.h"

// Stepper Isolated bench wiring. EN is TMC ENN (active low): held HIGH until
// the driver mode is verified, then LOW.
constexpr uint8_t STEP_PIN = 7, DIR_PIN = 8, EN_PIN = 9;
constexpr uint8_t UART_RX_PIN = 18, UART_TX_PIN = 17, LED_PIN = 48;
constexpr gpio_num_t ENC_A = GPIO_NUM_4, ENC_B = GPIO_NUM_5, ENC_Z = GPIO_NUM_6;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;   // R100 on the fitted module
constexpr uint32_t DRIVER_BAUD = 115200, STEP_HIGH_US = 3, DIR_SETUP_US = 20;
constexpr uint32_t MAX_PULSE_HZ = 20000;
constexpr uint16_t DEFAULT_CURRENT_MA = 360, DEFAULT_MRES = 16;
constexpr int64_t U256 = 256;
constexpr size_t MAX_ACTIONS = 128;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);
pcnt_unit_handle_t pcnt = nullptr;

enum ActionKind : uint8_t {
  HOLD = 0, MOVE = 1, STAIR = 2, STEP_RESPONSE = 3, VELOCITY = 4, REVERSAL = 5
};

struct PlanAction {
  uint8_t kind;
  int32_t microsteps;
  uint32_t speedMilliFs;
  uint32_t dwellMs;
  uint16_t repeats;
  uint16_t currentMa;
  uint16_t mres;
  uint16_t sampleHz;
  char tag[24];
};

PlanAction plan[MAX_ACTIONS];
size_t planCount = 0, actionIndex = 0;
bool driverReady = false, encoderReady = false, running = false, abortRequested = false;
bool lastDirection = true, directionKnown = false;
int64_t commandU256 = 0;
uint16_t currentMa = DEFAULT_CURRENT_MA, mres = DEFAULT_MRES;
int64_t stepU256 = U256 / DEFAULT_MRES;
volatile int32_t indexCount = 0;
uint32_t sampleSeq = 0, samplePeriodUs = 500, nextSampleUs = 0;
String inputBuffer;

void IRAM_ATTR onIndex() { indexCount++; }

int32_t encoderCount() {
  int value = 0;
  if (pcnt) pcnt_unit_get_count(pcnt, &value);
  return value;
}

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
  // accum_count only carries the count across the limits if the limits are
  // watch points; without these the count silently restarts at 0 every
  // 30000 counts (1.5 shaft revolutions on the 20000 CPR encoder).
  pcnt_unit_add_watch_point(pcnt, 30000);
  pcnt_unit_add_watch_point(pcnt, -30000);
  pcnt_unit_enable(pcnt); pcnt_unit_clear_count(pcnt); pcnt_unit_start(pcnt);
  pinMode(ENC_Z, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(ENC_Z), onIndex, RISING);
  return true;
}

template <typename Reader>
bool readTwice(Reader read, uint32_t &value) {
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    const uint32_t a = read(), b = read();
    if (a == b && a != 0xFFFFFFFFUL) { value = a; return true; }
    delay(5);
  }
  return false;
}

bool verifyDriverMode(bool report) {
  uint32_t gconf = 0, chopconf = 0, status = 0, iholdIrun = 0;
  if (!readTwice([] { return driver.GCONF(); }, gconf) ||
      !readTwice([] { return driver.CHOPCONF(); }, chopconf) ||
      !readTwice([] { return driver.DRV_STATUS(); }, status) ||
      !readTwice([] { return driver.IHOLD_IRUN(); }, iholdIrun)) return false;
  const bool spread = (gconf & (1UL << 2)) != 0;
  const bool intpol = (chopconf & (1UL << 28)) != 0;
  const bool stealth = (status & (1UL << 30)) != 0;
  const uint8_t ihold = iholdIrun & 0x1F, irun = (iholdIrun >> 8) & 0x1F;
  if (report) {
    Serial.printf("# DRIVER,gconf=0x%08lX,chopconf=0x%08lX,"
                  "drv_status=0x%08lX,spreadcycle=%u,intpol=%u,stealth=%u,"
                  "ihold=%u,irun=%u,coolstep=0\n",
                  (unsigned long)gconf, (unsigned long)chopconf,
                  (unsigned long)status, spread, intpol, stealth, ihold, irun);
  }
  return spread && !intpol && !stealth && ihold == irun;
}

bool applyMres(uint16_t value) {
  if (value < 1 || value > 256 || (value & (value - 1))) return false;
  uint8_t code = 8;
  for (uint16_t v = value; v > 1; v >>= 1) --code;
  constexpr uint32_t MASK = 0x0F000000UL;
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    uint32_t chop = 0;
    if (!readTwice([] { return driver.CHOPCONF(); }, chop)) continue;
    driver.CHOPCONF((chop & ~MASK) | ((uint32_t)code << 24));
    uint32_t back = 0;
    if (readTwice([] { return driver.CHOPCONF(); }, back) &&
        ((back & MASK) >> 24) == code && !(back & (1UL << 28))) {
      mres = value; stepU256 = U256 / value; return true;
    }
    delay(5);
  }
  return false;
}

void setCurrent(uint16_t milliamps) {
  driver.rms_current(milliamps);
  driver.ihold(driver.irun());  // no standstill-current reduction
  currentMa = milliamps;
}

bool configureDriver() {
  digitalWrite(EN_PIN, HIGH);
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(100);
  // Chip ID first: version 0x21 is a TMC2209 answering; 0x00/0xFF means the
  // UART gets no reply (VM off, PDN_UART wiring, or no common ground).
  const uint32_t ioin = driver.IOIN();
  Serial.printf("# DRIVER_ID,ioin=0x%08lX,version=0x%02X,expected=0x21\n",
                (unsigned long)ioin, (unsigned)((ioin >> 24) & 0xFF));
  driver.begin();
  driver.toff(5); driver.blank_time(24);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  setCurrent(DEFAULT_CURRENT_MA);
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);
  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);       // no velocity switchover to StealthChop
  driver.intpol(false);     // MicroPlyer explicitly disabled
  driver.semin(0);          // CoolStep disabled
  driver.TCOOLTHRS(0);      // StallGuard/CoolStep threshold disabled
  driver.pwm_autoscale(false);
  uint8_t connection = 2;
  for (uint8_t attempt = 0; attempt < 5 && connection; ++attempt) {
    connection = driver.test_connection();
    if (connection) delay(50);
  }
  if (connection) return false;
  if (!applyMres(DEFAULT_MRES)) return false;
  const bool ok = verifyDriverMode(true);
  if (ok) digitalWrite(EN_PIN, LOW);
  return ok;
}

void emitSample(int eventCode = 0) {
  Serial.printf("D,%lu,%llu,%u,%d,%lld,%ld,%ld,%u,%u\n",
                (unsigned long)sampleSeq++,
                (unsigned long long)esp_timer_get_time(),
                (unsigned)actionIndex, eventCode, (long long)commandU256,
                (long)encoderCount(), (long)indexCount, currentMa, mres);
}

bool pollAbort() {
  while (Serial.available()) {
    const char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      inputBuffer.trim();
      if (inputBuffer.equalsIgnoreCase("ABORT")) {
        abortRequested = true;
        Serial.println("# ABORT_ACK");
      }
      inputBuffer = "";
    } else if (inputBuffer.length() < 180) inputBuffer += c;
  }
  return abortRequested;
}

void serviceSample() {
  const uint32_t now = micros();
  if ((int32_t)(now - nextSampleUs) >= 0) {
    emitSample(0);
    nextSampleUs = now + samplePeriodUs;
  }
}

bool waitUs(uint32_t durationUs) {
  const uint32_t end = micros() + durationUs;
  while ((int32_t)(micros() - end) < 0) {
    serviceSample();
    if (pollAbort()) return false;
    delayMicroseconds(20);
  }
  return true;
}

bool dwellMs(uint32_t durationMs) {
  const uint64_t end = esp_timer_get_time() + (uint64_t)durationMs * 1000ULL;
  while (esp_timer_get_time() < end) {
    serviceSample();
    if (pollAbort()) return false;
    delayMicroseconds(20);
  }
  return true;
}

bool pulseOne(bool positive, bool markEdge) {
  if (!directionKnown || positive != lastDirection) {
    digitalWrite(DIR_PIN, positive ? HIGH : LOW);
    delayMicroseconds(DIR_SETUP_US);
    lastDirection = positive; directionKnown = true;
  }
  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(STEP_HIGH_US);
  digitalWrite(STEP_PIN, LOW);
  commandU256 += positive ? stepU256 : -stepU256;
  if (markEdge) emitSample(positive ? 1 : -1);
  return !pollAbort();
}

bool pulseTrain(uint32_t count, bool positive, uint32_t periodUs,
                bool markEdges = false) {
  if (periodUs <= STEP_HIGH_US || periodUs < 1000000UL / MAX_PULSE_HZ)
    return false;
  uint32_t target = micros();
  for (uint32_t i = 0; i < count; ++i) {
    if (!pulseOne(positive, markEdges)) return false;
    target += periodUs;
    while ((int32_t)(micros() - target) < 0) {
      serviceSample();
      if (pollAbort()) return false;
      delayMicroseconds(10);
    }
  }
  return true;
}

uint32_t periodFor(const PlanAction &a) {
  const uint64_t pulseMilliHz = (uint64_t)a.speedMilliFs * a.mres;
  if (!pulseMilliHz) return 0;
  return (uint32_t)(1000000000ULL / pulseMilliHz);
}

bool executeAction(const PlanAction &a) {
  if (a.currentMa != currentMa) setCurrent(a.currentMa);
  if (a.mres != mres && !applyMres(a.mres)) {
    Serial.println("# ERR,reason=mres_write_failed"); return false;
  }
  if (!verifyDriverMode(false)) {
    Serial.println("# ERR,reason=driver_mode_changed"); return false;
  }
  samplePeriodUs = 1000000UL / a.sampleHz;
  nextSampleUs = micros();
  const uint32_t count = abs(a.microsteps);
  const bool positive = a.microsteps >= 0;
  const uint32_t periodUs = periodFor(a);

  if (a.kind == HOLD) return dwellMs(a.dwellMs);
  if (a.kind == MOVE || a.kind == VELOCITY) {
    if (!periodUs || !pulseTrain(count, positive, periodUs)) return false;
    return dwellMs(a.dwellMs);
  }
  if (a.kind == STAIR) {
    for (uint32_t i = 0; i < count; ++i) {
      if (!pulseOne(positive, true) || !dwellMs(a.dwellMs)) return false;
    }
    return true;
  }
  if (a.kind == STEP_RESPONSE) {
    for (uint16_t rep = 0; rep < a.repeats; ++rep) {
      for (uint32_t i = 0; i < count; ++i)
        if (!pulseOne(positive, true)) return false;
      if (!dwellMs(a.dwellMs)) return false;
      for (uint32_t i = 0; i < count; ++i)
        if (!pulseOne(!positive, true)) return false;
      if (!dwellMs(a.dwellMs)) return false;
    }
    return true;
  }
  if (a.kind == REVERSAL) {
    if (!periodUs) return false;
    for (uint16_t rep = 0; rep < a.repeats; ++rep) {
      if (!pulseTrain(count, true, periodUs, true)) return false;
      if (!pulseTrain(2 * count, false, periodUs, true)) return false;
      if (!pulseTrain(count, true, periodUs, true)) return false;
    }
    return dwellMs(a.dwellMs);
  }
  Serial.println("# ERR,reason=unknown_action");
  return false;
}

void runPlan() {
  if (!driverReady || !encoderReady || !planCount) {
    Serial.println("# ERR,reason=not_ready_or_empty"); return;
  }
  if (!verifyDriverMode(true)) {
    Serial.println("# ERR,reason=driver_mode_not_confirmed"); return;
  }
  running = true; abortRequested = false; sampleSeq = 0;
  digitalWrite(LED_PIN, HIGH);
  Serial.printf("# RUN_START,actions=%u,command_u256=%lld,encoder=%ld\n",
                (unsigned)planCount, (long long)commandU256, (long)encoderCount());
  bool ok = true;
  for (actionIndex = 0; actionIndex < planCount && !abortRequested; ++actionIndex) {
    const PlanAction &a = plan[actionIndex];
    Serial.printf("# ACTION_BEGIN,index=%u,tag=%s,kind=%u,current_ma=%u,mres=%u\n",
                  (unsigned)actionIndex, a.tag, a.kind, a.currentMa, a.mres);
    ok = executeAction(a);
    Serial.printf("# ACTION_END,index=%u,tag=%s,command_u256=%lld,encoder=%ld,ok=%u\n",
                  (unsigned)actionIndex, a.tag, (long long)commandU256,
                  (long)encoderCount(), ok);
    if (!ok) break;
  }
  digitalWrite(LED_PIN, LOW); running = false;
  const char *result = abortRequested ? "ABORTED" : (ok ? "COMPLETE" : "FAILED");
  Serial.printf("# RUN_%s,command_u256=%lld,encoder=%ld\n",
                result, (long long)commandU256, (long)encoderCount());
}

void handleCommand(String line) {
  line.trim();
  if (!line.length()) return;
  if (line.equalsIgnoreCase("PING")) {
    Serial.printf("# PONG,fw=motor_identification,ready=%u,actions=%u\n",
                  driverReady && encoderReady, (unsigned)planCount);
    return;
  }
  if (line.equalsIgnoreCase("DRIVER")) {
    if (!verifyDriverMode(true))
      Serial.println("# ERR,reason=driver_mode_not_confirmed");
    return;
  }
  if (line.equalsIgnoreCase("STATUS")) {
    Serial.printf("# STATUS,command_u256=%lld,encoder=%ld,index=%ld,current_ma=%u,"
                  "mres=%u,actions=%u,running=%u\n", (long long)commandU256,
                  (long)encoderCount(), (long)indexCount, currentMa, mres,
                  (unsigned)planCount, running);
    return;
  }
  if (line.equalsIgnoreCase("ZERO")) {
    if (running) { Serial.println("# ERR,reason=running"); return; }
    commandU256 = 0; indexCount = 0;
    if (pcnt) pcnt_unit_clear_count(pcnt);
    Serial.println("# ZERO_OK"); return;
  }
  if (line.equalsIgnoreCase("PLAN CLEAR")) {
    if (running) { Serial.println("# ERR,reason=running"); return; }
    planCount = 0; Serial.println("# PLAN_CLEARED"); return;
  }
  if (line.startsWith("PLAN ADD ")) {
    if (running || planCount >= MAX_ACTIONS) {
      Serial.println("# ERR,reason=running_or_plan_full"); return;
    }
    PlanAction a = {};
    char tag[24] = "";
    unsigned kind, repeats, current, resolution, sampleRate;
    long steps;
    unsigned long speedMilli, dwell;
    const int n = sscanf(line.c_str() + 9, "%u %ld %lu %lu %u %u %u %u %23s",
                         &kind, &steps, &speedMilli, &dwell, &repeats, &current,
                         &resolution, &sampleRate, tag);
    if (n != 9 || kind > REVERSAL || resolution < 1 || resolution > 256 ||
        (resolution & (resolution - 1)) || sampleRate < 10 || sampleRate > 5000 ||
        current < 1 || current > 2000 || repeats < 1) {
      Serial.println("# ERR,reason=bad_plan_fields"); return;
    }
    a.kind = kind; a.microsteps = steps; a.speedMilliFs = speedMilli;
    a.dwellMs = dwell; a.repeats = repeats; a.currentMa = current;
    a.mres = resolution; a.sampleHz = sampleRate;
    strncpy(a.tag, tag, sizeof(a.tag) - 1);
    plan[planCount++] = a;
    Serial.printf("# PLAN_ADD,index=%u,tag=%s\n", (unsigned)(planCount - 1), a.tag);
    return;
  }
  if (line.equalsIgnoreCase("REINIT")) {
    // Retry the driver setup without a reset, e.g. after VM was switched on.
    if (running) { Serial.println("# ERR,reason=running"); return; }
    driverReady = configureDriver();
    Serial.printf("# REINIT_DONE,ready=%u\n", driverReady && encoderReady);
    return;
  }
  if (line.equalsIgnoreCase("PLAN RUN")) { runPlan(); return; }
  if (line.equalsIgnoreCase("ABORT")) {
    abortRequested = true; Serial.println("# ABORT_ACK"); return;
  }
  Serial.println("# ERR,reason=unknown_command");
}

void setup() {
  pinMode(EN_PIN, OUTPUT); digitalWrite(EN_PIN, HIGH);
  pinMode(STEP_PIN, OUTPUT); pinMode(DIR_PIN, OUTPUT); pinMode(LED_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW); digitalWrite(DIR_PIN, HIGH); digitalWrite(LED_PIN, LOW);
  Serial.begin(115200);
  const uint32_t start = millis();
  while (!Serial && millis() - start < 3000) delay(10);
  delay(200);
  Serial.println("# BOOT,fw=motor_identification");
  Serial.printf("# PINS,step=%u,dir=%u,uart_rx=%u,uart_tx=%u,enc_a=%d,enc_b=%d,enc_z=%d\n",
                STEP_PIN, DIR_PIN, UART_RX_PIN, UART_TX_PIN,
                (int)ENC_A, (int)ENC_B, (int)ENC_Z);
  encoderReady = setupEncoder();
  if (!encoderReady) Serial.println("# ERR,reason=pcnt_init_failed");
  driverReady = configureDriver();
  if (!driverReady) Serial.println("# ERR,reason=driver_config_failed");
  Serial.printf("# READY,driver=%u,encoder=%u,current_ma=%u,mres=%u\n",
                driverReady, encoderReady, currentMa, mres);
  Serial.println("# IDLE,upload plan; no motion before PLAN RUN");
}

void loop() {
  while (!running && Serial.available()) {
    const char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      if (inputBuffer.length()) {
        String line = inputBuffer; inputBuffer = ""; handleCommand(line);
      }
    } else if (inputBuffer.length() < 200) inputBuffer += c;
  }
}

// Combined encoder monitor stream and conservative TMC2209 motor test.
//
// Encoder: A=GPIO4, B=GPIO5, Z=GPIO6 (hardware x4 quadrature decode).
// TMC2209: STEP=GPIO7, DIR=GPIO8, EN=GPIO9, TX=GPIO17 via 1k,
//          RX=GPIO18 direct, UART address 0.
//
// The first INFO request from encoder_monitor.py schedules one finite move:
// 3200 microsteps at 1/16 microstepping over 60 s. This is one revolution at
// 1 rpm for a 200-full-step/revolution motor. EN is disabled before and after.

#include <Arduino.h>
#include <TMCStepper.h>
#include "driver/gpio.h"
#include "driver/pulse_cnt.h"
#include "esp_timer.h"

// Encoder inputs.
static const gpio_num_t PIN_A = GPIO_NUM_4;
static const gpio_num_t PIN_B = GPIO_NUM_5;
static const gpio_num_t PIN_Z = GPIO_NUM_6;

// Motor driver pins and configuration.
constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.10F;  // R100 sense resistors fitted to this module.
constexpr uint16_t CURRENT_RMS_MA = 400;
constexpr uint16_t MICROSTEPS = 16;
constexpr uint32_t MOVE_STEPS = 3200;
constexpr uint32_t STEP_INTERVAL_US = 18750;  // 53.333 Hz; 60 s total.
constexpr uint32_t START_DELAY_US = 3000000;

// Encoder stream configuration.
constexpr bool USE_INTERNAL_PULLUPS = true;
constexpr uint32_t GLITCH_FILTER_NS = 1000;
constexpr int PCNT_LIMIT = 30000;
constexpr uint32_t MIN_RATE_HZ = 10;
constexpr uint32_t MAX_RATE_HZ = 5000;

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);
static pcnt_unit_handle_t pcnt_unit = nullptr;

static volatile uint32_t z_count = 0;
static volatile int32_t z_pos = 0;
static volatile bool z_latched = false;

static uint32_t rate_hz = 1000;
static uint64_t period_us = 1000;
static uint64_t next_sample_us = 0;
static uint32_t seq = 0;
static bool streaming = true;

static bool motor_ready = false;
static bool control_pins_verified = false;
static uint8_t control_pin_test_bits = 0;
static uint32_t control_ioin_low = 0;
static uint32_t control_ioin_step_high = 0;
static uint32_t control_ioin_dir_high = 0;
static bool move_scheduled = false;
static bool move_started = false;
static bool moving = false;
static uint32_t move_steps_done = 0;
static uint64_t move_start_us = 0;
static uint64_t next_step_us = 0;

static void IRAM_ATTR onIndex() {
  int count = 0;
  pcnt_unit_get_count(pcnt_unit, &count);
  z_pos = count;
  z_count++;
  z_latched = true;
}

static void setupPcnt() {
  pcnt_unit_config_t unit_cfg = {};
  unit_cfg.low_limit = -PCNT_LIMIT;
  unit_cfg.high_limit = PCNT_LIMIT;
  unit_cfg.flags.accum_count = 1;
  ESP_ERROR_CHECK(pcnt_new_unit(&unit_cfg, &pcnt_unit));

  pcnt_glitch_filter_config_t filter_cfg = {};
  filter_cfg.max_glitch_ns = GLITCH_FILTER_NS;
  ESP_ERROR_CHECK(pcnt_unit_set_glitch_filter(pcnt_unit, &filter_cfg));

  pcnt_chan_config_t cha_cfg = {};
  cha_cfg.edge_gpio_num = PIN_A;
  cha_cfg.level_gpio_num = PIN_B;
  pcnt_channel_handle_t cha = nullptr;
  ESP_ERROR_CHECK(pcnt_new_channel(pcnt_unit, &cha_cfg, &cha));

  pcnt_chan_config_t chb_cfg = {};
  chb_cfg.edge_gpio_num = PIN_B;
  chb_cfg.level_gpio_num = PIN_A;
  pcnt_channel_handle_t chb = nullptr;
  ESP_ERROR_CHECK(pcnt_new_channel(pcnt_unit, &chb_cfg, &chb));

  ESP_ERROR_CHECK(pcnt_channel_set_edge_action(
      cha, PCNT_CHANNEL_EDGE_ACTION_DECREASE, PCNT_CHANNEL_EDGE_ACTION_INCREASE));
  ESP_ERROR_CHECK(pcnt_channel_set_level_action(
      cha, PCNT_CHANNEL_LEVEL_ACTION_KEEP, PCNT_CHANNEL_LEVEL_ACTION_INVERSE));
  ESP_ERROR_CHECK(pcnt_channel_set_edge_action(
      chb, PCNT_CHANNEL_EDGE_ACTION_INCREASE, PCNT_CHANNEL_EDGE_ACTION_DECREASE));
  ESP_ERROR_CHECK(pcnt_channel_set_level_action(
      chb, PCNT_CHANNEL_LEVEL_ACTION_KEEP, PCNT_CHANNEL_LEVEL_ACTION_INVERSE));

  ESP_ERROR_CHECK(pcnt_unit_add_watch_point(pcnt_unit, PCNT_LIMIT));
  ESP_ERROR_CHECK(pcnt_unit_add_watch_point(pcnt_unit, -PCNT_LIMIT));
  ESP_ERROR_CHECK(pcnt_unit_enable(pcnt_unit));
  ESP_ERROR_CHECK(pcnt_unit_clear_count(pcnt_unit));
  ESP_ERROR_CHECK(pcnt_unit_start(pcnt_unit));

  if (USE_INTERNAL_PULLUPS) {
    gpio_pullup_en(PIN_A);
    gpio_pullup_en(PIN_B);
  }
}

static bool configureMotor() {
  digitalWrite(EN_PIN, HIGH);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);

  DriverSerial.begin(115200, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();
  driver.toff(5);
  driver.I_scale_analog(false);
  driver.internal_Rsense(false);
  driver.mstep_reg_select(true);
  driver.pdn_disable(true);
  // SpreadCycle is intentionally used for first-motion bring-up: it provides
  // deterministic current regulation and more starting torque than an
  // uncalibrated StealthChop setup at very low speed.
  driver.en_spreadCycle(true);
  driver.intpol(false);
  driver.microsteps(MICROSTEPS);
  driver.rms_current(CURRENT_RMS_MA);
  driver.ihold(driver.irun());
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);
  driver.pwm_autoscale(true);
  driver.pwm_autograd(true);
  delay(50);

  const uint32_t ioin = driver.IOIN();
  const uint8_t version = static_cast<uint8_t>(ioin >> 24);
  const uint8_t mres = driver.mres();
  const uint16_t current = driver.rms_current();
  const uint32_t drv_status = driver.DRV_STATUS();

  // Verify the signals at the TMC's own IOIN register while ENN remains HIGH.
  // The STEP rising edge is harmless here because the power stage is disabled.
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);
  delay(100);
  const uint32_t ioin_both_low = driver.IOIN();
  const bool step_low_seen = (ioin_both_low & (1UL << 7)) == 0;
  const bool dir_low_seen = (ioin_both_low & (1UL << 9)) == 0;
  digitalWrite(STEP_PIN, HIGH);
  delay(100);
  const uint32_t ioin_step_high = driver.IOIN();
  const bool step_high_seen = (ioin_step_high & (1UL << 7)) != 0;
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, HIGH);
  delay(100);
  const uint32_t ioin_dir_high = driver.IOIN();
  const bool dir_high_seen = (ioin_dir_high & (1UL << 9)) != 0;
  control_ioin_low = ioin_both_low;
  control_ioin_step_high = ioin_step_high;
  control_ioin_dir_high = ioin_dir_high;
  digitalWrite(DIR_PIN, LOW);
  const bool control_pins_ok = step_low_seen && step_high_seen &&
                               dir_low_seen && dir_high_seen;
  control_pins_verified = control_pins_ok;
  control_pin_test_bits = (step_low_seen ? 1u : 0u) |
                          (step_high_seen ? 2u : 0u) |
                          (dir_low_seen ? 4u : 0u) |
                          (dir_high_seen ? 8u : 0u);

  Serial.printf("#MOTOR_CONFIG,version=0x%02X,current_requested_mA=%u,"
                "current_readback_mA=%u,microsteps=%u,mres_code=%u,"
                "drv_status=0x%08lX,en=disabled\n",
                version, CURRENT_RMS_MA, current, MICROSTEPS, mres,
                static_cast<unsigned long>(drv_status));
  Serial.printf("#CONTROL_PIN_TEST,step_low=%u,step_high=%u,dir_low=%u,"
                "dir_high=%u,result=%s,en=disabled,ioin_low=0x%08lX,"
                "ioin_step_high=0x%08lX,ioin_dir_high=0x%08lX\n",
                step_low_seen ? 1 : 0, step_high_seen ? 1 : 0,
                dir_low_seen ? 1 : 0, dir_high_seen ? 1 : 0,
                control_pins_ok ? "PASS" : "FAIL",
                static_cast<unsigned long>(ioin_both_low),
                static_cast<unsigned long>(ioin_step_high),
                static_cast<unsigned long>(ioin_dir_high));

  const bool uart_ok = version == 0x21 &&
                       ioin != 0 && ioin != 0xFFFFFFFFUL;
  const bool config_ok = mres == 4 && current >= 350 && current <= 450;
  if (!uart_ok || !config_ok || !control_pins_ok) {
    Serial.printf("#MOTOR_REFUSED,uart_ok=%u,config_ok=%u,control_pins_ok=%u\n",
                  uart_ok ? 1 : 0, config_ok ? 1 : 0,
                  control_pins_ok ? 1 : 0);
    return false;
  }
  return true;
}

static void printInfo() {
  Serial.printf("#INFO,fw=stepper_encoder_slow_test,pin_a=%d,pin_b=%d,pin_z=%d,"
                "rate_hz=%lu,pullups=%d,glitch_ns=%lu,decode=x4,streaming=%d,"
                "motor_ready=%d,control_pins_ok=%d,control_pin_bits=0x%X,"
                "ioin_low=0x%08lX,ioin_step_high=0x%08lX,ioin_dir_high=0x%08lX,"
                "current_rms_mA=%u,microsteps=%u,move_steps=%lu,"
                "step_interval_us=%lu\n",
                static_cast<int>(PIN_A), static_cast<int>(PIN_B),
                static_cast<int>(PIN_Z), static_cast<unsigned long>(rate_hz),
                USE_INTERNAL_PULLUPS ? 1 : 0,
                static_cast<unsigned long>(GLITCH_FILTER_NS), streaming ? 1 : 0,
                motor_ready ? 1 : 0, control_pins_verified ? 1 : 0,
                control_pin_test_bits,
                static_cast<unsigned long>(control_ioin_low),
                static_cast<unsigned long>(control_ioin_step_high),
                static_cast<unsigned long>(control_ioin_dir_high),
                CURRENT_RMS_MA, MICROSTEPS,
                static_cast<unsigned long>(MOVE_STEPS),
                static_cast<unsigned long>(STEP_INTERVAL_US));
}

static void setRate(long hz) {
  if (hz < static_cast<long>(MIN_RATE_HZ)) hz = MIN_RATE_HZ;
  if (hz > static_cast<long>(MAX_RATE_HZ)) hz = MAX_RATE_HZ;
  rate_hz = static_cast<uint32_t>(hz);
  period_us = 1000000ULL / rate_hz;
  next_sample_us = esp_timer_get_time();
  Serial.printf("#ACK,rate_hz=%lu\n", static_cast<unsigned long>(rate_hz));
}

static void zeroCounts() {
  pcnt_unit_clear_count(pcnt_unit);
  noInterrupts();
  z_count = 0;
  z_pos = 0;
  z_latched = false;
  interrupts();
  Serial.println("#ACK,zero");
}

static void stopMotor(const char *reason) {
  digitalWrite(EN_PIN, HIGH);
  digitalWrite(STEP_PIN, LOW);
  if (moving || move_scheduled) {
    Serial.printf("#MOTOR_STOP,reason=%s,steps=%lu\n", reason,
                  static_cast<unsigned long>(move_steps_done));
  }
  moving = false;
  move_scheduled = false;
}

static void scheduleMove() {
  if (!motor_ready || move_started || move_scheduled) return;
  move_scheduled = true;
  move_start_us = esp_timer_get_time() + START_DELAY_US;
  Serial.println("#MOTOR_SCHEDULED,start_in_s=3,duration_s=60,direction=HIGH");
}

static void handleCommand(const String &cmd) {
  if (cmd.length() == 0) return;
  const char command = toupper(cmd[0]);
  if (command == 'R') {
    setRate(cmd.substring(1).toInt());
  } else if (command == 'Z') {
    zeroCounts();
  } else if (command == 'S') {
    streaming = cmd.length() > 1 && cmd[1] == '1';
    next_sample_us = esp_timer_get_time();
    Serial.printf("#ACK,streaming=%d\n", streaming ? 1 : 0);
    if (!streaming) stopMotor("stream_stopped");
  } else if (command == 'I') {
    printInfo();
    scheduleMove();
  } else if (command == 'X') {
    stopMotor("abort_command");
  } else {
    Serial.printf("#ERR,unknown_command=%s\n", cmd.c_str());
  }
}

static void pollSerial() {
  static String buffer;
  while (Serial.available()) {
    const char ch = static_cast<char>(Serial.read());
    if (ch == '\n' || ch == '\r') {
      buffer.trim();
      handleCommand(buffer);
      buffer = "";
    } else if (buffer.length() < 32) {
      buffer += ch;
    }
  }
}

static void updateMotor(uint64_t now) {
  if (move_scheduled && !moving && static_cast<int64_t>(now - move_start_us) >= 0) {
    move_scheduled = false;
    move_started = true;
    moving = true;
    move_steps_done = 0;
    digitalWrite(DIR_PIN, HIGH);
    digitalWrite(STEP_PIN, LOW);
    digitalWrite(EN_PIN, LOW);
    delay(20);
    const uint32_t enabled_ioin = driver.IOIN();
    const uint32_t enabled_status = driver.DRV_STATUS();
    Serial.printf("#MOTOR_ENABLED,ioin=0x%08lX,enn=%u,gstat=0x%08lX,"
                  "drv_status=0x%08lX,ihold_irun=0x%08lX,chopconf=0x%08lX\n",
                  static_cast<unsigned long>(enabled_ioin), driver.enn() ? 1 : 0,
                  static_cast<unsigned long>(driver.GSTAT()),
                  static_cast<unsigned long>(enabled_status),
                  static_cast<unsigned long>(driver.IHOLD_IRUN()),
                  static_cast<unsigned long>(driver.CHOPCONF()));
    next_step_us = esp_timer_get_time();
    Serial.println("#MOTOR_START,current_rms_mA=400,microsteps=16,target_steps=3200");
  }

  if (!moving || static_cast<int64_t>(now - next_step_us) < 0) return;

  digitalWrite(STEP_PIN, HIGH);
  delayMicroseconds(4);
  digitalWrite(STEP_PIN, LOW);
  move_steps_done++;
  next_step_us += STEP_INTERVAL_US;

  if (move_steps_done == 64) {
    const uint32_t status = driver.DRV_STATUS();
    Serial.printf("#MOTOR_RUNNING_STATUS,drv_status=0x%08lX,ot=%u,otpw=%u,"
                  "s2ga=%u,s2gb=%u,s2vsa=%u,s2vsb=%u,ola=%u,olb=%u,cs=%u\n",
                  static_cast<unsigned long>(status), driver.ot() ? 1 : 0,
                  driver.otpw() ? 1 : 0, driver.s2ga() ? 1 : 0,
                  driver.s2gb() ? 1 : 0, driver.s2vsa() ? 1 : 0,
                  driver.s2vsb() ? 1 : 0, driver.ola() ? 1 : 0,
                  driver.olb() ? 1 : 0, driver.cs_actual());
  }

  if (move_steps_done >= MOVE_STEPS) {
    digitalWrite(EN_PIN, HIGH);
    moving = false;
    Serial.printf("#MOTOR_COMPLETE,steps=%lu,en=disabled\n",
                  static_cast<unsigned long>(move_steps_done));
  }
}

static void streamEncoder(uint64_t now) {
  if (!streaming || static_cast<int64_t>(now - next_sample_us) < 0) return;

  next_sample_us += period_us;
  if (static_cast<int64_t>(now - next_sample_us) >
      static_cast<int64_t>(10 * period_us)) {
    next_sample_us = now + period_us;
  }

  int count = 0;
  pcnt_unit_get_count(pcnt_unit, &count);
  uint32_t levels = (gpio_get_level(PIN_A) ? 1u : 0u) |
                    (gpio_get_level(PIN_B) ? 2u : 0u) |
                    (gpio_get_level(PIN_Z) ? 4u : 0u);

  noInterrupts();
  if (z_latched) levels |= 8u;
  z_latched = false;
  const uint32_t zc = z_count;
  const int32_t zp = z_pos;
  interrupts();

  char line[96];
  const int length = snprintf(
      line, sizeof(line), "D,%lu,%llu,%d,%lu,%lu,%ld\n",
      static_cast<unsigned long>(seq++),
      static_cast<unsigned long long>(now), count,
      static_cast<unsigned long>(levels), static_cast<unsigned long>(zc),
      static_cast<long>(zp));
  Serial.write(reinterpret_cast<const uint8_t *>(line), length);
}

void setup() {
  pinMode(EN_PIN, OUTPUT);
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);

  Serial.begin(921600);
  delay(200);

  setupPcnt();
  pinMode(PIN_Z, USE_INTERNAL_PULLUPS ? INPUT_PULLUP : INPUT);
  attachInterrupt(digitalPinToInterrupt(PIN_Z), onIndex, RISING);

  motor_ready = configureMotor();
  setRate(rate_hz);
  printInfo();
  Serial.println(motor_ready
                     ? "#READY,waiting_for_encoder_monitor"
                     : "#NOT_READY,motor_will_remain_disabled");
}

void loop() {
  pollSerial();
  const uint64_t now = esp_timer_get_time();
  updateMotor(now);
  streamEncoder(now);
}

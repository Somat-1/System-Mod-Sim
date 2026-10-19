// ESP32-S3 quadrature encoder streamer (Motor Identification v1, 2026-10-02)
//
// A  -> GPIO4, B -> GPIO5, Z (index) -> GPIO6.
//
// A/B are decoded in hardware by the PCNT peripheral in x4 mode, so counts are
// never missed regardless of loop timing. Z is caught by a GPIO interrupt that
// latches the PCNT count at the index edge.
//
// Output (USB CDC, one line per sample, default 1000 Hz):
//   D,<seq>,<t_us>,<count>,<levels>,<z_count>,<z_pos>
//     levels  bit0 = A, bit1 = B, bit2 = Z (instantaneous),
//             bit3 = Z edge seen since previous sample (latched)
//     z_count number of Z rising edges since boot / last zero
//     z_pos   PCNT count captured at the most recent Z edge
// Lines starting with '#' are info/acks.
//
// Commands (newline-terminated):
//   R<hz>   set sample rate, 10..5000 Hz      e.g. R2000
//   Z       zero count and index statistics
//   S0/S1   stop / start streaming
//   I       print configuration

#include <Arduino.h>
#include "driver/gpio.h"
#include "driver/pulse_cnt.h"
#include "esp_timer.h"

static const gpio_num_t PIN_A = GPIO_NUM_4;
static const gpio_num_t PIN_B = GPIO_NUM_5;
static const gpio_num_t PIN_Z = GPIO_NUM_6;

// Enable the ESP32 internal pull-ups (~45 kOhm). Needed for open-collector
// encoders without external pull-ups; harmless for push-pull outputs.
static const bool USE_INTERNAL_PULLUPS = true;

// PCNT glitch filter: pulses shorter than this are ignored (max ~12 us).
static const uint32_t GLITCH_FILTER_NS = 1000;

static const int PCNT_LIMIT = 30000;  // hardware wraps here; driver accumulates
static const uint32_t MIN_RATE_HZ = 10;
static const uint32_t MAX_RATE_HZ = 5000;

static pcnt_unit_handle_t pcnt_unit = nullptr;

static volatile uint32_t z_count = 0;
static volatile int32_t z_pos = 0;
static volatile bool z_latched = false;

static uint32_t rate_hz = 1000;
static uint64_t period_us = 1000;
static uint64_t next_sample_us = 0;
static uint32_t seq = 0;
static bool streaming = true;

static void IRAM_ATTR onIndex() {
  int c = 0;
  pcnt_unit_get_count(pcnt_unit, &c);
  z_pos = c;
  z_count = z_count + 1;
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

  // Standard x4 quadrature decoding (ESP-IDF rotary encoder example).
  ESP_ERROR_CHECK(pcnt_channel_set_edge_action(cha, PCNT_CHANNEL_EDGE_ACTION_DECREASE, PCNT_CHANNEL_EDGE_ACTION_INCREASE));
  ESP_ERROR_CHECK(pcnt_channel_set_level_action(cha, PCNT_CHANNEL_LEVEL_ACTION_KEEP, PCNT_CHANNEL_LEVEL_ACTION_INVERSE));
  ESP_ERROR_CHECK(pcnt_channel_set_edge_action(chb, PCNT_CHANNEL_EDGE_ACTION_INCREASE, PCNT_CHANNEL_EDGE_ACTION_DECREASE));
  ESP_ERROR_CHECK(pcnt_channel_set_level_action(chb, PCNT_CHANNEL_LEVEL_ACTION_KEEP, PCNT_CHANNEL_LEVEL_ACTION_INVERSE));

  // Watch points at the limits let the driver extend the count beyond 16 bits.
  ESP_ERROR_CHECK(pcnt_unit_add_watch_point(pcnt_unit, PCNT_LIMIT));
  ESP_ERROR_CHECK(pcnt_unit_add_watch_point(pcnt_unit, -PCNT_LIMIT));

  ESP_ERROR_CHECK(pcnt_unit_enable(pcnt_unit));
  ESP_ERROR_CHECK(pcnt_unit_clear_count(pcnt_unit));
  ESP_ERROR_CHECK(pcnt_unit_start(pcnt_unit));

  // PCNT routes A/B through the GPIO matrix with input enabled; set pulls here.
  if (USE_INTERNAL_PULLUPS) {
    gpio_pullup_en(PIN_A);
    gpio_pullup_en(PIN_B);
  }
}

static void printInfo() {
  Serial.printf("#INFO,fw=esp32s3_encoder_stream_v1,pin_a=%d,pin_b=%d,pin_z=%d,rate_hz=%lu,"
                "pullups=%d,glitch_ns=%lu,decode=x4,streaming=%d\n",
                (int)PIN_A, (int)PIN_B, (int)PIN_Z, (unsigned long)rate_hz,
                USE_INTERNAL_PULLUPS ? 1 : 0, (unsigned long)GLITCH_FILTER_NS, streaming ? 1 : 0);
}

static void setRate(long hz) {
  if (hz < (long)MIN_RATE_HZ) hz = MIN_RATE_HZ;
  if (hz > (long)MAX_RATE_HZ) hz = MAX_RATE_HZ;
  rate_hz = (uint32_t)hz;
  period_us = 1000000ULL / rate_hz;
  next_sample_us = esp_timer_get_time();
  Serial.printf("#ACK,rate_hz=%lu\n", (unsigned long)rate_hz);
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

static void handleCommand(const String &cmd) {
  if (cmd.length() == 0) return;
  char c = toupper(cmd[0]);
  if (c == 'R') {
    setRate(cmd.substring(1).toInt());
  } else if (c == 'Z') {
    zeroCounts();
  } else if (c == 'S') {
    streaming = cmd.length() > 1 && cmd[1] == '1';
    next_sample_us = esp_timer_get_time();
    Serial.printf("#ACK,streaming=%d\n", streaming ? 1 : 0);
  } else if (c == 'I') {
    printInfo();
  } else {
    Serial.printf("#ERR,unknown_command=%s\n", cmd.c_str());
  }
}

static void pollSerial() {
  static String buf;
  while (Serial.available()) {
    char ch = (char)Serial.read();
    if (ch == '\n' || ch == '\r') {
      buf.trim();
      handleCommand(buf);
      buf = "";
    } else if (buf.length() < 32) {
      buf += ch;
    }
  }
}

void setup() {
  Serial.begin(921600);  // ignored by USB CDC, kept for clarity
  delay(200);

  setupPcnt();

  pinMode(PIN_Z, USE_INTERNAL_PULLUPS ? INPUT_PULLUP : INPUT);
  attachInterrupt(digitalPinToInterrupt(PIN_Z), onIndex, RISING);

  setRate(rate_hz);
  printInfo();
}

void loop() {
  pollSerial();
  if (!streaming) return;

  uint64_t now = esp_timer_get_time();
  if ((int64_t)(now - next_sample_us) < 0) return;

  // If we fell far behind (e.g. host stalled the USB link), resync instead of bursting.
  next_sample_us += period_us;
  if ((int64_t)(now - next_sample_us) > (int64_t)(10 * period_us)) next_sample_us = now + period_us;

  int count = 0;
  pcnt_unit_get_count(pcnt_unit, &count);
  uint32_t levels = (gpio_get_level(PIN_A) ? 1u : 0u) | (gpio_get_level(PIN_B) ? 2u : 0u) |
                    (gpio_get_level(PIN_Z) ? 4u : 0u);

  noInterrupts();
  if (z_latched) levels |= 8u;
  z_latched = false;
  uint32_t zc = z_count;
  int32_t zp = z_pos;
  interrupts();

  char line[96];
  int n = snprintf(line, sizeof(line), "D,%lu,%llu,%d,%lu,%lu,%ld\n", (unsigned long)seq++,
                   (unsigned long long)now, count, (unsigned long)levels, (unsigned long)zc, (long)zp);
  Serial.write((const uint8_t *)line, n);
}

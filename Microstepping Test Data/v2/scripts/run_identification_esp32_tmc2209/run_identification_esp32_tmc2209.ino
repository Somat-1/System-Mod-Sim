/*
 * Rev 2 microstepping identification runner: ESP32-S3 + TMC2209.
 *
 * Framework: Arduino-ESP32 3.x (ESP-IDF 5.x)
 * Libraries: TMCStepper
 * Timing: ESP-IDF RMT TX driver
 *
 * This runner owns the timing-critical portions of the campaign:
 *   Block 0, conditioning, A1, A2, B, E, final Block 0.
 * Blocks C and D are executed by run_identification_dedicated_controller.py.
 *
 * Safety model:
 *   - EN is held disabled until an explicit RUN command is received.
 *   - Every requested move is checked against the known campaign envelope.
 *   - ABORT is sampled between RMT batches and during every dwell.
 *   - Trigger is forced low and EN high on every normal or aborted exit.
 *
 * Serial commands at 115200 baud:
 *   CHECK  validate and print the campaign without enabling the motor
 *   RUN    execute all 12 MRES/current combinations
 *   ABORT  stop at the next safe batch/dwell boundary
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <driver/rmt_tx.h>
#include <esp_err.h>
#include <esp_timer.h>

namespace rig {

constexpr uint8_t STEP_PIN = 1;
constexpr uint8_t DIR_PIN = 2;
constexpr uint8_t EN_PIN = 5;
constexpr uint8_t UART_TX_PIN = 4;
constexpr uint8_t UART_RX_PIN = 6;
constexpr uint8_t TRIG_OUT_PIN = 7;
constexpr uint8_t TRIG_ECHO_PIN = 15;

constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.03F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;
constexpr uint16_t MOTOR_FULL_STEPS_PER_REV = 200;

constexpr uint32_t RMT_RESOLUTION_HZ = 1000000;
constexpr uint32_t RMT_HIGH_TICKS = 2;
constexpr size_t RMT_MAX_SYMBOLS = 64;
constexpr uint32_t DIR_SETUP_US = 5;

constexpr float BURST_FULL_STEPS_S = 250.0F;
constexpr float CONDITIONING_FULL_STEPS_S = 150.0F;
constexpr uint32_t DWELL_LADDER_MS = 400;
constexpr uint32_t LOOP_DWELL_MS = 300;
constexpr uint32_t DOUBLET_DWELL_MS = 1000;
constexpr uint8_t LADDER_REPEATS = 25;
constexpr uint8_t LOOP_REPEATS = 10;
constexpr uint8_t DOUBLET_REPEATS = 20;

// Largest planned excursion: A1, N=32, 25 positive moves at MRES 1/1.
// Units are 1/16 full steps: 25 * 32 * 16 = 12800.
constexpr int64_t MAX_ABS_POSITION_U16 = 12800;

constexpr bool USE_TRIGGER_ECHO = true;

struct CurrentLevel {
  const char *name;
  uint16_t setRmsMa;
  uint16_t measuredRmsMa;
};

constexpr CurrentLevel CURRENT_LEVELS[] = {
    {"I_lo", 360, 355},
    {"I_mid", 600, 556},
    {"I_hi", 750, 715},
};

constexpr uint16_t MRES_VALUES[] = {16, 4, 2, 1};
constexpr int16_t N_VALUES[] = {1, 2, 4, 8, 16, 32};
constexpr int16_t REFERENCE_MOVES[] = {
    16, -16, 4, -4, 1, -1, -16, 16, -4, 4, -1, 1};
constexpr int16_t NEST_DESC[] = {
    32, -16, 8, -4, 2, -1, 1, -2, 4, -8, 16, -32};
constexpr int16_t NEST_ASYM[] = {8, -3, 2, -5, 6, -8, 4, -4};
constexpr int16_t NEST_MINOR[] = {
    64, -16, 2, -2, -16, 2, -2, -16, 2, -2, -16};

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);

rmt_channel_handle_t rmtChannel = nullptr;
rmt_encoder_handle_t rmtEncoder = nullptr;
rmt_symbol_word_t rmtSymbols[RMT_MAX_SYMBOLS];

volatile bool aborted = false;
bool driverConfigured = false;
int64_t positionU16 = 0;
uint16_t activeMres = 16;
const CurrentLevel *activeCurrent = nullptr;
uint32_t runIndex = 0;
char activeBlock[48] = "IDLE";

void printCsvHeader() {
  Serial.println(
      "timestamp_us,event,run_index,current,mres,block,label,"
      "direction,pulses,pulse_rate_hz,period_ticks,position_u16,"
      "rmt_done_us,trigger_echo");
}

void logEvent(const char *event, const char *label, int direction = 0,
              uint32_t pulses = 0, float pulseRateHz = 0.0F,
              uint32_t periodTicks = 0, int64_t timestampUs = -1,
              int64_t rmtDoneUs = -1) {
  if (timestampUs < 0) {
    timestampUs = esp_timer_get_time();
  }
  const int echo = USE_TRIGGER_ECHO ? digitalRead(TRIG_ECHO_PIN) : -1;
  Serial.printf(
      "%lld,%s,%lu,%s,%u,%s,%s,%d,%lu,%.9g,%lu,%lld,%lld,%d\n",
      static_cast<long long>(timestampUs), event,
      static_cast<unsigned long>(runIndex),
      activeCurrent ? activeCurrent->name : "NA", activeMres, activeBlock,
      label ? label : "", direction, static_cast<unsigned long>(pulses),
      static_cast<double>(pulseRateHz),
      static_cast<unsigned long>(periodTicks),
      static_cast<long long>(positionU16),
      static_cast<long long>(rmtDoneUs), echo);
}

void pollAbort() {
  while (Serial.available() > 0) {
    String command = Serial.readStringUntil('\n');
    command.trim();
    command.toUpperCase();
    if (command == "ABORT") {
      aborted = true;
      logEvent("ABORT_REQUEST", "serial");
    }
  }
}

bool dwellMs(uint32_t durationMs, const char *label) {
  const int64_t startUs = esp_timer_get_time();
  const int64_t deadlineUs = startUs + static_cast<int64_t>(durationMs) * 1000;
  logEvent("DWELL_START", label, 0, 0, 0.0F, 0, startUs);
  while (!aborted && esp_timer_get_time() < deadlineUs) {
    pollAbort();
    delay(2);
  }
  logEvent("DWELL_END", label);
  return !aborted;
}

void setTrigger(bool high, const char *label) {
  digitalWrite(TRIG_OUT_PIN, high ? HIGH : LOW);
  delayMicroseconds(5);
  logEvent(high ? "BLOCK_TRIGGER_HIGH" : "BLOCK_TRIGGER_LOW", label);
}

bool beginBlock(const char *name) {
  if (aborted) {
    return false;
  }
  snprintf(activeBlock, sizeof(activeBlock), "%s", name);
  setTrigger(true, name);
  return true;
}

void endBlock() {
  setTrigger(false, activeBlock);
  snprintf(activeBlock, sizeof(activeBlock), "%s", "IDLE");
}

bool initialiseRmt() {
  rmt_tx_channel_config_t channelConfig = {};
  channelConfig.gpio_num = static_cast<gpio_num_t>(STEP_PIN);
  channelConfig.clk_src = RMT_CLK_SRC_DEFAULT;
  channelConfig.resolution_hz = RMT_RESOLUTION_HZ;
  channelConfig.mem_block_symbols = RMT_MAX_SYMBOLS;
  channelConfig.trans_queue_depth = 2;
  channelConfig.flags.invert_out = false;
  channelConfig.flags.with_dma = false;
  channelConfig.flags.io_loop_back = false;
  channelConfig.flags.io_od_mode = false;

  esp_err_t error = rmt_new_tx_channel(&channelConfig, &rmtChannel);
  if (error != ESP_OK) {
    Serial.printf("RMT channel creation failed: %s\n", esp_err_to_name(error));
    return false;
  }

  rmt_copy_encoder_config_t encoderConfig = {};
  error = rmt_new_copy_encoder(&encoderConfig, &rmtEncoder);
  if (error != ESP_OK) {
    Serial.printf("RMT encoder creation failed: %s\n", esp_err_to_name(error));
    return false;
  }

  error = rmt_enable(rmtChannel);
  if (error != ESP_OK) {
    Serial.printf("RMT enable failed: %s\n", esp_err_to_name(error));
    return false;
  }
  return true;
}

bool configureDriver() {
  DriverSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.en_spreadCycle(true);
  driver.intpol(false);
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  const uint8_t connection = driver.test_connection();
  Serial.printf("TMC2209 connection test: %u (0 means OK)\n", connection);
  if (connection != 0) {
    return false;
  }
  driverConfigured = true;
  return true;
}

bool applyRunConfiguration(const CurrentLevel &current, uint16_t mres) {
  if (positionU16 != 0) {
    Serial.println("Refusing MRES/current change away from the run origin.");
    aborted = true;
    return false;
  }
  if (!(mres == 16 || mres == 4 || mres == 2 || mres == 1)) {
    Serial.println("Unsupported MRES.");
    aborted = true;
    return false;
  }

  activeCurrent = &current;
  activeMres = mres;
  driver.rms_current(current.setRmsMa);
  const uint8_t configuredIrun = driver.irun();
  driver.ihold(configuredIrun);
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);
  driver.microsteps(mres);
  driver.en_spreadCycle(true);
  driver.intpol(false);

  digitalWrite(EN_PIN, LOW);
  delay(10);
  logEvent("RUN_CONFIG", "SpreadCycle_intpol_off_hold_equals_run");
  Serial.printf(
      "# current=%s set_rms_mA=%u measured_rms_mA=%u IRUN=%u IHOLD=%u "
      "MRES=%u readback_MRES=%u\n",
      current.name, current.setRmsMa, current.measuredRmsMa, configuredIrun,
      driver.ihold(), mres, driver.microsteps());
  return true;
}

bool updatePosition(int direction, uint32_t pulseCount) {
  const int64_t unitsPerPulse = 16 / activeMres;
  const int64_t candidate =
      positionU16 + static_cast<int64_t>(direction) *
                        static_cast<int64_t>(pulseCount) * unitsPerPulse;
  if (llabs(candidate) > MAX_ABS_POSITION_U16) {
    Serial.printf(
        "Travel guard rejected move: candidate=%lld u16, limit=%lld u16\n",
        static_cast<long long>(candidate),
        static_cast<long long>(MAX_ABS_POSITION_U16));
    aborted = true;
    return false;
  }
  positionU16 = candidate;
  return true;
}

bool pulseBatch(int32_t signedPulses, float pulseRateHz, const char *label) {
  if (aborted || signedPulses == 0) {
    return !aborted;
  }
  const uint32_t pulseCount =
      static_cast<uint32_t>(signedPulses > 0 ? signedPulses : -signedPulses);
  if (pulseCount > RMT_MAX_SYMBOLS) {
    Serial.println("Pulse batch exceeds the fixed RMT symbol buffer.");
    aborted = true;
    return false;
  }
  if (!(pulseRateHz > 0.0F)) {
    Serial.println("Pulse rate must be positive.");
    aborted = true;
    return false;
  }

  const int direction = signedPulses > 0 ? 1 : -1;
  const uint32_t periodTicks =
      static_cast<uint32_t>(llround(RMT_RESOLUTION_HZ / pulseRateHz));
  if (periodTicks <= RMT_HIGH_TICKS || periodTicks > 32767) {
    Serial.printf("RMT period out of range: %lu ticks\n",
                  static_cast<unsigned long>(periodTicks));
    aborted = true;
    return false;
  }

  if (!updatePosition(direction, pulseCount)) {
    return false;
  }

  digitalWrite(DIR_PIN, direction > 0 ? HIGH : LOW);
  delayMicroseconds(DIR_SETUP_US);
  for (uint32_t index = 0; index < pulseCount; ++index) {
    rmtSymbols[index].level0 = 1;
    rmtSymbols[index].duration0 = RMT_HIGH_TICKS;
    rmtSymbols[index].level1 = 0;
    rmtSymbols[index].duration1 = periodTicks - RMT_HIGH_TICKS;
  }

  rmt_transmit_config_t transmitConfig = {};
  transmitConfig.loop_count = 0;
  transmitConfig.flags.eot_level = 0;
  transmitConfig.flags.queue_nonblocking = false;

  const int64_t submitUs = esp_timer_get_time();
  esp_err_t error =
      rmt_transmit(rmtChannel, rmtEncoder, rmtSymbols,
                   pulseCount * sizeof(rmt_symbol_word_t), &transmitConfig);
  if (error == ESP_OK) {
    error = rmt_tx_wait_all_done(rmtChannel, -1);
  }
  const int64_t doneUs = esp_timer_get_time();
  if (error != ESP_OK) {
    Serial.printf("RMT transmission failed: %s\n", esp_err_to_name(error));
    aborted = true;
    return false;
  }

  // submit_us + k*period_ticks reconstructs every commanded rising edge.
  logEvent("PULSE_BATCH", label, direction, pulseCount, pulseRateHz,
           periodTicks, submitUs, doneUs);
  pollAbort();
  return !aborted;
}

bool burst(int32_t pulses, const char *label) {
  return pulseBatch(pulses, BURST_FULL_STEPS_S * activeMres, label);
}

bool physicalMove(int32_t fullSteps, float fullStepRate, const char *label) {
  return pulseBatch(fullSteps * static_cast<int32_t>(activeMres),
                    fullStepRate * activeMres, label);
}

bool runReference(const char *blockName) {
  if (!beginBlock(blockName)) return false;
  if (!dwellMs(2000, "lead_in")) return false;
  for (size_t index = 0;
       index < sizeof(REFERENCE_MOVES) / sizeof(REFERENCE_MOVES[0]); ++index) {
    char label[24];
    snprintf(label, sizeof(label), "reference_%u",
             static_cast<unsigned>(index + 1));
    if (!burst(REFERENCE_MOVES[index], label)) return false;
    if (!dwellMs(1000, label)) return false;
  }
  if (!dwellMs(2000, "tail")) return false;
  endBlock();
  return true;
}

bool runConditioning(const char *targetBlock) {
  char blockName[48];
  snprintf(blockName, sizeof(blockName), "COND_BEFORE_%s", targetBlock);
  if (!beginBlock(blockName)) return false;
  if (!physicalMove(4, CONDITIONING_FULL_STEPS_S, "conditioning_plus"))
    return false;
  if (!physicalMove(-4, CONDITIONING_FULL_STEPS_S, "conditioning_minus"))
    return false;
  if (!dwellMs(2000, "conditioning_settle")) return false;
  endBlock();
  return true;
}

bool runA1() {
  if (!beginBlock("A1")) return false;
  for (const int16_t n : N_VALUES) {
    char label[24];
    snprintf(label, sizeof(label), "N%d_positive", n);
    for (uint8_t repeat = 0; repeat < LADDER_REPEATS; ++repeat) {
      if (!burst(n, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
    }
    snprintf(label, sizeof(label), "N%d_negative", n);
    for (uint8_t repeat = 0; repeat < LADDER_REPEATS; ++repeat) {
      if (!burst(-n, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
    }
  }
  endBlock();
  return true;
}

bool runA2() {
  if (!beginBlock("A2")) return false;
  for (const int16_t n : N_VALUES) {
    char label[24];
    snprintf(label, sizeof(label), "N%d_alternating", n);
    for (uint8_t repeat = 0; repeat < LADDER_REPEATS; ++repeat) {
      if (!burst(n, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
      if (!burst(-n, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
    }
  }
  endBlock();
  return true;
}

template <size_t N>
bool runLoopPattern(const char *label, const int16_t (&pattern)[N]) {
  for (uint8_t repeat = 0; repeat < LOOP_REPEATS; ++repeat) {
    for (size_t index = 0; index < N; ++index) {
      if (!burst(pattern[index], label) || !dwellMs(LOOP_DWELL_MS, label))
        return false;
    }
  }
  return true;
}

bool runB() {
  if (!beginBlock("B")) return false;
  if (!runLoopPattern("descending", NEST_DESC)) return false;
  if (!runLoopPattern("asymmetric", NEST_ASYM)) return false;
  if (!runLoopPattern("minor", NEST_MINOR)) return false;
  endBlock();
  return true;
}

bool runE() {
  if (!beginBlock("E")) return false;
  constexpr int16_t doubletN[] = {1, 2, 4, 8, 16};
  for (const int16_t n : doubletN) {
    char label[24];
    snprintf(label, sizeof(label), "N%d_doublet", n);
    for (uint8_t repeat = 0; repeat < DOUBLET_REPEATS; ++repeat) {
      if (!burst(n, label) || !burst(-n, label) ||
          !dwellMs(DOUBLET_DWELL_MS, label))
        return false;
    }
  }
  endBlock();
  return true;
}

bool checkAtOrigin(const char *where) {
  if (positionU16 == 0) {
    logEvent("ORIGIN_CHECK_OK", where);
    return true;
  }
  Serial.printf("Origin check failed at %s: %lld u16\n", where,
                static_cast<long long>(positionU16));
  aborted = true;
  return false;
}

bool runOneConfiguration() {
  if (!runReference("BLOCK_0_START") || !checkAtOrigin("BLOCK_0_START"))
    return false;
  if (!runConditioning("A1") || !runA1() || !checkAtOrigin("A1"))
    return false;
  if (!runConditioning("A2") || !runA2() || !checkAtOrigin("A2"))
    return false;
  if (!runConditioning("B") || !runB() || !checkAtOrigin("B"))
    return false;
  if (!runConditioning("E") || !runE() || !checkAtOrigin("E"))
    return false;
  if (!runReference("BLOCK_0_END") || !checkAtOrigin("BLOCK_0_END"))
    return false;
  return true;
}

void safeStop(const char *reason) {
  digitalWrite(TRIG_OUT_PIN, LOW);
  digitalWrite(EN_PIN, HIGH);
  snprintf(activeBlock, sizeof(activeBlock), "%s", "IDLE");
  logEvent("SAFE_STOP", reason);
}

void printPlan() {
  Serial.println("# Validated ESP32 campaign plan");
  Serial.println("# backends: Block 0 + conditioning + A1 + A2 + B + E");
  Serial.println("# MRES: 16,4,2,1");
  Serial.println("# currents RMS set/readback mA: 360/355,600/556,750/715");
  Serial.println("# executions: 4 MRES x 3 currents = 12");
  Serial.println("# all blocks and every execution are net zero");
  Serial.println("# maximum planned excursion: 12800 units of 1/16 full step");
}

void runCampaign() {
  aborted = false;
  positionU16 = 0;
  runIndex = 0;
  printCsvHeader();
  logEvent("CAMPAIGN_START", "RUN");

  for (const uint16_t mres : MRES_VALUES) {
    for (const CurrentLevel &current : CURRENT_LEVELS) {
      ++runIndex;
      if (!applyRunConfiguration(current, mres) || !runOneConfiguration()) {
        safeStop(aborted ? "ABORTED" : "FAILED");
        return;
      }
      logEvent("RUN_COMPLETE", "origin_verified");
    }
  }
  safeStop("CAMPAIGN_COMPLETE");
}

}  // namespace rig

void setup() {
  using namespace rig;
  pinMode(EN_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  pinMode(TRIG_OUT_PIN, OUTPUT);
  pinMode(TRIG_ECHO_PIN, INPUT_PULLDOWN);
  digitalWrite(EN_PIN, HIGH);
  digitalWrite(DIR_PIN, LOW);
  digitalWrite(TRIG_OUT_PIN, LOW);

  Serial.begin(CONSOLE_BAUD);
  Serial.setTimeout(50);
  const uint32_t waitStart = millis();
  while (!Serial && millis() - waitStart < 3000) {
    delay(10);
  }
  Serial.println("\n# ESP32-S3/TMC2209 identification runner");

  if (!initialiseRmt() || !configureDriver()) {
    safeStop("INITIALISATION_FAILED");
    Serial.println("# Correct the hardware connection before issuing RUN.");
    return;
  }
  printPlan();
  Serial.println("# Enter CHECK or RUN. RUN enables physical motion.");
}

void loop() {
  using namespace rig;
  if (!Serial.available()) {
    delay(10);
    return;
  }
  String command = Serial.readStringUntil('\n');
  command.trim();
  command.toUpperCase();
  if (command == "CHECK") {
    printPlan();
    Serial.println("# CHECK complete; motor remains disabled.");
  } else if (command == "RUN") {
    if (!driverConfigured) {
      Serial.println("# RUN rejected: TMC2209 initialisation did not pass.");
    } else {
      runCampaign();
      Serial.println("# Enter CHECK or RUN.");
    }
  } else if (command == "ABORT") {
    aborted = true;
    safeStop("ABORT_WHILE_IDLE");
  } else if (command.length() > 0) {
    Serial.println("# Unknown command. Use CHECK, RUN, or ABORT.");
  }
}

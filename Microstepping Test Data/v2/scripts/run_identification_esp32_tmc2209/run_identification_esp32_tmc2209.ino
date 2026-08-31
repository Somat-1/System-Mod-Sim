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
 *   VISIBLE move +100 full steps, dwell, and return -100 full steps
 *   QUICK12 representative A1/A2/B/E diagnostic at 1/2 and full step
 *   RECOVER38 one-time recovery of the known interrupted-run offset
 *   RUN    execute all 12 MRES/current combinations
 *   ABORT  stop at the next safe batch/dwell boundary
 */

#include <Arduino.h>
#include <TMCStepper.h>
#include <driver/rmt_tx.h>
#include <esp_err.h>
#include <esp_timer.h>

namespace rig {

// ESP32-S3 <-> TMC2209 wiring.
// VIO -> 3V3; logic GND -> common logic/motor GND; CLK left floating.
// MS1=GND and MS2=GND select DRIVER_ADDRESS=0 below.
// EN requires an external 10 kohm pull-up to 3V3 so the driver remains
// disabled before the ESP32 firmware configures this output.
constexpr uint8_t STEP_PIN = 6;
constexpr uint8_t DIR_PIN = 7;
constexpr uint8_t EN_PIN = 5;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t DIAG_PIN = 4;

// Independent measurement-system synchronization pair. These are not
// TMC2209 pins; GPIO 1/2 became available after correcting STEP and DIR.
constexpr uint8_t TRIG_OUT_PIN = 1;
constexpr uint8_t TRIG_ECHO_PIN = 2;

constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.03F;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint32_t CONSOLE_BAUD = 115200;
constexpr uint16_t MOTOR_FULL_STEPS_PER_REV = 200;

constexpr uint32_t RMT_RESOLUTION_HZ = 1000000;
constexpr uint32_t RMT_HIGH_TICKS = 2;
constexpr size_t RMT_MAX_SYMBOLS = 64;
// ESP32-S3 has 48 RMT symbols per hardware memory block. Two blocks avoid
// an exact-fit boundary for the largest 64-pulse campaign batch.
constexpr size_t RMT_MEM_BLOCK_SYMBOLS = 96;
constexpr uint32_t DIR_SETUP_US = 5;

constexpr float BURST_FULL_STEPS_S = 250.0F;
constexpr float CONDITIONING_FULL_STEPS_S = 150.0F;
constexpr uint32_t DWELL_LADDER_MS = 400;
constexpr uint32_t LOOP_DWELL_MS = 300;
constexpr uint32_t DOUBLET_DWELL_MS = 1000;
constexpr uint8_t LADDER_REPEATS = 25;
constexpr uint8_t LOOP_REPEATS = 10;
constexpr uint8_t DOUBLET_REPEATS = 20;
constexpr float MARKER_FULL_STEPS_S = 50.0F;
constexpr uint32_t MARKER_DWELL_MS = 500;
constexpr uint32_t EXPECTED_RUNS = 12;
constexpr uint32_t MEASURED_PHASES_PER_RUN = 6;
constexpr uint32_t STEP_CONDITIONS_PER_RUN = 20;  // A1 6 + A2 6 + B 3 + E 5

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
uint32_t completedRuns = 0;
uint32_t completedMeasuredPhases = 0;
uint32_t completedStepConditions = 0;
uint8_t quickLogicalPulseScale = 1;
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
  channelConfig.mem_block_symbols = RMT_MEM_BLOCK_SYMBOLS;
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
  driver.toff(5);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.en_spreadCycle(true);
  driver.intpol(false);
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  uint8_t connection = 2;
  for (uint8_t attempt = 1; attempt <= 10; ++attempt) {
    connection = driver.test_connection();
    Serial.printf("TMC2209 connection test attempt %u: %u (0 means OK)\n",
                  attempt, connection);
    if (connection == 0) break;
    delay(200);
  }
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
  // Write TMC2209 CHOPCONF.MRES bits 27:24 by explicit raw mask. The mixed-type
  // bitfield in TMCStepper 0.7.3 loses bit 3 on this ESP32 toolchain: codes
  // 0..7 work, but full-step code 8 is emitted as code 0 (1/256).
  const uint8_t requestedMresCode =
      mres == 16 ? 4 : mres == 4 ? 6 : mres == 2 ? 7 : 8;
  driver.en_spreadCycle(true);
  driver.intpol(false);
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  uint32_t chopconf = driver.CHOPCONF();
  chopconf = (chopconf & ~MRES_MASK) |
             (static_cast<uint32_t>(requestedMresCode) << 24);
  driver.CHOPCONF(chopconf);

  const uint32_t chopconfReadback = driver.CHOPCONF();
  const uint8_t readbackMresCode =
      static_cast<uint8_t>((chopconfReadback & MRES_MASK) >> 24);
  const uint16_t readbackMres =
      readbackMresCode <= 8 ? (256U >> readbackMresCode) : 0;
  if (readbackMresCode != requestedMresCode || readbackMres != mres) {
    Serial.printf("MRES configuration failed: requested=%u code=%u "
                  "readback=%u code=%u\n", mres, requestedMresCode,
                  readbackMres, readbackMresCode);
    digitalWrite(EN_PIN, HIGH);
    aborted = true;
    return false;
  }

  digitalWrite(EN_PIN, LOW);
  delay(10);
  logEvent("RUN_CONFIG", "SpreadCycle_intpol_off_hold_equals_run");
  Serial.printf(
      "# current=%s set_rms_mA=%u measured_rms_mA=%u IRUN=%u IHOLD=%u "
      "MRES=%u MRES_code=%u readback_MRES=%u readback_code=%u\n",
      current.name, current.setRmsMa, current.measuredRmsMa, configuredIrun,
      driver.ihold(), mres, requestedMresCode, readbackMres,
      readbackMresCode);
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

  // RMT ticks reconstruct every relative rising-edge time exactly.
  // submitUs/doneUs bracket the transaction in the console clock domain.
  logEvent("PULSE_BATCH", label, direction, pulseCount, pulseRateHz,
           periodTicks, submitUs, doneUs);
  pollAbort();
  return !aborted;
}

bool burst(int32_t pulses, const char *label) {
  return pulseBatch(pulses, BURST_FULL_STEPS_S * activeMres, label);
}

bool checkAtOrigin(const char *where);
void safeStop(const char *reason);

bool physicalMove(int32_t fullSteps, float fullStepRate, const char *label) {
  int32_t remaining = fullSteps * static_cast<int32_t>(activeMres);
  while (remaining != 0) {
    const int32_t batch = remaining > 0
                              ? min(remaining, static_cast<int32_t>(RMT_MAX_SYMBOLS))
                              : max(remaining, -static_cast<int32_t>(RMT_MAX_SYMBOLS));
    if (!pulseBatch(batch, fullStepRate * activeMres, label)) return false;
    remaining -= batch;
  }
  return true;
}

// Markers are deliberately slow, net-zero physical moves. Configuration
// markers move positive first and encode run number by amplitude. Phase
// markers move negative first and encode the upcoming phase by amplitude.
bool runMarker(const char *name, int32_t amplitudeFullSteps,
               bool positiveFirst) {
  char blockName[48];
  snprintf(blockName, sizeof(blockName), "MARKER_%s", name);
  if (!beginBlock(blockName)) return false;
  logEvent("MARKER_START", name);
  const int32_t signedAmplitude = positiveFirst ? amplitudeFullSteps
                                                 : -amplitudeFullSteps;
  if (!dwellMs(MARKER_DWELL_MS, "marker_lead") ||
      !physicalMove(signedAmplitude, MARKER_FULL_STEPS_S, "marker_out") ||
      !dwellMs(MARKER_DWELL_MS, "marker_endpoint") ||
      !physicalMove(-signedAmplitude, MARKER_FULL_STEPS_S, "marker_return") ||
      !dwellMs(MARKER_DWELL_MS, "marker_returned"))
    return false;
  logEvent("MARKER_COMPLETE", name);
  endBlock();
  return checkAtOrigin(name);
}

bool runCampaignStartSignature() {
  if (!beginBlock("MARKER_CAMPAIGN_START")) return false;
  logEvent("MARKER_START", "CAMPAIGN_START_40_20_40");
  const int16_t signature[] = {40, -40, 20, -20, 40, -40};
  for (const int16_t move : signature) {
    if (!physicalMove(move, MARKER_FULL_STEPS_S, "campaign_signature") ||
        !dwellMs(MARKER_DWELL_MS, "campaign_signature"))
      return false;
  }
  logEvent("MARKER_COMPLETE", "CAMPAIGN_START_40_20_40");
  endBlock();
  return checkAtOrigin("CAMPAIGN_START_SIGNATURE");
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
    ++completedStepConditions;
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
    ++completedStepConditions;
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
  ++completedStepConditions;
  if (!runLoopPattern("asymmetric", NEST_ASYM)) return false;
  ++completedStepConditions;
  if (!runLoopPattern("minor", NEST_MINOR)) return false;
  ++completedStepConditions;
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
    ++completedStepConditions;
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
  if (!runMarker("BLOCK_0_START", 10, false)) return false;
  if (!runReference("BLOCK_0_START") || !checkAtOrigin("BLOCK_0_START"))
    return false;
  ++completedMeasuredPhases;
  if (!runMarker("A1", 20, false)) return false;
  if (!runConditioning("A1") || !runA1() || !checkAtOrigin("A1"))
    return false;
  ++completedMeasuredPhases;
  if (!runMarker("A2", 30, false)) return false;
  if (!runConditioning("A2") || !runA2() || !checkAtOrigin("A2"))
    return false;
  ++completedMeasuredPhases;
  if (!runMarker("B", 40, false)) return false;
  if (!runConditioning("B") || !runB() || !checkAtOrigin("B"))
    return false;
  ++completedMeasuredPhases;
  if (!runMarker("E", 50, false)) return false;
  if (!runConditioning("E") || !runE() || !checkAtOrigin("E"))
    return false;
  ++completedMeasuredPhases;
  if (!runMarker("BLOCK_0_END", 60, false)) return false;
  if (!runReference("BLOCK_0_END") || !checkAtOrigin("BLOCK_0_END"))
    return false;
  ++completedMeasuredPhases;
  return true;
}

// Short IDS-validation sequence. It deliberately uses the same burst rates and
// dwell times as the long campaign, but only representative subsets and fewer
// repeats. QUICK12 runs it at MRES=2 and MRES=1, where one STEP pulse represents
// 2.5 um and 5 um nominal travel respectively.
bool runQuickA1() {
  constexpr int16_t values[] = {1, 4, 16};
  if (!beginBlock("QUICK_A1")) return false;
  for (const int16_t n : values) {
    char label[24];
    snprintf(label, sizeof(label), "N%d_positive", n);
    for (uint8_t repeat = 0; repeat < 5; ++repeat)
      if (!burst(n * quickLogicalPulseScale, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
    snprintf(label, sizeof(label), "N%d_negative", n);
    for (uint8_t repeat = 0; repeat < 5; ++repeat)
      if (!burst(-n * quickLogicalPulseScale, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
  }
  endBlock();
  return checkAtOrigin("QUICK_A1");
}

bool runQuickA2() {
  constexpr int16_t values[] = {1, 4, 16};
  if (!beginBlock("QUICK_A2")) return false;
  for (const int16_t n : values) {
    char label[24];
    snprintf(label, sizeof(label), "N%d_alternating", n);
    for (uint8_t repeat = 0; repeat < 5; ++repeat) {
      if (!burst(n * quickLogicalPulseScale, label) || !dwellMs(DWELL_LADDER_MS, label) ||
          !burst(-n * quickLogicalPulseScale, label) || !dwellMs(DWELL_LADDER_MS, label)) return false;
    }
  }
  endBlock();
  return checkAtOrigin("QUICK_A2");
}

bool runQuickB() {
  if (!beginBlock("QUICK_B")) return false;
  for (size_t index = 0; index < sizeof(NEST_DESC) / sizeof(NEST_DESC[0]); ++index)
    if (!burst(NEST_DESC[index] * quickLogicalPulseScale, "descending") ||
        !dwellMs(LOOP_DWELL_MS, "descending")) return false;
  endBlock();
  return checkAtOrigin("QUICK_B");
}

bool runQuickE() {
  if (!beginBlock("QUICK_E")) return false;
  constexpr int16_t values[] = {1, 4, 16};
  for (const int16_t n : values) {
    char label[24];
    snprintf(label, sizeof(label), "N%d_doublet", n);
    for (uint8_t repeat = 0; repeat < 5; ++repeat)
      if (!burst(n * quickLogicalPulseScale, label) ||
          !burst(-n * quickLogicalPulseScale, label) ||
          !dwellMs(DOUBLET_DWELL_MS, label)) return false;
  }
  endBlock();
  return checkAtOrigin("QUICK_E");
}

void runQuick12Diagnostic() {
  aborted = false;
  positionU16 = 0;
  runIndex = 0;
  printCsvHeader();
  logEvent("QUICK12_START", "MRES2_then_MRES1");
  constexpr uint16_t logicalMres[] = {2, 1};
  // Native full-step code 8 does not persist on the attached module. Verify
  // half-step once, then create whole-step endpoints with two contiguous,
  // verified half-step pulses. This is explicit in the serial labels.
  if (!applyRunConfiguration(CURRENT_LEVELS[1], 2)) {
    safeStop("QUICK12_PREFLIGHT_FAILED");
    return;
  }
  digitalWrite(EN_PIN, HIGH);
  logEvent("QUICK12_PREFLIGHT_OK", "HALF_STEP_NATIVE_FULL_ENDPOINT_PAIRED");
  for (const uint16_t logicalResolution : logicalMres) {
    ++runIndex;
    quickLogicalPulseScale = logicalResolution == 1 ? 2 : 1;
    if (!applyRunConfiguration(CURRENT_LEVELS[1], 2)) {
      safeStop(aborted ? "QUICK12_ABORTED" : "QUICK12_FAILED");
      return;
    }
    logEvent("QUICK12_CONFIG_START", logicalResolution == 2
                 ? "HALF_STEP_NATIVE" : "FULL_STEP_ENDPOINT_PAIRED_HALF");
    if (!runMarker(logicalResolution == 2 ? "HALF_STEP_START" : "FULL_ENDPOINT_START",
                   logicalResolution == 2 ? 20 : 30, true) ||
        !runQuickA1() || !runMarker("QUICK_A2", 10, false) ||
        !runQuickA2() || !runMarker("QUICK_B", 15, false) ||
        !runQuickB() || !runMarker("QUICK_E", 20, false) ||
        !runQuickE() || !checkAtOrigin("QUICK_CONFIG_COMPLETE")) {
      safeStop(aborted ? "QUICK12_ABORTED" : "QUICK12_FAILED");
      return;
    }
    logEvent("QUICK12_CONFIG_COMPLETE", logicalResolution == 2
                 ? "HALF_STEP_NATIVE" : "FULL_STEP_ENDPOINT_PAIRED_HALF");
  }
  quickLogicalPulseScale = 1;
  logEvent("QUICK12_COMPLETE", "origin_verified");
  safeStop("QUICK12_COMPLETE");
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
  Serial.println("# per execution: Block 0 start, A1, A2, B, E, Block 0 end");
  Serial.println("# A1/A2 N: 1,2,4,8,16,32 pulses; B: all 3 loop patterns");
  Serial.println("# E N: 1,2,4,8,16 pulses (all protocol-defined E sizes)");
  Serial.println("# campaign marker: +40,-40,+20,-20,+40,-40 full steps");
  Serial.println("# configuration marker: positive-first, 20 + 5*run_index full steps");
  Serial.println("# phase markers: negative-first 10,20,30,40,50,60 full steps");
  Serial.println("# C and D are not ESP tests; run them with the dedicated-controller runner");
  Serial.println("# all blocks and every execution are net zero");
  Serial.println("# maximum planned excursion: 12800 units of 1/16 full step");
}

void runCampaign() {
  aborted = false;
  positionU16 = 0;
  runIndex = 0;
  completedRuns = 0;
  completedMeasuredPhases = 0;
  completedStepConditions = 0;
  printCsvHeader();
  logEvent("CAMPAIGN_START", "RUN");

  for (const uint16_t mres : MRES_VALUES) {
    for (const CurrentLevel &current : CURRENT_LEVELS) {
      ++runIndex;
      if (!applyRunConfiguration(current, mres) ||
          (runIndex == 1 && !runCampaignStartSignature())) {
        safeStop(aborted ? "ABORTED" : "FAILED");
        return;
      }
      char configurationMarker[24];
      snprintf(configurationMarker, sizeof(configurationMarker), "CONFIG_%02lu",
               static_cast<unsigned long>(runIndex));
      if (!runMarker(configurationMarker, 20 + 5 * runIndex, true) ||
          !runOneConfiguration()) {
        safeStop(aborted ? "ABORTED" : "FAILED");
        return;
      }
      ++completedRuns;
      logEvent("RUN_COMPLETE", "origin_verified");
    }
  }
  if (completedRuns != EXPECTED_RUNS ||
      completedMeasuredPhases != EXPECTED_RUNS * MEASURED_PHASES_PER_RUN ||
      completedStepConditions != EXPECTED_RUNS * STEP_CONDITIONS_PER_RUN) {
    logEvent("CAMPAIGN_COUNT_FAILED", "incomplete_execution_count");
    safeStop("FAILED");
    return;
  }
  logEvent("CAMPAIGN_COUNTS_OK", "12_runs_72_phases_240_step_conditions");
  safeStop("CAMPAIGN_COMPLETE");
}

void recoverInterruptedOffset38() {
  aborted = false;
  positionU16 = 0;
  runIndex = 0;
  printCsvHeader();
  if (!applyRunConfiguration(CURRENT_LEVELS[1], 16)) {
    safeStop("RECOVERY_CONFIG_FAILED");
    return;
  }
  // The preceding logged abort occurred at +38 units of 1/16 full step.
  positionU16 = 38;
  if (!beginBlock("RECOVER_ABORT_38") ||
      !pulseBatch(-38, MARKER_FULL_STEPS_S * activeMres, "return_to_origin") ||
      !dwellMs(1000, "recovery_settle")) {
    safeStop("RECOVERY_FAILED");
    return;
  }
  endBlock();
  if (!checkAtOrigin("RECOVER_ABORT_38")) {
    safeStop("RECOVERY_ORIGIN_FAILED");
    return;
  }
  safeStop("RECOVERY_COMPLETE");
}

bool visibleLeg(int direction, const char *label) {
  // At 1/16 MRES, 100 full steps are 1600 STEP pulses. The fixed RMT
  // buffer holds 64 pulses, so each leg is emitted as 25 contiguous,
  // individually logged batches. Every rising edge remains reconstructible
  // from submit_us, pulse_rate_hz and the batch pulse count.
  constexpr uint32_t FULL_STEPS = 100;
  constexpr uint32_t PULSES_PER_BATCH = 64;
  constexpr uint32_t BATCH_COUNT =
      FULL_STEPS * 16 / PULSES_PER_BATCH;
  constexpr float FULL_STEP_RATE_HZ = 50.0F;
  for (uint32_t batch = 0; batch < BATCH_COUNT; ++batch) {
    if (!pulseBatch(direction * static_cast<int32_t>(PULSES_PER_BATCH),
                    FULL_STEP_RATE_HZ * activeMres, label)) {
      return false;
    }
  }
  return true;
}

void runVisibleMotion() {
  aborted = false;
  if (positionU16 != 0) {
    Serial.printf("# VISIBLE rejected: commanded position is %lld u16, not zero.\n",
                  static_cast<long long>(positionU16));
    safeStop("VISIBLE_NOT_AT_ORIGIN");
    return;
  }

  runIndex = 1;
  printCsvHeader();
  logEvent("VISIBLE_START", "plus_100_then_return");
  if (!applyRunConfiguration(CURRENT_LEVELS[1], 16) ||
      !beginBlock("VISIBLE_100") ||
      !dwellMs(1000, "origin_dwell") ||
      !visibleLeg(+1, "forward_100_full_steps") ||
      !dwellMs(2000, "positive_endpoint") ||
      !visibleLeg(-1, "return_100_full_steps") ||
      !dwellMs(1000, "returned_dwell")) {
    safeStop(aborted ? "VISIBLE_ABORTED" : "VISIBLE_FAILED");
    return;
  }
  endBlock();
  if (!checkAtOrigin("VISIBLE_100")) {
    safeStop("VISIBLE_ORIGIN_FAILED");
    return;
  }
  logEvent("VISIBLE_COMPLETE", "origin_verified");
  safeStop("VISIBLE_COMPLETE");
}

}  // namespace rig

void setup() {
  using namespace rig;
  pinMode(EN_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  pinMode(DIAG_PIN, INPUT);
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
  Serial.println("# Enter CHECK, VISIBLE, QUICK12, RECOVER38, or RUN. Motion commands enable the motor.");
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
  } else if (command == "VISIBLE") {
    if (!driverConfigured) {
      Serial.println("# VISIBLE rejected: TMC2209 initialisation did not pass.");
    } else {
      runVisibleMotion();
      Serial.println("# Enter CHECK, VISIBLE, or RUN.");
    }
  } else if (command == "QUICK12") {
    if (!driverConfigured) {
      Serial.println("# QUICK12 rejected: TMC2209 initialisation did not pass.");
    } else {
      runQuick12Diagnostic();
      Serial.println("# Enter CHECK, VISIBLE, QUICK12, or RUN.");
    }
  } else if (command == "RUN") {
    if (!driverConfigured) {
      Serial.println("# RUN rejected: TMC2209 initialisation did not pass.");
    } else {
      runCampaign();
      Serial.println("# Enter CHECK or RUN.");
    }
  } else if (command == "RECOVER38") {
    if (!driverConfigured) {
      Serial.println("# RECOVER38 rejected: TMC2209 initialisation did not pass.");
    } else {
      recoverInterruptedOffset38();
      Serial.println("# Recovery command completed; do not issue RECOVER38 again.");
    }
  } else if (command == "ABORT") {
    aborted = true;
    safeStop("ABORT_WHILE_IDLE");
  } else if (command.length() > 0) {
    Serial.println("# Unknown command. Use CHECK, VISIBLE, QUICK12, RECOVER38, RUN, or ABORT.");
  }
}

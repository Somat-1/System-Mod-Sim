/*
 * TEMPORARY diagnostic sketch -- not part of the v4 campaign.
 *
 * Communication-only probe: opens the TMC2209 UART link and repeatedly polls
 * test_connection()/version()/GCONF/IOIN/DRV_STATUS. It never writes a
 * configuration register and never drives STEP or DIR, so it is safe to run
 * with the motor supply off -- the expected result in that state is a clean
 * "no response" report rather than a garbled one, since the TMC2209's logic
 * (and therefore its UART transceiver) is powered from the same VM rail as
 * the motor, not from the ESP's 3V3.
 *
 * STEP GPIO5, DIR GPIO6 are configured as outputs held LOW throughout but are
 * otherwise unused here. UART ESP_RX GPIO18 / ESP_TX GPIO17, matching the v4
 * campaign firmware exactly.
 *
 * Prints one CSV line per poll. Commands over USB CDC: STATUS (dump current
 * state on demand), ABORT (stop the polling loop; POLL resumes it).
 */

#include <Arduino.h>
#include <TMCStepper.h>

constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.03F;
constexpr uint32_t DRIVER_BAUD = 115200;

constexpr uint32_t POLL_INTERVAL_MS = 1000;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

bool polling = true;
uint32_t pollIndex = 0;

void colour(uint8_t red, uint8_t green, uint8_t blue) {
  rgbLedWrite(LED_PIN, red, green, blue);
}

void printHeader() {
  Serial.println(
      "timestamp_us,poll_index,uart_ok,version_hex,crc_error,gconf_hex,"
      "ioin_hex,drv_status_hex,detail");
}

// test_connection() returns 0 only when the read succeeds and the register
// contents are self-consistent; version() is read separately because a UART
// short or missing VM can still clock out a fixed pattern that passes a
// weaker check.
void pollOnce() {
  ++pollIndex;
  const uint8_t connectionResult = driver.test_connection();
  const uint8_t version = driver.version();
  const bool crcError = driver.CRCerror;
  const uint32_t gconf = driver.GCONF();
  const bool gconfCrc = driver.CRCerror;
  const uint32_t ioin = driver.IOIN();
  const uint32_t drvStatus = driver.DRV_STATUS();
  const bool uartOk = connectionResult == 0 && version == 0x21 && !gconfCrc;

  const char *detail = uartOk ? "responding"
                       : (connectionResult != 0)
                           ? "test_connection_failed;check_VM_power_and_wiring"
                           : "responded_but_unexpected_readback";

  Serial.printf("%llu,%lu,%u,0x%02X,%u,0x%08lX,0x%08lX,0x%08lX,%s\n",
                static_cast<unsigned long long>(esp_timer_get_time()),
                static_cast<unsigned long>(pollIndex),
                static_cast<unsigned>(uartOk), version,
                static_cast<unsigned>(crcError),
                static_cast<unsigned long>(gconf),
                static_cast<unsigned long>(ioin),
                static_cast<unsigned long>(drvStatus), detail);

  colour(uartOk ? 0 : 24, uartOk ? 22 : 0, 0);
}

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);

  colour(24, 24, 24);  // One white boot blink.
  delay(1000);
  colour(0, 0, 0);

  Serial.begin(115200);
  Serial.setTimeout(20);
  const uint32_t serialStart = millis();
  while (!Serial && millis() - serialStart < 3000) delay(10);
  // The S3's native USB-CDC only flushes once a host attaches, so anything
  // printed during re-enumeration is lost. Give the capture time to open.
  delay(4000);

  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();

  Serial.println("# TMP_V4_UART_PROBE (temporary diagnostic; no motion, no config writes)");
  Serial.println("# STEP=GPIO5(idle_low) DIR=GPIO6(idle_low) UART_RX=GPIO18 UART_TX=GPIO17");
  Serial.println("# EXPECT_NO_RESPONSE_IF_TMC_VM_UNPOWERED; COMMANDS=STATUS,ABORT,POLL");
  printHeader();
}

void loop() {
  if (Serial.available()) {
    String command = Serial.readStringUntil('\n');
    command.trim();
    command.toUpperCase();
    if (command == "ABORT") {
      polling = false;
      Serial.println("# POLLING_STOPPED");
    } else if (command == "POLL") {
      polling = true;
      Serial.println("# POLLING_RESUMED");
    } else if (command == "STATUS") {
      Serial.printf("# STATUS,polling=%u,poll_index=%lu\n", polling,
                    static_cast<unsigned long>(pollIndex));
    } else if (command.length() > 0) {
      Serial.printf("# UNKNOWN_COMMAND=%s\n", command.c_str());
    }
  }

  if (polling) {
    pollOnce();
    delay(POLL_INTERVAL_MS);
  } else {
    delay(20);
  }
}

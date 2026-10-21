// No-motion ESP32-S3/TMC2209 UART communication check.
// Pinout: STEP=7, DIR=8, EN=9, UART TX=17 through 1k, UART RX=18.
// EN is held HIGH and STEP is held LOW for the entire test.

#include <Arduino.h>
#include <TMCStepper.h>

constexpr uint8_t STEP_PIN = 7;
constexpr uint8_t DIR_PIN = 8;
constexpr uint8_t EN_PIN = 9;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr float R_SENSE_OHM = 0.11F;  // Not used to set current in this check.

HardwareSerial DriverSerial(1);
TMC2209Stepper driver(&DriverSerial, R_SENSE_OHM, DRIVER_ADDRESS);

static void printHex32(const char *name, uint32_t value) {
  Serial.printf("# %-10s = 0x%08lX\n", name, static_cast<unsigned long>(value));
}

void setup() {
  // Put all motion-related pins in a safe state before starting either UART.
  pinMode(EN_PIN, OUTPUT);
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(EN_PIN, HIGH);   // TMC ENN is active LOW: outputs disabled.
  digitalWrite(STEP_PIN, LOW);  // Never changed by this sketch.
  digitalWrite(DIR_PIN, LOW);

  Serial.begin(115200);
  const uint32_t waitStart = millis();
  while (!Serial && millis() - waitStart < 3000) delay(10);

  Serial.println();
  Serial.println("# TMC2209 no-motion UART diagnostic");
  Serial.println("# EN=GPIO9 HIGH (disabled), STEP=GPIO7 LOW, DIR=GPIO8 LOW");
  Serial.println("# UART RX=GPIO18, TX=GPIO17 through 1k, address=0");

  DriverSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  delay(100);

  // Read first, before changing any TMC configuration.
  const uint32_t ioinBefore = driver.IOIN();
  const uint32_t gstat = driver.GSTAT();
  const uint32_t chopconf = driver.CHOPCONF();
  const uint32_t drvStatus = driver.DRV_STATUS();
  const uint8_t ifcntBefore = driver.IFCNT();

  printHex32("IOIN", ioinBefore);
  printHex32("GSTAT", gstat);
  printHex32("CHOPCONF", chopconf);
  printHex32("DRV_STATUS", drvStatus);
  Serial.printf("# IFCNT      = %u\n", ifcntBefore);
  Serial.printf("# VERSION    = 0x%02X (TMC2209 expected: 0x21)\n",
                static_cast<unsigned>((ioinBefore >> 24) & 0xFF));
  Serial.printf("# IOIN pins  ENN=%u MS1=%u MS2=%u STEP=%u DIR=%u\n",
                driver.enn(), driver.ms1(), driver.ms2(),
                driver.step(), driver.dir());

  // driver.begin() only writes GCONF fields needed for UART/register control.
  // With ENN held HIGH and no STEP edges this cannot command motor motion.
  driver.begin();
  delay(50);
  const uint8_t ifcntAfter = driver.IFCNT();
  const uint32_t gconfAfter = driver.GCONF();
  const uint32_t ioinAfter = driver.IOIN();
  const uint8_t writesAccepted = static_cast<uint8_t>(ifcntAfter - ifcntBefore);

  printHex32("GCONF", gconfAfter);
  Serial.printf("# IFCNT after = %u (accepted writes: %u)\n",
                ifcntAfter, writesAccepted);

  const bool readsGood = ((ioinBefore >> 24) & 0xFF) == 0x21 &&
                         ioinBefore != 0xFFFFFFFFUL &&
                         ioinAfter != 0xFFFFFFFFUL;
  const bool writeGood = writesAccepted >= 2;
  const bool safePins = digitalRead(EN_PIN) == HIGH &&
                        digitalRead(STEP_PIN) == LOW;

  if (readsGood && writeGood && safePins) {
    Serial.println("# RESULT: PASS - TMC2209 UART reads and writes are working; motor stayed disabled.");
  } else {
    Serial.printf("# RESULT: FAIL - reads=%s writes=%s safe_pins=%s\n",
                  readsGood ? "OK" : "BAD",
                  writeGood ? "OK" : "BAD",
                  safePins ? "OK" : "BAD");
    Serial.println("# Check VM power, common GND, PDN/UART bridge, 1k TX resistor, and MS1/MS2=GND.");
  }
}

void loop() {
  // Continuously enforce the no-motion state.
  digitalWrite(EN_PIN, HIGH);
  digitalWrite(STEP_PIN, LOW);
  delay(100);
}

# Hardware identification runners

## Campaign split

The timing-critical and controller-paced measurements are deliberately split:

| Runner | Assigned measurements |
|---|---|
| ESP32-S3 + TMC2209 | Block 0, conditioning, A1, A2, B, E, final Block 0 |
| EVO dedicated controller, axis X | Block 0, conditioning, C, D, final Block 0 |

Each runner automatically covers four MRES settings by three current levels:
12 configurations per runner. One ESP32 `RUN` command and one Python
invocation therefore cover the full divided campaign. All configurations
begin at the same working origin; neither runner homes between configurations.

## Current mapping

| Level | TMC command RMS | TMC measured RMS | Dedicated `SC` peak |
|---|---:|---:|---:|
| I_lo | 360 mA | approximately 355 mA | 502 mA |
| I_mid | 600 mA | 556 mA | 786 mA |
| I_hi | 750 mA | 715 mA | 1011 mA |

Run and hold current are equal. The TMC runner writes `IHOLD = IRUN`,
`IHOLDDELAY = 0`, and `TPOWERDOWN = 0`.

## ESP32-S3 and TMC2209

Sketch:
`scripts/run_identification_esp32_tmc2209/run_identification_esp32_tmc2209.ino`

Requirements:

- Arduino-ESP32 3.x, which supplies the ESP-IDF 5 RMT TX driver.
- TMCStepper.
- ESP32-S3 target.

Wiring:

| Signal | GPIO |
|---|---:|
| STEP | 1 |
| DIR | 2 |
| EN | 5 |
| UART TX through 1 kohm | 4 |
| UART RX direct to bridged PDN/UART | 6 |
| TRIG_OUT | 7 |
| TRIG_ECHO, optional loopback from GPIO 7 | 15 |

TMC configuration is fixed at address 0, `R_SENSE = 0.03 ohm`,
SpreadCycle, UART MRES selection, and MicroPlyer disabled. MS1 and MS2 remain
grounded. DIAG is unused.

The RMT clock is 1 MHz. Each logged `PULSE_BATCH` row records direction,
pulse count, pulse rate, exact RMT period ticks, enqueue time, and completion
time. Relative STEP-edge timing is reconstructed exactly from the RMT ticks.
`TRIG_OUT` remains high for the duration of each block.

### Compile and upload

Use the already validated Arduino toolchain. With Arduino CLI and the usual
Espressif core package, the equivalent workflow is:

```powershell
arduino-cli compile --fqbn esp32:esp32:esp32s3 .\scripts\run_identification_esp32_tmc2209
arduino-cli upload --fqbn esp32:esp32:esp32s3 --port COM_PORT .\scripts\run_identification_esp32_tmc2209
```

Replace `COM_PORT` with the board port. Open the console at 115200 baud.

1. Send `CHECK`. This validates and prints the plan while EN remains disabled.
2. Confirm at least the planned four-revolution positive travel envelope from
   the fixed working origin.
3. Start interferometer acquisition and serial-log capture.
4. Send `RUN`.
5. Send `ABORT` if required. It is honored at the next RMT-batch or dwell
   boundary; trigger is driven low and EN is driven high.

The sketch produces CSV-formatted console output. Capture the complete serial
stream because it is the command/input record for later model comparison.

## EVO dedicated controller

Runner: `scripts/run_identification_dedicated_controller.py`

The script uses 115200 baud, 8N1, CR termination, axis X, and the documented
`DM/SM/SS/SC/MH/MA/SP/DP/DS/ME/MO/SO/CO` commands. It verifies `DM 0`,
uses absolute six-decimal `MA` targets, and polls `DS` until status 1.

Install pyserial only for live execution:

```powershell
python -m pip install pyserial
```

### Dry run

```powershell
python .\scripts\run_identification_dedicated_controller.py --dry-run
```

The dry run does not open a serial port or wait in real time. It generates all
12 configurations and writes the complete command/event CSV under
`data/hardware_runs`.

### Required preflight

Before live execution:

1. Confirm axis X and `DM 0`.
2. Confirm position units are revolutions using the separately supervised
   `SP X 0.000000`, `MR X 1.000000` interferometer test.
3. Confirm the `SM` direction flag with a small supervised move.
4. Determine the home direction/maximum steps, sensor mask, overtravel, and
   safe homing speed parameters.
5. Select the post-home working position.
6. Measure available positive and negative travel from that working origin.
   The positive allowance must cover the longest D trajectory, approximately
   6.1 motor revolutions including acceleration/deceleration distance.
7. Connect Status output 1 to the acquisition trigger input. The runner uses
   `SO 32` at block start and `CO 32` at block end.

### Live invocation

Supply the verified rig values explicitly:

```powershell
python .\scripts\run_identification_dedicated_controller.py --execute `
  --port COM_PORT `
  --confirm-position-units REVOLUTIONS `
  --direction DIRECTION_FLAG `
  --working-position-rev WORKING_REV `
  --positive-limit-rev POSITIVE_MARGIN_REV `
  --negative-limit-rev NEGATIVE_MARGIN_REV `
  --home-max-steps HOME_MAX_STEPS `
  --home-sensor-mask SENSOR_MASK `
  --home-overtravel OVERTRAVEL `
  --home-min-speed HOME_MIN_SPEED `
  --home-max-speed HOME_MAX_SPEED `
  --home-accel HOME_ACCEL `
  --home-decel HOME_DECEL `
  --home-ramp-type HOME_RAMP_TYPE
```

The homing speed values use the controller's documented 0.01-rad/s and
0.01-rad/s-squared units. If the controller has already been homed during the
same uninterrupted session, replace all home arguments with `--skip-home`.
The runner still moves to the working position and applies `SP X 0.000000`.

The three slow D rates, 0.125, 0.375, and 1.25 full steps/s, are software
paced as individual absolute microstep moves. Every acknowledgment is
timestamped. The remaining five rates use `ACCEL = DECEL = 628` and
`RAMPTYPE = 1`; the log marks the first 0.5 s as discarded identification
data.

Ctrl+C requests cancellation. Whether the run completes or fails, the script
attempts `CO 32` followed by `MO X` and records the shutdown result.

## Session order

1. Home once at the beginning of a dedicated-controller session.
2. Move to the common working position and define it with `SP X 0.000000`.
3. Run the dedicated-controller campaign without re-homing.
4. Establish the same physical working origin for the ESP32 rig and verify its
   travel envelope.
5. Run `CHECK`, then `RUN`, while capturing the ESP32 serial log.
6. Do not home between MRES/current configurations. Re-home only after a
   deliberate relocation, then repeat Block 0 before accepting new data.

# Hardware identification runners

## Campaign split

The timing-critical and controller-paced measurements are deliberately split:

| Runner | Assigned measurements |
|---|---|
| ESP32-S3 + TMC2209 | Block 0, conditioning, A1, A2, B, E, final Block 0 |
| EVO dedicated controller, axis X | Block 0, conditioning, C, D, final Block 0 |

The active dedicated-controller campaign covers MRES 1/4, 1/2, and 1/1 at
two peak-current levels: six configurations in one Python invocation. The
older ESP32 runner retains its original matrix. All dedicated-controller
configurations begin at the same working origin and do not re-home.

## Active dedicated-controller current mapping

| Level | Relative setting | Dedicated `SC` peak |
|---|---:|---:|
| I_50pct | 50% | 200 mA |
| I_100pct | 100% | 400 mA |

Run and hold current are equal. The controller manual defines `SC` in mA
Ipeak; 400 mA is therefore entered directly as `SC X 400 400`.

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
controller commands. It verifies `DM 0` and polls `DS` until status 1. The
standard homed workflow uses absolute `MA`; current-position mode captures
`DP X` and uses signed relative `MR` commands without `SP`.

Install pyserial only for live execution:

```powershell
python -m pip install pyserial
```

### Dry run

```powershell
python .\scripts\run_identification_dedicated_controller.py --dry-run
```

The dry run does not open a serial port or wait in real time. It generates all
six configurations and writes the complete command/event CSV under
`data/hardware_runs`.

### Short diagnostic first

Runner: `scripts/run_dedicated_controller_diagnostic.py`

The diagnostic uses MRES 1/4 at `SC X 400 400` and lasts approximately
92.3 seconds, excluding initialization and serial overhead. It exercises the reference
fingerprint, conditioning, both shortened creep approaches, the fastest
software-paced plateau (1.25 full steps/s), the fastest controller-paced
plateau (200 full steps/s), and the data-visible marker signatures.

Dry-run command:

```powershell
python .\scripts\run_dedicated_controller_diagnostic.py --dry-run
```

The diagnostic has the verified local rig settings embedded: COM5, axis X
(J19), direction 0, MRES 1/4, and 400 mA peak. It performs no homing and sends
no positioning move during initialization. It captures the displayed starting
position with `DP X`, then sends every move as a signed relative `MR` command.

Run it with:

```powershell
python .\scripts\run_dedicated_controller_diagnostic.py --execute
```

The script initializes the controller and then waits automatically. Start IDS
acquisition at that prompt and press Enter to begin the first marker. Retain
both the IDS CSV and controller CSV. After a successful diagnostic, leave the
controller powered and do not move the stage. The subsequent full-run command
can use `--use-current-position-as-origin` to retain the same relative-motion
workflow without homing or `SP`.

### Data-visible separator signatures

Every marker makes a rapid negative leap, dwells for 1.0 second, returns by
the same positive amount, and settles for 0.5 second. Configuration markers
are 68, 72, 76, 80, 84, and 88 full steps for runs 1 through 6. Test markers
identify the following block by amplitude:

| Following test | Marker amplitude (full steps) |
|---|---:|
| Conditioning before C | 12 |
| Conditioning before D | 16 |
| C | 20 |
| D 0.125 | 24 |
| D 0.375 | 28 |
| D 1.25 | 32 |
| D 3.5 | 36 |
| D 9.5 | 40 |
| D 27.5 | 44 |
| D 70 | 48 |
| D 200 | 52 |
| Final Block 0 reference | 56 |

Each marker also has its own `SO 32`/`CO 32` trigger interval and a
`MARKER_SIGNATURE` row in the controller CSV. The planned six-configuration
motion time is 2552.1 seconds (42 min 32 s); allow approximately 46 minutes
in acquisition for serial polling and command overhead.

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
  --home-ramp-type HOME_RAMP_TYPE `
  --wait-for-acquisition
```

The homing speed values use the controller's documented 0.01-rad/s and
0.01-rad/s-squared units. If the controller has already been homed during the
same uninterrupted session, replace all home arguments with `--skip-home`.
The runner still moves to the working position and applies `SP X 0.000000`.
After the diagnostic has already established that origin, prefer
`--reuse-working-origin`; it performs neither operation and first verifies the
displayed position is zero.

The three slow D rates, 0.125, 0.375, and 1.25 full steps/s, are software
paced as individual absolute microstep moves. Every acknowledgment is
timestamped. The remaining five rates use `ACCEL = DECEL = 628` and
`RAMPTYPE = 1`; the log marks the first 0.5 s as discarded identification
data.

Ctrl+C requests cancellation. Whether the run completes or fails, the script
attempts `CO 32` followed by `MO X` and records the shutdown result.

## Session order

1. Place the stage at a starting position with sufficient travel both ways.
2. Run the diagnostic; it captures `DP X` and uses only relative `MR` moves.
3. Record and review the approximately 92-second diagnostic.
4. Leave the controller powered and the stage untouched.
5. For the later full run, use `--use-current-position-as-origin` and start IDS
   acquisition at the prompt.
6. Do not relocate the stage between MRES/current configurations. After a
   deliberate relocation, repeat the diagnostic before accepting data.

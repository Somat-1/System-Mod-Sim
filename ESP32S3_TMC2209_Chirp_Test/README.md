# ESP32-S3 + TMC2209 bounded chirp test

This folder contains a standalone Arduino sketch for a Dewesoft accelerometer
resonance survey. It performs:

1. 2 s stationary lead-in representing 0 Hz.
2. Linear up-sweep from 1 Hz to 1 kHz in 30 s.
3. 0.5 s stationary turnaround.
4. Linear down-sweep from 1 kHz to 1 Hz in 30 s.
5. 2 s stationary tail.

The digital trigger is high for the complete 64.5 s sequence. The first
nonzero frequency is 1 Hz because a literal 0 Hz oscillation has an infinite
period.

## Why the motion is bounded

This is not a one-direction STEP-frequency ramp. The firmware alternates one
positive microstep with one negative microstep. Commanded position therefore
toggles only between the starting point and +1 configured microstep and
returns to the origin after every cycle and at every exit.

Defaults are:

- MRES: 1/4 step.
- Current: 360 mA RMS.
- 2 mm screw pitch and 200 full steps/revolution.
- Theoretical position increment: 2.5 micrometres.
- Maximum STEP transition rate: 2 kHz at a 1 kHz excitation frequency.
- SpreadCycle enabled; MicroPlyer disabled.

The TMC2209 and ESP32-S3 have ample timing margin at this rate. The important
limitation is mechanical, not digital timing: an individual microstep may not
produce its theoretical displacement because of detent torque, friction and
load compliance. Consequently this test is suitable for a commanded-position
to acceleration resonance survey. It is not by itself a calibrated
force-to-acceleration FRF. For a quantitative physical FRF, also measure the
actual stage motion or applied force/current.

## Resonance exclusion band

The sketch suppresses STEP commands from 120 to 230 Hz on both sweeps. This
band brackets the motor-dominated first resonance predicted by the current
model. Dewesoft data inside the band must be marked as deliberately
unexcited, not interpreted as a low response.

After a low-amplitude physical pilot, edit these constants if the observed
motor resonance lies elsewhere:

```cpp
constexpr bool RESONANCE_NOTCH_ENABLED = true;
constexpr float NOTCH_LOW_HZ = 120.0F;
constexpr float NOTCH_HIGH_HZ = 230.0F;
```

Do not disable the notch until the pilot confirms that resonance motion is
controlled.

## Wiring

The GPIO definitions are copied from the currently validated v2 firmware,
not the older table in `HARDWARE_RUNNERS.md`:

| Signal | ESP32-S3 GPIO |
|---|---:|
| STEP | 6 |
| DIR | 7 |
| EN | 5 |
| TMC UART TX | 17 |
| TMC UART RX | 18 |
| Dewesoft trigger output | 1 |
| Optional trigger echo input | 2 |

TMC2209 address is 0 and `R_SENSE` is 0.03 ohm. Maintain a common ground
between the ESP32, driver, motor supply and Dewesoft digital input. Confirm
the Dewesoft input voltage compatibility before connecting GPIO 1.

## Compile and upload

Requirements:

- Arduino-ESP32 3.x.
- TMCStepper.
- ESP32-S3 board target.

From the `Sytem Mod & Sim` directory:

```powershell
arduino-cli compile --fqbn esp32:esp32:esp32s3 .\ESP32S3_TMC2209_Chirp_Test
arduino-cli upload --fqbn esp32:esp32:esp32s3 --port COM_PORT .\ESP32S3_TMC2209_Chirp_Test
```

Replace `COM_PORT` with the ESP32 port and open its serial console at 115200
baud.

## Test procedure

1. Centre the stage and remove anything that could collide.
2. Start the serial console and send `CHECK`. The motor remains disabled.
3. Confirm the reported TMC connection and MRES readback are OK.
4. Start Dewesoft acquisition, including GPIO 1 as the digital trigger.
5. Send `RUN`.
6. Send `ABORT` if unexpected motion or noise occurs. The firmware returns
   the single outstanding microstep where possible, lowers the trigger and
   disables the motor.
7. Keep the serial output with the Dewesoft file. Its trigger edge and sweep
   segment timestamps define the commanded frequency law for processing.

Perform the first run with the mechanism unloaded and observed directly. A
bounded command cannot prevent every mechanical resonance amplification, so
the physical emergency stop and travel limits remain authoritative.

# ESP32-S3 + TMC2209 step-size sweep

## Purpose

For each microstep resolution `MRES` in `{1, 2, 4, 8, 16, 32, 64}`: step
the motor one microstep forward, wait, one microstep back, wait,
repeating for ~35 s, bracketed by a large out-and-back marker move before
and after each block. Checks whether a single commanded microstep at each
resolution produces the physically-correct step size, independent of any
adaptive current-control feature that could distort it.

## Hardware/firmware conventions (reused, not reinvented)

Pin map, TMCStepper instantiation, UART setup, and the driver
configuration sequence are copied verbatim from the two existing sketches
in this repo:

- `../../../ESP32S3_TMC2209_Chirp_Test/ESP32S3_TMC2209_Chirp_Test.ino`
- `../../v2/scripts/run_identification_esp32_tmc2209/run_identification_esp32_tmc2209.ino`

Both already run the driver in **SpreadCycle-only mode**
(`driver.en_spreadCycle(true)`, `driver.intpol(false)`, constant
`IHOLD=IRUN`) and never call `pwm_autoscale()`/`TCOOLTHRS()`/`PWMCONF()`
(StealthChop) or any CoolStep register — this sketch reproduces that
exactly, so **StealthChop and CoolStep are never enabled** and cannot
perturb microstep positioning during the sweep, per your request.

`MRES` is set via the same raw `CHOPCONF` bits-27:24 write + readback
verification both existing sketches use (not `driver.microsteps()`,
which both existing sketches' comments note is unreliable in TMCStepper
0.7.3 on this ESP32 toolchain), generalized here to cover all seven
requested values (`mresCode = 8 - log2(mres)`).

Step pulses use the chirp test's simple bit-banged
`digitalWrite`/`delayMicroseconds` approach rather than the RMT-batch
approach from the other sketch — appropriate here since steps are
seconds apart, not a fast burst.

## What's new here

**StallGuard/DRV_STATUS logging** — neither existing sketch reads these
registers at all (confirmed by search: zero hits for `SG_RESULT`,
`DRV_STATUS`, `stallguard` anywhere in this repo). This sketch logs, once
per single-microstep step and once per marker move: `SG_RESULT` (raw
StallGuard load value), the raw `DRV_STATUS` register (hex), and its
individually decoded flags (`stst`, `otpw`, `ot`, `s2ga`, `s2gb`, `ola`,
`olb`, `stallguard`, `cs_actual`) via the TMCStepper library's per-bit
accessors.

**Marker sizing** — fixed at `MARKER_FULL_STEPS = 20` full-step-equivalents
(not 20 microsteps), so a marker's physical travel is identical across
every block regardless of MRES — always a large, easily recognizable
jump next to the (deliberately shrinking, as MRES rises) single-microstep
zig-zags under test.

**Position accounting** — a signed accumulator in units of 1/64 full step
(the finest unit needed to represent any MRES in {1,...,64} as an
integer), with a hard travel guard and an origin check after every block,
following the safety pattern already established by both existing
sketches (position-bound guard + `checkAtOrigin`-equivalent + unconditional
`EN_PIN` disable on any exit path).

## Files

- `esp32_tmc2209_stepsize_sweep/esp32_tmc2209_stepsize_sweep.ino` — the
  sketch (folder name matches the `.ino` filename per Arduino IDE
  convention). Serial commands: `CHECK`, `RUN`, `ABORT` (same convention
  as both existing sketches). CSV log over Serial at 115200 baud, header
  printed once at startup.

## Status

**Written and manually reviewed, not yet compiled or flashed** — no
Arduino toolchain was available in this environment to build it. Format
strings and register-accessor calls were checked by hand against the
TMCStepper API used in the two existing sketches. Before running on real
hardware: compile in the Arduino IDE/arduino-cli with the ESP32 board
package and TMCStepper library installed, and verify on the bench with a
short single-MRES test before committing to the full 7-block sweep
(~7 × (2 × marker + 35 s) ≈ 4-5 minutes total).

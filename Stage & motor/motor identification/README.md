# Motor identification

Rotor-side identification of the TMC2209-driven 200-step hybrid motor using
the encoder on the motor shaft. This campaign separates static microstep error,
detent torque, friction hysteresis, step dynamics, and speed-dependent lag.

This is **not** a MicroPlyer or StealthChop comparison. Both are excluded:

- TMC2209 `en_spreadCycle=1` (SpreadCycle only);
- `TPWMTHRS=0` (no speed-dependent return to StealthChop);
- `intpol=0` (MicroPlyer disabled; only native MRES table entries);
- CoolStep disabled and `IHOLD=IRUN`.

The firmware reads `GCONF`, `CHOPCONF`, `DRV_STATUS`, and `IHOLD_IRUN` back. It
refuses to run if SpreadCycle, interpolation-off, inactive StealthChop, and
equal hold/run current are not all confirmed. The host checks the same startup
line independently.

## Contents

```
motor identification/
|-- README.md
|-- scripts/
|   |-- sequence_design.py      experiment definition and validation
|   |-- render_sequence.py      PNG/SVG campaign overview
|   `-- run_campaign.py         upload, mode check, recording, abort
|-- sequence/
|   |-- motor_identification_plan.csv
|   |-- motor_identification_summary.json
|   `-- motor_identification_sequence.png / .svg
|-- firmware/motor_identification/motor_identification.ino
`-- recordings/                 created by the host runner
```

## Default campaign

The generated plan is 37.1 minutes before operator overhead and contains 83
actions. It starts and ends at the same commanded electrical phase.

| Block | Time | Measurement |
|---|---:|---|
| Warm-up/reference | 4.10 min | opening standstill and six bidirectional run-in pairs |
| Electrical alignment | 0.26 min | one electrical cycle out and back |
| Static map | 16.00 min | all 3,200 native MRES=16 positions over one shaft revolution, then the identical positions in reverse |
| Current sweep | 3.27 min | four electrical cycles out/back at 50, 75, and 100% of the nominal test current |
| Dynamics | 2.43 min | 12 full-step and one-microstep out/back responses in both directions at four rotor angles |
| Velocity | 4.90 min | 0.25 to 400 full-step/s plateaus, each path immediately reversed |
| Reversal | 5.00 min | 1 to 32 microstep triangular reversals at 0.2 and 1 Hz |
| Closing reference | 1.14 min | repeated two-cycle staircase and final standstill |

The default nominal test current is **360 mA RMS** because that is the current
used by the present rig firmware. It is not a claim about the motor rating.
Change it before generating the plan if needed; never exceed the motor or
driver limits.

Regenerate the executable plan and its figure:

```powershell
python scripts/sequence_design.py
python scripts/render_sequence.py

# Example with another tested operating current and no friction block
python scripts/sequence_design.py --nominal-current-ma 500 --no-reversal
python scripts/render_sequence.py --nominal-current-ma 500 --no-reversal
```

## Hardware assumptions

Run this campaign with the motor decoupled on the bench, or with a known load.
The encoder must measure the motor shaft directly. A stage/screw encoder mixes
motor error with coupling, screw, nut, bearing, and guideway error and cannot
support the separation below.

The pin map matches [`../CONNECTIONS.md`](../CONNECTIONS.md): STEP 5, DIR 6,
TMC UART 17/18, and encoder A/B/Z on 4/7/15. Confirm the encoder's electrical
levels before connection. The firmware assumes 0.11 ohm sense resistors.

Before a measurement:

1. Centre or decouple anything whose travel is limited. The warm-up and speed
   blocks rotate repeatedly; there are deliberately no stage soft limits here.
2. Fit the shaft encoder rigidly and load its previously measured eccentricity
   correction in the later analysis.
3. Record motor/driver identity, supply voltage, nominal test current, encoder
   CPR, load condition, and temperature sensor details.
4. Flash the firmware, then send `ZERO` only at the reference position.

## Run

Compile with the same ESP32-S3 board options used elsewhere in this repository:

```powershell
arduino-cli compile --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" `
  firmware/motor_identification
arduino-cli upload --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" `
  --port COM11 firmware/motor_identification
```

Validate the plan without hardware, then run it:

```powershell
python scripts/run_campaign.py --dry-run
python scripts/run_campaign.py --port COM11
```

Immediately before motion, the runner asks for confirmation that the motor is
decoupled or has a known load, the encoder is on the shaft, and rotation is
safe. During a run, `Ctrl+C` sends `ABORT`. Use `--yes` only in an already
interlocked automated setup.

Telemetry is recorded at 500 Hz for settled staircases, 1 kHz for reversals,
and 2 kHz for step responses and velocity plateaus. Encoder counting itself is
performed in hardware. Data rows are:

```
D,seq,t_us,action,event,command_u256,encoder,index_count,current_ma,mres
```

`command_u256` is commanded phase in 1/256-full-step units even when MRES
changes. `event=+1/-1` marks the exact direction of a commanded step edge;
ordinary periodic samples use zero. Lines beginning with `#` contain the mode
readback and action boundaries and remain valid CSV comments.

## Separation in analysis

First subtract the independently characterised encoder eccentricity. For the
200-step, 50-electrical-cycle motor, inspect harmonics at:

| Source | Mechanical-revolution orders |
|---|---|
| residual encoder mounting | 1/rev and 2/rev |
| current waveform mismatch/distortion | 50/rev, 100/rev, 200/rev |
| detent torque | primarily 200/rev |
| friction | forward/reverse offset and hysteresis rather than a fixed order |

For the static map, interpolate the forward and reverse passes onto identical
command positions. Their mean is the repeatable static angle error; half their
difference is directional hysteresis. In the current sweep, a detent component
in angle scales approximately as `1/current`, whereas phase-current waveform
angle error is approximately current-independent. Fit all current levels
together rather than interpreting any one harmonic as uniquely causal.

Use the marked step edges for natural frequency and damping fits. Use only the
constant-rate interior of each velocity action for mean lag and ripple; discard
the first and last electrical cycle. The closing two-cycle staircase and final
standstill quantify thermal drift and repeatability of the opening result.

The scripts generate and execute measurements; the physical campaign has not
yet been run on hardware. Verify direction, UART readback, encoder polarity,
and thermal limits in a short bring-up run before committing to the full plan.

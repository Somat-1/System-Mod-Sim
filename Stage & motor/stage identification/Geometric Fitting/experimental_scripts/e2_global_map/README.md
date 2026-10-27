# E2 global-map measurement

This folder contains the generated target design, executable measurement plan,
sequence renders, terminal logger, and dedicated ESP32-S3 firmware for
Experiment E2. E2 measures the direction-averaged geometric error over the
full ISO validation span for model B. Homing and end-to-end warm-up remain a
separate workflow.

## Precondition

Before arming the firmware, the stage must be homed, warmed, energized at the
ISO v2/E1 driver settings, physically parked at `x_ref = 61.35 mm`, and being
recorded by the interferometer. Pressing `p` explicitly assigns that physical
position to the firmware coordinate. Do not arm elsewhere.

## Generate and render

```bash
python3 e2_sequence.py
python3 render_e2_sequence.py
```

The generator reads the executed E0 target list from
`../../../iso230-2/sequence/iso_targets_seed2302_rom122.695.csv`. It divides
2–120 mm into 400 equal strata and selects one 64-microstep-grid target from
each, excluding every point within ±64 microsteps of an E0 target. Seeds start
at 5200 and increment until all 50 screw-phase bins contain at least five
targets. The accepted deterministic design uses seed **5315**.

Generated outputs:

- `sequence/e2_targets_seeded_grid64.csv` — 400 targets and screw-phase bins;
- `sequence/e2_measurement_plan.csv` — complete controller action plan;
- `sequence/e2_plan_summary.json` — config, hash, seed, timing, and checks;
- `firmware/e2_global_map/e2_config.h` — compiled settings and target array;
- `rendered/e2_global_map_sequence.png` — full trajectory and time budget;
- `rendered/e2_target_design_and_rmt.png` — coverage, phase balance, mirrored
  timing, and representative RMT command.

The nominal sequence is 38.6 minutes: 400 upward stops, the turnaround
reference, 400 downward stops, and seven reference records in total.

## Endpoint safety adaptation

The supplied procedure requests a 3 mm approach/overrun around a 2–120 mm
target span, i.e. turn points at −1 and 123 mm. Those positions are outside the
proven E1 soft envelope of 2.000–120.695 mm and the measured travel. The
executable sequence therefore uses **2.000 mm** as the lower approach and
**120.695 mm** as the upper turn point while retaining the requested 2–120 mm
sampling span. This reduced endpoint clearance is recorded in the JSON summary,
both figures, and every firmware session header. Do not widen the limits without
a separate travel/switch validation.

## Firmware

`firmware/e2_global_map/e2_global_map.ino` preserves the E1 command path:

- 1 MHz ESP32 RMT STEP generation with six rotating symbol buffers;
- trapezoidal or triangular timing precomputed before each transmission;
- 128-pulse RMT transactions queued without CPU pulse timing;
- no TMC UART or serial telemetry during a pulse train;
- MSCNT and driver status read only after each record window;
- switch, soft-limit, incomplete-move, RMT, and driver-fault handling;
- safe park through the fixed-side reference approach when the coordinate is
  still trustworthy.

Compile and upload using the same ESP32-S3 board package as E1:

```bash
arduino-cli compile \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" \
  firmware/e2_global_map

arduino-cli compile --upload -p /dev/cu.usbmodem3101 \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" \
  firmware/e2_global_map
```

Then run `python3 e2_terminal.py [port]`. Commands are:

- `p` — confirm the precondition and arm the coordinate at x_ref;
- `r` — execute INIT → DATUM → UP → TURN_REF → DOWN → END_REF → PARK;
- `s` — status; `t` — target summary; `0` — release the motor and invalidate
  the coordinate; Ctrl+C — request abort; Ctrl+] — close the terminal.

At the turnaround, the controller estimates the down pass, end references, and
summary margin. A projected finish beyond 58 minutes skips the down pass, marks
the data set incomplete for direction averaging, records END_REF, and parks.

## Logging boundaries

The terminal stores the complete timestamped serial stream in `data/`. The
session header logs the config hash, seed, all 400 targets, screw-phase bins,
motion settings, soft limits, and endpoint adaptation. Every `# STOP` line
contains state, map/reference kind, pass, direction, target index, integer
microsteps, millimetres, screw-phase bin, move start, arrival, record start,
record end, and post-window MSCNT.

As with E1, the ESP32 cannot read the IDS interferometer or environmental
sensors. The IDS/host analysis remains responsible for stop-window mean/std/
slope/sample count, saturation detection, measured-reference drift checks,
air-temperature validity, and the final `stops.csv`, `events.csv`, `env.csv`,
and `summary.json` analysis products. Stops are flagged after acquisition; none
are repeated inline because repetition would break the mirrored timing.

## Acceptance summary

Before model-B fitting, require at least 95% valid stops, paired directions for
at least 95% of targets, no watchdog abort, consistent reference visits,
recording time at or below 60 minutes, and air-temperature change below 0.5 K.

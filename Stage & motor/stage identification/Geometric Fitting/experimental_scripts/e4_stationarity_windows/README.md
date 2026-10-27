# E4 stationarity-window measurement

This package contains the deterministic target design, complete action plan,
ESP32-S3 firmware, lossless terminal logger, and academic renders for Experiment
E4. E4 estimates screw-periodic harmonics locally at seven positions to test the
stationarity assumption in model C.

## Generate and render

From this directory:

```bash
python3 e4_sequence.py
python3 render_e4_sequence.py
```

The generator produces:

- `sequence/e4_targets_grid192.csv` — all 357 target coordinates and phases;
- `sequence/e4_measurement_plan.csv` — the complete controller action plan;
- `sequence/e4_plan_summary.json` — settings, hash, checks, counts, and timing;
- `firmware/e4_stationarity_windows/e4_config.h` — generated firmware constants;
- `rendered/e4_stationarity_windows_sequence.png` — trajectory and time budget;
- `rendered/e4_window_design_and_rmt.png` — spatial, phase, traversal, and RMT detail.

The seven windows are 7–13, 24–30, 41–47, 58–64, 75–81, 92–98, and
109–115 mm. Each contains 51 targets at 0.12 mm spacing. Its first 50 targets
cover every screw-phase bin once; the endpoint repeats the starting phase after
exactly three screw revolutions.

## Precondition and sequence

E4 does not home or warm the stage. Before pressing `p`, complete the established
homing and end-to-end warm-up, leave the motor energised, start the IDS and
environmental recordings, and physically park at `x_ref = 61.35 mm`. Pressing
`p` assigns that physical point to the controller coordinate; do not arm at any
other position.

Command `r` executes:

`INIT → DATUM×3 → W1↑↓ → REF → … → W7↑↓ → END_REF×3 → PARK`

There are 714 analysis stops and 12 references, or 726 recorded stops in total.
The nominal duration is 36.0 minutes. The fixed reference is approached from
58.35 mm. A projected-finish check after W4 stops the remaining windows if the
run would exceed 58 minutes, but still attempts end references and a safe park.

## Firmware and logger

`firmware/e4_stationarity_windows/e4_stationarity_windows.ino` preserves the
E1 command path and hardware contract:

- ESP32-S3, TMC2209, SpreadCycle, 400 mA run/hold;
- 1 MHz RMT STEP output with six rotating 128-symbol buffers;
- precomputed trapezoidal/triangular pulse timing and no serial/UART work inside
  a pulse train;
- post-record MSCNT/driver telemetry, switch and soft-limit checks, controlled
  abort, incomplete-motion detection, and conditional safe parking.

Compile with the same board package and options used for E1:

```bash
arduino-cli compile \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" \
  firmware/e4_stationarity_windows
```

Run `python3 e4_terminal.py [port]`. Commands are `p` to arm at x_ref, `r` to
run, `s` for status, `t` for the target summary, `0` to release the motor and
invalidate the coordinate, Ctrl+C to request a controlled abort, and Ctrl+] to
close the terminal. Raw timestamped logs are written under `data/`.

The ESP32 cannot read IDS or environmental values. Host-side processing remains
responsible for stop means/std/slopes/sample counts, saturation and settling
flags, reference drift correction, environmental validity, harmonic fits, and
the final stationarity verdict. The full procedure is preserved at
`../../instructions/E4_stationarity_windows_procedure.md`.


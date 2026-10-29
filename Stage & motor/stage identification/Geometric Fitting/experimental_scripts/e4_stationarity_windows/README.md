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
run, `h` to home both ends (the step_timing_test routine: datum at the negative
switch), check the travel holds the soft limits, centre at x_ref from 58.35 mm
and arm in place of `p`, `g` for `h` then `r`, `s` for status, `t` for the target summary, `0` to release the motor and
invalidate the coordinate, Ctrl+C to request a controlled abort, and Ctrl+] to
close the terminal. Raw timestamped logs are written under `data/`.

The ESP32 cannot read IDS or environmental values. Host-side processing remains
responsible for stop means/std/slopes/sample counts, saturation and settling
flags, reference drift correction, environmental validity, harmonic fits, and
the final stationarity verdict. The full procedure is preserved at
`../../instructions/E4_stationarity_windows_procedure.md`.

## Completed-run processing

The completed 27 October 2026 run is processed reproducibly from:

- `data/E4_Experimental_Run.csv` — the 1 kHz TwinCAT/IDS capture;
- `data/e4_terminal_20261027_155838.log` — controller events and stop windows;
- `e4_run.py` — scope/log clock alignment and stage-coordinate reconstruction;
- `render_e4_findings.py` — validity checks, E4 and E1 harmonic refits, complex
  profile-model selection, tables, and publication figures.

Because the local Homebrew Python 3.14 installation needs its matching Expat
library at runtime, regenerate all findings from this directory with:

```bash
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib \
  /opt/homebrew/bin/python3 render_e4_findings.py
```

Generated data products are under `analysis/`:

- `e4_validity_table.csv` — acquisition and step-count checks;
- `e4_t1_window_harmonics.csv` — E4 and consistently refitted E1 windows;
- `e4_t2_profile_models.csv` — P0–P3 comparison by χ² and BIC;
- `e4_findings_summary.json` — machine-readable conclusions and conventions.

The raw-window evidence figure and four appendix figures are written as 320 dpi
PNG files under `rendered/experimental_findings/`. The three presentation-ready
results figures are grouped under `rendered/presentation/`:

- `e4_figure1_geometry_fold.png` — geometry fold, windows told apart by line style;
- `e4_figure1_geometry_fold.gif` — the same fold animated: the highlight steps
  through the seven windows from motor end to support bearing (~17 s loop);
- `e4_figure2_harmonic_profile.png` — k = 1/2 amplitude and phase profiles;
- `e4_figure3_reversal.png` — reversal profile and direct co-variation test.

The processing currently places the motor end at 0.000 mm and the support
bearing at 122.770 mm, inferred from the logged homed travel. These two physical
locations must be confirmed on the stage before their labels or physical
profile-model parameters are treated as final.


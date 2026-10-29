# E1 spectral-windows measurement

This folder contains the generated sequence, sequence render, terminal and
dedicated ESP32-S3 firmware for Experiment E1. E1 is measurement-only: homing
and end-to-end warm-up are completed beforehand by the separate warm-up
workflow.

## Precondition

Before arming this firmware, the stage must be:

1. homed and warmed endstop-to-endstop;
2. parked physically at `x_ref = 61.35 mm`;
3. still energized at the ISO v2 driver settings; and
4. being recorded by the interferometer.

The sketch cannot read the IDS interferometer. Pressing `p` is the operator's
explicit confirmation that the physical stage and the firmware's absolute
coordinate agree. Do not arm at any other position.

## Generate and render

```bash
python3 e1_sequence.py
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib python3 render_e1_sequence.py
```

The generator creates a plain 0.12 mm grid: 192 microsteps, exactly three
64-microstep motor cycles. Consequently, every target within a window has the
same motor-cycle phase. The motor error becomes one constant window offset and
is removed by the existing per-window offset subtraction. No seed or random
jitter is used.

The resulting 101 points per 12 mm window resolve approximately through the
8th screw harmonic; the analysis only interprets harmonics through the 5th.

It writes:

- `sequence/e1_targets_grid192.csv` — 101 integer-microstep targets in each
  of W1, W2 and W3;
- `sequence/e1_measurement_plan.csv` — the measurement-only controller plan;
- `sequence/e1_plan_summary.json` — full settings, configuration hash, counts,
  duration and validation checks;
- `firmware/e1_spectral_windows/e1_config.h` — the same targets and settings
  compiled into the sketch.

The overview is `rendered/e1_spectral_windows_sequence.png`. A second render,
`rendered/e1_window_rmt_detail.png`, expands one window into its pre-approach,
two 101-target passes, overrun/reversal, 33 mm total path, one 0.12 mm velocity
profile, per-stop timing and the 128 + 64 RMT transaction split. The plan
contains 606 analysis stops instead of 726. With the current 2.0 s settling
time it is estimated at 33.8 minutes including a two-minute pilot-review
budget, leaving about 26 minutes against the 60-minute limit.

## Firmware

`firmware/e1_spectral_windows/e1_spectral_windows.ino` uses the ESP32 RMT
peripheral for every motion. Its pulse-profile generator is the same buffered
1 MHz RMT implementation used for the 2/5/8 mm/s bandwidth test:

- trapezoidal/triangular timing is precomputed into RMT symbols;
- chunks are queued without per-pulse CPU timing;
- no TMC UART access or serial telemetry occurs inside a pulse train;
- MSCNT and driver status are read only after each record window.

Compile and upload:

```bash
arduino-cli compile \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" \
  firmware/e1_spectral_windows

arduino-cli compile --upload -p /dev/cu.usbmodem3101 \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" \
  firmware/e1_spectral_windows
```

Then run `python3 e1_terminal.py [port]`. The terminal stores the complete
session stream in `data/`.

Firmware commands:

- `p` — confirm the precondition and assign the current physical position to
  x_ref;
- `r` — execute INIT → DATUM → PILOT → operator gate → W1 → REF → W2 → REF →
  W3 → END_REF → PARK;
- `h` — home both ends (the step_timing_test routine: datum at the negative
  switch), check the travel holds the soft limits, centre at x_ref from
  58.35 mm and arm; replaces the manual `p` precondition;
- `g` — `h`, then `r` with the pilot gate continuing by itself;
- `c` / `a` — continue or abort at the pilot gate;
- `s` — status; `t` — target summary; `0` — release the motor and invalidate
  the coordinate; Ctrl+C — abort.

The time watchdog runs after W2. If the projected finish exceeds 60 minutes,
the firmware retains the W1/W2 data, skips W3, records END_REF, parks and marks
the session incomplete. A switch, soft-limit or RMT fault disables automatic
parking because the absolute position can no longer be trusted.

## Logging boundaries

Each `# STOP` record contains state, window, pass, direction, target index,
integer-microstep target, move-start, arrival, record-start, record-end and
post-window MSCNT. `# TARGET` lines put the complete target list in the session
header. The host/analysis layer remains responsible for IDS signal checks,
measured position/reference watchdogs, environmental data, and stop-window
mean/std/slope/sample-count calculations.

Set `sync_gpio` in `e1_config.json` only after a dedicated line is wired to the
Beckhoff. Its default is `-1`, so synchronization uses the event clock map.

## Completed-run findings

The publication findings are regenerated from the aligned 1 kHz IDS capture and
ESP log with:

```bash
python3 render_e1_findings.py
```

The analysis removes the Step A gain (`−25.0866 ppm`) and a direction-averaged
offset per window, folds on absolute commanded screw angle, tests harmonics
1–8 against a 1% family-wise false-alarm amplitude, and writes:

- `rendered/experimental_findings/e1_findings_figure1_raw_windows.png` — corrected raw evidence;
- `rendered/experimental_findings/e1_findings_figure2_folded_revolution.png` — folded direction
  average and reversal;
- `rendered/experimental_findings/e1_findings_figure3_harmonic_spectrum.png` — amplitude spectrum;
- `rendered/experimental_findings/e1_findings_table_harmonics.png` — main-text Table 1
  (harmonics 1–4), including the k = 1 W1/W2 amplitude-ratio contrast;
- `rendered/experimental_findings/e1_findings_appendix_table_harmonics.png` — full
  harmonic table (harmonics 1–8) for the appendix;
- `rendered/experimental_findings/e1_findings_appendix_validity.png` — references, record quality,
  and harmonic-fit residuals;
- `analysis/e1_findings_summary.json` and `analysis/e1_harmonics_table.csv` —
  machine-readable results and method metadata.

Measured observations use markers only; continuous lines denote fitted models
or reference/acceptance levels. Window colours are consistent throughout.
The folded figure and amplitude spectrum are the key main-text figures. The
validity figure remains in the appendix, but panel (d) should be cited in the
main text: coherent position-dependent structure remains after the order-2
harmonic fit, so the residual is not white noise. The repeat-phase σ used for
the false-alarm threshold therefore includes unmodelled slow structure and
makes the harmonic-selection test conservative.

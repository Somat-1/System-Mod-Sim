# Geometric Fitting

This project isolates the position-dependent geometric component of the stage's
ISO 230-2 v3 error. The first completed stage is the linear gain fit described
in `instructions/G_step_A_linear_processing.md`.

## Structure

- `instructions/` — preserved analysis specification supplied with the task.
- `data/` — source inventory and data-availability notes; source data remain in
  the existing `iso230-2` project to avoid duplication.
- `scripts/` — reproducible processing and rendering code.
- `experimental_scripts/` — hardware-experiment sequence definitions, separated
  by experiment; E1, E2, and E4 are in `e1_spectral_windows/`,
  `e2_global_map/`, and `e4_stationarity_windows/` respectively.
- `results/` — machine-readable fit results and per-target table.
- `rendered/` — A4 portrait PNG summary. No PDF is generated.

## Run

From this directory:

```bash
python3 scripts/process_iso_v3_step_a.py
```

On the current Homebrew installation, Python 3.14 requires its matching Expat
library to be selected explicitly:

```bash
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib python3 scripts/process_iso_v3_step_a.py
```

The current analysis reads:

- `../iso230-2/analysis/iso230_2_v3_results.json`
- `../iso230-2/data/terminal_20261021_183115.log`

## Experimental sequences

E1, E2, and E4 are measurement-only workflows. Their generated plans, ESP32-S3
firmware, terminal loggers, and academic sequence renders are self-contained in
`experimental_scripts/e1_spectral_windows/` and
`experimental_scripts/e2_global_map/`, and
`experimental_scripts/e4_stationarity_windows/`. The supplied E2 design
specification and the compiled E4 specification are stored as
`instructions/E2_global_map_procedure.md` and
`instructions/E4_stationarity_windows_procedure.md`.

The raw `../ISO230-2_v3.csv` capture is absent from the workspace. Therefore,
the present run exactly reconstructs the per-target weighted fit, direction
fits, leave-one-target-out sensitivity, error-power shares, reversal trend and
predicted ISO metrics from the saved aggregate statistics. It does not claim a
new drift-corrected fit, per-run slopes, or raw-window settling analysis.

## Important v3 adaptation

The attached instructions contain a v2 note referring to 221 windows. The
actual v3 controller log ends with `ISO_DONE,runs=5,stops=240`: 24 targets × 2
directions × 5 runs. All 240 ISO stops are represented by the saved statistics,
with five samples per direction at every target.


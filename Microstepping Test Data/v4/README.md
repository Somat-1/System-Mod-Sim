# v4 — settling-error characterization

## Purpose

v3's D-block campaign varies commanded **rate** at a fixed move structure
and studies *tracking* — the position error while the move is happening.
v4 holds rate fixed and varies commanded **distance** instead, and looks
at what happens *after* each move ends: overshoot, ringing, and creep
back to the final position. Same hardware (EVO dedicated controller,
protocol in `../v2/docs/Stepper Motor Controller Command list.pdf`), same
axis (X), same proven runner infrastructure as
`../v3/scripts/run_identification_dedicated_controller.py`.

## Experiment design

- **Test distances** (`TEST_DISTANCES_FULL_STEPS`): `1, 2, 4, 8, 16, 32,
  64, 128, 256, 512, 1000` full steps — a doubling sweep from
  single-microstep scale up to several revolutions. 1 full step = 1.8° =
  10 µm of stage travel (`L=2 mm` lead pitch). The small-distance end is
  the primary interest (settling error is proportionally largest there);
  the large-distance points are included for context/completeness.
- **Move speed**: fixed at `SETTLE_MOVE_FULL_STEPS_S = 250` full-steps/s,
  burst-style (`SS` with `MIN=MAX`, no ramp — same convention as v3's
  reference moves), for every distance. Holding speed fixed makes
  distance the sole independent variable.
- **Per-distance test** (`run_settling_test`): out-and-back from the
  shared origin — move `+D`, dwell `SETTLE_DWELL_S = 30 s` (records
  settling), move `-D`, dwell `30 s` again (records return settling),
  then an `assert_origin` check. Both directions are exercised every
  time, so direction-dependent settling (e.g. backlash-like asymmetry)
  is visible without a separate test.
- **Dwell duration (30 s)**: long enough to capture both the fast
  structural ringdown (~183–211 Hz, ζ≈0.2–0.5 per
  `../v3/Parameter Optimization/calibration_bracketing/step4c_ringdown_fit.py`
  — decays within ~100 ms) and slow LuGre presliding/creep relaxation
  (v3's `C` block used 60 s `creep_record` dwells for the same reason;
  halved here since this campaign repeats the dwell many more times).
- **Segmentation**: every distance test is preceded by a `run_marker`
  call with a unique, monotonically increasing amplitude — the exact
  device v3 uses for the same purpose. Two independent ways to segment
  the resulting IDS trace afterward: (1) the CSV log's own `block`/
  `label` columns, exact; (2) for pure-trace analysis with no CSV
  cross-reference, a marker's short (~1 s) reverse dwell contrasts
  sharply with every settling block's much longer (30 s) dwell.
- **Configuration sweep**: `run_campaign` loops all 3 MRES values × 2
  current levels (6 configurations total, ~12 min each ≈ 72 min for the
  full campaign), identical structure to v3.

## Files

- `scripts/run_settling_dedicated_controller.py` — the controller
  driver. `--dry-run` validates the whole sequence and produces a CSV
  log with no hardware; `--execute --port COMx ...` runs it for real
  (same CLI conventions as v3: homing options, travel-limit guards,
  `DS`-based `wait_ready` polling, `assert_origin` checks after every
  block, safe shutdown on any failure or Ctrl-C).
- `scripts/plot_planned_sequence.py` — reads the most recent `--dry-run`
  log and reconstructs the ideal commanded-position-vs-time preview for
  one representative configuration (`run_index=1`), with block shading,
  entirely analytically (no hardware, no dynamics model — this is a
  preview of what will be *commanded*, not a prediction of how the stage
  will actually respond).
- `rendered_assets/planned_settling_sequence_preview.png` — that preview.
- `data/hardware_runs/` — dry-run and (eventually) live CSV logs, same
  schema as v3 (`CsvEventLog.FIELDNAMES`).

## Status

Sequence designed, implemented, and validated in `--dry-run` only. Not
yet run against real hardware. Before `--execute`: confirm the 30 s
per-distance dwell is enough once a first live run is inspected, confirm
travel limits/homing arguments for the actual bench setup, and decide
whether `REPEATS_PER_DISTANCE` (currently 1) should be increased for
repeatability statistics.

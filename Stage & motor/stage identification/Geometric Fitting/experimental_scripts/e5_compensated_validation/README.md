# E5: compensated validation (ISO 230-2, paired)

Does the frozen geometric and reversal model, applied as a **command
pre-correction**, reduce the stage's positioning error on an independent
session, and by as much as the offline replay predicts? The procedure is
`E5_validation_procedure.md` (supplied 2 Nov). This folder holds:

- the frozen model;
- the command table;
- the session plan;
- the renders of both the sequence and the planned analysis.

Status: **sequence generated and validated, renders done. No firmware yet, no data.**

## Files

| Path | Contents |
|---|---|
| `e5_config.json` | the session definition; the SHA-256 of this file is recorded in the summary |
| `e5_geometry.py` | a standalone evaluator of the frozen model (no dependency on the E2 pipeline) |
| `e5_freeze_model.py` | refits E2's selected model and writes `model/geom_model_E2_frozen.json` |
| `e5_make_table.py` | targets, compensation law, plan, and the procedure's §4.3 checks |
| `render_e5_sequence.py` | the three figures in `rendered/` |
| `model/` | the frozen model (sha256 `98ebd4d24c77…`) and E2's fitted points |
| `sequence/e5_command_table.csv` | every stop in execution order, with the columns from §4.1 plus `recording`, `x_nom_mm` and `G_rel_um` |
| `sequence/e5_measurement_plan.csv` | the full action list (moves, stops, state changes, the gate) |
| `sequence/e5_b3_targets.csv` | the 32 off-grid targets and their assigned electrical offsets |
| `sequence/e5_plan_summary.json` | check results, timing, hashes, seed, and spot checks |
| `sequence/e5_expected_metrics.json` | the offline-replay prediction for U and C |

```sh
export DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib
/opt/homebrew/bin/python3 e5_freeze_model.py        # reproduces the same hash
/opt/homebrew/bin/python3 e5_make_table.py          # --check-only to write nothing
/opt/homebrew/bin/python3 render_e5_sequence.py
```

## The frozen model

E2's selected model is M_C: H = 4, K_s = 2, K_C = 8, λ = 1, with drift
variant D2. It is refitted from the E2 capture through
`render_e2_findings.py`, and the script asserts that the refit reproduces E2's
residual and β₀ exactly. All its numbers are then stored in JSON. Before the
file is written, the standalone evaluator is checked against E2's `predict()`
on a dense grid, and the two differ by 0.0 µm.

- RMS residual 0.369 µm; cross-validation error 0.434 µm.
- Gain −29.81 ppm.
- β₀ = −1.757 µm.
- Domain 2.280–119.920 mm.

The fit uses E2 data only. Nothing from E3 enters it.

The datum a₀ is the model's intercept. The correction acts on Ĝ − a₀.

## Compensation law

    x_c = x* − [Ĝ(x*, θ(u_nom)) − a₀] − (d/2)·β₀,     u_cmd = round(x_c · 1600)

- Ĝ is evaluated at the nominal target and its screw angle, as §2 specifies.
- Pre-approach and overrun positions come from the nominal targets.
- Uncompensated (U) rows command `u_nom` unchanged.

All 14 checks pass (`e5_plan_summary.json`):
- **Correction size:** corrections span −4.38 to +4.38 µm, well inside the 6 µm limit.
- **Quantisation:** |q| ≤ 0.310 µm, with an RMS of 0.164 µm.
- **Sign check:** up minus down is between 1.250 and 1.875 µm, which is +1.757 ± 0.625 µm as required. The up-approach is sent higher.
- **U runs:** they are identical, command for command, to the E0 v4 plan.
- **Spot checks:** three targets (5.0, 63.12 and 117.5 mm) are re-derived step by step in the summary. For example, at 5 mm Ĝ − a₀ = +2.860 µm, so the up command is 4.998019 mm = 7997 µsteps and the down command is 7994 µsteps.

## Where this departs from the procedure, and why

1. **Two IDS recordings instead of one 60-minute window.** The procedure
   estimates about 3 s per stop. E0 v4's real cycle is **5.05 min per ISO run**,
   about 6.3 s per stop, set by 2 s settle, 0.5 s record and the move between
   targets about 4.9 mm apart. As a single window the session takes
   **62.2 min**, against a budget of 52 min (60 minus an 8-minute margin).
   Cutting B3 to one repetition, which the procedure lists as the first thing
   to drop, still does not fit. The plan therefore splits after B1c:

   | Recording | Contents | Duration |
   |---|---|---|
   | R1 | B1a, B2 (10 runs), B1b after run 5, B1c | 51.3 min |
   | R2 | B1d, B3 (2 reps), B1e | 10.8 min |

   Both recordings share one homing. At the split the stage holds at x_ref
   behind an operator gate, so the step-count origin and θ carry over. B2 is
   left untouched, as the procedure demands. Each recording has its own
   reference anchors for drift correction. Setting
   `split_recording_before_b3: false` shows the single-window plan failing its
   time check. R1's 51.3 min is only 0.7 min under budget; the E2 timing model
   was accurate to 0.1 min (38.6 min estimated, 38.7 min actual).

2. **E0 v4's per-run reference visit is kept.** Every ISO run opens with E0
   v4's own x_ref visit, so the U runs match E0 v4 exactly. That gives drift
   anchors at every run in addition to B1a, B1b and B1c.

3. **B3 offsets are permuted, not placed in x order.** The 32 offsets −32,
   −30, …, +30 cover the electrical cycle once. If target k took offset
   o_k = −32 + 2k in x order, electrical phase would track position, and so
   track the slow error s(x). A seeded permutation (seed 5500) assigns them
   instead.
   - B3 spans 5.0–117.5 mm, the ISO span, so no target relies on spline
     extrapolation.
   - Every B3 target is more than one grid point away from every E0, E2 and E3
     target.

4. **The ISO targets are not on the 64-µstep grid.** §1 assumes they move
   "only a few µsteps off the electrical grid after correction". In fact E0
   v4's targets are already at scattered electrical offsets (23541 µsteps is
   offset 21, for example), so the ISO stops sample phase too, though
   unevenly. B3 is still the clean test because its offsets are uniform and
   independent of x (figure 2d).

## The prediction does not meet the acceptance criteria in §6

I replayed the planned session offline. For each target and direction, the
predicted compensated deviation is E0 v4's measured deviation plus the
correction the table actually commands, with quantisation included. The
result is in `sequence/e5_expected_metrics.json` and figure 3.

| | E | M | A | mean B |
|---|---|---|---|---|
| U (E0 v4, measured) | 8.52 | 6.88 | 9.25 | −1.82 |
| C (replay) | **2.44** | **2.04** | 3.56 | −0.07 |
| §6 pass line | ≤ 1.96 | ≤ 1.735 | | ≤ 0.3 in magnitude |

On E0 v4, the frozen model is predicted to **fail M_C and E_C** while still
cutting M by a factor of 3.4. The "offline prediction" of M 1.335 µm and
E 1.560 µm that the pass lines are built on cannot be reproduced from
anything in the repository:

- E2's own validation on v3 gives M 1.22 µm and E 1.48 µm;
- the E3-repeat findings give M 1.96 µm and E 2.30 µm on v4, with a
  posterior M of 2.03 µm (5th–95th percentile 1.86–2.22 µm).

The baseline the replay predicts on is the same stage and targets as E5's U
runs. So the criteria should be settled before the run, as §6 itself says:
either re-base them on the v4 replay, or record where 1.335 and 1.560 µm come
from.

## Before running

- **Firmware.** The firmware is not written yet. It should be the E0 v4 RMT
  state machine executing `e5_command_table.csv` row by row (§4.2). It needs
  to:
  - echo `row` and `u_cmd` at every stop;
  - print the model and config hashes and the B3 seed at sequence start;
  - honour the gate at the recording split.
- **Homing phase.** The E3-repeat analysis found a homing reference shift of
  about 2° between sessions (δ₁ ≈ −1.99°, confirmed by MSCNT). That is exactly
  the "k = 1 remains, phase shifted" failure mode in §7. Logging MSCNT at both
  home switches, as §3.2 asks, is what makes it diagnosable.
- **Bench dry run.** Do the dry run from §4.3, with the stage disconnected or
  over a short range, and check that the firmware's `u_cmd` echoes match the
  table.

## Figures

- `rendered/e5_sequence_windows.png`: the whole session on one time axis,
  with each window bracketed and labelled with the validation it serves:
  - B2, the paired U vs C ISO runs, for §5.2–5.4;
  - B3, off-grid, for §5.5;
  - the B1a–e drift anchors, for §5.1;
  - the two IDS recordings and the pause between them.
- `rendered/e5_command_correction.png`: what the C runs command at each ISO
  target, for both directions. The up-minus-down gap is the reversal
  correction.
- `rendered/e5_expected_result.png`: the replay prediction, showing U against
  C along the travel, and E, M and |mean B| against the §6 pass lines.

The numeric pre-run checks (sign, quantisation, electrical-offset coverage,
spot checks) are in `sequence/e5_plan_summary.json` rather than in a figure.

# What this initial-testing pass actually established

Synthesis across `sensitivity_dynamic_fit.py`, `model_parameter_optimization.py`
round 1, and round 2. Written as a capstone, not a new experiment.

## 1. The pipeline itself is validated

Real measured target (`MRES_1/fitted_blocks.npz`, an actual v6 EL5101
capture) → real 15-state nonlinear Rev 4.2 LuGre model → real
`scipy.optimize.least_squares` run, twice, with the second warm-started
from the first. Cost dropped monotonically and reproducibly
(0.904 → 0.690 → 0.576 µm RMS-to-measured). This is not a demonstration —
it is a working fitting pipeline that measurably improves against real
hardware data.

## 2. What is NOT yet a calibration result

`sigma0_way`, `Fc_way`, `Fs_way` do **not** have trustworthy fitted values
from this data:

- `Fc_way`/`Fs_way` moved in exact lockstep (their ratio stayed pinned at
  1.5, baseline's own ratio) across two rounds, both times running to
  whichever upper bound was set (5×, then 20× baseline) with cost still
  falling. This is the signature of one effective degree of freedom
  ("way-port friction level"), not two independently-identifiable
  parameters — the classic collinearity `calibration_bracketing/
  SENSITIVITY_ANALYSIS.md` already warns about for this same parameter
  pair, now reproduced independently on different data (a v6 oscillation
  block instead of v3 D-rate legs).
- `sigma0_way` is not stable either: its round-1 optimum (7.32e6) was not
  its round-2 optimum (6.86e6) — it moved *because* Fc/Fs got more room,
  meaning its apparent best value is conditional on wherever Fc/Fs happen
  to be capped, not a property of the data on its own.
- The values the search wants (`Fc_way`→60 N, `Fs_way`→90 N) are
  physically implausible for this port — 100–200× the *nut* port's own
  baseline (0.5 N) on the same rig.

**Conclusion: this 3-parameter subspace is the wrong subspace.** Chasing
it further (a third bounds-widening round) would very likely repeat the
same pattern at a still-larger, still-implausible scale rather than find
a real optimum.

## 3. What IS a solid, reproducible finding

- **The baseline model overshoots.** On both step directions, baseline
  simulated peaks exceed measured by 4–5 µm (forward 16.7 vs 11.3–11.7;
  return −6.7 vs −1.9 to −2.2). Increasing way-port friction moves both
  peaks toward measured, on both polarities, non-trivially (round 2 got
  within ~0.2–1.1 µm of measured on both) — so the *direction* of the fix
  (more damping/friction needed) is real, even though the specific
  parameters carrying it are not identifiable from this alone.
- **A ~0.65–0.67 µm plateau offset that friction tuning structurally
  cannot touch.** Measured settles at 9.326 µm / −0.025 µm against a
  10.0/0.0 µm command; baseline AND both optimized rounds settle at
  *exactly* 10.000/0.000 µm regardless of `sigma0/Fc/Fs_way` — because at
  MRES 1 every commanded position sits exactly on this model's detent
  minimum, so friction parameters have zero steady-state effect there.
  This was checked directly (not inferred) and held unchanged, to the
  fourth decimal, across both rounds. It must come from something outside
  this model's friction ports entirely — sensor/encoder scale calibration,
  lead-screw tolerance, or an unmodeled backlash/asymmetry — and it puts a
  hard floor under achievable RMS no matter how well the transient gets
  matched.
- **Noise cannot substitute for either finding above.** Verified
  numerically (not just asserted): adding the block's own real sensor
  noise floor (0.022 µm std) to the simulated trace left cost unchanged;
  larger noise only made it worse. Both gaps above are systematic bias,
  not measurement scatter.

## 4. Where the model's own documentation already pointed

`Rev 4.2/model_parameters.json`'s `sigma0_way_caveat` flags `sigma1_way`
(776.7 N·s/m) as **not recomputed and inconsistent with the target-zeta=0.7
rule** since `sigma0_way` was last revised. `sigma1_way` is bristle
*damping* — mechanically the right lever for "the transient rings down
faster than predicted," which is exactly what we're seeing — rather than
`Fc_way`/`Fs_way`, which set static/Coulomb friction *levels*, not decay
rate. This wasn't something this session discovered independently; the
model's own notes already said this parameter was stale. The optimizer
reaching for Fc/Fs instead (the only knobs it was given) and running them
to an implausible extreme is consistent with it compensating for a
missing damping correction it wasn't allowed to make directly.

## Bottom line

Three things are now established with reasonable confidence: (1) the
fitting pipeline works and produces real, verifiable improvement; (2) the
specific three friction-level parameters tested are not independently
identifiable from a single oscillation block and shouldn't be trusted as
calibrated; (3) there is a genuine, friction-independent ~0.65 µm
plateau-offset problem that no amount of parameter fitting in this
subspace will resolve. None of this yet constitutes a calibrated model —
it constitutes a diagnosis of what would need to change (free `sigma1_way`
instead of `Fc_way`/`Fs_way`; separately resolve the plateau offset) before
a calibration attempt on this data would be meaningful. That diagnosis,
not a parameter table, is this initial-testing pass's actual output.

# Rev 4.2 way-port parameter optimization — initial run

**Script**: `model_parameter_optimization.py` | **Raw output**:
`model_optimization_results.json` | **Plot**: `plots/model_fit_overlay.png`

Fits three Rev 4.2 LuGre way-port parameters (`sigma0_way`, `Fc_way`,
`Fs_way`) against the MRES 1 oscillation block's first 3 of 15 cycles
(6 edges, 5.78 s), using the same measured target (`MRES_1/fitted_blocks.npz`)
as `sensitivity_dynamic_fit.py`. `sigma1_way`, `sigma2_way`, `vs_way` held at
Rev 4.2 baseline. Ran for real: `scipy.optimize.least_squares`, 25 accepted
function evaluations, 608 s wall time.

## Why only 3 parameters / 3 cycles

One model evaluation (a 6-edge piecewise Radau integration of the full
15-state Rev 4.2 model) costs ~10 s even after loosening solver tolerances
from the `generate_stepping_rev42.py` convention (`rtol` 1e-6→1e-4, `atol`
tightened only on the bristle states). A 6-parameter fit over all 15 cycles
was timed first (49 s/eval at tight tolerance) and projected to run for
hours — this run is deliberately the largest version that finishes as an
actual test rather than a job started and abandoned.

## Result

| | baseline | optimized |
|---|---|---|
| RMS-to-measured | 0.9041 µm | 0.6897 µm (**−23.7%**) |
| `sigma0_way` (N/m) | 3.45e6 | 7.319e6 (+112%) |
| `Fc_way` (N) | 3.0 | **15.0** (+400%, = upper bound) |
| `Fs_way` (N) | 4.5 | **22.5** (+400%, = upper bound) |

`least_squares` reported `success=True` (`xtol` termination), **but treat
that with caution**: `Fc_way` and `Fs_way` both landed exactly on the
artificial upper bound I set (5× baseline) — `xtol` triggers when the
optimizer can no longer take a useful step, which also happens when a
parameter is pinned against a bound with the true unconstrained gradient
still pointing outward. This is not evidence of a true interior optimum for
those two; it means **the search ran out of room**, not that it converged.
Only `sigma0_way` landed at a genuine interior value.

Both `Fc_way` and `Fs_way` hit the *same* relative bound (exactly 5×)
together, which is the same kind of signature `calibration_bracketing/
SENSITIVITY_ANALYSIS.md` flags as `Fs_way`/`sigma0_way` collinearity —
worth an explicit correlation check before trusting the two independently
for this segment.

## What the fit actually improved

Checked directly (not just from the aggregate RMS number): the settled
plateau level and the peak transient overshoot after the first edge.

| | measured | baseline sim | optimized sim |
|---|---|---|---|
| settled plateau (forward, µm) | 9.326 | 10.0000 | 10.0000 |
| settled plateau (return, µm) | −0.025 | 0.0000 | 0.0000 |
| peak overshoot after edge (µm) | 11.74 | 16.72 | 13.73 |

**The fit is entirely about overshoot, not plateau level.** Both parameter
sets settle at *exactly* the commanded target regardless of friction
values — at MRES 1 every commanded position sits exactly on a detent
minimum (`sin(4·N_r·θ)=0` there by construction), so nothing in the
friction ports moves the steady state at all here. Raising `sigma0_way`/
`Fc_way`/`Fs_way` damps the post-edge ringing: baseline overshoots to
16.72 µm against a measured 11.74 µm peak, optimized brings that to
13.73 µm — real progress, but still ~2 µm over measured, consistent with
the fit not having reached a true optimum (still bound-limited on 2 of 3
parameters) rather than having found the best possible match.

## A limitation this run exposed, not fixed

The measured trace settles **short of the ideal target by ~0.67 µm**
(9.326 µm vs. the commanded 10.0 µm) on the forward plateau, while the
return plateau lands almost exactly on 0 (−0.025 µm). No combination of
`sigma0/Fc/Fs_way` can produce that: as shown above, this model's steady
state always lands exactly on the commanded position at MRES 1 positions,
independent of friction parameters. That forward-only ~0.67 µm shortfall
must come from something this run didn't touch — sensor/encoder scale
calibration, lead-screw calibration tolerance, or a backlash/asymmetry
this model doesn't represent — not from friction tuning. Worth resolving
*before* spending more optimizer time on this block, since it puts a
floor under achievable RMS no matter how well the transient is matched.

## Round 2 — widened bounds, warm-started from round 1

User observed the optimized trace still didn't quite reach the measured
overshoot level and asked about adding noise to close the gap. Checked
that numerically first (see chat/session log): adding noise at the real
sensor-floor level (0.022 µm std, measured from this block's own quiet
plateau) left the cost essentially unchanged (0.690→0.690 µm), and larger
noise only ever made it worse (1.0 µm std → 1.216 µm) — RMS error can't be
improved by adding independent noise to a biased signal, so the actual fix
was to let the bound-limited round-1 optimizer keep going, not add jitter.

Re-ran with bounds widened to 0.05×–20× the *original* Rev 4.2 baseline
(round 1 used 0.2×–5×) and warm-started from round 1's result. 17 accepted
`nfev` / 38 real model evaluations, 379 s.

| | baseline | round 1 | round 2 |
|---|---|---|---|
| RMS-to-measured | 0.9041 µm | 0.6897 µm (−23.7%) | **0.5762 µm (−36.3%)** |
| `sigma0_way` (N/m) | 3.45e6 | 7.319e6 | 6.864e6 |
| `Fc_way` (N) | 3.0 | 15.0 (bound) | **60.0 (new bound)** |
| `Fs_way` (N) | 4.5 | 22.5 (bound) | **90.0 (new bound)** |

**`Fc_way`/`Fs_way` pinned at the new bound again**, and moved together at
*exactly* the same ratio both times (22.5/15 = 90/60 = 1.5, matching the
baseline's own 4.5/3 = 1.5 ratio exactly). This is the important result of
round 2, more than the RMS number: within this 3-parameter subspace, cost
keeps falling as long as "way-port friction level" scales up, with no
interior optimum found even after a 4x bound expansion — and `sigma0_way`
moved *down* between rounds (7.32e6→6.86e6) as Fc/Fs got room to grow
further, meaning it isn't stable on its own either; its apparent optimum
depends on wherever Fc/Fs happens to be capped. That is a textbook
identifiability degeneracy signature, not three independently-calibrated
numbers.

Checked what actually improved, the same way as round 1 (peak transient
excursion after each edge, not just the aggregate RMS):

| | measured | baseline | round 1 | round 2 |
|---|---|---|---|---|
| forward peak (µm) | 11.3–11.7 | 16.72 | 13.73 | **12.38** |
| return trough (µm) | −1.9 to −2.2 | −6.72 to −6.74 | *(not checked)* | **−2.38 to −2.40** |
| forward plateau (µm) | 9.326 | 10.000 | 10.000 | 10.000 (unchanged) |

Genuinely good news: round 2 brought **both** the forward overshoot and
the return undershoot to within ~0.2–1.1 µm of measured (down from
baseline's 5+ µm misses on both), confirming increased way-port friction
is pushing the transient shape in the right direction on both polarities,
not just one. The ~0.65–0.67 µm plateau offset is exactly unchanged, as
predicted — still a friction-independent, unmodeled effect.

**Recommendation: stop widening bounds here.** `Fc_way = 60 N`,
`Fs_way = 90 N` are implausible for a linear-guide friction port on this
stage — compare the *nut* port's own baseline (`Fc_nut = 0.5 N`), and note
`model_parameters.json`'s own `sigma0_way_caveat` already flags
`sigma1_way` (776.7 N·s/m) as **not recomputed and inconsistent with the
target-zeta=0.7 rule** since `sigma0_way` was last revised — i.e. the
model's own documentation says the damping parameter, not the Coulomb/
static levels, is the one known to be stale. A third bounds-widening round
would likely just repeat this pattern at a still-larger, still-implausible
scale. The more physically motivated next experiment is to free `sigma1_way`
(and hold Fc_way/Fs_way at baseline or a modest value) — `sigma1_way`
directly controls bristle damping rate, which is mechanically the right
lever for "transient rings down faster than the model predicts," rather
than static/Coulomb friction magnitude.

## Honest bottom line

- The pipeline works end-to-end: real measured target, real 15-state
  nonlinear model, real gradient-based optimizer, real (if modest) cost
  reduction — this was an actual optimization run, not a demo.
- The specific numbers (`Fc_way`→15, `Fs_way`→22.5) are **not** a
  trustworthy calibration result — they're bound-limited, not converged,
  and possibly collinear with each other.
- Next steps before trusting any number here: widen (or remove) the 5×
  bound and see where `Fc_way`/`Fs_way` actually want to land; check their
  correlation directly; and separately investigate the ~0.67 µm plateau
  offset, since it's an unmodeled effect the optimizer cannot absorb no
  matter how it's run.
- Ramp legs were not attempted — extrapolating `one_eval_s` (~10 s for 6
  edges) to a 25 mm leg's edge count (hundreds to low thousands at fine
  MRES) puts a single evaluation at tens of minutes to hours; that needs a
  cheaper per-edge integration (e.g. adaptive segment merging, or a coarser
  model) before it's worth attempting.

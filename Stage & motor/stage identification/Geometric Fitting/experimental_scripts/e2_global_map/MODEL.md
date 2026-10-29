# The E2 geometric model

What `render_e2_findings.py` fits, why each term is there, and where every
number in it comes from. Companion to `README.md`, which covers the
measurement itself rather than the model.

Two figures draw the model rather than describe it, both written by
`render_e2_model_shape.py`:

| Figure | Shows |
|---|---|
| `rendered/model/e2_model_exploded_components.png` | one strip per term along the travel, all on the same µm scale, with a two-turn detail column — what the model is made of and how big each part is where |
| `rendered/model/e2_model_xtheta_table.png` | the periodic part over position × screw angle, with per-turn slices at 10, 50 and 78 mm — the compensation table itself |

---

## 1. What the model is for

E2 measures the direction-averaged positioning error over the full travel and
turns it into a continuous map

```
G(x, theta)   [um]
```

that predicts the error at any commanded position. Two things motivate its
shape:

- **The error is not a single smooth curve of position.** It carries a strong
  component that repeats once and twice per screw revolution — a 2 mm spatial
  period — riding on a slower drift along the travel.
- **That periodic component is not stationary.** Its amplitude and phase change
  with position: the once-per-revolution term runs from 1.07 µm to 2.31 µm
  across the travel. A model with fixed harmonic coefficients cannot represent
  that, which is why the harmonics carry position-dependent coefficients.

The model deliberately covers only the *direction-averaged* error. Reversal is
a separate quantity with a separate model (§11).

---

## 2. The response being fitted

The chain from capture to response, in `record_quality()` and `preprocess()`:

| Step | What happens |
|---|---|
| 1 | Each stop's position is the mean of the 0.5 s record window that opens 2.0 s after arrival (`settle_s`, `record_s`), ~454 samples at 1 kHz. |
| 2 | `raw = (measured_mm − commanded_mm) × 1000` → error in µm, signed positive when the stage overshoots the command. |
| 3 | Drift correction subtracted (§10), from the reference visits. |
| 4 | Up and down visits to the same target are paired by plan index; 400 pairs survive. |
| 5 | **`ybar = (e_up + e_down) / 2`**, then its mean is removed. This is the response `y`. |
| 6 | `B = e_up − e_down` is kept aside as the reversal. |

Averaging the two directions is what makes `y` *geometric*: reversal enters
`e_up` and `e_down` with opposite sign, so it cancels in the mean and survives
only in `B`. Everything the model fits is therefore error that depends on where
the stage is, not on how it got there.

The mean of `ybar` is removed because only differences along the travel are
identifiable — a constant offset is absorbed by the datum and carries no
information about geometry.

---

## 3. The two coordinates

### 3.1 Position `x`

Commanded position in mm. 400 targets spanning 2.28–119.92 mm, median spacing
0.28 mm.

The design (`e2_sequence.py`) divides 2–120 mm into 400 equal strata and takes
one target from each, so coverage is uniform by construction rather than by
luck. Targets within ±64 microsteps of an E0/ISO target are excluded, so the
E2 map and the ISO v3 session remain statistically independent and the ISO
cross-check in §12 is a genuine out-of-sample test.

### 3.2 Screw phase `theta`

```
usteps = round(x × 1600)
theta  = 2*pi * (usteps mod 3200) / 3200          [rad]
```

Where those two integers come from:

| Quantity | Value | Origin |
|---|---:|---|
| Screw lead | 2.0 mm | the ball screw: one revolution advances the nut 2 mm |
| Full steps per revolution | 200 | 1.8° stepper |
| Microstep resolution `MRES` | 16 | driver setting, `run_current_ma` 400, `intpol` off |
| **Microsteps per mm** | **1600** | 200 × 16 / 2 mm |
| **Microsteps per screw revolution** | **3200** | 1600 × 2 mm |

So `theta` is the angular position of the screw within one revolution, and a
harmonic at `k·theta` is an error that repeats `k` times per revolution — i.e.
with spatial period `2/k` mm. `k = 1` is 2 mm, `k = 2` is 1 mm.

### 3.3 The 64-microstep grid, and why it matters

Every target is placed on a multiple of 64 microsteps. 64 microsteps is
exactly one *electrical* cycle of the motor (4 full steps × 16 microsteps).

The motor has its own error that repeats every electrical cycle, and it is not
small. Putting every target on the same electrical phase makes that error a
**constant** across the whole dataset, where it is absorbed by the offset
instead of masquerading as position-dependent geometry. Without this the
screw and motor contributions would be confounded.

The design additionally bins `theta` into 50 bins and requires at least 5
targets in every bin (`phase_min_count`), searching seeds from 5200 upward
until that holds — seed **5315** is the accepted design. This guarantees the
harmonic coefficients are identifiable at every part of the travel rather than
only where the sampling happened to land.

---

## 4. The model equation

```
y(x, theta) =  a0
             + a1 * (x - xbar)                        gain
             + s(x)                                   slow trend
             + sum over k = 1,2   [ A_k(x) cos(k theta) + B_k(x) sin(k theta) ]
             + sum over k = 3..H  [ a_k    cos(k theta) + b_k    sin(k theta) ]
```

with `xbar` the mean of the fitted `x`. Built in `prepare_design()`; the
column layout for the selected model is

| Columns | Block | Count |
|---|---|---:|
| 0 | offset `a0` | 1 |
| 1 | gain `a1` | 1 |
| 2–5 | slow trend `s(x)` | 4 |
| 6–29 | `A_1(x)`, `B_1(x)` | 12 + 12 |
| 30–53 | `A_2(x)`, `B_2(x)` | 12 + 12 |
| 54–57 | stationary `k = 3, 4` | 4 |
| | **total** | **58** |

58 columns, 39.24 effective parameters after the penalty (§7).

---

## 5. What each term does

### 5.1 Offset `a0`

Carries nothing physical. `ybar` is mean-centred, so `a0` sits at ~0 and exists
only so the other terms are measured about the right level.

### 5.2 Gain `a1` — the linear scale error

A single slope over the whole travel, reported in ppm (`a1 × 1000`, since `a1`
is µm/mm). Currently **−29.81 ± 0.58 ppm**.

Physically this is everything that makes the stage's length scale differ from
the interferometer's, lumped together: screw lead error, thermal expansion of
the screw relative to the air path, the interferometer's wavelength
compensation, and cosine error from any residual misalignment. E2 cannot
separate those — none of temperature, pressure, compensation state or
alignment angle were recorded — so the gain is reported as one number and
attributed no further.

It is fitted rather than assumed because it is large: −29.81 ppm over 117.6 mm
of travel is 3.5 µm, comparable to everything else in the model combined.

### 5.3 Slow trend `s(x)` — the non-linear part of the drift

A cubic B-spline in `x` with `K_s` interior knots, **projected orthogonal to
the constant and linear columns** before use. The projection is the important
part:

```python
q          = [1, x - xbar]
projection = pinv(q) @ basis
projected  = basis - q @ projection
_, s, vt   = svd(projected)
transform  = vt[:rank].T             # rank = singular values > 1e-10 * s[0]
```

A raw B-spline basis can represent a constant and a straight line perfectly
well, so without this the slow trend and the gain would be fighting over the
same signal and neither would be interpretable. After projection, `s(x)` spans
only what a straight line *cannot* represent, so:

- the gain `a1` is the true whole-travel slope, not an artefact of knot
  placement, and
- `s(x)` is by construction the curvature that remains.

The SVD step then drops any direction left numerically empty by the
projection, which is why `K_s = 2` gives 6 raw basis columns but only **4**
columns in the design.

Interpretation: the bending of the error along the travel — the screw's own
pitch variation, support-bearing effects near the ends, anything that varies
over tens of mm rather than per revolution.

### 5.4 The varying harmonics, `k = 1` and `k = 2`

These are the heart of the model. Instead of two numbers per harmonic, each
coefficient is itself a cubic B-spline in `x`:

```
A_k(x) = sum_j alpha_kj  N_j(x)
B_k(x) = sum_j beta_kj   N_j(x)
```

so the harmonic's amplitude `sqrt(A^2 + B^2)` and phase `atan2(-B, A)` both
vary continuously along the travel. `harmonic_profile()` evaluates exactly
this, with uncertainties propagated through the delta method.

Fitted ranges:

| | amplitude min | amplitude max |
|---|---:|---:|
| `k = 1`, once per revolution (2 mm) | 1.07 µm | 2.31 µm |
| `k = 2`, twice per revolution (1 mm) | 0.16 µm | 0.70 µm |

A factor of 2.2 across the travel for `k = 1`. That variation is the thing a
stationary model gets wrong, and it is why `M_C` beats `M_B` (§9.4).

Physically: the screw's periodic error is not the same everywhere. Nut
preload, screw straightness and support conditions all change along the
travel, so the per-revolution signature changes with them.

**Only `k = 1` and `k = 2` are allowed to vary.** They dominate the periodic
content, and letting higher harmonics vary too would cost 24 columns each for
signal that the data cannot support.

### 5.5 The stationary harmonics, `k = 3` to `H`

Two columns each, constant along the travel. They catch the remaining
periodic content — a non-sinusoidal per-revolution shape needs more than the
fundamental to describe it — without claiming any position dependence.

With `H = 4`, that is `k = 3` (0.667 mm period) and `k = 4` (0.5 mm).

---

## 6. The B-spline bases

`raw_bspline(x, interior, domain)`:

- **degree 3** (cubic). Twice differentiable, which is what makes a
  second-difference penalty meaningful.
- **interior knots equally spaced** over the domain. No data-adaptive knot
  placement — with a penalty doing the regularising, knot placement stops
  mattering, which is exactly why the penalised form is preferred (§7).
- **clamped** at both ends (the end knots are repeated `degree + 1` times), so
  the basis spans the full domain.
- **`domain` = (min x, max x) of the data being fitted**, and `x` is clipped to
  it before evaluation. Evaluating outside the fitted range therefore holds the
  edge value rather than extrapolating a cubic — the one thing that keeps
  predictions finite when the ISO cross-check asks for targets near 2 mm and
  117.5 mm.
- column count = `interior + degree + 1`. For `K_C = 8`: 12 columns.

---

## 7. Regularisation — the penalty

`coefficient_penalty()` builds a **second-difference penalty** on each of the
four coefficient-spline blocks (`A_1`, `B_1`, `A_2`, `B_2`): 4 blocks × 10 rows
= 40 penalty rows over the 58 columns.

The fit minimises

```
|| y - X beta ||^2  +  lambda * || D beta ||^2
```

where `D` takes second differences between *neighbouring B-spline
coefficients*. Penalising second differences penalises curvature in `A_k(x)`
and `B_k(x)`: `lambda -> 0` lets them wiggle freely, `lambda -> infinity`
drives them to straight lines, and a straight line of zero slope is the
stationary model. So `lambda` interpolates continuously between `M_C` and
`M_B`.

**Only the C1/C2 blocks are penalised.** The offset, gain, slow trend and
stationary harmonics are untouched — the penalty's job is to control how fast
the screw harmonics are allowed to change along the travel and nothing else.
This has a consequence that bit during development and is now guarded against
in §9.3.

**Effective degrees of freedom** replace the raw column count:

```
edf = trace( X (X'X + lambda P'P)^-1 X' )
```

58 columns, `edf = 39.24`. The error variance uses `n - edf`, and BIC uses
`edf`, so a penalised fit is charged for the flexibility it actually uses.

### Why penalised rather than choosing a knot count

Measured directly (see `TEMP_E0_E2_update_notes.md` §8): fixed knots force a
trade between interior accuracy and edge stability. `K_C = 3` already costs 3×
`K_C = 1` on a held-out end block; `K_C = 4` and `K_C = 5` are not merely worse
but unusable there (401 µm and 22 532 µm). The penalty removes the trade — a
rich `K_C = 8` basis with `lambda` chosen by CV gets the interior accuracy of
`K_C = 5` *and* better edge behaviour than `K_C = 1`.

---

## 8. Estimation

`fit_linear()`:

| | |
|---|---|
| Solver | least squares on the augmented system `[X; sqrt(lambda) P]` vs `[y; 0]` |
| Covariance | `sigma^2 (X'X + lambda P'P)^-1` |
| `sigma^2` | `RSS / (n - edf)` |
| Conditioning | `cond` of the column-normalised *augmented* matrix — the system actually solved. Currently **2.06** (it was 5.89 before the penalty). |
| Robustness | a Huber refit is run as a **diagnostic only**. The reported coefficients are the least-squares ones; the Huber fit reports how far any coefficient moves, in units of its own standard error. Currently max shift **0.23 sigma**, 0 coefficients over 1 sigma — no target is driving the fit. |

---

## 9. Model selection

### 9.1 The three families

| Family | What it is | Role |
|---|---|---|
| `M_B` | stationary harmonics, no varying coefficients | the null the varying model must beat |
| `M_C_fixed` | varying coefficients on a fixed knot grid, `K_C` = 1…5 | diagnostic only; it is what panel (a) of appendix A2 plots |
| `M_C` | varying coefficients, rich `K_C = 8` basis, penalised | **the selectable family** |

Grid: `H` ∈ {2, 3, 4}, `K_s` ∈ 0…8, and for `M_C` `lambda` over nine values
from 10⁻² to 10².

### 9.2 Repeated blocked cross-validation

`cv_assignments()` splits the travel into 10 contiguous blocks and repeats that
over **five shifted sets of block edges** (offsets 0, 0.2, 0.4, 0.6, 0.8 of a
block width). Blocks, not random folds, because neighbouring targets are 0.28
mm apart and a random split would leak: a held-out point would sit between two
training points a third of a millimetre away and be trivially interpolable.

Only interior folds are scored. The two end folds ask the model to predict
*beyond* its training support, which is extrapolation, not the quantity being
selected on — but they are recorded, and §9.3 uses them.

The repetition is not cosmetic. One partition yields eight interior scores and
a standard error of ±0.044 µm — wide enough that the one-SE rule admitted
`K_C = 1`, the worst varying-coefficient model on every other measure. Five
partitions give ±0.012 µm and the rule discriminates.

### 9.3 The one-SE rule, and the stability guard

```
best     = lowest CV
eligible = CV <= best_CV + best_SE
stable   = max(end_low, end_high) <= 5 * own CV          <- the guard
choose   = fewest effective parameters among stable
```

Cost is *effective* parameters, so penalised and fixed-knot models are compared
on the flexibility they use rather than on basis size.

**The guard exists because of a real failure.** Without it the rule selected
`H = 4, K_s = 7, lambda = 3.2` — respectable interior CV, 39.2 edf, and an
end-fold error of **440 µm**. The penalty covers the C1/C2 curves but *not*
the slow trend, so a large `K_s` buys flexibility the interior folds never
charge for. Passing models sit near 3× their interior score; failures run 60×
to 20 000×, so `CV_END_FACTOR = 5` is nowhere near a borderline cut.

### 9.4 `M_C` over `M_B`

Two conditions, both required:

1. `M_C`'s CV beats `M_B`'s by more than the larger of the two standard errors
   — currently 0.148 µm against a 0.015 µm requirement, so by 10×;
2. `M_C` also wins on the independent E1/E4 windows (`V2_M_C_better`).

### 9.5 Selected

| Parameter | Value |
|---|---|
| `H` maximum harmonic | 4 |
| `K_s` slow-trend knots | 2 |
| `K_C` coefficient knots | 8 (penalised) |
| `lambda` | 1.0 |
| columns | 58 |
| effective parameters | 39.24 |
| in-sample residual rms | 0.369 µm |

---

## 10. Drift variants D0 / D1 / D2

The run takes ~39 minutes and the stage drifts over that time. Seven reference
visits to a fixed position (`x_ref = 61.35 mm`) measure it: 3 DATUM at the
start, 1 TURN_REF at the turnaround, 3 END_REF at the end.

| Variant | Correction |
|---|---|
| `D0` | none |
| `D1` | linear in time between the DATUM and END_REF anchors |
| `D2` | piecewise linear through all three anchor groups, including the turnaround |

**`D2` is the default** and what the reported model is fitted on. `D0` and
`D1` are refitted in `appendix_drift()` as a sensitivity check — if a
conclusion moved between them it would not be a conclusion.

Drift is corrected in *time*, before the direction averaging, because it is a
property of when a stop happened, not where.

---

## 11. What the model deliberately excludes

| Excluded | Why, and where it lives instead |
|---|---|
| **Reversal** | Direction-dependent, not geometric. Fitted separately by `reversal_fit()` on `B = e_up − e_down`, with its own offset, slope and `k = 1, 2` harmonics. Currently −1.757 µm (D2). |
| **Electrical-cycle error** | Designed out by the 64-microstep grid (§3.3); constant across the dataset, absorbed by the offset. |
| **Thermal/time dependence** | Removed by the drift correction, not modelled as a term. |
| **Higher varying harmonics** | `k ≥ 3` are stationary; the data cannot support 24 more columns each. |

---

## 12. How the model is checked

| Check | What it tests | Current |
|---|---|---|
| **V2** — E1/E4 windows | Ten windows from two *other* experiments, never used in fitting. Residual rms, and the agreement of the model's `C1`/`C2` against each window's own measured spectrum. | window rms 0.378 µm; `\|C1\|+\|C2\|` error 0.245 µm; 92.2 % of the gap to each window's local floor closed |
| **V3** — ISO v3 cross-check | Applies `G` (and then the reversal) to the ISO session's per-target means, which are independent by design (§3.1). | E 7.855 → 1.482 µm; M 6.463 → 1.218 µm |
| Blocked CV | interpolation error | 0.434 ± 0.012 µm |
| End folds | extrapolation stability | 1.21 / 0.63 µm |

The `mean_external_local_floor_um` of 0.306 µm is the useful reference point:
it is the rms each window achieves fitting *itself*, so it is the best any
global map could do there. The model reaches 0.378 µm against that floor.

---

## 13. Every numeric constant and where it comes from

### Physical / design

| Constant | Value | Origin |
|---|---:|---|
| `lead_mm` | 2.0 | ball screw |
| `steps_per_rev` | 200 | 1.8° stepper |
| `mres` | 16 | driver setting |
| `ustep_per_mm` | 1600 | derived: 200 × 16 / 2 |
| `screw_rev_usteps` | 3200 | derived: 1600 × 2 |
| `motor_cycle_usteps` | 64 | 4 full steps × 16 — one electrical cycle |
| `target_count` | 400 | one per stratum over 2–120 mm |
| `phase_min_count` | 5 | minimum targets per screw-phase bin, 50 bins |
| seed | 5315 | first seed ≥ 5200 satisfying the phase-bin constraint |
| `x_ref_mm` | 61.35 | the reference position, mid-travel |
| `velocity_mm_s`, `acceleration_mm_s2` | 2.0, 40.0 | motion profile |
| `settle_s`, `record_s` | 2.0, 0.5 | dwell structure |

### Model

| Constant | Value | Origin |
|---|---:|---|
| spline degree | 3 | cubic; needed for a second-difference penalty to mean curvature |
| `PEN_KNOTS` | 8 | deliberately rich — `lambda` regularises, so the count stops mattering |
| `PEN_ORDER` | 2 | second differences → penalise curvature of `A_k(x)`, `B_k(x)` |
| `PEN_LAMBDAS` | 10⁻²…10², 9 values | selected `lambda = 1.0` is interior, so the grid is not binding |
| `H`, `K_s`, `K_C`, `lambda` | 4, 2, 8, 1.0 | chosen by §9 |

### Selection

| Constant | Value | Origin |
|---|---:|---|
| CV blocks | 10 | contiguous, to defeat 0.28 mm spacing leakage |
| `CV_OFFSETS` | 0, .2, .4, .6, .8 | five partitions; SE ±0.044 → ±0.012 µm |
| `CV_MIN_TEST` | 5 | a shifted end block can be narrow |
| `CV_END_FACTOR` | 5.0 | passes near 3×, failures ≥60× — a wide margin |

### Quality thresholds

| Constant | Value | Origin |
|---|---:|---|
| `STD_LIMIT_UM` | 0.5 | per-record standard deviation; a stop noisier than this is not settled |
| `SLOPE_EFFECT_LIMIT_UM` | 0.25 | drift *within* a record window, over its 0.5 s |
| `SATURATION_INCREMENT_NM` | 8333 | the AquadB link ceiling: 1 nm per increment at a 120 ns clock |
| `reference_drift_warning_um` | 0.5 | DATUM-to-END_REF drift over the run |
| `MOTOR_END_MM`, `SUPPORT_BEARING_MM` | 0.0, 122.770 | mechanical ends, shaded on every figure |

---

## 14. Diagnostics worth watching

| Diagnostic | Meaning | Current | Concern if |
|---|---|---:|---|
| `condition_number_scaled` | conditioning of the solved system | 2.06 | large — columns fighting |
| `max_abs_slow_vs_C_column_correlation` | overlap between `s(x)` and the harmonic columns | 0.138 | approaching 1 — the slow trend and the harmonics become unidentifiable |
| `huber_max_shift_sigma` | largest coefficient move under robust loss | 0.231 | > 1 — a few targets driving the fit |
| `huber_coefficients_over_1sigma` | how many do | 0 | > 0 |
| `median_record_mean_noise_um` | noise on each stop's mean | 0.0003 µm | — it is ~1000× below the residual, so the model is limited by structure, not measurement |

That last row is the one that frames everything: with a 0.369 µm residual
against a 0.0003 µm measurement noise, every µm in this model is real stage
behaviour.

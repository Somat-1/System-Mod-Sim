# Step A — Linear Geometric Error from the ISO 230-2 Baseline

**Goal:** extract the linear (gain) component of the geometric error G(x) from the existing ISO 230-2 baseline data. Quantify how much of the error it explains, and prepare the residual for the next steps (periodogram → B → C → D).

**Input:** the ISO baseline session (27 targets, 5 bidirectional runs, 5–135 mm).
**Model:** ȳ(x) = a₀ + a₁·x
**Scope:** position-only. Friction is *not* modelled here. The reversal is only characterized.

---

## 0. Outputs this step must produce

| Output | Used for |
|---|---|
| a₁ ± CI (µm/mm and ppm) | gain correction; comparison with screw class; thermal attribution |
| a₁ per run and per direction | stability check; flag for friction asymmetry varying with position |
| residual r(x) = ȳ − (a₀ + a₁x) | periodogram, starting point for model B |
| share of error power explained by the gain | thesis error budget |
| predicted ISO metrics after gain correction | first bar of the improvement chart |
| reversal Bᵢ vs x (characterized, not modelled) | input for the friction work |

---

## 0b. Adaptations for the ISO230-2_v2 data set (2 mm/s run)

The steps below were written for the planned 27-target sequence. The executed v2 run differs, so apply these adaptations:

| Item | Planned | v2 data | Action |
|---|---|---|---|
| Targets | 27, 5–135 mm | **m = 24**, ≈ 2–120 mm | use the actual target list; x_c = mean of the 24 targets; span = P_max − P_min |
| Datum | x_ref datum pair after warm-up | **global mean d** over all stops (your formula (1)) | fine for Step A: d only shifts a₀. State the convention; keep it identical for E1–E5 |
| Reference value | rounded commanded position | **nominal Pᵢ** | at 1/16 µstep (0.625 µm/µstep) rounding is up to ±0.31 µm per target, a fixed per-target offset that G would fit as "geometry". **Check that every Pᵢ is an integer multiple of 0.625 µm**; if not, use the rounded commanded position |
| Recorded windows | 24 × 2 × 5 = **240** | **221** windows in the AquadB check | find which stops are missing; use the actual nᵢ↑, nᵢ↓ per target in §4 (variance formula already uses n per target). Exclude targets with n < 3 in one direction, or flag them |
| Stop timing | common clock | **clock map** (r = 0.999960, c = 43 725 ms) | add per-window checks: window std and slope (§1), and confirm the window contains no saturated AquadB sample. Already shown: 0/221 saturated, minimum margin 1.14 s ✔ |
| Drift / count-loss reference | x_ref visits | the **61.35 mm reference visits** between runs | use them for the drift correction in §3 **and** as a count-loss check: AquadB saturation during moves showed up to ~0.8 µm travel mismatch in the v1 test, the same order as B |
| Hold current | settle after TPOWERDOWN | 400 mA run; hold setting to confirm | if IHOLD < IRUN, check the 0.5 s record window starts after the TPOWERDOWN + IHOLDDELAY transition; if IHOLD = IRUN there is no issue |
| Temperature | air + stage logged | environmental temperature via the IDS | thermal attribution in §8 limited to air temperature; state that the stage temperature wasn't logged |
| Validation role | — | ISO 230-2 §5.2: acceptance targets **shall differ** from compensation sampling points | supports using this data set for **Step A + validation only**, and fitting B/C on E2/E3 |

Reported v2 metrics (for reference): B = 0.98 µm, R = 2.43 µm, E = 6.61 µm, M = 6.36 µm, A = 7.95 µm.
M ≫ B, so **the direction-averaged (geometric) error dominates, and the reversal is about 1/6 of it**. This is the expected outcome for Step A, and it backs the argument that the friction model alone can't compensate this stage.

---

## 1. Load and check validity

Files from the session folder: `stops.csv`, `events.csv`, `env.csv`, `summary.json`.

Checks, before any fitting:
1. **Settling:** per stop, |slope| × t_record < ½ of the expected repeatability. ≥ 95 % of stops must pass. List the failures.
2. **Watchdog:** no aborts or lost-step flags during ISO_CYCLE.
3. **Drift:** plot x_ref against time over the session. Note the total drift and the drift rate.
4. **Environment:** ΔT_air and ΔT_stage over the session; mean stage temperature T̄_stage (needed in §8).

Record everything in a short validity table. It goes in the thesis appendix.

---

## 2. Deviation per stop

For each recorded ISO stop (run j, direction ↑/↓, target i):

$$x_{ij}^{\uparrow/\downarrow} = (x_{\text{IFM}} - x_{\text{IFM,datum}}) - (x_{\text{cmd}} - x_{\text{cmd,datum}})$$

- Use the **commanded position after rounding to microsteps**, not the nominal target.
- Units: **µm** for deviations, **mm** for positions.
- Sign convention: positive means the stage went *further* than commanded. Keep it fixed and state it.

---

## 3. Drift correction (do both, report both)

- **Uncorrected:** use x_ij as is.
- **Corrected:** interpolate the x_ref readings linearly in time and subtract the interpolated drift at each stop's timestamp.

Carry both versions through §4–§7. If a₁ differs between them by more than its CI, drift is significant and must be discussed. Otherwise, report the corrected version as primary.

---

## 4. Per-target statistics

For each target i, using n = 5 runs per direction:

| Quantity | Formula |
|---|---|
| mean per direction | x̄ᵢ↑, x̄ᵢ↓ |
| std per direction | sᵢ↑, sᵢ↓ (n − 1 normalization) |
| **direction-averaged mean** | **ȳᵢ = ½(x̄ᵢ↑ + x̄ᵢ↓)** |
| its variance | var(ȳᵢ) = ¼(sᵢ↑²/n + sᵢ↓²/n) |
| reversal | Bᵢ = x̄ᵢ↑ − x̄ᵢ↓ |

ȳᵢ is the input for G. Bᵢ is only characterized (§9).

---

## 5. Weighted linear fit

**Centre the position** to decorrelate offset and slope: x̃ = x − x_c, with x_c = mean of the target positions (~70 mm).

$$\bar y_i = a_0' + a_1\,\tilde x_i \qquad (a_0 = a_0' - a_1 x_c)$$

- **Weights:** wᵢ = 1 / var(ȳᵢ).
- **Fit:** weighted least squares (`numpy.linalg.lstsq` on √w-scaled rows, or `statsmodels.WLS`).
- **Important: scale the covariance by the reduced χ².** The linear model ignores the periodic and slow error, which act as extra "noise" around the line. var(ȳᵢ) only contains repeatability, so the raw WLS confidence interval would be **far too small**:

$$\chi^2_\nu = \frac{\sum w_i r_i^2}{N-2}, \qquad \text{Cov}_{\text{reported}} = \chi^2_\nu \cdot (\Phi^T W \Phi)^{-1}$$

  Report χ²_ν itself. A value ≫ 1 is expected, and is direct evidence that the linear model is incomplete (motivation for B).
- **Cross-check:** ordinary least squares (unweighted). a₁ should agree within the CI.

**Report:**
- a₁ in µm/mm and **ppm** (1 µm/mm = 1000 ppm), with the 95 % CI;
- the total gain error over the span: a₁ × (P_max − P_min) [µm] (≈ 118 mm for v2);
- χ²_ν.

---

## 6. Robustness of a₁

| Check | How | What it tells you |
|---|---|---|
| **Per run** | fit a₁ on each run's direction-averaged data separately (5 values) | run-to-run stability; drift or thermal trend within the session if a₁ changes monotonically with run index |
| **Per direction** | fit a₁↑ on x̄ᵢ↑ and a₁↓ on x̄ᵢ↓ separately | **a₁↑ ≈ a₁↓ expected.** A significant difference means the reversal changes linearly along the travel → friction (or preload) varying with position. Flag it for the friction work; don't model it here |
| **Leave-one-target-out** | refit 27× without one target | influence of single points (e.g. end targets); any target that moves a₁ by > ½ CI is suspicious |

---

## 7. Residual and error-power share

$$r_i = \bar y_i - (a_0 + a_1 x_i)$$

Compute:
- rms(ȳ) before correction, rms(r) after correction;
- **share of the direction-averaged error power explained by the gain:** 1 − Σrᵢ² / Σ(ȳᵢ − mean ȳ)²;
- the same share relative to the **total** error budget (including reversal and repeatability), so it can be compared with the earlier ~42 %.

Plot ȳᵢ with the fitted line, and rᵢ against x. The residual should show visible structure (periodic and slow). That's the input for the next steps.

---

## 8. Physical interpretation of a₁

Split a₁ into its possible causes. Report each as an estimate, not as a fit:

| Cause | Estimate | Note |
|---|---|---|
| **Thermal expansion of the screw** | α·(T̄_stage − 20 °C), α ≈ 11.5 ppm/K for steel | ISO 230-2 refers to 20 °C |
| **Interferometer wavelength** | residual after environmental compensation (≈ 1 ppm/K air, ≈ 0.3 ppm/hPa) | check that compensation was active; otherwise correct with Edlén |
| **Cosine error** | ≈ −θ²/2 for beam misalignment θ | typically < 1 ppm with good alignment |
| **Lead deviation of the screw** | the rest: a₁ − thermal − wavelength − cosine | compare a₁ × 300 mm with the ISO 3408 e_p / v_300p tolerance of the screw's accuracy class |

**Conclusion to state:** whether the gain is within the screw's specification, and how much is thermal (which determines whether G needs a temperature term later).

---

## 9. Characterize the reversal (no modelling)

- Plot Bᵢ against x with its uncertainty: σ(Bᵢ) = √(sᵢ↑²/n + sᵢ↓²/n).
- Report its mean, spread and trend: a linear fit of Bᵢ against x, slope ± CI.
- Note any visible periodic pattern in Bᵢ. That would be friction varying per revolution (see the periodogram in §11).

These numbers are handed over to the friction work unchanged.

---

## 10. Predicted ISO metrics after gain correction

Recompute the ISO 230-2 quantities with the linear correction applied to the data:

$$\bar x_{i\uparrow/\downarrow}^{\,\text{corr}} = \bar x_{i\uparrow/\downarrow} - (a_0 + a_1 x_i)$$

| Metric | Expected change |
|---|---|
| E↑, E↓, E, M (systematic deviations) | **decrease** |
| A (accuracy) | decreases |
| B (reversal), R (repeatability), sᵢ | **unchanged by construction.** If they change, there's a bug |

This "predicted after gain correction" row is the first step in the improvement chart. It's confirmed experimentally later by the G-compensated ISO run (experiment E5).

---

## 11. First look at the periodic content (preliminary)

- Lomb–Scargle periodogram of rᵢ: `astropy.timeseries.LombScargle(x, r, dy)`, spatial frequency 1/118 to ~1/2 mm⁻¹ (above that, 24 points carry no information).
- With 24 points, expect at best a hint of a peak at 1/L (= 0.5 mm⁻¹). Treat it as preliminary; the clean periodogram comes from experiment E1.
- **Fold** rᵢ onto the lead phase: φᵢ = (xᵢ mod L)/L, and plot rᵢ against φᵢ. A sinusoid-like pattern shows the periodic error directly, even with sparse data.
- Do the same for Bᵢ (friction-related periodicity, see §9).

---

## 12. Script outline

```python
# ── STEP A: linear geometric error from the ISO baseline ───────────────
import numpy as np, pandas as pd
from astropy.timeseries import LombScargle

L = 2.0                                         # screw lead [mm]

# 1–3: load, validate, deviations (µm), optional drift correction
st = pd.read_csv("stops.csv")
iso = st[st.state == "ISO_CYCLE"].copy()
# v2 convention: deviation vs reference position, minus the global datum d (your formula (1))
P = iso.cmd_pos_mm                                # use ROUNDED commanded position (k · 0.625 µm)
assert np.allclose((P * 1600) % 1, 0, atol=1e-6), "targets not on the 1/16 µstep grid"
iso["dev"] = (iso.x_mean - P) * 1e3               # µm
iso["dev"] -= iso["dev"].mean()                   # global datum d
# iso["dev"] -= drift_interp(iso.t_record_start)  # corrected variant, from the 61.35 mm visits

# 4: per-target statistics
g = iso.groupby(["target_i", "direction"]).dev.agg(["mean", "std", "count"]).unstack()
x = iso.groupby("target_i").cmd_pos_mm.first().values
m_up, m_dn = g[("mean", "up")].values, g[("mean", "down")].values
s_up, s_dn = g[("std", "up")].values,  g[("std", "down")].values
n = g[("count", "up")].values
y   = 0.5 * (m_up + m_dn)
var = 0.25 * (s_up**2 / n + s_dn**2 / n)
B   = m_up - m_dn

# 5: weighted, centred linear fit with χ²-scaled covariance
xc = x.mean()
Phi = np.column_stack([np.ones_like(x), x - xc])
w = 1 / var
A = Phi.T @ (w[:, None] * Phi)
beta = np.linalg.solve(A, Phi.T @ (w * y))
r = y - Phi @ beta
chi2nu = np.sum(w * r**2) / (len(x) - 2)
cov = chi2nu * np.linalg.inv(A)
a1, a1_se = beta[1], np.sqrt(cov[1, 1])          # µm/mm
print(f"a1 = {a1*1e3:.1f} ± {1.96*a1_se*1e3:.1f} ppm   (chi2_nu = {chi2nu:.1f})")

# 7: explained share of direction-averaged error power
share = 1 - np.sum(r**2) / np.sum((y - y.mean())**2)

# 11: preliminary periodogram + folding
f = np.linspace(1/140, 0.6, 2000)
p = LombScargle(x, r, np.sqrt(var)).power(f)
phase = (x % L) / L
```

---

## 13. Deliverables checklist

- [ ] Validity table (settling, drift, watchdog, environment)
- [ ] a₁ ± CI in ppm, χ²_ν, OLS cross-check, drift-corrected vs uncorrected
- [ ] a₁ per run, per direction, leave-one-out
- [ ] Explained share of error power (direction-averaged and total)
- [ ] Physical attribution table (thermal, wavelength, cosine, lead vs screw class)
- [ ] Bᵢ vs x with trend
- [ ] Predicted ISO metrics after gain correction (B and R unchanged)
- [ ] Preliminary periodogram + folded residual (for rᵢ and Bᵢ)
- [ ] Plots: ȳ with fit; residual vs x; folded residual; Bᵢ vs x; x_ref vs time

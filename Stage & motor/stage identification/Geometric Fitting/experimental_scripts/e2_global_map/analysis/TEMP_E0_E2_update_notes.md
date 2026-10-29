# Temp note — the corrected E0, what it does to E2, and what v4 would do

Working note, not a deliverable. Delete when the content has been folded into
the E2 README / TABLES.

Covers three things:

1. why E0 changed and by how much,
2. what the E2 re-run produced,
3. what would happen if E0 and E2 were built on ISO v4 instead of v3.

---

## 1. Why E0 changed

The v3 ISO evaluation took each 0.5 s averaging window from the board log and
mapped it onto the interferometer trace with a fitted clock. The clock was
fitted against a stop detector that smoothed velocity with a **centred
51-sample boxcar**, which reports rest about 34 ms late. The fit absorbed that
as a ~38 ms offset on the whole board timeline, and because the firmware starts
the next move at the instant the record window closes, there was no margin: all
240 windows closed about 38 ms **into the following move**.

In those 38 ms the stage travels ~23 µm, always continuing in the current pass
direction, so each window mean is dragged by:

| | shift |
|---|---:|
| `xbar_up` | −0.617 ± 0.021 µm |
| `xbar_dn` | +0.617 ± 0.028 µm |
| direction average `ybar` | ≤ 0.050 µm, rms 0.017 µm, no position dependence (+0.065 ± 0.216 ppm) |

`B = xbar_up − xbar_dn` subtracts two biases of opposite sign, so they add:
B shifts by +1.235 µm. `(xbar_up + xbar_dn)/2` sums them, so they cancel — which
is why the gain, the residual and M looked fine all along.

Full diagnosis: `iso230-2/analysis/ISO230_2_V3_REVERSAL_WINDOWING_POSTMORTEM.md`.
Evidence figure: `iso230-2/analysis/plots/temp/v3_record_windows.png`.

The correction is not a reinterpretation — keeping the published pipeline and
merely dropping the **last 40 ms** of each logged window reproduces the blind
result to three decimals on B and R.

---

## 2. E0 before and after

E0 is now built from `iso230_2_v3_results_blind.json` (holds found in the trace,
averaging window placed inside the measured hold, board clock not used), run
through the unchanged `process_iso_v3_step_a.py`. Outputs live in
`Geometric Fitting/results_v3_corrected/`.

| Quantity | published | corrected | note |
|---|---:|---:|---|
| G slope | −25.09 ± 13.72 ppm | **−27.43 ± 14.22 ppm** | both inside each other's CI |
| reduced χ² | 153.8 | 168.6 | |
| residual rms | 1.108 µm | 1.171 µm | |
| residual peak-to-peak | 4.042 µm | 4.137 µm | |
| E↑ / E↓ | 6.116 / 6.611 | 6.192 / 6.735 | per-direction ranges, barely move |
| M | 6.364 µm | 6.463 µm | |
| **E** | 6.611 µm | **7.855 µm** | |
| **A** | 7.949 µm | **9.151 µm** | |
| **R** | 2.427 µm | **3.617 µm** | |
| **B̄** | −0.365 µm | **−1.601 µm** | |
| **B** (max\|Bᵢ\|) | 0.982 µm | **2.264 µm** | |
| B trend | +1.66 ± 4.32 ppm | +2.05 ± 4.19 ppm | still consistent with zero |

**The slope move is a weighting effect, not a signal one.** Feeding the
corrected ȳ through the *published* weights gives −25.01 ppm; feeding the
published ȳ through the *corrected* weights gives −27.65 ppm. Clipping the
window changes each target's s↑/s↓ slightly, and with χ²_ν ≈ 150 the weighted
fit is sensitive to that. No conclusion depends on it.

F1 and F2 have been re-rendered. `results/` and `E0/TABLES.md` still carry the
old numbers, so **F1 currently disagrees with T1 and T2**.

---

## 3. The E2 re-run

`render_e2_findings.py` was repointed at the corrected inputs (two path
constants, now env-overridable via `E2_STEP_A` / `E2_ISO_RESULTS`). No method,
model or plotting code was touched.

Of 291 leaf values in `e2_findings_summary.json`, **262 are identical**. That is
the expected result: E2's own capture was never affected, so its model, drift
variants, cross-validation, reversal fit and validity table are unchanged.

### Unchanged

| | value |
|---|---|
| selected model | M_C, H = 4, K_s = 2, K_C = 1 |
| CV RMSE | 0.487 µm (vs M_B 0.597 µm) |
| residual rms | 0.397 µm |
| E2 gain | −29.89 ± 0.61 ppm |
| reversal D0 / D1 / D2 | −1.614 / −1.758 / −1.757 µm |
| V2 window validation | raw 1.220 → M_B 0.525 → M_C 0.399 µm |
| validity table | identical, all checks unchanged |

### Changed

Only the Step A reference and the v3 cross-check block:

| | before | after |
|---|---:|---:|
| `gain.Step_A_ppm` | −25.09 | **−27.43** |
| `gain.Step_A_ci95_ppm` | 13.72 | 14.22 |
| raw E / M / A / R | 6.611 / 6.364 / 7.949 / 2.427 | **7.855 / 6.463 / 9.151 / 3.617** |
| after G: E / M / A | 1.963 / 1.129 / 3.311 | **3.317 / 1.335 / 4.971** |
| after G+B: E / A | 2.304 / 4.030 | **1.560 / 3.353** |

Figures regenerated: all eleven. Only **figure 3(b)** and
**appendix e0-electrical-phase** change in content; the rest differ only by a
1–2 px glyph-positioning shift from a different matplotlib build (3.10.9 here),
same typeface, same layout.

---

## 4. New E2 findings

### 4.1 E2's reversal model now transfers to the v3 session

This is the substantive change. E2 predicts a reversal of β₀ = −1.757 µm (D2).
Removing it from the v3 per-target means leaves:

| v3 source | mean up−down | max | rms |
|---|---:|---:|---:|
| as published | +1.392 µm | 2.118 | 1.438 |
| **corrected** | **+0.157 µm** | **0.900** | **0.399** |

Against the published v3, E2's reversal appeared to over-correct by 1.39 µm —
which would have read as "the reversal does not transfer between sessions".
Against the corrected v3 the residual is 0.157 µm, about nine times smaller and
comparable to E2's own 0.397 µm model residual. **The apparent session-to-session
reversal mismatch was the windowing defect, not physics.**

### 4.2 The predicted ISO improvement is larger than reported

| | before | after |
|---|---|---|
| E, raw → after G+B | 6.611 → 2.304 µm (−65 %) | **7.855 → 1.560 µm (−80 %)** |
| M, raw → after G | 6.364 → 1.129 µm | **6.463 → 1.335 µm** |
| mean-deviation rms | — | 1.414 → **0.344 µm** |

### 4.3 E2 and E0 now agree better on the gain

| | E2 | Step A | gap |
|---|---:|---:|---:|
| before | −29.89 ± 0.61 ppm | −25.09 ± 13.72 ppm | 4.80 ppm |
| after | −29.89 ± 0.61 ppm | **−27.43 ± 14.22 ppm** | **2.46 ppm** |

Both were always inside Step A's CI; the agreement is simply tighter now.

### 4.4 Not a bug (correcting an earlier suspicion)

`V3_ISO_metrics_um.G_plus_reversal.R` equals the raw R rather than being
recomputed. That is deliberate — `render_e2_findings.py:554-559` overwrites
`R_up`/`R_down`/`R` at the G and G+B levels with the measured values, with the
comment that compensation changes the predicted systematic terms E/M/A while
repeatability and the measured scatters do not. I had flagged this as a possible
R/A transposition; it is not.

---

## 5. What if E0 and E2 used ISO v4 instead?

### 5.1 The two sessions measure the same stage

| | value |
|---|---|
| corr(ȳ_v3, ȳ_v4) | **+0.965** |
| corr(residual_v3, residual_v4) | +0.790 |
| per-target Bᵢ correlation | +0.871 |
| ȳ difference | 0.399 µm rms, mean 0.000 |
| ȳ difference vs position | +5.03 ± 4.71 ppm |

### 5.2 E0 on v4 — the primary estimator breaks

Same `process_iso_v3_step_a.py`, `--aggregate iso230_2_v4_results_blind.json`:

| Estimator | v3 published | v3 corrected | **v4** |
|---|---:|---:|---:|
| weighted (primary) | −25.09 ± 13.72 | −27.43 ± 14.22 | **+3.55 ± 17.50** |
| OLS cross-check | −24.17 | −24.04 | −19.01 |
| Theil–Sen | — | −20.29 | −15.60 |
| leave-one-target-out | −28.07 … −20.67 | −30.27 … −21.85 | **−18.15 … +11.10** |
| reduced χ² | 153.8 | 168.6 | **714.2** |
| residual rms | 1.108 | 1.171 | 1.564 µm |
| residual peak-to-peak | 4.042 | 4.137 | **7.260 µm** |
| direction-average power explained | 34.3 % | 31.5 % | **−7.5 %** |

The weighted fit on v4 returns a *positive* slope and a **negative** explained
power share — i.e. applying it makes the direction average worse. With
χ²_ν = 714 the weights are near-degenerate and a couple of low-s targets
dominate. The robust estimators are well behaved and agree with v3 to within
~5 ppm, which matches the +5.03 ± 4.71 ppm trend in the ȳ difference.

**Conclusion: do not move E0's primary estimator onto v4 as written.** Either
keep v3 as the E0 baseline, or change the estimator (OLS or Theil–Sen, or
χ²-scaled weights with a variance floor) before using v4.

### 5.3 E2 validated against v4 instead of v3

Re-run with `E2_ISO_RESULTS=iso230_2_v4_results_blind.json`:

| | v3 corrected | v4 |
|---|---:|---:|
| raw E / M / A / R | 7.855 / 6.463 / 9.151 / 3.617 | 8.520 / 6.880 / 9.250 / 3.620 |
| after G: E / M | 3.317 / 1.335 | 3.994 / 1.925 |
| after G+B: E / A | 1.560 / 3.353 | 2.408 / 3.544 |
| residual up−down after G+B | +0.157 µm (rms 0.399) | **−0.061 µm (rms 0.418)** |
| mean-deviation rms after G+B | 0.344 µm | 0.496 µm |

v4 is the **better** test of E2's reversal model — the residual mean drops from
+0.157 to −0.061 µm — while the predicted improvement is slightly weaker
(E 8.52 → 2.41 µm, −72 %). Both sessions support the model; neither contradicts
it once v3 is read correctly.

**E2 does not need v4.** Swapping it in strengthens the reversal claim and
softens the E-improvement claim, both modestly. The reason to run the
cross-check against both is that it is cheap and shows the result is not an
artefact of one session.

---

## 6. State on disk

Changed:

- `Geometric Fitting/results_v3_corrected/` — new, Step A from the corrected v3
- `Geometric Fitting/E0/rendered/F1_*.png`, `F2_*.png` — re-rendered
- `experimental_scripts/e2_global_map/render_e2_findings.py` — two input paths,
  now env-overridable; nothing else
- `experimental_scripts/e2_global_map/analysis/*` and `rendered/*` — regenerated

Not changed:

- `iso230-2/analysis/iso230_2_evaluate.py` — still has the 51-sample detector
- `iso230-2/analysis/iso230_2_v3_results.json` — still the defective table
- `Geometric Fitting/results/` and `E0/TABLES.md` — still the old numbers
- E1 and E4 — see below

## 7. E1 and E4

Both subtract the Step A gain and then remove one offset per window, so only the
within-window tilt from the 2.34 ppm slope change survives:

| | windows | tilt per window | leakage into a 2 mm harmonic |
|---|---|---:|---:|
| E1 | 3 × 12 mm | 28 nm p-p | 1.5 nm |
| E4 | 7 × 6 mm | 14 nm p-p | 1.5 nm |

Against E1's own combined detection floor of 0.03 µm and its 0.1 µm
"not detected below" threshold, this is at or under the floor. No harmonic
amplitude or model selection moves. The only stale item is the caption in
`render_e1_findings.py:301`, which prints the Step A gain as −25.09 ppm.

**Neither needs re-rendering for correctness.**

## 8. The K_C refit — penalised splines, not K_C = 2–3

The stiffness in 60–78 mm was real. K_C = 1 was never the best varying-
coefficient model; it was only the cheapest one the one-SE rule let through.

### 8.1 Why K_C = 1 got selected

The single blocked-CV partition scores eight interior folds, and their spread
set the standard error the rule compares against: ±0.044 µm. That band was
wide enough to admit K_C = 1 even though it was the *worst* M_C variant on CV.
Repeating the CV over five shifted sets of block edges measures the same
quantity more precisely — SE drops to ±0.012 µm — and K_C ≤ 2 is then excluded
everywhere.

### 8.2 The comparison (H = 3, K_s = 2 family, repeated CV)

| model | size | CV | end-fold low | end-fold high | window rms | window C1+C2 | zone C1+C2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| K_C = 1 | 28 | 0.4820 | 1.57 | 1.04 | 0.4172 | 0.3191 | 0.3433 |
| K_C = 2 | 32 | 0.4679 | 2.30 | 0.87 | 0.4052 | 0.2874 | 0.3000 |
| K_C = 3 | 36 | 0.4482 | 5.27 | 0.78 | 0.3943 | 0.2326 | 0.1676 |
| K_C = 5 | 44 | 0.4513 | **22 532** | 3.24 | 0.3933 | 0.2314 | 0.1532 |
| **penalised, K_C = 8, λ = 0.32** | 42.6 edf | **0.4454** | **1.27** | **0.67** | **0.3928** | 0.2316 | 0.1579 |

Two things decide it:

- **Fixed knots trade interior accuracy against edge stability.** The end-fold
  columns are a held-out *end* block, so the model has to predict beyond its
  training support. K_C = 3 already costs 3× K_C = 1 there; K_C = 4 and 5 are
  not merely worse but unusable (401 µm and 22 532 µm).
- **The penalty removes the trade.** A rich K_C = 8 basis with a second-
  difference penalty gets the interior accuracy of K_C = 5 *and* better edge
  behaviour than K_C = 1, at fewer effective parameters than K_C = 5.

So: **penalised splines**, not K_C = 2–3.

### 8.3 What was implemented

In `render_e2_findings.py` only:

- `coefficient_penalty()` — second-difference penalty on the C1/C2 coefficient
  spline blocks. The offset, gain, slow trend and stationary higher harmonics
  are untouched.
- `fit_linear(..., penalty=λ)` — augmented least squares; effective degrees of
  freedom from the hat-matrix trace; BIC and the error variance both use edf.
  With `penalty=None` the arithmetic is identical to before.
- `cv_assignments()` — blocked CV repeated over five shifted edge sets.
- `blocked_cv()` — three families now: `M_B` (stationary), `M_C_fixed`
  (fixed-knot, kept as the diagnostic behind panel (a)) and `M_C` (the
  penalised family, which is what gets selected). λ replaces K_C as the
  flexibility knob.
- `one_se()` — unchanged in spirit: cheapest model within one SE of the best.
  Cost is now *effective* parameters, and a stability guard rejects any model
  whose held-out end block exceeds `CV_END_FACTOR` = 5× its interior score.

**The stability guard is not cosmetic.** Without it the rule selected
H = 4, K_s = 7, λ = 3.2 — 39.2 edf, good interior CV, and an end-fold error of
**440 µm**. The penalty covers the C1/C2 curves but not the slow trend, so a
large K_s buys flexibility the interior folds never charge for. Passing models
sit near 3× their interior score; failures run 60× to 20 000×, so the cut is
not delicate.

### 8.4 Selected model and what moved

`M_C`: H = 4, K_s = 2, K_C = 8, **λ = 1.0**, 39.2 effective parameters (58
columns). M_B is unchanged at H = 2, K_s = 2.

| | before (K_C = 1) | after (penalised) |
|---|---:|---:|
| interior CV | 0.4872 ± 0.0436 | **0.4340 ± 0.0123** |
| end-fold, low / high | 1.385 / 1.738 | **1.214 / 0.629** |
| in-sample residual rms | 0.3966 µm | **0.3685 µm** |
| scaled condition number | 5.89 | **2.06** |
| E2 gain | −29.895 ppm | −29.806 ppm |
| independent-window rms | 0.3992 | **0.3776 µm** |
| independent-window \|C1\|+\|C2\| error | 0.3141 | **0.2446 µm** |
| …the same, 60–78 mm only | 0.3322 | **0.1957 µm** |
| gap closed vs the local floor | 0.896 | **0.922** |
| ISO prediction, E after G+B | 1.560 | 1.482 µm |
| ISO prediction, M after G | 1.335 | 1.218 µm |

The worst single case was the E4 window at 78 mm, right in the zone: its
harmonic error falls from 0.586 µm to 0.263 µm.

Two windows get slightly worse on C1+C2 — E4 @ 95 (0.108 → 0.147) and
E4 @ 112 (0.413 → 0.432). Both are outside the zone and the net is strongly
positive.

### 8.5 Caveats

- `max_abs_slow_vs_C_column_correlation` rises from 0.096 to 0.138. Still small,
  but the richer coefficient basis overlaps the slow trend slightly more.
- Figure 2 panel (d) now auto-scales to roughly −150…+50° because the k = 2
  phase band widens sharply past ~115 mm. That region is inside the endpoint
  shading, so it is flagged, but the axis is stretched by it.
- λ was chosen from a 9-point grid spanning 10⁻²…10². The optimum (1.0) is
  interior to the grid, so the grid is not binding.

## 9. Open items

- `iso230_2_evaluate.py` still contains the defect; anything re-run through it
  will reproduce the wrong B/R/E/A.
- `E0/TABLES.md` T1/T2 and `results/` need updating, or F1 and the tables will
  keep disagreeing.
- `E0/TABLES.md` T4 lists three checks as "NOT RUN — raw capture unavailable to
  the analysis". That is no longer true; the capture was used for all of the
  above.
- `render_command_recorded_overlay.py` and `analyze_open_loop_diagnostic.py`
  were never audited for the same windowing assumption.

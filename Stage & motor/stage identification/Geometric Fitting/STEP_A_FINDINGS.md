# Step A — findings

Linear geometric identification from the ISO 230-2 v3 baseline. Compiled
2026-10-27 from `results/step_a_results.json`,
`results/per_target_statistics.csv` and `rendered/iso_v3_step_a_a4.png`,
produced by `scripts/process_iso_v3_step_a.py`.

**Status: performed, with three of the specified analyses not run** because the
raw capture was unavailable to the analysis machine. See
[§7](#7-what-the-specification-asked-for-and-did-not-get).

---

## 1. Headline

The direction-averaged geometric error carries a linear gain term of

$$G_\mathrm{linear}(x) = 0.17835 - 0.0250866\,(x - 61.1013)\ \mu\mathrm{m}$$

| Quantity | Value |
|---|---|
| Gain slope a₁ | **−0.025087 ± 0.013715 µm/mm** (95 %) |
| In ppm | **−25.09 ± 13.72 ppm** |
| Over the 112.500 mm span | **−2.822 µm** |
| Significance | 3.8 σ (a₁ / SE = 0.025087 / 0.006613) |
| Reduced χ² | **153.8** |
| Share of direction-averaged error power | **34.3 %** |
| Share of total budget (incl. reversal, repeatability) | 33.4 % |

The sign convention is that positive deviation means measured travel exceeded
command, so the stage runs **short** by 25 ppm — it under-travels progressively
along +x.

## 2. The finding that matters most

**The linear term is the smaller half of the problem.** Removing it takes the
direction-averaged RMS from 1.367 µm to 1.108 µm. The residual that remains is
not noise:

| Residual after gain removal | Value |
|---|---|
| RMS | 1.108 µm |
| Peak-to-peak | **4.042 µm** |
| Most negative | −2.481 µm at 112.0 mm |
| Most positive | +1.561 µm at 73.6 mm |

The residual's 4.042 µm peak-to-peak **exceeds the 2.822 µm that the gain
contributes over the whole span**. A linear model is therefore not a useful
standalone description of this stage's geometry, and the reduced χ² of 153.8
says so formally: the target-to-target structure is √153.8 ≈ 12× larger than
the directional-mean repeatability predicts. Panel (b) of the render shows the
structure plainly — it is smooth and position-locked, not scatter.

This is the quantitative justification for E1.

## 3. Validity of the underlying session

| Check | Result |
|---|---|
| ISO stops recorded | 240 (24 targets × 2 directions × 5 runs) — complete |
| `ISO_DONE` marker | present; 5 runs, 29.9 min |
| Abort / watchdog / fault / lost-step | none |
| Targets on the 1/16-microstep grid | all |
| Record window | 500 ms, uniform |
| Settle | 2000–2002 ms |
| Reference stops | 8 |
| Driver | MRES 16, IRUN = IHOLD = 11 (400 mA), TPOWERDOWN 0, IHOLDDELAY 0, intpol 0 |

IRUN equals IHOLD with both delays zero, so there was **no run-to-hold current
transition inside any record window** — the measured settled position is not
contaminated by a current step. That is worth keeping: it is the condition E1
must also run under, and `e1_config.json` does preserve it (400/400, 0, 0).

Dataset: 24 targets, 5.000–117.500 mm, span 112.500 mm, mean spacing 4.89 mm.
Note this differs from the Step A specification's stated 27 targets over
5–135 mm; the v3 session as actually run is what was analysed.

## 4. Robustness of the slope

| Variant | Slope (ppm) | 95 % CI (ppm) |
|---|---:|---:|
| Weighted, direction-averaged (primary) | −25.09 | ±13.72 |
| OLS cross-check | −24.17 | ±14.24 |
| Up direction only | −26.79 | ±12.52 |
| Down direction only | −34.56 | ±11.61 |
| Leave-one-target-out range | −28.07 … −20.67 | — |

The OLS cross-check agrees with the weighted fit to 0.9 ppm, and every
leave-one-target-out slope falls inside the primary 95 % interval, so no single
target drives the result.

**Directional difference is 7.76 ppm and is not established** — the two
intervals overlap substantially. Do not read a direction-dependent gain from
this.

One caution when quoting the direction fits: they are weighted by their own
per-direction variances (median s↑ = 0.229 µm, s↓ = 0.324 µm), which differ
enough that the two direction slopes average to −30.7 ppm rather than the
−25.09 ppm of the direction-averaged fit. That is a weighting artefact, not a
contradiction — OLS on the same data averages consistently at −24.17 ppm. Quote
the primary fit, not the mean of the two directions.

## 5. Reversal (characterised, not modelled)

| Quantity | Value |
|---|---|
| Mean B | −0.365 µm |
| SD of B | 0.369 µm |
| Max \|B\| | 0.982 µm |
| Trend of B with position | +1.66 ± 4.32 ppm |

The reversal trend is **indistinguishable from zero**, so backlash does not
measurably vary along the travel. Mean reversal is negative and small relative
to the geometric error — roughly a quarter of the residual's peak-to-peak. This
is input for the friction work, not for G.

## 6. Predicted effect of applying the gain correction

| Metric | Measured | After linear correction |
|---|---:|---:|
| E | 6.611 µm | **4.203 µm** |
| M | 6.364 µm | **4.042 µm** |
| A | 7.949 µm | **5.682 µm** |
| B | 0.982 µm | 0.982 µm |
| R | 2.427 µm | 2.427 µm |
| median s | 0.280 µm | 0.280 µm |
| max s | 0.524 µm | 0.524 µm |

B and R are unchanged by construction — a position-only correction cannot touch
reversal or repeatability. This is the correct behaviour and a useful check
that the correction was applied as a pure function of position.

E improves by 36 %, A by 29 %. Both remain far above the 0.28 µm repeatability
floor, which is the same point as §2: most of the error is still on the table
after the gain is removed.

## 7. What the specification asked for and did not get

`instructions/G_step_A_linear_processing.md` §13 lists nine deliverables. Six
were produced; three were not, all for the same reason.

| Deliverable | Status |
|---|---|
| Validity table | partial — sequence, grid, watchdog and driver done |
| a₁ ± CI, χ²_ν, OLS cross-check | done |
| a₁ per direction, leave-one-out | done |
| Explained share of error power | done |
| B_i vs x with trend | done |
| Predicted ISO metrics | done |
| Preliminary periodogram + folded residual | done |
| Plots | done — panels (a)–(g); `x_ref` vs time not possible |
| **Drift-corrected fit (§3, "do both, report both")** | **not run** |
| **a₁ per run (§6)** | **not run** |
| **Physical attribution table (§8)** | **not run** |

The first two were blocked by the absence of `ISO230-2_v3.csv` on the machine
that ran the analysis; the third by missing environment records (air and stage
temperature, pressure, interferometer compensation state) and the screw accuracy
class.

**The raw capture is not actually missing.** It is present in this workspace as
`stage identification/ISO230-2_v3.csv`, 43.8 MB, dated 2026-10-22, untracked by
git — which is why the Mac could not see it. The drift-corrected fit, per-run
slopes, per-stop settling checks and in-window saturation checks are all
runnable without new measurements.

That matters for one specific reason. §3 of the specification says that if a₁
differs between the drift-corrected and uncorrected fits by more than its CI,
drift is significant and must be discussed, and that **the corrected version is
the one to report as primary**. The number quoted throughout this document is
the uncorrected fit. Until the drift-corrected fit is run, −25.09 ppm is
provisional.

Physical attribution remains genuinely blocked. −25 ppm is within the range that
thermal expansion alone could produce — for steel at ~11 ppm/K, a little over
2 K of uncompensated stage temperature offset would do it — so the slope cannot
yet be apportioned among thermal expansion, interferometer scale error, cosine
misalignment and true screw lead deviation. Without that split it is not safe to
compensate the gain in firmware, because the thermal part is not a fixed
property of the machine.

## 8. Periodic content — screening only

The largest Lomb–Scargle peak sits at a **3.28 mm spatial period**
(0.3045 cycles/mm, normalised power 0.54).

Treat this as a hint and nothing more:

- 24 targets over 112.5 mm is a mean spacing of 4.89 mm, so the mean-sampling
  Nyquist period is about 9.8 mm. A reported peak at 3.28 mm is well inside the
  aliased region. Irregular spacing lets Lomb–Scargle probe below that limit in
  principle, but not reliably with 24 points.
- 3.28 mm is **not** a harmonic of the nominal 2 mm lead (that family is 2.000,
  1.000, 0.667 mm …), so it does not correspond to any obvious screw order.
- Panel (e), the residual folded onto the 2 mm lead, shows a trend across phase
  rather than a clean closed cycle — consistent with the dominant residual
  structure being longer-period than the lead, not screw-synchronous.

No harmonic count H for model B should be set from this. That is exactly what E1
is designed to deliver, at 0.12 mm spacing with 50 distinct screw phases per
window.

## 9. Conclusions

1. A real linear gain of **−25.09 ± 13.72 ppm** is present, significant at
   3.8 σ and robust to target removal and to the choice of estimator.
2. It explains **34.3 %** of direction-averaged spatial error power. Correcting
   it would take E from 6.61 to 4.20 µm and A from 7.95 to 5.68 µm.
3. **The remaining 66 % is structured, not random** — 4.04 µm peak-to-peak,
   χ²_ν = 153.8. The linear model is inadequate on its own, which is the result
   that justifies proceeding to E1 and model B.
4. Reversal is small (max 0.982 µm) and has no measurable position trend.
5. The slope is **provisional** until the drift-corrected fit is run, and
   **not yet attributable** to a physical cause, so it should not be programmed
   as a compensation constant.

## 10. Recommended next actions

1. Re-run `scripts/process_iso_v3_step_a.py` on this machine, where
   `ISO230-2_v3.csv` is present, to obtain the drift-corrected fit, per-run
   slopes and per-stop settling checks. Compare a₁ against the uncorrected
   −25.09 ppm using the §3 criterion.
2. Commit `ISO230-2_v3.csv` (LFS is already configured for `*.csv`) so the
   analysis is reproducible from a clean clone and the Mac stops reporting it
   missing.
3. Record air and stage temperature for E1 so the thermal share of any gain can
   be separated, which Step A could not do.
4. Run E1 to fix the harmonic count H. Do not set H from the 3.28 mm screening
   peak.

---

## Source files

| File | Contents |
|---|---|
| `results/step_a_results.json` | all fit coefficients, robustness, ISO metrics |
| `results/per_target_statistics.csv` | 24 rows: means, SDs, B, fit, residual |
| `results/STEP_A_SUMMARY.md` | the generated short summary |
| `rendered/iso_v3_step_a_a4.png` | panels (a)–(g), A4 portrait |
| `scripts/process_iso_v3_step_a.py` | the processing code |
| `instructions/G_step_A_linear_processing.md` | the specification |

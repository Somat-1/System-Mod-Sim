# E0 — tables

The quantities that read better as numbers than as marks. Source:
`../results/step_a_results.json` and `../results/per_target_statistics.csv`.
Sign convention: positive deviation means measured travel exceeded command.

---

## T1. The gain slope and its variants

| Estimator | Slope (µm/mm) | Slope (ppm) | 95 % CI (ppm) | SE (µm/mm) | χ²_ν |
|---|---:|---:|---:|---:|---:|
| **Weighted, direction-averaged (primary)** | **−0.025087** | **−25.09** | ±13.72 | 0.006613 | 153.8 |
| OLS cross-check | −0.024170 | −24.17 | ±14.24 | 0.006868 | — |
| Forward only (↑) | −0.026794 | −26.79 | ±12.52 | 0.006036 | 120.9 |
| Reverse only (↓) | −0.034558 | −34.56 | ±11.61 | 0.005598 | 84.3 |
| Leave-one-target-out, min | −0.028068 | −28.07 | — | — | — |
| Leave-one-target-out, max | −0.020673 | −20.67 | — | — | — |

Fitted form, position centred at x_c = 61.1013 mm to decorrelate offset
from slope:

$$G_\mathrm{linear}(x) = 0.17835 - 0.0250866\,(x - 61.1013)\ \mu\mathrm{m}$$

| Derived | Value |
|---|---:|
| Significance (a₁ / SE) | 3.79 σ |
| CI as a fraction of the estimate | ±54.7 % |
| Gain over the 112.500 mm span | −2.822 µm |
| Intercept at x = 0 | +1.711 µm |
| Forward − reverse slope difference | 7.76 ppm (CIs overlap; not established) |

Two cautions:

- **χ²_ν = 153.8.** The fit is formally a poor description; the CI above is
  χ²-scaled to account for it. Target-to-target structure is √153.8 ≈ 12×
  larger than directional-mean repeatability predicts.
- **Do not average the two direction slopes.** They carry different
  per-direction weights (median s↑ = 0.229 µm, s↓ = 0.324 µm), so their mean is
  −30.7 ppm, not the −25.09 ppm of the primary fit. OLS on the same data
  averages consistently at −24.17 ppm. Quote the primary row.

---

## T2. ISO 230-2 metrics, measured and predicted after gain correction

| Metric | Meaning | Measured (µm) | After gain (µm) | Change |
|---|---|---:|---:|---:|
| E↑ | Systematic deviation, forward | 6.116 | 3.918 | −35.9 % |
| E↓ | Systematic deviation, reverse | 6.611 | 4.203 | −36.4 % |
| **E** | **Systematic deviation, bidirectional** | **6.611** | **4.203** | **−36.4 %** |
| **M** | **Mean bidirectional range** | **6.364** | **4.042** | **−36.5 %** |
| **A** | **Accuracy** | **7.949** | **5.682** | **−28.5 %** |
| **B** | **Reversal (max \|B_i\|)** | **0.982** | **0.982** | **0** |
| B̄ | Mean reversal | −0.365 | −0.365 | 0 |
| R↑ | Repeatability, forward | 2.098 | 2.098 | 0 |
| R↓ | Repeatability, reverse | 1.850 | 1.850 | 0 |
| **R** | **Repeatability, bidirectional** | **2.427** | **2.427** | **0** |
| s median | Per-target SD, median | 0.280 | 0.280 | 0 |
| s max | Per-target SD, worst | 0.524 | 0.524 | 0 |

B, R and s are unchanged by construction — a position-only correction cannot
touch reversal or repeatability. That they come through identical is the check
that the correction was applied as a pure function of position.

Rendered as [F2_metric_improvement.png](rendered/F2_metric_improvement.png).

---

## T3. Error power explained by the gain

| Quantity | Value |
|---|---:|
| Direction-averaged RMS, before | 1.3669 µm |
| Direction-averaged RMS, after | 1.1083 µm |
| **Share of direction-averaged error power explained** | **34.3 %** |
| Share of the total budget (incl. reversal and repeatability) | 33.4 % |

What survives the gain:

| Residual r_i = ȳ_i − G_linear(x_i) | Value |
|---|---:|
| RMS | 1.108 µm |
| **Peak-to-peak** | **4.042 µm** |
| Most negative | −2.481 µm at 112.0 mm |
| Most positive | +1.561 µm at 73.6 mm |
| Mean error bar on ȳ_i (±2σ) | ±0.189 µm |

The residual's 4.042 µm peak-to-peak **exceeds the 2.822 µm the gain
contributes over the whole span**, and its 1.108 µm RMS is roughly 12× the
0.095 µm (1σ) uncertainty on each point. The remaining error is structure, not
scatter. That is the case for E1.

---

## T4. Validity of the underlying session

| Check | Result | Status |
|---|---|---|
| ISO stops recorded | 240 = 24 targets × 2 directions × 5 runs | PASS |
| Completion marker | `ISO_DONE,runs=5,stops=240`, 29.9 min | PASS |
| Abort / watchdog / fault / lost step | none found | PASS |
| Targets on the 1/16-microstep grid | all | PASS |
| Record window | 500 ms, uniform | PASS |
| Settle window | 2000–2002 ms | PASS |
| Reference stops | 8 | PASS |
| Runs present | 0, 1, 2, 3, 4 | PASS |
| Directions present | both | PASS |
| Per-stop settling test | raw capture unavailable to the analysis | **NOT RUN** |
| x_ref drift correction | raw capture unavailable to the analysis | **NOT RUN** |
| Per-run slopes | raw capture unavailable to the analysis | **NOT RUN** |
| Environment (air/stage T, pressure, compensation) | never recorded | **NOT AVAILABLE** |

### Driver state during the session

| Register | Value | Why it matters |
|---|---:|---|
| MRES | 16 | 1 microstep = 625 nm at 100 full steps/mm |
| intpol | 0 | no MicroPlyer interpolation |
| IRUN | 11 (400 mA) | — |
| IHOLD | 11 (400 mA) | **equal to IRUN** |
| TPOWERDOWN | 0 | — |
| IHOLDDELAY | 0 | — |

IRUN equals IHOLD with both delays zero, so **no run-to-hold current transition
occurred inside any record window**. The settled position is not contaminated
by a current step. E1 preserves this condition (`e1_config.json`: 400/400 mA,
TPOWERDOWN 0, IHOLDDELAY 0).

### Dataset

| Property | Value |
|---|---:|
| Targets | 24 |
| Range | 5.000 – 117.500 mm |
| Span | 112.500 mm |
| Mean target spacing | 4.89 mm |
| Runs per direction | 5 |

Note this differs from the Step A specification's stated 27 targets over
5–135 mm. The v3 session as actually run is what was analysed.

---

## T5. Preliminary periodic content — screening only

| Quantity | Value |
|---|---:|
| Largest Lomb–Scargle peak | 3.28 mm spatial period |
| Frequency | 0.3045 cycles/mm |
| Normalised power | 0.54 |
| Mean-sampling Nyquist period | ≈ 9.8 mm |

**Do not set the harmonic count H from this.** The peak sits well inside the
aliased region for 24 points at 4.89 mm mean spacing, and 3.28 mm is not a
harmonic of the nominal 2 mm lead (that family is 2.000, 1.000, 0.667 mm …).
Fixing H is E1's job.

# Floor versus optical-table SEIS noise, 2–14 minutes

## Scope and method

The two 20 kHz recordings were compared over their common 2–14 minute interval. The operator-confirmed 26–30 minute disturbance in `TvalSeis.csv` does not overlap this interval and is therefore excluded by construction.

No accelerometer sensitivity or full-scale range is present in the project, so all results are reported as recorded channel voltage (mV), not m/s². Baseline noise is the median of locally demeaned, non-overlapping 1 s RMS values; no frequency filter was applied.

## Baseline noise comparison

| Axis | Floor median RMS (mV) | Table median RMS (mV) | Table/floor | Floor 95th pct (mV) | Table 95th pct (mV) |
|---|---:|---:|---:|---:|---:|
| X | 1.35374 | 0.212044 | 0.1566 | 1.55943 | 0.336823 |
| Y | 1.33895 | 0.254412 | 0.19 | 1.50344 | 0.385853 |
| Z | 2.31774 | 0.362895 | 0.1566 | 2.52993 | 0.559557 |
| Vector √(X²+Y²+Z²) | 3.00613 | 0.497024 | 0.1653 | 3.2238 | 0.71609 |

A ratio below 1 means the optical-table recording has the lower baseline noise.

## Integrity and clipping checks

| Recording | Axis | Raw min (mV) | Raw max (mV) | Samples at exact extrema | Longest identical run | Non-finite samples | Digital clipping flag |
|---|---|---:|---:|---:|---:|---:|---|
| Floor | X | -12.71075 | 13.10551 | 2 (1.3889e-05%) | 3 | 0 | none detected |
| Floor | Y | -11.02793 | 13.86547 | 2 (1.3889e-05%) | 3 | 0 | none detected |
| Floor | Z | -19.41973 | 18.97961 | 2 (1.3889e-05%) | 2 | 0 | none detected |
| Optical table | X | -6.232798 | 5.548298 | 2 (1.3889e-05%) | 3 | 0 | none detected |
| Optical table | Y | -4.111171 | 6.742716 | 2 (1.3889e-05%) | 3 | 0 | none detected |
| Optical table | Z | -9.512901 | 7.039011 | 2 (1.3889e-05%) | 3 | 0 | none detected |

No digital clipping or flat-top signature was detected by the exact-extrema and repeated-code checks. Analog sensor saturation cannot be ruled out without the sensor full-scale range or calibration.

## Sampling and drift checks

| Recording | Samples | Fs (Hz) | Min Δt (s) | Max Δt (s) | Irregular Δt count |
|---|---:|---:|---:|---:|---:|
| Floor | 14400000 | 20000 | 5e-05 | 5e-05 | 0 |
| Optical table | 14400000 | 20000 | 5e-05 | 5e-05 | 0 |

| Recording | Axis | Linear mean slope (mV/min) | Last–first 60 s shift (mV) | |shift| / median RMS |
|---|---|---:|---:|---:|
| Floor | X | -8.23817e-06 | -0.00122709 | 0.0009064 |
| Floor | Y | 1.20119e-05 | 0.00110522 | 0.0008254 |
| Floor | Z | -3.7366e-06 | 0.00156057 | 0.0006733 |
| Optical table | X | 8.15591e-07 | 0.000198026 | 0.0009339 |
| Optical table | Y | 3.46304e-06 | -0.000753082 | 0.00296 |
| Optical table | Z | -4.63112e-06 | 0.000869067 | 0.002395 |

The drift table operates on 1 s channel means. It is intentionally kept separate from AC RMS so slow offset drift is not mislabeled as broadband noise.

### RMS non-stationarity

| Recording | Axis | Linear RMS slope (mV/min) | Last–first 60 s median RMS (mV) |
|---|---|---:|---:|
| Floor | X | 0.0182249 | 0.242081 |
| Floor | Y | 0.00929414 | 0.150092 |
| Floor | Z | -0.0226182 | -0.191394 |
| Optical table | X | 0.0012849 | 0.0159764 |
| Optical table | Y | 0.000462013 | 0.0295921 |
| Optical table | Z | -0.00413366 | 0.00827281 |

This second table describes change in vibration/noise amplitude over the segment, rather than electrical zero-offset drift.

## Outputs

- `plots/seis_time_domain_rms_02_14min.png`: overlaid time histories and RMS trends.
- `plots/seis_time_domain_rms_02_14min.pdf`: publication-ready vector version.
- `seis_time_rms_02_14min.csv`: 1 s means and locally demeaned RMS values.
- `cache/seis_time_rms_02_14min.npz`: compact plotting cache.
- `exclusion_windows.csv`: retained 26–30 minute disturbance definition.

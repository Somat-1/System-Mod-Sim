# Step A result — ISO 230-2 v3 linear geometric fit

## Result

The uncorrected, direction-averaged weighted fit is

\[
G_\mathrm{linear}(x) = 0.17835 - 0.0250866(x - 61.1013)\ \mu\mathrm{m},
\]

where position is in millimetres. The slope is **−0.02509 ± 0.01372
µm/mm** at 95% confidence, or **−25.09 ± 13.72 ppm**. Across the measured
112.500 mm span, this corresponds to **−2.822 µm** of linear gain error.

The OLS cross-check gives −24.17 ± 14.24 ppm and therefore agrees with the
weighted estimate. The reduced χ² is 153.8: target-to-target structure is far
larger than directional-mean repeatability predicts, so a linear model alone is
not an adequate description of the geometry.

## Effect of linear correction

| Quantity | Measured | Predicted after gain correction |
|---|---:|---:|
| Direction-averaged RMS about its mean | 1.367 µm | 1.108 µm |
| E | 6.611 µm | 4.203 µm |
| M | 6.364 µm | 4.042 µm |
| A | 7.949 µm | 5.682 µm |
| B | 0.982 µm | 0.982 µm |
| R | 2.427 µm | 2.427 µm |

The fitted gain explains 34.3% of direction-averaged spatial error power and
33.4% of the reconstructed total power including reversal and repeatability.
B and R remain unchanged, as required for a position-only correction.

## Robustness and reversal

- Directional slopes are −26.79 ppm (up) and −34.56 ppm (down). Their
  confidence intervals overlap; the data do not establish a direction-dependent
  gain.
- Leave-one-target-out slopes range from −28.07 to −20.67 ppm, all inside the
  main fit's 95% interval.
- Mean reversal is −0.365 µm; maximum absolute reversal is 0.982 µm.
- The reversal trend is +1.66 ± 4.32 ppm and is not distinguishable from zero.

The largest preliminary Lomb–Scargle peak corresponds to a 3.28 mm spatial
period. This is screening evidence only: 24 irregular targets are too sparse
to assign a physical period confidently, and the nominal 2 mm-lead fold does
not by itself establish a screw-periodic term.

## Validity and limitation

The controller log records all 240 intended ISO stops (24 targets × two
directions × five runs), an `ISO_DONE` marker, integer 1/16-microstep target
commands, and no abort, fault, watchdog or lost-step marker. IRUN equals IHOLD
(11), with TPOWERDOWN and IHOLDDELAY both zero, so there was no run-to-hold
current transition in the record windows.

`ISO230-2_v3.csv` is absent. Consequently, drift-corrected fitting, direct
x_ref drift estimation, per-stop settling tests, saturation checks within each
record window, and per-run slopes were not recomputed. Environmental records
and the screw accuracy class are also unavailable, so the −25 ppm slope cannot
yet be apportioned defensibly among thermal expansion, interferometer scale,
cosine alignment and screw lead deviation.

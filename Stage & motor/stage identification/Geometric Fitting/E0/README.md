# E0 — the ISO 230-2 v3 baseline, visualised

E0 is not a new experiment. It is the presentation layer for the existing ISO
230-2 v3 session and the Step A fit: the figures and tables that carry the
result into the write-up, and the argument for why E1 is needed.

```
E0/
├── render_e0_figures.py   regenerates both figures from ../results/
├── TABLES.md              T1 slope variants · T2 ISO metrics · T3 power share
│                          T4 validity · T5 periodic screening
└── rendered/
    ├── F1_deviation_residual_reversal.png
    └── F2_metric_improvement.png
```

## Regenerate

```bash
python render_e0_figures.py
```

It reads `../results/per_target_statistics.csv` and
`../results/step_a_results.json` and writes both figures. Nothing is recomputed
— the fit coefficients come from `process_iso_v3_step_a.py`, so the figures
cannot drift away from the numbers in `TABLES.md`.

## Figure 1 — the core result

Three panels on one shared position axis, which is the point: the reader sees
that (b) has structure while (c) does not.

| Panel | Shows |
|---|---|
| (a) | x̄_i↑ and x̄_i↓ with ±2s (run-to-run) shading, the direction average ȳ_i in black, and G_linear(x) fitted to those black points |
| (b) | residual r_i = ȳ_i − G_linear(x_i) with ±2σ-of-the-mean error bars |
| (c) | reversal B_i with ±2σ-of-the-mean, the mean (−0.365 µm) and the fitted trend (+1.66 ± 4.32 ppm) |

Markers are not joined. The targets are ~5 mm apart with an expected 2 mm
periodic term between them, so a line across that gap would assert smoothness
the data cannot support. Note the error bars differ by panel: (a) is ±2s, the
run-to-run spread, while (b) and (c) are ±2σ of the mean, smaller by √n.

The reading: a −25 ppm gain is real and visible in (a), but (b) shows a 4.04 µm
peak-to-peak residual against ±0.19 µm error bars — structure roughly 12× the
measurement uncertainty, which is what χ²_ν = 154 means. Meanwhile (c) is flat
within its error bars, so reversal carries no position dependence and friction
stays a separate model from G.

## Figure 2 — before and after the gain

Grouped bars for E, M, A, B and R, measured against predicted after a
position-only gain correction. E and M fall by ~36 %, A by 29 %; B and R are
identical by construction. This is the first bar of the improvement chart that
E1/E2 and then the friction work extend.

## Caveat carried from Step A

The slope is **provisional**. The drift-corrected fit the specification asks
for as primary was never run, because `ISO230-2_v3.csv` was not visible to the
machine that produced the results — it is present in this workspace. See
`../STEP_A_FINDINGS.md` §7.

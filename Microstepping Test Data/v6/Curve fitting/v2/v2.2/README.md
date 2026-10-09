# v2.2 — first 60 s of the idle tail: detrended noise and its spectrum

[v2](../README.md) characterizes the whole 153.65 s post-recording tail and
establishes what is in it: a boundary transient, slow wander, and a ~5-9 nm
sample-to-sample noise floor. v2.2 narrows to the **first 60 s** and pushes on
one question: with a low-order trend removed, what is the residual, how is it
distributed, and where does its power sit in frequency?

```
v2.2/
├── README.md                              this file
├── tail_60s_summary.json                  every number below, machine-readable
├── scripts/
│   └── analyze_tail_60s.py                the whole analysis
└── rendered_assets/
    ├── detrend_overlay.png                raw + fitted polynomial + residual + order selection
    ├── residual_distribution.png          residual with σ/p-p in the legend, and its histogram
    └── welch_psd.png                      Welch PSD + cumulative RMS
```

```powershell
cd "Microstepping Test Data/v6/Curve fitting/v2/v2.2/scripts"
python analyze_tail_60s.py
```

## The segment

| | |
|---|---|
| Tail boundary (from v2) | t = 1148.667 s |
| **Segment start** | **t = 1148.717 s** (boundary + 50 ms) |
| Segment end | t = 1208.716 s |
| Length | 60.000 s, 60 000 samples @ 1 kHz |
| Raw σ / peak-to-peak | 16.49 nm / 107.0 nm |

### Why 50 ms is dropped

The tail's opening contains a ~1.15 µm drop with no explanation available in
this dataset (v2's `large_outliers_gt_300nm`). It is not part of the noise
process and it dominates any statistic computed with it included, so it is
excluded — but the cut is **measured rather than assumed.** Growing the skip
and watching peak-to-peak settle:

| Skip | 0 ms | 5 ms | 10 ms | 50 ms | 100 ms | 500 ms | 1000 ms |
|---|---|---|---|---|---|---|---|
| p-p | 1271 nm | 275 nm | 219 nm | **107 nm** | 107 nm | 107 nm | 107 nm |

Peak-to-peak stops changing at 50 ms and is then identical out to 1 s, and the
last step-to-step jump above 100 nm is at 26 ms. 50 ms is where the transient
is fully out and nothing further is being thrown away.

## The subtracted polynomial

```
p(t) = -8.88226 - 0.685266·t + 0.0630326·t² - 0.000831767·t³
```

with **t in seconds from segment start** and **p in nanometres**. (Fitted by
least squares on a [−1, 1]-mapped domain for conditioning, then converted to
these plain-t coefficients; the segment is referenced to its own median, so
the constant term is a small offset rather than the ~−194 µm absolute stage
position.)

`detrend_overlay.png` shows it over the raw trace, the residual underneath,
and the order sweep that chose it.

### Why order 3

The order is selected, not assumed: sweep 0-8 and take the first order whose
successor buys less than 1% further reduction in residual σ.

| Order | 0 | 1 | 2 | **3** | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| Residual σ [nm] | 16.49 | 14.95 | 14.61 | **14.21** | 14.18 | 14.18 | 13.82 | 12.62 | 12.61 |
| Improvement | — | 9.3% | 2.3% | 2.7% | **0.2%** | 0.0% | 2.5% | 8.7% | 0.1% |

The improvement collapses to 0.2% at order 4, which is the knee. It then picks
up again at orders 6-7 — and that matters: **that is not a better trend, it is
the polynomial acquiring enough freedom to start tracing the slow wander
itself.** The wander is a result (v2 quantifies it), not a nuisance to be
fitted away; a residual with it removed would report a noise floor that is
partly an artifact of the detrend. Stopping at the first knee keeps the
subtraction a trend removal.

The cubic removes **13.8%** of the raw standard deviation — modest by
construction, and the point: what is left is still the wander plus the noise,
with only the smooth part of the trend gone.

## Residual statistics

`residual_distribution.png` — the residual time series with σ and peak-to-peak
annotated in the legend, ±1σ/±3σ direct-labelled, and the histogram beside it.

| | |
|---|---|
| **Standard deviation** | **14.21 nm** |
| **Peak-to-peak** | **102.4 nm** |
| Mean | ~0 (by construction) |
| Skew | +0.368 |
| Excess kurtosis | +0.273 |
| Within ±1σ | 69.8% (68.3% if normal) |
| Histogram bins | 116 (Freedman–Diaconis) |

The distribution is close to normal but not quite: a mild positive skew and a
slightly peaked centre, both visible against the overlaid Normal(0, σ) curve.
That is expected — **the residual is not a stationary white process.** It is
dominated by the slow wander, which the spectrum below makes explicit, and a
σ over 60 s of a wandering signal is a description of that wander, not an
encoder noise figure. v2's per-second statistics (~5-9 nm) remain the number
to quote for the noise floor.

## Welch PSD

`welch_psd.png` — Hann window, **16 384-sample segments (16.4 s)**, **50%
overlap**, 6 averaged segments, 0.0610 Hz resolution, per-segment constant
detrend.

Long segments are the right trade here: the feature of interest at the low end
is the wander, which lives below a few Hz and needs fine resolution to be
separated from DC at all. 16 384 still leaves 6 averages, which is enough to
keep the floor readable.

### Where the power is

| Band | RMS | Mean ASD |
|---|---:|---:|
| DC – 1 Hz | **11.92 nm** | 11.70 nm/√Hz |
| 1 – 10 Hz | 2.74 nm | 0.914 nm/√Hz |
| 10 – 100 Hz | 0.93 nm | 0.098 nm/√Hz |
| 100 – 300 Hz | 1.34 nm | 0.095 nm/√Hz |
| **300 – 400 Hz** | **4.34 nm** | 0.434 nm/√Hz |
| 400 – 500 Hz | 2.28 nm | 0.228 nm/√Hz |

Three features:

1. **The wander dominates completely.** Half the total RMS is below
   **0.08 Hz**, and the single lowest resolvable bin (0.0610 Hz) carries
   5.52 nm on its own. The PSD peaks at 0.122 Hz with 644 nm²/Hz — roughly
   **five orders of magnitude** above the broadband floor. This is the same
   phenomenon v2 attributes to thermal cycling, now localized in frequency.
2. **A genuine white floor between ~10 and 250 Hz**, median 0.0073 nm²/Hz =
   **0.085 nm/√Hz**. This is the part that behaves like sensor noise.
3. **A resonance at 357.9 Hz**, ~25× above the local floor — an independent
   confirmation of v2's 362.8 ± 3.0 Hz, found here on a different segment with
   a different detrend, agreeing to within 1.4%. The rise into Nyquist is also
   reproduced.

**The aliasing caveat from v2 applies unchanged.** Nyquist is 500 Hz and the
PSD climbs steeply into it — the signature of folding. 357.9 Hz is
indistinguishable from 642 Hz (or 1358 Hz…) at a 1 kHz sample rate. It is the
right number for designing a filter that runs *on this series*; it is not
evidence of a 358 Hz mechanical mode. Confirming that needs a faster capture.

### Parseval check, and why it is not exactly 1

Integrating the PSD gives **13.18 nm** against a time-domain σ of **14.21 nm**
— a ratio of 0.934. The 7% shortfall is expected and is not an error: Welch
removes each 16.4 s segment's own mean (`detrend="constant"`), which discards
power below ~0.06 Hz — and this signal has most of its variance exactly there.
The gap is a direct measure of how much of the residual lives below the
resolution the segment length can reach.

The cumulative-RMS panel shows this: the curve is already at ~6 nm in the
first plotted bin, reaches ~12 nm by 2 Hz, then runs essentially flat to
250 Hz before a visible step at the 358 Hz resonance.

# v6.1 — the standstill tail of v6's EXPrun, clipped and DC-removed

A derived dataset, not a new measurement. It takes [`../v6/EXPrun.csv`](../v6/EXPrun.csv)
— the full v6 MRES trajectory campaign — clips off the trailing stretch where
**nothing was commanded and the stage is standing still**, subtracts the DC
offset, and exports the result in a form that can go straight into a
frequency-domain analysis.

```
v6.1/
├── README.md                          this file
├── scripts/
│   ├── clip_standstill_tail.py        detection, clip, DC removal, export, render
│   └── welch_psd_tail.py              Welch PSD (Hann) over the whole clip
├── standstill_tail/
│   ├── standstill_tail.csv            t_s, pos_um, pos_um_ac  (152,651 rows)
│   ├── standstill_tail.npz            same arrays + clip metadata, as float64
│   ├── clip_summary.json              every clip number below, machine-readable
│   ├── standstill_tail.png            the time-domain render
│   └── standstill_tail.svg            same, vector
└── psd/
    ├── welch_psd_tail.png             the PSD plot
    ├── welch_psd_tail.svg             same, vector
    ├── welch_psd_tail.npz             f_hz, psd_nm2_per_hz, cumulative_rms_nm
    └── welch_summary.json             configuration + every PSD number below
```

Run from `scripts/`, in this order:

```powershell
python clip_standstill_tail.py
python welch_psd_tail.py
```

It reads `../../v6/EXPrun.csv` and reuses
[`../v6/analysis/scripts/ids_common.py`](../v6/analysis/scripts/ids_common.py)
for loading and motion detection, so the scale convention (**1 count = 1 nm**)
and the stillness test are the same ones the rest of v6 uses.

## The clip

The campaign's 24th and last trajectory leg decelerates into its final rest at
t = 1148.67 s, and the scope keeps recording to t = 1302.32 s. Nothing is
commanded in between — no marker, no oscillation, no leg. That makes it the
only segment in the capture with no settling inherited from a preceding move.

| | |
|---|---:|
| Source capture | 1302.317 s at 1 kHz, uniformly sampled, no gaps |
| Detected standstill start | t = 1148.667 s |
| Settling guard dropped | 1.0 s |
| **Exported clip** | **t = 1149.667 – 1302.317 s** |
| Duration | 152.650 s, **152,651 samples** |
| DC offset subtracted | **−194.321443 µm** |

The boundary is **detected, not assumed**: `ids_common.still_mask` (windowed
slope, the same test used elsewhere in v6) is scanned backward from the end of
the capture for the earliest sample after which the axis is continuously still
all the way to the last sample. It independently reproduces the t = 1148.667 s
that [`../v6/Curve fitting/v2`](../v6/Curve%20fitting/v2/README.md) found.

### Why the 1 s guard

The detected boundary's first ~4 ms carry a single ~1.15 µm settling transient
(four consecutive >300 nm sample-to-sample jumps, catalogued in
[`../v6/Curve fitting/v2/noise_summary.json`](../v6/Curve%20fitting/v2/noise_summary.json)).
Against a signal whose whole remaining range is 232 nm, that impulse is both an
amplitude outlier and — more to the point here — a broadband one: left in, it
smears energy across the entire spectrum and contaminates every bin. Dropping
1 s from the front removes it with room to spare, at a cost of 0.65 % of the
record.

It is removed, not hidden: panel 1 of the render marks the detected start and
the exported clip separately, so the gap between them is visible.

## What's in the exported columns

| Column | Meaning |
|---|---|
| `t_s` | time from the start of the **original capture** [s] — so a row is traceable back to `EXPrun.csv` |
| `pos_um` | raw position in the capture's own frame (referenced to its first sample) [µm] |
| `pos_um_ac` | `pos_um` minus the clip's mean — **this is the DC-removed series** [µm] |

`standstill_tail.npz` carries those three plus `t_rel_s` (time from the start of
the clip, i.e. starting at 0 — usually what you want as the time base),
`pos_nm_ac` (the DC-removed series in nm, since the signal is tens of nm), and
the clip metadata as scalars (`fs_hz`, `dt_s`, `dc_offset_um`, `clip_start_s`,
`clip_end_s`, `detected_start_s`, `guard_s`).

```python
import numpy as np
d = np.load("standstill_tail/standstill_tail.npz")
x, fs = d["pos_nm_ac"], float(d["fs_hz"])   # nm, 1000.0 Hz
```

## What the clipped signal looks like

`standstill_tail.png` / `.svg` — three stacked panels:

1. **Provenance.** The full 1302 s capture with the exported window shaded, so
   it is obvious which part of the campaign this came from and that the clip
   contains none of the motion.
2. **The clip, DC removed**, full 152.65 s. Plotted in **nm**, because ±50 nm
   drawn on a µm axis is a flat line at zero.
3. **One second of it at native 1 kHz resolution**, undecimated, with the
   sample markers drawn — the 1 nm counter quantisation and the character of
   the sample-to-sample noise are both visible here and nowhere else.

| Metric (DC-removed clip) | Value |
|---|---:|
| σ | 37.38 nm |
| Peak-to-peak | 232.0 nm |
| Mean \|Δsample\| | 6.81 nm |
| Median \|Δsample\| | 5.00 nm |
| Residual linear trend | +0.492 nm/s (+75.1 nm over the clip) |
| σ after also removing that linear trend | 30.45 nm |

**The DC offset is not the only low-frequency content.** Subtracting the mean
zeroes bin 0, but the panel-2 trace is visibly dominated by a slow wander of
±50–150 nm over tens of seconds — that is what the 37.38 nm σ is mostly
measuring, not sample-to-sample noise (6.81 nm). A linear detrend only takes it
from 37.38 to 30.45 nm, so the wander is not a straight line either. If the
interest is the higher-frequency content, expect to window or detrend further
before transforming, or that leakage from this wander will sit on top of
everything else.

## Spectral properties of the clip

| | |
|---|---:|
| Sample rate | 1000 Hz |
| Nyquist | 500 Hz |
| Record length | 152.650 s |
| Frequency resolution (1/T) | 0.00655 Hz |
| Amplitude quantisation | 1 nm (1 encoder count) |

**Caveat carried over from v6: the capture may be aliased.** The PSD of this
same segment rises steeply into Nyquist, the usual signature of folding, so a
feature measured at *f* here is indistinguishable from 1000 − *f*, 1000 + *f*,
… . A frequency read off this 1 kHz series is the right number for designing
something that *runs* on this 1 kHz series; it is not by itself evidence of a
mechanical mode at that frequency. Settling it needs a faster capture. See
[`../v6/Curve fitting/v2`](../v6/Curve%20fitting/v2/README.md).

## The Welch PSD

`welch_psd_tail.py` transforms the exported clip — the whole 152.65 s of it,
DC already removed, nothing further detrended.

| | |
|---|---|
| Window | **Hann** |
| `nperseg` | 16384 samples = 16.384 s |
| Overlap | 50 % (8192 samples) |
| Segments averaged | **17** |
| Resolution | 0.0610 Hz |
| `detrend` | `"constant"` (per segment) |
| Scaling | `density` → nm²/Hz |

`nperseg` is held at 16384 **to match
[`../v6/Curve fitting/v2/v2.2`](../v6/Curve%20fitting/v2/v2.2/README.md)**,
which ran the same estimator over the first 60 s of this tail. Same window,
same resolution, same per-segment detrend, so the two spectra differ because
of the data rather than the configuration. What the full tail buys at that
setting is averaging — **17 segments against v2.2's 6**, cutting the variance
of the estimate by roughly √(17/6) ≈ 1.7×. The alternatives are recorded in
`welch_summary.json` under `nperseg_trade`:

| `nperseg` | Segment | Δf | Averages |
|---:|---:|---:|---:|
| 4096 | 4.096 s | 0.244 Hz | 73 |
| 8192 | 8.192 s | 0.122 Hz | 36 |
| **16384** | **16.384 s** | **0.0610 Hz** | **17** |
| 32768 | 32.768 s | 0.0305 Hz | 8 |

### Results

`psd/welch_psd_tail.png` — the PSD on log-log axes, and beneath it the
cumulative RMS, which answers "which frequencies actually carry the motion"
in a way the PSD alone does not.

| | |
|---|---:|
| Peak above 10 Hz | **357.85 Hz**, 125× the floor |
| Broadband floor, 10–250 Hz | **0.082 nm/√Hz** (0.0067 nm²/Hz) |
| Wander peak | 0.0610 Hz — the lowest resolvable bin |
| Integrated RMS | 26.17 nm |
| 90 % / 99 % of RMS reached by | 0.244 Hz / 356 Hz |

| Band | RMS [nm] |
|---|---:|
| 0–1 Hz | 25.44 |
| 1–10 Hz | 3.64 |
| 10–100 Hz | 0.87 |
| 100–300 Hz | 1.29 |
| 300–400 Hz | 4.14 |
| 400–500 Hz | 2.22 |

Both headline numbers **reproduce v2.2 on 2.5× the record**: 357.85 Hz against
its 357.85 Hz (the same bin), and a 0.082 nm/√Hz floor against its 0.085.

This is a consistency check, **not an independent confirmation** — v2.2's 60 s
window is a subset of this clip, so the two estimates share data and the first
6 of these 17 segments largely are v2.2's. What it does establish is that the
peak and the floor are properties of the whole tail rather than of its opening
minute, and that neither moves when the averaging is nearly tripled. For a
genuinely independent check of the peak, the ring-down fits in
[`v2.3`](../v6/Curve%20fitting/v2/v2.3/README.md) are the ones to cite: they
use different data in a different motion state and land at 343.8–350.0 Hz.

### Two things to read carefully

**The 50 % point is not resolvable, and is not reported as a frequency.**
51 % of the total power sits in the two lowest bins (≤ 0.061 Hz). The wander
is slower than one 16.384 s segment, so this estimate can only say it is below
its own resolution — interpolating a "half the RMS below 0.024 Hz" out of the
first bin would be an artifact of the interpolation, not a measurement.
Resolving the wander itself needs longer segments, which this 152 s record
cannot supply many of. The quantiles in `welch_summary.json` carry a
`resolved` flag for exactly this reason.

**The Parseval ratio is 0.700, and that is expected.** The integrated spectrum
gives 26.17 nm RMS against the time-domain σ of 37.38 nm. `detrend="constant"`
removes each 16.384 s segment's own mean, which discards the wander *between*
segments — content below roughly 1/`nperseg`. The missing variance is that
wander, not an estimator error. If you want a number where the two agree,
detrend first, as v2.2 does (its ratio is 0.928).

## Related work on this same segment

This folder deliberately produces only the clipped, DC-removed series and a
look at it — the frequency analysis itself is not here. Some has already been
done on the same tail, and is worth reading before repeating it:

| Where | What it covers |
|---|---|
| [`../v6/Curve fitting/v2`](../v6/Curve%20fitting/v2/README.md) | noise-floor metrics, zoom cascade, interactive view, candidate smoothers; locates a peak at 362.8 ± 3.0 Hz |
| [`../v6/Curve fitting/v2/v2.2`](../v6/Curve%20fitting/v2/v2.2/README.md) | first 60 s only: polynomial detrend + Welch PSD; white floor 0.085 nm/√Hz over 10–250 Hz |
| [`../v6/Curve fitting/v2/v2.3`](../v6/Curve%20fitting/v2/v2.3/README.md) | is that peak a real mode? ring-down fits say yes — 343.8–350.0 Hz, ζ ≈ 0.022–0.029 |

## Plot conventions

Same as v6: measured series in `#1769aa`, annotation/region accents in
`#d94801`, recessive grid, long traces decimated with a min/max envelope so
extremes survive (panel 3 is undecimated by design).

# v2 — post-recording idle tail: a pure noise-floor measurement

`EXPrun.csv` runs 1302.3 s. The last commanded motion (leg #24 of 24, MRES 32
I-fast) ends at t=1147.06 s, but the capture keeps recording until
t=1302.3 s. That last **153.65 s (153,651 samples)** has no commanded motion
in it at all — nothing is supposed to move. That makes it the cleanest
available segment in this dataset for asking "how much does the measured
position vary when nothing is happening?" Every other quiet stretch in the
capture (between two commanded moves) inherits some settling from the move
before it; this one doesn't.

```
v2/
├── README.md                                this file
├── noise_summary.json                       every metric below, machine-readable
├── scripts/
│   ├── analyze_post_recording_tail.py       tail detection, noise metrics, cascade + interactive
│   └── filter_overlay_tail.py               the three candidate smoothers
├── rendered_assets/
│   ├── tail_zoom_cascade.png                full tail + 3 progressively zoomed panels
│   ├── tail_interactive.html                same segment, full resolution, self-directed zoom
│   └── filter_overlay_tail.png              raw + all three candidate smoothers
├── v2.2/                                    first 60 s only: polynomial detrend + Welch PSD
└── v2.3/                                    is the ~360 Hz peak real? spectra of moving segments
```

## [`v2.2/`](v2.2/README.md) — the first 60 s, detrended

A follow-on that narrows to the **first 60 s** of this tail (excluding the
opening transient below) and subtracts a low-order polynomial before measuring
anything, so the residual is closer to the stationary process that a standard
deviation, a histogram and a PSD all assume. It reports σ = 14.21 nm and
peak-to-peak = 102.4 nm after an order-3 detrend, and a Welch spectrum
(Hann, 16.4 s segments, 50% overlap) that localizes the wander below 0.08 Hz,
finds a white floor of 0.085 nm/√Hz between 10 and 250 Hz, and independently
reproduces the resonance this README locates at 362.8 Hz — at 357.9 Hz, on a
different segment. See [`v2.2/README.md`](v2.2/README.md).

## [`v2.3/`](v2.3/README.md) — is that ~360 Hz peak actually a mode?

Both this README and v2.2 measure the peak in the **same state**: at rest,
nothing commanded. v2.3 applies the same Welch treatment to the states this
segment cannot reach — the oscillation dwells mid-campaign, constant-velocity
stretches of the trajectory legs, and the free decay right after a commanded
microstep. The ring-down is decisive: a free decay rings at the system's own
damped natural frequency, and it lands at 343.8–350.0 Hz across all four MRES
values. So the peak **is a real mechanical mode** (ζ ≈ 0.022–0.029, Q ≈ 17–23,
decaying in 15–21 ms) rather than an artifact of the idle state — but it is
**not dominant**: while moving, commanded full-step-rate harmonics and aliased
microstep-rate lines sit three to four orders of magnitude above it. The
aliasing caveat above is unchanged, and v2.3 narrows the experiment that would
settle it to a few seconds of data. See [`v2.3/README.md`](v2.3/README.md).

```powershell
cd "Microstepping Test Data/v6/Curve fitting/v2/scripts"
python analyze_post_recording_tail.py
python filter_overlay_tail.py
```

## Finding the segment

The boundary is found the same way `ramp_extent`/`still_mask` are already
used elsewhere in v6: scan backward from the end of the capture for the
earliest sample after which the axis is continuously "still" (windowed slope
below threshold, `ids_common.still_mask`) all the way to the last sample —
not a fixed offset after the last leg, an actual detected boundary.

| | |
|---|---|
| Last commanded leg ends | t = 1147.064 s |
| Tail starts | t = 1148.667 s |
| Tail ends | t = 1302.317 s (end of capture) |
| Duration | 153.65 s, 153,651 samples |

## Visualizing it: full segment, progressively zoomed

`tail_zoom_cascade.png` is 4 stacked panels, each ~10-15× narrower than the
one above it: **full tail (153.65 s) → 10 s → 1 s → 0.1 s**, all centered on
the same instant (the temporal midpoint of the tail). The shaded band on each
panel marks exactly where the next narrower one zooms into.

`tail_interactive.html` is the same segment at **full sample resolution — no
decimation** (153,651 points renders smoothly in WebGL directly, so nothing
is thrown away the way the full-capture interactive viewer in `../../analysis/`
has to for its 1.3M-sample capture). Open it in a browser: drag to box-zoom
both axes, scroll to zoom, double-click to reset. This is the "open it and
zoom in yourself" view — the cascade's window choices aren't the only way to
look at it.

**Each cascade panel is normalized independently, to its own start level** —
position is shown as `y[n] - median(y[:25])`, so `y=0` always means "the level
this particular panel starts at," not a value shared across panels (marked
with a solid grey line). The median over the leading 25 samples matters
rather than the literal first sample: the full tail's first sample lands
*inside* the ~4-sample boundary transient described below and is itself the
outlier, so zeroing to it literally would put the whole settled trace ~250 nm
off zero. The y-axis on each panel is then centered on, and tightly fit to,
that panel's own data (using a 1st-99th percentile spread rather than the
true min/max, so that transient can't flatten the rest of the trace) —
enlarging real variation instead of showing it against a loose, shared scale.
**Each panel additionally marks its own local mean** (dotted line, value in
the title) since, because of the slow wander described below, the true local
center genuinely shifts from window to window. `tail_interactive.html` keeps
the whole-tail convention instead (normalized to the tail's single entry
point, one mean line) since it has no discrete "windows" to zero separately —
it's one continuous, self-directed zoom.

## Candidate smoothers on the same segment

`filter_overlay_tail.png` (from `scripts/filter_overlay_tail.py`) applies three
smoothers directly to the same raw tail samples — each computed independently
from the raw signal, so they're directly comparable rather than chained — and
overlays all three on the raw trace:

| Method | Parameters |
|---|---|
| Median despike | width 4 samples |
| Despike + centered Gaussian | window 10 samples, σ = W/6 ≈ 1.67 samples |
| Cubic B-spline (LSQ) | interior knots every ~750 samples |

It uses the **same four-window cascade as `tail_zoom_cascade.png`**
(full tail → 10 s → 1 s → 0.1 s, each panel zeroed to its own start level, the
shaded band marking where the next one zooms into), since at full-tail scale
the three methods are nearly indistinguishable and only separate once you
zoom in. Residual spread after each
(whole-tail std): raw 37.47 nm → despike 37.08 nm → despike+Gaussian
37.00 nm → B-spline 36.57 nm. The std barely moves because it is dominated by
the slow wander (below), which none of these remove — what they do remove is
peak-to-peak: 1271 nm → 276 → 258 → 176 nm. The per-window std in each panel's
legend is where the filtering actually shows: in the 0.1 s window the raw 6.2 nm
drops to 1.7 nm under despike+Gaussian.

### Where the Gaussian's parameters come from

They are measured, not assumed. `temp/psd_resonance_probe.py` detrends a steady
10 s block of the tail (mean + linear, so the thermal wander can't leak into the
high-frequency end), runs a Welch PSD, and locates the peak:

- **fc ≈ 363 Hz**, ~100× above the noise floor, and stable at **362.8 ± 3.0 Hz**
  across 9 of 10 independent 10 s blocks spanning the whole tail.
- The band originally assumed to hold the resonance (120–160 Hz) is flat floor:
  **0.25 nm²** against **19.11 nm²** in 300–400 Hz — roughly 75× less power.

Kernel geometry follows from that, with σ = W/6 so ±3σ spans the window and the
bell tapers to ~zero at its bounds. Two cycles of fc would give W ≈ 5.5, but at
**2.74 samples/cycle** the ripple is too marginally sampled for phase
cancellation to work — a 5-sample kernel leaves visible waviness. Widening to
**W = 10** moves the −3 dB point to **112 Hz**, which rejects the band outright
rather than relying on phase:

| W | −3 dB | gain at fc | 300–400 Hz power |
|---|---:|---:|---:|
| 5 | 225 Hz | 0.161 | 19.11 → 0.67 nm² |
| **10** | **112 Hz** | **0.001** | **19.11 → 0.000 nm²** |

The macro curve is untouched: the wander lives below 5 Hz, two decades under the
cutoff, so there is no corner-cutting.

**Caveat — aliasing.** Nyquist is 500 Hz and the PSD rises steeply into it, the
usual signature of folding. A measured 365 Hz is indistinguishable from 635 Hz
(or 1365 Hz…) at this sample rate. 365 Hz is the right number for designing a
filter that runs *on this 1 kHz series*; it is not evidence of a 365 Hz
mechanical mode. Confirming the true frequency needs a faster capture.

## Noise metrics

Three different questions, three different metrics — a single number here
would be misleading, because the tail contains three genuinely different
phenomena at three different timescales (see "what's actually in the tail"
below).

| Metric | Value | What it answers |
|---|---:|---|
| **Mean \|Δsample\|** (step-to-step) | **6.85 nm** | The literal "average variation in position" — mean absolute change from one 1 kHz sample to the next. Insensitive to slow drift by construction. |
| Median \|Δsample\| | 5.00 nm | Same, robust to the one large outlier (below). |
| Mean std over 153×1 s windows | 8.75 nm | Short-timescale noise, drift-free within each window. |
| **Minimum** std over those windows | **4.13 nm** | The quietest single second found — the closest analogue to the "encoder noise floor" convention used elsewhere in this repo (a comparable ~3 nm was found there over a 1 s window on different data). |
| Maximum std over those windows | 39.5 nm | The noisiest single second — see the outlier note below; this is the window containing it. |
| Whole-tail std (all 153.65 s) | 37.5 nm | Not a noise number — dominated by slow wander over the full duration (below), not by sample-to-sample noise. |
| Whole-tail peak-to-peak | 1.271 µm | Not a noise number either — dominated almost entirely by one transient (below). |

**If asked for a single number: ~5-9 nm** (median/mean \|Δsample\|, and the
typical 1 s-window std) is the steady-state noise-floor indicator. The
37.5 nm whole-tail std is real but answers a different question — see below.

## What's actually in the tail

Three distinct things, at three different timescales, are combined in the
raw numbers above:

1. **A single ~1.15 µm transient in the first 4 ms of the tail.**
   `large_outliers_gt_300nm` in `noise_summary.json` lists 4 consecutive
   step-to-step jumps >300 nm, all within t=1148.667-1148.671 s — the very
   first few samples after the tail boundary. This one event accounts for
   nearly all of the whole-tail peak-to-peak range (excluding the first 1 s
   drops ptp from 1.271 µm to 0.232 µm, a 5.5× reduction) but barely moves
   the mean \|Δsample\| (6.85 → 6.81 nm) or the whole-tail std (37.5 → 37.4 nm)
   — one sample among 150k+ can't move either average much. Read as a brief
   mechanical/electrical settling artifact right at the moment motion
   classification flips to "still," not as ongoing noise.

2. **Slow wander over the full 153.65 s**, visible in the top panel of
   `tail_zoom_cascade.png` as a loose ±0.1-0.2 µm meander with what looks
   like a slow, roughly periodic component (a broad rise-and-fall on a scale
   of tens of seconds — consistent with thermal cycling, though this dataset
   alone can't confirm that). A linear fit finds a slope of
   **+0.49 nm/s (+75 nm total over the tail)** — real but an order of
   magnitude too small to explain the full ±0.1-0.2 µm swing by itself, so
   the wander isn't purely linear drift. This is what the whole-tail std
   (37.5 nm) is actually measuring, and it is **not** the same thing as
   sample noise.

3. **A steady ~5-9 nm sample-to-sample noise floor**, visible once the view
   is zoomed enough (the 1 s and 0.1 s panels) that the slow wander looks
   flat and only the fast component remains. This holds consistently at
   both the 1 s and 0.1 s zoom levels (6.1 nm and 6.2 nm std respectively in
   the example window shown) and across the 153 one-second windows generally
   (4-9 nm typical, excluding the outlier window). **This is the number to
   use as a noise-contribution indicator.**

## Interpreting `noise_summary.json`

```
whole_tail                    -- raw statistics over all 153,651 samples (mixes 1+2+3 above)
step_to_step_variation        -- mean/median/max |y[n+1]-y[n]|  (answers "3" cleanly; max shows "1")
local_1s_windows              -- std per 1 s window, min/mean/max across 153 windows (answers "3", drift-free)
linear_drift                  -- slope + total drift implied by a linear fit (partial answer to "2")
large_outliers_gt_300nm       -- every step-to-step jump >300 nm, with its timestamp (identifies "1")
steady_state_excl_first_1s    -- whole-tail stats recomputed skipping the first 1 s (checks whether "1" was driving "2" -- it isn't)
```

# v6 curve fitting — how it works

Curve-fits `../EXPrun.csv` (the full v6 campaign capture), grouped into one
folder per MRES value. Every source block — the oscillation, each of the 6
trajectory legs — is fitted independently, so a transition in one block can
never bleed into the fit of its neighbour.

**Every block gets two fitted curves, not one.** The original design (still
present, labelled "reference") followed
`../../v3/temp_smoothed_curve_fit/fit_mres4_100pct.py`: discrete, step-like
motion gets a piecewise-constant fit (find each transition, take a robust
median of the settled plateau between edges, with a guard trimmed around
each edge specifically to exclude the transition's own overshoot); continuous
ramp motion gets a 9 ms despike + 20 ms Gaussian smoother heavy enough to
average the individual commanded microsteps out of the ramp entirely. Both
were built to recover the *idealized commanded shape* — flat plateaus, a
smooth ramp — and both do that well.

That turned out to be the wrong target. The actual use for this data is
fitting a simulation model's *dynamic response* — and a model gets tuned
against overshoot and discrete stepping, not against a curve that has already
had both filtered out. So every block now also gets a **"dynamic" fit**: a
much lighter despike-plus-Gaussian pass (see below for how light, and how
that's chosen) that removes sensor noise only, leaving the real transient
shape — overshoot after a commanded step, the individual microsteps riding on
a ramp — intact. That's the curve to fit a model against; "reference" is kept
alongside it only because it's still the cleanest read of the *commanded*
step size and mean ramp trend.

```
Curve fitting/
├── scripts/
│   └── fit_v6_curves.py   run this; writes into the four folders below
├── MRES_1/
├── MRES_4/
├── MRES_16/
└── MRES_32/
```

Run from `scripts/`:

```powershell
python fit_v6_curves.py
```

Each `MRES_*/` folder gets:

| File | What it is |
|---|---|
| `oscillation_fit_overlay.png` | the full 15-cycle, ~30 s block: measured vs. both fits |
| `oscillation_microstep_zoom.png` | the first 6 cycles, zoomed, with the seed edge times marked |
| `trajectory_fits_overlay.png` | all 6 trajectory legs (D/I × slow/mod/fast): measured vs. both fits, full extent |
| `trajectory_fits_apex_zoom.png` (+ `.svg`) | same 6 legs, zoomed on the direction reversal at the top of the ramp |
| `trajectory_fits_midramp_zoom.png` (+ `.svg`) | same 6 legs, zoomed on a plain constant-velocity stretch |
| `fitted_blocks.npz` | every block's raw arrays: time, measured, **reference**, **dynamic** |
| `fit_summary.json` | method + parameters (both fits) and file manifest for that MRES |

In every plot, dashed grey is "reference" and solid teal is "dynamic" — see
above for what each means. The two trajectory zoom PNGs also get an `.svg`
twin: vector, so it never pixelates no matter how far it's zoomed into after
the fact (the PNGs are still 300 DPI at 18×10in for anyone who'd rather not
open a vector viewer).

At full extent (`trajectory_fits_overlay.png` / `oscillation_fit_overlay.png`),
overshoot and individual steps are real but small relative to the whole
block, so both fits can still look close together there. The zoom views make
the difference obvious:

- **Apex zoom** — centered on the direction reversal, where "reference"
  (the heavy filter) is most likely to show lag or over-rounding that
  "dynamic" doesn't.
- **Mid-ramp zoom** — a window 40% of the way from ramp start to apex, clear
  of both the start transient and the turnaround, sized to a target number of
  *commanded steps* rather than a fixed duration (see below), at
  individual-sample resolution (plotted as dots, not just a line).

The mid-ramp zoom also makes the MRES difference visible in a way the
oscillation-block fits don't on their own: at **MRES 1** (10 µm full steps)
the raw trace is a visible staircase that "dynamic" tracks and "reference"
smooths away (see `MRES_1/trajectory_fits_midramp_zoom.png`); by **MRES 32**
individual 0.3 µm-nominal steps are far too fine to resolve at
trajectory-leg position scale even for "dynamic," and the raw trace already
looks like a smooth line before any fitting happens
(`MRES_32/trajectory_fits_midramp_zoom.png`) — consistent with what the
oscillation-block fits in the same folder already show about step fidelity
collapsing at fine MRES.

### Sizing the zoom window to commanded steps, not a fixed duration

A fixed time window is what made this indistinguishable at fine MRES in an
earlier version: at MRES 32 a commanded step is 32× smaller than at MRES 1,
so the same 1 s window that clearly shows ~4 steps at MRES 1 packs in ~128 at
MRES 32 — far too dense to resolve at any zoom. The mid-ramp window is sized
instead to `TARGET_STEPS_EACH_SIDE = 10` commanded steps each side of center
(`step_um / velocity_um_s`, where velocity is that leg's own commanded rate —
`RATES_FSPS[...] × FULL_STEP_UM`, independent of MRES since MRES only
subdivides a full step, it doesn't change the commanded feed rate), floored
at `MIN_SAMPLES_EACH_SIDE = 40` samples so a fast-rate + fine-MRES leg (whose
steps land under 1 ms apart, i.e. faster than the 1 kHz sample period) still
gets a window with something in it rather than collapsing to 2-3 samples.

## Why v6 needed a different seeding strategy than v3

v3's script fits one run whose host controller logged a `MOVE_ACK` timestamp
for every commanded microstep (`identification_controller_log.csv`). Those
timestamps are used to **seed** each expected transition; the fit then
**refines** each seed to the nearest real edge in the encoder trace, and
takes a **median** of the settled plateau between edges. Three steps, in that
order: seed, refine, estimate.

`EXPrun.csv` is a bare EL5101 encoder trace with no accompanying event log
(see `../analysis/scripts/plot_v6_exp_run.py`'s docstring — this is also why
that script has to infer trajectory-leg and oscillation-block boundaries
straight from the signal). There are no `MOVE_ACK` timestamps to seed from.
The refine and estimate steps carry over unchanged; only the seed step had to
be rebuilt from what the signal itself contains.

### Seeding the oscillation's 30 transitions (15 cycles, forward + return)

Two paths, chosen per MRES group at runtime and recorded in that group's
`fit_summary.json` under `oscillation.seed_method`:

**1. Detected directly (used for MRES 1 and 4 here).** Every marker (the
small out-and-back move that separates blocks, since there's no event log to
read boundaries from instead) shows up as a cluster of large single-sample
jumps lasting several hundred ms or more — the move itself takes that long at
`MARKER_RATE_MILLIHZ`. A single oscillation microstep, by contrast, is one
pulse at the much faster `OSCILLATION_RATE_MILLIHZ` and clusters to under
~80 ms. Where a run of ≥28 such short events, spaced 0.8-1.3 s apart (the
firmware's 1 s forward/return dwell), is found before the group's first
trajectory leg, those detected times are used **as** the seed times directly
— the closest available analogue of a logged `MOVE_ACK`.

**2. Derived from marker timing (used for MRES 16 and 32 here).** At finer
MRES the commanded step is small enough that its individual transitions don't
reliably clear a noise-safe jump threshold — the same limitation
`v5/analysis/scripts/plot_v5_microstepping.py` documents for its own MRES 32
block. Instead:

- Find the marker events in two windows, kept **separate on purpose**: a
  narrow 25 s window ending at the leg's own ramp-threshold crossing (for the
  marker immediately before this group's first leg — always close, observed
  13-20 s out) and a wider 60 s window further back (for the CONFIG and
  OSCILLATION markers that open the group). A single wide window was tried
  first and is fragile: at the slowest commanded trajectory rate a leg's own
  ramp moves as little as 0.275 µm/ms, mostly *below* the jump threshold, so
  it fragments into many small, irregularly-spaced events that can land at
  the same array positions the real markers would otherwise occupy and throw
  off fixed-offset indexing into the event list.
- `osc0` (the first forward transition) is the OSCILLATION marker's return
  move **finishing** (its detected start time, plus that event's own
  duration — a marker's return isn't instantaneous, it's a ~0.3-0.55 s ramp —
  plus the firmware's 500 ms settle) rather than when the return move merely
  *starts*. Missing the event's own duration was caught by checking one case
  (MRES 16) against the raw signal: without it, `osc0` landed while the axis
  was still actively ramping from −992 µm to −195 µm, and the resulting
  "oscillation" fit was one giant ~800 µm step — obviously not a MRES 16
  oscillation (commanded step 0.625 µm). See `oscillation_fit_overlay.png`
  for what the corrected version looks like instead.
- The period is `(next marker's onset − osc0) / 15`. Measured this way it
  comes out to 2.0035 s (MRES 16) and 2.0038 s (MRES 32) — matching the
  *directly detected* periods from MRES 1 and 4 (~2.006-2.012 s) to within
  0.3%, which is the cross-check that the derived path is finding the same
  real timing the detected path finds directly.

Either way, every seed time is then refined identically (see below) — the
seed only has to be close enough for the refinement step to find the real
edge, not exact.

### The "reference" fit: refining edges and estimating plateaus (`fit_steps`)

Unchanged from v3, applied to whichever seed times came out of the step
above:

1. Despike the block with a centered 9 ms median filter.
2. For each seed time, search ±0.35 s (wider than v3's ±0.12 s — v6's seeds
   are measured/derived, not logged, so allow more slack) and pick the point
   where `|despiked[i + span] − despiked[i − span]|` is largest (`span` = 6
   ms) — the point of steepest local change, i.e. the real transition.
3. Between consecutive refined edges, take the **median** of the despiked
   signal with a 0.2 s guard trimmed from each side (so the transition's own
   settling ring-down never enters a plateau's estimate). That median is the
   fitted value for the whole plateau — a genuine flat segment, not a rounded
   corner.

This is the fit that deliberately throws the overshoot away (that's what the
0.2 s guard is *for* — keeping the ring-down out of the plateau estimate) and
why it's "reference," not the primary curve, now.

### The "dynamic" fit: light enough to keep overshoot and individual steps

Same two-stage filter shape (median despike, then Gaussian) as v3's
continuous-block filter, just applied much more lightly, and applied to the
*whole* block rather than per-plateau — nothing gets flattened to a constant.

**Oscillation blocks** use a fixed light pass: 3 ms despike, 2 ms Gaussian
sigma (`OSC_DYNAMIC_MEDIAN_MS` / `OSC_DYNAMIC_SIGMA_MS`). Fixed because every
oscillation transition is well isolated (the firmware's 1 s forward/return
dwell puts them nowhere near each other), so there's no "how many steps are
crowded into this window" problem the way there is on a ramp.

**Trajectory legs** scale both parameters to that specific leg's own
commanded step interval (`ramp_dynamic_params`), because the interval spans
three orders of magnitude across this dataset — from 36 ms (MRES 1, slow) to
0.156 ms (MRES 32, fast), the latter already faster than the 1 ms sample
period. A single fixed light filter would either still blur MRES 1's steps or
do essentially nothing at MRES 32:

```
median_ms = clip(0.25 × step_interval_ms, 1.0, 9.0)
sigma_ms  = clip(0.15 × step_interval_ms, 0.5, 5.0)
```

At MRES 1 slow (36 ms interval) that's 9.0/5.0 ms — light relative to the
interval, and `MRES_1/trajectory_fits_midramp_zoom.png` shows the resulting
staircase clearly. At MRES 32 fast (0.14 ms interval) it floors out at
1.0/0.5 ms — essentially just despiking, because there is no step structure
left to preserve once steps arrive faster than samples do; the honest result
is a smooth line, not a fabricated staircase. Every leg's actual
`(median_ms, sigma_ms)` is in that panel's title and in `fit_summary.json`
under `trajectory_legs.dynamic.params_by_leg_ms`.

### What's genuinely different from v3, and what isn't

| | v3 (`fit_mres4_100pct.py`) | v6 "reference" | v6 "dynamic" |
|---|---|---|---|
| Step seed source | logged `MOVE_ACK` timestamps | detected transitions, or derived from marker timing | same seeds, used only for the apex-zoom marker line |
| Edge refinement | max local slope, ±0.12 s | max local slope, ±0.35 s (wider — seeds are measured/derived, not logged) | n/a (no per-segment edges; one continuous light filter) |
| Step/oscillation filter | median plateau, 0.04 s guard (overshoot excluded) | median plateau, 0.2 s guard (overshoot excluded) | despike 3 ms + Gaussian 2 ms, whole block (**overshoot kept**) |
| Ramp/continuous filter | 9 ms median + 20 ms Gaussian | same (idealized trend) | despike/Gaussian scaled per leg to ~0.25×/0.15× its own step interval, clipped to [1,9]/[0.5,5] ms (**steps kept**) |
| Per-block independence | yes | yes | yes |
| Grouping | one run, all blocks in one montage | one folder per MRES value | same folder, second curve in every plot |

v3 only ever produced the "reference" kind of fit — an idealized command
shape. v6's "reference" column carries that forward with a different seed
source (since there's no event log — see above), same filter parameters,
same purpose. "Dynamic" is new: not a v3 concept at all, added because
fitting a simulation model needs the actual transient response, not the
idealized one.

The seed-timing work above (detecting or deriving `osc0`/period) went wrong
twice during development before matching the directly-detected MRES 1/4
periods to within 0.3% — see the marker-duration and decoupled-window notes
in that section and in `fit_v6_curves.py`'s own comments. It still matters
for "dynamic": the oscillation block's extraction window is bounded by the
same seed times, and with only a light filter left to smooth it, that window
running even slightly into the next marker's own fast excursion is no longer
invisible the way it was under the heavy reference filter — it showed up
directly as a blown-out y-axis (±100 µm on a ~10 µm MRES 1 oscillation)
during development, fixed by trimming the window to end within one dwell
period past the last real transition instead of 1.5 dwell periods past it.

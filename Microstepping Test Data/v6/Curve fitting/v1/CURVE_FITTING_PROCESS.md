# v6 curve fitting — mathematical procedure

This documents the actual algorithm `scripts/fit_v6_curves.py` runs: the
filtering operators, the edge-detection scoring function, the plateau
estimator, and the formulas that derive every seed time and every adaptive
filter parameter from the signal itself. For *why* this design was chosen —
two fits instead of one, why v6's seeding had to be rebuilt from scratch, the
bugs found while building it — see `README.md` in this folder. This document
is the "what, exactly, gets computed"; that one is the "why."

## 1. Signal model and notation

The raw capture is a sequence of position samples

```
y[n],  n = 0 .. N-1,   sampled at  dt = 1/1000 s   (1 kHz)
t[n] = n · dt
```

Every quantity below is either a **filter** applied to (a slice of) `y`, an
**index** derived from `y` by thresholding or scoring, or a **parameter**
computed from the commanded MRES/rate rather than from `y` at all. The
pipeline composes these in a fixed order per block; nothing here is fit in
the least-squares sense — every "fit" is a deterministic filter or estimator
applied to the measured signal.

## 2. Core filtering operators

Two operators, and one composition of them, are reused throughout.

### 2.1 Median despike, M_w

Given a window duration `w` (ms), convert to a sample count and force it odd:

```
k  = round(w · 1e-3 / dt)
k' = k + (1 - k mod 2)          # smallest odd integer >= k
r  = (k' - 1) / 2
```

Then

```
M_w[y](n) = median( y[n-r], y[n-r+1], ..., y[n+r] )
```

with the sequence edge-replicated ("nearest") outside `[0, N-1]`. This is a
robust despike: single-sample sensor glitches are suppressed without
attenuating real transitions that span several samples.

### 2.2 Centered Gaussian smoother, G_σ

Given a sigma `σ` (ms), convert to samples and set a truncation radius:

```
σ_n = σ · 1e-3 / dt
R   = max(1, ceil(4 · σ_n))
g[k] = exp(-k² / (2σ_n²)) / Σ_{j=-R}^{R} exp(-j² / (2σ_n²)),   k = -R .. R
```

```
G_σ[y](n) = Σ_{k=-R}^{R}  y[n-k] · g[k]
```

`g` is symmetric, so this is a **zero-phase** filter — no time lag is
introduced, unlike a causal (one-sided) smoother. Edges are handled the same
"nearest" way as §2.1.

### 2.3 Combined smoothing pass, S_{w,σ}

```
S_{w,σ}[y] = G_σ( M_w[y] )
```

despike first, then smooth. Every "fit" in this pipeline that isn't a
piecewise-constant plateau is one evaluation of `S_{w,σ}` for some `(w, σ)` —
the entire design space is *which* `(w, σ)` to use, and that choice is where
"reference" and "dynamic" diverge (§4).

## 3. Seeding: recovering transition times with no event log

v3's equivalent pipeline seeds each expected step transition from a logged
`MOVE_ACK` timestamp. `EXPrun.csv` has no such log, so 30 seed times
`τ_0 .. τ_29` (15 forward-edges, 15 return-edges) are recovered from `y`
itself before any fitting happens.

### 3.1 Jump-event detector

```
d[n] = | y[n+1] - y[n] |
E = { n : d[n] > θ },   θ = 0.8 µm/sample
```

Consecutive members of `E` are grouped into **events**: `n_i, n_{i+1} ∈ E`
join the same event iff `t[n_{i+1}] - t[n_i] < 0.15 s`. Each event `j` has a
start time `t_j` and a duration `Δ_j = t[last] - t[first]` within the group.
A marker (fast, large out-and-back move) produces one long event
(`Δ_j` from ~80 ms to ~1.5 s, scaling with its own amplitude); a single
oscillation microstep produces one short event (`Δ_j < 80 ms`).

### 3.2 Case A — direct detection (used for MRES 1, 4 in this dataset)

Let `M = { j : Δ_j ≤ 80 ms }` (candidate microstep events). If a contiguous
run of `M`, sorted by time, has consecutive gaps `t_{j+1} - t_j ∈ [0.8, 1.3] s`
(the firmware's 1 s forward/return dwell) and length ≥ 28, take its first 30
members' start times directly:

```
τ_i = t_{run[i]},   i = 0 .. 29
```

### 3.3 Case B — derived from marker timing (used for MRES 16, 32)

At finer MRES the commanded step rarely clears `θ`, so §3.2 fails. Instead:

```
next_onset  = t_j  for the marker event j immediately preceding this
              group's first trajectory leg (searched in a 25 s window
              ending at that leg's own ramp-threshold crossing)

osc_return  = (t_j, Δ_j)  for the OSCILLATION marker's return-move event
              (searched separately, further back, among marker-duration
              events strictly before next_onset - 5 s)

osc0   = t_j + Δ_j + 0.5 s                    # move finishes, then settle
period = ( next_onset - osc0 ) / 15
```

```
τ_{2k}   = osc0 + k · period,           k = 0 .. 14      (forward edges)
τ_{2k+1} = osc0 + (k + 0.5) · period,   k = 0 .. 14      (return edges)
```

`t_j + Δ_j` (the move **finishing**, not starting) matters: using `t_j` alone
put `osc0` mid-ramp on the marker's own return move in one checked case
(MRES 16), and the derived period cross-checks the directly-detected periods
from §3.2 to within 0.3% only once the duration term is included.

## 4. The reference fit: piecewise-constant plateaus

Applied to the oscillation block, given seed times `τ_0..τ_29` from §3 and a
search half-width `h = 0.35 s`:

**Edge refinement.** Despike first: `ŷ = M_{9ms}[y]`. Let
`s = round(0.006 / dt)` samples (a 6 ms half-span). For each seed `τ_i`,
search the window `n ∈ [ idx(τ_i - h), idx(τ_i + h) ]` and score every
candidate:

```
score(n) = | ŷ[n+s] - ŷ[n-s] |
edge_i   = argmax_n  score(n)
```

`score` is a discrete central-difference slope estimate over a fixed 12 ms
baseline — maximizing it finds the steepest local change, i.e. the real
transition, even if the seed itself is off by up to `h`. (If `edge_i ≤ edge_{i-1}`,
it's bumped to `edge_{i-1} + 1` to keep edges strictly increasing.)

**Plateau estimate.** For each interval `(edge_i, edge_{i+1})` of width
`W = edge_{i+1} - edge_i` samples, with guard
`g = min( round(0.2/dt), ⌊W/4⌋ )`:

```
fitted[n] = median( ŷ[edge_i + g] .. ŷ[edge_{i+1} - g] )   for all n in the interval
```

trimming `g` samples off both sides of the interval before taking the
median excludes the transition's own overshoot/ring-down from the plateau
estimate by construction — this is a piecewise-constant, idealized-command
curve, and the overshoot exclusion is deliberate, not incidental.

Trajectory legs get the continuous analogue with no edges at all:
`reference = S_{9ms, 20ms}[y]` over the whole ramp — a single heavy smoothing
pass, no seeding involved.

## 5. The dynamic fit: transient-preserving

No edges, no plateaus — a single evaluation of `S_{w,σ}` (§2.3) over the
*whole* block, with `(w, σ)` chosen to be light relative to whatever
structure needs to survive.

**Oscillation blocks** use a fixed pair, since every transition is isolated
in time by a 1 s dwell on each side:

```
dynamic = S_{3ms, 2ms}[y]
```

**Trajectory legs** scale `(w, σ)` to that specific leg's own commanded step
interval, because the interval is not constant across this dataset:

```
Δx   = FULL_STEP_UM / mres                        # commanded step size, µm
v    = rate_fsps · FULL_STEP_UM                   # commanded velocity, µm/s
                                                   # (independent of MRES)
τ_step = 1000 · Δx / v                             # step interval, ms

w = clip( 0.25 · τ_step, 1.0, 9.0 )                # despike width, ms
σ = clip( 0.15 · τ_step, 0.5, 5.0 )                # Gaussian sigma, ms

dynamic = S_{w, σ}[y]
```

`τ_step` spans three orders of magnitude across the four MRES values and
three rates in this dataset — 36 ms (MRES 1, 27.5 full-steps/s) down to
0.156 ms (MRES 32, 200 full-steps/s), the latter already below the 1 ms
sample period — which is why `(w, σ)` is recomputed per leg rather than
fixed. At the floor (`w=1.0, σ=0.5`), the filter is doing essentially nothing
beyond despiking, and correctly renders a smooth line rather than a
fabricated staircase: there is no recoverable step structure once the
commanded interval is shorter than the sample period.

Every leg's resolved `(w, σ)` is written to `fit_summary.json` under
`trajectory_legs.dynamic.params_by_leg_ms`.

## 6. Zoom-window sizing (visualization, not a fit)

The two trajectory zoom plots (apex, mid-ramp) are windowed to a target
number of commanded steps rather than a fixed duration, using the same
`Δx, v` as §5:

```
half_window = max( κ · Δx / v,  n_min · dt ),   κ = 10,  n_min = 40
```

the `n_min · dt` floor exists so a fast-rate + fine-MRES leg (steps under
1 ms apart) still gets a window with a visible number of raw samples in it,
rather than collapsing to 2-3 points. A fixed-duration window (tried first)
is why an earlier version of this plot was indistinguishable at fine MRES:
the same 1 s window that shows ~4 steps at MRES 1 packs in ~128 at MRES 32.

## 7. Setup and execution

```
Curve fitting/
├── CURVE_FITTING_PROCESS.md   this document
├── README.md                  design rationale and debugging history
├── scripts/fit_v6_curves.py   implements §2-§6; run this
└── MRES_1/, MRES_4/, MRES_16/, MRES_32/
```

```powershell
cd "Microstepping Test Data/v6/Curve fitting/scripts"
python fit_v6_curves.py
```

`main()` loads `../../../EXPrun.csv` once, locates all 24 trajectory legs and
the whole-run "still" mask, then runs §3-§6 independently for each of the 4
MRES groups. Outputs per group: `oscillation_fit_overlay.png`,
`oscillation_microstep_zoom.png`, `trajectory_fits_overlay.png`,
`trajectory_fits_apex_zoom.png`/`.svg`, `trajectory_fits_midramp_zoom.png`/`.svg`,
`fitted_blocks.npz` (every block's `measured`/`reference`/`dynamic` arrays —
**fit a simulation model against `dynamic`**), and `fit_summary.json`.

## 8. Known limits

- **Sub-sample-period steps aren't recoverable** (§5): at MRES 32 fast,
  `τ_step < dt`, so there is no staircase in the raw signal for any filter to
  preserve.
- **Case B (§3.3) depends on marker timing being well-formed** — cross-checked
  against Case A to 0.3% agreement, but only 4 MRES groups' worth of
  validation exists in this dataset.
- **Both fits assume the extraction window itself is correct** (`ramp_extent`
  for legs, the seed-time span for the oscillation) — that boundary
  detection lives in `ids_common.py` / `plot_v6_exp_run.py`, documented in
  `../../analysis/README.md`, not in this script.

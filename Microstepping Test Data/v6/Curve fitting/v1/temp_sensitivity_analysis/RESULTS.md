# Dynamic-fit sensitivity — initial testing

**Script**: `sensitivity_dynamic_fit.py` | **Raw output**: `sensitivity_results.json`
**Plots**: `plots/sensitivity_overview.png`, `plots/floor_case_perturbation.png`

Same type of analysis as `v3/Parameter Optimization/calibration_bracketing/
step6_sensitivity_table.py` (cost = RMS(measured - fit), sensitivity =
`0.5*(|cost(+10%)-cost(base)| + |cost(-10%)-cost(base)|)/cost(base)`, 1%
noise floor), adapted to the much cheaper `smooth()` despike+Gaussian pass
in `fit_v6_curves.py` instead of a LuGre time-domain simulation. Runs in
under a second — no ODE solve, no re-parsing the raw capture, everything
comes from each MRES group's already-cached `fitted_blocks.npz` and
`fit_summary.json`.

**Scope** (initial testing, not the full sequence): four segments chosen to
contrast the two ends of the ramp-dynamic clip range plus an oscillation
comparison —

| segment | kind | base filter (median/sigma ms) |
|---|---|---|
| `osc_MRES1` | oscillation | 3.0 / 2.0 (fixed) |
| `osc_MRES32` | oscillation | 3.0 / 2.0 (fixed) |
| `leg_MRES1_Dmod` | ramp | 3.57 / 2.14 (unclipped) |
| `leg_MRES32_Dfast` | ramp | 1.0 / 0.5 (floored at the hard lower bound) |

## Results

| segment | param | base (ms) | base cost (µm) | sensitivity | verdict |
|---|---|---|---|---|---|
| osc_MRES1 | median_ms | 3.00 | 0.3050 | 0.0000 | insensitive |
| osc_MRES1 | sigma_ms | 2.00 | 0.3050 | 0.0378 | sensitive |
| osc_MRES32 | median_ms | 3.00 | 0.0531 | 0.0000 | insensitive |
| osc_MRES32 | sigma_ms | 2.00 | 0.0531 | 0.0023 | insensitive |
| leg_MRES1_Dmod | median_ms | 3.57 | 1.9631 | 0.0317 | sensitive |
| leg_MRES1_Dmod | sigma_ms | 2.14 | 1.9631 | 0.0349 | sensitive |
| leg_MRES32_Dfast | median_ms | 1.00 | 0.0618 | 0.0000 | insensitive |
| leg_MRES32_Dfast | sigma_ms | 0.50 | 0.0618 | **0.3166** | **sensitive** |

## Findings

1. **`median_ms` reads as zero-sensitivity everywhere tested — this is a
   discretization artifact, not evidence the despike window doesn't
   matter.** At 1 kHz, `median_filter`'s window is
   `round(median_ms / dt_ms)`, forced odd. A ±10% nudge around 3.0 ms
   (→ 2.7/3.3) or 1.0 ms (→ 0.9/1.1) rounds right back to the same integer
   window both directions, so the fit is bit-for-bit identical and the cost
   doesn't move at all — confirmed in `plots/floor_case_perturbation.png`
   (left panel is exactly flat at 0). The one case where `median_ms` *did*
   show up as sensitive, `leg_MRES1_Dmod` (base 3.57 ms), is because 3.57 ms
   sits close enough to a rounding/parity boundary that +10% still rounds to
   the same window as the base but −10% rounds down to a different one —
   i.e. even that "sensitivity" number is really measuring which side of an
   integer boundary the perturbation landed on, not a smooth dependence.
   **Takeaway: don't trust a ±10% sweep on `median_ms` at this sample rate
   to characterize its effect continuously — the honest statement is "this
   filter step only responds to changes on the order of ±1 sample
   (±0.5-1 ms)," and a sweep design that actually steps in whole samples
   would be more informative than a fractional-percent one.**

2. **The oscillation dynamic fit is far more sensitive to `sigma_ms` at
   MRES 1 than MRES 32, for a real (not artifactual) reason.** Same fixed
   3.0/2.0 ms filter in both cases, but MRES 1's ~10 µm commanded
   one-microstep oscillation is well above the noise floor (base cost
   0.305 µm), so a change in blur amount changes how much of the real step
   transition survives. MRES 32's commanded step is ~0.3 µm — already at or
   below what this fit can resolve (base cost 0.053 µm is mostly noise
   already) — so blurring it more or less barely changes anything. This
   matches `CURVE_FITTING_PROCESS.md`'s own note that MRES 32 trajectory
   fits are "honestly smooth" with nothing left to preserve; the same
   pattern shows up quantitatively here for the oscillation block.

3. **The one segment that stands out sharply is `leg_MRES32_Dfast`'s
   `sigma_ms`: 0.32 relative sensitivity, ~10x every other cell.** This is
   the ramp-dynamic filter floored at its hard lower bound
   (`RAMP_DYNAMIC_SIGMA_BOUNDS_MS = (0.5, 5.0)`). At `sigma_ms = 0.5 ms`
   (= 0.5 samples at 1 kHz), the Gaussian kernel is in the regime where it
   goes from "does almost nothing" to "measurably smooths" over a very
   small parameter range — exactly where a filter's response is steepest
   per unit change in σ. `plots/floor_case_perturbation.png` (right panel)
   shows the actual fit-vs-fit difference: a ±10% σ change moves the
   "dynamic" curve by up to ~0.09 µm, on the same order as the entire base
   RMS-to-measured cost (0.062 µm) for this segment — i.e. this isn't a
   relative-sensitivity artifact, the absolute effect is real and roughly
   as large as what's being fit in the first place. Two short near-zero
   stretches in that same panel line up with this leg's rest/dwell
   sections, where there's nothing for σ to blur regardless of its value.
   **Takeaway: MRES 32 fast being pinned at the ramp-dynamic filter's floor
   isn't a "settled, stop worrying about it" state — it's the single most
   parameter-sensitive combination tested. If a friction/tracking-error
   model ever gets fit against this leg's `dynamic` curve, the σ floor
   value deserves a deliberate choice, not just "whatever the clip
   produced."**

## Caveats

- `median_ms`'s apparent insensitivity is a sampling-rate artifact (see
  Finding 1) — re-running with perturbations sized in whole samples, not a
  fixed ±10%, would give a more honest picture of that parameter's true
  effect.
- Only 4 of the 4 MRES × 6-leg × (oscillation + legs) = 28 possible blocks
  were tested, chosen to contrast known extremes (bound-floored vs.
  mid-range vs. two MRES values on the same fixed oscillation filter). This
  is explicitly initial testing, not a claim about every block's behavior —
  in particular nothing here says anything about MRES 4 or MRES 16, or any
  of the `I-*` (individual-command) legs.
- Cost is RMS-to-*measured* position, i.e. it rewards a fit that hugs the
  noisy raw signal. That's the right metric for probing the filter's own
  sensitivity (this analysis's goal), but it is **not** the same thing as
  "how good the `dynamic` fit is for its actual intended purpose" (matching
  a simulation model's dynamic response) — a fit that hugs sensor noise
  tightest is not necessarily the best target for that.

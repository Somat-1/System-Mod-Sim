"""Curve-fit the v6 EXPrun capture, grouped by MRES.

Same method as `Microstepping Test Data/v3/temp_smoothed_curve_fit/
fit_mres4_100pct.py`: discrete step-like motion gets a piecewise-constant fit
(seed the expected transition times, refine each edge from the local encoder
slope, take a robust median per plateau); continuous ramp motion gets a
median despike plus centered Gaussian smoother. Each source block is still
fitted independently so no transition bleeds into a neighbour.

v3's version seeds its step transitions from `MOVE_ACK` timestamps in a host
controller's event log. v6 has no such log (`EXPrun.csv` is a bare EL5101
capture -- see plot_v6_exp_run.py's docstring), so the discrete blocks here
(each MRES's 15-cycle one-microstep oscillation) instead seed from timing
recovered directly out of the signal:

  1. Every marker (the small out-and-back move used to segment blocks, since
     there is no event log) shows up as a cluster of large single-sample
     jumps lasting 100-800+ ms -- the move itself takes that long at
     MARKER_RATE_MILLIHZ. A single oscillation microstep, by contrast, is one
     pulse at OSCILLATION_RATE_MILLIHZ and clusters to under ~80 ms.
  2. Where the oscillation's own 30 individual transitions clear the jump
     threshold (MRES 1 and 4 here), they are used AS the seed times directly
     -- the closest available analogue of MOVE_ACK.
  3. Where they do not (MRES 16 and 32 -- the commanded step is too small to
     reliably exceed a noise-safe jump threshold), the two markers preceding
     the oscillation (CONFIG, then OSCILLATION) and the one following it
     (before the group's first trajectory leg) are used instead: osc0 is the
     OSCILLATION marker's return plus its 500 ms settle, and the period is
     (next marker's onset - osc0) / 15.
Either way, every seed time is then refined the same way v3's are: snapped to
the nearest local maximum of |slope| within a small search radius, so small
seed error does not survive into the fit.

Trajectory legs are the continuous case: reuses `leg_rest_levels`'s
`ramp_extent` window (already the genuine full extent of the ramp, settled on
both sides).

REVISED: v3's filters (0.2 s guard-trimmed median plateaus for steps, 20 ms
Gaussian for ramps) are both tuned to produce an *idealized command* shape --
flat plateaus, a smooth ramp -- by design excluding exactly the things a
tracking-error/friction model needs to be fit against: the overshoot right
after each commanded step, and the discrete micro-step structure riding on
top of each ramp (visible once `trajectory_fits_midramp_zoom.png` is zoomed
enough -- see that plot's own git history). Once the actual goal became
matching a simulation model's dynamic response to this data, an idealized
target is the wrong target.

Every block therefore now gets TWO fitted curves, not one:
  - "reference": the original idealized fit (flat plateau / heavy Gaussian),
    kept because it is still the cleanest read of the *commanded* step size
    and the *mean* ramp trend.
  - "dynamic": a much lighter despike-plus-Gaussian pass (ramp: adaptive to
    that leg's own commanded step interval; oscillation: a fixed light pass
    sized to preserve a transition's overshoot/ring-down) that removes only
    sensor noise and keeps the real transient shape -- overshoot, individual
    steps -- intact. This is the one to fit a model against.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, median_filter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "analysis" / "scripts"))
from ids_common import load_ids, still_mask, ramp_extent, FULL_STEP_UM, apply_style as ids_apply_style
from plot_v6_exp_run import CSV, find_trajectories, MRES_VALUES, RATES_FSPS

OUT_ROOT = Path(__file__).resolve().parent.parent
CYCLES = 15
OSCILLATION_HALF_DWELL_S = 1.0
MARKER_SETTLE_S = 0.5
MEDIAN_WINDOW_MS = 9.0
GAUSSIAN_SIGMA_MS = 20.0
RAW_COLOR, FIT_COLOR, REF_COLOR = "#777777", "#007f86", "#9a9a9a"
SEED_COLOR = "#bd6b18"

# "Dynamic" pass: light enough to leave overshoot and individual steps
# intact, unlike the "reference" (idealized) MEDIAN_WINDOW_MS/GAUSSIAN_SIGMA_MS
# pass above, which is tuned to erase both on purpose.
OSC_DYNAMIC_MEDIAN_MS = 3.0
OSC_DYNAMIC_SIGMA_MS = 2.0
RAMP_DYNAMIC_MEDIAN_FRAC = 0.25   # of the commanded step interval
RAMP_DYNAMIC_SIGMA_FRAC = 0.15
RAMP_DYNAMIC_MEDIAN_BOUNDS_MS = (1.0, 9.0)
RAMP_DYNAMIC_SIGMA_BOUNDS_MS = (0.5, 5.0)


def ramp_dynamic_params(step_um, velocity_um_s):
    """(median_ms, sigma_ms) for a ramp's light pass, scaled to that leg's own
    commanded step interval so the filter is light relative to the structure
    it needs to preserve at every MRES/rate combination, not just MRES 1."""
    step_interval_ms = 1e3 * step_um / velocity_um_s
    median_ms = np.clip(RAMP_DYNAMIC_MEDIAN_FRAC * step_interval_ms,
                        *RAMP_DYNAMIC_MEDIAN_BOUNDS_MS)
    sigma_ms = np.clip(RAMP_DYNAMIC_SIGMA_FRAC * step_interval_ms,
                       *RAMP_DYNAMIC_SIGMA_BOUNDS_MS)
    return float(median_ms), float(sigma_ms)

# Trajectory zoom plots: rendered large and at high DPI, plus an SVG
# (vector, so it never pixelates no matter how far the viewer zooms into it)
# alongside every PNG -- the whole point of these plots is to be zoomed into
# further after the fact.
ZOOM_FIGSIZE = (18, 10)
ZOOM_DPI = 300
RATE_BY_LABEL = {"D-slow": RATES_FSPS[0], "D-mod": RATES_FSPS[1], "D-fast": RATES_FSPS[2],
                 "I-slow": RATES_FSPS[0], "I-mod": RATES_FSPS[1], "I-fast": RATES_FSPS[2]}


def savefig_png_and_svg(fig, out_dir, stem, dpi=ZOOM_DPI):
    fig.savefig(out_dir / f"{stem}.png", dpi=dpi)
    fig.savefig(out_dir / f"{stem}.svg")


def apply_style():
    ids_apply_style()


# ---------------------------------------------------------------- events ---

def big_jump_events(t, pos, dt_s, lo, hi, jump_threshold_um, cluster_gap_s=0.15):
    """Cluster consecutive large single-sample jumps within [lo, hi) into
    events. Returns (event_start_idx, event_duration_s) arrays, oldest first."""
    d = np.abs(np.diff(pos[lo:hi]))
    idx = np.flatnonzero(d > jump_threshold_um) + 1 + lo
    if idx.size == 0:
        return np.array([], dtype=int), np.array([])
    events, cur = [], [idx[0]]
    for i in idx[1:]:
        if t[i] - t[cur[-1]] < cluster_gap_s:
            cur.append(i)
        else:
            events.append(cur)
            cur = [i]
    events.append(cur)
    starts = np.array([c[0] for c in events])
    durs = np.array([(c[-1] - c[0]) * dt_s for c in events])
    return starts, durs


def find_oscillation_seed_times(t, pos, dt_s, leg1_start_idx, search_back_s=60.0):
    """Return (seed_times, method) for one MRES group's 30 oscillation
    transitions (15 forward + 15 return, alternating)."""
    lo = max(0, leg1_start_idx - int(search_back_s / dt_s))
    starts, durs = big_jump_events(t, pos, dt_s, lo, leg1_start_idx, jump_threshold_um=0.8)
    if starts.size == 0:
        raise RuntimeError("no marker/oscillation events found before this group's first leg")

    # Case 1: the oscillation's own transitions cleared the jump threshold --
    # a run of >=28 short (<=80ms) events spaced 0.8-1.3s apart.
    micro = starts[durs <= 0.080]
    if micro.size >= 28:
        gaps = np.diff(t[micro])
        good = np.concatenate([[True], (gaps >= 0.8) & (gaps <= 1.3)])
        # longest contiguous run of "good"-spaced micro events
        run_start = run_len = best_start = best_len = 0
        for i, g in enumerate(good):
            if g:
                if run_len == 0:
                    run_start = i
                run_len += 1
            else:
                run_len = 0
            if run_len > best_len:
                best_len, best_start = run_len, run_start
        run = micro[best_start:best_start + best_len]
        if run.size >= 28:
            return t[run[:30]], "detected (oscillation transitions cleared the jump threshold)"

    # Case 2: fall back to marker timing, from two DECOUPLED, narrow searches
    # rather than one 60 s window. A single wide window is fragile: a slow
    # trajectory leg's per-sample motion (as little as 0.275 um/ms at the
    # slowest commanded rate) mostly falls *below* the jump threshold, so its
    # ramp fragments into many small, irregularly-timed events instead of one
    # clean long one, and those fragments can land at the same array
    # positions the CONFIG/OSCILLATION/NEXT marker events would otherwise
    # occupy. Searching the CONFIG/OSCILLATION markers and the NEXT marker in
    # two separate, narrower windows keeps each search well clear of the
    # other's contamination sources.
    #
    # NEXT marker: always close to the leg (observed 13-20 s), so a 25 s
    # window ending at the leg's own threshold-crossing comfortably contains
    # it without reaching into the oscillation block itself.
    near_lo = max(0, leg1_start_idx - int(25.0 / dt_s))
    near_starts, near_durs = big_jump_events(t, pos, dt_s, near_lo, leg1_start_idx,
                                             jump_threshold_um=0.8)
    next_marker = near_starts[(near_durs > 0.080) & (near_durs <= 1.5)]
    if next_marker.size < 1:
        raise RuntimeError("could not find the marker preceding this group's first leg")
    next_marker_onset_t = t[next_marker[0]]

    # CONFIG/OSCILLATION markers: always the first marker-like pair-of-pairs
    # in the wide window, well clear of both the leg ahead and (for MRES 1)
    # the pre-campaign quiet stretch behind.
    far_mask = (durs > 0.080) & (durs <= 1.5) & (t[starts] < next_marker_onset_t - 5.0)
    osc_candidates, osc_candidate_durs = starts[far_mask], durs[far_mask]
    if osc_candidates.size < 2:
        raise RuntimeError("could not find the CONFIG/OSCILLATION marker pair")
    # osc_candidates[-1] is when the OSCILLATION marker's RETURN MOVE STARTS,
    # not when it finishes -- the move itself (its own event duration, ~0.3-
    # 0.55 s here) has to elapse before the settle timer even starts, or osc0
    # lands mid-ramp and the "oscillation" fit captures the marker's own
    # multi-hundred-um return instead of the few-hundred-nm-to-um commanded
    # step (this was checked against the raw signal for MRES 16: without the
    # event's own duration, osc0 landed at -215 um while the axis was still
    # actively ramping from -992 to -195).
    osc_marker_return_start_t = t[osc_candidates[-1]]
    osc_marker_return_duration_s = osc_candidate_durs[-1]
    osc0 = osc_marker_return_start_t + osc_marker_return_duration_s + MARKER_SETTLE_S
    period = (next_marker_onset_t - osc0) / CYCLES
    seeds = np.empty(2 * CYCLES)
    for k in range(CYCLES):
        seeds[2 * k] = osc0 + k * period
        seeds[2 * k + 1] = osc0 + (k + 0.5) * period
    return seeds, f"derived from marker timing (osc0={osc0:.3f}s, period={period:.4f}s)"


# ------------------------------------------------------------ fit: steps ---

def fit_steps(t, y, seed_times, dt_s, median_ms=MEDIAN_WINDOW_MS,
             search_radius_s=0.35):
    """Piecewise-constant fit: seed_times -> local-slope-refined edges ->
    robust median per plateau. Same algorithm as v3's fit_steps."""
    median_n = max(1, round(median_ms * 1e-3 / dt_s))
    median_n += 1 - median_n % 2
    clean = median_filter(y, size=median_n, mode="nearest")
    slope_span = max(1, round(0.006 / dt_s))

    edges = []
    for seed in seed_times:
        lo = max(slope_span, int(np.searchsorted(t, seed - search_radius_s)))
        hi = min(t.size - slope_span, int(np.searchsorted(t, seed + search_radius_s)))
        if hi <= lo:
            continue
        candidates = np.arange(lo, hi, dtype=int)
        score = np.abs(clean[candidates + slope_span] - clean[candidates - slope_span])
        edge = int(candidates[int(np.argmax(score))])
        if edges and edge <= edges[-1]:
            edge = edges[-1] + 1
        edges.append(edge)

    bounds = np.asarray([0, *edges, t.size], dtype=int)
    fitted = np.empty_like(y)
    for left, right in zip(bounds[:-1], bounds[1:]):
        width = right - left
        guard = min(round(0.20 / dt_s), max(0, width // 4))
        inner_left, inner_right = left + guard, right - guard
        plateau = clean[inner_left:inner_right]
        if plateau.size == 0:
            plateau = clean[left:right]
        fitted[left:right] = float(np.median(plateau))
    return fitted, np.asarray(edges)


# --------------------------------------------------------- fit: smoothed ---

def smooth(y, dt_s, median_ms=MEDIAN_WINDOW_MS, sigma_ms=GAUSSIAN_SIGMA_MS):
    """Centered median despike plus centered Gaussian filter. Same as v3."""
    median_n = max(1, round(median_ms * 1e-3 / dt_s))
    median_n += 1 - median_n % 2
    sigma_n = sigma_ms * 1e-3 / dt_s
    despiked = median_filter(y, size=median_n, mode="nearest")
    return gaussian_filter1d(despiked, sigma=sigma_n, mode="nearest",
                             radius=max(1, int(np.ceil(4 * sigma_n))))


# --------------------------------------------------------------- styling ---

def style(ax):
    ax.grid(True, color="#b6b6b6", alpha=.38, linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Time (s)", fontsize=8)
    ax.set_ylabel("Position [µm]", fontsize=8)


# ------------------------------------------------------------------ main ---

def process_group(t, pos, dt_s, mres, group_spans, still, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    leg1_idx = group_spans[0][0]

    seed_times, method = find_oscillation_seed_times(t, pos, dt_s, leg1_idx)
    print(f"MRES {mres}: oscillation seed method: {method}")

    # Start just before osc0, not a full second before: the marker's own
    # return ramp can still be in progress a second out (its move duration is
    # added into osc0 already, but only up to the settle point), and pulling
    # that tail into the window would land it in fit_steps's "before the
    # first edge" plateau, contaminating the very first segment.
    # End just under one dwell (OSCILLATION_HALF_DWELL_S = 1.0 s) past the
    # last seed time, not 1.5 s: beyond one dwell is the next marker's own
    # move, which the light "dynamic" filter below no longer smooths away.
    osc_lo = max(0, int(np.searchsorted(t, seed_times[0] - 0.15)))
    osc_hi = min(pos.size, int(np.searchsorted(t, seed_times[-1] + 0.85)))
    osc_t0 = t[osc_lo]
    osc_local_t = t[osc_lo:osc_hi] - osc_t0
    osc_measured = pos[osc_lo:osc_hi].copy()
    # Baseline from a short window ending right at osc0, not the start of the
    # (wider, for visual context) extraction window: seed_times[0] - 1.0 can
    # still be inside the preceding marker's own return ramp (it only has to
    # be well clear by the time the oscillation itself starts), and a median
    # taken there would zero against the wrong reference entirely.
    base_lo, base_hi = np.searchsorted(t, (seed_times[0] - 0.30, seed_times[0] - 0.02))
    baseline = float(np.median(pos[base_lo:base_hi])) if base_hi > base_lo else \
        float(np.median(osc_measured[:min(50, osc_measured.size)]))
    osc_measured -= baseline
    # "Reference": the idealized flat-plateau fit -- still the cleanest read
    # of the commanded step size, kept for that, not as the primary curve.
    osc_reference, osc_edges = fit_steps(osc_local_t, osc_measured,
                                         seed_times - osc_t0, dt_s)
    # "Dynamic": a light despike + light Gaussian pass over the WHOLE block
    # (not per-plateau, so nothing is flattened) -- this is what preserves
    # the overshoot/ring-down right after each commanded transition, which
    # `osc_reference` deliberately discards via its 0.2 s edge guard.
    osc_dynamic = smooth(osc_measured, dt_s, median_ms=OSC_DYNAMIC_MEDIAN_MS,
                         sigma_ms=OSC_DYNAMIC_SIGMA_MS)
    seed_local = seed_times - osc_t0

    # y-limits from the dynamic fit WITHIN THE CORE OSCILLATION REGION ONLY
    # (it can legitimately overshoot past the reference plateaus -- that's
    # the point), not the padded extraction window: osc_hi runs 1.5 s past
    # the last seed time for plotting context, which can reach into the next
    # marker's own fast excursion. That used to matter to `osc_reference`
    # (which had here it out, contaminates the plateau median a little); now
    # that osc_dynamic is only lightly smoothed (2-3 ms), the marker's own
    # multi-hundred-um move survives at near full amplitude if it's included,
    # and blew the axis out to +-100 um for a ~10 um MRES 1 oscillation.
    core = (osc_local_t >= seed_local[0] - 0.05) & (osc_local_t <= seed_local[-1] + 0.55)
    core_vals = osc_dynamic[core] if core.any() else osc_dynamic
    fit_span = float(np.ptp(core_vals))
    fit_mid = float(np.mean(core_vals))
    half_range = max(0.75 * fit_span, 0.05) * 1.6
    ylim = (fit_mid - half_range, fit_mid + half_range)

    # Oscillation overlay: full 30 s block, measured vs both fits.
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.plot(osc_local_t, osc_measured, color=RAW_COLOR, alpha=.5, linewidth=.5,
            label="Measured")
    ax.step(osc_local_t, osc_reference, where="post", color=REF_COLOR, linewidth=1.0,
            linestyle="--", alpha=.8, label="Reference (idealized plateau)")
    ax.plot(osc_local_t, osc_dynamic, color=FIT_COLOR, linewidth=1.3,
            label="Dynamic fit (overshoot preserved)")
    ax.set_ylim(*ylim)
    ax.set_title(f"MRES {mres} — one-microstep oscillation: measured vs fit",
                loc="left", fontsize=10)
    style(ax)
    ax.legend(loc="best", fontsize=8, framealpha=.9)
    fig.tight_layout()
    fig.savefig(out_dir / "oscillation_fit_overlay.png", dpi=160)
    plt.close(fig)

    # Zoom: first 6 cycles, with seed times marked.
    fig, ax = plt.subplots(figsize=(12, 4.2))
    zoom_end = seed_local[min(11, seed_local.size - 1)] + 0.6
    vis = osc_local_t <= zoom_end
    ax.plot(osc_local_t[vis], osc_measured[vis], color=RAW_COLOR, alpha=.7,
            linewidth=.6, marker=".", markersize=3, markeredgewidth=0,
            label="Measured (1 kHz)")
    ax.step(osc_local_t[vis], osc_reference[vis], where="post", color=REF_COLOR,
            linewidth=1.2, linestyle="--", alpha=.8,
            label="Reference (idealized plateau)")
    ax.plot(osc_local_t[vis], osc_dynamic[vis], color=FIT_COLOR, linewidth=1.8,
            label="Dynamic fit (overshoot preserved)")
    for s in seed_local[seed_local <= zoom_end]:
        ax.axvline(s, color=SEED_COLOR, alpha=.5, linewidth=.9, linestyle="--")
    ax.set_xlim(0, zoom_end)
    ax.set_ylim(*ylim)
    ax.set_title(f"MRES {mres} — first six cycles (orange dashed = seed edge times)",
                loc="left", fontsize=10)
    style(ax)
    ax.legend(loc="best", fontsize=8, framealpha=.9)
    fig.tight_layout()
    fig.savefig(out_dir / "oscillation_microstep_zoom.png", dpi=160)
    plt.close(fig)

    # Trajectory legs: 6 panels. Each gets the same two fits as the
    # oscillation: "reference" (the old heavy despike+Gaussian, an idealized
    # mean trend -- kept thin/dashed for context) and "dynamic" (a light,
    # per-leg-adaptive pass that keeps the individual commanded microsteps
    # riding on the ramp visible, per ramp_dynamic_params).
    leg_series = []
    apex_local_ts = []
    dynamic_params_by_label = {}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    labels = ["D-slow", "D-mod", "D-fast", "I-slow", "I-mod", "I-fast"]
    for ax, (s0, apex, e0), label in zip(axes.ravel(), group_spans, labels):
        a, b = ramp_extent(still, s0, e0, dt_s, pad_s=1.0)
        local_t = t[a:b] - t[a]
        measured = pos[a:b].copy()
        base = float(np.median(measured[:min(50, measured.size)]))
        measured -= base
        reference = smooth(measured, dt_s)
        step_um = FULL_STEP_UM / mres
        velocity_um_s = RATE_BY_LABEL[label] * FULL_STEP_UM
        dyn_median_ms, dyn_sigma_ms = ramp_dynamic_params(step_um, velocity_um_s)
        dynamic_params_by_label[label] = (dyn_median_ms, dyn_sigma_ms)
        dynamic = smooth(measured, dt_s, median_ms=dyn_median_ms, sigma_ms=dyn_sigma_ms)
        leg_series.append((label, local_t, measured, reference, dynamic))
        apex_local_ts.append(t[apex] - t[a])

        ax.plot(local_t, measured, color=RAW_COLOR, alpha=.5, linewidth=.5,
                label="Measured")
        ax.plot(local_t, reference, color=REF_COLOR, linewidth=.9, linestyle="--",
                alpha=.8, label="Reference (idealized trend)")
        ax.plot(local_t, dynamic, color=FIT_COLOR, linewidth=1.1,
                label="Dynamic fit (steps preserved)")
        ax.set_title(f"MRES {mres}  {label}", fontsize=9)
        style(ax)
    axes.ravel()[0].legend(loc="best", fontsize=7, framealpha=.9)
    fig.suptitle(f"MRES {mres} — 10 mm trajectory legs: measured and both fits",
                fontsize=12)
    fig.text(.5, .01,
            "dashed = idealized trend (heavy filter); solid = dynamic fit "
            "(light, per-leg-adaptive filter, individual steps preserved)",
            ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, .95))
    fig.savefig(out_dir / "trajectory_fits_overlay.png", dpi=160)
    plt.close(fig)

    # Two zoomed views per leg: at full extent the fit sits almost exactly on
    # top of the measured trace (the point of the filter working well), which
    # makes it impossible to see *how* well from the overlay above alone.
    #   - apex zoom: the direction reversal, where a smoothing filter is most
    #     likely to show lag or over-rounding if it has any.
    #   - mid-ramp zoom: a plain constant-velocity stretch (40% of the way
    #     from start to apex, clear of both the start transient and the
    #     turnaround), windowed to a target NUMBER OF COMMANDED STEPS rather
    #     than a fixed duration. A fixed time window is exactly why this was
    #     indistinguishable at fine MRES: at MRES 32 a commanded step is 32x
    #     smaller than at MRES 1, so the same 1 s window that clearly shows
    #     ~4 steps at MRES 1 packs in ~128 at MRES 32, far too dense to
    #     resolve individually at any zoom. Sizing the window to the
    #     commanded step size and the leg's own commanded rate keeps the
    #     number of visible transitions comparable across every MRES/rate
    #     combination instead.
    TARGET_STEPS_EACH_SIDE = 10
    MIN_SAMPLES_EACH_SIDE = 40   # floor so a fast+fine combo doesn't collapse
                                 # to a handful of samples with nothing to see
    for zoom_name, zoom_title, center_frac_of_apex, out_stem in (
            ("apex", "direction reversal", None, "trajectory_fits_apex_zoom"),
            ("mid-ramp", "constant-velocity stretch", 0.40, "trajectory_fits_midramp_zoom")):
        fig, axes = plt.subplots(2, 3, figsize=ZOOM_FIGSIZE)
        for ax, (label, local_t, measured, reference, dynamic), apex_t in zip(
                axes.ravel(), leg_series, apex_local_ts):
            center = apex_t if center_frac_of_apex is None else center_frac_of_apex * apex_t
            if center_frac_of_apex is None:
                half_win = min(0.8, 0.10 * apex_t)
            else:
                step_um = FULL_STEP_UM / mres
                velocity_um_s = RATE_BY_LABEL[label] * FULL_STEP_UM
                half_win = TARGET_STEPS_EACH_SIDE * step_um / velocity_um_s
            half_win = max(half_win, MIN_SAMPLES_EACH_SIDE * dt_s)
            vis = (local_t >= center - half_win) & (local_t <= center + half_win)
            ax.plot(local_t[vis], measured[vis], color=RAW_COLOR, alpha=.75,
                    linewidth=.7, marker=".", markersize=4.5,
                    markeredgewidth=0, label="Measured (1 kHz)")
            ax.plot(local_t[vis], reference[vis], color=REF_COLOR, linewidth=1.3,
                    linestyle="--", alpha=.8, label="Reference (idealized trend)")
            ax.plot(local_t[vis], dynamic[vis], color=FIT_COLOR, linewidth=2.2,
                    label="Dynamic fit (steps preserved)")
            if center_frac_of_apex is None:
                ax.axvline(apex_t, color=SEED_COLOR, alpha=.4, linewidth=.9,
                          linestyle="--")
            ax.set_xlim(center - half_win, center + half_win)
            n_vis = int(vis.sum())
            dyn_median_ms, dyn_sigma_ms = dynamic_params_by_label[label]
            ax.set_title(f"MRES {mres}  {label}  ({n_vis} samples; dynamic "
                        f"filter {dyn_median_ms:.1f}/{dyn_sigma_ms:.1f} ms)",
                        fontsize=9.5)
            style(ax)
        axes.ravel()[0].legend(loc="best", fontsize=8, framealpha=.9)
        fig.suptitle(f"MRES {mres} — trajectory legs zoomed on the {zoom_title}",
                    fontsize=13)
        window_note = ("fixed window" if center_frac_of_apex is None else
                       f"window sized to ~{TARGET_STEPS_EACH_SIDE} commanded "
                       "steps each side of center, floored at "
                       f"{MIN_SAMPLES_EACH_SIDE} samples each side")
        fig.text(.5, .01,
                "dynamic fit: median despike + Gaussian, both sized per-panel "
                "to a fraction of that leg's own commanded step interval "
                f"(title); dots are individual 1 kHz samples; {window_note}",
                ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .03, 1, .95))
        savefig_png_and_svg(fig, out_dir, out_stem)
        plt.close(fig)

    # Arrays: both fits saved for every block, "dynamic" is the one to fit a
    # simulation model against, "reference" is the idealized commanded shape.
    arrays = {
        "oscillation_time_s": osc_local_t, "oscillation_measured_um": osc_measured,
        "oscillation_reference_um": osc_reference, "oscillation_dynamic_um": osc_dynamic,
        "oscillation_seed_times_s": seed_local,
    }
    for label, local_t, measured, reference, dynamic in leg_series:
        arrays[f"leg_{label}_time_s"] = local_t
        arrays[f"leg_{label}_measured_um"] = measured
        arrays[f"leg_{label}_reference_um"] = reference
        arrays[f"leg_{label}_dynamic_um"] = dynamic
    np.savez_compressed(out_dir / "fitted_blocks.npz", **arrays)

    # Summary.
    summary = {
        "mres": mres,
        "source": str(CSV),
        "sample_period_ms": 1e3 * dt_s,
        "sample_rate_hz": 1.0 / dt_s,
        "oscillation": {
            "reference": {
                "method": "measurement-aligned piecewise-constant plateaus",
                "seed_method": method,
                "median_despike_window_ms": MEDIAN_WINDOW_MS,
                "edge_refinement": "maximum local encoder slope within 0.35 s of seed",
                "plateau_estimator": "median, 0.2 s guard from each edge",
            },
            "dynamic": {
                "method": "centered median despike plus Gaussian, light pass "
                          "over the whole block (overshoot/ring-down preserved)",
                "median_despike_window_ms": OSC_DYNAMIC_MEDIAN_MS,
                "gaussian_sigma_ms": OSC_DYNAMIC_SIGMA_MS,
            },
            "n_cycles": CYCLES,
        },
        "trajectory_legs": {
            "reference": {
                "method": "centered median despike plus Gaussian smoother "
                          "(idealized mean trend)",
                "median_window_ms": MEDIAN_WINDOW_MS,
                "gaussian_sigma_ms": GAUSSIAN_SIGMA_MS,
            },
            "dynamic": {
                "method": "centered median despike plus Gaussian, both "
                          "parameters scaled per leg to a fraction of that "
                          "leg's own commanded step interval (individual "
                          "microsteps preserved)",
                "median_frac_of_step_interval": RAMP_DYNAMIC_MEDIAN_FRAC,
                "sigma_frac_of_step_interval": RAMP_DYNAMIC_SIGMA_FRAC,
                "median_bounds_ms": list(RAMP_DYNAMIC_MEDIAN_BOUNDS_MS),
                "sigma_bounds_ms": list(RAMP_DYNAMIC_SIGMA_BOUNDS_MS),
                "params_by_leg_ms": {lbl: {"median_ms": m, "sigma_ms": s}
                                    for lbl, (m, s) in dynamic_params_by_label.items()},
            },
            "phase": "zero",
            "legs": labels,
        },
        "boundary_handling": "each block (oscillation, each leg) fitted independently",
        "outputs": {
            "oscillation_overlay": "oscillation_fit_overlay.png",
            "oscillation_zoom": "oscillation_microstep_zoom.png",
            "trajectory_overlay": "trajectory_fits_overlay.png",
            "trajectory_apex_zoom": "trajectory_fits_apex_zoom.png/.svg",
            "trajectory_midramp_zoom": "trajectory_fits_midramp_zoom.png/.svg",
            "arrays": "fitted_blocks.npz",
        },
    }
    (out_dir / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                               encoding="utf-8")
    print(f"MRES {mres}: wrote oscillation_fit_overlay.png, "
          f"oscillation_microstep_zoom.png, trajectory_fits_overlay.png, "
          f"fitted_blocks.npz, fit_summary.json -> {out_dir}")


def main():
    apply_style()
    t, pos, meta = load_ids(CSV)
    dt_s = meta["dt_s"]
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz")

    spans = find_trajectories(t, pos, dt_s)
    still = still_mask(pos, dt_s)
    if len(spans) != 6 * len(MRES_VALUES):
        raise RuntimeError(f"expected {6*len(MRES_VALUES)} legs, found {len(spans)}")

    for gi, mres in enumerate(MRES_VALUES):
        group_spans = spans[gi * 6:(gi + 1) * 6]
        out_dir = OUT_ROOT / f"MRES_{mres}"
        process_group(t, pos, dt_s, mres, group_spans, still, out_dir)


if __name__ == "__main__":
    main()

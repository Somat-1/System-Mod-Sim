"""v2: the post-recording idle tail — a pure noise-floor characterization.

`EXPrun.csv` runs 1302.3 s; the last commanded motion (leg #24 of 24, MRES 32
I-fast) ends at t=1147.06 s, but the capture keeps recording until t=1302.3 s.
That last ~154 s has no commanded motion in it at all -- nothing is supposed
to move -- which makes it the cleanest available segment for asking "how much
does the measured position vary when nothing is happening?" Every other quiet
stretch in this capture sits between two commanded moves and inherits some of
their settling; this one doesn't.

The boundary is found the same way `ramp_extent`/`still_mask` are already
used elsewhere in v6 (see ids_common.py): scan backward from the end of the
capture for the earliest sample after which the axis is continuously
"still" (windowed slope below threshold) all the way to the last sample.

Two deliverables:
  1. A static zoom cascade: the full tail, then progressively narrower
     windows (10 s, 1 s, 0.1 s) all centered on the same instant, each
     panel showing where the next one comes from.
  2. A self-contained interactive Plotly HTML of the same segment at full
     sample resolution (no decimation -- 153,651 points renders smoothly in
     WebGL directly) so the window/center choices above aren't the only way
     to look at it.

Noise metrics (§ compute_noise_metrics) are reported for the whole tail and
for a set of local windows, because the two answer different questions: the
whole-tail statistics can't distinguish genuine sample-to-sample noise from
slow drift over 154 s, while short local windows and the mean absolute
step-to-step difference isolate the high-frequency contribution specifically.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "analysis" / "scripts"))
from ids_common import load_ids, still_mask, apply_style, C_MEASURED, C_COMMANDED, C_INK, C_MUTED, C_GRID
from plot_v6_exp_run import CSV, find_trajectories

OUT_DIR = Path(__file__).resolve().parent.parent / "rendered_assets"
ZOOM_LEVELS_S = [None, 10.0, 1.0, 0.1]   # None = full tail
LOCAL_WINDOW_S = 1.0                      # for the windowed-std noise estimate
START_REF_SAMPLES = 25                    # leading samples behind the "starts at zero" reference


def window_start_level(y):
    """The settled position level at a window's start, as a median over its
    first few samples rather than its literal first sample.

    The full tail's first sample lands inside the ~4-sample transient that
    marks the boundary (see large_outliers_gt_300nm), and is itself the
    outlier: the axis settles ~250 nm away from it within 20 ms. Zeroing to
    that one sample therefore puts the whole settled trace ~250 nm off
    zero. A short median is robust to that handful of samples and gives the
    level the window actually starts at.
    """
    k = max(1, min(y.size, START_REF_SAMPLES))
    return float(np.median(y[:k]))


def find_post_recording_tail(t, pos, dt_s):
    """(start_idx, end_idx) of the longest continuously-still stretch running
    all the way to the end of the capture, starting the search after the
    last commanded trajectory leg."""
    spans = find_trajectories(t, pos, dt_s)
    last_e0 = spans[-1][2]
    still = still_mask(pos, dt_s)
    tail_still = still[last_e0:]
    rev = tail_still[::-1]
    not_still = np.flatnonzero(~rev)
    if not_still.size:
        start_idx = last_e0 + (tail_still.size - not_still[0])
    else:
        start_idx = last_e0
    return start_idx, pos.size, spans


def compute_noise_metrics(t, pos, dt_s, start_idx, end_idx):
    seg = pos[start_idx:end_idx]
    seg_t = t[start_idx:end_idx]
    diffs = np.abs(np.diff(seg))

    # Whole-tail statistics: mixes genuine noise with any slow drift.
    whole = {
        "n_samples": int(seg.size),
        "duration_s": float(seg_t[-1] - seg_t[0]),
        "mean_um": float(seg.mean()),
        "std_um": float(seg.std(ddof=1)),
        "ptp_um": float(np.ptp(seg)),
        "min_um": float(seg.min()),
        "max_um": float(seg.max()),
    }

    # Mean absolute sample-to-sample variation: the literal "average
    # variation in position" -- insensitive to slow drift by construction,
    # since it only ever looks one sample ahead.
    step_variation = {
        "mean_abs_diff_um": float(diffs.mean()),
        "median_abs_diff_um": float(np.median(diffs)),
        "max_abs_diff_um": float(diffs.max()),
    }

    # Local (1 s) windowed std: isolates noise from drift on a timescale
    # short enough that drift within one window is negligible, the same
    # "quietest window" convention used elsewhere in this repo for the
    # encoder noise floor.
    win_n = int(round(LOCAL_WINDOW_S / dt_s))
    n_windows = seg.size // win_n
    local_stds = np.array([seg[i * win_n:(i + 1) * win_n].std(ddof=1)
                           for i in range(n_windows)])
    local = {
        "window_s": LOCAL_WINDOW_S,
        "n_windows": int(n_windows),
        "mean_std_um": float(local_stds.mean()),
        "min_std_um": float(local_stds.min()),
        "max_std_um": float(local_stds.max()),
    }

    # Linear drift: is the axis still creeping after motion has stopped?
    slope_um_s, intercept_um = np.polyfit(seg_t - seg_t[0], seg, 1)
    drift = {
        "slope_um_per_s": float(slope_um_s),
        "total_drift_um_over_tail": float(slope_um_s * whole["duration_s"]),
    }

    # Large step-to-step outliers: whole_tail.ptp_um and
    # step_to_step_variation.max_abs_diff_um can both be dominated by a
    # single event rather than being representative of the steady noise
    # level -- report them explicitly rather than let them hide inside an
    # average. Threshold (300 nm) is well above the ~5-9 nm steady-state
    # scale found below but well below a genuine commanded step (>=150 nm
    # even at MRES 32).
    outlier_mask = diffs > 0.300
    outliers = [{"t_s": float(seg_t[i]), "diff_um": float(diffs[i])}
               for i in np.flatnonzero(outlier_mask)]

    # Steady-state: the same whole-tail statistics recomputed after skipping
    # the first 1 s, to check whether excluding it (a) removes the outliers
    # above and (b) brings std down toward the local-window figure. If it
    # doesn't move much, the elevated whole-tail std is drift/wander over the
    # full 154 s, not a settling transient -- worth knowing before quoting
    # "the" noise number.
    steady_mask = (seg_t - seg_t[0]) >= 1.0
    steady_seg = seg[steady_mask]
    steady = {
        "excludes_first_s": 1.0,
        "n_samples": int(steady_seg.size),
        "std_um": float(steady_seg.std(ddof=1)),
        "ptp_um": float(np.ptp(steady_seg)),
        "mean_abs_diff_um": float(np.abs(np.diff(steady_seg)).mean()),
    }

    return {"whole_tail": whole, "step_to_step_variation": step_variation,
            "local_1s_windows": local, "linear_drift": drift,
            "large_outliers_gt_300nm": outliers, "steady_state_excl_first_1s": steady}


def style(ax):
    ax.grid(True, color=C_GRID, alpha=.7, linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Time from start of capture [s]", fontsize=9)
    ax.set_ylabel("Relative position [nm]\n(rel. to this window's own start level)", fontsize=9)
    ax.tick_params(labelsize=8)


def plot_zoom_cascade(t, pos, start_idx, end_idx, metrics, out_path):
    apply_style()
    seg_t0, seg_t1 = t[start_idx], t[end_idx - 1]
    center = 0.5 * (seg_t0 + seg_t1)

    windows = []
    for lvl in ZOOM_LEVELS_S:
        if lvl is None:
            windows.append((seg_t0, seg_t1))
        else:
            windows.append((center - lvl / 2, center + lvl / 2))

    # Single column, 4 stacked panels -- full tail down to the narrowest zoom.
    fig, axes = plt.subplots(len(windows), 1, figsize=(12, 4.3 * len(windows)))

    for i, (ax, (lo, hi)) in enumerate(zip(axes, windows)):
        mask = (t >= lo) & (t <= hi)
        # Each window is normalized to ITS OWN first visible sample -- not a
        # shared reference -- so every panel starts at exactly 0 regardless
        # of where in the tail it sits.
        yy = (pos[mask] - window_start_level(pos[mask])) * 1000.0   # µm -> nm, zeroed to window start
        local_std = float(yy.std(ddof=1)) if yy.size > 1 else float("nan")
        local_mean = float(yy.mean()) if yy.size else float("nan")
        # Path simplification off: the wide panels compress the first
        # milliseconds into a fraction of a pixel, and the default
        # simplification thins that transient out of the rendered line.
        with plt.rc_context({"path.simplify": False}):
            ax.plot(t[mask], yy, color=C_MEASURED, lw=.8 if i == 0 else 1.1,
                    marker="." if i == len(windows) - 1 else None, markersize=3)
        ax.axhline(local_mean, color=C_COMMANDED, lw=1.1, ls=":", alpha=.85,
                  label=f"window mean ({local_mean:+.1f} nm)")
        ax.axhline(0.0, color=C_MUTED, lw=.7, ls="-", alpha=.5,
                  label="window start level")

        # Centered, tightly-fit y-range, from a robust (1st-99th percentile)
        # spread rather than the true min/max -- so the ~4-sample transient
        # right at the tail's own start (see large_outliers_gt_300nm) can't
        # single-handedly flatten the full-tail panel's view of the drift and
        # noise it's meant to make visible. 0 (the window's own start, always
        # its first plotted point) is forced into the range regardless, so
        # the panel visibly starts at zero even when that first sample is
        # more extreme than the rest of the window's 1st-99th percentile
        # spread -- only clip the tail of the outlier, never the start point.
        p_lo, p_hi = np.percentile(yy, [1, 99])
        y_lo, y_hi = min(p_lo, 0.0), max(p_hi, 0.0)
        y_mid, y_span = 0.5 * (y_lo + y_hi), max(y_hi - y_lo, 0.5)
        pad = 0.25 * y_span
        ax.set_ylim(y_mid - y_span / 2 - pad, y_mid + y_span / 2 + pad)

        width_label = "full tail" if ZOOM_LEVELS_S[i] is None else f"{ZOOM_LEVELS_S[i]:g} s window"
        ax.set_title(f"{width_label}  ·  mean = {local_mean:+.1f} nm  ·  "
                    f"std = {local_std:.1f} nm  ·  n = {yy.size}",
                    loc="left", fontsize=9.5)
        style(ax)
        ax.set_xlim(lo, hi)
        ax.legend(loc="upper right", fontsize=6.5, framealpha=.9)
        # Mark where the next (narrower) window sits, on this panel.
        if i + 1 < len(windows):
            nlo, nhi = windows[i + 1]
            ax.axvspan(nlo, nhi, color=C_COMMANDED, alpha=.12, lw=0)
            ax.axvline(nlo, color=C_COMMANDED, lw=.8, alpha=.5)
            ax.axvline(nhi, color=C_COMMANDED, lw=.8, alpha=.5)

    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_interactive(t, pos, start_idx, end_idx, metrics, out_path):
    seg_t = t[start_idx:end_idx]
    entry_val = float(pos[start_idx])
    # Same normalization as the cascade: relative to the tail's own first
    # sample, so 0 on this axis always means "back at the entry point" --
    # and in nm, not µm, to match.
    seg_y = (pos[start_idx:end_idx] - entry_val) * 1000.0
    mean_nm = float(seg_y.mean())          # this window's own mean (there's
    std_nm = metrics["whole_tail"]["std_um"] * 1000.0   # only one window here: the whole tail)

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=seg_t, y=seg_y, mode="lines",
        line=dict(color=C_MEASURED, width=1),
        name="Measured (EL5101), full sample resolution",
        hovertemplate="t=%{x:.4f} s<br>Δpos=%{y:.2f} nm<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=[seg_t[0], seg_t[-1]], y=[mean_nm, mean_nm], mode="lines",
        line=dict(color=C_COMMANDED, width=1, dash="dash"),
        name=f"Window mean ({mean_nm:+.1f} nm)",
    ))
    fig.add_trace(go.Scatter(
        x=[seg_t[0], seg_t[-1]], y=[0.0, 0.0], mode="lines",
        line=dict(color=C_MUTED, width=1),
        name="Tail entry point (t=%.3fs)" % seg_t[0],
    ))
    fig.add_trace(go.Scatter(
        x=[seg_t[0], seg_t[-1], seg_t[-1], seg_t[0]],
        y=[mean_nm - std_nm, mean_nm - std_nm, mean_nm + std_nm, mean_nm + std_nm],
        fill="toself", mode="none", fillcolor="rgba(217,72,1,0.08)",
        name=f"±1σ ({std_nm:.1f} nm)", hoverinfo="skip",
    ))

    fig.update_layout(
        title=dict(text=("EXPrun — post-recording idle tail, interactive"
                         f"<br><sub>{seg_y.size:,} samples, {seg_t[-1]-seg_t[0]:.1f} s, "
                         "no commanded motion, full resolution (not decimated), "
                         "normalized to the tail's entry point  ·  "
                         f"window mean {mean_nm:+.1f} nm  ·  std {std_nm:.1f} nm  ·  "
                         f"mean |Δsample| "
                         f"{metrics['step_to_step_variation']['mean_abs_diff_um']*1000:.2f} nm  ·  "
                         "drag on the chart to box-zoom both axes, scroll to zoom, "
                         "double-click to reset</sub>"),
                  x=0.01, xanchor="left"),
        dragmode="zoom",
        xaxis=dict(title="Time from start of capture [s]",
                  rangeslider=dict(visible=True), showspikes=True),
        yaxis=dict(title="Relative position [nm] (rel. to tail entry point)",
                  showspikes=True),
        template="plotly_white",
        hovermode="x unified",
        margin=dict(t=100),
    )
    fig.write_html(out_path, include_plotlyjs=True, full_html=True,
                   config=dict(scrollZoom=True, doubleClick="reset", displaylogo=False))


def main():
    t, pos, meta = load_ids(CSV)
    dt_s = meta["dt_s"]
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz")

    start_idx, end_idx, spans = find_post_recording_tail(t, pos, dt_s)
    print(f"last commanded leg ends at t={t[spans[-1][2]]:.3f}s")
    print(f"post-recording tail: t=[{t[start_idx]:.3f}, {t[end_idx-1]:.3f}]s "
          f"({t[end_idx-1]-t[start_idx]:.2f}s, {end_idx-start_idx:,} samples)")

    metrics = compute_noise_metrics(t, pos, dt_s, start_idx, end_idx)
    print(json.dumps(metrics, indent=2))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_zoom_cascade(t, pos, start_idx, end_idx, metrics, OUT_DIR / "tail_zoom_cascade.png")
    print("wrote tail_zoom_cascade.png")
    plot_interactive(t, pos, start_idx, end_idx, metrics, OUT_DIR / "tail_interactive.html")
    print("wrote tail_interactive.html")

    summary = {
        "source": str(CSV),
        "tail_start_s": float(t[start_idx]),
        "tail_end_s": float(t[end_idx - 1]),
        "last_commanded_leg_end_s": float(t[spans[-1][2]]),
        "sample_rate_hz": meta["fs_hz"],
        "metrics": metrics,
    }
    (OUT_DIR.parent / "noise_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("wrote ../noise_summary.json")


if __name__ == "__main__":
    main()

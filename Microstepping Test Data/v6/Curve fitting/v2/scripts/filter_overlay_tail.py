"""v2 filter overlay: three candidate smoothers applied to the post-recording
idle tail, rendered on top of each other and the raw trace for comparison.

Same segment as analyze_post_recording_tail.py (the 153.65 s post-recording
tail -- see that module's docstring for why it's the cleanest available
no-motion stretch in EXPrun.csv). Three methods are applied directly to the
raw samples of that segment:

  1. Median despike            M_4(raw)          -- running median, width 4 samples
  2. Despike + Gaussian        G_{10,1.67}(M_4(raw)) -- centered Gaussian on top of (1),
                                                       window 10 samples, sigma = W/6.
                                                       Sized from the tail's own measured
                                                       ~363 Hz resonance rather than
                                                       assumed -- see ../temp/.
  3. Cubic B-spline            LSQ fit directly to raw, interior knots spaced
                                every ~750 samples (requested range: 500-1000)

All three are computed independently from the same raw segment (2 does NOT
feed into 3) so each represents a distinct, directly comparable candidate
smoother -- not a sequential pipeline.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.interpolate import make_lsq_spline

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "analysis" / "scripts"))
from ids_common import load_ids, apply_style, C_MEASURED, C_INK, C_MUTED, C_GRID
from plot_v6_exp_run import CSV
from analyze_post_recording_tail import find_post_recording_tail, window_start_level

OUT_DIR = Path(__file__).resolve().parent.parent / "rendered_assets"

# --- filter parameters, as specified ---------------------------------------
MEDIAN_WIDTH_SAMPLES = 4
# Sized from the measured resonance, not assumed: the Welch PSD of this same
# idle tail peaks at fc ~= 363 Hz (../temp/psd_resonance_probe.py; stable to
# +/-3 Hz across 9 of 10 independent 10 s blocks, ~100x above the noise floor).
# At 2.74 samples/cycle the ripple is too marginally sampled for phase
# cancellation to work, so the window is widened past the 2-cycle rule until
# the cutoff simply rejects the band: W = 10 puts -3 dB at 112 Hz and removes
# 99.9% of fc, while staying two decades above the <5 Hz thermal wander.
GAUSSIAN_WINDOW_SAMPLES = 10          # centered kernel support
GAUSSIAN_SIGMA_SAMPLES = 10 / 6       # = W/6 ~= 1.667 samples (3-sigma rule)
SPLINE_KNOT_SPACING_SAMPLES = 750     # requested range: 500-1000 samples; 750 used here
SPLINE_DEGREE = 3                     # cubic

# Same window cascade as tail_zoom_cascade.png, so the two figures line up
# panel for panel. None = the full tail; the rest are centered on its midpoint.
ZOOM_LEVELS_S = [None, 10.0, 1.0, 0.1]

C_DESPIKE = "#8a4fd6"
C_SMOOTH = "#d94801"
C_SPLINE = "#0f9d58"


def median_despike(y, width):
    return median_filter(y, size=width, mode="nearest")


def despike_gaussian(y, width, win, sigma):
    despiked = median_despike(y, width)
    # truncate=3.0 IS the 3-sigma rule: cut at +/-3 sigma = +/-W/2, so the
    # kernel spans exactly the target window and tapers to ~zero at its bounds.
    return gaussian_filter1d(despiked, sigma=sigma, truncate=3.0, mode="nearest")


def cubic_bspline_fit(t, y, knot_spacing):
    k = SPLINE_DEGREE
    interior_idx = np.arange(knot_spacing, t.size - knot_spacing, knot_spacing)
    interior_knots = t[interior_idx]
    t_knots = np.r_[(t[0],) * (k + 1), interior_knots, (t[-1],) * (k + 1)]
    spl = make_lsq_spline(t, y, t_knots, k=k)
    return spl(t), interior_knots.size


def style(ax):
    ax.grid(True, color=C_GRID, alpha=.7, linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Time from start of capture [s]", fontsize=9)
    ax.set_ylabel("Relative position [nm]\n(rel. to tail's own start)", fontsize=9)
    ax.tick_params(labelsize=8)


def main():
    t, pos, meta = load_ids(CSV)
    dt_s = meta["dt_s"]
    fs_hz = meta["fs_hz"]
    start_idx, end_idx, _ = find_post_recording_tail(t, pos, dt_s)

    seg_t = t[start_idx:end_idx]
    seg = pos[start_idx:end_idx]           # µm, raw
    n = seg.size

    despiked = median_despike(seg, MEDIAN_WIDTH_SAMPLES)
    smoothed = despike_gaussian(seg, MEDIAN_WIDTH_SAMPLES,
                                GAUSSIAN_WINDOW_SAMPLES, GAUSSIAN_SIGMA_SAMPLES)
    spline_y, n_interior_knots = cubic_bspline_fit(seg_t, seg, SPLINE_KNOT_SPACING_SAMPLES)

    # Sanity check on the Gaussian's stated "140 Hz": the -3 dB point of a
    # Gaussian kernel with std sigma (samples) at sample rate fs is
    # fc = sqrt(ln2 / (2 pi^2 sigma_s^2)), sigma_s = sigma/fs.
    sigma_s = GAUSSIAN_SIGMA_SAMPLES / fs_hz
    fc_3db = np.sqrt(np.log(2) / (2 * np.pi**2 * sigma_s**2))
    print(f"Gaussian W={GAUSSIAN_WINDOW_SAMPLES}, sigma={GAUSSIAN_SIGMA_SAMPLES:.3f} "
          f"samples @ {fs_hz:.0f} Hz -> -3dB cutoff ~= {fc_3db:.1f} Hz "
          f"(sized to the measured ~363 Hz resonance; see ../temp/)")
    print(f"spline: {n_interior_knots} interior knots over {n} samples "
          f"(~{n / (n_interior_knots + 1):.0f} samples/knot)")

    for name, y in [("raw", seg), ("median despike", despiked),
                     ("despike+gaussian", smoothed), ("cubic B-spline", spline_y)]:
        print(f"  {name:>18s}: std={y.std(ddof=1)*1000:7.2f} nm  "
              f"ptp={np.ptp(y)*1000:8.2f} nm")

    center = 0.5 * (seg_t[0] + seg_t[-1])
    windows = []
    for lvl in ZOOM_LEVELS_S:
        if lvl is None:
            windows.append((seg_t[0], seg_t[-1]))
        else:
            windows.append((center - lvl / 2, center + lvl / 2))

    despike_label = f"median despike (w={MEDIAN_WIDTH_SAMPLES})"
    gauss_label = (f"despike+Gaussian (w={GAUSSIAN_WINDOW_SAMPLES}, "
                   f"$\\sigma$={GAUSSIAN_SIGMA_SAMPLES:.2f})")
    spline_label = f"cubic B-spline (knots~{SPLINE_KNOT_SPACING_SAMPLES} smp)"

    def render(traces, out_name):
        """One figure, four stacked panels (full tail, 10 s, 1 s, 0.1 s), with
        `traces` (list of (y, color, label)) overlaid on the raw samples in
        each. Same window cascade as tail_zoom_cascade.png."""
        apply_style()
        fig, axes = plt.subplots(len(windows), 1, figsize=(13, 4.3 * len(windows)))
        for i, (ax, (lo, hi)) in enumerate(zip(axes, windows)):
            m = (seg_t >= lo) & (seg_t <= hi)
            tt = seg_t[m]
            wide = i == 0
            # Each panel zeroed to its OWN start level (median of its leading
            # samples), matching tail_zoom_cascade.png -- so every panel
            # starts at 0 regardless of where in the tail it sits.
            ref = window_start_level(seg[m])

            def with_std(label, y_win):
                """Label plus that trace's std within THIS window -- the whole
                point of the comparison is how much spread each method leaves,
                and it differs by zoom level (slow wander dominates the wide
                panels, sample noise the narrow ones)."""
                s = y_win.std(ddof=1) * 1000.0 if y_win.size > 1 else float("nan")
                return f"{label}  ·  std = {s:.1f} nm"

            with plt.rc_context({"path.simplify": False}):
                ax.plot(tt, (seg[m] - ref) * 1000.0, color=C_MEASURED,
                        lw=.6 if wide else 1.0, alpha=.55,
                        label=with_std("raw", seg[m]), zorder=1)
                for j, (y, color, label) in enumerate(traces):
                    ax.plot(tt, (y[m] - ref) * 1000.0, color=color,
                            lw=(.8 + .2 * j) if wide else (1.3 + .3 * j),
                            alpha=.93, label=with_std(label, y[m]), zorder=2 + j)
            style(ax)
            ax.set_xlim(lo, hi)

            # Centered, tightly-fit y-range from a robust 1st-99th percentile
            # spread (the ~1.15 um transient at the tail's own start would
            # otherwise flatten the wide panel), always keeping 0 in view.
            yy_ref = (seg[m] - ref) * 1000.0
            p_lo, p_hi = np.percentile(yy_ref, [1, 99])
            y_lo, y_hi = min(p_lo, 0.0), max(p_hi, 0.0)
            y_mid, y_span = 0.5 * (y_lo + y_hi), max(y_hi - y_lo, 0.5)
            pad = 0.25 * y_span
            ax.set_ylim(y_mid - y_span / 2 - pad, y_mid + y_span / 2 + pad)

            width_label = ("full tail" if ZOOM_LEVELS_S[i] is None
                           else f"{ZOOM_LEVELS_S[i]:g} s window")
            ax.set_title(f"{width_label}  ·  n = {tt.size:,}",
                         loc="left", fontsize=9.5)
            ax.legend(loc="upper right", fontsize=7.5, framealpha=.92)
            # Mark where the next (narrower) panel zooms into.
            if i + 1 < len(windows):
                nlo, nhi = windows[i + 1]
                ax.axvspan(nlo, nhi, color=C_INK, alpha=.07, lw=0)
                ax.axvline(nlo, color=C_INK, lw=.7, alpha=.35)
                ax.axvline(nhi, color=C_INK, lw=.7, alpha=.35)

        fig.tight_layout()
        out_path = OUT_DIR / out_name
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print("wrote", out_path)

    render([(despiked, C_DESPIKE, despike_label),
            (smoothed, C_SMOOTH, gauss_label),
            (spline_y, C_SPLINE, spline_label)],
           "filter_overlay_tail.png")


if __name__ == "__main__":
    main()

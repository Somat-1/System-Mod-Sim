"""v2.2: the first 60 s of the idle tail — detrended noise and its spectrum.

`../..`'s v2 characterizes the whole 153.65 s post-recording tail. This
narrows to the **first 60 s** of it and asks a sharper question: with a
low-order trend taken out, what is the residual noise, how is it distributed,
and where does its power sit in frequency?

Three things make that a different analysis from v2's:

1. **The opening transient is excluded.** The tail's first samples contain a
   ~1.15 um drop (v2's `large_outliers_gt_300nm`) that has no mechanical or
   electrical explanation in this dataset and is not part of the noise
   process. `SKIP_S` drops it. The cut is measured rather than assumed: the
   segment's peak-to-peak falls 1271 -> 275 -> 219 -> 107 nm as the skip grows
   to 50 ms and then stops changing entirely (identical at 50, 100, 500 and
   1000 ms), and the last step-to-step jump above 100 nm is at 26 ms. 50 ms is
   the point where the transient is out and nothing further is being discarded.

2. **A low-order polynomial is subtracted.** v2 establishes that the tail
   carries slow wander an order of magnitude larger than a linear drift can
   explain, and that the whole-tail std (37.5 nm) is measuring that wander
   rather than noise. Subtracting a low-order polynomial removes the part of
   it that is a smooth trend over this window, so the residual is closer to
   the stationary process a std, a histogram and a PSD all assume.

   The order is selected, not assumed -- see `select_order`.

3. **The spectrum is estimated with long Welch segments.** `NPERSEG` = 16384
   samples = 16.4 s, Hann window, 50% overlap: 0.061 Hz resolution from 6
   averaged segments. Long segments are the right trade here because the
   feature of interest at the low end is the wander, which lives below a few
   Hz and needs fine resolution to separate from DC at all.

Run from this directory:

    python analyze_tail_60s.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy import signal, stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
V22_DIR = HERE.parent
V2_DIR = V22_DIR.parent
V6_DIR = V2_DIR.parent.parent

sys.path.insert(0, str(V6_DIR / "analysis" / "scripts"))
sys.path.insert(0, str(V2_DIR / "scripts"))

from ids_common import load_ids
from plot_v6_exp_run import CSV
from analyze_post_recording_tail import find_post_recording_tail

OUT_DIR = V22_DIR / "rendered_assets"

SKIP_S = 0.050          # drop the unexplained opening transient (see docstring)
SEGMENT_S = 60.0        # the window this whole version is about
MAX_ORDER = 8           # highest polynomial order evaluated by select_order
KNEE_IMPROVEMENT = 0.01 # <1% further reduction in residual std => stop

# Welch: long segments, Hann, 50% overlap.
NPERSEG = 16384
OVERLAP_FRAC = 0.5
WINDOW = "hann"

# Slots 1-3 of the dataviz reference categorical theme (aqua darkened one step
# to #199e70 so every mark clears 3:1 on white). Validated over all pairs:
# normal-vision, protanopia and deuteranopia OKLab dE all pass. Colour follows
# the entity: blue is the raw measurement, orange the fitted trend, aqua
# everything derived from the residual.
C_RAW = "#2563a6"
C_FIT = "#d55e00"
C_RESID = "#008b6d"
C_INK = "#17202a"
C_INK2 = "#46515c"
C_MUTED = "#77818a"
C_GRID = "#dde3e8"
C_BASELINE = "#aeb7bf"
C_SOFT_BLUE = "#eaf2f8"
C_SOFT_ORANGE = "#fff1e8"
LABEL_BBOX = dict(facecolor="white", edgecolor=C_GRID, linewidth=0.6,
                  alpha=0.94, pad=3.0)


def apply_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": C_BASELINE,
        "axes.labelcolor": C_INK2,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlepad": 9,
        "axes.labelsize": 10,
        "axes.grid": False,
        "grid.color": C_GRID,
        "grid.linewidth": 0.7,
        "grid.alpha": 0.8,
        "xtick.color": C_MUTED,
        "ytick.color": C_MUTED,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.frameon": False,
        "legend.edgecolor": C_GRID,
        "legend.fontsize": 9,
        "font.size": 10,
        "text.color": C_INK,
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def style(ax, xlabel=None, ylabel=None):
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(C_BASELINE)
    ax.grid(axis="y")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


# ------------------------------------------------------------- the segment ---

def extract_segment(t, pos, dt_s):
    """The 60 s window, in seconds-from-window-start and nanometres.

    Zeroed to its own median rather than its first sample: the point of this
    version is the *shape* of the trend and the spread around it, and a
    median-referenced signal keeps the fitted polynomial's constant term a
    meaningful small number instead of the ~-194 um absolute stage offset.
    """
    start, end, _ = find_post_recording_tail(t, pos, dt_s)
    lo = start + int(round(SKIP_S / dt_s))
    hi = lo + int(round(SEGMENT_S / dt_s))
    if hi > end:
        raise RuntimeError(f"tail is shorter than {SKIP_S + SEGMENT_S} s")
    seg_t = t[lo:hi] - t[lo]
    seg_y = (pos[lo:hi] - np.median(pos[lo:hi])) * 1000.0    # um -> nm
    return seg_t, seg_y, dict(tail_start_s=float(t[start]),
                              segment_start_s=float(t[lo]),
                              segment_end_s=float(t[hi - 1]),
                              n_samples=int(hi - lo))


# ---------------------------------------------------------------- detrend ---

def select_order(seg_t, seg_y):
    """Lowest polynomial order past which fitting harder stops paying.

    Sweep 0..MAX_ORDER and take the first order whose *successor* improves the
    residual std by less than KNEE_IMPROVEMENT. On this segment that is a
    cubic: the sweep reads 9.3%, 2.3%, 2.7%, then 0.2% and 0.0%.

    The stopping rule matters because the curve is not monotone in the useful
    sense -- improvement picks up again at orders 6-7 (2.5%, 8.7%). That is
    not a better trend, it is the polynomial acquiring enough freedom to start
    tracing the slow wander itself. Removing the wander is exactly what this
    analysis must NOT do: the wander is a result (v2 quantifies it), and a
    residual with it fitted away would report a noise floor that is partly an
    artifact of the detrend. Stopping at the first knee keeps the subtraction
    a trend removal rather than a fit to the signal of interest.
    """
    stds = []
    for order in range(MAX_ORDER + 1):
        resid = seg_y - Polynomial.fit(seg_t, seg_y, order)(seg_t)
        stds.append(float(np.std(resid)))
    chosen = MAX_ORDER
    for order in range(len(stds) - 1):
        if (stds[order] - stds[order + 1]) / stds[order] < KNEE_IMPROVEMENT:
            chosen = order
            break
    return chosen, stds


def poly_expression(coef, dp=6):
    """p(t) = a0 + a1 t + ... as a readable string, t in seconds, p in nm."""
    terms = []
    for i, c in enumerate(coef):
        if i == 0:
            terms.append(f"{c:+.{dp}g}")
        elif i == 1:
            terms.append(f"{c:+.{dp}g}·t")
        else:
            terms.append(f"{c:+.{dp}g}·t^{i}")
    return "p(t) = " + " ".join(terms)


# ------------------------------------------------------------------ plots ---

def plot_detrend_overlay(seg_t, seg_y, fit_y, resid, order, stds, expr, out_dir):
    """Show the detrend decision in one scan: input, residual, then order choice."""
    raw_sd, raw_ptp = float(np.std(seg_y)), float(np.ptp(seg_y))
    resid_sd, resid_ptp = float(np.std(resid)), float(np.ptp(resid))

    fig = plt.figure(figsize=(13.2, 7.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.10, 1.0],
                          left=0.075, right=0.975, top=0.855, bottom=0.115,
                          hspace=0.43, wspace=0.32)

    ax = fig.add_subplot(gs[0, :])
    ax.plot(seg_t, seg_y, color=C_RAW, linewidth=0.55, alpha=0.82,
            label="Measured position")
    ax.plot(seg_t, fit_y, color=C_FIT, linewidth=2.2,
            label=f"Order-{order} trend")
    ax.set_title("Measured position and fitted trend", loc="left")
    ax.legend(loc="upper left", ncol=2)
    ax.text(0.985, 0.06,
            f"raw  σ {raw_sd:.2f} nm   ·   peak-to-peak {raw_ptp:.1f} nm",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            color=C_INK2, bbox=LABEL_BBOX)
    style(ax, ylabel="Position relative to median [nm]")
    ax.set_xlim(0, seg_t[-1])

    ax = fig.add_subplot(gs[1, :2])
    ax.axhspan(-resid_sd, resid_sd, color=C_RESID, alpha=0.09, zorder=0,
               label="±1σ")
    ax.axhline(0.0, color=C_BASELINE, linewidth=1.0, zorder=1)
    ax.plot(seg_t, resid, color=C_RESID, linewidth=0.55, zorder=2,
            label="Detrended residual")
    ax.set_title("Residual after subtracting the trend", loc="left")
    ax.text(0.985, 0.06,
            f"σ {resid_sd:.2f} nm   ·   peak-to-peak {resid_ptp:.1f} nm",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            color=C_INK2, bbox=LABEL_BBOX)
    style(ax, xlabel="Time from segment start [s]", ylabel="Residual [nm]")
    ax.set_xlim(0, seg_t[-1])

    ax = fig.add_subplot(gs[1, 2])
    orders = np.arange(len(stds))
    ax.plot(orders, stds, color=C_MUTED, linewidth=1.5, marker="o",
            markersize=4.5, markerfacecolor="white", markeredgewidth=1.0,
            zorder=2)
    ax.plot(order, stds[order], marker="o", markersize=9, color=C_FIT,
            markeredgecolor="white", markeredgewidth=1.5, zorder=3)
    next_gain = 100 * (stds[order] - stds[order + 1]) / stds[order]
    ax.annotate(f"chosen: {order}\nnext gain: {next_gain:.1f}%",
                (order, stds[order]), textcoords="offset points",
                xytext=(10, 9), fontsize=8.5, color=C_FIT,
                fontweight="semibold")
    ax.set_xticks(orders)
    ax.set_title("Polynomial order", loc="left")
    style(ax, xlabel="Order", ylabel="Residual σ [nm]")

    fig.suptitle("Idle-tail detrending", x=0.075, ha="left", y=0.965,
                 fontsize=16, fontweight="semibold")
    fig.text(0.075, 0.915,
             "First 60 s after the opening transient · cubic is the first-knee choice",
             ha="left", fontsize=9.5, color=C_MUTED)
    fig.text(0.075, 0.035, f"Subtracted trend:  {expr}", ha="left",
             fontsize=8.5, color=C_INK2, family="DejaVu Sans Mono")
    fig.savefig(out_dir / "detrend_overlay.png")
    plt.close(fig)


def plot_distribution(seg_t, resid, out_dir):
    sd, ptp = float(np.std(resid)), float(np.ptp(resid))
    mean = float(np.mean(resid))
    skew = float(stats.skew(resid))
    kurt = float(stats.kurtosis(resid))
    within = 100.0 * np.mean(np.abs(resid - mean) <= sd)

    fig = plt.figure(figsize=(13.2, 5.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.15, 1.0], wspace=0.25,
                          left=0.07, right=0.975, top=0.80, bottom=0.15)

    ax = fig.add_subplot(gs[0])
    ax.axhspan(-sd, sd, color=C_RESID, alpha=0.10, zorder=0,
               label=f"±1σ  ({sd:.1f} nm)")
    for level in (-3 * sd, 3 * sd):
        ax.axhline(level, color=C_MUTED, linewidth=0.9,
                   linestyle=(0, (3, 3)), zorder=1)
    ax.axhline(0.0, color=C_BASELINE, linewidth=1.0, zorder=1)
    ax.plot(seg_t, resid, color=C_RESID, linewidth=0.55, zorder=2)
    ax.set_title("Residual over time", loc="left")
    ax.legend(loc="upper left")
    ax.text(0.985, 0.05, f"peak-to-peak  {ptp:.1f} nm   ·   dashed  ±3σ",
            transform=ax.transAxes, fontsize=8.7, color=C_INK2,
            ha="right", va="bottom", bbox=LABEL_BBOX)
    style(ax, xlabel="Time from segment start [s]", ylabel="Residual [nm]")
    ax.set_xlim(0, seg_t[-1])

    ax = fig.add_subplot(gs[1])
    bins = np.histogram_bin_edges(resid, bins="fd")
    ax.hist(resid, bins=bins, color=C_RESID, alpha=0.78, edgecolor="white",
            linewidth=0.35, density=True, zorder=2)
    xs = np.linspace(resid.min(), resid.max(), 512)
    ax.plot(xs, stats.norm.pdf(xs, mean, sd), color=C_INK,
            linewidth=1.8, linestyle=(0, (5, 2)), zorder=3,
            label="Normal reference")
    ax.axvline(mean, color=C_BASELINE, linewidth=1.0, zorder=1)
    ax.set_title("Distribution", loc="left")
    ax.legend(loc="upper right")
    ax.text(0.04, 0.05,
            f"σ  {sd:.2f} nm\n"
            f"skew  {skew:+.3f}\n"
            f"excess kurtosis  {kurt:+.3f}\n"
            f"inside ±1σ  {within:.1f}%",
            transform=ax.transAxes, fontsize=8.8, color=C_INK2,
            ha="left", va="bottom", linespacing=1.45, bbox=LABEL_BBOX)
    style(ax, xlabel="Residual [nm]", ylabel="Probability density")

    fig.suptitle("Residual spread and distribution", x=0.07, ha="left",
                 y=0.955, fontsize=16, fontweight="semibold")
    fig.text(0.07, 0.885,
             "Same order-3 detrended series in both panels · shading shows the central ±1σ band",
             ha="left", fontsize=9.5, color=C_MUTED)
    fig.savefig(out_dir / "residual_distribution.png")
    plt.close(fig)
    return dict(std_nm=sd, peak_to_peak_nm=ptp, mean_nm=mean, skew=skew,
                excess_kurtosis=kurt, within_1sigma_pct=float(within),
                n_bins=int(len(bins) - 1))


def plot_psd(f, pxx, resid, fs, n_segments, out_dir):
    integrated = float(np.sqrt(np.trapezoid(pxx, f)))
    cumulative = np.sqrt(np.cumsum(pxx) * (f[1] - f[0]))
    half = float(np.interp(0.5 * cumulative[-1], cumulative, f))

    fig = plt.figure(figsize=(13.2, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.38, 1.0], hspace=0.38,
                          left=0.078, right=0.975, top=0.84, bottom=0.09)

    ax = fig.add_subplot(gs[0])
    ax.axvspan(330, 370, color=C_SOFT_ORANGE, zorder=0,
               label="330–370 Hz mode band")
    ax.loglog(f[1:], pxx[1:], color=C_RESID, linewidth=0.9, zorder=3,
              label="Welch PSD")
    ax.axvline(362.8, color=C_FIT, linewidth=1.35, linestyle=(0, (4, 3)),
               zorder=4, label="v2 peak: 362.8 Hz")
    ax.set_title("Power spectral density", loc="left")
    ax.legend(loc="lower left", ncol=3)
    style(ax, ylabel="PSD [nm²/Hz]")
    ax.set_xlim(f[1], fs / 2)
    ax.grid(True, which="major", axis="both")
    ax.grid(True, which="minor", axis="both", linewidth=0.45, alpha=0.42)

    ax = fig.add_subplot(gs[1])
    ax.semilogx(f[1:], cumulative[1:], color=C_RAW, linewidth=1.8, zorder=3,
                label="Cumulative RMS")
    ax.axhline(np.std(resid), color=C_INK2, linewidth=1.0,
               linestyle=(0, (5, 3)), zorder=2,
               label=f"Time-domain σ: {np.std(resid):.2f} nm")
    ax.axvline(half, color=C_FIT, linewidth=1.0, linestyle=(0, (3, 3)),
               zorder=2)
    ax.scatter([half], [0.5 * cumulative[-1]], s=30, color=C_FIT,
               edgecolor="white", linewidth=0.8, zorder=4)
    ax.annotate(f"50% of RMS below {half:.2f} Hz",
                (half, 0.5 * cumulative[-1]), textcoords="offset points",
                xytext=(10, -22), fontsize=9, color=C_FIT,
                fontweight="semibold", bbox=LABEL_BBOX)
    ax.set_title("Cumulative RMS — where the power accumulates", loc="left")
    ax.legend(loc="lower right", ncol=2)
    style(ax, xlabel="Frequency [Hz]", ylabel="Cumulative RMS [nm]")
    ax.set_xlim(f[1], fs / 2)
    ax.grid(True, which="major", axis="both")
    ax.grid(True, which="minor", axis="x", linewidth=0.45, alpha=0.42)

    fig.suptitle("Spectrum of the detrended idle tail", x=0.078, ha="left",
                 y=0.965, fontsize=16, fontweight="semibold")
    fig.text(0.078, 0.91,
             f"Welch: {NPERSEG:,} samples · {n_segments} averages · Δf {fs / NPERSEG:.4f} Hz · "
             f"integrated RMS {integrated:.2f} nm (time-domain σ {np.std(resid):.2f} nm)",
             ha="left", fontsize=9.5, color=C_MUTED)
    fig.savefig(out_dir / "welch_psd.png")
    plt.close(fig)
    return integrated, half

def main():
    # Windows consoles default to cp1252, which cannot encode the Greek and
    # maths glyphs in the progress output below.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t, pos, meta = load_ids(CSV)
    dt_s = meta["dt_s"]
    fs = 1.0 / dt_s
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {fs:.0f} Hz")

    seg_t, seg_y, seg_meta = extract_segment(t, pos, dt_s)
    print(f"segment: t={seg_meta['segment_start_s']:.3f}.."
          f"{seg_meta['segment_end_s']:.3f} s "
          f"({seg_meta['n_samples']:,} samples), "
          f"raw sd={np.std(seg_y):.2f} nm, p-p={np.ptp(seg_y):.1f} nm")

    order, stds = select_order(seg_t, seg_y)
    poly = Polynomial.fit(seg_t, seg_y, order)
    coef = poly.convert().coef
    fit_y = poly(seg_t)
    resid = seg_y - fit_y
    expr = poly_expression(coef)
    print(f"order {order} selected; {expr}")

    plot_detrend_overlay(seg_t, seg_y, fit_y, resid, order, stds, expr, OUT_DIR)
    dist = plot_distribution(seg_t, resid, OUT_DIR)

    noverlap = int(NPERSEG * OVERLAP_FRAC)
    f, pxx = signal.welch(resid, fs=fs, window=WINDOW, nperseg=NPERSEG,
                          noverlap=noverlap, detrend="constant",
                          return_onesided=True, scaling="density")
    n_segments = 1 + (resid.size - NPERSEG) // (NPERSEG - noverlap)
    integrated, half_f = plot_psd(f, pxx, resid, fs, n_segments, OUT_DIR)
    print(f"Welch: {n_segments} segments, Δf={fs / NPERSEG:.4f} Hz; "
          f"integrated {integrated:.2f} nm vs sd {np.std(resid):.2f} nm; "
          f"half the RMS below {half_f:.2f} Hz")

    df = float(f[1] - f[0])

    def band_rms(lo_hz, hi_hz):
        m = (f >= lo_hz) & (f < hi_hz)
        return float(np.sqrt(np.sum(pxx[m]) * df))

    bands = {f"{lo:g}-{hi:g}_Hz": band_rms(lo, hi) for lo, hi in
             ((0, 1), (1, 10), (10, 100), (100, 300), (300, 400), (400, 500))}
    floor_m = (f >= 10) & (f < 250)
    # The global argmax is the 0.12 Hz wander bin, which is a trend artifact
    # rather than a mechanical feature; the peak that matters for the
    # mechanics is the largest one ABOVE the wander, reported separately.
    hf = f > 10.0
    peak = int(np.flatnonzero(hf)[int(np.argmax(pxx[hf]))])
    wander_peak = int(np.argmax(pxx[1:])) + 1
    summary = {
        "source": str(CSV),
        "segment": {
            **seg_meta,
            "skip_after_tail_boundary_s": SKIP_S,
            "skip_rationale": "excludes the ~1.15 um unexplained opening "
                              "transient; peak-to-peak stops changing past "
                              "50 ms (1271 -> 275 -> 219 -> 107 nm) and the "
                              "last >100 nm step-to-step jump is at 26 ms",
            "duration_s": SEGMENT_S,
            "sample_rate_hz": fs,
            "raw_std_nm": float(np.std(seg_y)),
            "raw_peak_to_peak_nm": float(np.ptp(seg_y)),
        },
        "detrend": {
            "order": order,
            "expression_nm_vs_s": expr,
            "coefficients_ascending": [float(c) for c in coef],
            "selection_rule": f"first order whose successor improves residual "
                              f"std by < {KNEE_IMPROVEMENT:.0%}",
            "residual_std_by_order_nm": stds,
            "std_reduction_pct": float(100 * (1 - np.std(resid) / np.std(seg_y))),
        },
        "residual": dist,
        "welch_psd": {
            "window": WINDOW,
            "nperseg": NPERSEG,
            "segment_length_s": NPERSEG / fs,
            "overlap_fraction": OVERLAP_FRAC,
            "n_segments_averaged": int(n_segments),
            "frequency_resolution_hz": float(fs / NPERSEG),
            "nyquist_hz": float(fs / 2),
            "detrend_per_segment": "constant",
            "integrated_rms_nm": integrated,
            "parseval_ratio": float(integrated / np.std(resid)),
            "peak_frequency_above_10hz_hz": float(f[peak]),
            "peak_psd_above_10hz_nm2_per_hz": float(pxx[peak]),
            "wander_peak_frequency_hz": float(f[wander_peak]),
            "wander_peak_psd_nm2_per_hz": float(pxx[wander_peak]),
            "first_bin_rms_nm": float(np.sqrt(pxx[1] * df)),
            "band_rms_nm": bands,
            "white_floor_10_250hz": {
                "median_psd_nm2_per_hz": float(np.median(pxx[floor_m])),
                "asd_nm_per_sqrt_hz": float(np.sqrt(np.median(pxx[floor_m]))),
            },
            "half_rms_below_hz": float(half_f),
        },
        "outputs": {
            "detrend_overlay": "rendered_assets/detrend_overlay.png",
            "distribution": "rendered_assets/residual_distribution.png",
            "psd": "rendered_assets/welch_psd.png",
        },
    }
    (V22_DIR / "tail_60s_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"wrote 3 figures -> {OUT_DIR}")
    print(f"wrote tail_60s_summary.json -> {V22_DIR}")


if __name__ == "__main__":
    main()

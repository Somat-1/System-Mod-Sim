#!/usr/bin/env python3
"""E0 — the ISO 230-2 v3 baseline, drawn.

Step A produced the numbers; this draws them. Figure 1 is the argument in one
view: the standard ISO deviation chart with the fitted gain on it, the residual
that survives the gain, and the reversal. Sharing one position axis is the
point - the reader sees that the residual has structure while the reversal does
not, which is why G and friction stay separate models.

Figure 2 is the improvement chart's first bar: the ISO metrics before and after
a position-only gain correction.

Reads  ../results/per_target_statistics.csv and ../results/step_a_results.json.
Writes rendered/F1_deviation_residual_reversal.png
       rendered/F2_metric_improvement.png
"""
from pathlib import Path
import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, INK3, SURFACE, style   # noqa

RESULTS = Path(os.environ.get("E0_RESULTS", _HERE.parent / "results"))
OUT = Path(os.environ.get("E0_OUT", _HERE / "rendered"))

UP, DN = "#5f8bc5", "#c88e5c"      # the two directions
FIT, MEAN = "#99483e", "#356449"   # fitted gain, fitted trend / after-state
AVG = "#161616"                    # the direction average, which the fit uses

MINUS = "−"                   # typographic minus, not a hyphen


def num(value, places=2):
    """Format with a real minus sign, for anything a reader sees."""
    return f"{value:.{places}f}".replace("-", MINUS)


def thesis_style():
    """House style, then the page: white paper and a LaTeX-matching serif."""
    style()
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "figure.facecolor": "#ffffff",
        "axes.facecolor": "#ffffff",
        "savefig.facecolor": "#ffffff",
    })


def load():
    with (RESULTS / "per_target_statistics.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    col = lambda k: np.array([float(r[k]) for r in rows])   # noqa: E731
    d = {k: col(k) for k in
         ("cmd_mm", "xbar_up_um", "xbar_dn_um", "s_up_um", "s_dn_um",
          "ybar_um", "var_ybar_um2", "B_um", "fit_um", "residual_um")}
    d["n_up"] = np.array([float(r["n_up"]) for r in rows])
    d["n_dn"] = np.array([float(r["n_dn"]) for r in rows])
    # var(B) = s_up^2/n_up + s_dn^2/n_dn, which is 4 x var(ybar) since
    # ybar = (xbar_up + xbar_dn)/2.
    d["sigma_ybar"] = np.sqrt(d["var_ybar_um2"])
    d["sigma_B"] = 2.0 * d["sigma_ybar"]
    with (RESULTS / "step_a_results.json").open(encoding="utf-8") as fh:
        res = json.load(fh)
    return d, res


def line(fit, x):
    return fit["intercept_centered_um"] + fit["slope_um_per_mm"] * (
        x - fit["x_center_mm"])


def frame(ax, lo, hi, below=0.0, above=0.0):
    """Set the y-limits so a strip of the given height - as a fraction of the
    final axis range - stays clear below and above the data, for a legend or a
    note to sit in without covering a marker, a line or a band."""
    rng = (hi - lo) / (1.0 - below - above)
    ax.set_ylim(lo - below * rng, hi + above * rng)


def figure1(d, res, out):
    fit = res["fit_uncorrected"]
    rev = res["reversal"]["fit"]
    x = d["cmd_mm"]
    xs = np.linspace(0, 122, 200)

    thesis_style()
    fig, (ax, bx, cx) = plt.subplots(
        3, 1, figsize=(8.6, 9.2), sharex=True,
        gridspec_kw=dict(hspace=0.15, left=0.105, right=0.975,
                         top=0.965, bottom=0.062))

    # (a) the standard ISO chart, with the gain on it -----------------------
    # Markers only: the targets are ~5 mm apart and a 2 mm periodic term is
    # expected between them, so any line or band drawn across that gap asserts
    # smoothness the data cannot support.
    for y, s, c, m, lab in (
            (d["xbar_up_um"], d["s_up_um"], UP, "o",
             "forward  $\\bar{x}_i\\uparrow$"),
            (d["xbar_dn_um"], d["s_dn_um"], DN, "s",
             "reverse  $\\bar{x}_i\\downarrow$")):
        ax.fill_between(x, y - 2 * s, y + 2 * s, color=c, alpha=0.10, lw=0)
        ax.plot(x, y, m, color=c, ms=4.2, label=lab)
    ax.plot(x, d["ybar_um"], "D", color=AVG, ms=2.9,
            label="direction average $\\bar{y}_i$  (the fit's data)")
    ax.plot(xs, line(fit, xs), color=FIT, lw=1.8,
            label=f"$G_\\mathrm{{linear}}$: {num(fit['slope_ppm'], 1)} ± "
                  f"{fit['slope_ci95_ppm']:.1f} ppm (95 %)")
    ax.plot([], [], "s", color=INK3, alpha=0.35, ms=7,
            label="shaded: ±2$s$ (run-to-run)")
    ax.axhline(0, color=INK3, lw=0.7, zorder=0)
    ax.set_ylabel("deviation (µm)")
    ax.set_title("(a)  ISO 230-2 deviation with the fitted gain",
                 size=10, loc="left")
    ax.legend(frameon=False, fontsize=8.0, ncol=3, loc="lower left",
              handletextpad=0.5, columnspacing=1.3)
    frame(ax,
          float(np.min([d["xbar_up_um"] - 2 * d["s_up_um"],
                        d["xbar_dn_um"] - 2 * d["s_dn_um"], d["ybar_um"]])),
          float(np.max([d["xbar_up_um"] + 2 * d["s_up_um"],
                        d["xbar_dn_um"] + 2 * d["s_dn_um"], d["ybar_um"]])),
          below=0.21, above=0.04)

    # (b) what the gain leaves behind ---------------------------------------
    bx.axhline(0, color=INK2, lw=1.0, zorder=1)
    bx.errorbar(x, d["residual_um"], yerr=2 * d["sigma_ybar"], fmt="o",
                color=FIT, ms=4.0, elinewidth=1.1, capsize=2.6,
                mfc="white", mec=FIT, mew=1.2,
                label="$r_i$ ± 2$\\sigma$ of the mean")
    rms = float(np.sqrt(np.mean(d["residual_um"] ** 2)))
    bx.axhspan(-rms, rms, color=INK3, alpha=0.12, lw=0,
               label=f"±rms  {rms:.3f} µm")
    bx.set_ylabel("residual $r_i$ (µm)")
    bx.set_title("(b)  Residual after the gain — structured, not noise",
                 size=10, loc="left")
    bx.legend(frameon=False, fontsize=8.0, ncol=2, loc="lower left")
    frame(bx, float(np.min(d["residual_um"] - 2 * d["sigma_ybar"])),
          float(np.max(d["residual_um"] + 2 * d["sigma_ybar"])),
          below=0.14, above=0.04)

    # (c) the reversal, for the friction hand-off ---------------------------
    cx.axhline(0, color=INK3, lw=0.7, zorder=0)
    cx.errorbar(x, d["B_um"], yerr=2 * d["sigma_B"], fmt="o", color=MEAN,
                ms=4.0, elinewidth=1.1, capsize=2.6, mfc="white", mec=MEAN,
                mew=1.2, label="$B_i$ ± 2$\\sigma$ of the mean")
    cx.axhline(res["reversal"]["mean_um"], color=INK2, lw=1.2, ls=(0, (5, 3)),
               label=f"mean  {num(res['reversal']['mean_um'], 3)} µm")
    cx.plot(xs, line(rev, xs), color=FIT, lw=1.6,
            label=f"trend  {rev['slope_ppm']:+.2f} ± "
                  f"{rev['slope_ci95_ppm']:.2f} ppm")
    cx.set_ylabel("$B_i$ = forward − reverse (µm)")
    cx.set_xlabel("commanded position (mm)")
    cx.set_title("(c)  Reversal", size=10, loc="left")
    cx.legend(frameon=False, fontsize=8.0, ncol=3, loc="upper left")
    cx.annotate("$B_i$ < 0: the stage lags the command in the "
                "direction of travel.", (0.995, 0.022),
                xycoords="axes fraction", ha="right", va="bottom",
                size=8.0, color=INK2)
    frame(cx, float(np.min(d["B_um"] - 2 * d["sigma_B"])),
          float(np.max(d["B_um"] + 2 * d["sigma_B"])),
          below=0.09, above=0.13)

    cx.set_xlim(0, 122)
    for a in (ax, bx, cx):
        a.grid(True, axis="y")
        a.grid(False, axis="x")

    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "F1_deviation_residual_reversal.png")
    plt.close(fig)


def figure2(res, out):
    keys = ["E", "M", "A", "B", "R"]
    before = [res["iso_metrics_measured_um"][k] for k in keys]
    after = [res["iso_metrics_predicted_after_gain_um"][k] for k in keys]
    xs = np.arange(len(keys))
    w = 0.36

    thesis_style()
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    fig.subplots_adjust(left=0.098, right=0.975, top=0.90, bottom=0.135)
    ax.set_axisbelow(True)
    ax.bar(xs - w / 2, before, w, color=INK3, lw=0, label="measured")
    ax.bar(xs + w / 2, after, w, color=MEAN, lw=0,
           label="after gain correction")
    for i, (b, a) in enumerate(zip(before, after)):
        ax.text(i - w / 2, b + 0.13, f"{b:.2f}", ha="center", size=8.0,
                color=INK2)
        ax.text(i + w / 2, a + 0.13, f"{a:.2f}", ha="center", size=8.0,
                color=MEAN)
        if b - a > 0.05:
            ax.text(i, max(b, a) + 0.72, f"−{100*(b-a)/b:.0f}%", ha="center",
                    size=8.4, color=FIT)
    ax.set_xticks(xs)
    ax.set_xticklabels(keys, size=11)
    ax.set_ylabel("µm")
    ax.set_ylim(0, max(before) * 1.26)
    ax.set_title("ISO 230-2 metrics before and after a position-only gain "
                 "correction", size=10, loc="left")
    ax.legend(frameon=False, fontsize=8.6, loc="upper right")
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")

    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "F2_metric_improvement.png")
    plt.close(fig)


def main() -> int:
    if not (RESULTS / "step_a_results.json").is_file():
        print(f"missing Step A results under {RESULTS}")
        return 2
    d, res = load()
    figure1(d, res, OUT)
    figure2(res, OUT)
    fit = res["fit_uncorrected"]
    print(f"targets {len(d['cmd_mm'])}  span "
          f"{d['cmd_mm'].min():.1f}-{d['cmd_mm'].max():.1f} mm")
    print(f"gain {fit['slope_ppm']:.2f} +/- {fit['slope_ci95_ppm']:.2f} ppm, "
          f"chi2_nu {fit['reduced_chi2']:.1f}")
    print(f"residual rms {np.sqrt(np.mean(d['residual_um']**2)):.3f} um, "
          f"pk-pk {np.ptp(d['residual_um']):.3f} um")
    print(f"mean |2 sigma(ybar)| {2*d['sigma_ybar'].mean():.3f} um")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

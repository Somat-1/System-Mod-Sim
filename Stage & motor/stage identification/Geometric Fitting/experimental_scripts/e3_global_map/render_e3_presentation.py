#!/usr/bin/env python3
"""Paper figures for E3.

Three figures, in the form a reader needs rather than the form the analysis
produced:

  Fig 1  transfer to unseen positions - the result
  Fig 2  which harmonics reproduce, with uncertainty
  Fig 3  reversal, binned (optional; one table row may do instead)

Two choices here are deliberate and worth stating. The model is drawn as a
continuous curve on a fine x grid, never as a line through the measured
points: the targets are irregularly spaced, so joining them draws a sawtooth
that belongs to the sampling and not to the stage. And no harmonic phase is
plotted unless its amplitude is resolved, because the phase of a harmonic
consistent with zero is not defined.

Titles live in the caption, not on the axes.
"""
from __future__ import annotations

import csv
import json
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
E2 = HERE.parent / "e2_global_map"
E5 = HERE.parent / "e5_compensated_validation"
for path in (E2, E5, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import render_e2_findings as e2f                          # noqa: E402
import e2_run                                             # noqa: E402
import e3_run                                             # noqa: E402
from e5_geometry import FrozenModel                       # noqa: E402
from process_e3 import reduce_session, M_C, own_fit       # noqa: E402

OUT = HERE / "rendered" / "presentation"
SUMMARY = HERE / "analysis" / "e3_findings_summary.json"

C_E2, C_E3 = "#176d8c", "#b4473f"
C_INK, C_MUTED, C_GRID = "#2d3338", "#566069", "#dde1e4"
C_EDGE = "#eceff1"                 # weakly constrained ends
WINDOWS = ((4.0, 16.0), (74.0, 86.0))
EDGE_LOW, EDGE_HIGH = 5.0, 110.0
CHAR_INCH = 0.0455


def caption(fig, left, text, bottom=0.013, size=6.9) -> None:
    usable = fig.get_figwidth() * (1.0 - left) - 0.12
    width = max(40, int(usable / CHAR_INCH))
    lines = textwrap.wrap(" ".join(text.split()), width=width)
    fig.text(left, bottom, "\n".join(lines), fontsize=size, color=C_MUTED,
             va="bottom", linespacing=1.5)


def letter(ax, text, dx=-0.085, dy=1.06) -> None:
    ax.text(dx, dy, text, transform=ax.transAxes, fontsize=9.5, weight="bold",
            color=C_INK, ha="left", va="top")


def tidy(ax, grid=True) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(color=C_GRID, lw=0.5)


def theta_of_x(x_mm: np.ndarray) -> np.ndarray:
    """Screw phase at a position, the same map the plan uses."""
    return 2 * np.pi * np.mod(x_mm * 1600.0, 3200.0) / 3200.0


def edge_shading(ax, label: bool = False) -> None:
    lo, hi = ax.get_xlim()
    if lo < EDGE_LOW:
        ax.axvspan(lo, EDGE_LOW, color=C_EDGE, zorder=0)
    if hi > EDGE_HIGH:
        ax.axvspan(EDGE_HIGH, hi, color=C_EDGE, zorder=0)
    if label:
        ax.text(EDGE_HIGH + 0.6, ax.get_ylim()[0], " weakly\n constrained",
                fontsize=5.9, color=C_MUTED, va="bottom", linespacing=1.3)


# ------------------------------------------------------------------ fig 1

def figure_transfer(d3: dict, summary: dict, transfer: list[dict]) -> None:
    model = FrozenModel()
    x, ybar = d3["x"], d3["D2"]["ybar"]
    inside = model.in_domain(x)
    xi, yi = x[inside], ybar[inside]
    predicted_at = model.G_rel(xi, theta_of_x(xi))
    # One centring constant for both, so the curve and the dots share a datum.
    offset_model = predicted_at.mean()
    yi = yi - yi.mean()
    predicted_at = predicted_at - offset_model
    residual = yi - predicted_at

    fig = plt.figure(figsize=(7.8, 7.3))
    gs = fig.add_gridspec(5, 2, height_ratios=[0.30, 2.25, 0.80, 2.25, 0.80],
                          left=0.090, right=0.960, top=0.960, bottom=0.185,
                          hspace=0.26, wspace=0.44)

    bar = fig.add_subplot(gs[0, 0])
    bar.hlines(0, 2, 120.7, color="#9aa3a9", lw=2.0, capstyle="butt")
    for i, (lo, hi) in enumerate(WINDOWS):
        bar.hlines(0, lo, hi, color=C_E3, lw=5.0, capstyle="butt")
        bar.text((lo + hi) / 2, 0.55, "AB"[i], ha="center", va="bottom",
                 fontsize=6.6, weight="bold", color=C_E3)
    bar.set(xlim=(0, 122), ylim=(-1.1, 1.6), yticks=[])
    bar.set_xticks([2, 40, 80, 120])
    bar.tick_params(axis="x", labelsize=6.2, pad=1.5, length=2)
    bar.text(123, 0, "mm", fontsize=6.2, color=C_MUTED, va="center")
    bar.spines[:].set_visible(False)

    for i, (lo, hi) in enumerate(WINDOWS):
        row = 1 + 2 * i
        ax = fig.add_subplot(gs[row, 0])
        fine = np.linspace(lo, hi, 4000)
        curve = model.G_rel(fine, theta_of_x(fine)) - offset_model
        ax.plot(fine, curve, lw=0.9, color=C_E2, zorder=2,
                label="E2 model (continuous)")
        window = (xi >= lo) & (xi <= hi)
        ax.plot(xi[window], yi[window], ls="none", marker="o", ms=3.0,
                color=C_E3, zorder=3, label="E3 measured")
        ax.set(xlim=(lo, hi), ylabel=r"$\bar{y}$ [µm]")
        ax.tick_params(labelbottom=False)
        letter(ax, "a" if i == 0 else "")
        if i == 0:
            ax.legend(frameon=False, fontsize=6.6, ncol=1, loc="upper right",
                      handlelength=1.5, borderpad=0.1, labelspacing=0.25)
        ax.text(0.012, 0.955, f"window {'AB'[i]}", transform=ax.transAxes,
                va="top", fontsize=6.8, weight="bold", color=C_MUTED)
        tidy(ax)

        ax = fig.add_subplot(gs[row + 1, 0])
        ax.axhline(0, color="#9aa3a9", lw=0.7)
        ax.plot(xi[window], residual[window], ls="none", marker="o", ms=2.4,
                color=C_MUTED)
        ax.set(xlim=(lo, hi), ylim=(-1.9, 1.9), ylabel="resid\n[µm]")
        ax.set_yticks([-1, 0, 1])
        # The residual is offset within each window: the model's slow term and
        # E3's differ locally even though both are centred over the travel.
        offset = float(residual[window].mean())
        ax.axhline(offset, color=C_E3, lw=0.8, ls=(0, (3, 2)))
        ax.text(0.985, 0.92, f"mean {offset:+.2f} µm", transform=ax.transAxes,
                ha="right", va="top", fontsize=6.2, color=C_E3)
        if i == len(WINDOWS) - 1:
            ax.set_xlabel("position [mm]")
        tidy(ax)

    ax = fig.add_subplot(gs[1:3, 1])
    span = 1.06 * max(np.abs(predicted_at).max(), np.abs(yi).max())
    points = ax.scatter(predicted_at, yi, c=xi, cmap="viridis", s=8,
                        linewidths=0, alpha=0.9)
    ax.plot([-span, span], [-span, span], lw=0.8, color="#4b5257",
            ls=(0, (3, 2)), zorder=1)
    ax.set(xlim=(-span, span), ylim=(-span, span), aspect="equal",
           xlabel="model prediction [µm]", ylabel="E3 measured [µm]")
    bar2 = fig.colorbar(points, ax=ax, fraction=0.042, pad=0.025)
    bar2.set_label("position [mm]", fontsize=6.8, labelpad=1.5)
    bar2.ax.tick_params(labelsize=6.4)
    ax.text(0.05, 0.955, f"r = {summary['frozen_E2_prediction_on_E3']['correlation']:.3f}",
            transform=ax.transAxes, va="top", fontsize=7.4)
    letter(ax, "b")
    tidy(ax)

    ax = fig.add_subplot(gs[3:5, 1])
    names = [row["name"] for row in transfer]
    values = [row["value"] for row in transfer]
    colours = [row["colour"] for row in transfer]
    y = np.arange(len(names))[::-1]
    ax.barh(y, values, height=0.60, color=colours)
    floor = summary["fits"]["E3_M_C"]["rms_um"]
    ax.axvline(floor, color=C_INK, lw=1.0, ls=(0, (3, 2)), zorder=3)
    ax.text(floor + 0.03, y[0] + 0.62, f"same-day limit\n{floor:.3f} µm",
            fontsize=6.5, color=C_INK, va="top", linespacing=1.35)
    for yi_, value in zip(y, values):
        ax.text(value - 0.035, yi_, f"{value:.3f}", ha="right", va="center",
                fontsize=6.6, color="white", weight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.6)
    ax.set(xlabel="E3 residual [µm rms]", xlim=(0, max(values) * 1.12),
           ylim=(y.min() - 0.7, y.max() + 1.25))
    letter(ax, "c")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=C_GRID, lw=0.5)

    f = summary["frozen_E2_prediction_on_E3"]
    caption(fig, 0.095,
            f"E2's geometric model, frozen before this session was captured, predicting "
            f"{f['targets_in_domain']} positions drawn from an independent seed. (a) Two 12 mm "
            f"windows; the curve is the model evaluated continuously, the dots are measured means, "
            f"and the strip beneath each window is the residual. (b) Every in-domain target, "
            f"coloured by position; the largest departures sit at the ends of the travel. (c) The "
            f"model takes the error from {f['measured_rms_um']:.3f} µm rms to "
            f"{f['residual_rms_um']:.3f} µm, a factor of {f['reduction_factor']:.2f}. The "
            f"stationary model M_B reaches only {transfer[1]['value']:.3f} µm, so the "
            f"position-varying coefficients carry real structure; truncating to two harmonics "
            f"costs nothing. The dashed line is the same-day repeatability limit — the residual of "
            f"E3's own fit to its own data. The residual carries a local offset in each window: "
            f"the model's slow term and E3's differ by a few tenths of a micron over a 12 mm span, "
            f"which the full-travel view hides.")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e3_fig1_transfer.png", dpi=320)
    plt.close(fig)


# ------------------------------------------------------------------ fig 2

def figure_harmonics(model2: dict, model3: dict, grid: np.ndarray,
                     summary: dict) -> None:
    fig = plt.figure(figsize=(7.45, 6.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[2.2, 2.2, 1.55],
                          left=0.095, right=0.975, top=0.945, bottom=0.185,
                          hspace=0.46, wspace=0.28)

    panels = "abcd"
    for row, k in enumerate((1, 2)):
        a = e2f.harmonic_profile(model2, grid, k)
        b = e2f.harmonic_profile(model3, grid, k)

        ax = fig.add_subplot(gs[row, 0])
        for profile, colour, label in ((a, C_E2, "E2"), (b, C_E3, "E3")):
            ax.fill_between(grid, profile["amplitude"] - 2 * profile["amplitude_sigma"],
                            profile["amplitude"] + 2 * profile["amplitude_sigma"],
                            color=colour, alpha=0.20, lw=0)
            ax.plot(grid, profile["amplitude"], lw=1.2, color=colour, label=label)
        ax.set(xlim=(grid[0], grid[-1]), xlabel="position [mm]",
               ylabel=f"$A_{k}$ [µm]")
        edge_shading(ax, label=(row == 0))
        letter(ax, panels[2 * row])
        if row == 0:
            ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper right",
                      bbox_to_anchor=(0.80, 1.02), borderpad=0.1,
                      handlelength=1.5)
        ax.text(0.015, 0.955, f"k = {k}", transform=ax.transAxes, va="top",
                fontsize=7.4, weight="bold", color=C_INK)
        tidy(ax)

        ax = fig.add_subplot(gs[row, 1])
        reference = e2f.phase_reference(a)
        for profile, colour, label in ((a, C_E2, "E2"), (b, C_E3, "E3")):
            phase = np.degrees([e2f.nearest_unwrap(p, reference)
                                for p in profile["phase"]])
            sigma = np.degrees(profile["phase_sigma"])
            ax.fill_between(grid, phase - 2 * sigma, phase + 2 * sigma,
                            color=colour, alpha=0.20, lw=0)
            ax.plot(grid, phase, lw=1.2, color=colour, label=label)
        ax.set(xlim=(grid[0], grid[-1]), xlabel="position [mm]",
               ylabel=rf"$\phi_{k}$ [deg]")
        edge_shading(ax)
        letter(ax, panels[2 * row + 1])
        ax.text(0.015, 0.955, f"k = {k}", transform=ax.transAxes, va="top",
                fontsize=7.4, weight="bold", color=C_INK)
        tidy(ax)

    ax = fig.add_subplot(gs[2, :])
    rows = []
    for k in (3, 4):
        for model, colour, label in ((model2, C_E2, "E2"), (model3, C_E3, "E3")):
            profile = e2f.harmonic_profile(model, grid, k)
            rows.append((f"k = {k}\n{label}", float(np.mean(profile["amplitude"])),
                         2 * float(np.mean(profile["amplitude_sigma"])), colour))
    positions = np.arange(len(rows))
    ax.bar(positions, [row[1] for row in rows], 0.58,
           yerr=[row[2] for row in rows], color=[row[3] for row in rows],
           error_kw=dict(lw=0.9, ecolor=C_INK, capsize=2.5))
    ax.set_xticks(positions)
    ax.set_xticklabels([row[0] for row in rows], fontsize=6.8)
    ax.set_ylabel("mean amplitude [µm]")
    ax.set_ylim(0, max(row[1] + row[2] for row in rows) * 1.25)
    letter(ax, "e", dx=-0.050, dy=1.14)
    ax.text(0.985, 0.94, "error bars ±2σ", transform=ax.transAxes, ha="right",
            va="top", fontsize=6.6, color=C_MUTED)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=C_GRID, lw=0.5)

    h = summary["harmonic_comparison"]
    sel = summary.get("own_selection", {}).get("M_C", {})
    caption(fig, 0.095,
            f"Screw-harmonic amplitude and phase along the travel, both sessions, with ±2σ bands "
            f"from the fitted covariance. (a–d) The first two harmonics reproduce: amplitudes agree "
            f"to {abs(h[0]['amplitude_ratio'] - 1) * 100:.1f}% and "
            f"{abs(h[1]['amplitude_ratio'] - 1) * 100:.1f}%, and the phase curves overlap within "
            f"their bands across the interior. The shaded ends, below {EDGE_LOW:.0f} mm and above "
            f"{EDGE_HIGH:.0f} mm, are weakly constrained by the spline and are where k = 1 diverges. "
            f"(e) Harmonics 3 and 4 are resolved in E2 — {h[2]['E2_mean_amplitude_um']:.3f} and "
            f"{h[3]['E2_mean_amplitude_um']:.3f} µm, both above zero at 2σ — and absent in E3, "
            f"where {h[2]['E3_mean_amplitude_um']:.3f} and {h[3]['E3_mean_amplitude_um']:.3f} µm "
            f"are consistent with zero. That is the failure of reproducibility: present one day, "
            f"gone the next. Their phases are undefined wherever the amplitude is not resolved and "
            f"are not plotted. E3's own model selection, run on this capture alone, chooses "
            f"H = {sel.get('H', 2)} where E2 chose H = 4.")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e3_fig2_harmonics.png", dpi=320)
    plt.close(fig)


# ------------------------------------------------------------------ fig 3

def figure_reversal(d2: dict, d3: dict, summary: dict) -> None:
    fig = plt.figure(figsize=(7.45, 3.6))
    gs = fig.add_gridspec(1, 1, left=0.095, right=0.975, top=0.945,
                          bottom=0.315)
    ax = fig.add_subplot(gs[0, 0])

    edges = np.arange(2.0, 122.0, 8.0)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for data, colour, label in ((d2, C_E2, "E2"), (d3, C_E3, "E3")):
        means, errors, keep = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (data["x"] >= lo) & (data["x"] < hi)
            if sel.sum() < 3:
                continue
            values = data["D2"]["B"][sel]
            keep.append(0.5 * (lo + hi))
            means.append(values.mean())
            errors.append(2 * values.std(ddof=1) / np.sqrt(values.size))
        ax.errorbar(keep, means, yerr=errors, marker="o", ms=3.6, lw=1.0,
                    color=colour, capsize=2.5, elinewidth=0.9, label=label)
        ax.axhline(summary["reversal"][label]["beta0_um"], color=colour, lw=0.9,
                   ls=(0, (3, 2)))
    ax.set(xlabel="position [mm]", ylabel=r"$B = e_\uparrow - e_\downarrow$ [µm]")
    ax.margins(y=0.20)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower left",
              borderpad=0.2, handlelength=1.6)
    ax.text(0.985, 0.055, "points: 8 mm bin means ±2σ      dashed: fitted β₀",
            transform=ax.transAxes, ha="right", fontsize=6.6, color=C_MUTED)
    tidy(ax)

    r2, r3 = summary["reversal"]["E2"], summary["reversal"]["E3"]
    delta = r3["beta0_um"] - r2["beta0_um"]
    sigma = np.hypot(r2["beta0_sigma_um"], r3["beta0_sigma_um"])
    caption(fig, 0.095,
            f"Reversal, binned in 8 mm steps so the per-target scatter does not obscure the trend. "
            f"The constant part moved from {r2['beta0_um']:.3f} ± {r2['beta0_sigma_um']:.3f} µm to "
            f"{r3['beta0_um']:.3f} ± {r3['beta0_sigma_um']:.3f} µm, a change of {abs(delta):.3f} µm "
            f"or {abs(delta) / sigma:.1f}σ — small, but far outside the fitting uncertainty. The "
            f"position dependence also changed sign, from {r2['slope_um_per_100mm']:+.3f} to "
            f"{r3['slope_um_per_100mm']:+.3f} µm per 100 mm, so this is a change in how reversal "
            f"varies along the screw and not only in its level.")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e3_fig3_reversal.png", dpi=320)
    plt.close(fig)


# -------------------------------------------------------------- appendix

def figure_independence(d2: dict, d3: dict) -> None:
    fig = plt.figure(figsize=(3.6, 3.0))
    ax = fig.add_axes([0.185, 0.275, 0.785, 0.70])
    bins = np.arange(0, 0.62, 0.04)
    for data, colour, label in ((d2, C_E2, "E2"), (d3, C_E3, "E3")):
        ax.hist(np.diff(data["usteps"]) / 1600.0, bins=bins, histtype="step",
                lw=1.2, color=colour, label=label)
    ax.set(xlabel="step between consecutive targets [mm]", ylabel="count")
    ax.legend(frameon=False, fontsize=7)
    tidy(ax)
    caption(fig, 0.06,
            "The two target sets are drawn from the same stratified design but "
            "independent seeds; no E3 target coincides with an E2 one.",
            bottom=0.025, size=6.2)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e3_appendix_independence.png", dpi=320)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    e2f.style()
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))

    d2 = reduce_session(e2_run, "E2")
    d3 = reduce_session(e3_run, "E3")

    f = summary["frozen_E2_prediction_on_E3"]
    # The names carry commas, so this file is read as CSV rather than split.
    with (HERE / "analysis" / "e3_transfer_table.csv").open(
            newline="", encoding="utf-8") as fh:
        table = {row["name"]: float(row["e3_residual_um"])
                 for row in csv.DictReader(fh)}
    transfer = [
        {"name": "uncorrected", "value": f["measured_rms_um"], "colour": "#6b7880"},
        {"name": "M_B  stationary", "value": table["M_B  H2, stationary"],
         "colour": "#9aa3a9"},
        {"name": "M_C  H = 4 (selected)", "value": table["M_C  H4, varying (selected)"],
         "colour": C_E2},
        {"name": "H = 2, varying", "value": table["H2, varying"], "colour": C_E2},
    ]

    model2 = own_fit(d2, M_C)
    model3 = own_fit(d3, M_C)
    grid = np.linspace(max(d2["x"].min(), d3["x"].min()),
                       min(d2["x"].max(), d3["x"].max()), 500)

    figure_transfer(d3, summary, transfer)
    figure_harmonics(model2, model3, grid, summary)
    figure_reversal(d2, d3, summary)
    figure_independence(d2, d3)
    print(f"wrote 4 figures to {OUT.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Render the geometry fit that every later block subtracts.

The fit is a constant plus the 1/rev and 2/rev mechanical orders, least-squares
fitted to the *direction-pair mean* of the static map (actions 15/16) so that
friction hysteresis cannot bias it. No linear term is included, so nothing but
true geometry orders is removed.

Once fitted it is never refitted. The current sweep, dynamics and velocity
blocks all evaluate this same curve at absolute shaft angle. The right-hand
panel shows why that matters: over the 28.8 deg the current sweep spans, the
curve is very nearly a straight line, so refitting it there would absorb a
local offset and slope that belong to the motor, not to the encoder.

Writes analysis/staircase/geometry_fit.png.
"""
from pathlib import Path
import json
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from staircase_figures import (                                    # noqa: E402
    CACHE, ROOT, DPC, APC, C2C, REV_U256, ELEC_DEG, INK, INK2, INK3,
    SURFACE, GRID, BAND, style, design, harmonic_fit, paired, load_static,
)
from current_sweep_figures import titleblock3                      # noqa: E402

OUT = Path(os.environ.get("STAIRCASE_OUT", ROOT / "analysis" / "staircase"))

ONE = "#2a78d6"          # categorical slot 1 - the 1/rev term
TWO = "#eb6834"          # slot 2 - the 2/rev term
CURRENT_SWEEP_DEG = 28.8  # 16 full steps, the span of each current level
V1_ECCENTRICITY_UM = 159.52


def fit():
    d = load_static(CACHE / "static.csv")
    p = paired(d, 15, 16)
    ideal = p["cmd"] * C2C
    mean = ((p["ef"] - ideal) + (p["er"] - ideal)) / 2.0
    angle = p["cmd"] / REV_U256 * 360.0
    coef, names, _, geo = harmonic_fit(angle, mean, (1, 2), linear=False)
    return angle, mean, coef, names, geo


def component(angle, coef, names, order):
    A, nm = design(angle, (1, 2), linear=False)
    out = np.zeros_like(np.asarray(angle, dtype=float))
    for part in ("cos", "sin"):
        k = nm.index(f"{part}{order}")
        out += A[:, k] * coef[names.index(f"{part}{order}")]
    return out


def main():
    if not (CACHE / "static.csv").exists():
        raise SystemExit(f"missing {CACHE / 'static.csv'} - see "
                         "staircase_figures.py for the extraction command")
    OUT.mkdir(parents=True, exist_ok=True)
    style()

    angle, mean, coef, names, geo = fit()
    g = np.linspace(0, 360, 2001)
    A, _ = design(g, (1, 2), linear=False)
    total = A @ coef
    c1 = component(g, coef, names, 1)
    c2 = component(g, coef, names, 2)
    const = float(coef[names.index("const")])

    fitted_at_data = design(angle, (1, 2), linear=False)[0] @ coef
    resid = mean - fitted_at_data
    var_expl = 100 * (1 - np.var(resid) / np.var(mean))

    # Behaviour over the window the current sweep actually spans.
    w = np.linspace(0, CURRENT_SWEEP_DEG, 400)
    wfit = design(w, (1, 2), linear=False)[0] @ coef
    line = np.polyval(np.polyfit(w, wfit, 1), w)
    max_dev = float(np.max(np.abs(wfit - line)))

    fig = plt.figure(figsize=(11.0, 8.0))
    ax = fig.add_axes([0.075, 0.435, 0.845, 0.275])
    bx = fig.add_axes([0.075, 0.085, 0.385, 0.215])
    cx = fig.add_axes([0.535, 0.085, 0.385, 0.215])

    # --- the fit over one revolution -----------------------------------
    ax.grid(False)
    ax.grid(True, axis="y")
    ax.axhline(0, color=INK3, lw=0.8, zorder=2)
    ax.axvspan(0, CURRENT_SWEEP_DEG, color=BAND, alpha=0.20, lw=0, zorder=1)
    ax.plot(g, (total - const) * DPC, color=INK, lw=2.6, zorder=5,
            label=f"geometry fit — amplitude "
                  f"{geo[1]['amp'] * DPC:.4f}° at phase "
                  f"{geo[1]['phase_deg']:.0f}°")
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.set_xlabel("shaft angle (degrees)")
    ax.set_ylabel("correction (degrees)")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=1,
              handlelength=2.0, borderaxespad=0.0, fontsize=8.6)
    sec = ax.secondary_yaxis("right", functions=(lambda v: v / DPC,
                                                 lambda v: v * DPC))
    sec.set_ylabel("encoder counts", color=INK2)
    sec.spines["right"].set_color(GRID)
    ax.annotate(f"shaded: the {CURRENT_SWEEP_DEG:g}° each current level "
                f"spans (panel at right)",
                xy=(CURRENT_SWEEP_DEG + 6, 0.94), xycoords=("data", "axes fraction"),
                color=INK2, size=8.4, va="top")
    ax.annotate(f"the 1/rev term alone is this curve to within "
                f"{geo[2]['amp'] / geo[1]['amp'] * 100:.1f}%;\n"
                "the 2/rev remainder is below left",
                xy=(0.985, 0.06), xycoords="axes fraction", ha="right",
                va="bottom", color=INK2, size=8.4)

    # --- the 2/rev term on its own scale --------------------------------
    bx.grid(False)
    bx.grid(True, axis="y")
    bx.axhline(0, color=INK3, lw=0.8)
    bx.plot(g, c2 * APC, color=TWO, lw=2.0)
    bx.set_xlim(0, 360)
    bx.set_xticks(np.arange(0, 361, 90))
    bx.set_xlabel("shaft angle (degrees)")
    bx.set_ylabel("2/rev term (arcsec)")
    bx.set_title(f"2/rev is {geo[2]['amp'] / geo[1]['amp'] * 100:.1f}% of 1/rev "
                 "- the correction is effectively pure eccentricity",
                 size=9.2, loc="left", pad=8)

    # --- the fit across the current-sweep window -------------------------
    cx.grid(False)
    cx.grid(True, axis="y")
    cx.plot(w, (wfit - wfit[0]) * DPC, color=INK, lw=2.4, label="geometry fit")
    cx.plot(w, (line - wfit[0]) * DPC, color=TWO, lw=1.4, ls=(0, (5, 3)),
            label="best straight line")
    cx.set_xlim(0, CURRENT_SWEEP_DEG)
    cx.set_xlabel("shaft angle (degrees)")
    cx.set_ylabel("correction (degrees)")
    cx.legend(fontsize=8.0, loc="lower left")
    cx.set_title(f"over {CURRENT_SWEEP_DEG:g}° it departs from a straight "
                 f"line by only {max_dev:.2f} counts",
                 size=9.2, loc="left", pad=8)

    titleblock3(
        fig,
        "The geometry fit subtracted by every block",
        f"Constant + 1/rev + 2/rev, fitted once to the direction-pair mean of "
        f"the static map ({len(angle):,} commanded positions) and then "
        f"evaluated at absolute shaft angle wherever it is needed. No linear "
        f"term, so only true geometry orders are removed.",
        f"Amplitude {geo[1]['amp'] * DPC:.4f}° "
        f"({geo[1]['amp']:.2f} counts) at phase "
        f"{geo[1]['phase_deg']:.1f}°, accounting for "
        f"{var_expl:.1f}% of the variance in the direction-pair mean. The "
        f"2/rev term is negligible, so the correction is a single sinusoid: "
        f"an eccentric encoder disc.", width=118)

    fig.savefig(OUT / "geometry_fit.png")
    plt.close(fig)

    metrics = {
        "source": "static map, actions 15/16, direction-pair mean",
        "model": "const + cos/sin at 1 and 2 cycles per revolution; "
                 "no linear term; never refitted",
        "positions_fitted": int(len(angle)),
        "const_counts": const,
        "coefficients_counts": {n: float(c) for n, c in zip(names, coef)},
        "one_per_rev": {"amp_counts": geo[1]["amp"],
                        "amp_deg": geo[1]["amp"] * DPC,
                        "amp_arcsec": geo[1]["amp"] * APC,
                        "phase_deg": geo[1]["phase_deg"]},
        "two_per_rev": {"amp_counts": geo[2]["amp"],
                        "amp_deg": geo[2]["amp"] * DPC,
                        "amp_arcsec": geo[2]["amp"] * APC,
                        "phase_deg": geo[2]["phase_deg"],
                        "percent_of_one_per_rev":
                            geo[2]["amp"] / geo[1]["amp"] * 100},
        "variance_explained_pct": float(var_expl),
        "residual_rms_counts": float(np.sqrt(np.mean(resid ** 2))),
        "current_sweep_window": {
            "span_deg": CURRENT_SWEEP_DEG,
            "fit_change_over_window_counts": float(wfit[-1] - wfit[0]),
            "max_departure_from_straight_line_counts": max_dev},
        "implied_scale_radius_mm":
            V1_ECCENTRICITY_UM / (geo[1]["amp"] * DPC * np.pi / 180) / 1000,
    }
    (OUT / "geometry_fit.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    np.savetxt(OUT / "geometry_fit_curve.csv",
               np.column_stack([g, (total - const) * DPC, c1 * DPC, c2 * DPC]),
               delimiter=",", comments="",
               header="shaft_deg,total_deg,one_per_rev_deg,two_per_rev_deg",
               fmt="%.6f")

    print("coefficients (encoder counts):")
    for n, c in zip(names, coef):
        print(f"  {n:<8} {c:+10.4f}")
    print(f"1/rev      {geo[1]['amp']:.3f} counts = {geo[1]['amp'] * DPC:.4f} deg "
          f"= {geo[1]['amp'] * APC:.1f} arcsec at {geo[1]['phase_deg']:.1f} deg")
    print(f"2/rev      {geo[2]['amp']:.3f} counts = {geo[2]['amp'] * APC:.1f} arcsec "
          f"({geo[2]['amp'] / geo[1]['amp'] * 100:.1f}% of 1/rev)")
    print(f"variance   {var_expl:.2f}% explained, residual rms "
          f"{metrics['residual_rms_counts']:.3f} counts")
    print(f"over {CURRENT_SWEEP_DEG:g} deg: changes "
          f"{metrics['current_sweep_window']['fit_change_over_window_counts']:+.2f} "
          f"counts, departs from a straight line by {max_dev:.3f} counts")
    print(f"implied scale radius {metrics['implied_scale_radius_mm']:.1f} mm")
    print(f"wrote {OUT / 'geometry_fit.png'}")


if __name__ == "__main__":
    main()

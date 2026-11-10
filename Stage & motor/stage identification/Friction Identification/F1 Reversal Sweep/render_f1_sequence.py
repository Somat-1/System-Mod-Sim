#!/usr/bin/env python3
"""Render the paper-width F1 reversal-sweep sequence overview (procedure §8)."""
from __future__ import annotations

import json
import re
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from f1_sequence import CONFIG, plan, strokes, summarise

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "rendered" / "f1_sequence.png"
GEOM = HERE.parents[1] / "Geometric Fitting" / "experimental_scripts" / "e5_compensated_validation"

INK = "#20262b"
MUTED = "#59636b"
GRID = "#dfe3e6"
GROUP_SHADE = {"A": "#eef0f2", "B": "#e3e7ea", "C": "#eef0f2"}


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.edgecolor": "#7d858c",
        "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c",
        "ytick.color": "#41474c",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def tidy(ax, grid: bool = True) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(color=GRID, lw=0.5, zorder=0)


def spans(segs, key):
    out = []
    for s in segs:
        k = key(s)
        if out and out[-1][0] == k:
            out[-1][2] = s.t1
        else:
            out.append([k, s.t0, s.t1])
    return [(k, a / 60, b / 60) for k, a, b in out]


def add_ribbon(ax, segs, c, total):
    """Two tiers: position groups on top, sub-blocks below; widths are time."""
    cols = c["colors"]
    top = ax.inset_axes([0, 1.10, 1, .065])
    low = ax.inset_axes([0, 1.025, 1, .065])
    pos = c["positions_mm"]
    for g, a, b in spans(segs, lambda s: s.group):
        is_ref = g.startswith("R")
        top.barh(0, b - a, left=a, height=.9, color=cols["R"] if is_ref else "#4a545c",
                 edgecolor="white", lw=.6)
        label = g if is_ref else f"{g} · {pos[g]:g} mm"
        top.text((a + b) / 2, 0, label, ha="center", va="center", color="white",
                 weight="bold", fontsize=7)
    for (g, sub), a, b in spans(segs, lambda s: (s.group, s.sub)):
        colour = cols["R"] if sub == "R" else cols[sub]
        low.barh(0, b - a, left=a, height=.9, color=colour, edgecolor="white", lw=.6)
        if sub != "R" and b - a > .55:
            low.text((a + b) / 2, 0, sub, ha="center", va="center", color="white",
                     weight="bold", fontsize=7)
    for r in (top, low):
        r.set(xlim=(0, total), ylim=(-.5, .5), xticks=[], yticks=[])
        r.spines[:].set_visible(False)


def panel_trajectory(ax, segs, c, total):
    cols = c["colors"]
    setup, transit, ref, measured = [], [], [], {}
    for s in segs:
        line = [(s.t0 / 60, s.u0 - s.ref_um), (s.t1 / 60, s.u1 - s.ref_um)]
        if s.kind == "setup":
            setup.append(line)
        elif s.kind == "transit":
            transit.append(line)
        elif s.kind == "ref":
            ref.append(line)
        else:
            measured.setdefault(s.sub, []).append(line)
    ax.add_collection(LineCollection(setup, colors=cols["setup"], linewidths=.65, zorder=1))
    ax.add_collection(LineCollection(transit, colors=cols["setup"], linewidths=.9,
                                     linestyles=(0, (3, 2)), zorder=1))
    ax.add_collection(LineCollection(ref, colors=cols["R"], linewidths=1.0,
                                     linestyles=(0, (4, 2)), zorder=3))
    for sub, lines in measured.items():
        thin = sub in ("P", "Q")
        ax.add_collection(LineCollection(lines, colors=cols[sub], linewidths=.35 if thin else 1.35,
                                         alpha=.35 if thin else 1, zorder=3))

    # P and Q read as oscillation: light envelope fill with a thin outline.
    for sub in ("P", "Q"):
        loops = {}
        for s in segs:
            if s.sub == sub and s.kind == "measure":
                key = (s.group, re.sub(r" (to start|c\d+ (up|down)|return)( steps)?$", "", s.tag))
                a = loops.setdefault(key, [s.t0, s.t1, 0.0])
                a[1] = s.t1
                a[2] = max(a[2], abs(s.u1 - s.ref_um))
        # One envelope per group: each loop holds its half-stroke across its span.
        for g in {k[0] for k in loops}:
            rows = sorted(v for k, v in loops.items() if k[0] == g)
            env_t = np.array([[a, b] for a, b, _ in rows]).ravel() / 60
            env_a = np.repeat([h for _, _, h in rows], 2)
            for sign in (1, -1):
                ax.fill_between(env_t, 0, sign * env_a, color=cols[sub], alpha=.05, lw=0,
                                zorder=2)
                ax.plot(env_t, sign * env_a, color=cols[sub], lw=.55, alpha=.8, zorder=3)

    # Transit labels at the start of each new position.
    pos = c["positions_mm"]
    for g, a, _ in spans(segs, lambda s: s.group):
        if g in pos:
            ax.text(a + .25, -3.2e4, rf"$\rightarrow$ {pos[g]:g} mm", fontsize=7, color=MUTED,
                    ha="left", va="top", bbox=dict(fc="white", ec="none", pad=.4), zorder=5)

    for g, a, b in spans(segs, lambda s: s.group):
        if g in GROUP_SHADE:
            ax.axvspan(a, b, color=GROUP_SHADE[g], alpha=.35, lw=0, zorder=-2)
    ax.axhspan(-1, 1, color="#f4f5f6", zorder=-1)
    ax.set_yscale("symlog", linthresh=1.0, linscale=.55)
    ticks = [-1e5, -1e4, -1e3, -100, -10, -1, 1, 10, 100, 1e3, 1e4, 1e5]
    labels = ["−100 mm", "−10 mm", "−1 mm", "−100 µm", "−10 µm", "−1 µm",
              "+1 µm", "+10 µm", "+100 µm", "+1 mm", "+10 mm", "+100 mm"]
    ax.set_yticks(ticks, labels)
    ax.set(xlim=(0, total), ylim=(-1.6e5, 1.6e5),
           ylabel=r"commanded position relative to $x_p$", xlabel="elapsed session time [min]")
    ax.set_title("(a)", loc="left", pad=38)
    tidy(ax)
    add_ribbon(ax, segs, c, total)


def panel_amplitudes(ax, c):
    cols = c["colors"]
    s90, dmin, q = c["f0b_dependent"]["s90_um"], c["f0b_dependent"]["d_min_um"], c["ustep_um"]
    ys = {"V": 3, "P": 2, "M": 1, "Q": 0}
    # V: 100 um ramp; the slow, fine-resolution part spans 0 ... 2 s90.
    ax.hlines(ys["V"], .039, c["v_ramp_um"], color=cols["V"], lw=3.5, alpha=.25)
    ax.hlines(ys["V"], .039, c["v_slow_factor_of_s90"] * s90, color=cols["V"], lw=3.5, alpha=.8)
    ax.text(c["v_ramp_um"] * 1.35, ys["V"], "100 µm ramp, slow to 2$s_{90}$", color=MUTED,
            fontsize=7, va="center")
    grid = strokes(c)
    ax.hlines(ys["P"], min(grid), max(grid), color=cols["P"], lw=3.5, alpha=.22)
    ax.scatter(grid, np.full(len(grid), ys["P"]), s=18, facecolor="white",
               edgecolor=cols["P"], lw=.9, zorder=3)
    ax.text(.045, ys["P"] + .2, f"{len(grid)} strokes · ×3 positions", color=MUTED,
            fontsize=7, va="bottom")
    so, si = c["m_outer_factor"] * s90, c["m_inner_factor"] * s90
    ax.scatter([si, so], [ys["M"]] * 2, s=22, marker="s", facecolor="white",
               edgecolor=cols["M"], lw=1, zorder=3)
    ax.plot([si, si, so, so], [ys["M"] + .14, ys["M"] + .24, ys["M"] + .24, ys["M"] + .14],
            color=cols["M"], lw=.8)
    ax.text(np.sqrt(si * so), ys["M"] + .3, "inner inside outer", color=MUTED, fontsize=7,
            ha="center", va="bottom")
    ax.scatter([s90], [ys["Q"]], s=24, marker="D", facecolor="white", edgecolor=cols["Q"],
               lw=1, zorder=3)
    ax.text(s90 / 1.4, ys["Q"], "half-cycle 1.3 / 4 / 12 s", color=MUTED, fontsize=7,
            ha="right", va="center", bbox=dict(fc="white", ec="none", pad=.3))

    ax.axvspan(s90 / 2, 2 * s90, color="#9aa7ad", alpha=.14, lw=0, zorder=0)
    guides = [
        (dmin, "$d_{\\min}$\n5 nm (F0b)", MUTED, 3.5),
        (q, "1 µstep\n39 nm", MUTED, 3.5),
        (c["beta_0_um"], rf"$\beta_0$ = {c['beta_0_um']:.2f} µm", "#9b4d48", 3.5),
        (c["electrical_cycle_um"], "electrical cycle\n40 µm", MUTED, 3.5),
    ]
    for value, label, color, y in guides:
        ax.axvline(value, color=color, lw=.7, ls=(0, (2, 2)), alpha=.82, zorder=0)
        ax.annotate(label, xy=(value, y), xytext=(4, 0), textcoords="offset points",
                    ha="left", va="bottom", color=color, fontsize=7, linespacing=.9)
    ax.text(s90, 2.5, "$s_{90}$ band (F0b)", color="#5d6b72", fontsize=7, ha="center",
            va="center")
    ax.set_xscale("log")
    ax.set(xlim=(.0025, 6500), ylim=(-.45, 4.1), yticks=list(ys.values()),
           yticklabels=["V  virgin", "P  loops", "M  nested", "Q  rate"],
           xlabel="commanded amplitude [µm]")
    ax.set_title("(b)", loc="left", pad=4)
    for tick, name in zip(ax.get_yticklabels(), ys):
        tick.set_color(cols[name])
        tick.set_fontweight("bold")
    tidy(ax)


def first_harmonic_amplitude(x_mm):
    """A1(x) from the frozen E2/E5 geometric model (spline-modulated first harmonic)."""
    sys.path.insert(0, str(GEOM))
    from e5_geometry import FrozenModel, bspline_basis
    m = FrozenModel()
    names, beta = m.data["state"]["names"], np.asarray(m.data["beta"])
    a = beta[[i for i, n in enumerate(names) if n.startswith("a1_")]]
    b = beta[[i for i, n in enumerate(names) if n.startswith("b1_")]]
    x = np.clip(x_mm, *m.domain)
    cb = bspline_basis(x, int(m.Kc), m.domain)
    return np.hypot(cb @ a, cb @ b), m.domain


def panel_map(ax, c):
    cols = c["colors"]
    travel, (lo, hi) = 120.0, c["constrained_span_mm"]
    ax.add_patch(plt.Rectangle((0, -.18), travel, .36, fc="white", ec="#7d858c", lw=.8))
    for a, b in ((0, lo), (hi, travel)):
        ax.add_patch(plt.Rectangle((a, -.18), b - a, .36, fc="#d9dde0", ec="none"))
    for name, (a, b) in c["features_mm"].items():
        if a == b:
            ax.plot([a, a], [-.18, .18], color=INK, lw=1)
        else:
            ax.add_patch(plt.Rectangle((a, -.18), b - a, .36, fc="none", ec=INK,
                                       hatch="////", lw=.6))
        ax.text((a + b) / 2, -.28, name, ha="center", va="top", fontsize=7, color=MUTED)
    for r in c["reference_mm"]:
        ax.plot(r, 0, "o", ms=6, mfc="white", mec=cols["R"], mew=1.2, zorder=4)
        ax.text(r, -.28, f"ref {r:g}", ha="center", va="top", fontsize=7, color=cols["R"])
    for name, x in c["positions_mm"].items():
        theta = (round(x * 1600) % 3200) / 3200 * 360
        ax.plot(x, 0, "o", ms=7, color=INK, zorder=5)
        ax.text(x, .27, f"$x_{name}$ = {x:g} mm\n" + rf"$\theta \approx {theta:.0f}^\circ$",
                ha="center", va="bottom", fontsize=7, color=INK, linespacing=1.0)
    # Faint first-harmonic amplitude above the bar (own scale, right axis).
    xs = np.linspace(0, travel, 600)
    amp, dom = first_harmonic_amplitude(xs)
    inside = (xs >= dom[0]) & (xs <= dom[1])
    strip = ax.inset_axes([0, 1.06, 1, .55], sharex=ax)
    grey = "#8c979e"
    strip.plot(xs[inside], amp[inside], color=grey, lw=.9)
    for x in c["positions_mm"].values():
        a1, _ = first_harmonic_amplitude(np.array([x]))
        strip.plot(x, a1[0], "o", ms=3.2, color=INK)
    strip.set_ylim(0, 1.2 * amp[inside].max())
    strip.set_yticks([0, round(amp[inside].max(), 1)])
    strip.tick_params(labelsize=7, colors=grey, length=2, labelbottom=False, bottom=False)
    strip.set_ylabel("$A_1$ [µm]", fontsize=7, color=grey, labelpad=2)
    strip.spines[["top", "right", "bottom"]].set_visible(False)
    strip.spines["left"].set_color("#c3c9cd")
    ax.set(xlim=(-2, 122), ylim=(-.75, .95), yticks=[], xlabel="stage position [mm]",
           xticks=np.arange(0, 121, 20))
    ax.spines[["top", "right", "left"]].set_visible(False)
    strip.set_title("(c)", loc="left", pad=2)


def caption(fig, c, summary):
    text = (
        "F1 command sequence. (a) Commanded position relative to the current test position "
        "on a symmetric logarithmic scale; the ribbon marks the three test positions and, "
        "within each, the virgin curves (V), symmetric loops (P), nested loops (M) and, at "
        "the middle position, the rate check (Q). Each virgin ramp reverses the preceding "
        "3 mm wipe. Reference pairs (purple dashed) are drawn relative to the test position "
        "of the adjacent group, so R1–R3 reach different heights; wipes are grey and "
        "inter-position transits grey dashed. (b) Amplitudes covered by each "
        "sub-block; the loop strokes are dense around the pre-sliding length s₉₀ measured "
        "in F0b and extend to 3 mm, bracketing the expected saturation lost motion β₀. "
        "(c) Test positions along the travel at different screw angles and first-harmonic "
        "wobble amplitudes A₁ (grey line, frozen geometric model), inside the "
        "well-constrained span and away from known local features. "
        f"Placeholder defaults s₉₀ = {c['f0b_dependent']['s90_um']:g} µm and d_min = "
        f"{c['f0b_dependent']['d_min_um'] * 1000:g} nm until F0b is processed. "
        f"{summary['segments']} segments, {summary['duration_min']:.1f} min."
    )
    fig.text(.115, .012, textwrap.fill(text, width=132), fontsize=7, color=MUTED,
             va="bottom", linespacing=1.35)


def main() -> None:
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    p, r2 = plan(c)
    summary = summarise(p, r2)
    total = p.now / 60
    style()
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.3, 1.2, .5], left=.115, right=.955,
                          top=.905, bottom=.19, hspace=.62)
    panel_trajectory(fig.add_subplot(gs[0]), p.segs, c, total)
    panel_amplitudes(fig.add_subplot(gs[1]), c)
    panel_map(fig.add_subplot(gs[2]), c)
    caption(fig, c, summary)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=320, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    print(f"{OUTPUT}\n{summary['segments']} segments, {total:.1f} min, R2 {'in' if r2 else 'out'}")


if __name__ == "__main__":
    main()

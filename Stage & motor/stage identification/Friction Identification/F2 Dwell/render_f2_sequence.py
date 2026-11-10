#!/usr/bin/env python3
"""Render the paper-width F2 dwell sequence overview (procedure §8).

    python3 render_f2_sequence.py                   # the variant in f2_config.json
    python3 render_f2_sequence.py --variant reduced # procedure §4.3

writes rendered/f2_sequence.png (or f2_sequence_reduced.png).

The figure is drawn from the timed plan built by f2_sequence.py, so it cannot
drift from the table the board replays: both come from f2_config.json and the
two builders cross-check each other. Time at rest, not amplitude, is the
variable here, so the dwell is the thing the eye has to be able to measure:
every dwell is a flat bar on the zero line whose length is its duration and
whose shade is its duration, and the post-dwell ramp separates the two modes
by direction alone.
"""
from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.cm import ScalarMappable

from f2_sequence import CONFIG, blocks_of, build, checks, ramp_um, summarise, trial_order

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered"

INK = "#20262b"
MUTED = "#59636b"
GRID = "#dfe3e6"
BLOCK_SHADE = {"B": "#e6eaec", "A": "#eff1f3", "L": "#e6eaec"}
DWELL_EDGE = "#6b7478"


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
        ax.grid(color=GRID, lw=.5, zorder=0)


def dwell_scale(c: dict):
    """Sequential green, light = shortest dwell, dark = longest, on a log norm.

    A mid stop is given explicitly: the two end colours in f2_config.json
    interpolated on their own put every dwell below ~10 s in the palest
    quarter, where 0.1 s and 2 s are not tellable apart on paper. The light
    end is a visible green rather than near-white, so that a filled 0.1 s
    marker in (b) cannot be mistaken for an open one.
    """
    cmap = LinearSegmentedColormap.from_list(
        "dwell", [c["colors"]["dwell_light"], "#5c9f73", c["colors"]["dwell_dark"]])
    return cmap, LogNorm(vmin=.08, vmax=400.0)


def spans(segs, key):
    """Contiguous runs of equal key, as (key, t0_min, t1_min)."""
    out = []
    for s in segs:
        k = key(s)
        if out and out[-1][0] == k:
            out[-1][2] = s.t1
        else:
            out.append([k, s.t0, s.t1])
    return [(k, a / 60, b / 60) for k, a, b in out]


# -- (a) ribbon + trajectory ----------------------------------------------
def add_ribbon(ax, segs, c, total) -> None:
    """Two tiers: blocks on top, repetitions inside each block below."""
    cols = c["colors"]
    pos = c["positions_mm"]
    top = ax.inset_axes([0, 1.105, 1, .062])
    low = ax.inset_axes([0, 1.030, 1, .062])
    longest = max(s.td_s for s in segs)
    for b, a, z in spans(segs, lambda s: s.block):
        ref = b.startswith("R")
        top.barh(0, z - a, left=a, height=.9, color=cols["R"] if ref else "#49535b",
                 edgecolor="white", lw=.6)
        if b == "L":
            label = f"L · {longest:g} s"
        elif b in pos:
            label = f"{b} · {pos[b]:g} mm"
        else:
            label = b
        top.text((a + z) / 2, 0, label, ha="center", va="center", color="white",
                 weight="bold", fontsize=7)
    for (b, rep), a, z in spans(segs, lambda s: (s.block, s.rep)):
        if not rep:                     # the reference blocks have no repetitions
            low.barh(0, z - a, left=a, height=.9, color="#d6dadd", edgecolor="white", lw=.6)
            continue
        low.barh(0, z - a, left=a, height=.9, color="#8e9aa1" if rep % 2 else "#73808a",
                 edgecolor="white", lw=.6)
        if z - a > 1.2:
            low.text((a + z) / 2, 0, f"rep {rep}", ha="center", va="center", color="white",
                     weight="bold", fontsize=6.6)
    for r in (top, low):
        r.set(xlim=(0, total), ylim=(-.5, .5), xticks=[], yticks=[])
        r.spines[:].set_visible(False)
    low.text(-.006, 0, "reps", transform=low.transAxes, ha="right", va="center",
             fontsize=6.6, color=MUTED)


def dwell_key(ax, c) -> None:
    """Small colour key, in the band above the trial ramps, which is empty."""
    cmap, norm = dwell_scale(c)
    cax = ax.inset_axes([.44, .845, .17, .027])
    bar = ax.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                             orientation="horizontal")
    bar.set_ticks([.1, 1, 10, 100], labels=["0.1", "1", "10", "100"])
    bar.ax.minorticks_off()
    bar.ax.tick_params(labelsize=6.4, length=1.6, pad=1.4, colors=MUTED)
    bar.outline.set(edgecolor=DWELL_EDGE, linewidth=.5)
    cax.set_title("dwell time [s]", fontsize=6.8, color=MUTED, pad=2.2)


def panel_trajectory(ax, segs, c, total) -> None:
    cols = c["colors"]
    cmap, norm = dwell_scale(c)
    setup, transit, ref, ramps = [], [], [], []
    for s in segs:
        line = [(s.t0 / 60, s.u0 - s.ref_um), (s.t1 / 60, s.u1 - s.ref_um)]
        if s.kind == "transit":
            transit.append(line)
        elif s.kind == "ref":
            ref.append(line)
        elif s.kind == "measure":
            ramps.append(line)
        elif s.kind in ("setup", "settle", "hold"):
            setup.append(line)
    ax.add_collection(LineCollection(setup, colors=cols["setup"], linewidths=.7, zorder=2))
    ax.add_collection(LineCollection(transit, colors=cols["transit"], linewidths=.9,
                                     linestyles=(0, (3, 2)), zorder=2))
    ax.add_collection(LineCollection(ref, colors=cols["R"], linewidths=1.0,
                                     linestyles=(0, (4, 2)), zorder=4))

    # The dwells: flat bars on the zero line, length and shade both the duration.
    # A thin dark underlay gives the palest shades an outline on white paper.
    bars = [[(s.t0 / 60, 0.0), (s.t1 / 60, 0.0)] for s in segs if s.kind == "dwell"]
    tds = [s.td_s for s in segs if s.kind == "dwell"]
    ax.add_collection(LineCollection(bars, colors=DWELL_EDGE, linewidths=3.6, zorder=5))
    ax.add_collection(LineCollection(bars, colors=cmap(norm(tds)), linewidths=2.4,
                                     zorder=6))
    # The post-dwell ramps last, so the mode is readable over the dwell bars.
    ax.add_collection(LineCollection(ramps, colors=cols["ramp"], linewidths=1.15, zorder=7))

    for b, a, z in spans(segs, lambda s: s.block):
        if b in BLOCK_SHADE:
            ax.axvspan(a, z, color=BLOCK_SHADE[b], alpha=.6, lw=0, zorder=-2)
    ax.axhspan(-1, 1, color="#f4f5f6", zorder=-1)

    # Where the trajectory moves to another test position, say so once.
    pos = c["positions_mm"]
    seen = None
    for b, a, _z in spans(segs, lambda s: s.block):
        x = pos.get(b if b in pos else "B")
        if b.startswith("R") or x == seen:
            continue
        seen = x
        ax.text(a + .3, -3.0e4, rf"$\rightarrow$ {x:g} mm", fontsize=7, color=MUTED,
                ha="left", va="top", bbox=dict(fc="white", ec="none", pad=.4), zorder=8)

    ax.set_yscale("symlog", linthresh=1.0, linscale=.55)
    ticks = [-1e4, -1e3, -100, -10, -1, 1, 10, 100, 1e3, 1e4]
    labels = ["−10 mm", "−1 mm", "−100 µm", "−10 µm", "−1 µm",
              "+1 µm", "+10 µm", "+100 µm", "+1 mm", "+10 mm"]
    ax.set_yticks(ticks, labels)
    ax.set(xlim=(0, total), ylim=(-6.5e4, 6.5e4),
           ylabel=r"commanded position relative to $x_p$",
           xlabel="elapsed session time [min]")
    ax.set_title("(a)", loc="left", pad=38)
    tidy(ax)
    legend(ax, c)
    dwell_key(ax, c)
    add_ribbon(ax, segs, c, total)


def legend(ax, c) -> None:
    cols = c["colors"]
    handles = [
        plt.Line2D([], [], color=cols["setup"], lw=.9, label="wipe, settle, hold"),
        plt.Line2D([], [], color=cols["transit"], lw=.9, ls=(0, (3, 2)), label="transit"),
        plt.Line2D([], [], color=cols["R"], lw=1.0, ls=(0, (4, 2)), label="reference pair"),
        plt.Line2D([], [], color="#5c9f73", lw=2.4, label="dwell"),
        plt.Line2D([], [], color=cols["ramp"], lw=1.15, label="post-dwell ramp"),
    ]
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(.065, .995),
                    frameon=True, fontsize=6.8, handlelength=1.6, borderpad=.4,
                    labelspacing=.28, ncols=1)
    leg.get_frame().set(edgecolor="#d3d8db", facecolor="white", linewidth=.5)
    leg.set_zorder(9)


# -- (b) dwell coverage ----------------------------------------------------
def panel_coverage(ax, c, blocks, order) -> None:
    cmap, norm = dwell_scale(c)
    rows, labels = [], []
    for b in blocks:
        for mode in b["modes"]:
            long = len(b["dwells"]) == 1 and b["dwells"][0]["s"] >= 300
            name = "long dwell" if long else ("reversal" if mode == "R" else "continuation")
            rows.append((b["name"], mode))
            labels.append(f"{b['name']} · {name}")
    y = {key: len(rows) - 1 - i for i, key in enumerate(rows)}

    lo, hi = c["operational_dwell_range_s"]
    ax.axvspan(lo, hi, color="#9aa7ad", alpha=.13, lw=0, zorder=0)
    ax.text(.055, -1.16, f"shaded: expected operational dwell range, {lo:g}–{hi:g} s "
            "(to be agreed)", ha="left", va="center", fontsize=6.8, color="#5d6b72")

    # 1, 2 and 2.5 s sit within a fifth of a decade, so the labels descend a
    # staircase rather than sharing a line.
    guides = [(c["earlier_dwells_s"]["f0b_probe_lo"], "F0b probe", 2),
              (c["earlier_dwells_s"]["f1_virgin_settle"], "F1 settle", 1),
              (c["earlier_dwells_s"]["e5_iso_stop"], "E5 stop", 0),
              (c["earlier_dwells_s"]["f0b_probe_hi"], "F0b probe", 0)]
    top = len(rows) - .55
    for value, label, tier in guides:
        ax.axvline(value, color=MUTED, lw=.7, ls=(0, (1.6, 1.8)), alpha=.85, zorder=1)
        ax.annotate(label, xy=(value, top + .22 + .38 * tier), xytext=(3.2, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=6.8,
                    color=MUTED)
    ax.text(.055, top + 1.52, "dotted: dwells already in F0b, F1 and E5",
            ha="left", va="center", fontsize=6.8, color=MUTED, style="italic")

    counts: dict[tuple, int] = {}
    for name, _rep, td, mode in order:
        n = counts.get((name, mode, td), 0)
        counts[(name, mode, td)] = n + 1
        shade = cmap(norm(td))
        ax.scatter(td, y[(name, mode)] + .13 * n, s=26, zorder=4,
                   facecolor=shade if mode == "R" else "white",
                   edgecolor=DWELL_EDGE, linewidth=.8)
    for (name, mode), row in y.items():
        tds = sorted({td for n, _r, td, m in order if (n, m) == (name, mode)})
        ax.plot([min(tds), max(tds)], [row, row], color="#c9d0d4", lw=.8, zorder=2)

    ax.scatter([], [], s=26, facecolor="#5c9f73", edgecolor=DWELL_EDGE, lw=.8,
               label="reversal (R)")
    ax.scatter([], [], s=26, facecolor="white", edgecolor=DWELL_EDGE, lw=.8,
               label="continuation (C)")
    leg = ax.legend(loc="lower right", bbox_to_anchor=(1.0, -.02), fontsize=6.8,
                    frameon=True, borderpad=.4, labelspacing=.3, handletextpad=.3,
                    columnspacing=1.0, ncols=2)
    leg.get_frame().set(edgecolor="#d3d8db", facecolor="white", linewidth=.5)

    ax.set_xscale("log")
    ax.set(xlim=(.05, 500), ylim=(-1.40, top + 1.80), yticks=list(y.values()),
           yticklabels=labels, xlabel="dwell time $T_d$ [s]")
    ax.set_title("(b)", loc="left", pad=4)
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")
        tick.set_color(INK)
    tidy(ax)
    ax.grid(axis="y", visible=False)


# -- (c) trial anatomy ----------------------------------------------------
def break_mark(ax, x, yc, h) -> None:
    """A // gap in a bar: white sliver, then two slanted strokes."""
    w = .085
    ax.add_patch(plt.Rectangle((x - w, yc - h / 2), 2 * w, h, fc="white", ec="none",
                               zorder=6))
    for dx in (-w * .55, w * .55):
        ax.plot([x + dx - .05, x + dx + .05], [yc - h / 2, yc + h / 2],
                color=DWELL_EDGE, lw=.8, zorder=7, solid_capstyle="butt")


def bar(ax, x0, x1, yc, h, color, edge=DWELL_EDGE, lw=.6, z=5):
    ax.add_patch(plt.Rectangle((x0, yc - h / 2), x1 - x0, h, fc=color, ec=edge, lw=lw,
                               zorder=z))


def panel_anatomy(ax, c, short_td, long_td) -> None:
    cmap, norm = dwell_scale(c)
    cols = c["colors"]
    span = ramp_um(c)
    ramp_s = span / c["ramp_rate_um_s"]
    wipe_s = c["wipe_mm"] * 1000.0 / (c["transit_mm_s"] * 1000.0)

    X_WIPE = (-4.55, -3.35)
    X_SETTLE = (-3.35, -3.00)
    X_LONG = (-3.00, 0.0)
    X_SHORT = (-1.00, 0.0)
    X_RAMP = (0.0, 1.10)
    X_HOLD = (1.10, 1.50)
    H = .44          # tall enough to carry the long dwell's label inside the bar

    ax.plot([0, 0], [-1.42, 1.55], color=INK, lw=.7, ls=(0, (1.6, 1.8)), zorder=3)
    ax.text(.06, 1.62, "end of the dwell ($t = 0$)", fontsize=7, color=INK, ha="left",
            va="center")

    for mode, yc in (("R", .85), ("C", -.85)):
        direction = -1 if mode == "R" else 1
        bar(ax, *X_WIPE, yc, H, cols["setup"])
        bar(ax, *X_SETTLE, yc, H, "#d6dadd")
        bar(ax, *X_LONG, yc, H, cmap(norm(long_td)))
        break_mark(ax, -2.55, yc, H)
        bar(ax, *X_SHORT, yc, H * .46, cmap(norm(short_td)), lw=.5, z=8)
        end = yc + direction * .46
        ax.plot(X_RAMP, [yc, end], color=cols["ramp"], lw=1.6, solid_capstyle="round",
                zorder=6)
        ax.plot(X_HOLD, [end] * 2, color=cols["setup"], lw=1.3, zorder=6)
        ax.text(X_HOLD[1] + .10, end + .12,
                "down: reversal (R)" if mode == "R" else "up: continuation (C)",
                fontsize=7.2, color=cols["ramp"], weight="bold", ha="left", va="center")

    # Part names above the top (R) timeline, durations below the bottom (C) one.
    for x, text in ((np.mean(X_WIPE), f"{c['wipe_mm']:g} mm wipe from below"),
                    (-1.95, "dwell")):
        ax.text(x, .85 + H / 2 + .09, text, fontsize=7, color=MUTED, ha="center",
                va="bottom")
    base = -.85 - H / 2
    for x, text in ((np.mean(X_WIPE), f"{wipe_s:.1f} s"),
                    (np.mean(X_SETTLE), f"settle {c['pre_dwell_settle_s']:g} s")):
        ax.text(x, base - .14, text, fontsize=7, color=MUTED, ha="center", va="top")
    # The two drawn lengths: the long one named on the bar itself (there is no
    # room for it below, under the settle label), the short one on a leader.
    ax.text(-1.95, .85, f"$T_d$ = {long_td:g} s", fontsize=7, color="white",
            weight="bold", ha="center", va="center", zorder=9)
    ax.plot([X_SHORT[0]] * 2, [base, base - .41], color=MUTED, lw=.6, zorder=4)
    ax.text(X_SHORT[0], base - .44, f"$T_d$ = {short_td:g} s", fontsize=7, color=MUTED,
            ha="center", va="top")
    ax.text(.78, base - .44,
            f"{span:.1f} µm at {c['ramp_rate_um_s']:g} µm/s, {ramp_s:.0f} s",
            fontsize=7, color=MUTED, ha="center", va="top")
    ax.text(1.70, base - .44, f"hold {c['post_ramp_hold_s']:g} s",
            fontsize=7, color=MUTED, ha="left", va="top")

    # What each part is for, bracketed under the whole panel.
    for (x0, x1), text in ((X_LONG, "creep measured here"),
                           (X_RAMP, "response measured here")):
        y = -2.15
        ax.plot([x0, x0, x1, x1], [y + .10, y, y, y + .10], color=INK, lw=.7)
        ax.text((x0 + x1) / 2, y - .08, text, fontsize=7.2, color=INK, ha="center",
                va="top", style="italic")

    ax.set(xlim=(-4.85, 3.45), ylim=(-2.60, 1.80), xticks=[], yticks=[.85, -.85],
           yticklabels=["R", "C"])
    ax.set_xlabel("time, aligned at the end of the dwell — not to scale; "
                  "// marks the elided part of the long dwell", labelpad=4)
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")
        tick.set_fontsize(9)
        tick.set_color(INK)
    ax.set_title("(c)", loc="left", pad=4)
    ax.spines[:].set_visible(False)
    ax.tick_params(axis="y", length=0)


# -- caption ---------------------------------------------------------------
def caption(fig, c, summary, short_td, long_td) -> None:
    f0b = c["f0b_dependent"]
    text = (
        "F2 command sequence. (a) Commanded position relative to the current test "
        "position on a symmetric logarithmic scale; each trial is a 3 mm wipe from "
        "below (grey), a dwell on the zero line (green, shaded and sized by its "
        f"duration) and a {summary['ramp_um']:.1f} µm post-dwell ramp, down (reversal, R) "
        "or up (continuation, C), so the mode is readable from direction alone. The "
        "ribbon marks the blocks and, inside each, the repetitions; the order of dwell "
        "times and modes is shuffled per repetition, which is why short and long flat "
        "segments alternate irregularly. Reference pairs at 15 and 105 mm (purple "
        "dashed) are drawn against the test position of the adjacent block, so R1 and R3 "
        "reach ±45 mm. (b) Dwell times covered at each position, one marker per planned "
        "trial stacked by repetition, filled for reversal and open for continuation; "
        "dotted lines mark the dwells already present in F0b, F1 and E5, and the shaded "
        "band the dwell range expected in service. (c) Anatomy of one trial per mode, "
        f"aligned at the end of the dwell, with the dwell drawn at {short_td:g} s and "
        f"{long_td:g} s: creep is measured during the dwell, the friction response during "
        "the ramp. The ramp length is 2s₉₀ and is provisional — placeholder defaults "
        f"s₉₀ = {f0b['s90_um']:g} µm and d_min = {f0b['d_min_um'] * 1000:g} nm until F0b "
        f"is processed. {summary['trials']} trials, {summary['segments']} segments, "
        f"{summary['duration_min']:.1f} min of a {c['gate_min']:.0f} min gate, "
        f"{summary['dwelling_min']:.1f} min of it at rest."
    )
    fig.text(.105, .012, textwrap.fill(text, width=132), fontsize=7, color=MUTED,
             va="bottom", linespacing=1.35)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variant", choices=("full", "reduced"), default=None)
    args = ap.parse_args()

    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    variant = args.variant or c["variant"]
    blocks, references = blocks_of(c, variant)
    order = trial_order(c, blocks)
    p = build(c, blocks, references)
    summary = summarise(p, variant, blocks, order, checks(p, blocks, order))
    total = p.now / 60.0

    tds = sorted({td for _b, _r, td, _m in order})
    short_td = next((t for t in tds if t >= 2.0), tds[0])
    long_td = max(t for t in tds if t <= 120.0) if any(t <= 120 for t in tds) else max(tds)

    style()
    fig = plt.figure(figsize=(7.2, 7.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.15, 1.25, 1.05], left=.105, right=.965,
                          top=.905, bottom=.20, hspace=.72)
    panel_trajectory(fig.add_subplot(gs[0]), p.segs, c, total)
    panel_coverage(fig.add_subplot(gs[1]), c, blocks, order)
    panel_anatomy(fig.add_subplot(gs[2]), c, short_td, long_td)
    caption(fig, c, summary, short_td, long_td)

    OUT.mkdir(parents=True, exist_ok=True)
    name = "f2_sequence.png" if variant == "full" else f"f2_sequence_{variant}.png"
    fig.savefig(OUT / name, dpi=320, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    bad = [k for k, v in summary["checks"].items() if not v]
    print(f"{OUT / name}\n{summary['trials']} trials, {summary['segments']} segments, "
          f"{total:.1f} min, dwelling {summary['dwelling_min']:.1f} min, "
          f"variant {variant}")
    print("checks: " + ("all pass" if not bad else "FAILED " + ", ".join(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())

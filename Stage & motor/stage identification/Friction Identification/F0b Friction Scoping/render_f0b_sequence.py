#!/usr/bin/env python3
"""Render the paper-width F0b friction-scoping sequence overview."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f0b_config.json"
OUTPUT = HERE / "rendered" / "f0b_sequence.png"

INK = "#20262b"
MUTED = "#59636b"
GRID = "#dfe3e6"
SETUP = "#b8bec2"
REFERENCE = "#73578f"


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


def block_bounds(config: dict) -> dict[str, tuple[float, float]]:
    out, start = {}, 0.0
    for block in config["blocks"]:
        end = start + block["duration_min"]
        out[block["name"]] = (start, end)
        start = end
    return out


def build_segments(config: dict) -> list[dict]:
    """Build a schematic trace retaining analysed/setup identity.

    Timing within a block is proportional and illustrative; block boundaries,
    ramp directions, wipe sides and the V-block speeds follow the procedure.
    """
    bounds = block_bounds(config)
    seg: list[dict] = []

    def add(block, t0, t1, x0, x1, measured=True, dashed=False):
        seg.append(dict(block=block, t0=t0, t1=t1, x0=x0, x1=x1,
                        measured=measured, dashed=dashed))

    # Reference pairs are dashed so they remain distinct from D in greyscale.
    for name in ("R1", "R2", "R3"):
        a, b = bounds[name]
        u = np.array([0.0, 0.24, 0.50, 0.74, 1.0])
        t = a + u * (b - a)
        x = np.array([0.0, -45000.0, 0.0, 45000.0, 0.0])
        for i in range(4):
            add(name, t[i], t[i + 1], x[i], x[i + 1], dashed=True)

    # N: three measured holds spanning the travel. Each is approached from
    # below; the 9 min of holds plus 2.7 min of moves/approaches fill the block.
    a, b = bounds["N"]
    now, previous = a, 0.0
    for target, hold_min in ((-45000.0, 2.0), (0.0, 5.0), (45000.0, 2.0)):
        add("N", now, now + .82, previous, target - 3000, False)
        add("N", now + .82, now + .90, target - 3000, target, False)
        add("N", now + .90, now + .90 + hold_min, target, target)
        now += .90 + hold_min
        previous = target
    assert np.isclose(now, b)

    # L: one setup wipe, then increasing and decreasing stroke lists.
    a, b = bounds["L"]
    add("L", a, a + .08, 20000, -3000, False)
    add("L", a + .08, a + .16, -3000, 0, False)
    strokes = np.asarray(config["ladder_strokes_um"])
    order = np.r_[strokes, strokes[::-1]]
    cell = (b - a - .16) / len(order)
    for i, stroke in enumerate(order):
        t0 = a + .16 + i * cell
        points = (0.0, stroke, 0.0, stroke, 0.0)
        for j in range(4):
            add("L", t0 + j * cell / 4, t0 + (j + 1) * cell / 4,
                points[j], points[j + 1])

    # S1/S2 reverse the wipe direction; S3/S4 continue it. Repeat once.
    a, b = bounds["S"]
    designs = ((-3000, -20), (3000, 20), (-3000, 20), (3000, -20)) * 2
    cell = (b - a) / len(designs)
    previous = 0.0
    for i, (wipe, target) in enumerate(designs):
        t0 = a + i * cell
        add("S", t0, t0 + .035 * cell, previous, wipe, False)
        add("S", t0 + .035 * cell, t0 + .085 * cell, wipe, 0, False)
        add("S", t0 + .085 * cell, t0 + cell, 0, target)
        previous = float(target)

    # D: setup wipe from below, analysed dwell, then 5 um reversal ramp (ABBA).
    a, b = bounds["D"]
    dwell_s = (1, 60, 60, 1)
    cell = (b - a) / 4
    previous = 0.0
    for i, dwell in enumerate(dwell_s):
        t0 = a + i * cell
        prep_end = t0 + .17 * cell
        add("D", t0, t0 + .07 * cell, previous, -3000, False)
        add("D", t0 + .07 * cell, prep_end, -3000, 0, False)
        dwell_fraction = .06 if dwell == 1 else .62
        dwell_end = prep_end + dwell_fraction * cell
        add("D", prep_end, dwell_end, 0, 0)
        add("D", dwell_end, t0 + cell, 0, -5)
        previous = -5.0

    # V: two fast cycles, one 100 s up + 100 s down cycle, two fast cycles.
    a, b = bounds["V"]
    fast_cycle = 2.0 / 60.0
    slow_half = 100.0 / 60.0
    used = 4 * fast_cycle + 2 * slow_half
    now = a + max((b - a - used) / 2, 0)
    for _ in range(2):
        add("V", now, now + fast_cycle / 2, 0, 1000)
        add("V", now + fast_cycle / 2, now + fast_cycle, 1000, 0)
        now += fast_cycle
    add("V", now, now + slow_half, 0, 1000)
    add("V", now + slow_half, now + 2 * slow_half, 1000, 0)
    now += 2 * slow_half
    for _ in range(2):
        add("V", now, now + fast_cycle / 2, 0, 1000)
        add("V", now + fast_cycle / 2, now + fast_cycle, 1000, 0)
        now += fast_cycle
    return seg


def add_ribbon(ax, config: dict) -> None:
    ribbon = ax.inset_axes([0, 1.025, 1, .072], transform=ax.transAxes)
    start = 0.0
    for block in config["blocks"]:
        width = block["duration_min"]
        ribbon.barh(0, width, left=start, height=.9, color=block["color"],
                    edgecolor="white", linewidth=.6)
        ribbon.text(start + width / 2, 0, block["name"], ha="center", va="center",
                    color="white", weight="bold", fontsize=7)
        start += width
    ribbon.set(xlim=(0, config["duration_min"]), ylim=(-.5, .5), xticks=[], yticks=[])
    ribbon.spines[:].set_visible(False)


def panel_trajectory(ax, config: dict) -> None:
    colors = {b["name"]: b["color"] for b in config["blocks"]}
    setup_lines, coloured, reference = [], {}, []
    for s in build_segments(config):
        line = [(s["t0"], s["x0"]), (s["t1"], s["x1"])]
        if not s["measured"]:
            setup_lines.append(line)
        elif s["dashed"]:
            reference.append(line)
        else:
            coloured.setdefault(s["block"], []).append(line)
    ax.add_collection(LineCollection(setup_lines, colors=SETUP, linewidths=.65, zorder=1))
    ax.add_collection(LineCollection(reference, colors=REFERENCE, linewidths=1.0,
                                     linestyles=(0, (4, 2)), zorder=3))
    for block, lines in coloured.items():
        # The ladder contains 104 half-cycles in a narrow time span. Keep the
        # individual oscillations visible rather than letting their strokes
        # merge into an apparently solid polygon.
        width = .42 if block == "L" else 1.35
        alpha = .58 if block == "L" else 1.0
        ax.add_collection(LineCollection(lines, colors=colors[block], linewidths=width,
                                         alpha=alpha, zorder=3))

    # Light ladder envelope and a thin boundary: the envelope is the design
    # message at full-paper scale, while the fine lines still show that it is
    # built from oscillatory loops.
    l0, l1 = block_bounds(config)["L"]
    strokes = np.asarray(config["ladder_strokes_um"])
    order = np.r_[strokes, strokes[::-1]]
    cell = (l1 - l0 - .16) / len(order)
    envelope_t = l0 + .16 + (np.arange(order.size) + .5) * cell
    ax.fill_between(envelope_t, 0, order, color=colors["L"], alpha=.10, zorder=2)
    ax.plot(envelope_t, order, color=colors["L"], lw=.8, zorder=3)

    ax.axhspan(-1, 1, color="#f4f5f6", zorder=-1)
    ax.axhline(45000, color=REFERENCE, lw=.6, ls=(0, (2, 2)), alpha=.55)
    ax.axhline(-45000, color=REFERENCE, lw=.6, ls=(0, (2, 2)), alpha=.55)
    ax.set_yscale("symlog", linthresh=1.0, linscale=.55)
    ticks = [-10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000]
    labels = ["−10 mm", "−1 mm", "−100 µm", "−10 µm", "−1 µm",
              "+1 µm", "+10 µm", "+100 µm", "+1 mm", "+10 mm"]
    ax.set_yticks(ticks, labels)
    ax.set(xlim=(0, config["duration_min"]), ylim=(-65000, 65000),
           ylabel=r"commanded position relative to $x_0$", xlabel="elapsed session time [min]")
    ax.set_title("(a) Commanded trajectory", loc="left", pad=25)
    tidy(ax)
    add_ribbon(ax, config)

    n0, _ = block_bounds(config)["N"]
    for x, y, label, va in ((n0 + 1.90, -45000, "N-end1", "bottom"),
                            (n0 + 6.30, 0, "N-mid", "bottom"),
                            (n0 + 10.70, 45000, "N-end2", "top")):
        ax.text(x, y, label, color=colors["N"], fontsize=7, weight="bold",
                ha="center", va=va, bbox=dict(fc="white", ec="none", pad=.3))
    right = config["duration_min"] - .35
    ax.text(right, 45000, "reference pair · 105 mm", ha="right", va="bottom",
            fontsize=7, color=REFERENCE, bbox=dict(fc="white", ec="none", pad=.5))
    ax.text(right, -45000, "reference pair · 15 mm", ha="right", va="top",
            fontsize=7, color=REFERENCE, bbox=dict(fc="white", ec="none", pad=.5))


def panel_amplitudes(ax, config: dict) -> None:
    colors = {b["name"]: b["color"] for b in config["blocks"]}
    ys = {"L": 3, "S": 2, "D": 1, "V": 0}
    strokes = np.asarray(config["ladder_strokes_um"])
    ax.hlines(ys["L"], strokes.min(), strokes.max(), color=colors["L"],
              lw=3.5, alpha=.22, zorder=1)
    ax.scatter(strokes, np.full(strokes.size, ys["L"]), s=18, marker="o",
               facecolor="white", edgecolor=colors["L"], lw=.9, zorder=3)
    ax.hlines(ys["S"], .039, 20, color=colors["S"], lw=3.5, alpha=.28)
    ax.vlines([.039, 20], ys["S"] - .12, ys["S"] + .12, color=colors["S"], lw=1.1)
    ax.plot(5, ys["D"], marker="|", ms=11, mew=1.6, color=colors["D"])
    ax.plot(1000, ys["V"], marker="|", ms=11, mew=1.6, color=colors["V"])

    guides = [
        (.002, "interferometer floor\n" + r"$\approx 2$ nm", MUTED),
        (.039, "1 µstep\n39 nm", MUTED),
        (config["beta_0_um"], rf"$\beta_0$ = {config['beta_0_um']:.3f} µm", "#9b4d48"),
        (40, "electrical cycle\n40 µm", MUTED),
    ]
    for value, label, color in guides:
        ax.axvline(value, color=color, lw=.7, ls=(0, (2, 2)), alpha=.82, zorder=0)
        ax.annotate(label, xy=(value, 3.48), xytext=(4, 0),
                    textcoords="offset points", ha="left", va="bottom",
                    color=color, fontsize=7, linespacing=.9)

    ax.text(.055, 3.14, "13 commanded strokes", color=MUTED, fontsize=7, va="bottom")
    ax.text(.055, 2.14, "20 µm ramp · 0.25 µm/s", color=MUTED, fontsize=7, va="bottom")
    ax.text(6.0, 1.0, "5 µm after 1 s / 60 s", color=MUTED, fontsize=7, va="center")
    ax.text(850, 0.0, "1 mm at 0.01 / 1 mm/s", color=MUTED, fontsize=7,
            ha="right", va="center")

    ax.set_xscale("log")
    ax.set(xlim=(.0011, 6500), ylim=(-.42, 3.96),
           yticks=list(ys.values()),
           yticklabels=["L  ladder", "S  ramps", "D  dwell", "V  speed"],
           xlabel="commanded amplitude [µm]")
    ax.set_title("(b) Amplitude coverage", loc="left", pad=6)
    for tick, name in zip(ax.get_yticklabels(), ys):
        tick.set_color(colors[name])
        tick.set_fontweight("bold")
    tidy(ax)


def caption(fig, config: dict) -> None:
    rest_min = next(item["duration_min"] for item in config["time_budget"]
                    if item["name"] == "holds, dwells and settles")
    rest_pct = round(100 * rest_min / config["duration_min"])
    text = (
        "F0b command sequence. (a) Commanded position relative to x₀ = 60 mm on a "
        "symmetric logarithmic scale; the ribbon marks reference pairs R1–R3, noise "
        "floor N at 15, 60 and 105 mm, reversal ladder L, slow ramps S, dwell probe D "
        "and speed probe V. "
        "Measurement segments are coloured; positioning moves and memory wipes are grey; "
        "V follows the fast–slow–fast order. "
        "(b) Amplitudes covered by each measurement block, from one 39 nm microstep to "
        "a 3 mm stroke; the ladder brackets the expected saturation lost motion "
        f"(β₀ = {config['beta_0_um']:.3f} µm). {config['segments']} segments, "
        f"{config['duration_min']:.1f} min, {config['gate_min'] - config['duration_min']:.1f} min "
        f"inside the gate, of which {rest_pct}% is at rest by design."
    )
    fig.text(.115, .025, textwrap.fill(text, width=130), fontsize=7, color=MUTED,
             va="bottom", linespacing=1.35)


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    assert np.isclose(sum(b["duration_min"] for b in config["blocks"]),
                      config["duration_min"])
    assert np.isclose(sum(b["duration_min"] for b in config["time_budget"]),
                      config["duration_min"])
    style()
    fig = plt.figure(figsize=(7.0, 6.25))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.25, 1.25], left=.115, right=.975,
                          top=.93, bottom=.20, hspace=.36)
    panel_trajectory(fig.add_subplot(gs[0]), config)
    panel_amplitudes(fig.add_subplot(gs[1]), config)
    caption(fig, config)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=320, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()

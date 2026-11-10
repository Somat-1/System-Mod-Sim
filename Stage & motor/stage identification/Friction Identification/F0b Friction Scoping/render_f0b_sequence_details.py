#!/usr/bin/env python3
"""Render detailed, predominantly visual views of the F0b L/S/D/V blocks."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f0b_config.json"
OUT = HERE / "rendered" / "sequence_details"

INK = "#20262b"
MUTED = "#59636b"
GRID = "#dfe3e6"
SETUP = "#b8bec2"


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
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def finish(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(which="major", axis="both", color=GRID, lw=.55, zorder=0)
    ax.grid(which="minor", axis="y", color=GRID, lw=.35, alpha=.45, zorder=0)
    ax.tick_params(which="both", direction="out")


def colours(config: dict) -> dict[str, str]:
    return {b["name"]: b["color"] for b in config["blocks"]}


def panel_l(ax, config: dict, panel: str = "") -> None:
    color = colours(config)["L"]
    strokes = np.asarray(config["ladder_strokes_um"])
    order = np.r_[strokes, strokes[::-1]]

    for i, stroke in enumerate(order):
        x = i + np.linspace(0, 1, 5)
        y = np.array([0, stroke, 0, stroke, 0])
        ax.plot(x, y, color=color, lw=.78, alpha=.88, zorder=3)

    ax.axvline(13, color="#aeb4b9", lw=.7, ls=(0, (3, 2)))
    box = dict(fc="white", ec="none", pad=.8, alpha=.9)
    ax.text(.24, .94, r"pass 1 · increasing $\rightarrow$", transform=ax.transAxes,
            ha="center", va="top", color=color, fontsize=7, bbox=box)
    ax.text(.76, .94, r"pass 2 · decreasing $\rightarrow$", transform=ax.transAxes,
            ha="center", va="top", color=color, fontsize=7, bbox=box)

    ax.set_yscale("symlog", linthresh=.039, linscale=.7)
    ax.set_yticks([0, .039, 1, 10, 100, 1000, 3000],
                  ["0", "39 nm", "1 µm", "10 µm", "100 µm", "1 mm", "3 mm"])
    ticks = [.5, 3.5, 6.5, 9.5, 12.5, 13.5, 16.5, 19.5, 22.5, 25.5]
    ax.set_xticks(ticks, ["1", "4", "7", "10", "13", "13", "10", "7", "4", "1"])
    ax.set(xlim=(0, 26), ylim=(-.02, 5000), xlabel="stroke number within pass",
           ylabel=r"commanded displacement from $x_0$")
    ax.set_title(f"{panel}L · reversal ladder", loc="left", pad=5)
    finish(ax)


def panel_s(ax, config: dict, panel: str = "") -> None:
    color = colours(config)["S"]
    designs = ((-3000, -20, "S1", "rev."), (3000, 20, "S2", "rev."),
               (-3000, 20, "S3", "cont."), (3000, -20, "S4", "cont.")) * 2
    ax.axvspan(0, 4, color=color, alpha=.025, zorder=-1)
    ax.axvspan(4, 8, color=color, alpha=.055, zorder=-1)
    previous = 0.0
    for i, (wipe, target, name, mode) in enumerate(designs):
        ax.plot([i, i + .12, i + .24], [previous, wipe, 0], color=SETUP, lw=.9, zorder=2)
        ax.annotate("", xy=(i + .96, target), xytext=(i + .24, 0),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.35), zorder=3)
        ax.text(i + .5, .965, f"{name}\n{mode}", transform=ax.get_xaxis_transform(),
                color=color, ha="center", va="top", fontsize=7, linespacing=.9,
                weight="bold" if name else "normal")
        previous = float(target)
    for x in range(1, 8):
        ax.axvline(x, color=GRID, lw=.5)
    ax.axvline(4, color="#9aa2a8", lw=.8, ls=(0, (3, 2)))
    ax.text(.24, .055, "repetition 1", transform=ax.transAxes, ha="center",
            va="bottom", color=MUTED, fontsize=7)
    ax.text(.76, .055, "repetition 2", transform=ax.transAxes, ha="center",
            va="bottom", color=MUTED, fontsize=7)

    ax.set_yscale("symlog", linthresh=1, linscale=.55)
    ax.set_yticks([-3000, -20, 0, 20, 3000], ["−3 mm", "−20 µm", "0", "+20 µm", "+3 mm"])
    ax.set_xticks(np.arange(8) + .5, ["S1", "S2", "S3", "S4"] * 2)
    ax.set(xlim=(0, 8), ylim=(-7000, 7000), xlabel="ramp order",
           ylabel=r"commanded displacement from $x_0$")
    ax.set_title(f"{panel}S · slow ramps", loc="left", pad=5)
    finish(ax)


def panel_d(ax, config: dict, panel: str = "") -> None:
    color = colours(config)["D"]
    dwells = (1, 60, 60, 1)
    now, previous = 0.0, 0.0
    centres, boundaries, starts = [], [], []
    for trial, dwell in enumerate(dwells):
        start = now
        starts.append(start)
        # Eight seconds of schematic positioning, then the literal dwell and
        # 20 s measured 5 um reversal ramp from the procedure.
        ax.plot([now, now + 4, now + 8], [previous, -3000, 0], color=SETUP, lw=.9)
        now += 8
        ax.plot([now, now + dwell], [0, 0], color=color, lw=2.2, solid_capstyle="butt")
        ax.plot([now + dwell, now + dwell + 20], [0, -5], color=color, lw=1.35)
        now += dwell + 20
        centres.append((start + now) / 2)
        boundaries.append(now)
        previous = -5.0
        if trial % 2:
            ax.axvspan(start, now, color=color, alpha=.035, zorder=-1)
    for x in boundaries[:-1]:
        ax.axvline(x, color=GRID, lw=.55)
    for i, (x, dwell) in enumerate(zip(centres, dwells), start=1):
        ax.text(x, .96, f"D{i}\n{dwell} s", transform=ax.get_xaxis_transform(),
                ha="center", va="top", color=color, fontsize=7, weight="bold",
                linespacing=.9)
    ax.text(.99, .05, r"wipe  $\rightarrow$  hold  $\rightarrow$  5 µm reversal",
            transform=ax.transAxes,
            ha="right", va="bottom", color=MUTED, fontsize=7,
            bbox=dict(fc="white", ec="none", pad=.6, alpha=.88))

    ax.set_yscale("symlog", linthresh=1, linscale=.55)
    ax.set_yticks([-3000, -5, 0, 5], ["−3 mm", "−5 µm", "0", "+5 µm"])
    ax.set_xticks(np.arange(0, now + 1, 30))
    ax.set(xlim=(0, now), ylim=(-6000, 10), xlabel="elapsed block time [s]",
           ylabel=r"commanded displacement from $x_0$")
    ax.set_title(f"{panel}D · dwell probe", loc="left", pad=5)
    finish(ax)


def panel_v(ax, config: dict, panel: str = "") -> None:
    color = colours(config)["V"]
    points = [(0.0, 0.0)]
    now = 0.0
    for _ in range(2):
        points.extend([(now + 1, 1000), (now + 2, 0)])
        now += 2
    slow_start = now
    points.extend([(now + 100, 1000), (now + 200, 0)])
    now += 200
    slow_end = now
    for _ in range(2):
        points.extend([(now + 1, 1000), (now + 2, 0)])
        now += 2
    p = np.asarray(points)
    ax.plot(p[:, 0], p[:, 1], color=color, lw=1.35, marker="o", ms=2.5, zorder=3)
    ax.axvspan(0, slow_start, color=color, alpha=.08)
    ax.axvspan(slow_start, slow_end, color=color, alpha=.035)
    ax.axvspan(slow_end, now, color=color, alpha=.08)
    ax.text(2, 1060, "fast ×2", ha="center", va="bottom", color=color, fontsize=7,
            weight="bold")
    ax.text((slow_start + slow_end) / 2, 1060, "slow ×1", ha="center", va="bottom",
            color=color, fontsize=7, weight="bold")
    ax.text(now - 2, 1060, "fast ×2", ha="center", va="bottom", color=color,
            fontsize=7, weight="bold")
    ax.text((slow_start + slow_end) / 2, 510, "100 s", ha="center", va="bottom",
            color=MUTED, fontsize=7)

    ax.set(xlim=(-2, now + 2), ylim=(-40, 1190), xlabel="elapsed block time [s]",
           ylabel=r"commanded displacement from $x_0$",
           yticks=[0, 500, 1000], yticklabels=["0", "0.5 mm", "1 mm"])
    ax.set_title(f"{panel}V · speed probe", loc="left", pad=5)
    finish(ax)


PANELS = (("L", panel_l, "f0b_L_ladder.png"),
          ("S", panel_s, "f0b_S_slow_ramps.png"),
          ("D", panel_d, "f0b_D_dwell_probe.png"),
          ("V", panel_v, "f0b_V_speed_probe.png"))


def save_figure(fig, png: Path) -> None:
    fig.savefig(png, dpi=320, bbox_inches="tight", pad_inches=.04)
    fig.savefig(png.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.04)


def render_l_individual(config: dict) -> None:
    """Publication version of L with the two passes separated for legibility."""
    color = colours(config)["L"]
    strokes = np.asarray(config["ladder_strokes_um"])
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.35), sharex=True, sharey=True)
    fig.subplots_adjust(left=.14, right=.98, top=.88, bottom=.16, hspace=.16)
    fig.suptitle("L · reversal ladder", x=.14, y=.965, ha="left", fontsize=10,
                 weight="bold", color=INK)
    fig.text(.98, .945, r"setup: 3 mm wipe from below $\rightarrow x_0$",
             ha="right", va="top", fontsize=7, color=MUTED)

    for ax, order, title in zip(axes, (strokes, strokes[::-1]),
                                ("(a) Increasing pass", "(b) Decreasing pass")):
        for i, stroke in enumerate(order):
            x = i + np.linspace(.06, .94, 5)
            y = np.array([0, stroke, 0, stroke, 0])
            ax.plot(x, y, color=color, lw=.9, zorder=3)
            # One hollow marker per commanded amplitude; the W-shaped line
            # itself carries the four half-cycles without marker clutter.
            ax.plot(x[1], stroke, marker="o", ms=3.0, mfc="white", mec=color,
                    mew=.8, zorder=4)
        ax.set_yscale("symlog", linthresh=.039, linscale=.7)
        ax.set_yticks([0, .039, 1, 10, 100, 1000, 3000],
                      ["0", "39 nm", "1 µm", "10 µm", "100 µm", "1 mm", "3 mm"])
        ax.set_ylim(-.02, 5000)
        ax.set_title(title, loc="left", pad=4, fontsize=8.5, color=INK)
        finish(ax)

    axes[-1].set_xticks(np.arange(13) + .5, [str(i) for i in range(1, 14)])
    axes[-1].set_xlim(0, 13)
    axes[-1].set_xlabel("stroke number within pass")
    fig.supylabel(r"commanded displacement from $x_0$", x=.025, fontsize=8.5,
                  color="#2d3338")
    fig.text(.98, .035,
             r"each stroke:  $x_0 \rightarrow x_0+s \rightarrow x_0 \rightarrow x_0+s \rightarrow x_0$",
             ha="right", va="bottom", fontsize=7, color=MUTED)
    save_figure(fig, OUT / "f0b_L_ladder.png")
    plt.close(fig)


def render_individual(config: dict) -> None:
    render_l_individual(config)
    for key, draw, filename in PANELS:
        if key == "L":
            continue
        fig, ax = plt.subplots(figsize=(7.0, 2.85))
        fig.subplots_adjust(left=.13, right=.98, top=.90, bottom=.21)
        draw(ax, config)
        save_figure(fig, OUT / filename)
        plt.close(fig)


def render_combined(config: dict) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(7.0, 9.6))
    fig.subplots_adjust(left=.13, right=.98, top=.975, bottom=.07, hspace=.68)
    for label, (_, draw, _) , ax in zip(("(a) ", "(b) ", "(c) ", "(d) "), PANELS, axes):
        draw(ax, config, panel=label)
    save_figure(fig, OUT / "f0b_sequence_details.png")
    plt.close(fig)


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    render_individual(config)
    render_combined(config)
    for path in sorted(OUT.iterdir()):
        if path.suffix in (".png", ".pdf"):
            print(path)


if __name__ == "__main__":
    main()

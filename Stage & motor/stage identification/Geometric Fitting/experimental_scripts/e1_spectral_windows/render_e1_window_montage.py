#!/usr/bin/env python3
"""The three E1 windows side by side, each cropped to its own 12 mm span.

Same scale in every panel - 12 mm of travel, 101 targets up and 101 back down
on a 0.12 mm grid - so the three windows can be compared directly. The
pre-approach and the overrun run off the top and bottom of each panel by
construction: they are outside the window, which is the point of them.

Writes rendered/e1_window_montage.png.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import e1_run
from e1_run import WINDOWS, OVERRUN_MM

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parents[3] / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, style   # noqa: E402

OUT = Path(os.environ.get("E1_OUT", _HERE / "rendered"))
COLOURS = {"W1": "#4c78a8", "W2": "#598a76", "W3": "#9d6770"}
USTEP_PER_MM = 1600
MINUS = "−"


def thesis_style():
    style()
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "figure.facecolor": "#ffffff",
        "axes.facecolor": "#ffffff",
        "savefig.facecolor": "#ffffff",
    })


def main():
    thesis_style()
    run = e1_run.load()
    hold = run.hold()
    st = run.stops

    span = max(run.events.when(w, "state_enter")[0] - run.events.when(w, "move_start", "lower_pre")[0]
               for w, _, _ in WINDOWS)
    ends = []
    for name, _, _ in WINDOWS:
        t0 = run.events.when(name, "move_start", "lower_pre")[0]
        m = st.window == name
        ends.append(st.record_end_s[m].max() - t0)
    width = max(ends) / 60.0 + 0.2

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.6), sharex=True)
    for ax, (name, a, b) in zip(axes, WINDOWS):
        colour = COLOURS[name]
        t0 = run.events.when(name, "move_start", "lower_pre")[0]
        t1 = t0 + width * 60.0
        i, j = run.index_at([t0, t1])
        tm = (run.t[i:j] - t0) / 60.0
        ax.plot(tm, run.pos[i:j], lw=1.1, color=colour, rasterized=True)

        m = st.window == name
        up, dn = m & (st.direction == 1), m & (st.direction == -1)
        for mask, mk in ((up, "o"), (dn, "s")):
            ax.plot((st.record_start_s[mask] - t0) / 60.0, hold[mask], mk,
                    ms=2.4, mfc="#ffffff", mec=INK, mew=0.6, zorder=4)

        t_turn = (run.events.when(name, "arrive", "upper_overrun")[0] - t0) / 60.0
        t_pre = (run.events.when(name, "arrive", "lower_pre")[0] - t0) / 60.0
        ax.axvline(t_turn, color=INK2, lw=0.8, ls=":", zorder=2)
        ax.plot([t_pre], [a], "v", ms=7, mfc="#ffffff", mec=INK2, mew=1.1,
                clip_on=False, zorder=6)
        ax.plot([t_turn], [b], "^", ms=7, mfc="#ffffff", mec=INK2, mew=1.1,
                clip_on=False, zorder=6)

        ax.set_ylim(a, b)
        ax.set_yticks(np.arange(a, b + 0.1, 2.0))
        ax.set_ylabel("position [mm]")
        ax.set_title(f"{name}: {a:.0f}–{b:.0f} mm — 101 targets up, "
                     f"101 down, 0.12 mm apart", fontsize=11, loc="left")
        ax.grid(axis="y", lw=0.5)

        rt = ax.twinx()
        rt.set_ylim(a * USTEP_PER_MM, b * USTEP_PER_MM)
        rt.set_yticks(np.arange(a, b + 0.1, 2.0) * USTEP_PER_MM)
        rt.set_ylabel("microsteps", fontsize=9)
        rt.tick_params(labelsize=8)
        rt.grid(False)

    axes[-1].set_xlabel("time since the window's pre-approach began [min]")
    axes[-1].set_xlim(0, width)
    handles = [
        plt.Line2D([], [], ls="none", marker="o", ms=5, mfc="#ffffff",
                   mec=INK, mew=0.8, label="hold, upward pass"),
        plt.Line2D([], [], ls="none", marker="s", ms=5, mfc="#ffffff",
                   mec=INK, mew=0.8, label="hold, downward pass"),
        plt.Line2D([], [], ls="none", marker="v", ms=7, mfc="#ffffff",
                   mec=INK2, mew=1.1,
                   label=f"pre-approach, {OVERRUN_MM:.0f} mm below the crop"),
        plt.Line2D([], [], ls="none", marker="^", ms=7, mfc="#ffffff",
                   mec=INK2, mew=1.1,
                   label=f"overrun and turnaround, {OVERRUN_MM:.0f} mm above"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 0.962))
    fig.suptitle("E1 windows, each cropped to its own 12 mm span",
                 fontsize=12.5, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "e1_window_montage.png", dpi=200)
    print("wrote", OUT / "e1_window_montage.png")


if __name__ == "__main__":
    main()

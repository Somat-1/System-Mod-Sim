#!/usr/bin/env python3
"""Render publication-style figures from the executable E3 plan.

E3 repeats E2's design at an independent set of 400 targets. The stratified
draw makes the target-to-target step deliberately irregular, so the overview
carries two zoom panels that resolve individual stops: the question the plan
has to answer in a figure is how far apart consecutive targets actually are.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "rendered"
OVERVIEW = OUT / "e3_global_map_sequence.png"

COL = {
    "INIT": "#69747c", "DATUM": "#73578f", "UP": "#1b6f8a", "TURN_REF": "#a16a2d",
    "DOWN": "#a85b45", "END_REF": "#73578f", "PARK": "#4d5961", "DONE": "#4d5961",
    "record": "#67527c", "settle": "#c59b50", "motion": "#697b86",
    "small": "#b4473f", "grid": "#dde1e4",
}
# The fine zoom is placed on a single-grid-point gap flanked by typical ones,
# so the figure shows the smallest and the median step at the same scale.
FINE_FIRST, FINE_LAST = 110, 118


def load():
    config = json.loads((ROOT / "e3_config.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "sequence/e3_plan_summary.json").read_text(encoding="utf-8"))
    with (ROOT / "sequence/e3_measurement_plan.csv").open(newline="", encoding="utf-8") as fh:
        plan = list(csv.DictReader(fh))
    return config, summary, plan


def move_time(distance_mm: float, c: dict) -> float:
    distance = abs(distance_mm)
    velocity, acceleration = c["velocity_mm_s"], c["acceleration_mm_s2"]
    ramp_distance = velocity * velocity / acceleration
    if distance <= ramp_distance:
        return 2.0 * math.sqrt(distance / acceleration)
    return 2.0 * velocity / acceleration + (distance - ramp_distance) / velocity


def schedule(plan: list[dict], c: dict):
    """Replay the plan on the same timing model the generator used."""
    now, position, state = 0.0, c["x_ref_mm"], "INIT"
    segments, blocks, stops = [], [], []
    block_start = 0.0
    for row in plan:
        if row["action"] == "state_enter":
            if now > block_start:
                blocks.append((block_start, now, state))
            state, block_start = row["state"], now
            continue
        if row["target_mm"]:
            target = float(row["target_mm"])
            end = now + move_time(target - position, c)
            segments.append((now, end, position, target, state, "move"))
            now, position = end, target
        settle = float(row["settle_s"] or 0)
        record = float(row["record_s"] or 0)
        hold_start = now
        if settle:
            segments.append((now, now + settle, position, position, state, "settle"))
            now += settle
        if record:
            segments.append((now, now + record, position, position, state, "record"))
            now += record
        if row["record"] == "True":
            stops.append((hold_start, now, position, state,
                          int(row["index"]) if row["index"] else 0,
                          row["direction"]))
    if now > block_start:
        blocks.append((block_start, now, state))
    return segments, blocks, stops, now


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.edgecolor": "#7d858c",
        "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c",
        "ytick.color": "#41474c",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def staircase(stops: list, lo_index: int, hi_index: int, direction: str):
    """The commanded position against session time over a run of stops."""
    chosen = [s for s in stops
              if s[5] == direction and lo_index <= s[4] <= hi_index]
    time, position = [], []
    for start, end, pos, _state, _index, _dir in chosen:
        time += [start, end]
        position += [pos, pos]
    return np.array(time), np.array(position), chosen


def render_overview(c, summary, segments, blocks, stops, total) -> None:
    total_min = total / 60.0
    fig = plt.figure(figsize=(7.45, 8.7))
    gs = fig.add_gridspec(4, 6, height_ratios=[0.62, 2.45, 1.75, 0.80],
                          left=0.095, right=0.975, top=0.895, bottom=0.068,
                          hspace=0.52, wspace=1.05)
    fig.text(0.105, 0.978, "Experiment E3 — independent global map", fontsize=16,
             weight="bold", color="#20262b", va="top")
    fig.text(0.105, 0.942,
             "400 new stratified targets on the 64-µstep grid · mirrored up/down passes · "
             "step size deliberately irregular",
             fontsize=8.5, color="#566069", va="top")

    ax0 = fig.add_subplot(gs[0, :])
    for start, end, state in blocks:
        width = (end - start) / 60.0
        ax0.barh(0, width, left=start / 60.0, height=0.40,
                 color=COL.get(state, "#7a7f84"), edgecolor="white", linewidth=0.55)
        if width > 0.9:
            label = state.replace("END_REF", "END").replace("TURN_REF", "TURN")
            ax0.text((start + end) / 120.0, 0, label, ha="center", va="center",
                     fontsize=6.1, color="white", weight="bold")
    ax0.set(xlim=(0, total_min), ylim=(-0.34, 0.34), yticks=[])
    ax0.set_title("(a) Controller-state allocation", loc="left", pad=6)
    ax0.tick_params(axis="x", bottom=False, labelbottom=False)
    ax0.spines[:].set_visible(False)

    ax = fig.add_subplot(gs[1, :], sharex=ax0)
    for state in ["DATUM", "UP", "TURN_REF", "DOWN", "END_REF", "PARK"]:
        lines = [[(a / 60.0, x0), (b / 60.0, x1)]
                 for a, b, x0, x1, st, _ in segments if st == state]
        if lines:
            ax.add_collection(LineCollection(lines, colors=COL[state], linewidths=1.0))
    ax.axhline(c["x_ref_mm"], color=COL["DATUM"], lw=0.7, ls=(0, (3, 2)))
    ax.text(0.5, c["x_ref_mm"] + 2.2, r"$x_{ref}=61.35$ mm", color=COL["DATUM"], fontsize=7)

    coarse_t, coarse_x, _ = staircase(stops, FINE_FIRST - 25, FINE_LAST + 25, "up")
    ax.axvspan(coarse_t[0] / 60.0, coarse_t[-1] / 60.0, color=COL["UP"], alpha=0.14, lw=0)
    ax.annotate("(c)", xy=(coarse_t[0] / 60.0, coarse_x.max() + 8), fontsize=7.4,
                color=COL["UP"], weight="bold", ha="center")

    ax.set(xlim=(0, total_min), ylim=(0, 124), ylabel="commanded position [mm]")
    ax.set_title("(b) Complete commanded trajectory", loc="left", pad=6)
    ax.grid(axis="x", color=COL["grid"], lw=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    box = (f"targets  {summary['targets']} (seed {summary['seed_used']})\n"
           f"analysis stops  {summary['analysis_stops']}\n"
           f"reference records  {summary['reference_stops']}\n"
           f"span  {summary['target_min_mm_actual']:.2f}–"
           f"{summary['target_max_mm_actual']:.2f} mm\n"
           f"step  {summary['spacing_min_mm']:.2f}–{summary['spacing_max_mm']:.2f} mm "
           f"(mean {summary['spacing_mean_mm']:.2f})\n"
           f"nominal duration  {total_min:.1f} min\n"
           f"margin to the 58 min gate  "
           f"{c['watchdog_projection_min'] - total_min:.1f} min")
    # Top left: the up pass climbs away from this corner, so it stays clear.
    ax.text(0.015, 0.97, box, transform=ax.transAxes, ha="left", va="top",
            fontsize=7.2, linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.48", fc="white", ec="#aeb4b9", lw=0.7))

    # (c) a coarse zoom: enough stops to show the staircase as a recording sees it.
    ax = fig.add_subplot(gs[2, 0:3])
    ax.plot(coarse_t - coarse_t[0], coarse_x, color=COL["UP"], lw=1.0,
            solid_joinstyle="miter")
    fine_t, fine_x, fine_stops = staircase(stops, FINE_FIRST, FINE_LAST, "up")
    ax.axvspan(fine_t[0] - coarse_t[0], fine_t[-1] - coarse_t[0],
               color=COL["small"], alpha=0.13, lw=0)
    ax.annotate("(d)", xy=((fine_t[0] + fine_t[-1]) / 2 - coarse_t[0], fine_x.max() + 1.6),
                fontsize=7.4, color=COL["small"], weight="bold", ha="center")
    ax.set(xlabel="time within the up pass [s]", ylabel="position [mm]")
    ax.set_title(f"(c) Up pass, targets {FINE_FIRST - 25}–{FINE_LAST + 25}",
                 loc="left", pad=6)
    ax.grid(color=COL["grid"], lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.035, 0.955,
            f"each tread = {c['settle_s']:.1f} s settle + {c['record_s']:.1f} s record",
            transform=ax.transAxes, va="top", fontsize=6.8, color="#566069")

    # (d) the fine zoom: individual risers, labelled in mm and in grid points.
    # Every riser is labelled, the single-grid-point one included: at 40 µm
    # against a 2.4 mm span it is nearly invisible, which is the point.
    ax = fig.add_subplot(gs[2, 3:6])
    base = fine_x[0]
    ax.plot(fine_t - fine_t[0], (fine_x - base) * 1000.0, color=COL["UP"], lw=1.25,
            solid_joinstyle="miter")
    for j in range(len(fine_stops) - 1):
        x_lo = (fine_stops[j][2] - base) * 1000.0
        x_hi = (fine_stops[j + 1][2] - base) * 1000.0
        grid_points = int(round((x_hi - x_lo) / 40.0))
        when = fine_stops[j + 1][0] - fine_t[0]
        smallest = grid_points == 1
        # A 1 gp riser is thinner than the line, so its label clears the trace.
        ax.text(when + (-0.4 if smallest else 0.5),
                (x_lo + x_hi) / 2.0 + (190.0 if smallest else 0.0),
                f"{grid_points * 0.04:.2f} mm · {grid_points} gp",
                fontsize=6.1, va="bottom" if smallest else "center",
                ha="right" if smallest else "left",
                color=COL["small"] if smallest else "#4b5257",
                weight="bold" if smallest else "normal",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.8)
                if smallest else None)
    ax.set(xlabel="time within the up pass [s]",
           ylabel=f"position − {base:.2f} mm [µm]")
    ax.set_title(f"(d) Nine stops, targets {FINE_FIRST}–{FINE_LAST}",
                 loc="left", pad=6)
    span = fine_t[-1] - fine_t[0]
    ax.set_xlim(-span * 0.03, span * 1.34)
    ax.margins(y=0.12)
    ax.grid(color=COL["grid"], lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    up_positions = np.array([s[2] for s in stops if s[5] == "up" and s[3] == "UP"])
    single = int((np.round(np.diff(up_positions) / 0.04).astype(int) == 1).sum())
    ax.text(0.035, 0.965,
            "gp = grid point = 64 µsteps = 0.040 mm\n"
            "     = one motor electrical cycle\n"
            f"the 1 gp step is one of {single} in the 400-target list",
            transform=ax.transAxes, va="top", fontsize=6.4, color="#566069",
            linespacing=1.45)

    ax2 = fig.add_subplot(gs[3, 0:6], sharex=ax0)
    n_records = summary["analysis_stops"] + summary["reference_stops"]
    components = [
        ("motion + transits", total - n_records * (c["settle_s"] + c["record_s"]), COL["motion"]),
        ("settle windows", n_records * c["settle_s"], COL["settle"]),
        ("record windows", n_records * c["record_s"], COL["record"]),
    ]
    left = 0.0
    for label, seconds, color in components:
        width = seconds / 60.0
        ax2.barh(0, width, left=left, height=0.48, color=color, edgecolor="white")
        if width > 2.0:
            ax2.text(left + width / 2.0, 0, f"{label}\n{width:.1f} min", ha="center",
                     va="center", fontsize=6.7, color="white", weight="bold")
        left += width
    ax2.set(yticks=[], xlabel="elapsed experiment time [min]", ylim=(-0.42, 0.42))
    ax2.set_title("(e) Recording-time budget", loc="left", pad=6)
    ax2.spines[["left", "right", "top"]].set_visible(False)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OVERVIEW, dpi=320)
    plt.close(fig)


def main() -> int:
    config, summary, plan = load()
    style()
    segments, blocks, stops, total = schedule(plan, config)
    if abs(total / 60.0 - summary["estimated_recording_min"]) > 1e-6:
        raise RuntimeError("render schedule disagrees with generated plan summary")
    if len(stops) != summary["analysis_stops"] + summary["reference_stops"]:
        raise RuntimeError("render replay recovered the wrong number of recorded stops")
    render_overview(config, summary, segments, blocks, stops, total)
    print(f"wrote {OVERVIEW.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

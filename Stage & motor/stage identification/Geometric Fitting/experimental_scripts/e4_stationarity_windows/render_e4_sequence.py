#!/usr/bin/env python3
"""Render publication-style figures from the executable E4 plan."""
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
OVERVIEW = OUT / "e4_stationarity_windows_sequence.png"
DETAIL = OUT / "e4_window_design_and_rmt.png"

WINDOW_COLORS = [
    "#1b6f8a", "#2b7f72", "#6a8052", "#a27832",
    "#a85b45", "#875b83", "#586a91",
]
COL = {
    "INIT": "#69747c", "DATUM": "#73578f", "REF": "#a16a2d",
    "END_REF": "#73578f", "PARK": "#4d5961", "record": "#67527c",
    "settle": "#c59b50", "motion": "#697b86",
}


def load():
    config = json.loads((ROOT / "e4_config.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "sequence/e4_plan_summary.json").read_text(encoding="utf-8"))
    with (ROOT / "sequence/e4_measurement_plan.csv").open(newline="", encoding="utf-8") as fh:
        plan = list(csv.DictReader(fh))
    with (ROOT / "sequence/e4_targets_grid192.csv").open(newline="", encoding="utf-8") as fh:
        targets = list(csv.DictReader(fh))
    return config, summary, plan, targets


def move_time(distance_mm: float, c: dict) -> float:
    distance = abs(distance_mm)
    velocity, acceleration = c["velocity_mm_s"], c["acceleration_mm_s2"]
    ramp_distance = velocity * velocity / acceleration
    if distance <= ramp_distance:
        return 2.0 * math.sqrt(distance / acceleration)
    return 2.0 * velocity / acceleration + (distance - ramp_distance) / velocity


def schedule(plan: list[dict], c: dict):
    now, position, state = 0.0, c["x_ref_mm"], "INIT"
    segments, blocks = [], []
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
        if settle:
            segments.append((now, now + settle, position, position, state, "settle"))
            now += settle
        if record:
            segments.append((now, now + record, position, position, state, "record"))
            now += record
    if now > block_start:
        blocks.append((block_start, now, state))
    return segments, blocks, now


def state_color(state: str) -> str:
    if state.startswith("W") and state[1:].isdigit():
        return WINDOW_COLORS[int(state[1:]) - 1]
    return COL.get(state, "#7a7f84")


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


def render_overview(c: dict, summary: dict, segments: list, blocks: list, total: float) -> None:
    total_min = total / 60.0
    fig = plt.figure(figsize=(7.45, 6.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.78, 3.15, 0.95],
                          left=0.105, right=0.97, top=0.86, bottom=0.10, hspace=0.31)
    fig.text(0.105, 0.975, "Experiment E4 — stationarity windows", fontsize=16,
             weight="bold", color="#20262b", va="top")
    fig.text(0.105, 0.928,
             "Seven local screw-revolution tests · 50 unique phases per window · bidirectional stop–settle–record",
             fontsize=8.5, color="#566069", va="top")

    ax0 = fig.add_subplot(gs[0])
    for start, end, state in blocks:
        width = (end - start) / 60.0
        ax0.barh(0, width, left=start / 60.0, height=0.40,
                 color=state_color(state), edgecolor="white", linewidth=0.55)
        if width > 0.38:
            label = state.replace("END_REF", "END")
            ax0.text((start + end) / 120.0, 0, label, ha="center", va="center",
                     fontsize=6.1, color="white", weight="bold")
    ax0.set(xlim=(0, total_min), ylim=(-0.34, 0.34), yticks=[])
    ax0.set_title("(a) Controller-state allocation", loc="left", pad=6)
    ax0.tick_params(axis="x", bottom=False, labelbottom=False)
    ax0.spines[:].set_visible(False)

    ax = fig.add_subplot(gs[1], sharex=ax0)
    states = ["DATUM", "W1", "REF", "W2", "W3", "W4", "W5", "W6", "W7", "END_REF", "PARK"]
    for state in states:
        lines = [[(a / 60.0, x0), (b / 60.0, x1)]
                 for a, b, x0, x1, st, _ in segments if st == state]
        if lines:
            ax.add_collection(LineCollection(lines, colors=state_color(state), linewidths=1.0))
    ax.axhline(c["x_ref_mm"], color=COL["DATUM"], lw=0.7, ls=(0, (3, 2)))
    ax.text(0.35, c["x_ref_mm"] + 2.0, r"$x_{ref}=61.35$ mm",
            color=COL["DATUM"], fontsize=7)
    for i, centre in enumerate(c["window_centres_mm"]):
        ax.axhspan(centre - 3.0, centre + 3.0, color=WINDOW_COLORS[i], alpha=0.035)
    ax.set(xlim=(0, total_min), ylim=(0, 121), ylabel="commanded position [mm]")
    ax.set_title("(b) Complete commanded trajectory", loc="left", pad=6)
    ax.grid(axis="x", color="#dde1e4", lw=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    box = (f"windows  {summary['windows']} × 6 mm\n"
           f"targets  {summary['targets_per_window']} per window\n"
           f"analysis stops  {summary['analysis_stops']}\n"
           f"reference records  {summary['reference_stops']}\n"
           f"nominal duration  {total_min:.1f} min")
    ax.text(0.985, 0.97, box, transform=ax.transAxes, ha="right", va="top",
            fontsize=7.4, linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.48", fc="white", ec="#aeb4b9", lw=0.7))

    ax2 = fig.add_subplot(gs[2], sharex=ax0)
    n_records = summary["total_record_stops"]
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
    ax2.axvline(c["watchdog_projection_min"], color="#b4473f", lw=0.9, ls="--")
    ax2.text(c["watchdog_projection_min"] - 0.25, -0.37, "58 min projected-finish gate",
             ha="right", va="top", fontsize=6.5, color="#8d3934")
    ax2.set(yticks=[], xlabel="elapsed experiment time [min]", ylim=(-0.42, 0.42))
    ax2.set_title("(c) Recording-time budget", loc="left", pad=6)
    ax2.spines[["left", "right", "top"]].set_visible(False)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OVERVIEW, dpi=320)
    plt.close(fig)


def render_detail(c: dict, summary: dict, targets: list[dict]) -> None:
    grouped = {}
    for row in targets:
        grouped.setdefault(row["window"], []).append(row)

    fig, axes = plt.subplots(2, 2, figsize=(8.27, 7.25))
    fig.subplots_adjust(left=0.09, right=0.965, top=0.875, bottom=0.10,
                        hspace=0.43, wspace=0.30)
    fig.text(0.09, 0.975, "E4 local phase coverage and command fidelity", fontsize=16,
             weight="bold", color="#20262b", va="top")
    fig.text(0.09, 0.938,
             "Executable design: 6 mm = 3 screw revolutions = 50 intervals on the 0.12 mm grid",
             fontsize=8.6, color="#566069", va="top")

    ax = axes[0, 0]
    ax.hlines(0, c["soft_min_mm"], c["soft_max_mm"], color="#8d979e", lw=1.0)
    for i, centre in enumerate(c["window_centres_mm"]):
        rows = grouped[f"W{i + 1}"]
        x = np.array([float(row["target_mm"]) for row in rows])
        ax.axvspan(x[0], x[-1], ymin=0.31, ymax=0.69, color=WINDOW_COLORS[i], alpha=0.78)
        ax.vlines(x, -0.13, 0.13, color=WINDOW_COLORS[i], lw=0.34, alpha=0.9)
        ax.text(centre, 0.34, f"W{i + 1}\n{centre:.0f}", ha="center", va="bottom",
                fontsize=7, color=WINDOW_COLORS[i], weight="bold")
    ax.plot(c["x_ref_mm"], -0.31, marker="D", ms=4.2, color=COL["DATUM"])
    ax.text(c["x_ref_mm"] + 1.6, -0.31, r"$x_{ref}$", va="center", fontsize=7,
            color=COL["DATUM"])
    ax.set(xlim=(0, 122), ylim=(-0.46, 0.66), yticks=[], xlabel="commanded position [mm]")
    ax.set_title("(a) Seven evenly distributed local windows", loc="left", pad=7)
    ax.spines[["left", "right", "top"]].set_visible(False)

    ax = axes[0, 1]
    indices = np.arange(1, 51)
    for i in range(7):
        phase = np.array(summary["phase_bins_by_window"][f"W{i + 1}"][:50])
        ax.scatter(indices, np.full(50, i + 1), c=phase, cmap="twilight_shifted",
                   vmin=0, vmax=49, s=15, marker="s", linewidths=0)
        endpoint_phase = summary["phase_bins_by_window"][f"W{i + 1}"][-1]
        ax.scatter([51], [i + 1], c=[endpoint_phase], cmap="twilight_shifted",
                   vmin=0, vmax=49, s=24, marker="D", edgecolor="#30363b", linewidth=0.35)
    ax.axvline(50.5, color="#7e858a", lw=0.7, ls=(0, (3, 2)))
    ax.text(50.7, 7.42, "repeated\nendpoint", fontsize=6.5, color="#60676c", ha="center")
    ax.set(xlim=(0, 52), ylim=(0.4, 7.75), yticks=np.arange(1, 8),
           yticklabels=[f"W{i}" for i in range(1, 8)], xlabel="target index")
    ax.set_title("(b) All 50 screw-phase bins once per window", loc="left", pad=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e3e5e7", lw=0.45)

    ax = axes[1, 0]
    lower, upper = -3.0, 9.0
    up = np.linspace(0, 6, 51)
    order_x = np.concatenate(([lower], up, [upper], up[::-1]))
    order_n = np.arange(order_x.size)
    ax.plot(order_n, order_x, color=WINDOW_COLORS[0], lw=1.15)
    ax.scatter(np.arange(1, 52), up, s=6, color="#1b6f8a", label="ascending records")
    ax.scatter(np.arange(53, 104), up[::-1], s=6, color="#a85b45", label="descending records")
    ax.scatter([0, 52], [lower, upper], marker="x", s=28, color="#444b50", label="approach / overrun")
    ax.axhspan(0, 6, color="#1b6f8a", alpha=0.035)
    ax.set(xlim=(-2, 105), ylim=(-3.8, 9.8), xlabel="command order within one window",
           ylabel="position relative to window start [mm]")
    ax.set_title("(c) Identical bidirectional traversal in each window", loc="left", pad=7)
    ax.legend(frameon=False, fontsize=6.7, loc="upper center", ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#e1e3e5", lw=0.5)

    ax = axes[1, 1]
    pulses = c["grid_spacing_usteps"]
    distance = c["grid_spacing_mm"]
    velocity, acceleration = c["velocity_mm_s"], c["acceleration_mm_s2"]
    if distance <= velocity * velocity / acceleration:
        t_acc = math.sqrt(distance / acceleration)
        t_cruise = 0.0
        peak = math.sqrt(distance * acceleration)
    else:
        t_acc = velocity / acceleration
        t_cruise = (distance - velocity * velocity / acceleration) / velocity
        peak = velocity
    t_end = 2.0 * t_acc + t_cruise
    time = np.linspace(0, t_end, 500)
    speed = np.where(time <= t_acc, acceleration * time,
            np.where(time <= t_acc + t_cruise, peak, acceleration * (t_end - time)))
    position = np.zeros_like(time)
    for j in range(1, len(time)):
        position[j] = position[j - 1] + 0.5 * (speed[j - 1] + speed[j]) * (time[j] - time[j - 1])
    ax.fill_between(time * 1000, speed, color="#1b6f8a", alpha=0.15)
    ax.plot(time * 1000, speed, color="#1b6f8a", lw=1.6)
    axr = ax.twinx()
    axr.plot(time * 1000, position, color="#67527c", lw=1.05, ls=":")
    axr.set_ylabel("displacement [mm]", color="#67527c")
    axr.set_ylim(0, distance * 1.12)
    ax.text(0.03, 0.96,
            "adjacent target move  0.120 mm = 192 STEP pulses\n"
            "RMT transactions  128 + 64 symbols\n"
            f"then SETTLE {c['settle_s']:.2f} s + RECORD {c['record_s']:.2f} s",
            transform=ax.transAxes, va="top", fontsize=7.1, linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.38", fc="white", ec="#b8bdc1", lw=0.6))
    ax.set(xlabel="time from move command [ms]", ylabel="commanded velocity [mm/s]")
    ax.set_title("(d) Representative buffered RMT move", loc="left", pad=7)
    ax.spines["top"].set_visible(False)
    axr.spines["top"].set_visible(False)
    ax.grid(color="#e1e3e5", lw=0.5)

    fig.text(0.09, 0.035,
             r"Per-window inference: $e_w(\theta)=a_{0,w}+\sum_h[c_{h,w}\cos(h\theta)+s_{h,w}\sin(h\theta)]$.  "
             "Compare amplitude and unwrapped absolute phase across W1–W7; references monitor drift and slip.",
             fontsize=7.1, color="#5b6268")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(DETAIL, dpi=320)
    plt.close(fig)


def main() -> int:
    config, summary, plan, targets = load()
    style()
    segments, blocks, total = schedule(plan, config)
    if abs(total / 60.0 - summary["estimated_recording_min"]) > 1e-6:
        raise RuntimeError("render schedule disagrees with generated plan summary")
    render_overview(config, summary, segments, blocks, total)
    render_detail(config, summary, targets)
    print(f"wrote {OVERVIEW.relative_to(ROOT)}")
    print(f"wrote {DETAIL.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

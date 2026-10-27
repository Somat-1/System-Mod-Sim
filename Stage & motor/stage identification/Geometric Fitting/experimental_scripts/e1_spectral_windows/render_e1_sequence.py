#!/usr/bin/env python3
"""Render the generated E1 measurement sequence as an A4 landscape PNG."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

from e1_sequence import (
    DEFAULT_CONFIG,
    build_measurement_plan,
    generate_targets,
    load_config,
    move_time_s,
    validate_config,
    validate_plan,
)


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "rendered" / "e1_spectral_windows_sequence.png"
DETAIL_OUTPUT = ROOT / "rendered" / "e1_window_rmt_detail.png"

COLORS = {
    "INIT": "#68717a",
    "DATUM": "#7a6f91",
    "PILOT": "#a9783f",
    "W1": "#4c78a8",
    "W2": "#598a76",
    "W3": "#9d6770",
    "REF": "#817765",
    "END_REF": "#817765",
    "PARK": "#68717a",
    "DONE": "#68717a",
}


def schedule(rows, config):
    time_s = 0.0
    position = config["x_ref_mm"]
    segments = []
    blocks = []
    active_state = None
    block_start = 0.0
    block_label = ""
    ref_number = 0

    def change_block(label, state):
        nonlocal active_state, block_start, block_label
        if active_state is not None and time_s > block_start:
            blocks.append((block_start, time_s, block_label, active_state))
        active_state = state
        block_start = time_s
        block_label = label

    for row in rows:
        if row.action == "state_enter":
            label = row.state
            if row.state == "REF":
                ref_number += 1
                label = f"REF {ref_number}"
            change_block(label, row.state)
            continue
        if row.target_mm is not None:
            duration = move_time_s(row.target_mm - position, config)
            if duration:
                segments.append((time_s, time_s + duration, position,
                                 row.target_mm, row.state, "move"))
                time_s += duration
            position = row.target_mm
        if row.action == "stop":
            settle_end = time_s + row.settle_s
            segments.append((time_s, settle_end, position, position,
                             row.state, "settle"))
            time_s = settle_end
            record_end = time_s + row.record_s
            segments.append((time_s, record_end, position, position,
                             row.state, "record"))
            time_s = record_end
        elif row.action == "operator_gate":
            pause_end = time_s + config["pilot_pause_budget_s"]
            segments.append((time_s, pause_end, position, position,
                             "PILOT", "operator"))
            time_s = pause_end
    if active_state is not None and time_s > block_start:
        blocks.append((block_start, time_s, block_label, active_state))
    return segments, blocks, time_s


def render_window_detail(config, targets):
    """Explain one window path and one adjacent-target RMT command."""
    count = len(targets["W1"])
    spacing = config["nominal_spacing_mm"]
    window_span = config["windows_mm"][0][1] - config["windows_mm"][0][0]
    overrun = config["overrun_mm"]
    adjacent_moves = 2 * (count - 1)
    long_moves = 3  # pre->first, last->overrun, overrun->last
    path_mm = adjacent_moves * spacing + long_moves * overrun
    move_012_s = move_time_s(spacing, config)
    move_3_s = move_time_s(overrun, config)
    motion_s = adjacent_moves * move_012_s + long_moves * move_3_s
    dwell_s = 2 * count * (config["settle_s"] + config["record_s"])

    fig = plt.figure(figsize=(11.69, 8.27), facecolor="#f5f3ef")
    outer = fig.add_gridspec(
        2, 2, height_ratios=[1.52, 1], width_ratios=[1.18, 1],
        left=0.07, right=0.975, top=0.84, bottom=0.13, hspace=0.34, wspace=0.25,
    )
    fig.text(0.07, 0.985, "Inside one E1 spectral window",
             fontsize=20, weight="bold", color="#20242a", va="top")
    fig.text(
        0.07, 0.935,
        "Regular 0.12 mm target grid · bidirectional measurement · RMT-generated STEP pulses",
        fontsize=10, color="#4c5560", va="top",
    )

    # Complete path, parameterized by cumulative distance travelled.
    path_ax = fig.add_subplot(outer[0, :])
    path_ax.set_facecolor("#fffdfa")
    up_positions = [i * spacing for i in range(count)]
    up_travel = [overrun + x for x in up_positions]
    down_positions = list(reversed(up_positions))
    down_travel = [overrun + window_span + 2 * overrun +
                   (window_span - x) for x in down_positions]
    path_ax.axhspan(0, window_span, color=COLORS["W1"], alpha=0.07)
    path_ax.plot([0, overrun], [-overrun, 0], color="#65717d", lw=2.0)
    path_ax.plot(up_travel, up_positions, color=COLORS["W1"], lw=2.2,
                 label="upward pass")
    path_ax.plot([overrun + window_span, overrun + window_span + overrun],
                 [window_span, window_span + overrun], color="#65717d", lw=2.0)
    path_ax.plot([overrun + window_span + overrun,
                  overrun + window_span + 2 * overrun],
                 [window_span + overrun, window_span], color="#65717d", lw=2.0)
    path_ax.plot(down_travel, down_positions, color=COLORS["W3"], lw=2.2,
                 label="downward pass")
    path_ax.scatter(up_travel, up_positions, s=8, color=COLORS["W1"], zorder=3)
    path_ax.scatter(down_travel, down_positions, s=8, color=COLORS["W3"], zorder=3)
    path_ax.axhline(0, color="#858d94", lw=0.8, ls="--")
    path_ax.axhline(window_span, color="#858d94", lw=0.8, ls="--")
    path_ax.axhline(-overrun, color="#aeb3b7", lw=0.7, ls=":")
    path_ax.axhline(window_span + overrun, color="#aeb3b7", lw=0.7, ls=":")
    path_ax.annotate("start: pre-approach\n−3 mm", xy=(0, -overrun), xytext=(1.0, -1.1),
                     arrowprops=dict(arrowstyle="->", color="#59636f"),
                     fontsize=8.5, color="#3f4851")
    path_ax.annotate("reverse outside measured range\n+3 mm overrun", 
                     xy=(overrun + window_span + overrun, window_span + overrun),
                     xytext=(22.2, 14.4), ha="center",
                     arrowprops=dict(arrowstyle="->", color="#59636f"),
                     fontsize=8.5, color="#3f4851")
    path_ax.annotate("finish at first target", xy=(path_mm, 0), xytext=(29.0, 2.2),
                     arrowprops=dict(arrowstyle="->", color=COLORS["W3"]),
                     fontsize=8.5, color=COLORS["W3"], ha="center")
    path_ax.text(3.35, 11.25, "101 recorded targets (up)", color=COLORS["W1"],
                 fontsize=9, weight="bold")
    path_ax.text(24.1, 5.5, "101 recorded targets (down)", color=COLORS["W3"],
                 fontsize=9, weight="bold", rotation=-35)
    path_ax.set_xlim(-0.4, path_mm + 0.5)
    path_ax.set_ylim(-4.4, window_span + overrun + 1.4)
    path_ax.set_xlabel("cumulative distance travelled from the pre-approach point [mm]")
    path_ax.set_ylabel("position relative to window start [mm]")
    path_ax.set_title("(a) Complete path through one 12 mm window", loc="left",
                      fontsize=12, weight="bold", pad=8)
    path_ax.grid(color="#d9dcde", lw=0.6)
    path_ax.legend(loc="lower center", ncol=2, frameon=False, fontsize=8.5)
    path_ax.text(
        0.985, 0.055,
        f"window span  {window_span:.0f} mm\n"
        f"total path  {path_mm:.0f} mm\n"
        f"recorded stops  {2 * count}\n"
        f"0.12 mm moves  {adjacent_moves}\n"
        f"3 mm moves  {long_moves}",
        transform=path_ax.transAxes, ha="right", va="bottom", fontsize=8.5,
        color="#30373d", bbox=dict(boxstyle="round,pad=0.45", fc="#f4f1ea",
                                    ec="#c9c5bb", lw=0.8),
    )

    # Ideal continuous profile for one adjacent target move.
    profile_ax = fig.add_subplot(outer[1, 0])
    profile_ax.set_facecolor("#fffdfa")
    velocity = config["velocity_mm_s"]
    accel = config["acceleration_mm_s2"]
    t_acc = velocity / accel
    d_acc = 0.5 * accel * t_acc**2
    t_cruise = (spacing - 2 * d_acc) / velocity
    t_total = 2 * t_acc + t_cruise
    ts = [t_total * i / 400 for i in range(401)]
    vs, xs = [], []
    for t in ts:
        if t <= t_acc:
            vs.append(accel * t)
            xs.append(0.5 * accel * t * t)
        elif t <= t_acc + t_cruise:
            vs.append(velocity)
            xs.append(d_acc + velocity * (t - t_acc))
        else:
            td = t - t_acc - t_cruise
            vs.append(velocity - accel * td)
            xs.append(d_acc + velocity * t_cruise + velocity * td - 0.5 * accel * td * td)
    time_ms = [t * 1000 for t in ts]
    profile_ax.fill_between(time_ms, vs, color=COLORS["W1"], alpha=0.16)
    profile_ax.plot(time_ms, vs, color=COLORS["W1"], lw=2.1)
    for boundary in (t_acc * 1000, (t_acc + t_cruise) * 1000, t_total * 1000):
        profile_ax.axvline(boundary, color="#8c949b", lw=0.8, ls="--")
    profile_ax.text(25, 1.05, "80 pulses\nramp up", ha="center", va="center",
                    fontsize=8.5, color=COLORS["W1"], weight="bold")
    profile_ax.text(55, 1.78, "32 pulses\nconstant", ha="center", va="center",
                    fontsize=8, color="#3e765f", weight="bold")
    profile_ax.text(85, 1.05, "80 pulses\nramp down", ha="center", va="center",
                    fontsize=8.5, color=COLORS["W3"], weight="bold")
    profile_ax.set_xlim(0, t_total * 1000)
    profile_ax.set_ylim(0, 2.25)
    profile_ax.set_xlabel("time from move command [ms]")
    profile_ax.set_ylabel("commanded velocity [mm/s]", color=COLORS["W1"])
    profile_ax.set_title("(b) One adjacent-target move: 0.12 mm", loc="left",
                         fontsize=12, weight="bold", pad=8)
    profile_ax.grid(color="#d9dcde", lw=0.6)
    distance_ax = profile_ax.twinx()
    distance_ax.plot(time_ms, xs, color="#544766", lw=1.5, ls=":")
    distance_ax.set_ylim(0, spacing * 1.12)
    distance_ax.set_ylabel("displacement [mm]", color="#544766")
    profile_ax.text(
        0.02, 0.96,
        "50 ms acceleration  +  10 ms cruise  +  50 ms deceleration",
        transform=profile_ax.transAxes, va="top", fontsize=8.2, color="#374047",
    )

    # One stop and the exact RMT buffer split.
    detail_grid = outer[1, 1].subgridspec(2, 1, hspace=1.05)
    stop_ax = fig.add_subplot(detail_grid[0])
    stop_ax.set_facecolor("#fffdfa")
    durations = [move_012_s, config["settle_s"], config["record_s"]]
    labels = ["", "SETTLE\n2.00 s", "RECORD\n0.50 s"]
    colors = [COLORS["W1"], "#b39a63", "#735a9a"]
    left = 0.0
    for duration, label, color in zip(durations, labels, colors):
        stop_ax.barh(0, duration, left=left, height=0.55, color=color,
                     edgecolor="white", linewidth=1.0)
        stop_ax.text(left + duration / 2, 0, label, ha="center", va="center",
                     color="white", fontsize=8.2, weight="bold")
        left += duration
    stop_ax.set_xlim(0, sum(durations))
    stop_ax.set_ylim(-0.42, 0.62)
    stop_ax.set_yticks([])
    stop_ax.annotate("MOVE  0.11 s", xy=(move_012_s / 2, 0.23),
                     xytext=(0.34, 0.47), ha="center", fontsize=7.8,
                     color=COLORS["W1"], weight="bold",
                     arrowprops=dict(arrowstyle="->", color=COLORS["W1"], lw=0.8))
    stop_ax.set_title("(c) One recorded stop: 2.61 s", loc="left",
                      fontsize=11.5, weight="bold", pad=7)
    stop_ax.spines[["left", "right", "top"]].set_visible(False)

    buffer_ax = fig.add_subplot(detail_grid[1])
    buffer_ax.set_facecolor("#fffdfa")
    buffer_ax.barh(0, 80, left=0, height=0.52, color=COLORS["W1"])
    buffer_ax.barh(0, 32, left=80, height=0.52, color="#2d8b70")
    buffer_ax.barh(0, 80, left=112, height=0.52, color=COLORS["W3"])
    buffer_ax.axvline(config["rmt_chunk_pulses"], color="#20242a", lw=1.5, ls="--")
    buffer_ax.text(40, 0, "80 accel", ha="center", va="center", color="white",
                   fontsize=8, weight="bold")
    buffer_ax.text(96, 0, "32\ncruise", ha="center", va="center", color="white",
                   fontsize=7.3, weight="bold")
    buffer_ax.text(152, 0, "80 decel", ha="center", va="center", color="white",
                   fontsize=8, weight="bold")
    buffer_ax.text(64, 0.46, "RMT transaction 1 · 128 symbols", ha="center",
                   fontsize=7.8, color="#3e474f")
    buffer_ax.text(160, 0.46, "transaction 2 · 64 symbols", ha="center",
                   fontsize=7.8, color="#3e474f")
    buffer_ax.set_xlim(0, config["grid_spacing_usteps"])
    buffer_ax.set_ylim(-0.45, 0.72)
    buffer_ax.set_yticks([])
    buffer_ax.set_xlabel("STEP pulse / microstep index")
    buffer_ax.set_title("(d) One logical RMT train: 192 STEP pulses", loc="left",
                        fontsize=11.5, weight="bold", pad=10)
    buffer_ax.spines[["left", "right", "top"]].set_visible(False)

    fig.text(
        0.07, 0.035,
        f"One window: {motion_s:.1f} s commanded motion + {dwell_s:.1f} s settling/recording"
        f" = {(motion_s + dwell_s) / 60:.2f} min, excluding the transit from x_ref to the pre-approach point.",
        fontsize=8.5, color="#59616a",
    )
    DETAIL_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(DETAIL_OUTPUT, dpi=220, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> int:
    config, _ = load_config(DEFAULT_CONFIG)
    validate_config(config)
    targets = generate_targets(config)
    rows = build_measurement_plan(config, targets)
    summary = validate_plan(config, targets, rows)
    segments, blocks, total_s = schedule(rows, config)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.edgecolor": "#8a9096",
        "axes.labelcolor": "#2f3438",
        "xtick.color": "#454b50",
        "ytick.color": "#454b50",
        "legend.fontsize": 7,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })
    duration_min = total_s / 60.0
    fig = plt.figure(figsize=(7.45, 5.25), facecolor="white")
    grid = fig.add_gridspec(
        2, 1, height_ratios=[0.82, 3.15],
        left=0.105, right=0.97, top=0.945, bottom=0.115, hspace=0.20,
    )
    timeline = fig.add_subplot(grid[0])
    timeline.set_xlim(0, duration_min)
    timeline.set_ylim(-0.40, 0.76)
    for start, end, label, state in blocks:
        if end <= start:
            continue
        x, width = start / 60.0, (end - start) / 60.0
        timeline.barh(
            0.15, width, left=x, height=0.38,
            color=COLORS.get(state, "#777777"), edgecolor="white", linewidth=0.6,
        )
        short_label = label.replace("REF ", "R").replace("END_REF", "END")
        if width >= 1.0:
            timeline.text(
                x + width / 2, 0.15, short_label,
                ha="center", va="center", fontsize=6.5,
                color="white", weight="bold", clip_on=True,
            )
        elif state == "DATUM":
            timeline.plot([x + width / 2, x + width / 2], [-0.05, -0.13],
                          color="#555b60", lw=0.55)
            timeline.text(x + width / 2 + 0.12, -0.16, short_label,
                          ha="left", va="top", fontsize=6.2, color="#3f454a")
        else:
            timeline.plot([x + width / 2, x + width / 2], [0.35, 0.48],
                          color="#555b60", lw=0.55)
            timeline.text(x + width / 2, 0.51, short_label, ha="center", va="bottom",
                          fontsize=6.2, color="#3f454a")
    timeline.vlines(0, -0.34, 0.34, color=COLORS["INIT"], lw=0.8)
    timeline.text(0.18, -0.32, "INIT", ha="left", va="center",
                  fontsize=6.2, color="#3f454a")
    timeline.vlines(duration_min, -0.34, 0.34, color=COLORS["PARK"], lw=0.8)
    timeline.text(duration_min - 0.18, -0.18, "PARK / DONE", ha="right", va="center",
                  fontsize=6.2, color="#3f454a")
    timeline.set_title("(a) Time-aligned sequence states", loc="left", pad=7)
    timeline.set_yticks([])
    timeline.tick_params(axis="x", bottom=False, labelbottom=False)
    timeline.spines[["left", "right", "top", "bottom"]].set_visible(False)

    ax = fig.add_subplot(grid[1], sharex=timeline)
    ax.set_facecolor("white")
    line_styles = {
        "DATUM": (0, (2.5, 1.5)),
        "PILOT": (0, (1.2, 1.2)),
        "W1": "solid",
        "W2": (0, (4, 1.6)),
        "W3": (0, (5, 1.5, 1.2, 1.5)),
        "REF": (0, (2.5, 1.5)),
        "END_REF": (0, (2.5, 1.5)),
    }
    for state in line_styles:
        lines = [[(a / 60.0, x0), (b / 60.0, x1)]
                 for a, b, x0, x1, segment_state, _ in segments
                 if segment_state == state]
        if lines:
            ax.add_collection(LineCollection(
                lines, colors=COLORS[state], linewidths=1.05,
                linestyles=line_styles[state], alpha=0.98,
            ))

    for low, high in config["windows_mm"]:
        ax.hlines((low, high), 0, duration_min, color="#c9cdd1",
                  lw=0.5, ls=(0, (2, 3)), zorder=0)
    ax.axhline(config["x_ref_mm"], color=COLORS["REF"], lw=0.65,
               ls=(0, (3, 2)), alpha=0.8)
    ax.text(10.3, config["x_ref_mm"] + 1.4, "x_ref",
            color=COLORS["REF"], fontsize=6.5, ha="center")

    ax.axvline(duration_min, color="#777d82", lw=0.65, ls=(0, (2, 2)))

    metrics = f"Measurement duration: {duration_min:.1f} min"
    ax.text(0.018, 0.972, metrics, transform=ax.transAxes,
            ha="left", va="top", fontsize=7,
            color="#33383d", linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.45", facecolor="white",
                      edgecolor="#aeb3b7", linewidth=0.6))

    legend_handles = [
        Line2D([0], [0], color=COLORS["PILOT"], lw=1.2,
               ls=line_styles["PILOT"], label="pilot"),
        Line2D([0], [0], color=COLORS["W1"], lw=1.2,
               ls=line_styles["W1"], label="W1"),
        Line2D([0], [0], color=COLORS["W2"], lw=1.2,
               ls=line_styles["W2"], label="W2"),
        Line2D([0], [0], color=COLORS["W3"], lw=1.2,
               ls=line_styles["W3"], label="W3"),
        Line2D([0], [0], color=COLORS["REF"], lw=1.2,
               ls=line_styles["REF"], label="reference/datum"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", ncol=2,
              frameon=True, fancybox=False, framealpha=1.0,
              edgecolor="#aeb3b7", facecolor="white",
              borderpad=0.5, handlelength=2.4, columnspacing=1.2)

    ax.set_xlim(0, duration_min)
    x_ticks = list(range(0, 31, 5)) + [duration_min]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(value) for value in range(0, 31, 5)] +
                       [f"{duration_min:.1f}"])
    ax.set_ylim(0, 123)
    ax.set_xlabel("elapsed recording time [min]")
    ax.set_ylabel("commanded position [mm]")
    ax.set_title("(b) Commanded trajectory", loc="left", pad=6)
    ax.grid(axis="x", color="#e0e2e4", lw=0.5)
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right"]].set_visible(False)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=320, facecolor="white")
    plt.close(fig)
    render_window_detail(config, targets)
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"wrote {DETAIL_OUTPUT.relative_to(ROOT)}")
    print(json.dumps({"estimated_min": total_s / 60.0,
                      "overview": str(OUTPUT), "detail": str(DETAIL_OUTPUT)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

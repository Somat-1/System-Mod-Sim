#!/usr/bin/env python3
"""Interactive ISO230-2 v3 raw, command, and corrected-trajectory viewer."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import render_command_recorded_overlay as overlay
from render_raw_measurement import read_scope

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / "ISO230-2_v3.csv"
TARGETS = HERE.parent / "sequence" / "iso_targets_seed2302_rom122.695.csv"


def build_plan():
    with TARGETS.open(newline="", encoding="utf-8") as handle:
        targets = np.asarray([float(row["cmd_mm"]) for row in csv.DictReader(handle)])

    xref, pre, overrun = 61.35, 58.35, 3.0
    speed, accel, dwell, sample_dt = 2.0, 40.0, 2.5, 0.02
    times, positions, stops = [0.0], [xref], []
    now, position = 0.0, xref

    def move(target):
        nonlocal now, position
        distance = abs(target - position)
        if distance == 0:
            return
        sign = 1.0 if target > position else -1.0
        critical = speed * speed / accel
        if distance <= critical:
            t_acc = np.sqrt(distance / accel)
            duration = 2.0 * t_acc
            peak = accel * t_acc
        else:
            t_acc = speed / accel
            duration = distance / speed + speed / accel
            peak = speed
        local = np.r_[np.arange(sample_dt, duration, sample_dt), duration]
        if distance <= critical:
            travelled = np.where(
                local <= t_acc,
                0.5 * accel * local**2,
                distance - 0.5 * accel * (duration - local)**2,
            )
        else:
            d_acc = 0.5 * accel * t_acc**2
            flat_end = duration - t_acc
            travelled = np.where(
                local <= t_acc,
                0.5 * accel * local**2,
                np.where(
                    local < flat_end,
                    d_acc + peak * (local - t_acc),
                    distance - 0.5 * accel * (duration - local)**2,
                ),
            )
        times.extend((now + local).tolist())
        positions.extend((position + sign * travelled).tolist())
        now += duration
        position = float(target)

    def hold(duration):
        nonlocal now
        now += duration
        times.append(now)
        positions.append(position)

    def record(state, run, direction, target_i):
        stops.append((now + dwell / 2.0, state, run, direction, target_i, position))
        hold(dwell)

    def visit_reference(state, run):
        move(pre)
        move(xref)
        record(state, run, 0, -1)

    for k in range(3):
        visit_reference("DATUM", k)
    for run in range(5):
        visit_reference("ISO_REF", run)
        move(targets[0] - overrun)
        for i, target in enumerate(targets, 1):
            move(target)
            record("ISO", run, +1, i)
        move(targets[-1] + overrun)
        for i, target in reversed(list(enumerate(targets, 1))):
            move(target)
            record("ISO", run, -1, i)
        move(targets[0] - overrun)
    for k in range(3):
        visit_reference("END_REF", k)
    move(pre)
    move(xref)
    return np.asarray(times), np.asarray(positions), stops


def main():
    scope = read_scope(SOURCE)
    plateaus = overlay.detect_plateaus(scope)
    plan_time, plan_position, stops = build_plan()
    command_stop = np.asarray([stop[5] for stop in stops])
    board_stop = np.asarray([stop[0] for stop in stops])

    datum = overlay.estimate_counter_datum(command_stop, plateaus.position_mm)
    plateau_absolute = plateaus.position_mm + datum
    alignment = overlay.align_sequences(command_stop, plateau_absolute)
    command_i = np.asarray([i for i, _ in alignment.pairs], dtype=int)
    plateau_i = np.asarray([j for _, j in alignment.pairs], dtype=int)
    rate, offset = overlay.fit_affine_clock(
        board_stop[command_i], plateaus.centre_s[plateau_i]
    )
    corrected_time = rate * plan_time + offset
    raw_absolute = scope.displacement_mm + datum

    overlay.style()
    fig, axes = plt.subplots(3, 1, figsize=(16.2, 9.5))
    raw_ax, command_ax, corrected_ax = axes
    fig.canvas.manager.set_window_title(
        "ISO230-2 v3 - complete commanded and recorded trajectory"
    )

    raw_ax.plot(scope.time_s / 60.0, raw_absolute, color="#155f8a", lw=0.68,
                label="raw interferometer", rasterized=True)
    raw_ax.set_title("(a) Complete raw interferometer recording", loc="left",
                     fontweight="bold")
    raw_ax.set_ylabel("Recorded position [mm]")
    raw_ax.set_xlim(0, scope.time_s[-1] / 60.0)
    raw_ax.legend(loc="upper right")

    command_ax.plot(plan_time / 60.0, plan_position, color="#4d4d4d", lw=0.76,
                    label="controller command", rasterized=True)
    command_ax.scatter(board_stop / 60.0, command_stop, s=6, color="#d55e00",
                       linewidths=0, label="recorded-stop command", zorder=3)
    command_ax.set_title("(b) Complete controller command from the flashed v3 plan",
                         loc="left", fontweight="bold")
    command_ax.set_ylabel("Intended position [mm]")
    command_ax.set_xlabel("Nominal controller sequence time [min]")
    command_ax.set_xlim(0, plan_time[-1] / 60.0)
    command_ax.legend(loc="upper right", ncol=2)

    corrected_ax.plot(scope.time_s / 60.0, raw_absolute, color="#155f8a", lw=0.76,
                      label="recorded position", rasterized=True)
    valid = (corrected_time >= 0) & (corrected_time <= scope.time_s[-1])
    corrected_ax.plot(corrected_time[valid] / 60.0, plan_position[valid],
                      color="#4d4d4d", lw=0.76, alpha=0.86,
                      linestyle=(0, (3, 2)), rasterized=True,
                      label="intended position - clock corrected")
    corrected_ax.scatter(plateaus.centre_s[plateau_i] / 60.0,
                         plateau_absolute[plateau_i], s=7, color="#009e73",
                         linewidths=0, label="paired stable arrival", zorder=4)
    corrected_ax.set_title(
        "(c) Clock-corrected command over the complete raw recording",
        loc="left", fontweight="bold")
    corrected_ax.set_ylabel("Position [mm]")
    corrected_ax.set_xlabel("Scope elapsed time [min]")
    corrected_ax.set_xlim(0, scope.time_s[-1] / 60.0)
    corrected_ax.legend(loc="upper right", ncol=2)

    for axis in axes:
        axis.set_ylim(-3, 123)
        axis.grid(True, color="0.9", lw=0.55)

    home_limits = [axis.get_xlim() + axis.get_ylim() for axis in axes]

    def rescale_visible(axis):
        left, right = sorted(axis.get_xlim())
        if axis is command_ax:
            xdata, ydata = plan_time / 60.0, plan_position
        else:
            xdata, ydata = scope.time_s / 60.0, raw_absolute
        use = (xdata >= left) & (xdata <= right)
        if np.any(use):
            low, high = float(np.min(ydata[use])), float(np.max(ydata[use]))
            pad = max(0.002, 0.08 * (high - low))
            axis.set_ylim(low - pad, high + pad)

    def on_key(event):
        key = (event.key or "").lower()
        if key == "r":
            for axis, limits in zip(axes, home_limits):
                axis.set_xlim(limits[0], limits[1])
                axis.set_ylim(limits[2], limits[3])
        elif key == "y":
            if event.inaxes in axes:
                rescale_visible(event.inaxes)
            else:
                for axis in axes:
                    rescale_visible(axis)
        elif key == "g":
            for axis in axes:
                axis.grid(not axis.xaxis._major_tick_kw.get("gridOn", False),
                          color="0.9", lw=0.55)
        else:
            return
        fig.canvas.draw_idle()

    def on_scroll(event):
        axis = event.inaxes
        if axis not in axes or event.xdata is None:
            return
        left, right = axis.get_xlim()
        factor = 0.75 if event.button == "up" else 1.0 / 0.75
        centre = event.xdata
        axis.set_xlim(centre - (centre - left) * factor,
                      centre + (right - centre) * factor)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("scroll_event", on_scroll)
    fig.suptitle("ISO 230-2 v3 - complete commanded and recorded trajectory",
                 x=0.067, y=0.985, ha="left", fontsize=15, fontweight="bold")
    fig.text(0.067, 0.015,
             "Toolbar: pan / box zoom / home   |   scroll: horizontal zoom   |   "
             "Y: fit vertical range to the visible window   |   R: reset   |   G: grid",
             ha="left", va="bottom", fontsize=8.5, color="0.28")
    fig.subplots_adjust(left=0.07, right=0.985, top=0.92, bottom=0.09, hspace=0.48)
    plt.show()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render the raw 120 ns AquadB v2/v5/v8 bandwidth experiment.

This first-stage renderer is intentionally descriptive. It unwraps the EL5101
UINT32 counter exactly, converts counts using the nominal 1 nm/count scale,
and identifies the twelve long experiment transfers from the raw trace. It
does not filter, interpolate, resample, or calculate a verdict.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import numpy as np
from scipy.signal import welch

from render_raw_measurement import read_scope


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE.parent.parent / "RMT_BW_Test.csv"
DEFAULT_LOG = HERE.parent / "data" / "terminal_20261023_155256.log"
DEFAULT_OUTPUT = HERE / "plots" / "aquadb_120ns"
MOTOR_SCRIPTS = HERE.parent.parent.parent / "motor identification" / "scripts"
if str(MOTOR_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(MOTOR_SCRIPTS))
from staircase_figures import INK, INK2, INK3, SURFACE, style  # noqa: E402


PROFILE_INFO = (
    ("v2", 2.0, "#3568a8"),
    ("v5", 5.0, "#b56732"),
    ("v8", 8.0, "#39705a"),
)
TRANSFER_THRESHOLD_MM = 35.0
MOTION_THRESHOLD_MM_S = 0.10
DETECTION_WINDOW_SAMPLES = 101
AQUADB_CLOCK_NS = 120.0
AQUADB_LIMIT_MM_S = 1e3 / AQUADB_CLOCK_NS
WELCH_NPERSEG = 2048
STEADY_TOLERANCE_FRACTION = 0.02
MIN_STEADY_SAMPLES = 4096


def moving_mean(values: np.ndarray, samples: int) -> np.ndarray:
    """Return a centred mean used only to locate motion boundaries."""
    kernel = np.ones(samples, dtype=float) / samples
    return np.convolve(values.astype(float), kernel, mode="same")


def find_experiment_transfers(pos_count: np.ndarray) -> list[tuple[int, int]]:
    """Recover the twelve v2/v5/v8 transfers from the measurement trace.

    The recording also contains setup/homing motion before the experiment and
    one v2 park transfer afterwards. The experiment transfers are the final 13
    long isolated motions excluding that last park.
    """
    # Preserve direction in the detection signal. The v2 pre-positioning move
    # reverses directly into the first experiment transfer; an absolute-value
    # mean would merge those two legs and hide their shared turning point.
    velocity_mm_s = moving_mean(
        np.diff(pos_count).astype(float), DETECTION_WINDOW_SAMPLES) / 1e3
    moving = np.abs(velocity_mm_s) > MOTION_THRESHOLD_MM_S
    long_runs: list[tuple[int, int]] = []
    i = 0
    while i < len(moving):
        if not moving[i]:
            i += 1
            continue
        j = i + 1
        direction = np.sign(velocity_mm_s[i])
        while (j < len(moving) and moving[j]
               and np.sign(velocity_mm_s[j]) == direction):
            j += 1
        travel_mm = abs(pos_count[j] - pos_count[i]) / 1e6
        if travel_mm > TRANSFER_THRESHOLD_MM:
            long_runs.append((i, j))
        i = j

    if len(long_runs) < 13:
        raise RuntimeError(
            "expected at least 13 long motions (12 experiment + park), "
            f"found {len(long_runs)}")
    transfers = long_runs[-13:-1]
    if len(transfers) != 12:
        raise RuntimeError(f"expected 12 experiment transfers, found {len(transfers)}")
    return transfers


def validate_controller_log(path: Path) -> dict[str, int]:
    """Confirm the profile/cycle structure recorded by the controller."""
    text = path.read_text(encoding="utf-8", errors="replace")
    counts = {}
    for name, _, _ in PROFILE_INFO:
        counts[name] = sum(
            f"# MOVE,{name},{cycle},{leg}," in text
            for cycle in (1, 2)
            for leg in ("low_to_pre", "pre_to_xref", "xref_to_low"))
    if counts != {"v2": 6, "v5": 6, "v8": 6}:
        raise RuntimeError(
            f"controller log does not contain the expected profile design: {counts}")
    return counts


def profile_groups(transfers: list[tuple[int, int]]) -> list[dict]:
    groups = []
    for i, (name, velocity, colour) in enumerate(PROFILE_INFO):
        runs = transfers[4 * i:4 * i + 4]
        groups.append({
            "name": name,
            "velocity": velocity,
            "colour": colour,
            "runs": runs,
            "start": runs[0][0],
            "end": runs[-1][1],
        })
    return groups


def panel_label(ax, text: str) -> None:
    ax.text(0.012, 0.955, text, transform=ax.transAxes, va="top",
            fontsize=9.2, fontweight="bold", color=INK,
            bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 1.8})


def raw_experiment_figure(scope, groups: list[dict]):
    """Publication-style overview plus one unfiltered panel per profile."""
    t = scope.time_s
    x = scope.displacement_mm - np.min(scope.displacement_mm)
    experiment_start = groups[0]["start"]
    experiment_end = groups[-1]["end"]

    fig = plt.figure(figsize=(10.6, 7.6))

    ax = fig.add_axes([0.080, 0.590, 0.900, 0.340])
    ax.plot(t, x, color=INK2, lw=0.62, rasterized=True,
            label="raw unwrapped counter")
    ax.axvspan(t[0], t[experiment_start], color=INK3, alpha=0.08, lw=0)
    ax.axvspan(t[experiment_end], t[-1], color=INK3, alpha=0.08, lw=0)
    for group in groups:
        a, b = group["start"], group["end"]
        ax.axvspan(t[a], t[b], color=group["colour"], alpha=0.12, lw=0)
        ax.text(0.5 * (t[a] + t[b]), 1.025, group["name"],
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=9.0, fontweight="bold", color=group["colour"])
    ax.text(0.5 * (t[0] + t[experiment_start]), 1.025, "setup",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=8.0, color=INK3)
    ax.text(0.5 * (t[experiment_end] + t[-1]), 1.025, "park",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=8.0, color=INK3)
    panel_label(ax, "(a)")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(-4, 126)
    ax.set_xlabel("elapsed recording time, $t$ (s)")
    ax.set_ylabel("relative displacement, $x$ (mm)")
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.grid(axis="y")
    ax.grid(axis="x", visible=False)

    left, right, gap = 0.080, 0.980, 0.032
    width = (right - left - 2 * gap) / 3
    for i, group in enumerate(groups):
        bx = fig.add_axes([left + i * (width + gap), 0.135, width, 0.320])
        a, b = group["start"], group["end"]
        local_t = t[a:b] - t[a]
        local_x = x[a:b]
        bx.plot(local_t, local_x, color=group["colour"], lw=0.72,
                rasterized=True)
        bx.set_xlim(0, local_t[-1])
        bx.set_ylim(-4, 66)
        bx.set_xlabel("time within profile (s)")
        if i == 0:
            bx.set_ylabel("relative displacement (mm)")
        else:
            bx.tick_params(labelleft=False)
        bx.grid(axis="y")
        bx.grid(axis="x", visible=False)
        panel_label(bx, f"({chr(ord('b') + i)})")
        bx.text(0.13, 1.035, group["name"], transform=bx.transAxes,
                va="bottom", fontsize=10.0, fontweight="bold",
                color=group["colour"])
        bx.text(0.25, 1.035,
                f"commanded $v$ = {group['velocity']:g} mm s$^{{-1}}$",
                transform=bx.transAxes, va="bottom", fontsize=8.2,
                color=group["colour"], alpha=0.82)

    legend_handles = [
        Patch(facecolor=p[2], edgecolor="none", alpha=0.45,
              label=f"{p[0]}  ({p[1]:g} mm s$^{{-1}}$)")
        for p in PROFILE_INFO
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.53, 0.030), frameon=False,
               handlelength=1.7, columnspacing=2.2)
    return fig


def velocity_profiles_figure(scope, groups: list[dict]):
    """Recorded velocity for each complete profile, including every dwell."""
    dt_s = scope.sample_ms / 1000.0
    velocity_mm_s = np.diff(scope.counter_unwrapped).astype(float) / 1e6 / dt_s
    velocity_mean = moving_mean(velocity_mm_s, DETECTION_WINDOW_SAMPLES)

    fig, axes = plt.subplots(3, 1, figsize=(10.6, 9.6))
    fig.subplots_adjust(left=0.105, right=0.985, top=0.940, bottom=0.075,
                        hspace=0.46)

    for i, (ax, group) in enumerate(zip(axes, groups)):
        start = group["start"]
        # The next profile begins immediately after the fourth low-position
        # dwell. The final v8 dwell ends when the 2 mm/s park move begins.
        window_end = (groups[i + 1]["start"] if i + 1 < len(groups)
                      else min(len(velocity_mm_s), group["end"] + 10_000))
        sample = np.arange(start, window_end)
        local_t = (sample - start) * dt_s

        ax.plot(local_t, velocity_mm_s[start:window_end],
                color=group["colour"], lw=0.30, alpha=0.22,
                rasterized=True, label="raw 1 ms samples")
        ax.plot(local_t, velocity_mean[start:window_end],
                color=group["colour"], lw=1.05,
                label="101 ms moving mean")

        command = group["velocity"]
        ax.axhline(command, color=group["colour"], lw=0.75,
                   ls=(0, (5, 3)), alpha=0.65)
        ax.axhline(-command, color=group["colour"], lw=0.75,
                   ls=(0, (5, 3)), alpha=0.65)
        for sign in (-1.0, 1.0):
            ax.axhline(sign * AQUADB_LIMIT_MM_S, color=INK, lw=0.72,
                       ls=(0, (2, 2)), alpha=0.60)
        ax.axhline(0, color=INK3, lw=0.65)

        dwell_intervals = []
        for run_index, (_, run_end) in enumerate(group["runs"]):
            if run_index + 1 < len(group["runs"]):
                next_start = group["runs"][run_index + 1][0]
            else:
                next_start = window_end
            if next_start > run_end:
                dwell_intervals.append((run_end, next_start))

        for dwell_index, (dwell_start, dwell_end) in enumerate(dwell_intervals):
            t0 = (dwell_start - start) * dt_s
            t1 = (dwell_end - start) * dt_s
            duration = (dwell_end - dwell_start) * dt_s
            ax.axvspan(t0, t1, color=INK3, alpha=0.10, lw=0)
            # Dwell traces sit at zero, so a common upper annotation row keeps
            # every arrow clear of both data and command/cap reference lines.
            label_y = 0.82
            ax.annotate("", xy=(t0, label_y), xytext=(t1, label_y),
                        xycoords=("data", "axes fraction"),
                        textcoords=("data", "axes fraction"),
                        arrowprops={"arrowstyle": "<->", "color": INK3,
                                    "lw": 0.65})
            ax.text(0.5 * (t0 + t1), label_y + 0.045,
                    f"{duration:.3f} s dwell",
                    transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=7.2, color=INK2)

        limit = max(9.1, command * 1.14)
        ax.set_xlim(0, local_t[-1])
        ax.set_ylim(-limit, limit)
        ax.set_ylabel("recorded velocity (mm s$^{-1}$)")
        ax.set_axisbelow(True)
        ax.grid(True, axis="both")
        panel_label(ax, f"({chr(ord('a') + i)})")
        ax.text(0.075, 1.025, group["name"], transform=ax.transAxes,
                va="bottom", fontsize=10.0, fontweight="bold",
                color=group["colour"])
        ax.text(0.125, 1.025,
                f"commanded $|v|$ = {command:g} mm s$^{{-1}}$",
                transform=ax.transAxes, va="bottom", fontsize=8.2,
                color=group["colour"], alpha=0.82)
        if i == 0:
            ax.legend(loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.02),
                      handlelength=2.1, columnspacing=1.4)
        if i == len(groups) - 1:
            ax.text(0.995, AQUADB_LIMIT_MM_S,
                    "+8.333 mm s$^{-1}$ AquadB cap",
                    transform=ax.get_yaxis_transform(), ha="right",
                    va="bottom", fontsize=7.2, color=INK2)
            ax.text(0.995, -AQUADB_LIMIT_MM_S,
                    "−8.333 mm s$^{-1}$ AquadB cap",
                    transform=ax.get_yaxis_transform(), ha="right",
                    va="top", fontsize=7.2, color=INK2)

    axes[-1].set_xlabel("time within profile (s)")
    return fig


def profile_window_end(groups: list[dict], index: int, sample_count: int) -> int:
    if index + 1 < len(groups):
        return groups[index + 1]["start"]
    return min(sample_count, groups[index]["end"] + 10_000)


def controller_handoff_gaps_ms(path: Path) -> list[dict]:
    """Read the low_to_pre -> pre_to_xref software gaps from the RMT log."""
    pattern = re.compile(
        r"# MOVE,(v[258]),([12]),(low_to_pre|pre_to_xref),"
        r"t0=(\d+),t1=(\d+)")
    moves: dict[tuple[str, int, str], tuple[int, int]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            moves[(match.group(1), int(match.group(2)), match.group(3))] = (
                int(match.group(4)), int(match.group(5)))
    gaps = []
    for name, _, _ in PROFILE_INFO:
        for cycle in (1, 2):
            first = moves[(name, cycle, "low_to_pre")]
            second = moves[(name, cycle, "pre_to_xref")]
            gaps.append({"profile": name, "cycle": cycle,
                         "gap_ms": second[0] - first[1]})
    return gaps


def montage_metrics(scope, groups: list[dict], log_path: Path) -> dict:
    """Compute dwell, saturation, and pre-position handoff measurements."""
    dt_s = scope.sample_ms / 1000.0
    pos_mm = scope.counter_unwrapped.astype(float) / 1e6
    velocity = np.diff(scope.counter_unwrapped).astype(float) / 1e6 / dt_s

    dwell_rows = []
    for group_index, group in enumerate(groups):
        window_end = profile_window_end(groups, group_index, len(velocity))
        for run_index, (_, run_end) in enumerate(group["runs"]):
            next_start = (group["runs"][run_index + 1][0]
                          if run_index + 1 < len(group["runs"])
                          else window_end)
            # Exclude 0.5 s at each boundary so transition-detection uncertainty
            # cannot leak moving samples into the settled-position estimate.
            settled = pos_mm[run_end + 500:next_start - 500]
            target = 61.35 if run_index % 2 == 0 else 2.0
            dwell_rows.append({
                "profile": group["name"],
                "cycle": run_index // 2 + 1,
                "target": "xref" if run_index % 2 == 0 else "low",
                "target_mm": target,
                "mean_raw_mm": float(np.mean(settled)),
                "sd_nm": float(np.std(settled) * 1e6),
                "samples": int(len(settled)),
                "detected_dwell_s": float((next_start - run_end) * dt_s),
            })

    # The scope counter and commanded axis do not share an absolute zero.
    # Register them at xref, then retain the low/xref span error.
    datum = float(np.mean([
        row["mean_raw_mm"] - row["target_mm"]
        for row in dwell_rows if row["target"] == "xref"
    ]))
    for row in dwell_rows:
        row["error_um"] = float(
            1000.0 * (row["mean_raw_mm"] - datum - row["target_mm"]))

    saturation = []
    for group in groups:
        transfers = []
        for transfer_index, (start, end) in enumerate(group["runs"], start=1):
            samples = np.abs(velocity[start:end])
            count = int(np.count_nonzero(samples >= AQUADB_LIMIT_MM_S - 5e-4))
            transfers.append({
                "transfer": transfer_index,
                "saturated_samples": count,
                "motion_samples": int(len(samples)),
                "saturated_percent": float(100.0 * count / len(samples)),
            })
        saturation.append({"profile": group["name"], "transfers": transfers})

    handoff_velocity = moving_mean(velocity, 11)
    handoffs = []
    for group in groups:
        for cycle, run_index in enumerate((0, 2), start=1):
            start, end = group["runs"][run_index]
            progress = pos_mm[start:end] - pos_mm[start]
            span = float(np.max(progress))
            candidates = np.flatnonzero(
                (progress > span - 4.0) & (progress < span - 2.0)) + start
            centre = int(candidates[np.argmin(handoff_velocity[candidates])])
            threshold = 0.10 * group["velocity"]
            left = centre
            while left > start and handoff_velocity[left] < threshold:
                left -= 1
            right = centre
            while right < end - 1 and handoff_velocity[right] < threshold:
                right += 1
            handoffs.append({
                "profile": group["name"], "cycle": cycle,
                "centre_index": centre,
                "minimum_velocity_mm_s": float(handoff_velocity[centre]),
                "below_10_percent_ms": float((right - left) * dt_s * 1000.0),
            })

    return {
        "dwell_datum": "mean measured xref dwell registered to 61.35 mm",
        "dwell_rows": dwell_rows,
        "saturation_threshold_mm_s": AQUADB_LIMIT_MM_S,
        "saturation": saturation,
        "handoff_velocity_mean_ms": 11.0,
        "handoff_rows": handoffs,
        "controller_handoff_gaps": controller_handoff_gaps_ms(log_path),
    }


def analysis_montage_figure(scope, groups: list[dict], metrics: dict):
    """Dwell accuracy, AquadB saturation, and RMT command handoff montage."""
    fig = plt.figure(figsize=(11.0, 8.8))
    ax = fig.add_axes([0.080, 0.585, 0.480, 0.350])
    bx = fig.add_axes([0.650, 0.585, 0.325, 0.350])
    cx = fig.add_axes([0.080, 0.095, 0.895, 0.350])

    # (a) Position errors during the settled portion of every dwell.
    for profile_index, (name, _, colour) in enumerate(PROFILE_INFO):
        rows = [r for r in metrics["dwell_rows"] if r["profile"] == name]
        for target, marker, offset, label in (
                ("xref", "o", -0.09, "$x_{ref}=61.35$ mm"),
                ("low", "s", 0.09, "$x_{low}=2.00$ mm")):
            selected = [r for r in rows if r["target"] == target]
            xs = profile_index + offset + np.array([-0.025, 0.025])
            ys = np.array([r["error_um"] for r in selected])
            yerr = np.array([r["sd_nm"] / 1000.0 for r in selected])
            ax.errorbar(xs, ys, yerr=yerr, ls="none", marker=marker, ms=5.8,
                        mfc=colour if target == "xref" else SURFACE,
                        mec=colour, mew=1.0, ecolor=colour, elinewidth=0.7,
                        capsize=2, label=label if profile_index == 0 else None)
    ax.axhline(0, color=INK3, lw=0.7)
    ax.set_xticks(range(3), [p[0] for p in PROFILE_INFO])
    ax.set_xlim(-0.45, 2.45)
    ax.set_ylabel("settled position error (µm)")
    ax.set_title("Dwell position offset, referenced at $x_{ref}$", loc="left")
    ax.grid(True)
    ax.legend(loc="center right")
    ax.text(0.985, 0.055,
            "ESP log: low stop = 2.070 mm\n"
            "(+70 µm from the 2.000 mm target)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.6,
            color=INK2)
    panel_label(ax, "(a)")

    # (b) Saturation occupancy for every retained transfer.
    for profile_index, ((name, _, colour), item) in enumerate(
            zip(PROFILE_INFO, metrics["saturation"])):
        values = np.array([r["saturated_percent"] for r in item["transfers"]])
        bx.bar(profile_index, np.mean(values), width=0.58, color=colour,
               alpha=0.72, lw=0)
        bx.plot(profile_index + np.array([-0.12, -0.04, 0.04, 0.12]), values,
                "o", ms=4.8, mfc=colour, mec=SURFACE, mew=0.65)
        count = sum(r["saturated_samples"] for r in item["transfers"])
        total = sum(r["motion_samples"] for r in item["transfers"])
        bx.text(profile_index, max(np.mean(values), 0) + 0.55,
                f"{count:,} / {total:,}", ha="center", va="bottom",
                fontsize=7.6, color=INK2)
    bx.set_xticks(range(3), [p[0] for p in PROFILE_INFO])
    bx.set_ylim(0, 15.2)
    bx.set_ylabel("motion samples at cap (%)")
    bx.set_title("AquadB saturation at 8.333 mm s$^{-1}$", loc="left")
    bx.grid(True, axis="y")
    panel_label(bx, "(b)")

    # (c) Physical response at the RMT transaction handoff.
    dt_s = scope.sample_ms / 1000.0
    velocity = np.diff(scope.counter_unwrapped).astype(float) / 1e6 / dt_s
    handoff_velocity = moving_mean(velocity, 11)
    half_before, half_after = 500, 2000
    for name, command, colour in PROFILE_INFO:
        rows = [r for r in metrics["handoff_rows"] if r["profile"] == name]
        traces = []
        for row in rows:
            centre = row["centre_index"]
            trace = handoff_velocity[centre - half_before:centre + half_after]
            traces.append(trace / command)
            time = np.arange(-half_before, half_after) * dt_s
            cx.plot(time, trace / command, color=colour, lw=0.65, alpha=0.38)
        cx.plot(time, np.mean(traces, axis=0), color=colour, lw=1.35,
                label=name)
    cx.axvline(0, color=INK, lw=0.75, ls=(0, (3, 3)))
    cx.axhline(0, color=INK3, lw=0.65)
    cx.axhline(1, color=INK3, lw=0.65, ls=(0, (4, 3)))
    cx.axvspan(-0.015, 0.045, color=INK3, alpha=0.10, lw=0)
    cx.set_xlim(-0.5, 2.0)
    cx.set_ylim(-0.08, 1.13)
    cx.set_xlabel("time from measured velocity minimum at 58.35 mm (s)")
    cx.set_ylabel("recorded velocity / commanded velocity")
    cx.set_title("RMT handoff: decelerate to pre-position, then issue final 3 mm move",
                 loc="left")
    cx.grid(True)
    cx.legend(loc="upper right", ncol=3)
    panel_label(cx, "(c)")

    gap_values = [r["gap_ms"] for r in metrics["controller_handoff_gaps"]]
    summaries = []
    for name, _, _ in PROFILE_INFO:
        values = [r["below_10_percent_ms"] for r in metrics["handoff_rows"]
                  if r["profile"] == name]
        summaries.append(f"{name}: {np.median(values):.0f} ms")
    cx.text(0.012, 0.075,
            f"ESP log gap: {min(gap_values)} ms for all six handoffs\n"
            "measured time below 10% of command (11 ms mean): "
            + " · ".join(summaries)
            + "\nRMT completion confirms emitted pulses, not mechanical settling",
            transform=cx.transAxes, ha="left", va="bottom", fontsize=7.8,
            color=INK2,
            bbox={"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.88,
                  "pad": 3.0})
    return fig


def steady_velocity_segments(scope, groups: list[dict]) -> list[dict]:
    """Select stationary cruise intervals without ramps, dwells, or restarts."""
    dt_s = scope.sample_ms / 1000.0
    velocity = np.diff(scope.counter_unwrapped).astype(float) / 1e6 / dt_s
    detector = moving_mean(velocity, DETECTION_WINDOW_SAMPLES)
    segments = []
    for group in groups:
        command = group["velocity"]
        tolerance = max(0.05, STEADY_TOLERANCE_FRACTION * command)
        start, end = group["start"], group["end"]
        selected = np.flatnonzero(
            np.abs(np.abs(detector[start:end]) - command) <= tolerance) + start
        if not len(selected):
            continue
        edges = np.flatnonzero(np.diff(selected) > 1) + 1
        profile_segments = []
        for indices in np.split(selected, edges):
            if len(indices) < MIN_STEADY_SAMPLES:
                continue
            a, b = int(indices[0]), int(indices[-1] + 1)
            values = velocity[a:b]
            profile_segments.append({
                "profile": group["name"],
                "start_index": a,
                "end_index": b,
                "start_time_s": float(scope.time_s[a]),
                "end_time_s": float(scope.time_s[b]),
                "duration_s": float((b - a) * dt_s),
                "direction": "positive" if np.mean(values) > 0 else "negative",
                "mean_velocity_mm_s": float(np.mean(values)),
                "sd_velocity_mm_s": float(np.std(values)),
                "values": values,
            })
        if len(profile_segments) != 4:
            raise RuntimeError(
                f"expected four steady cruise segments for {group['name']}, "
                f"found {len(profile_segments)}")
        segments.extend(profile_segments)
    return segments


def pwelch_velocity_figure(scope, groups: list[dict]):
    """Welch PSD of the four stationary cruise intervals at each velocity."""
    fs_hz = 1000.0 / scope.sample_ms
    segments = steady_velocity_segments(scope, groups)
    fig, axes = plt.subplots(3, 1, figsize=(10.6, 9.2), sharex=True)
    fig.subplots_adjust(left=0.110, right=0.985, top=0.955, bottom=0.085,
                        hspace=0.38)
    public_segments = []

    for panel_index, (ax, (name, command, colour)) in enumerate(
            zip(axes, PROFILE_INFO)):
        rows = [s for s in segments if s["profile"] == name]
        spectra = []
        for segment_index, row in enumerate(rows, start=1):
            frequency, psd = welch(
                row["values"], fs=fs_hz, window="hann",
                nperseg=WELCH_NPERSEG, noverlap=WELCH_NPERSEG // 2,
                detrend="constant", scaling="density", average="mean")
            spectra.append(psd)
            direction = "+" if row["direction"] == "positive" else "−"
            ax.loglog(frequency[1:], psd[1:], color=colour, lw=0.65,
                      alpha=0.27,
                      label="individual cruise segments" if segment_index == 1
                      else None)
            public_segments.append({
                k: v for k, v in row.items() if k != "values"
            } | {"segment": segment_index, "direction_symbol": direction})

        mean_psd = np.mean(spectra, axis=0)
        ax.loglog(frequency[1:], mean_psd[1:], color=colour, lw=1.35,
                  label="equal-segment mean PSD")
        ax.set_xlim(frequency[1], frequency[-1])
        ax.set_ylabel("PSD of velocity fluctuation\n"
                      "[(mm s$^{-1}$)$^2$ Hz$^{-1}$]")
        ax.grid(True, which="both")
        panel_label(ax, f"({chr(ord('a') + panel_index)})")
        ax.text(0.075, 0.955, name, transform=ax.transAxes, va="top",
                fontsize=10.0, fontweight="bold", color=colour)
        durations = [r["duration_s"] for r in rows]
        ax.text(0.125, 0.955,
                f"$|v|={command:g}$ mm s$^{{-1}}$ · four cruise segments · "
                f"{min(durations):.2f}–{max(durations):.2f} s",
                transform=ax.transAxes, va="top", fontsize=8.0,
                color=colour, alpha=0.85)
        if panel_index == 0:
            ax.legend(loc="lower left")

    axes[-1].set_xlabel("frequency (Hz)")
    axes[-1].text(
        0.995, 0.045,
        f"Hann · {WELCH_NPERSEG} samples ({WELCH_NPERSEG / fs_hz:.3f} s) · "
        f"50% overlap · $\\Delta f={fs_hz / WELCH_NPERSEG:.3f}$ Hz · "
        "constant detrend",
        transform=axes[-1].transAxes, ha="right", va="bottom",
        fontsize=7.5, color=INK2,
        bbox={"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.88,
              "pad": 2.5})
    settings = {
        "input": "raw 1 ms recorded velocity, constant mean removed per window",
        "segmentation": (
            "contiguous samples whose 101 ms mean absolute velocity is within "
            f"{100 * STEADY_TOLERANCE_FRACTION:.1f}% of command; intervals "
            f"shorter than {MIN_STEADY_SAMPLES / fs_hz:.3f} s excluded"
        ),
        "sample_rate_hz": fs_hz,
        "window": "Hann",
        "nperseg": WELCH_NPERSEG,
        "window_duration_s": WELCH_NPERSEG / fs_hz,
        "overlap_samples": WELCH_NPERSEG // 2,
        "overlap_percent": 50.0,
        "frequency_resolution_hz": fs_hz / WELCH_NPERSEG,
        "detrend": "constant",
        "scaling": "density",
        "aggregation": "arithmetic mean of four per-segment PSD estimates",
        "segments": public_segments,
    }
    return fig, settings


def write_manifest(output: Path, source: Path, log: Path, scope,
                   groups: list[dict]) -> Path:
    manifest = {
        "stage": "raw experiment rendering",
        "source": str(source.resolve()),
        "controller_log": str(log.resolve()),
        "sample_count": int(len(scope.time_s)),
        "sample_interval_ms": float(scope.sample_ms),
        "duration_s": float(scope.time_s[-1]),
        "start_time": scope.start_local_text,
        "end_time": scope.end_local_text,
        "channel": scope.channel,
        "counter_rollovers": int(scope.rollover_count),
        "processing": [
            "exact UINT32 counter unwrapping",
            "nominal conversion at 1 count = 1 nm",
            "constant offset for relative displacement",
            "101 ms moving mean used only for transfer-boundary detection",
        ],
        "profile_blocks": [
            {
                "name": g["name"],
                "command_velocity_mm_s": g["velocity"],
                "transfer_count": len(g["runs"]),
                "start_time_s": float(scope.time_s[g["start"]]),
                "end_time_s": float(scope.time_s[g["end"]]),
            }
            for g in groups
        ],
    }
    path = output / "raw_experiment_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--show", action="store_true",
        help="open the raw figure in an interactive Matplotlib window")
    args = parser.parse_args()

    scope = read_scope(args.source)
    if scope.sample_ms != 1.0:
        raise RuntimeError(f"expected 1 ms sampling, found {scope.sample_ms:g} ms")
    validate_controller_log(args.log)
    transfers = find_experiment_transfers(scope.counter_unwrapped)
    groups = profile_groups(transfers)
    montage_data = montage_metrics(scope, groups, args.log)

    args.output.mkdir(parents=True, exist_ok=True)
    style()
    fig = raw_experiment_figure(scope, groups)
    fig.savefig(args.output / "01_experiment.png")
    velocity_fig = velocity_profiles_figure(scope, groups)
    velocity_fig.savefig(args.output / "02_velocity_profiles.png")
    montage_fig = analysis_montage_figure(scope, groups, montage_data)
    montage_fig.savefig(args.output / "03_dwell_saturation_handoff_montage.png")
    montage_metrics_path = args.output / "03_montage_metrics.json"
    montage_metrics_path.write_text(
        json.dumps(montage_data, indent=2) + "\n", encoding="utf-8")
    pwelch_fig, pwelch_settings = pwelch_velocity_figure(scope, groups)
    pwelch_fig.savefig(args.output / "04_velocity_segment_pwelch.png")
    pwelch_metrics_path = args.output / "04_velocity_segment_pwelch.json"
    pwelch_metrics_path.write_text(
        json.dumps(pwelch_settings, indent=2) + "\n", encoding="utf-8")
    manifest = write_manifest(args.output, args.source, args.log, scope, groups)

    print(f"source: {args.source}")
    print(f"samples: {len(scope.time_s):,} at {scope.sample_ms:g} ms; "
          f"duration: {scope.time_s[-1]:.3f} s")
    for group in groups:
        print(f"{group['name']}: {len(group['runs'])} transfers, "
              f"{scope.time_s[group['start']]:.3f}--"
              f"{scope.time_s[group['end']]:.3f} s")
    print(f"wrote: {(args.output / '01_experiment.png').resolve()}")
    print(f"wrote: {(args.output / '02_velocity_profiles.png').resolve()}")
    print("wrote: "
          f"{(args.output / '03_dwell_saturation_handoff_montage.png').resolve()}")
    print(f"metrics: {montage_metrics_path.resolve()}")
    print(f"wrote: {(args.output / '04_velocity_segment_pwelch.png').resolve()}")
    print(f"metrics: {pwelch_metrics_path.resolve()}")
    print(f"manifest: {manifest.resolve()}")
    if args.show:
        print("interactive window open; close it to end the process")
        plt.show()
    else:
        plt.close(fig)
        plt.close(velocity_fig)
        plt.close(montage_fig)
        plt.close(pwelch_fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

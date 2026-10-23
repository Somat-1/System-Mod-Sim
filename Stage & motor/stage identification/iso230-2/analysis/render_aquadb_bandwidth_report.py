#!/usr/bin/env python3
"""Process and render the 120 ns AquadB bandwidth experiment.

The matching terminal log was not retained, so the twelve long diagnostic
transfers are recovered from the scope trace. Their order is fixed by the
generated test: four transfers at each of v2, v5, and v8. The first three in
each group are isolated low/xref transfers; the fourth is joined to the
profile-end homing sequence and is excluded from closure statistics.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from render_raw_measurement import read_scope


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE.parent.parent / "PLC_Bandwidth_test.csv"
DEFAULT_OUTPUT = HERE / "plots" / "aquadb_120ns"
MOTOR_SCRIPTS = HERE.parent.parent.parent / "motor identification" / "scripts"
if str(MOTOR_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(MOTOR_SCRIPTS))
from staircase_figures import INK, INK2, INK3, SURFACE, style  # noqa: E402


COMMAND_DISTANCE_MM = 59.350
ACCEL_MM_S2 = 40.0
AQUADB_CLOCK_NS = 120.0
CEILING_MM_S = 1e3 / AQUADB_CLOCK_NS
CEILING_NM_SAMPLE = 1e6 / AQUADB_CLOCK_NS
PASS_TOLERANCE_MM = 0.010
PROFILE_INFO = (
    ("v2", 2.0, "#5f8bc5"),
    ("v5", 5.0, "#c18453"),
    ("v8", 8.0, "#4f7862"),
)


def smooth(a: np.ndarray, n: int = 101) -> np.ndarray:
    return np.convolve(a.astype(float), np.ones(n) / n, mode="same")


def find_long_transfers(pos_count: np.ndarray) -> list[tuple[int, int]]:
    """Find the twelve diagnostic transfers from 101 ms mean speed."""
    speed = smooth(np.abs(np.diff(pos_count)), 101) / 1e3
    moving = speed > 0.03
    runs: list[tuple[int, int]] = []
    i = 0
    while i < len(moving):
        if not moving[i]:
            i += 1
            continue
        j = i + 1
        while j < len(moving) and moving[j]:
            j += 1
        if abs(pos_count[j] - pos_count[i]) / 1e6 > 35.0:
            runs.append((i, j))
        i = j
    if len(runs) != 12:
        raise RuntimeError(
            f"expected 12 long diagnostic transfers, found {len(runs)}"
        )
    return runs


def plateau(pos: np.ndarray, edge: int, before: bool) -> tuple[float, float]:
    """Mean count and scatter on a 1.5 s settled window beside a move."""
    if before:
        z = pos[max(0, edge - 2000):max(0, edge - 500)]
    else:
        z = pos[min(len(pos), edge + 500):min(len(pos), edge + 2000)]
    if len(z) < 100:
        raise RuntimeError("settled endpoint window is incomplete")
    return float(np.mean(z)), float(np.std(z))


def measure(scope, runs):
    pos = scope.counter_unwrapped
    step = np.diff(pos).astype(float)
    rows = []
    groups = []
    for profile_index, (name, velocity, colour) in enumerate(PROFILE_INFO):
        group = runs[4 * profile_index:4 * profile_index + 4]
        groups.append({"name": name, "velocity": velocity, "colour": colour,
                       "runs": group, "block": (group[0][0], group[-1][1])})
        for run_index, (a, b) in enumerate(group[:3], start=1):
            next_start = group[run_index][0]
            x0, s0 = plateau(pos, a, True)
            x1, s1 = plateau(pos, b, False)
            recorded = abs(x1 - x0) / 1e6
            duration = scope.time_s[b] - scope.time_s[a]
            z = np.abs(step[a:b])
            rows.append({
                "profile": name,
                "run": run_index,
                "command_velocity_mm_s": velocity,
                "acceleration_mm_s2": ACCEL_MM_S2,
                "command_distance_mm": COMMAND_DISTANCE_MM,
                "recorded_distance_mm": recorded,
                "error_um": 1000.0 * (recorded - COMMAND_DISTANCE_MM),
                "capture_percent": 100.0 * recorded / COMMAND_DISTANCE_MM,
                "detected_motion_s": duration,
                "detected_dwell_s": scope.time_s[next_start] - scope.time_s[b],
                "command_distance_over_time_mm_s": COMMAND_DISTANCE_MM / duration,
                "recorded_distance_over_time_mm_s": recorded / duration,
                "ceiling_occupancy_percent": 100.0 * float(
                    np.mean(z >= np.floor(CEILING_NM_SAMPLE))),
                "p99_increment_nm": float(np.quantile(z, 0.99)),
                "start_settled_sd_nm": s0,
                "end_settled_sd_nm": s1,
                "start_index": a,
                "end_index": b,
                "next_start_index": next_start,
            })
    return rows, groups


def summaries(rows):
    out = []
    for name, velocity, _ in PROFILE_INFO:
        rr = [r for r in rows if r["profile"] == name]
        distance = np.array([r["recorded_distance_mm"] for r in rr])
        item = {
            "profile": name,
            "command_velocity_mm_s": velocity,
            "acceleration_mm_s2": ACCEL_MM_S2,
            "n_transfers": len(rr),
            "mean_recorded_distance_mm": float(np.mean(distance)),
            "min_recorded_distance_mm": float(np.min(distance)),
            "max_recorded_distance_mm": float(np.max(distance)),
            "mean_error_um": float(1000.0 * (np.mean(distance) -
                                              COMMAND_DISTANCE_MM)),
            "repeatability_peak_to_peak_um": float(1000.0 * np.ptp(distance)),
            "mean_capture_percent": float(100.0 * np.mean(distance) /
                                          COMMAND_DISTANCE_MM),
            "mean_detected_motion_s": float(np.mean(
                [r["detected_motion_s"] for r in rr])),
            "mean_detected_dwell_s": float(np.mean(
                [r["detected_dwell_s"] for r in rr])),
            "mean_command_distance_over_time_mm_s": float(np.mean(
                [r["command_distance_over_time_mm_s"] for r in rr])),
            "mean_recorded_distance_over_time_mm_s": float(np.mean(
                [r["recorded_distance_over_time_mm_s"] for r in rr])),
            "mean_ceiling_occupancy_percent": float(np.mean(
                [r["ceiling_occupancy_percent"] for r in rr])),
            "max_settled_sd_nm": float(max(
                max(r["start_settled_sd_nm"], r["end_settled_sd_nm"])
                for r in rr)),
        }
        item["closure_pass"] = bool(abs(item["mean_error_um"]) <=
                                    PASS_TOLERANCE_MM * 1000.0)
        out.append(item)
    return out


def panel_label(ax, text):
    ax.text(0.012, 0.965, text, transform=ax.transAxes, va="top",
            fontsize=9.5, fontweight="bold", color=INK)


def figure_experiment(scope, groups):
    pos = scope.displacement_mm - np.min(scope.displacement_mm)
    t = scope.time_s
    fig = plt.figure(figsize=(7.6, 5.8))

    ax = fig.add_axes([0.093, 0.590, 0.890, 0.360])
    ax.plot(t, pos, color=INK, lw=0.65, rasterized=True)
    for g in groups:
        a, b = g["block"]
        ax.axvspan(t[a], t[b], color=g["colour"], alpha=0.16, lw=0)
        ax.text(0.5 * (t[a] + t[b]), 65.0, g["name"], ha="center",
                va="center", fontsize=9.5, fontweight="bold",
                color=g["colour"])
    panel_label(ax, "(a)")
    ax.set_xlim(0, t[-1])
    ax.set_ylim(-3, 70)
    ax.set_xlabel("$t$ (s)")
    ax.set_ylabel("relative $x$ (mm)")
    ax.grid(False)

    for i, g in enumerate(groups):
        bx = fig.add_axes([0.093 + 0.3005 * i, 0.098, 0.271, 0.330])
        a, b = g["block"]
        seg = pos[a:b]
        bx.plot(np.arange(len(seg)) * scope.sample_ms / 1000.0, seg,
                color=g["colour"], lw=0.72, rasterized=True)
        bx.set_xlim(0, len(seg) * scope.sample_ms / 1000.0)
        bx.set_ylim(-3, 65)
        bx.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        bx.set_xlabel("$t$ (s)")
        if i == 0:
            bx.set_ylabel("relative $x$ (mm)")
        else:
            bx.set_yticklabels([])
        bx.grid(False)
        panel_label(bx, f"({chr(ord('b') + i)})")
        bx.text(0.12, 0.96, g["name"], transform=bx.transAxes, va="top",
                fontsize=9.5, fontweight="bold", color=g["colour"])
        bx.text(0.12, 0.84, f"$v$ = {g['velocity']:g} mm s$^{{-1}}$",
                transform=bx.transAxes, va="top", fontsize=7.8, color=INK2)
    return fig


def figure_velocity(scope, rows):
    step = np.abs(np.diff(scope.counter_unwrapped)).astype(float)
    fig = plt.figure(figsize=(7.6, 3.45))
    left, right, bottom, top, gap = 0.105, 0.985, 0.175, 0.945, 0.027
    width = (right - left - 2 * gap) / 3
    for i, (name, velocity, colour) in enumerate(PROFILE_INFO):
        ax = fig.add_axes([left + i * (width + gap), bottom, width, top - bottom])
        r = next(z for z in rows if z["profile"] == name and z["run"] == 1)
        a, b = r["start_index"], r["end_index"]
        next_start = r["next_start_index"]
        end = min(len(step), next_start + 1000)
        raw_speed = step[a:end] / 1e3
        mean_speed = smooth(raw_speed, 101)
        time = np.arange(len(raw_speed)) * scope.sample_ms / 1000.0
        raw_line, = ax.plot(time, raw_speed, color=colour, lw=0.28,
                            alpha=0.28, rasterized=True, label="raw 1 ms")
        mean_line, = ax.plot(time, mean_speed, color=colour, lw=0.85,
                             rasterized=True, label="101 ms mean")
        ax.axhline(CEILING_MM_S, color=INK3, lw=0.8, ls=(0, (4, 3)))
        ceiling_hits = np.flatnonzero(
            raw_speed >= np.floor(CEILING_NM_SAMPLE) / 1e3)[::10]
        ax.plot(time[ceiling_hits],
                np.full(len(ceiling_hits), CEILING_MM_S),
                ls="none", marker="|", ms=1.6, mew=0.28,
                color=INK, alpha=0.45, zorder=5, rasterized=True)
        ax.axhline(velocity, color=INK2, lw=0.65, ls=(0, (1.5, 2.5)))
        move_end_s = (b - a) * scope.sample_ms / 1000.0
        dwell_end_s = (next_start - a) * scope.sample_ms / 1000.0
        ax.axvspan(move_end_s, dwell_end_s, color=INK3, alpha=0.10, lw=0)
        ax.annotate("", xy=(move_end_s, 0.55), xytext=(dwell_end_s, 0.55),
                    arrowprops=dict(arrowstyle="<->", color=INK, lw=0.75))
        ax.text(0.5 * (move_end_s + dwell_end_s), 0.76,
                f"{r['detected_dwell_s']:.1f} s dwell", ha="center",
                va="bottom", fontsize=7.2, color=INK)
        ax.set_xlim(0, time[-1])
        ax.set_ylim(0, 9.15)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        if i:
            ax.set_yticklabels([])
        ax.grid(False)
        panel_label(ax, f"({chr(ord('a') + i)})")
        ax.text(0.13, 0.965, name, transform=ax.transAxes, va="top",
                fontsize=9.5, fontweight="bold", color=colour)
        if i == 0:
            ax.text(0.98, CEILING_MM_S + 0.10, "8.333 mm s$^{-1}$",
                    transform=ax.get_yaxis_transform(), ha="right",
                    fontsize=7.2, color=INK3)
            ax.legend(handles=[raw_line, mean_line], loc="upper right",
                      bbox_to_anchor=(0.98, 0.86), fontsize=7.0,
                      handlelength=1.7, borderaxespad=0)
    fig.text((left + right) / 2, 0.045, "$t$ from detected move start (s)",
             ha="center")
    fig.text(0.022, (bottom + top) / 2, "recorded $|\\dot{x}|$ (mm s$^{-1}$)",
             rotation="vertical", va="center")
    return fig


def figure_verdict(scope, rows):
    pos = scope.counter_unwrapped
    step = np.abs(np.diff(pos)).astype(float)
    fig = plt.figure(figsize=(7.6, 6.15))
    ax = fig.add_axes([0.095, 0.610, 0.380, 0.335])
    bx = fig.add_axes([0.605, 0.610, 0.380, 0.335])
    cx = fig.add_axes([0.095, 0.095, 0.890, 0.405])

    bins = np.arange(0, 8501, 100)
    for name, _, colour in PROFILE_INFO:
        z = [step[r["start_index"]:r["end_index"]]
             for r in rows if r["profile"] == name]
        ax.hist(np.concatenate(z), bins=bins, histtype="step", lw=1.15,
                color=colour, label=name)
    ax.axvline(CEILING_NM_SAMPLE, color=INK, lw=0.8, ls=(0, (4, 3)))
    ax.set_yscale("log")
    ax.set_xlim(0, 8550)
    ax.set_xlabel("displacement per 1 ms sample (nm)")
    ax.set_ylabel("samples")
    ax.grid(False)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.01), ncol=3,
              borderaxespad=0, handlelength=1.7, columnspacing=1.0)
    panel_label(ax, "(a)")

    offsets = (-0.13, 0.0, 0.13)
    for i, (name, _, colour) in enumerate(PROFILE_INFO):
        for off, r in zip(offsets, [z for z in rows if z["profile"] == name]):
            bx.plot(i + off, r["recorded_distance_mm"], "o", ms=5.8,
                    mfc=colour, mec=SURFACE, mew=0.75)
    bx.axhline(COMMAND_DISTANCE_MM, color=INK, lw=0.85, ls=(0, (5, 3)),
               label="commanded span")
    bx.set_xlim(-0.45, 2.45)
    bx.set_ylim(38, 62)
    bx.set_xticks(range(3), [p[0] for p in PROFILE_INFO])
    bx.set_ylabel("recorded transfer (mm)")
    bx.grid(False)
    bx.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01),
              borderaxespad=0, handlelength=2.0)
    panel_label(bx, "(b)")

    for name, _, colour in PROFILE_INFO:
        r = next(z for z in rows if z["profile"] == name and z["run"] == 1)
        a, b = r["start_index"], r["end_index"]
        start, _ = plateau(pos, a, True)
        finish = min(len(pos), b + 4000)
        after = np.mean(pos[b + 500:min(len(pos), b + 2000)])
        direction = np.sign(after - start)
        first = max(0, a - 1000)
        y = direction * (pos[first:finish] - start) / 1e6
        tt = (np.arange(first, finish) - a) * scope.sample_ms / 1000.0
        cx.plot(tt, y, color=colour, lw=0.85, label=name, rasterized=True)
    cx.axhline(COMMAND_DISTANCE_MM, color=INK, lw=0.85, ls=(0, (5, 3)),
               label="commanded span")
    cx.set_xlim(-1, 42)
    cx.set_ylim(-2, 63)
    cx.set_xlabel("$t$ from detected move start (s)")
    cx.set_ylabel("recorded travel (mm)")
    cx.grid(False)
    cx.legend(loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=4,
              borderaxespad=0, handlelength=2.0, columnspacing=1.3)
    panel_label(cx, "(c)")
    return fig


def figure_saturation_impact(scope, rows):
    """Analogue of 04_v3_saturation_impact for v2, v5, and v8."""
    pos = scope.counter_unwrapped
    d = np.diff(pos).astype(float)
    dt = scope.sample_ms / 1000.0
    zoom_n = 160
    telem_ms = 20
    fig = plt.figure(figsize=(7.6, 8.0))

    left, right = 0.100, 0.985
    gap = 0.028
    width = (right - left - 2 * gap) / 3
    top_y, top_h = 0.700, 0.255
    zoom_y, zoom_h = 0.385, 0.225

    for i, (name, velocity, colour) in enumerate(PROFILE_INFO):
        r = next(z for z in rows if z["profile"] == name and z["run"] == 1)
        a, b = r["start_index"], r["end_index"]
        next_start = r["next_start_index"]
        direction = np.sign(np.mean(pos[b + 500:b + 1500]) -
                            np.mean(pos[a - 1500:a - 500]))
        if direction == 0:
            direction = 1.0

        # Top row: every raw sample through the full following dwell.
        ax = fig.add_axes([left + i * (width + gap), top_y, width, top_h])
        end = min(len(d), next_start + 800)
        seg = direction * d[a:end]
        tt = np.arange(len(seg)) * dt
        move_end = (b - a) * dt
        dwell_end = (next_start - a) * dt
        ax.axvspan(move_end, dwell_end, color=INK3, alpha=0.13, lw=0)
        ax.axvspan(move_end + 0.5, min(move_end + 2.0, dwell_end),
                   color="#356449", alpha=0.20, lw=0)
        zoom_start = a + 2000
        z0 = (zoom_start - a) * dt
        ax.axvspan(z0, z0 + zoom_n * dt, color="#99483e", alpha=0.24, lw=0)
        ax.plot(tt, seg, color=colour, lw=0.36, rasterized=True)
        for sign in (-1, 1):
            ax.axhline(sign * CEILING_NM_SAMPLE, color=INK, lw=0.8,
                       ls=(0, (4, 3)))
        ax.axhline(0, color=INK3, lw=0.6)
        ax.set_xlim(tt[0], tt[-1])
        ax.set_ylim(-9200, 9200)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        ax.set_xlabel("$t$ (s)", fontsize=8.5)
        if i == 0:
            ax.set_ylabel("displacement per 1 ms (nm)", fontsize=8.5)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=7.5)
        ax.grid(False)
        panel_label(ax, f"({chr(ord('a') + i)})")
        ax.text(0.17, 0.965, name, transform=ax.transAxes, va="top",
                fontsize=9.5, fontweight="bold", color=colour)
        ax.text(0.32, 0.96, f"$v$ = {velocity:g} mm s$^{{-1}}$",
                transform=ax.transAxes, va="top", fontsize=7.3, color=INK2)
        ax.text(0.5 * (move_end + dwell_end), -8200,
                f"dwell {r['detected_dwell_s']:.1f} s",
                ha="center", va="bottom", fontsize=7.2, color=INK2)
        if i == 0:
            ax.text(0.98, CEILING_NM_SAMPLE + 350,
                    "+8.333 mm s$^{-1}$", transform=ax.get_yaxis_transform(),
                    ha="right", fontsize=6.8, color=INK)

        # Middle row: sample-level view of the highlighted 160 ms strip.
        bx = fig.add_axes([left + i * (width + gap), zoom_y, width, zoom_h])
        z = direction * d[zoom_start:zoom_start + zoom_n]
        tz = np.arange(len(z))
        bx.plot(tz, z, color=colour, lw=0.9, marker="o", ms=1.6)
        for sign in (-1, 1):
            bx.axhline(sign * CEILING_NM_SAMPLE, color=INK, lw=0.8,
                       ls=(0, (4, 3)))
        bx.axhline(0, color=INK3, lw=0.6)
        first = int(np.argmin(z[:telem_ms]))
        for xv in range(first, zoom_n, telem_ms):
            bx.axvline(xv, color=INK3, lw=0.55, ls=(0, (2, 3)))
        bx.set_xlim(0, zoom_n - 1)
        bx.set_ylim(-9200, 9200)
        bx.set_xlabel("time within move (ms)", fontsize=8.5)
        if i == 0:
            bx.set_ylabel("displacement per 1 ms (nm)", fontsize=8.5)
        else:
            bx.set_yticklabels([])
        bx.tick_params(labelsize=7.5)
        bx.grid(False)
        panel_label(bx, f"({chr(ord('d') + i)})")
        if i == 0:
            bx.text(0.02, 0.04, "dotted lines 20 ms apart",
                    transform=bx.transAxes, fontsize=6.8, color=INK2)

    # Bottom left: raw ceiling occupancy for every retained transfer.
    gx = fig.add_axes([0.100, 0.075, 0.395, 0.220])
    means = []
    for i, (name, _, colour) in enumerate(PROFILE_INFO):
        rr = [r for r in rows if r["profile"] == name]
        values = np.array([r["ceiling_occupancy_percent"] for r in rr])
        means.append(values.mean())
        gx.bar(i, values.mean(), width=0.58, color=colour, alpha=0.72, lw=0)
        gx.plot(i + np.array([-0.10, 0.0, 0.10]), values, "o", ms=4.6,
                mfc=colour, mec=SURFACE, mew=0.65)
    gx.set_xticks(range(3), [p[0] for p in PROFILE_INFO])
    gx.set_ylim(0, 50)
    gx.set_ylabel("samples at ceiling (%)", fontsize=8.5)
    gx.tick_params(labelsize=7.5)
    gx.grid(False)
    panel_label(gx, "(g)")

    # Bottom right: the impact on recovered distance.
    hx = fig.add_axes([0.610, 0.075, 0.375, 0.220])
    for i, (name, _, colour) in enumerate(PROFILE_INFO):
        rr = [r for r in rows if r["profile"] == name]
        values = np.array([r["recorded_distance_mm"] for r in rr])
        hx.plot(i + np.array([-0.10, 0.0, 0.10]), values, "o", ms=5.4,
                mfc=colour, mec=SURFACE, mew=0.70)
    hx.axhline(COMMAND_DISTANCE_MM, color=INK, lw=0.85, ls=(0, (5, 3)))
    hx.set_xticks(range(3), [p[0] for p in PROFILE_INFO])
    hx.set_xlim(-0.45, 2.45)
    hx.set_ylim(38, 61.5)
    hx.set_ylabel("recorded transfer (mm)", fontsize=8.5)
    hx.tick_params(labelsize=7.5)
    hx.grid(False)
    panel_label(hx, "(h)")
    hx.text(0.98, COMMAND_DISTANCE_MM + 0.55, "59.350 mm command",
            transform=hx.get_yaxis_transform(), ha="right",
            fontsize=6.8, color=INK)
    return fig


def write_metrics(out: Path, source: Path, scope, rows, summary):
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "aquadb_120ns_transfer_metrics.csv"
    public = [k for k in rows[0]
              if k not in ("start_index", "end_index", "next_start_index")]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=public)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in public})

    supported = max(s["command_velocity_mm_s"] for s in summary
                    if s["closure_pass"])
    failed = [s["command_velocity_mm_s"] for s in summary
              if not s["closure_pass"]]
    result = {
        "source": str(source.resolve()),
        "sample_count": int(len(scope.time_s)),
        "sample_interval_ms": scope.sample_ms,
        "duration_s": float(scope.time_s[-1]),
        "start_time": scope.start_local_text,
        "end_time": scope.end_local_text,
        "counter_rollovers": scope.rollover_count,
        "aquadb_clock_ns": AQUADB_CLOCK_NS,
        "nominal_ceiling_mm_s": CEILING_MM_S,
        "command_distance_mm": COMMAND_DISTANCE_MM,
        "closure_tolerance_um": PASS_TOLERANCE_MM * 1000.0,
        "boundary_source": "scope-inferred; matching terminal log unavailable",
        "included_transfers": "first three isolated transfers per profile",
        "excluded_transfers": "fourth transfer per profile, joined to homing",
        "profile_summary": summary,
        "demonstrated_supported_velocity_mm_s": supported,
        "failed_test_velocities_mm_s": failed,
        "verdict": (f"closure demonstrated through {supported:g} mm/s; "
                    + (f"failed at {min(failed):g} mm/s" if failed else
                       "no tested profile failed")),
        "causal_caveat": (
            "No matching terminal log or independent stage encoder was supplied. "
            "The scope proves a recorded-distance shortfall, but cannot by itself "
            "separate EL5101 count loss from mechanical following error."
        ),
    }
    json_path = out / "aquadb_120ns_results.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return csv_path, json_path, result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    scope = read_scope(args.source)
    if scope.sample_ms != 1.0:
        raise RuntimeError(f"expected 1 ms sampling, found {scope.sample_ms:g} ms")
    runs = find_long_transfers(scope.counter_unwrapped)
    rows, groups = measure(scope, runs)
    summary = summaries(rows)
    csv_path, json_path, result = write_metrics(
        args.output, args.source, scope, rows, summary)

    style()
    figures = (
        ("01_experiment", figure_experiment(scope, groups)),
        ("02_velocity_profiles", figure_velocity(scope, rows)),
        ("03_bandwidth_verdict", figure_verdict(scope, rows)),
        ("04_saturation_impact_by_velocity",
         figure_saturation_impact(scope, rows)),
    )
    for stem, fig in figures:
        fig.savefig(args.output / f"{stem}.png")
        fig.savefig(args.output / f"{stem}.pdf")
        plt.close(fig)

    print(f"source: {args.source}")
    print(f"samples: {len(scope.time_s):,} at {scope.sample_ms:g} ms; "
          f"duration: {scope.time_s[-1]:.3f} s")
    print("profile  n  recorded mm       error um    capture %  clip %  status")
    for s in summary:
        print(f"{s['profile']:7} {s['n_transfers']:1d}  "
              f"{s['mean_recorded_distance_mm']:11.6f}  "
              f"{s['mean_error_um']:+10.1f}  "
              f"{s['mean_capture_percent']:9.3f}  "
              f"{s['mean_ceiling_occupancy_percent']:6.1f}  "
              f"{'PASS' if s['closure_pass'] else 'FAIL'}")
    print(f"verdict: {result['verdict']}")
    print("note: boundaries inferred from scope; no matching terminal log")
    print(f"wrote: {args.output.resolve()}")
    print(f"metrics: {csv_path.name}, {json_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

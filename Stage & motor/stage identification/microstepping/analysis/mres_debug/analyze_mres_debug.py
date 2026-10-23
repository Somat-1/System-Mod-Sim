#!/usr/bin/env python3
"""Analyze the 2026-10-22 MRES settling matrix.

The ESP32 board clock is aligned to the independent TwinCAT encoder export by
minimizing the robust command/encoder mismatch. Every reported displacement is
the difference between final 1 s means at consecutive stops; no datum or drift
correction is applied.
"""
from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parents[1]
DATA = CAMPAIGN / "data"
RAW = DATA / "MRESdebug.csv"
LOG = DATA / "terminal_20261022_171207.log"
PLAN = CAMPAIGN / "sequence" / "settling_plan.csv"
PLOTS = HERE / "plots"

COUNTS_PER_MM = 1_000_000.0
U256_PER_MM = 25_600.0
U32, U31 = 2**32, 2**31
SETTLED_S = 1.0
STOP_GUARD_S = 0.10
MRES_ORDER = [16, 32, 64, 256]
MRES_COLORS = {16: "#6f6f6b", 32: "#5f7f9f", 64: "#a87852", 256: "#557f6a"}
PHASE_COLORS = {"M1": "#7b708c", "M2": "#a38a62", "M3": "#66816d"}
DIST_ORDER = ["1/16 FS", "1/4 FS", "1 FS", "4 FS",
              "16 FS", "64 FS", "5 mm", "20 mm"]
DIST_UM = {"1/16 FS": 0.625, "1/4 FS": 2.5, "1 FS": 10.0,
           "4 FS": 40.0, "16 FS": 160.0, "64 FS": 640.0,
           "5 mm": 5000.0, "20 mm": 20000.0}


@dataclass
class Scope:
    t: np.ndarray
    x_mm: np.ndarray
    sample_ms: float
    start_epoch_s: float
    start_text: str
    end_text: str
    rollovers: int


@dataclass
class LogData:
    board_t: np.ndarray
    command_mm: np.ndarray
    row: np.ndarray
    host_elapsed: np.ndarray
    stop_board_t: dict[int, float]
    mres_events: list[tuple[float, int]]
    run_complete: bool
    mres_check_pass: bool


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "serif", "font.size": 8.5,
        "axes.labelsize": 9, "axes.titlesize": 9.5,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": "#666662",
        "axes.linewidth": 0.7, "grid.color": "#c7c7c2",
        "grid.linestyle": "--", "grid.linewidth": 0.5,
        "grid.alpha": 0.45, "savefig.dpi": 300, "figure.dpi": 130,
    })


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(0.012, 0.985, f"({label})", transform=ax.transAxes,
            ha="left", va="top", weight="bold", fontsize=9.5,
            color="#333330")


def _filetime_to_unix(ticks: int) -> float:
    return (ticks - 116_444_736_000_000_000) / 10_000_000.0


def read_scope(path: Path) -> Scope:
    raw = path.read_bytes()
    match = re.search(rb"(?m)^0\t-?\d+\s*$", raw)
    if not match:
        raise ValueError(f"No TwinCAT sample rows found in {path}")
    header = raw[:match.start()].decode("latin-1", errors="replace")
    meta: dict[str, list[str]] = {}
    for line in header.splitlines():
        parts = line.split("\t")
        if parts and parts[0].strip():
            meta[parts[0].strip()] = [p.strip() for p in parts[1:]]
    body = raw[match.start():].replace(b"EOF", b"")
    values = np.fromstring(body.replace(b"\t", b" "), sep=" ", dtype=np.int64)
    table = values.reshape(-1, 2)
    index, counter = table[:, 0], table[:, 1]
    delta_raw = np.diff(counter)
    delta = (delta_raw + U31) % U32 - U31
    unwrapped = np.empty(counter.size, dtype=np.int64)
    unwrapped[0] = counter[0]
    unwrapped[1:] = counter[0] + np.cumsum(delta, dtype=np.int64)
    sample_ms = float(meta.get("SampleTime[ms]", ["1"])[0])
    start = meta["Starttime of export"]
    end = meta.get("Endtime of export", [])
    return Scope(
        t=(index - index[0]).astype(float) * sample_ms / 1000.0,
        x_mm=(unwrapped - unwrapped[0]) / COUNTS_PER_MM,
        sample_ms=sample_ms,
        start_epoch_s=_filetime_to_unix(int(start[0])),
        start_text=" ".join(start[1:]), end_text=" ".join(end[1:]),
        rollovers=int(np.count_nonzero(np.abs(delta_raw) > U31)))


def read_log(path: Path, scope_start_epoch_s: float) -> LogData:
    bt, cmd, rows, host = [], [], [], []
    stops: dict[int, float] = {}
    events: list[tuple[float, int]] = []
    last_bt = float("nan")
    complete = checked = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s+(.*)$", line)
        if not match:
            continue
        host_s, payload = float(match.group(1)), match.group(2)
        if payload.startswith("# S,"):
            fields = payload.split(",")
            if len(fields) >= 5:
                last_bt = float(fields[1]) / 1000.0
                bt.append(last_bt)
                cmd.append(float(fields[2]) / U256_PER_MM)
                rows.append(int(fields[4]))
                host.append(host_s - scope_start_epoch_s)
        elif payload.startswith("# STOP,") and math.isfinite(last_bt):
            stops[int(payload.split(",")[1])] = last_bt
        elif payload.startswith("# MRES,") and math.isfinite(last_bt):
            events.append((last_bt, int(payload.split(",")[1])))
        elif payload.startswith("# RUN_COMPLETE"):
            complete = True
        elif payload == "# MRES_CHECK,pass":
            checked = True
    return LogData(np.asarray(bt), np.asarray(cmd), np.asarray(rows, int),
                   np.asarray(host), stops, events, complete, checked)


def read_plan(path: Path) -> list[dict[str, object]]:
    out = []
    with path.open(newline="", encoding="utf-8") as fh:
        for row_number, raw in enumerate(csv.DictReader(fh)):
            row: dict[str, object] = dict(raw)
            for key in ("block", "mres", "d_u256", "direction", "index",
                        "measured", "start_u256", "target_u256"):
                row[key] = int(raw[key])
            for key in ("speed_fs", "dwell_s", "t_start", "t_end"):
                row[key] = float(raw[key])
            # STOP records use the zero-based CSV row; index is only the
            # repeat number within one forward/reverse walk.
            row["_row"] = row_number
            out.append(row)
    return out


def robust_mad(values: np.ndarray) -> float:
    med = np.median(values)
    return float(np.median(np.abs(values - med)))


def alignment_cost(scope: Scope, log: LogData, offset_s: float) -> float:
    take = np.arange(0, log.board_t.size, 10)
    ts = log.board_t[take] + offset_s
    valid = (ts >= scope.t[0]) & (ts <= scope.t[-1])
    actual = np.interp(ts[valid], scope.t, scope.x_mm)
    residual = actual - log.command_mm[take][valid]
    return robust_mad(residual)


def align_log(scope: Scope, log: LogData) -> tuple[float, float, float]:
    hint = float(np.median(log.host_elapsed - log.board_t))
    coarse = np.arange(hint - 20.0, hint + 20.0001, 0.02)
    costs = np.array([alignment_cost(scope, log, z) for z in coarse])
    z0 = float(coarse[int(np.argmin(costs))])
    fine = np.arange(z0 - 0.04, z0 + 0.0401, 0.0005)
    fine_cost = np.array([alignment_cost(scope, log, z) for z in fine])
    offset = float(fine[int(np.argmin(fine_cost))])
    take = np.arange(0, log.board_t.size, 10)
    ts = log.board_t[take] + offset
    valid = (ts >= scope.t[0]) & (ts <= scope.t[-1])
    actual = np.interp(ts[valid], scope.t, scope.x_mm)
    command = log.command_mm[take][valid]
    datum = float(np.median(actual - command))
    residual = actual - command - datum
    return offset, datum, float(np.sqrt(np.mean(residual**2)))


def settled(scope: Scope, stop_t: float, window_s: float
            ) -> tuple[float, float]:
    # The STOP record follows the last 20 Hz telemetry sample. Keep 100 ms
    # clear of the next move so that serial/sample phase uncertainty cannot
    # leak its first encoder samples into the settled mean.
    end_t = stop_t - STOP_GUARD_S
    use = (scope.t >= end_t - window_s) & (scope.t < end_t)
    expected = window_s * 1000.0 / scope.sample_ms
    if np.count_nonzero(use) < 0.9 * expected:
        raise ValueError(f"Incomplete settled window at {stop_t:.3f} s")
    return float(np.mean(scope.x_mm[use])), float(np.std(scope.x_mm[use], ddof=1))


def extract_measurements(scope: Scope, log: LogData,
                         plan: list[dict[str, object]], offset_s: float
                         ) -> list[dict[str, object]]:
    positions: dict[int, tuple[float, float, float]] = {}
    dwell_by_row = {int(p["_row"]): float(p["dwell_s"]) for p in plan}
    window_by_row = {
        row: min(SETTLED_S, max(0.2, dwell - 2.0 * STOP_GUARD_S))
        for row, dwell in dwell_by_row.items()
    }
    for row, board_t in log.stop_board_t.items():
        stop_t = board_t + offset_s
        mean, sd = settled(scope, stop_t, window_by_row[row])
        positions[row] = (stop_t, mean, sd)
    measurements = []
    for p in plan:
        if int(p["measured"]) != 1:
            continue
        row, previous = int(p["_row"]), int(p["_row"]) - 1
        if row not in positions or previous not in positions:
            raise KeyError(f"Missing stop for measured row {row} or prior row")
        stop_t, after, after_sd = positions[row]
        _, before, before_sd = positions[previous]
        direction = int(p["direction"])
        command_um = direction * int(p["d_u256"]) / U256_PER_MM * 1000.0
        recorded_um = (after - before) * 1000.0
        error_um = recorded_um - command_um
        measurements.append({
            "row": row, "block": int(p["block"]), "mres": int(p["mres"]),
            "distance": str(p["d_label"]), "direction": direction,
            "approach": ("reversal" if str(p["note"]) == "reversal"
                         else "same-direction"),
            "command_um": command_um, "recorded_um": recorded_um,
            "error_um": error_um,
            "direction_normalized_error_um": direction * error_um,
            "stop_scope_s": stop_t, "settled_before_mm": before,
            "settled_after_mm": after, "before_sd_nm": before_sd * 1e6,
            "after_sd_nm": after_sd * 1e6,
            "before_window_s": window_by_row[previous],
            "after_window_s": window_by_row[row],
        })
    return measurements


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(measurements: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[int, str, int, str], list[dict[str, object]]] = defaultdict(list)
    for row in measurements:
        key = (int(row["mres"]), str(row["distance"]), int(row["direction"]),
               str(row["approach"]))
        groups[key].append(row)
    out = []
    ordered = sorted(groups.items(), key=lambda x: (
        MRES_ORDER.index(x[0][0]), DIST_ORDER.index(x[0][1]), x[0][2], x[0][3]))
    for key, rows in ordered:
        err = np.array([float(r["error_um"]) for r in rows])
        norm = np.array([float(r["direction_normalized_error_um"]) for r in rows])
        rec = np.array([float(r["recorded_um"]) for r in rows])
        out.append({
            "mres": key[0], "distance": key[1], "direction": key[2],
            "approach": key[3], "n": len(rows),
            "command_um": float(rows[0]["command_um"]),
            "mean_recorded_um": float(np.mean(rec)),
            "mean_error_um": float(np.mean(err)),
            "mean_direction_normalized_error_um": float(np.mean(norm)),
            "sd_error_um": (float(np.std(err, ddof=1)) if len(err) > 1
                            else float("nan")),
            "rms_error_um": float(np.sqrt(np.mean(err**2))),
            "mean_abs_error_um": float(np.mean(np.abs(err))),
        })
    return out


def trajectory(scope: Scope, log: LogData, plan_row: dict[str, object],
               previous_stop_bt: float, stop_bt: float, offset_s: float,
               settled_before_mm: float
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    row = int(plan_row["_row"])
    use = ((log.row == row) & (log.board_t >= previous_stop_bt)
           & (log.board_t <= stop_bt - STOP_GUARD_S))
    board_t, command = log.board_t[use], log.command_mm[use]
    if board_t.size < 2:
        return np.array([]), np.array([]), np.array([])
    scope_t = board_t + offset_s
    actual = np.interp(scope_t, scope.t, scope.x_mm)
    direction = int(plan_row["direction"])
    actual_progress = direction * (actual - settled_before_mm) * 1000.0
    command_progress = direction * (
        command - int(plan_row["start_u256"]) / U256_PER_MM) * 1000.0
    return scope_t - scope_t[0], actual_progress, command_progress


def build_trajectories(scope: Scope, log: LogData,
                       plan: list[dict[str, object]],
                       measurements: list[dict[str, object]], offset_s: float
                       ) -> dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    by_row = {int(r["row"]): r for r in measurements}
    out = {}
    for p in plan:
        row = int(p["_row"])
        if row not in by_row:
            continue
        out[row] = trajectory(
            scope, log, p, log.stop_board_t[row - 1], log.stop_board_t[row],
            offset_s, float(by_row[row]["settled_before_mm"]))
    return out


def common_median(traces: list[tuple[np.ndarray, np.ndarray]], n: int = 500
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    traces = [(t, y) for t, y in traces if t.size > 1]
    end = min(float(t[-1]) for t, _ in traces)
    grid = np.linspace(0.0, end, n)
    values = np.vstack([np.interp(grid, t, y) for t, y in traces])
    return grid, np.median(values, axis=0), np.percentile(values, [25, 75], axis=0)


def plot_trajectory_montage(measurements, trajectories) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(7.5, 8.6), constrained_layout=True)
    for i, (ax, dist, letter) in enumerate(
            zip(axes.flat, DIST_ORDER, "abcdefgh")):
        intended = DIST_UM[dist]
        for mres in MRES_ORDER:
            rows = [r for r in measurements if r["distance"] == dist
                    and r["mres"] == mres
                    and r["approach"] == "same-direction"]
            traces = [(trajectories[int(r["row"])][0],
                       trajectories[int(r["row"])][1])
                      for r in rows if int(r["row"]) in trajectories]
            if not traces:
                continue
            t, med, band = common_median(traces)
            ax.fill_between(t, band[0], band[1], color=MRES_COLORS[mres],
                            alpha=0.10, linewidth=0)
            ax.plot(t, med, color=MRES_COLORS[mres], lw=1.15,
                    label=f"1/{mres}")
        ax.axhline(intended, color="#3e3e3b", ls=(0, (3, 2)), lw=0.8,
                   label=("intended" if dist == DIST_ORDER[0] else None))
        ax.set_title(f"{dist}  ({intended:g} µm)", pad=4)
        if i >= 6:
            ax.set_xlabel("time from move start [s]")
        if i % 2 == 0:
            ax.set_ylabel("direction-normalized travel [µm]")
        ax.grid(True)
        panel_label(ax, letter)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=5, loc="outside upper center",
               frameon=False, handlelength=2.2, columnspacing=1.2)
    fig.savefig(PLOTS / "01_trajectory_montage.png", bbox_inches="tight")
    plt.close(fig)


def plot_endpoint_overlays(scope: Scope, measurements) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(7.5, 8.6),
                             constrained_layout=True, sharex=True)
    for i, (ax, dist, letter) in enumerate(
            zip(axes.flat, DIST_ORDER, "abcdefgh")):
        intended = DIST_UM[dist]
        for mres in MRES_ORDER:
            traces = []
            rows = [r for r in measurements if r["distance"] == dist
                    and r["mres"] == mres
                    and r["approach"] == "same-direction"]
            for row in rows:
                stop = float(row["stop_scope_s"])
                endpoint = stop - STOP_GUARD_S
                use = (scope.t >= endpoint - 1.5) & (scope.t < endpoint)
                elapsed = scope.t[use] - (endpoint - 1.5)
                direction = int(row["direction"])
                error = (direction * (
                    scope.x_mm[use] - float(row["settled_before_mm"])) * 1000.0
                    - intended)
                traces.append((elapsed, error))
                ax.plot(elapsed - 1.5, error, color=MRES_COLORS[mres],
                        alpha=0.08, lw=0.5)
            if traces:
                t, med, _ = common_median(traces, n=300)
                ax.plot(t - 1.5, med, color=MRES_COLORS[mres], lw=1.35,
                        label=f"1/{mres}")
        ax.axhline(0, color="#3e3e3b", ls=(0, (3, 2)), lw=0.8)
        ax.set_title(f"{dist}  ({intended:g} µm)", pad=4)
        if i >= 6:
            ax.set_xlabel("time before dwell end [s]")
        if i % 2 == 0:
            ax.set_ylabel("position − intended endpoint [µm]")
        ax.grid(True)
        panel_label(ax, letter)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="outside upper center",
               frameon=False, handlelength=2.2, columnspacing=1.4)
    fig.savefig(PLOTS / "02_endpoint_offset_overlays.png", bbox_inches="tight")
    plt.close(fig)


def plot_tracking_error_montage(measurements, trajectories) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(7.5, 8.6),
                             constrained_layout=True)
    for i, (ax, dist, letter) in enumerate(
            zip(axes.flat, DIST_ORDER, "abcdefgh")):
        for mres in MRES_ORDER:
            rows = [r for r in measurements if r["distance"] == dist
                    and r["mres"] == mres
                    and r["approach"] == "same-direction"]
            traces = []
            for row in rows:
                key = int(row["row"])
                if key not in trajectories:
                    continue
                t, actual, command = trajectories[key]
                traces.append((t, actual - command))
            if traces:
                t, med, band = common_median(traces)
                ax.fill_between(t, band[0], band[1], color=MRES_COLORS[mres],
                                alpha=0.10, linewidth=0)
                ax.plot(t, med, color=MRES_COLORS[mres], lw=1.15,
                        label=f"1/{mres}")
        ax.axhline(0, color="#3e3e3b", ls=(0, (3, 2)), lw=0.8)
        ax.set_title(f"{dist}  ({DIST_UM[dist]:g} µm)", pad=4)
        if i >= 6:
            ax.set_xlabel("time from move start [s]")
        if i % 2 == 0:
            ax.set_ylabel("encoder − instantaneous command [µm]")
        ax.grid(True)
        panel_label(ax, letter)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="outside upper center",
               frameon=False, handlelength=2.2, columnspacing=1.4)
    fig.savefig(PLOTS / "02_command_tracking_error_overlays.png",
                bbox_inches="tight")
    plt.close(fig)


def aggregate_same_direction(measurements):
    result = {}
    for mres in MRES_ORDER:
        for dist in DIST_ORDER:
            rows = [r for r in measurements if r["mres"] == mres
                    and r["distance"] == dist
                    and r["approach"] == "same-direction"]
            err = np.array([float(r["error_um"]) for r in rows])
            result[(mres, dist)] = {
                "n": int(err.size),
                "mean_abs": float(np.mean(np.abs(err))),
                "rms": float(np.sqrt(np.mean(err**2))),
                "median_abs": float(np.median(np.abs(err))),
            }
    return result


def plot_offset_summary(measurements) -> None:
    agg = aggregate_same_direction(measurements)
    x = np.array([DIST_UM[d] for d in DIST_ORDER])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.5, 3.25),
                                 constrained_layout=True)
    for mres in MRES_ORDER:
        y = np.array([agg[(mres, d)]["mean_abs"] for d in DIST_ORDER])
        ax.plot(x, y, marker="o", ms=3.7, lw=1.2,
                color=MRES_COLORS[mres], label=f"1/{mres}")
        bx.plot(x, 100.0 * y / x, marker="o", ms=3.7, lw=1.2,
                color=MRES_COLORS[mres], label=f"1/{mres}")
    for axis in (ax, bx):
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both")
        axis.set_xlabel("commanded distance [µm]")
    ax.set_ylabel("mean absolute landing error [µm]")
    bx.set_ylabel("mean absolute landing error [% command]")
    ax.legend(frameon=False, ncol=2, handlelength=2.0)
    panel_label(ax, "a")
    panel_label(bx, "b")
    fig.savefig(PLOTS / "03_offset_summary.png", bbox_inches="tight")
    plt.close(fig)


def plot_directional_montage(measurements) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.4),
                             constrained_layout=True, sharex=True, sharey=True)
    x = np.arange(len(DIST_ORDER))
    for ax, mres, letter in zip(axes.flat, MRES_ORDER, "abcd"):
        for direction, marker, label, color in (
                (1, "o", "forward", "#5f7f9f"),
                (-1, "s", "reverse", "#a87852")):
            means, sds = [], []
            for dist in DIST_ORDER:
                rows = [r for r in measurements if r["mres"] == mres
                        and r["distance"] == dist
                        and r["direction"] == direction
                        and r["approach"] == "same-direction"]
                values = np.array([
                    float(r["direction_normalized_error_um"]) for r in rows])
                means.append(float(np.mean(values)))
                sds.append(float(np.std(values, ddof=1))
                           if values.size > 1 else 0.0)
            ax.errorbar(x, means, yerr=sds, marker=marker, ms=3.5,
                        lw=1.0, capsize=2, color=color, label=label)
        ax.axhline(0, color="#555552", lw=0.75, ls=(0, (3, 2)))
        ax.set_title(f"MRES 1/{mres}")
        ax.set_xticks(x, DIST_ORDER, rotation=35, ha="right")
        ax.set_ylabel("direction-normalized landing error [µm]")
        ax.grid(True, axis="y")
        panel_label(ax, letter)
    axes.flat[0].legend(frameon=False)
    fig.savefig(PLOTS / "04_directional_offset_montage.png",
                bbox_inches="tight")
    plt.close(fig)


def plot_alignment(scope, log, offset_s, datum_mm, rms_mm) -> None:
    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(7.5, 5.0), constrained_layout=True,
        gridspec_kw={"height_ratios": [1.55, 1]})
    take = np.arange(0, log.board_t.size, 10)
    ts = log.board_t[take] + offset_s
    ax.plot(scope.t[::100] / 60.0, scope.x_mm[::100] - datum_mm,
            color="#73736f", lw=0.7, label="encoder")
    ax.plot(ts / 60.0, log.command_mm[take], color="#5f7f9f", lw=0.65,
            alpha=0.9, label="command")
    ax.set_ylabel("position [mm]")
    ax.set_xlabel("scope time [min]")
    ax.grid(True)
    ax.legend(frameon=False, ncol=2)
    panel_label(ax, "a")
    valid = (ts >= scope.t[0]) & (ts <= scope.t[-1])
    residual = (np.interp(ts[valid], scope.t, scope.x_mm) - datum_mm
                - log.command_mm[take][valid]) * 1000.0
    bx.plot(ts[valid] / 60.0, residual, color="#557f6a", lw=0.65)
    bx.axhline(0, color="#555552", lw=0.75, ls=(0, (3, 2)))
    bx.set_ylabel("encoder − command [µm]")
    bx.set_xlabel("scope time [min]")
    bx.grid(True)
    bx.text(0.99, 0.94, f"alignment RMS = {rms_mm * 1000:.2f} µm",
            transform=bx.transAxes, ha="right", va="top", color="#555552")
    panel_label(bx, "b")
    fig.savefig(PLOTS / "05_recording_and_alignment.png", bbox_inches="tight")
    plt.close(fig)


def plot_experimental_sequence(scope: Scope, log: LogData,
                               plan: list[dict[str, object]],
                               offset_s: float, datum_mm: float) -> None:
    """Full recording and an event-level map of the campaign phases."""
    phase_spans = {}
    block_spans = {}
    for mres in MRES_ORDER:
        block_rows = [p for p in plan if int(p["mres"]) == mres]
        first_row = int(block_rows[0]["_row"])
        last_row = int(block_rows[-1]["_row"])
        block_start = float(np.min(log.board_t[log.row == first_row]) + offset_s)
        block_end = float(log.stop_board_t[last_row] + offset_s)
        block_spans[mres] = (block_start, block_end)
        for phase in ("M1", "M2", "M3"):
            rows = [p for p in block_rows if str(p["phase"]) == phase]
            first = int(rows[0]["_row"])
            last = int(rows[-1]["_row"])
            start = float(np.min(log.board_t[log.row == first]) + offset_s)
            end = float(log.stop_board_t[last] + offset_s)
            phase_spans[(mres, phase)] = (start, end)

    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(8.0, 5.8), constrained_layout=True,
        gridspec_kw={"height_ratios": [1.6, 1.0]})

    # Panel a: the independent stage encoder and the controller command.
    ax.plot(scope.t[::100] / 60.0, scope.x_mm[::100] - datum_mm,
            color="#8b8b86", lw=0.7, label="stage encoder", zorder=2)
    for mres in MRES_ORDER:
        start, end = block_spans[mres]
        use = ((log.board_t + offset_s >= start)
               & (log.board_t + offset_s <= end))
        ax.plot((log.board_t[use] + offset_s) / 60.0,
                log.command_mm[use], color=MRES_COLORS[mres], lw=0.75,
                label=f"command, 1/{mres}", zorder=3)
        ax.axvspan(start / 60.0, end / 60.0, color=MRES_COLORS[mres],
                   alpha=0.045, linewidth=0, zorder=0)
        ax.text((start + end) / 120.0, 20.8, f"MRES 1/{mres}",
                ha="center", va="top", color=MRES_COLORS[mres],
                weight="bold", fontsize=8.5)
        for phase in ("M1", "M2", "M3"):
            p0, p1 = phase_spans[(mres, phase)]
            ax.axvline(p0 / 60.0, color=PHASE_COLORS[phase],
                       lw=0.55, alpha=0.55)
    ax.set_xlim(0, scope.t[-1] / 60.0)
    ax.set_ylim(-21.5, 22.5)
    ax.set_ylabel("position about stage centre [mm]")
    ax.grid(True, axis="y")
    ax.legend(frameon=False, ncol=5, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), handlelength=2.1,
              columnspacing=1.1)
    panel_label(ax, "a")

    # Panel b: every executed move, classified by campaign phase.
    phase_handles = {}
    for p in plan:
        if str(p["kind"]) != "move":
            continue
        row = int(p["_row"])
        row_samples = log.board_t[log.row == row]
        if not row_samples.size or row not in log.stop_board_t:
            continue
        t_mid = (float(row_samples[0]) + float(log.stop_board_t[row])) / 2.0
        t_mid = (t_mid + offset_s) / 60.0
        magnitude_um = abs(int(p["target_u256"]) - int(p["start_u256"])) / 25.6
        phase = str(p["phase"])
        if phase == "M1":
            marker, size, alpha, face = "D", 19, 0.9, PHASE_COLORS[phase]
            label = "M1  set MRES + marker"
        elif phase == "M2":
            marker, size, alpha, face = "s", 15, 0.75, PHASE_COLORS[phase]
            label = "M2  conditioning"
        elif int(p["measured"]) == 1:
            marker, size, alpha, face = "o", 13, 0.72, PHASE_COLORS[phase]
            label = "M3  measured displacement"
        else:
            marker, size, alpha, face = "o", 10, 0.45, "none"
            label = "M3  positioning move"
        edge = PHASE_COLORS[phase]
        artist = bx.scatter(t_mid, magnitude_um, s=size, marker=marker,
                            facecolors=face, edgecolors=edge, linewidths=0.7,
                            alpha=alpha, zorder=3)
        phase_handles.setdefault(label, artist)

    for mres in MRES_ORDER:
        start, end = block_spans[mres]
        bx.axvspan(start / 60.0, end / 60.0, color=MRES_COLORS[mres],
                   alpha=0.045, linewidth=0, zorder=0)
        for phase in ("M1", "M2", "M3"):
            p0, p1 = phase_spans[(mres, phase)]
            if phase != "M3":
                bx.axvspan(p0 / 60.0, p1 / 60.0,
                           color=PHASE_COLORS[phase], alpha=0.07,
                           linewidth=0, zorder=1)
            bx.axvline(p0 / 60.0, color=PHASE_COLORS[phase],
                       lw=0.55, alpha=0.55)

    y_values = np.array([DIST_UM[d] for d in DIST_ORDER])
    for y in y_values:
        bx.axhline(y, color="#bdbdb8", lw=0.42, ls=(0, (2, 3)), zorder=0)
    bx.set_yscale("log")
    bx.set_yticks(y_values, ["0.625", "2.5", "10", "40", "160",
                             "640", "5,000", "20,000"])
    bx.set_xlim(0, scope.t[-1] / 60.0)
    bx.set_ylim(0.4, 30000)
    bx.set_xlabel("scope time [min]")
    bx.set_ylabel("commanded move magnitude [µm]")
    bx.grid(True, axis="x")
    bx.legend(list(phase_handles.values()), list(phase_handles),
              frameon=False, ncol=2, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), handletextpad=0.5,
              columnspacing=1.4)
    panel_label(bx, "b")
    fig.savefig(PLOTS / "00_experimental_sequence.png", bbox_inches="tight")
    plt.close(fig)


def fit_gain(measurements, mres):
    rows = [r for r in measurements if r["mres"] == mres
            and r["approach"] == "same-direction"]
    x = np.array([float(r["command_um"]) for r in rows])
    y = np.array([float(r["recorded_um"]) for r in rows])
    design = np.column_stack([x, np.ones_like(x)])
    gain, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - design @ np.array([gain, intercept])
    return {
        "gain": float(gain),
        "gain_error_ppm": float((gain - 1.0) * 1e6),
        "intercept_um": float(intercept),
        "residual_rms_um": float(np.sqrt(np.mean(residual**2))),
        "n": int(x.size),
    }


def achieved_command_speeds(log: LogData, plan: list[dict[str, object]]
                            ) -> dict[str, dict[str, float]]:
    values: dict[tuple[int, str], list[float]] = defaultdict(list)
    for p in plan:
        distance = str(p["d_label"])
        if (int(p["measured"]) != 1 or str(p["note"]) == "reversal"
                or distance not in ("5 mm", "20 mm")):
            continue
        row = int(p["_row"])
        use = log.row == row
        t = log.board_t[use]
        direction = int(p["direction"])
        progress = direction * (
            log.command_mm[use] - int(p["start_u256"]) / U256_PER_MM) * 1000.0
        magnitude = DIST_UM[distance]
        i10 = np.where(progress >= 0.1 * magnitude)[0]
        i90 = np.where(progress >= 0.9 * magnitude)[0]
        if i10.size and i90.size and t[i90[0]] > t[i10[0]]:
            speed = (0.8 * magnitude / 1000.0) / (t[i90[0]] - t[i10[0]])
            values[(int(p["mres"]), distance)].append(float(speed))
    result = {}
    for mres in MRES_ORDER:
        per_distance = {
            distance: float(np.median(values[(mres, distance)]))
            for distance in ("5 mm", "20 mm")
        }
        both = values[(mres, "5 mm")] + values[(mres, "20 mm")]
        result[f"1/{mres}"] = {
            "median_mm_s": float(np.median(both)),
            **{f"{d}_median_mm_s": v for d, v in per_distance.items()},
        }
    return result


def main() -> None:
    configure_style()
    PLOTS.mkdir(parents=True, exist_ok=True)
    scope = read_scope(RAW)
    log = read_log(LOG, scope.start_epoch_s)
    plan = read_plan(PLAN)
    if not log.run_complete or not log.mres_check_pass:
        raise RuntimeError("Selected terminal log is not a verified complete run")
    offset_s, datum_mm, alignment_rms_mm = align_log(scope, log)
    measurements = extract_measurements(scope, log, plan, offset_s)
    if len(measurements) != 240:
        raise RuntimeError(f"Expected 240 measurements, found {len(measurements)}")
    cell_summary = summarize(measurements)
    trajectories = build_trajectories(
        scope, log, plan, measurements, offset_s)

    write_csv(HERE / "mres_debug_measurements.csv", measurements)
    write_csv(HERE / "mres_debug_cell_summary.csv", cell_summary)
    plot_trajectory_montage(measurements, trajectories)
    plot_tracking_error_montage(measurements, trajectories)
    plot_endpoint_overlays(scope, measurements)
    plot_offset_summary(measurements)
    plot_directional_montage(measurements)
    plot_alignment(scope, log, offset_s, datum_mm, alignment_rms_mm)
    plot_experimental_sequence(scope, log, plan, offset_s, datum_mm)

    same = [r for r in measurements if r["approach"] == "same-direction"]
    absolute = np.array([abs(float(r["error_um"])) for r in same])
    agg = aggregate_same_direction(measurements)
    result = {
        "source": str(RAW.relative_to(CAMPAIGN)),
        "terminal_log": str(LOG.relative_to(CAMPAIGN)),
        "scope": {
            "start": scope.start_text, "end": scope.end_text,
            "duration_s": float(scope.t[-1]),
            "sample_interval_ms": scope.sample_ms,
            "samples": int(scope.t.size), "rollovers": scope.rollovers,
        },
        "run": {
            "rows": len(plan), "measured_displacements": len(measurements),
            "mres_check": "pass", "run_complete": True,
            "mres": MRES_ORDER, "speed_mm_s": 0.78,
            "dwell_s": 3.0, "settled_window_s": SETTLED_S,
            "stop_guard_s": STOP_GUARD_S,
            "current_ma_rms": 400, "chopper": "SpreadCycle",
            "stealthchop": False, "interpolation": False,
        },
        "alignment": {
            "scope_time_equals_board_time_plus_s": offset_s,
            "position_datum_mm": datum_mm,
            "trajectory_rms_mm": alignment_rms_mm,
        },
        "same_direction": {
            "n": len(same),
            "median_abs_error_um": float(np.median(absolute)),
            "p95_abs_error_um": float(np.percentile(absolute, 95)),
            "max_abs_error_um": float(np.max(absolute)),
        },
        "gain_fits": {f"1/{m}": fit_gain(measurements, m)
                       for m in MRES_ORDER},
        "achieved_logical_command_velocity": achieved_command_speeds(log, plan),
        "mean_abs_error_um": {
            f"1/{m}": {d: agg[(m, d)]["mean_abs"] for d in DIST_ORDER}
            for m in MRES_ORDER
        },
    }
    (HERE / "mres_debug_results.json").write_text(
        json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

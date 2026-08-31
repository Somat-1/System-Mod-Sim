#!/usr/bin/env python3
"""Splice and render the reliable v3 dedicated-controller identification."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw_local"
PROCESSED_DIR = ROOT / "data" / "processed_local"
SPLICE_DIR = ROOT / "data" / "splices_local"
ASSET_DIR = ROOT / "rendered_assets"
IDS_PATH = RAW_DIR / "SteppingSequenceID.csv"
LOG_PATH = RAW_DIR / "identification_controller_log.csv"
INDEX_PATH = ROOT / "data" / "splice_index.csv"
SUMMARY_PATH = ROOT / "data" / "reliable_identification_summary.json"

FILETIME_UNIX_EPOCH = 116444736000000000
COUNTS_TO_NM = 1.0
LEAD_PITCH_M = 2.0e-3
MOTOR_FULL_STEPS_PER_REV = 200
CONFIG_ORDER = ((4, "I_50pct"), (4, "I_100pct"),
                (2, "I_50pct"), (2, "I_100pct"),
                (1, "I_50pct"), (1, "I_100pct"))
D_RATES = ("0.125", "0.375", "1.25", "3.5",
           "9.5", "27.5", "70", "200")

MEASURED_COLOR = "#136f63"
COMMAND_COLOR = "#d1495b"
MEASURED_LABEL = "Measured (IDS)"
COMMAND_LABEL = "Commanded"

# Dedicated-controller SC peak current per relative current level (see
# HARDWARE_RUNNERS.md): I_50pct -> 200 mA, I_100pct -> 400 mA.
CURRENT_LABELS = {"I_50pct": "200 mA (50%)", "I_100pct": "400 mA (100%)"}
CURRENT_SHORT = {"I_50pct": "50% I", "I_100pct": "100% I"}

THOUSANDS_FORMATTER = FuncFormatter(lambda value, _pos: f"{value:,.0f}")


def format_current(current):
    return CURRENT_LABELS.get(current, current)


def run_title(run):
    return (f"Run {run['run_index']} — MRES 1/{run['mres']}, "
            f"{format_current(run['current'])}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_ids(path: Path):
    start_filetime = None
    sample_period_ms = None
    data_line = None
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            fields = line.rstrip("\r\n").split("\t")
            if fields and fields[0] == "Starttime of export":
                start_filetime = int(fields[1])
            elif fields and fields[0] == "SampleTime[ms]":
                sample_period_ms = float(fields[1])
            if re.match(r"^\d+\t\d+\s*$", line):
                data_line = line_number
                break
    if start_filetime is None or sample_period_ms is None or data_line is None:
        raise RuntimeError("IDS export metadata or numeric data were not found")

    numeric = np.loadtxt(
        path, delimiter="\t", skiprows=data_line, usecols=(0, 1),
        dtype=np.uint64, comments="EOF",
    )
    raw = numeric[:, 1].astype(np.uint32)
    delta = np.diff(raw.astype(np.int64))
    delta[delta > 2**31] -= 2**32
    delta[delta < -(2**31)] += 2**32
    position_counts = np.empty(raw.size, dtype=np.int64)
    position_counts[0] = 0
    position_counts[1:] = np.cumsum(delta, dtype=np.int64)
    time_s = np.arange(raw.size, dtype=np.float64) * sample_period_ms * 1e-3
    start_epoch_s = (start_filetime - FILETIME_UNIX_EPOCH) / 1.0e7
    return {
        "time_s": time_s,
        "raw_uint32": raw,
        "position_nm": position_counts.astype(np.float64) * COUNTS_TO_NM,
        "sample_period_ms": sample_period_ms,
        "start_epoch_s": start_epoch_s,
        "start_utc": datetime.fromtimestamp(
            start_epoch_s, timezone.utc
        ).isoformat(),
    }


def parse_controller_log(path: Path, ids_start_epoch_s: float):
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            instant = datetime.fromisoformat(row["utc"])
            row["epoch_s"] = instant.timestamp()
            row["ids_time_s"] = row["epoch_s"] - ids_start_epoch_s
            row["run_index_int"] = int(row["run_index"] or 0)
            row["mres_int"] = int(row["mres"] or 0)
            rows.append(row)
    return rows


def event(rows, name, run_index=None):
    matches = [row for row in rows if row["event"] == name and
               (run_index is None or row["run_index_int"] == run_index)]
    if not matches:
        raise RuntimeError(f"Missing controller event {name}, run={run_index}")
    return matches[0]


def sample_bounds(start_s, end_s, sample_period_s, sample_count):
    start = max(0, int(np.ceil(start_s / sample_period_s)))
    end = min(sample_count, int(np.floor(end_s / sample_period_s)) + 1)
    if end <= start:
        raise RuntimeError(f"Empty splice at {start_s:.6f}..{end_s:.6f} s")
    return start, end


def slug(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def block_median(time_s, values, width=20):
    usable = values.size - values.size % width
    if usable < width:
        return time_s, values
    return (time_s[:usable].reshape(-1, width).mean(axis=1),
            np.median(values[:usable].reshape(-1, width), axis=1))


def command_series(rows, start_s, end_s):
    moves = [row for row in rows if row["event"] == "MOVE_ACK" and
             row["ideal_position_rev"] and
             start_s <= row["ids_time_s"] <= end_s]
    times = [start_s]
    positions_um = [0.0]
    for row in moves:
        times.append(row["ids_time_s"])
        positions_um.append(float(row["ideal_position_rev"]) *
                            LEAD_PITCH_M * 1.0e6)
    times.append(end_s)
    positions_um.append(positions_um[-1])
    return np.asarray(times), np.asarray(positions_um)


def pair_blocks(rows):
    blocks = []
    opened = None
    for row in rows:
        if row["event"] == "BLOCK_START":
            if opened is not None:
                raise RuntimeError("Nested/unclosed BLOCK_START encountered")
            opened = row
        elif row["event"] == "BLOCK_END":
            if opened is None:
                raise RuntimeError("BLOCK_END without BLOCK_START")
            if row["block"] != opened["block"]:
                raise RuntimeError(
                    f"Block mismatch {opened['block']} -> {row['block']}"
                )
            blocks.append({
                "run_index": opened["run_index_int"],
                "current": opened["current"],
                "mres": opened["mres_int"],
                "block": opened["block"],
                "start_s": opened["ids_time_s"],
                "end_s": row["ids_time_s"],
            })
            opened = None
    if opened is not None:
        raise RuntimeError(f"Unclosed block {opened['block']}")
    return blocks


def estimate_polarity(ids, rows, campaign_start_s, campaign_end_s):
    moves = [row for row in rows if row["event"] == "MOVE_ACK" and
             row["ideal_position_rev"] and
             campaign_start_s <= row["ids_time_s"] <= campaign_end_s]
    move_t = np.asarray([campaign_start_s] +
                        [row["ids_time_s"] for row in moves])
    move_q = np.asarray([0.0] +
                        [float(row["ideal_position_rev"]) for row in moves])
    mask = (ids["time_s"] >= campaign_start_s) & (
        ids["time_s"] <= campaign_end_s)
    sample_indices = np.flatnonzero(mask)[::20]
    sample_t = ids["time_s"][sample_indices]
    target_index = np.searchsorted(move_t, sample_t, side="right") - 1
    command_rev = move_q[np.maximum(target_index, 0)]
    measured_nm = ids["position_nm"][sample_indices]
    active = np.abs(command_rev) > 1.0e-9
    x = command_rev[active] - np.mean(command_rev[active])
    y = measured_nm[active] - np.mean(measured_nm[active])
    slope = float(np.dot(x, y) / np.dot(x, x))
    return (1 if slope >= 0.0 else -1), slope


def build_splices(ids, rows, campaign_start_s, campaign_end_s, polarity):
    sample_period_s = ids["sample_period_ms"] * 1.0e-3
    aligned_nm = polarity * ids["position_nm"]
    campaign_start_index, campaign_end_index = sample_bounds(
        campaign_start_s, campaign_end_s, sample_period_s, aligned_nm.size
    )
    campaign_baseline = float(np.median(
        aligned_nm[campaign_start_index:campaign_start_index + 200]
    ))
    campaign_t = ids["time_s"][campaign_start_index:campaign_end_index]
    campaign_y = (
        aligned_nm[campaign_start_index:campaign_end_index] - campaign_baseline
    )
    np.savez_compressed(
        PROCESSED_DIR / "campaign_timeseries.npz",
        time_from_ids_start_s=campaign_t,
        time_from_campaign_start_s=campaign_t - campaign_start_s,
        position_aligned_nm=campaign_y,
        polarity=np.asarray(polarity),
    )

    runs = []
    for run_index in range(1, 7):
        start_row = event(rows, "RUN_CONFIG", run_index)
        end_row = event(rows, "RUN_COMPLETE", run_index)
        start_s, end_s = start_row["ids_time_s"], end_row["ids_time_s"]
        first, last = sample_bounds(
            start_s, end_s, sample_period_s, aligned_nm.size
        )
        baseline = float(np.median(aligned_nm[first:first + 200]))
        local_t = ids["time_s"][first:last] - start_s
        local_y = aligned_nm[first:last] - baseline
        cmd_t, cmd_um = command_series(rows, start_s, end_s)
        np.savez_compressed(
            PROCESSED_DIR / f"run_{run_index:02d}_mres_{start_row['mres_int']}_"
                            f"{slug(start_row['current'])}.npz",
            time_s=local_t,
            position_aligned_nm=local_y,
            command_time_s=cmd_t - start_s,
            command_position_um=cmd_um,
            run_index=np.asarray(run_index),
            mres=np.asarray(start_row["mres_int"]),
            current=np.asarray(start_row["current"]),
        )
        runs.append({
            "run_index": run_index,
            "current": start_row["current"],
            "mres": start_row["mres_int"],
            "start_s": start_s,
            "end_s": end_s,
            "start_sample": first,
            "end_sample_exclusive": last,
            "duration_s": end_s - start_s,
            "position_range_um": float(np.ptp(local_y) * 1.0e-3),
            "end_residual_um": float(np.median(local_y[-200:]) * 1.0e-3),
        })

    block_records = []
    for sequence, block in enumerate(pair_blocks(rows), start=1):
        first, last = sample_bounds(
            block["start_s"], block["end_s"], sample_period_s,
            aligned_nm.size,
        )
        baseline = float(np.median(aligned_nm[first:first + 100]))
        local_t = ids["time_s"][first:last] - block["start_s"]
        local_y = aligned_nm[first:last] - baseline
        cmd_t, cmd_um = command_series(
            rows, block["start_s"], block["end_s"]
        )
        filename = (
            f"run_{block['run_index']:02d}_{sequence:03d}_"
            f"{slug(block['block'])}.npz"
        )
        np.savez_compressed(
            SPLICE_DIR / filename,
            time_s=local_t,
            position_aligned_nm=local_y,
            command_time_s=cmd_t - block["start_s"],
            command_position_um=cmd_um,
            run_index=np.asarray(block["run_index"]),
            mres=np.asarray(block["mres"]),
            current=np.asarray(block["current"]),
            block=np.asarray(block["block"]),
        )
        record = dict(block)
        record.update({
            "sequence": sequence,
            "start_sample": first,
            "end_sample_exclusive": last,
            "duration_s": block["end_s"] - block["start_s"],
            "samples": last - first,
            "position_range_um": float(np.ptp(local_y) * 1.0e-3),
            "end_residual_um": float(np.median(local_y[-100:]) * 1.0e-3),
            "local_npz": str((SPLICE_DIR / filename).relative_to(ROOT)),
        })
        block_records.append(record)

    with INDEX_PATH.open("w", encoding="utf-8", newline="") as handle:
        fields = list(block_records[0])
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(block_records)
    return aligned_nm, runs, block_records


def downsample_median(x, y, bin_samples=20):
    count = (len(x) // bin_samples) * bin_samples
    if count < bin_samples:
        return x, y
    return (x[:count].reshape(-1, bin_samples).mean(axis=1),
            np.median(y[:count].reshape(-1, bin_samples), axis=1))


def measured_window(ids, aligned_nm, start_s, end_s, bin_samples=20):
    first, last = sample_bounds(
        start_s, end_s, ids["sample_period_ms"] * 1.0e-3,
        aligned_nm.size,
    )
    baseline = float(np.median(aligned_nm[first:first + 100]))
    x = ids["time_s"][first:last] - start_s
    y = (aligned_nm[first:last] - baseline) * 1.0e-3
    return downsample_median(x, y, bin_samples)


def style_axis(ax, xlabel="Time (s)", ylabel="Position (µm)"):
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_formatter(THOUSANDS_FORMATTER)


def add_shared_legend(fig, ax, ncol=2):
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=ncol, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.995, 1.0))


def plot_campaign(ids, aligned_nm, runs, campaign_start_s, campaign_end_s):
    path = ASSET_DIR / "campaign_overview.png"
    x, y = measured_window(ids, aligned_nm, campaign_start_s,
                           campaign_end_s, bin_samples=100)
    fig, ax = plt.subplots(figsize=(15, 5.5), constrained_layout=True)
    ax.plot(x / 60.0, y, color=MEASURED_COLOR, lw=0.9)
    colors = ("#e8f1f2", "#f5e6cc")
    for run in runs:
        left = (run["start_s"] - campaign_start_s) / 60.0
        right = (run["end_s"] - campaign_start_s) / 60.0
        ax.axvspan(left, right, color=colors[(run["run_index"] - 1) % 2],
                   alpha=0.6, zorder=-1)
        ax.text((left + right) / 2.0, 0.97,
                f"R{run['run_index']}  1/{run['mres']}, "
                f"{CURRENT_SHORT.get(run['current'], run['current'])}",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=8.5)
    ax.set_title("Reliable v3 identification campaign — IDS position",
                 fontsize=14, fontweight="semibold")
    style_axis(ax, xlabel="Time from campaign start (min)")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_run_montage(ids, aligned_nm, rows, runs):
    path = ASSET_DIR / "configuration_montage.png"
    fig, axes = plt.subplots(3, 2, figsize=(15, 11), sharey=True,
                             constrained_layout=True)
    for ax, run in zip(axes.flat, runs):
        x, y = measured_window(ids, aligned_nm, run["start_s"],
                               run["end_s"], bin_samples=30)
        cmd_t, cmd_y = command_series(rows, run["start_s"], run["end_s"])
        ax.plot(x, y, color=MEASURED_COLOR, lw=0.75, label=MEASURED_LABEL)
        ax.step(cmd_t - run["start_s"], cmd_y, where="post",
                color=COMMAND_COLOR, lw=0.75, alpha=0.85, label=COMMAND_LABEL)
        ax.set_title(run_title(run), fontsize=10.5)
        style_axis(ax)
    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")
    add_shared_legend(fig, axes[0, 0])
    fig.suptitle("Commanded vs. measured position — six MRES/current "
                 "configurations", fontsize=15, fontweight="semibold")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def find_block(block_records, run_index, pattern):
    matches = [b for b in block_records if b["run_index"] == run_index and
               re.search(pattern, b["block"])]
    if not matches:
        raise RuntimeError(f"No block matching {pattern!r}, run {run_index}")
    return matches[0]


def plot_selected_montage(ids, aligned_nm, rows, runs, block_records,
                          pattern, filename, title):
    path = ASSET_DIR / filename
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharey=True,
                             constrained_layout=True)
    for ax, run in zip(axes.flat, runs):
        block = find_block(block_records, run["run_index"], pattern)
        x, y = measured_window(ids, aligned_nm, block["start_s"],
                               block["end_s"], bin_samples=10)
        cmd_t, cmd_y = command_series(rows, block["start_s"], block["end_s"])
        ax.plot(x, y, color=MEASURED_COLOR, lw=0.9, label=MEASURED_LABEL)
        ax.step(cmd_t - block["start_s"], cmd_y, where="post",
                color=COMMAND_COLOR, lw=0.85, label=COMMAND_LABEL)
        ax.set_title(run_title(run), fontsize=10.5)
        style_axis(ax)
    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")
    add_shared_legend(fig, axes[0, 0])
    fig.suptitle(title, fontsize=15, fontweight="semibold")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_reference_repeatability(ids, aligned_nm, runs, block_records):
    path = ASSET_DIR / "reference_repeatability_montage.png"
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharey="row",
                             constrained_layout=True)
    for ax, run in zip(axes.flat, runs):
        for pattern, label, color in (
                (r"^BLOCK_0_START$", "Block start", MEASURED_COLOR),
                (r"^BLOCK_0_END$", "Block end", COMMAND_COLOR)):
            block = find_block(block_records, run["run_index"], pattern)
            x, y = measured_window(ids, aligned_nm, block["start_s"],
                                   block["end_s"], bin_samples=10)
            ax.plot(x, y, color=color, lw=0.9, label=label)
        ax.set_title(run_title(run), fontsize=10.5)
        style_axis(ax)
    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")
    add_shared_legend(fig, axes[0, 0])
    fig.suptitle("Reference block repeatability — run start vs. run end",
                 fontsize=15, fontweight="semibold")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_velocity_montages(ids, aligned_nm, rows, runs, block_records):
    output_dir = ASSET_DIR / "velocity_montages"
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for run in runs:
        path = output_dir / (f"run_{run['run_index']:02d}_mres_{run['mres']}_"
                             f"{slug(run['current'])}.png")
        fig, axes = plt.subplots(2, 4, figsize=(16, 7), constrained_layout=True)
        for ax, rate in zip(axes.flat, D_RATES):
            block = find_block(block_records, run["run_index"],
                               rf"^D_{re.escape(rate)}$")
            x, y = measured_window(ids, aligned_nm, block["start_s"],
                                   block["end_s"], bin_samples=5)
            cmd_t, cmd_y = command_series(rows, block["start_s"], block["end_s"])
            ax.plot(x, y, color=MEASURED_COLOR, lw=0.85, label=MEASURED_LABEL)
            ax.step(cmd_t - block["start_s"], cmd_y, where="post",
                    color=COMMAND_COLOR, lw=0.8, label=COMMAND_LABEL)
            ax.set_title(f"{rate} full-steps/s", fontsize=10.5)
            style_axis(ax)
        for ax in axes[:, 1:].flat:
            ax.set_ylabel("")
        add_shared_legend(fig, axes[0, 0])
        fig.suptitle(f"Velocity-plateau sequence — {run_title(run)}",
                     fontsize=14, fontweight="semibold")
        fig.savefig(path, dpi=170)
        plt.close(fig)
        paths.append(path)
    return paths


def main():
    for directory in (PROCESSED_DIR, SPLICE_DIR, ASSET_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    ids = parse_ids(IDS_PATH)
    rows = parse_controller_log(LOG_PATH, ids["start_epoch_s"])
    campaign_start = event(rows, "CAMPAIGN_START")["ids_time_s"]
    campaign_end = event(rows, "CAMPAIGN_COMPLETE")["ids_time_s"]
    errors = [row for row in rows if row["event"] == "ERROR"]
    if errors:
        raise RuntimeError(f"Controller log contains {len(errors)} ERROR events")
    polarity, fitted_nm_per_rev = estimate_polarity(
        ids, rows, campaign_start, campaign_end
    )
    aligned_nm, runs, blocks = build_splices(
        ids, rows, campaign_start, campaign_end, polarity
    )
    assets = [
        plot_campaign(ids, aligned_nm, runs, campaign_start, campaign_end),
        plot_run_montage(ids, aligned_nm, rows, runs),
        plot_selected_montage(ids, aligned_nm, rows, runs, blocks,
                              r"^C$", "creep_c_montage.png",
                              "Condition C — slow/creep motion"),
        plot_reference_repeatability(ids, aligned_nm, runs, blocks),
    ]
    assets.extend(plot_velocity_montages(ids, aligned_nm, rows, runs, blocks))
    complete_events = Counter(row["event"] for row in rows)
    summary = {
        "dataset_status": "primary reliable v3 identification result",
        "ids_source": str(IDS_PATH.relative_to(ROOT)),
        "controller_log": str(LOG_PATH.relative_to(ROOT)),
        "source_sha256": {IDS_PATH.name: sha256(IDS_PATH),
                          LOG_PATH.name: sha256(LOG_PATH)},
        "ids_start_utc": ids["start_utc"],
        "sample_period_ms": ids["sample_period_ms"],
        "sample_count": int(ids["time_s"].size),
        "campaign_start_s_from_ids": campaign_start,
        "campaign_end_s_from_ids": campaign_end,
        "campaign_duration_s": campaign_end - campaign_start,
        "campaign_duration_min": (campaign_end - campaign_start) / 60.0,
        "encoder_scale_nm_per_count": COUNTS_TO_NM,
        "lead_pitch_mm_per_rev": LEAD_PITCH_M * 1.0e3,
        "measurement_polarity_applied": polarity,
        "fitted_absolute_nm_per_commanded_rev": abs(fitted_nm_per_rev),
        "completed_runs": complete_events["RUN_COMPLETE"],
        "spliced_blocks": len(blocks),
        "controller_error_events": len(errors),
        "runs": runs,
        "rendered_assets": [str(path.relative_to(ROOT)) for path in assets],
    }
    with SUMMARY_PATH.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(f"Processed {ids['time_s'].size:,} IDS samples")
    print(f"Campaign: {summary['campaign_duration_min']:.2f} min, "
          f"{len(runs)} runs, {len(blocks)} blocks")
    print(f"Measurement polarity: {polarity:+d}")
    print(f"Index: {INDEX_PATH}")
    print(f"Summary: {SUMMARY_PATH}")
    print(f"Rendered assets: {len(assets)}")


if __name__ == "__main__":
    main()

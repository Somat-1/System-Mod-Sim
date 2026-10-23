#!/usr/bin/env python3
"""Compare table and floor SEIS recordings in the 2-14 minute interval.

The channels are acceleration-sensor outputs, but no sensitivity calibration
exists in the repository. Results therefore remain in recorded mV rather than
claiming a physical acceleration unit. The script streams the raw CSV files,
captures uniformly decimated recorded samples for the long-duration waveform
display, computes locally demeaned 1 s RMS trends, checks data integrity and
drift, and writes a numerical report.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TABLE = ROOT / "data" / "2026-10-15_seis" / "TvalSeis.csv"
FLOOR = ROOT / "data" / "2026-10-15_seis_floor" / "TvalSeisfloor.csv"
ANALYSIS = ROOT / "analysis"
PLOTS = ANALYSIS / "plots"
CACHE = ANALYSIS / "cache"
START_S = 2 * 60.0
END_S = 14 * 60.0
DISPLAY_BIN_S = 0.050
RAW_DISPLAY_STEP_S = 0.005
RMS_BIN_S = 1.0
READ_BLOCK = 32 * 1024 * 1024

TABLE_C = ("#4f79a7", "#b87555", "#4f8873")
FLOOR_C = "#696863"
INK = "#161616"
MUTED = "#66635f"
GRID = "#dfded9"


def csv_blocks(path: Path):
    """Yield complete float64 CSV row blocks without splitting lines."""
    with path.open("rb") as stream:
        header = stream.readline().decode("ascii").strip()
        if header != "s,mV,mV,mV":
            raise ValueError(f"unexpected header in {path.name}: {header!r}")
        tail = b""
        while True:
            chunk = stream.read(READ_BLOCK)
            if not chunk:
                break
            chunk = tail + chunk
            cut = chunk.rfind(b"\n") + 1
            tail, complete = chunk[cut:], chunk[:cut]
            if complete:
                yield np.loadtxt(io.BytesIO(complete), delimiter=",",
                                 dtype=np.float64)
        if tail.strip():
            yield np.loadtxt(io.BytesIO(tail), delimiter=",", dtype=np.float64)


def reduce_bins(carry: np.ndarray, values: np.ndarray, width: int,
                mode: str):
    combined = np.vstack((carry, values)) if carry.size else values
    complete = len(combined) // width * width
    if not complete:
        return carry if carry.size else combined.copy(), ()
    blocks = combined[:complete].reshape(-1, width, 3)
    rest = combined[complete:].copy()
    if mode == "display":
        return rest, (np.nanmin(blocks, axis=1), np.nanmax(blocks, axis=1),
                      np.nanmean(blocks, axis=1))
    means = np.nanmean(blocks, axis=1)
    rms = np.sqrt(np.nanmean((blocks - means[:, None, :]) ** 2, axis=1))
    return rest, (means, rms)


def repeated_run(values: np.ndarray, previous_value: float | None,
                 previous_run: int) -> tuple[int, float, int]:
    """Return longest identical-value run and trailing state."""
    if not len(values):
        return previous_run, previous_value, previous_run
    change = np.flatnonzero(np.diff(values) != 0) + 1
    lengths = np.diff(np.r_[0, change, len(values)])
    longest = int(lengths.max())
    leading = int(change[0]) if len(change) else len(values)
    if previous_value is not None and values[0] == previous_value:
        longest = max(longest, previous_run + leading)
    trailing = int(len(values) - change[-1]) if len(change) else len(values)
    if not len(change) and previous_value is not None and values[0] == previous_value:
        trailing += previous_run
    return longest, float(values[-1]), trailing


def analyse_recording(label: str, path: Path, start_s: float, end_s: float):
    display_min, display_max, display_mean = [], [], []
    raw_display_time, raw_display_value = [], []
    rms_mean, rms_value = [], []
    display_carry = np.empty((0, 3), dtype=np.float64)
    rms_carry = np.empty((0, 3), dtype=np.float64)
    raw_min = np.full(3, np.inf)
    raw_max = np.full(3, -np.inf)
    min_count = np.zeros(3, dtype=np.int64)
    max_count = np.zeros(3, dtype=np.int64)
    longest_run = np.ones(3, dtype=np.int64)
    previous_value = [None, None, None]
    previous_run = np.zeros(3, dtype=np.int64)
    nonfinite = np.zeros(3, dtype=np.int64)
    sample_count = 0
    fs = None
    expected_dt = None
    previous_time = None
    dt_min, dt_max, irregular_dt = np.inf, -np.inf, 0

    for number, table in enumerate(csv_blocks(path), 1):
        table = np.atleast_2d(table)
        time = table[:, 0]
        if fs is None and len(time) > 2:
            expected_dt = float(np.median(np.diff(time[:min(len(time), 20000)])))
            fs = 1.0 / expected_dt
            display_width = int(round(DISPLAY_BIN_S * fs))
            raw_display_stride = int(round(RAW_DISPLAY_STEP_S * fs))
            rms_width = int(round(RMS_BIN_S * fs))
        if time[-1] < start_s:
            continue
        if time[0] >= end_s:
            break
        keep = (time >= start_s) & (time < end_s)
        if not np.any(keep):
            continue
        selected_t = time[keep]
        values = table[keep, 1:4]

        joined_t = (np.r_[previous_time, selected_t]
                    if previous_time is not None else selected_t)
        if len(joined_t) > 1:
            delta = np.diff(joined_t)
            dt_min = min(dt_min, float(delta.min()))
            dt_max = max(dt_max, float(delta.max()))
            irregular_dt += int(np.count_nonzero(
                np.abs(delta - expected_dt) > 0.25 * expected_dt))
        previous_time = float(selected_t[-1])

        finite = np.isfinite(values)
        nonfinite += np.count_nonzero(~finite, axis=0)
        if not np.all(finite):
            values = np.where(finite, values, np.nan)

        raw_offset = (-sample_count) % raw_display_stride
        raw_indices = np.arange(raw_offset, len(values), raw_display_stride)
        raw_display_time.append(selected_t[raw_indices])
        raw_display_value.append(values[raw_indices])
        sample_count += len(values)

        for axis in range(3):
            column = values[:, axis]
            valid = column[np.isfinite(column)]
            if not len(valid):
                continue
            chunk_min, chunk_max = float(valid.min()), float(valid.max())
            count_min = int(np.count_nonzero(valid == chunk_min))
            count_max = int(np.count_nonzero(valid == chunk_max))
            if chunk_min < raw_min[axis]:
                raw_min[axis], min_count[axis] = chunk_min, count_min
            elif chunk_min == raw_min[axis]:
                min_count[axis] += count_min
            if chunk_max > raw_max[axis]:
                raw_max[axis], max_count[axis] = chunk_max, count_max
            elif chunk_max == raw_max[axis]:
                max_count[axis] += count_max
            run, previous_value[axis], previous_run[axis] = repeated_run(
                valid, previous_value[axis], int(previous_run[axis]))
            longest_run[axis] = max(longest_run[axis], run)

        display_carry, reduced = reduce_bins(
            display_carry, values, display_width, "display")
        if reduced:
            lo, hi, mean = reduced
            display_min.append(lo); display_max.append(hi); display_mean.append(mean)
        rms_carry, reduced = reduce_bins(rms_carry, values, rms_width, "rms")
        if reduced:
            mean, rms = reduced
            rms_mean.append(mean); rms_value.append(rms)
        print(f"{label}: chunk {number}, {sample_count:,} selected samples",
              flush=True)

    if not sample_count or fs is None:
        raise ValueError(f"no samples selected from {path}")
    if len(display_carry):
        display_min.append(np.nanmin(display_carry, axis=0, keepdims=True))
        display_max.append(np.nanmax(display_carry, axis=0, keepdims=True))
        display_mean.append(np.nanmean(display_carry, axis=0, keepdims=True))
    if len(rms_carry):
        mean = np.nanmean(rms_carry, axis=0, keepdims=True)
        rms = np.sqrt(np.nanmean((rms_carry - mean) ** 2, axis=0,
                                 keepdims=True))
        rms_mean.append(mean); rms_value.append(rms)

    dmean = np.vstack(display_mean)
    rmean = np.vstack(rms_mean)
    rrms = np.vstack(rms_value)
    dtime = start_s + (np.arange(len(dmean)) + 0.5) * DISPLAY_BIN_S
    rtime = start_s + (np.arange(len(rmean)) + 0.5) * RMS_BIN_S
    baseline = np.nanmedian(rmean, axis=0)
    slopes = np.array([np.polyfit((rtime - start_s) / 60.0,
                                  rmean[:, i], 1)[0] for i in range(3)])
    head = np.nanmedian(rmean[:60], axis=0)
    tail = np.nanmedian(rmean[-60:], axis=0)
    rms_slopes = np.array([np.polyfit((rtime - start_s) / 60.0,
                                      rrms[:, i], 1)[0] for i in range(3)])
    rms_head = np.nanmedian(rrms[:60], axis=0)
    rms_tail = np.nanmedian(rrms[-60:], axis=0)
    vector_rms = np.sqrt(np.sum(rrms ** 2, axis=1))
    return dict(
        label=label, source=str(path), start_s=start_s, end_s=end_s,
        fs_hz=fs, sample_count=sample_count,
        raw_time_s=np.concatenate(raw_display_time),
        raw=np.vstack(raw_display_value),
        display_time_s=dtime, minimum=np.vstack(display_min),
        maximum=np.vstack(display_max), mean=dmean,
        rms_time_s=rtime, rms_mean=rmean, rms=rrms,
        vector_rms=vector_rms, baseline=baseline,
        raw_min=raw_min, raw_max=raw_max,
        min_count=min_count, max_count=max_count,
        longest_identical_run=longest_run, nonfinite=nonfinite,
        expected_dt_s=expected_dt, dt_min_s=dt_min, dt_max_s=dt_max,
        irregular_dt=irregular_dt, drift_slope_mV_per_min=slopes,
        end_shift_mV=tail - head, rms_slope_mV_per_min=rms_slopes,
        rms_end_shift_mV=rms_tail - rms_head)


def set_plot_style():
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": "#8f8d88", "axes.linewidth": 0.65,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "grid.color": GRID, "grid.linewidth": 0.55,
        "grid.linestyle": "--",
        "legend.frameon": False, "figure.facecolor": "white",
        "savefig.facecolor": "white", "lines.solid_capstyle": "round",
    })


def render_plot(table, floor, output: Path):
    set_plot_style()
    fig, axes = plt.subplots(3, 2, figsize=(14.5, 10), sharex="col",
                             constrained_layout=False)
    names = ("X", "Y", "Z")
    for axis, name in enumerate(names):
        left, right = axes[axis]
        for data, color, label in ((floor, FLOOR_C, "Floor"),
                                   (table, TABLE_C[axis], "Optical table")):
            x = data["raw_time_s"] / 60.0
            centre = data["baseline"][axis]
            raw = data["raw"][:, axis] - centre
            raw_alpha = 0.42 if label == "Floor" else 0.78
            raw_width = 0.45 if label == "Floor" else 0.55
            left.plot(x, raw, color=color, linewidth=raw_width,
                      alpha=raw_alpha, label=label, rasterized=True)
            rms_width = 1.20 if label == "Floor" else 1.05
            right.plot(data["rms_time_s"] / 60.0, data["rms"][:, axis],
                       color=color, linewidth=rms_width, alpha=0.92,
                       label=label)
        display_values = np.concatenate((
            floor["raw"][:, axis] - floor["baseline"][axis],
            table["raw"][:, axis] - table["baseline"][axis]))
        display_limit = 1.05 * np.nanpercentile(np.abs(display_values), 99.9)
        left.set_ylim(-display_limit, display_limit)
        left.set_ylabel(f"{name} relative output (mV)")
        right.set_ylabel(f"{name} AC RMS (mV)")
        left.grid(False)
        right.grid(True, linestyle="--", alpha=0.30)
        left.spines[["top", "right"]].set_visible(False)
        right.spines[["top", "right"]].set_visible(False)
        left.text(0.012, 0.965, f"({chr(ord('a') + axis)})",
                  transform=left.transAxes, ha="left", va="top",
                  fontsize=10, fontweight="bold", color=INK)
        right.text(0.012, 0.965, f"({chr(ord('d') + axis)})",
                   transform=right.transAxes, ha="left", va="top",
                   fontsize=10, fontweight="bold", color=INK)
    axes[0, 0].set_title("Raw samples (uniform 5 ms display decimation; y-range covers 99.9%)",
                         loc="left", color=INK, pad=8)
    axes[0, 1].set_title("Locally demeaned 1 s RMS trend", loc="left",
                         color=INK, pad=8)
    axes[-1, 0].set_xlabel("Recording time (min)")
    axes[-1, 1].set_xlabel("Recording time (min)")
    for ax in axes.ravel():
        ax.set_xlim(START_S / 60.0, END_S / 60.0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, loc="upper left",
                        bbox_to_anchor=(0.071, 0.985), ncol=2,
                        columnspacing=2.0, handlelength=3.5, fontsize=11)
    for line in legend.get_lines():
        line.set_linewidth(2.3)

    stats_ax = fig.add_axes([0.45, 0.900, 0.51, 0.078])
    stats_ax.axis("off")
    floor_med = np.nanmedian(floor["rms"], axis=0)
    table_med = np.nanmedian(table["rms"], axis=0)
    ratios = table_med / floor_med
    summary = [[name, f"{floor_med[i]:.4g}", f"{table_med[i]:.4g}",
                f"{ratios[i]:.3f}"] for i, name in enumerate(names)]
    stat_table = stats_ax.table(
        cellText=summary,
        colLabels=("Axis", "Floor RMS (mV)", "Table RMS (mV)", "Table/Floor"),
        cellLoc="center", colLoc="center", loc="center", edges="closed")
    stat_table.auto_set_font_size(False)
    stat_table.set_fontsize(7.5)
    stat_table.scale(1.0, 0.88)
    for (row, col), cell in stat_table.get_celld().items():
        cell.set_edgecolor("#b9b7b1")
        cell.set_linewidth(0.45)
        cell.set_facecolor("#f1f0ed" if row == 0 else "white")
        if row == 0 or col == 0:
            cell.get_text().set_fontweight("bold")
    stats_ax.set_title("Median locally demeaned 1 s RMS", fontsize=8.5,
                       fontweight="bold", pad=1.5, color=INK)

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.065, top=0.855,
                        wspace=0.20, hspace=0.18)
    fig.savefig(output, dpi=240)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def metrics(data):
    median = np.nanmedian(data["rms"], axis=0)
    p95 = np.nanpercentile(data["rms"], 95, axis=0)
    rail_fraction = ((data["min_count"] + data["max_count"])
                     / data["sample_count"])
    digital_clip = ((rail_fraction > 1e-4)
                    | (data["longest_identical_run"] >= 5))
    drift_ratio = np.abs(data["end_shift_mV"]) / median
    return dict(median=median, p95=p95, rail_fraction=rail_fraction,
                digital_clip=digital_clip, drift_ratio=drift_ratio,
                vector_median=float(np.nanmedian(data["vector_rms"])),
                vector_p95=float(np.nanpercentile(data["vector_rms"], 95)))


def write_trend_csv(table, floor, output: Path):
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("recording", "time_s", "mean_x_mV", "mean_y_mV",
                         "mean_z_mV", "rms_x_mV", "rms_y_mV", "rms_z_mV",
                         "vector_rms_mV"))
        for data in (floor, table):
            for i, time_s in enumerate(data["rms_time_s"]):
                writer.writerow((data["label"], f"{time_s:.6f}",
                                 *[f"{v:.10g}" for v in data["rms_mean"][i]],
                                 *[f"{v:.10g}" for v in data["rms"][i]],
                                 f"{data['vector_rms'][i]:.10g}"))


def write_exclusion_manifest(output: Path):
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("source", "start_s", "end_s", "status", "reason",
                         "detected_peak_s", "detected_peak_mV"))
        writer.writerow(("TvalSeis.csv", "1560.000", "1800.000", "active",
                         "operator-confirmed room entry and table perturbation",
                         "1672.725", "109.658"))


def save_cache(table, floor, output: Path):
    arrays = {}
    for prefix, data in (("table", table), ("floor", floor)):
        for key in ("raw_time_s", "raw", "display_time_s", "minimum", "maximum", "mean",
                    "rms_time_s", "rms_mean", "rms", "vector_rms"):
            arrays[f"{prefix}_{key}"] = data[key]
    np.savez_compressed(output, **arrays)


def write_report(table, floor, output: Path):
    tm, fm = metrics(table), metrics(floor)
    axes = ("X", "Y", "Z")
    lines = [
        "# Floor versus optical-table SEIS noise, 2–14 minutes", "",
        "## Scope and method", "",
        "The two 20 kHz recordings were compared over their common 2–14 minute interval. "
        "The operator-confirmed 26–30 minute disturbance in `TvalSeis.csv` does not overlap "
        "this interval and is therefore excluded by construction.", "",
        "No accelerometer sensitivity or full-scale range is present in the project, so all "
        "results are reported as recorded channel voltage (mV), not m/s². Baseline noise is "
        "the median of locally demeaned, non-overlapping 1 s RMS values; no frequency filter "
        "was applied.", "", "## Baseline noise comparison", "",
        "| Axis | Floor median RMS (mV) | Table median RMS (mV) | Table/floor | "
        "Floor 95th pct (mV) | Table 95th pct (mV) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for i, axis in enumerate(axes):
        ratio = tm["median"][i] / fm["median"][i]
        lines.append(f"| {axis} | {fm['median'][i]:.6g} | {tm['median'][i]:.6g} | "
                     f"{ratio:.4g} | {fm['p95'][i]:.6g} | {tm['p95'][i]:.6g} |")
    vector_ratio = tm["vector_median"] / fm["vector_median"]
    lines += [
        f"| Vector √(X²+Y²+Z²) | {fm['vector_median']:.6g} | "
        f"{tm['vector_median']:.6g} | {vector_ratio:.4g} | "
        f"{fm['vector_p95']:.6g} | {tm['vector_p95']:.6g} |", "",
        "A ratio below 1 means the optical-table recording has the lower baseline noise.", "",
        "## Integrity and clipping checks", "",
        "| Recording | Axis | Raw min (mV) | Raw max (mV) | Samples at exact extrema | "
        "Longest identical run | Non-finite samples | Digital clipping flag |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for data, stat in ((floor, fm), (table, tm)):
        for i, axis in enumerate(axes):
            count = int(data["min_count"][i] + data["max_count"][i])
            lines.append(f"| {data['label']} | {axis} | {data['raw_min'][i]:.7g} | "
                         f"{data['raw_max'][i]:.7g} | {count} "
                         f"({100 * stat['rail_fraction'][i]:.5g}%) | "
                         f"{data['longest_identical_run'][i]} | {data['nonfinite'][i]} | "
                         f"{'CHECK' if stat['digital_clip'][i] else 'none detected'} |")
    all_clear = not np.any(fm["digital_clip"]) and not np.any(tm["digital_clip"])
    clipping_text = ("No digital clipping or flat-top signature was detected by the exact-extrema "
                     "and repeated-code checks. " if all_clear else
                     "At least one digital clipping heuristic was triggered and should be inspected. ")
    lines += ["", clipping_text +
              "Analog sensor saturation cannot be ruled out without the sensor full-scale range "
              "or calibration.", "", "## Sampling and drift checks", "",
              "| Recording | Samples | Fs (Hz) | Min Δt (s) | Max Δt (s) | Irregular Δt count |",
              "|---|---:|---:|---:|---:|---:|",
              f"| {floor['label']} | {floor['sample_count']} | {floor['fs_hz']:.8g} | "
              f"{floor['dt_min_s']:.8g} | {floor['dt_max_s']:.8g} | {floor['irregular_dt']} |",
              f"| {table['label']} | {table['sample_count']} | {table['fs_hz']:.8g} | "
              f"{table['dt_min_s']:.8g} | {table['dt_max_s']:.8g} | {table['irregular_dt']} |", "",
              "| Recording | Axis | Linear mean slope (mV/min) | Last–first 60 s shift (mV) | "
              "|shift| / median RMS |", "|---|---|---:|---:|---:|" ]
    for data, stat in ((floor, fm), (table, tm)):
        for i, axis in enumerate(axes):
            lines.append(f"| {data['label']} | {axis} | "
                         f"{data['drift_slope_mV_per_min'][i]:.6g} | "
                         f"{data['end_shift_mV'][i]:.6g} | {stat['drift_ratio'][i]:.4g} |")
    lines += ["", "The drift table operates on 1 s channel means. It is intentionally kept "
              "separate from AC RMS so slow offset drift is not mislabeled as broadband noise.", "",
              "### RMS non-stationarity", "",
              "| Recording | Axis | Linear RMS slope (mV/min) | Last–first 60 s median RMS (mV) |",
              "|---|---|---:|---:|"]
    for data in (floor, table):
        for i, axis in enumerate(axes):
            lines.append(f"| {data['label']} | {axis} | "
                         f"{data['rms_slope_mV_per_min'][i]:.6g} | "
                         f"{data['rms_end_shift_mV'][i]:.6g} |")
    lines += ["", "This second table describes change in vibration/noise amplitude over the "
              "segment, rather than electrical zero-offset drift.", "",
              "## Outputs", "",
              "- `plots/seis_time_domain_rms_02_14min.png`: overlaid time histories and RMS trends.",
              "- `plots/seis_time_domain_rms_02_14min.pdf`: publication-ready vector version.",
              "- `seis_time_rms_02_14min.csv`: 1 s means and locally demeaned RMS values.",
              "- `cache/seis_time_rms_02_14min.npz`: compact plotting cache.",
              "- `exclusion_windows.csv`: retained 26–30 minute disturbance definition.", ""]
    output.write_text("\n".join(lines), encoding="utf-8")


def serial_summary(table, floor):
    result = {"interval_s": [START_S, END_S], "units": "mV",
              "rms_definition": "median locally demeaned non-overlapping 1 s RMS"}
    for key, data in (("floor", floor), ("table", table)):
        stat = metrics(data)
        result[key] = {
            "source": data["source"], "sample_count": int(data["sample_count"]),
            "fs_hz": float(data["fs_hz"]),
            "median_rms_mV": stat["median"].tolist(),
            "p95_rms_mV": stat["p95"].tolist(),
            "median_vector_rms_mV": stat["vector_median"],
            "drift_slope_mV_per_min": data["drift_slope_mV_per_min"].tolist(),
            "end_shift_mV": data["end_shift_mV"].tolist(),
            "rms_slope_mV_per_min": data["rms_slope_mV_per_min"].tolist(),
            "rms_end_shift_mV": data["rms_end_shift_mV"].tolist(),
            "digital_clipping_flag": stat["digital_clip"].tolist(),
            "nonfinite_samples": data["nonfinite"].tolist(),
            "irregular_dt_count": int(data["irregular_dt"]),
        }
    result["table_to_floor_median_rms"] = (
        np.array(result["table"]["median_rms_mV"])
        / np.array(result["floor"]["median_rms_mV"])).tolist()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-s", type=float, default=START_S)
    parser.add_argument("--end-s", type=float, default=END_S)
    args = parser.parse_args()
    if (args.start_s, args.end_s) != (START_S, END_S):
        raise ValueError("this report is intentionally fixed to the requested 2–14 min interval")

    PLOTS.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    table = analyse_recording("Optical table", TABLE, args.start_s, args.end_s)
    floor = analyse_recording("Floor", FLOOR, args.start_s, args.end_s)

    plot_path = PLOTS / "seis_time_domain_rms_02_14min.png"
    render_plot(table, floor, plot_path)
    write_trend_csv(table, floor, ANALYSIS / "seis_time_rms_02_14min.csv")
    write_report(table, floor, ANALYSIS / "SEIS_TIME_RMS_02_14MIN.md")
    save_cache(table, floor, CACHE / "seis_time_rms_02_14min.npz")
    write_exclusion_manifest(ANALYSIS / "exclusion_windows.csv")
    summary = serial_summary(table, floor)
    (ANALYSIS / "seis_time_rms_02_14min_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {plot_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

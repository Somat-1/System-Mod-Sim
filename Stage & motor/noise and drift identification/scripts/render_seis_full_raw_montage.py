#!/usr/bin/env python3
"""Render the complete table and floor recordings as a 3x2 raw montage.

Every displayed point is an actual recorded sample. Samples are selected at a
uniform 5 ms cadence to match the finite figure resolution; there is no bin
mean, envelope, interpolation, or filtering. Numerical analysis remains at the
original 20 kHz rate in analyze_seis_time_rms.py.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_seis_time_rms import (
    ANALYSIS, FLOOR, FLOOR_C, INK, MUTED, PLOTS, RAW_DISPLAY_STEP_S, TABLE,
    TABLE_C, csv_blocks, set_plot_style,
)


CACHE = ANALYSIS / "cache" / "seis_full_raw_5ms.npz"
OUTPUT = PLOTS / "seis_full_raw_montage.png"
DISPLAY_PERCENTILE = 99.9


def extract_recording(path: Path) -> dict:
    times, values = [], []
    stride = None
    sample_count = 0
    fs_hz = None
    last_time = None
    for number, block in enumerate(csv_blocks(path), 1):
        block = np.atleast_2d(block)
        time = block[:, 0]
        signal = block[:, 1:4]
        if stride is None:
            dt = float(np.median(np.diff(time[:min(len(time), 20000)])))
            fs_hz = 1.0 / dt
            stride = int(round(RAW_DISPLAY_STEP_S * fs_hz))
        offset = (-sample_count) % stride
        index = np.arange(offset, len(signal), stride)
        times.append(time[index])
        values.append(signal[index])
        sample_count += len(signal)
        last_time = float(time[-1])
        print(f"{path.name}: chunk {number}, {sample_count:,} samples", flush=True)
    raw_time = np.concatenate(times)
    raw = np.vstack(values)
    first_minute = raw_time < 60.0
    baseline = np.nanmedian(raw[first_minute], axis=0)
    return dict(time_s=raw_time, raw=raw, baseline=baseline,
                sample_count=sample_count, fs_hz=fs_hz,
                duration_s=last_time + 1.0 / fs_hz)


def write_cache(table: dict, floor: dict):
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    arrays = {}
    for prefix, data in (("table", table), ("floor", floor)):
        for key, value in data.items():
            arrays[f"{prefix}_{key}"] = np.asarray(value)
    np.savez_compressed(CACHE, **arrays)


def load_cache() -> tuple[dict, dict]:
    archive = np.load(CACHE)
    def load(prefix):
        return {key: archive[f"{prefix}_{key}"] for key in
                ("time_s", "raw", "baseline", "sample_count", "fs_hz",
                 "duration_s")}
    return load("table"), load("floor")


def render(table: dict, floor: dict):
    set_plot_style()
    fig, axes = plt.subplots(3, 2, figsize=(12.2, 8.7), sharex="col",
                             constrained_layout=False)
    names = ("X", "Y", "Z")
    recordings = (
        (table, "Optical table", TABLE_C),
        (floor, "Floor", (FLOOR_C, FLOOR_C, FLOOR_C)),
    )
    letters = iter("abcdef")
    for row, name in enumerate(names):
        for col, (data, recording, colours) in enumerate(recordings):
            ax = axes[row, col]
            trace = data["raw"][:, row] - data["baseline"][row]
            colour = colours[row]
            alpha = 0.78 if col == 0 else 0.58
            ax.plot(data["time_s"] / 60.0, trace, color=colour, lw=0.38,
                    alpha=alpha, rasterized=True)
            limit = 1.05 * np.nanpercentile(np.abs(trace), DISPLAY_PERCENTILE)
            ax.set_ylim(-limit, limit)
            ax.set_xlim(0, float(data["duration_s"]) / 60.0)
            if col == 0:
                ax.axvspan(26.0, min(30.0, float(data["duration_s"]) / 60.0),
                           color="#b8aaa1", alpha=0.14, lw=0)
            ax.set_ylabel(f"{name} relative output (mV)")
            ax.grid(False)
            ax.spines[["top", "right"]].set_visible(False)
            ax.text(0.012, 0.965, f"({next(letters)})", transform=ax.transAxes,
                    ha="left", va="top", fontsize=10, fontweight="bold",
                    color=INK)
            if row == 0:
                duration_min = float(data["duration_s"]) / 60.0
                ax.set_title(f"{recording}  ({duration_min:.2f} min)",
                             fontsize=11, fontweight="bold", pad=9, color=INK)
    axes[-1, 0].set_xlabel("Recording time (min)")
    axes[-1, 1].set_xlabel("Recording time (min)")
    fig.text(0.5, 0.018,
             "Actual recorded samples at uniform 5 ms display cadence; first-minute median removed; "
             "each y-range covers 99.9% of displayed samples. Shading marks the known 26–30 min disturbance.",
             ha="center", va="bottom", fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.082, right=0.985, bottom=0.080, top=0.94,
                        wspace=0.20, hspace=0.17)
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=240)
    fig.savefig(OUTPUT.with_suffix(".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="re-read the full CSV files even if the cache exists")
    args = parser.parse_args()
    if CACHE.exists() and not args.refresh:
        table, floor = load_cache()
    else:
        table = extract_recording(TABLE)
        floor = extract_recording(FLOOR)
        write_cache(table, floor)
    render(table, floor)
    print(f"Wrote {OUTPUT}")
    for name, data in (("table", table), ("floor", floor)):
        print(f"{name}: {int(data['sample_count']):,} source samples, "
              f"{len(data['time_s']):,} displayed samples, "
              f"{float(data['duration_s']) / 60.0:.6f} min")


if __name__ == "__main__":
    main()

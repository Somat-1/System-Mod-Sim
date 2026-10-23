#!/usr/bin/env python3
"""Compute and render Welch PSDs for the common 2-14 minute interval."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal.windows import hann

from analyze_seis_time_rms import (
    ANALYSIS, FLOOR, FLOOR_C, GRID, INK, MUTED, PLOTS, START_S, END_S,
    TABLE, TABLE_C, csv_blocks, set_plot_style,
)


SEGMENT_S = 10.0
OVERLAP_FRACTION = 0.50
CACHE = ANALYSIS / "cache" / "seis_pwelch_02_14min.npz"
SETTINGS = ANALYSIS / "seis_pwelch_02_14min_settings.json"
OUTPUT = PLOTS / "seis_pwelch_02_14min.png"


class WelchAccumulator:
    """Streaming one-sided Welch PSD with linear per-segment detrending."""

    def __init__(self, fs_hz: float, nperseg: int):
        self.fs_hz = float(fs_hz)
        self.nperseg = int(nperseg)
        self.noverlap = int(round(OVERLAP_FRACTION * self.nperseg))
        self.hop = self.nperseg - self.noverlap
        self.window = hann(self.nperseg, sym=False).astype(np.float64)
        self.scale = self.fs_hz * np.sum(self.window ** 2)
        self.trend_x = (np.arange(self.nperseg, dtype=np.float64)
                        - (self.nperseg - 1) / 2.0)
        self.trend_denominator = np.dot(self.trend_x, self.trend_x)
        self.total = np.zeros((self.nperseg // 2 + 1, 3), dtype=np.float64)
        self.buffer = np.empty((0, 3), dtype=np.float32)
        self.count = 0

    def feed(self, values: np.ndarray):
        if not len(values):
            return
        values = np.asarray(values, dtype=np.float32)
        self.buffer = (np.vstack((self.buffer, values))
                       if self.buffer.size else values.copy())
        while len(self.buffer) >= self.nperseg:
            segment = self.buffer[:self.nperseg].astype(np.float64)
            segment -= np.mean(segment, axis=0, keepdims=True)
            slope = self.trend_x @ segment / self.trend_denominator
            segment -= self.trend_x[:, None] * slope[None, :]
            spectrum = np.fft.rfft(segment * self.window[:, None], axis=0)
            density = np.abs(spectrum) ** 2 / self.scale
            density[1:-1] *= 2.0
            self.total += density
            self.count += 1
            self.buffer = self.buffer[self.hop:]

    def finish(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.count:
            raise ValueError("no complete Welch segments")
        frequency = np.fft.rfftfreq(self.nperseg, d=1.0 / self.fs_hz)
        return frequency, self.total / self.count


def compute(path: Path) -> dict:
    accumulator = None
    selected_samples = 0
    for block_number, block in enumerate(csv_blocks(path), 1):
        block = np.atleast_2d(block)
        time = block[:, 0]
        if accumulator is None and len(time) > 2:
            dt = float(np.median(np.diff(time[:min(len(time), 20000)])))
            fs_hz = 1.0 / dt
            nperseg = int(round(SEGMENT_S * fs_hz))
            accumulator = WelchAccumulator(fs_hz, nperseg)
        if time[-1] < START_S:
            continue
        if time[0] >= END_S:
            break
        keep = (time >= START_S) & (time < END_S)
        values = block[keep, 1:4]
        if len(values):
            accumulator.feed(values)
            selected_samples += len(values)
            print(f"{path.name}: chunk {block_number}, "
                  f"{selected_samples:,} samples, {accumulator.count} windows",
                  flush=True)
    frequency, psd = accumulator.finish()
    return dict(frequency_hz=frequency, psd_mV2_per_Hz=psd,
                fs_hz=accumulator.fs_hz, nperseg=accumulator.nperseg,
                noverlap=accumulator.noverlap, hop=accumulator.hop,
                segment_count=accumulator.count,
                selected_samples=selected_samples)


def save_results(table: dict, floor: dict):
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        CACHE,
        frequency_hz=table["frequency_hz"],
        table_psd_mV2_per_Hz=table["psd_mV2_per_Hz"],
        floor_psd_mV2_per_Hz=floor["psd_mV2_per_Hz"],
    )
    fs_hz = float(table["fs_hz"])
    nperseg = int(table["nperseg"])
    settings = {
        "analysis_interval_s": [START_S, END_S],
        "analysis_duration_s": END_S - START_S,
        "sampling_frequency_hz": fs_hz,
        "samples_per_channel_per_recording": int(table["selected_samples"]),
        "window": "periodic Hann (DFT-even)",
        "segment_length_s": SEGMENT_S,
        "segment_length_samples": nperseg,
        "overlap_fraction": OVERLAP_FRACTION,
        "overlap_samples": int(table["noverlap"]),
        "hop_s": float(table["hop"]) / fs_hz,
        "segment_count_per_recording": int(table["segment_count"]),
        "detrend": "linear, independently per segment and channel",
        "nfft": nperseg,
        "zero_padding": False,
        "frequency_resolution_hz": fs_hz / nperseg,
        "nyquist_frequency_hz": fs_hz / 2.0,
        "hann_enbw_hz": 1.5 * fs_hz / nperseg,
        "spectrum": "one-sided power spectral density",
        "scaling": "density (mV^2/Hz)",
        "averaging": "arithmetic mean of segment periodograms",
        "post_smoothing": False,
    }
    SETTINGS.write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    return settings


def render(table: dict, floor: dict, settings: dict):
    set_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.9), sharex=True,
                             sharey=True, constrained_layout=False)
    names = ("X", "Y", "Z")
    frequency = table["frequency_hz"]
    keep = frequency > 0
    for row, name in enumerate(names):
        ax = axes[row]
        ax.loglog(frequency[keep], floor["psd_mV2_per_Hz"][keep, row],
                  color=FLOOR_C, lw=1.0, alpha=0.86, label="Floor",
                  rasterized=True)
        ax.loglog(frequency[keep], table["psd_mV2_per_Hz"][keep, row],
                  color=TABLE_C[row], lw=0.90, alpha=0.92,
                  label="Optical table", rasterized=True)
        ax.set_xlim(settings["frequency_resolution_hz"],
                    settings["nyquist_frequency_hz"])
        if row == 0:
            ax.set_ylabel("PSD (mV²/Hz)")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_title(f"Axis {row + 1} ({name})", loc="left", fontsize=12,
                     fontweight="bold", pad=8, color=INK)
        ax.grid(True, which="major", linestyle="-", color="#cfcfca",
                linewidth=0.62, alpha=0.78)
        ax.grid(True, which="minor", linestyle="-", color="#e1e0dc",
                linewidth=0.42, alpha=0.72)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(0.012, 0.955, f"({chr(ord('a') + row)})",
                transform=ax.transAxes, ha="left", va="top", fontsize=10,
                fontweight="bold", color=INK)
        if row == 0:
            legend = ax.legend(loc="upper right", fontsize=10.5,
                               handlelength=2.8)
            for line in legend.get_lines():
                line.set_linewidth(2.0)
    fig.text(
        0.5, 0.018,
        f"Welch PSD: {settings['segment_length_s']:.0f} s periodic Hann windows, "
        f"50% overlap, linear detrending, Δf = {settings['frequency_resolution_hz']:.1f} Hz, "
        f"{settings['segment_count_per_recording']} averages; 2–14 min interval.",
        ha="center", va="bottom", fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.070, right=0.990, bottom=0.165, top=0.925,
                        wspace=0.045)
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=240)
    fig.savefig(OUTPUT.with_suffix(".pdf"))
    plt.close(fig)


def main():
    table = compute(TABLE)
    floor = compute(FLOOR)
    settings = save_results(table, floor)
    render(table, floor, settings)
    print(f"Wrote {OUTPUT}")
    print(json.dumps(settings, indent=2))


if __name__ == "__main__":
    main()

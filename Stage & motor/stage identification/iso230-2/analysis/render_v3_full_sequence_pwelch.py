#!/usr/bin/env python3
"""Render a Welch PSD of the complete ISO 230-2 v3 position sequence.

The complete trace is deliberately analysed rather than only the settled
windows.  Consequently, low-frequency power describes the programmed
move/dwell sequence as well as the instrument and stage.  The broadband panel
is split out so that weaker electrical/acquisition lines remain visible.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, welch

from render_raw_measurement import read_scope


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / "ISO230-2_v3.csv"
OUTPUT = HERE / "plots" / "iso230_2_v3_full_sequence_pwelch.png"
SETTINGS_OUTPUT = HERE / "iso230_2_v3_pwelch_settings.json"

NPERSEG = 2**16
OVERLAP = NPERSEG // 2

INK = "#333432"
MUTED = "#73746f"
BLUE = "#587fa3"
RUST = "#a76455"
GRID = "#b9bab5"
SURFACE = "#fbfbfa"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
        "font.size": 9.0,
        "axes.labelsize": 9.5,
        "axes.titlesize": 10.0,
        "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight",
    })


def finish_axis(ax: plt.Axes, *, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_linewidth(0.7)
    ax.tick_params(direction="out", length=3.2, width=0.7)
    if grid:
        ax.grid(True, which="major", color=GRID, alpha=0.48,
                linestyle="--", linewidth=0.55)
        ax.grid(True, which="minor", color=GRID, alpha=0.24,
                linestyle="--", linewidth=0.42)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.115, 1.04, label, transform=ax.transAxes, ha="left",
            va="bottom", fontsize=11.0, fontweight="bold", color=INK)


def select_low_frequency_peaks(freq: np.ndarray, psd: np.ndarray,
                               count: int = 6) -> list[int]:
    eligible = np.flatnonzero((freq >= 0.03) & (freq <= 5.0))
    spectrum_db = 10.0 * np.log10(np.maximum(psd[eligible], np.finfo(float).tiny))
    local, _ = find_peaks(spectrum_db, prominence=4.0, distance=4)
    candidates = eligible[local]
    candidates = candidates[np.argsort(psd[candidates])[::-1]]

    chosen: list[int] = []
    for candidate in candidates:
        # Suppress near-duplicate labels while retaining distinct harmonics.
        if all(abs(np.log(freq[candidate] / freq[old])) > 0.095 for old in chosen):
            chosen.append(int(candidate))
        if len(chosen) == count:
            break
    return sorted(chosen, key=lambda i: freq[i])


def nearest_line(freq: np.ndarray, psd: np.ndarray, nominal_hz: float,
                 half_width_hz: float = 0.35) -> tuple[float, float]:
    band = np.flatnonzero(np.abs(freq - nominal_hz) <= half_width_hz)
    if band.size == 0:
        return float("nan"), float("nan")
    index = int(band[np.argmax(psd[band])])
    return float(freq[index]), float(psd[index])


def main() -> int:
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)

    scope = read_scope(SOURCE)
    sample_rate_hz = 1000.0 / scope.sample_ms
    position_um = scope.displacement_mm * 1000.0
    freq, psd = welch(
        position_um,
        fs=sample_rate_hz,
        window="hann",
        nperseg=NPERSEG,
        noverlap=OVERLAP,
        detrend="constant",
        return_onesided=True,
        scaling="density",
        average="mean",
    )
    df_hz = sample_rate_hz / NPERSEG
    segment_count = 1 + (position_um.size - NPERSEG) // (NPERSEG - OVERLAP)
    low_peaks = select_low_frequency_peaks(freq, psd)
    line_peaks = [nearest_line(freq, psd, value)
                  for value in np.arange(50.0, 500.0, 50.0)]

    set_style()
    fig = plt.figure(figsize=(10.2, 7.1))
    grid = fig.add_gridspec(2, 2, height_ratios=(0.78, 1.25),
                            left=0.085, right=0.975, bottom=0.105, top=0.885,
                            hspace=0.48, wspace=0.28)
    ax_time = fig.add_subplot(grid[0, :])
    ax_low = fig.add_subplot(grid[1, 0])
    ax_high = fig.add_subplot(grid[1, 1])

    stride = max(1, position_um.size // 28000)
    ax_time.plot(scope.time_s[::stride] / 60.0,
                 scope.displacement_mm[::stride], color=MUTED,
                 linewidth=0.62, rasterized=True)
    ax_time.set_xlabel("elapsed time (min)")
    ax_time.set_ylabel("relative position (mm)")
    ax_time.set_title("Complete recorded position sequence", loc="left",
                      color=INK, pad=7)
    ax_time.margins(x=0)
    finish_axis(ax_time, grid=False)
    panel_label(ax_time, "(a)")

    positive = freq > 0.0
    ax_low.loglog(freq[positive], psd[positive], color=BLUE, linewidth=0.9)
    ax_low.set_xlim(df_hz, 5.0)
    ax_low.set_xlabel("frequency (Hz)")
    ax_low.set_ylabel(r"position PSD ($\mathrm{\mu m^2\,Hz^{-1}}$)")
    ax_low.set_title("Sequence and low-frequency content", loc="left", pad=7)
    finish_axis(ax_low)
    panel_label(ax_low, "(b)")

    offsets = [(4, 9), (4, -15), (4, 9), (4, -15), (4, 9), (4, -15)]
    for index, offset in zip(low_peaks, offsets):
        ax_low.plot(freq[index], psd[index], marker="o", markersize=2.8,
                    color=RUST, zorder=4)
        ax_low.annotate(f"{freq[index]:.3f}",
                        xy=(freq[index], psd[index]), xytext=offset,
                        textcoords="offset points", fontsize=7.0, color=INK,
                        ha="left", va="center",
                        arrowprops=dict(arrowstyle="-", color=RUST,
                                        lw=0.45, alpha=0.8))

    ax_high.loglog(freq[positive], psd[positive], color=BLUE, linewidth=0.82)
    ax_high.set_xlim(5.0, sample_rate_hz / 2.0)
    ax_high.set_xlabel("frequency (Hz)")
    ax_high.set_ylabel(r"position PSD ($\mathrm{\mu m^2\,Hz^{-1}}$)")
    ax_high.set_title("Broadband content", loc="left", pad=7)
    finish_axis(ax_high)
    panel_label(ax_high, "(c)")

    for nominal, (measured, value) in zip(np.arange(50.0, 500.0, 50.0),
                                           line_peaks):
        ax_high.axvline(nominal, color=RUST, linewidth=0.58,
                        linestyle=(0, (3, 3)), alpha=0.55, zorder=1)
        if np.isfinite(measured):
            ax_high.plot(measured, value, marker="o", markersize=2.2,
                         color=RUST, zorder=4)
    ax_high.text(0.975, 0.965, "50 Hz harmonic markers", transform=ax_high.transAxes,
                 ha="right", va="top", fontsize=7.5, color=RUST,
                 bbox=dict(boxstyle="round,pad=0.24", facecolor=SURFACE,
                           edgecolor=RUST, linewidth=0.45, alpha=0.92))

    fig.text(0.085, 0.950, "ISO 230-2 v3 | full-sequence Welch spectrum",
             ha="left", va="top", fontsize=14.0, fontweight="bold", color=INK)
    fig.text(0.085, 0.915,
             f"{position_um.size:,} samples at {sample_rate_hz:g} Hz; "
             f"Hann window, {NPERSEG / sample_rate_hz:.3f} s segments, "
             f"50% overlap, constant detrend; resolution {df_hz:.4f} Hz",
             ha="left", va="top", fontsize=8.5, color=MUTED)
    fig.text(0.085, 0.035,
             "Note: the low-frequency spectrum includes the deliberately non-stationary "
             "move/dwell programme; it is not a dwell-only vibration spectrum.",
             ha="left", va="bottom", fontsize=8.0, color=MUTED)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=300)
    plt.close(fig)

    settings = {
        "source": str(SOURCE),
        "output": str(OUTPUT),
        "samples": int(position_um.size),
        "duration_s": float(scope.time_s[-1] - scope.time_s[0]),
        "sample_rate_hz": sample_rate_hz,
        "signal": "relative encoder position",
        "signal_unit": "um",
        "window": "hann",
        "nperseg": NPERSEG,
        "segment_duration_s": NPERSEG / sample_rate_hz,
        "noverlap": OVERLAP,
        "overlap_percent": 50.0,
        "detrend": "constant (per segment)",
        "scaling": "density",
        "average": "mean",
        "frequency_resolution_hz": df_hz,
        "welch_segment_count": int(segment_count),
        "selected_low_frequency_peaks_hz": [float(freq[i]) for i in low_peaks],
        "harmonic_line_peaks_hz": [measured for measured, _ in line_peaks],
    }
    SETTINGS_OUTPUT.write_text(json.dumps(settings, indent=2) + "\n",
                               encoding="utf-8")

    print(f"source: {SOURCE}")
    print(f"samples: {position_um.size:,}; duration: {settings['duration_s']:.3f} s")
    print(f"fs: {sample_rate_hz:g} Hz; df: {df_hz:.6f} Hz; segments: {segment_count}")
    print("low-frequency peaks: " + ", ".join(
        f"{freq[i]:.4f} Hz" for i in low_peaks))
    print("50 Hz-family peaks: " + ", ".join(
        f"{measured:.3f} Hz" for measured, _ in line_peaks))
    print(f"wrote: {OUTPUT}")
    print(f"wrote: {SETTINGS_OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Render the expected Beckhoff 250–1000–250 Hz chirp command."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

TASK_DT_S = 250e-6
START_HZ = 250.0
END_HZ = 1000.0
SWEEP_S = 96.0
SETTLE_S = 5.0
AMPLITUDE_FULL_STEPS = 3.0
TERMINAL_MICROSTEPS = 64
SWEEP_RATE_HZ_S = (END_HZ - START_HZ) / SWEEP_S
TOTAL_S = 3.0 * SETTLE_S + 2.0 * SWEEP_S
T_UP_START = SETTLE_S
T_UP_END = T_UP_START + SWEEP_S
T_DOWN_START = T_UP_END + SETTLE_S
T_DOWN_END = T_DOWN_START + SWEEP_S

OUT_DIR = Path(__file__).resolve().parent / "Rendered"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def render_overview():
    t = np.linspace(0.0, TOTAL_S, 12000)
    frequency = np.zeros_like(t)
    amplitude = np.zeros_like(t)

    up = (t >= T_UP_START) & (t < T_UP_END)
    down = (t >= T_DOWN_START) & (t < T_DOWN_END)
    frequency[up] = START_HZ + SWEEP_RATE_HZ_S * (t[up] - T_UP_START)
    frequency[down] = END_HZ - SWEEP_RATE_HZ_S * (t[down] - T_DOWN_START)
    amplitude[up | down] = AMPLITUDE_FULL_STEPS

    fig, (ax_f, ax_x) = plt.subplots(
        2, 1, figsize=(14, 7.5), sharex=True, constrained_layout=True
    )
    ax_f.plot(t, frequency, color="#2463eb", linewidth=2.2)
    ax_f.fill_between(t, 0, frequency, color="#2463eb", alpha=0.12)
    ax_f.set_ylabel("Command frequency (Hz)")
    ax_f.set_ylim(0, 1050)
    ax_f.grid(True, alpha=0.25)
    ax_f.set_title("Beckhoff unnotched chirp — configured run overview")

    ax_x.fill_between(t, -amplitude, amplitude, where=amplitude > 0,
                      color="#0e9f6e", alpha=0.25, step="mid")
    ax_x.plot(t, amplitude, color="#0e9f6e", linewidth=1.4)
    ax_x.plot(t, -amplitude, color="#0e9f6e", linewidth=1.4)
    ax_x.axhline(0.0, color="#111827", linewidth=0.8)
    ax_x.set_ylabel("Command envelope\n(full steps)")
    ax_x.set_xlabel("Elapsed time (seconds)")
    ax_x.set_ylim(-3.5, 3.5)
    ax_x.grid(True, alpha=0.25)

    segments = (
        (0.0, T_UP_START, "Settle", "#6b7280"),
        (T_UP_START, T_UP_END, "Up: 250 → 1000 Hz", "#2463eb"),
        (T_UP_END, T_DOWN_START, "Settle", "#6b7280"),
        (T_DOWN_START, T_DOWN_END, "Down: 1000 → 250 Hz", "#9333ea"),
        (T_DOWN_END, TOTAL_S, "Settle", "#6b7280"),
    )
    for start, stop, label, color in segments:
        for ax in (ax_f, ax_x):
            ax.axvspan(start, stop, color=color, alpha=0.045)
        if stop - start > 8:
            ax_f.text((start + stop) / 2, 1020, label, ha="center", va="top",
                      fontsize=9.5, color=color, fontweight="bold")

    ax_x.text(
        0.01, 0.95,
        "±3 full-step peak envelope • no notch • no commanded frequency below 250 Hz",
        transform=ax_x.transAxes, va="top", fontsize=10
    )
    fig.savefig(OUT_DIR / "chirp_sequence_overview.png", dpi=180)
    plt.close(fig)


def sampled_chirp(local_t, start_hz, sweep_rate):
    phase = 2.0 * np.pi * (
        start_hz * local_t + 0.5 * sweep_rate * local_t * local_t
    )
    counts = np.rint(
        AMPLITUDE_FULL_STEPS * TERMINAL_MICROSTEPS * np.sin(phase)
    )
    return counts / TERMINAL_MICROSTEPS


def render_sampling_detail():
    cases = (
        ("Start of up-sweep", START_HZ, SWEEP_RATE_HZ_S, 0.012,
         "250 Hz: 16 command samples per period", "#2463eb"),
        ("Start of down-sweep", END_HZ, -SWEEP_RATE_HZ_S, 0.004,
         "1000 Hz: 4 command samples per period", "#9333ea"),
    )
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.5), constrained_layout=True)

    for ax, (title, start_hz, rate, duration, note, color) in zip(axes, cases):
        sample_t = np.arange(0.0, duration + 0.5 * TASK_DT_S, TASK_DT_S)
        sample_x = sampled_chirp(sample_t, start_hz, rate)
        fine_t = np.linspace(0.0, duration, 6000)
        fine_phase = 2.0 * np.pi * (
            start_hz * fine_t + 0.5 * rate * fine_t * fine_t
        )
        ideal_x = AMPLITUDE_FULL_STEPS * np.sin(fine_phase)

        ax.plot(fine_t * 1000.0, ideal_x, linestyle="--", color="#111827",
                alpha=0.65, linewidth=1.2, label="Ideal analytic sine")
        ax.step(sample_t * 1000.0, sample_x, where="post", color=color,
                linewidth=2.0, label="250 µs Beckhoff position command")
        ax.scatter(sample_t * 1000.0, sample_x, color=color, s=24, zorder=3)
        ax.axhline(0.0, color="#6b7280", linewidth=0.7)
        ax.set_ylim(-3.5, 3.5)
        ax.set_ylabel("Position (full steps)")
        ax.set_title(f"{title} — {note}", loc="left")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right")

    axes[-1].set_xlabel("Local sweep time (milliseconds)")
    fig.suptitle("Chirp waveform detail: analytic target versus cyclic setpoints", fontsize=14)
    fig.savefig(OUT_DIR / "chirp_sampling_detail.png", dpi=180)
    plt.close(fig)


def main():
    render_overview()
    render_sampling_detail()
    print(f"Rendered chirp previews to {OUT_DIR}")
    print(f"Total configured duration: {TOTAL_S:.1f} s")
    print(f"Samples per 1 kHz period: {1.0 / (TASK_DT_S * END_HZ):.1f}")


if __name__ == "__main__":
    main()
"""Render the expected Beckhoff oscillating-stepping command sequence."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

MRES_LEVELS = (1, 2, 4, 8, 16, 32, 64)
TERMINAL_MICROSTEPS = 64
CYCLES_PER_LEVEL = 15
ENDPOINT_DWELL_S = 5.0
LEAD_SETTLE_S = 5.0
INTER_LEVEL_SETTLE_S = 5.0
TAIL_SETTLE_S = 5.0
FULL_STEP_UM = 10.0

OUT_DIR = Path(__file__).resolve().parent / "Rendered"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def build_sequence():
    times = [0.0, LEAD_SETTLE_S]
    command_full_steps = [0.0, 0.0]
    level_spans = []
    t = LEAD_SETTLE_S

    for level_index, mres in enumerate(MRES_LEVELS):
        level_start = t
        amplitude = 1.0 / mres
        for _ in range(CYCLES_PER_LEVEL):
            times.append(t)
            command_full_steps.append(amplitude)
            t += ENDPOINT_DWELL_S
            times.append(t)
            command_full_steps.append(0.0)
            t += ENDPOINT_DWELL_S
            times.append(t)
            command_full_steps.append(0.0)
        level_spans.append((level_start, t, mres, amplitude))
        if level_index < len(MRES_LEVELS) - 1:
            t += INTER_LEVEL_SETTLE_S
            times.append(t)
            command_full_steps.append(0.0)

    t += TAIL_SETTLE_S
    times.append(t)
    command_full_steps.append(0.0)
    return np.asarray(times), np.asarray(command_full_steps), level_spans, t


def render_overview(times, command, spans, total_s):
    colors = plt.cm.tab10(np.linspace(0.0, 0.9, len(spans)))
    fig, (ax_cmd, ax_level) = plt.subplots(
        2, 1, figsize=(14, 7.5), sharex=True,
        gridspec_kw={"height_ratios": (3.0, 1.0)}, constrained_layout=True
    )

    for color, (start, stop, mres, _) in zip(colors, spans):
        ax_cmd.axvspan(start / 60.0, stop / 60.0, color=color, alpha=0.08)
        ax_level.broken_barh(
            [(start / 60.0, (stop - start) / 60.0)], (0.15, 0.7),
            facecolors=color, alpha=0.9
        )
        ax_level.text(
            (start + stop) / 120.0, 0.5, f"MRES {mres}",
            ha="center", va="center", color="white", fontweight="bold"
        )

    ax_cmd.step(times / 60.0, command * FULL_STEP_UM, where="post",
                color="#111827", linewidth=1.0)
    ax_cmd.set_title("Beckhoff oscillating-stepping run — complete command sequence")
    ax_cmd.set_ylabel("Commanded displacement (µm)")
    ax_cmd.set_ylim(-0.5, 10.8)
    ax_cmd.grid(True, alpha=0.25)
    ax_cmd.text(
        0.01, 0.95,
        "15 forward/return cycles per level • 5 s at each endpoint • no long trajectories",
        transform=ax_cmd.transAxes, va="top", fontsize=10
    )

    ax_level.set_ylim(0, 1)
    ax_level.set_yticks([])
    ax_level.set_xlabel("Elapsed time (minutes)")
    ax_level.set_xlim(0, total_s / 60.0)
    ax_level.set_title("Requested microstep-resolution blocks")
    ax_level.grid(True, axis="x", alpha=0.2)

    fig.savefig(OUT_DIR / "stepping_sequence_overview.png", dpi=180)
    plt.close(fig)


def main():
    times, command, spans, total_s = build_sequence()
    render_overview(times, command, spans, total_s)
    print(f"Rendered stepping previews to {OUT_DIR}")
    print(f"Total configured duration: {total_s:.1f} s ({total_s / 60.0:.2f} min)")


if __name__ == "__main__":
    main()
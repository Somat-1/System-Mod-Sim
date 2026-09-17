"""Preview/export the seeded irregular stepping command sequence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DEFAULT_SEED = 0x6D2B79F5
DEFAULT_RANDOM_MOVES = 16
TOTAL_RANGE_MM = 20.0
MIN_ABS_TARGET_MM = 0.5
FULL_STEP_MM = 0.010
MICROSTEPS_PER_FULL_STEP = 64
INITIAL_SETTLE_S = 5.0
DWELL_S = 5.0
MASK32 = 0xFFFFFFFF

OUT_DIR = Path(__file__).resolve().parent / "Rendered"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def xorshift32(state: int) -> int:
    state ^= (state << 13) & MASK32
    state ^= state >> 17
    state ^= (state << 5) & MASK32
    return state & MASK32


def build_commands(seed: int, random_moves: int, total_range_mm: float):
    if seed == 0:
        seed = DEFAULT_SEED
    if not 1 <= random_moves <= 30:
        raise ValueError("random_moves must be in [1, 30]")
    half_range = 0.5 * total_range_mm
    if half_range <= MIN_ABS_TARGET_MM:
        raise ValueError("half-range must exceed the minimum target magnitude")

    state = seed & MASK32
    target_offsets = []
    relative_counts = []
    for index in range(1, random_moves + 1):
        state = xorshift32(state)
        fraction = (state % 1_000_000) / 999_999.0
        magnitude = MIN_ABS_TARGET_MM + fraction * (half_range - MIN_ABS_TARGET_MM)
        requested = magnitude if index % 2 == 1 else -magnitude
        counts = int(np.rint(
            requested / FULL_STEP_MM * MICROSTEPS_PER_FULL_STEP
        ))
        quantized = counts * FULL_STEP_MM / MICROSTEPS_PER_FULL_STEP
        target_offsets.append(quantized)
        relative_counts.append(counts)

    target_offsets.append(0.0)  # final return to captured origin
    relative_counts.append(0)
    target_offsets = np.asarray(target_offsets, dtype=float)
    relative_counts = np.asarray(relative_counts, dtype=int)
    previous = np.concatenate(([0.0], target_offsets[:-1]))
    move_distances = target_offsets - previous
    nominal_times = INITIAL_SETTLE_S + DWELL_S * np.arange(len(target_offsets))
    return target_offsets, move_distances, relative_counts, nominal_times


def write_csv(seed, target_offsets, move_distances, relative_counts, nominal_times):
    path = OUT_DIR / f"irregular_commands_seed_{seed:08X}.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "command_index", "seed_hex", "target_offset_mm",
            "signed_move_distance_mm", "relative_terminal_counts",
            "nominal_command_time_s", "is_final_origin_return",
        ])
        for i, (target, move, counts, time_s) in enumerate(
            zip(target_offsets, move_distances, relative_counts, nominal_times), start=1
        ):
            writer.writerow([
                i, f"0x{seed:08X}", f"{target:.9f}", f"{move:.9f}",
                counts, f"{time_s:.3f}", int(i == len(target_offsets)),
            ])
    return path


def render(seed, target_offsets, move_distances, nominal_times, total_range_mm):
    indices = np.arange(1, len(target_offsets) + 1)
    colors = np.where(target_offsets > 0, "#2463eb",
                      np.where(target_offsets < 0, "#d94645", "#111827"))
    half_range = 0.5 * total_range_mm

    fig, (ax_target, ax_move) = plt.subplots(
        2, 1, figsize=(14, 8.5), sharex=True,
        gridspec_kw={"height_ratios": (2.1, 1.2)}, constrained_layout=True
    )
    ax_target.plot(indices, target_offsets, color="#6b7280", linewidth=1.3)
    ax_target.scatter(indices, target_offsets, c=colors, s=55, zorder=3)
    ax_target.axhline(0.0, color="#111827", linewidth=0.8)
    ax_target.axhline(half_range, color="#d97706", linestyle="--", linewidth=1.1)
    ax_target.axhline(-half_range, color="#d97706", linestyle="--", linewidth=1.1)
    ax_target.fill_between(indices, -half_range, half_range, color="#f59e0b", alpha=0.05)
    ax_target.set_ylabel("Target offset from origin (mm)")
    ax_target.set_ylim(-half_range * 1.12, half_range * 1.12)
    ax_target.grid(True, alpha=0.25)
    ax_target.set_title(
        f"Irregular stepping preview — seed 0x{seed:08X}, "
        f"{total_range_mm:g} mm total centered span"
    )
    ax_target.text(
        0.01, 0.96,
        "Alternating sides • 5 s dwell after arrival • final command returns to origin",
        transform=ax_target.transAxes, va="top", fontsize=10
    )

    ax_move.bar(indices, move_distances, color=colors, alpha=0.85)
    ax_move.axhline(0.0, color="#111827", linewidth=0.8)
    ax_move.set_ylabel("Signed move distance (mm)")
    ax_move.set_xlabel("Command index")
    ax_move.set_xticks(indices)
    ax_move.grid(True, axis="y", alpha=0.25)

    top = ax_target.secondary_xaxis("top")
    top.set_xticks(indices)
    top.set_xticklabels([f"{t:.0f}" for t in nominal_times], fontsize=8)
    top.set_xlabel("Nominal command time (s), excluding motion-to-target time")

    path = OUT_DIR / f"irregular_stepping_seed_{seed:08X}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=lambda value: int(value, 0), default=DEFAULT_SEED)
    parser.add_argument("--moves", type=int, default=DEFAULT_RANDOM_MOVES)
    parser.add_argument("--total-range-mm", type=float, default=TOTAL_RANGE_MM)
    args = parser.parse_args()

    targets, moves, counts, times = build_commands(
        args.seed, args.moves, args.total_range_mm
    )
    csv_path = write_csv(args.seed, targets, moves, counts, times)
    png_path = render(args.seed, targets, moves, times, args.total_range_mm)
    print(f"Rendered {png_path}")
    print(f"Wrote {csv_path}")
    print(f"Commands: {len(targets)} ({args.moves} irregular + final origin return)")
    print(f"Target range: {targets.min():.6f} to {targets.max():.6f} mm")


if __name__ == "__main__":
    main()
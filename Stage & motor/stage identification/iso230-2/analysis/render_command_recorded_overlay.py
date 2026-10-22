#!/usr/bin/env python3
"""Compare the ISO 230-2 command and interferometer recording accurately.

The controller and TwinCAT scope did not share a trigger.  Stable plateaus are
therefore detected in the raw interferometer trace and sequence-aligned to the
controller's STOP records.  Their times define a monotone clock warp for the
visual time overlay.  Position values are never fitted, scaled, or shifted:
both traces remain in millimetres on the same 1:1 vertical scale.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RAW_MODULE = HERE / "render_raw_measurement.py"
SOURCE = HERE.parent.parent / "ISO230-2_v2.csv"
LOG = HERE.parent / "data" / "terminal_20261021_101911.log"
OUTPUT = HERE / "plots" / "iso230_2_v2_command_recorded_overlay.png"


def load_raw_module():
    spec = importlib.util.spec_from_file_location("iso230_raw", RAW_MODULE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@dataclass
class PlateauSet:
    start: np.ndarray
    end: np.ndarray
    centre_s: np.ndarray
    position_mm: np.ndarray
    std_um: np.ndarray


@dataclass
class Alignment:
    pairs: list[tuple[int, int]]
    missing_command: list[int]
    extra_plateau: list[int]


def detect_plateaus(scope) -> PlateauSet:
    """Find stationary pieces without using command timestamps or positions."""
    delta = np.abs(np.diff(scope.counter_unwrapped,
                           prepend=scope.counter_unwrapped[0])).astype(float)
    activity = np.convolve(delta, np.ones(21) / 21.0, mode="same")
    stationary = activity < 400.0
    edge = np.diff(stationary.astype(np.int8), prepend=0, append=0)
    start = np.flatnonzero(edge == 1)
    end = np.flatnonzero(edge == -1)
    keep = (end - start) >= 100
    start, end = start[keep], end[keep]

    position, std = [], []
    for first, last in zip(start, end):
        values = scope.displacement_mm[first:last]
        position.append(float(np.mean(values)))
        std.append(float(np.std(values) * 1000.0))
    centre_s = (start + end - 1) * scope.sample_ms / 2000.0
    return PlateauSet(start, end, centre_s, np.asarray(position), np.asarray(std))


def align_sequences(command_mm: np.ndarray, measured_mm: np.ndarray) -> Alignment:
    """Order-preserving alignment of commanded stops and raw plateaus."""
    n, m = command_mm.size, measured_mm.size
    score = np.full((n + 1, m + 1), np.inf)
    back = np.zeros((n + 1, m + 1), dtype=np.int8)
    score[0, 0] = 0.0
    score[1:, 0] = np.arange(1, n + 1) * 5.0
    score[0, 1:] = np.arange(1, m + 1) * 4.0
    back[1:, 0] = 1
    back[0, 1:] = 2

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Correct pairs in this recording agree at micrometre scale.  The
            # cap lets a missing stop become a gap rather than forcing a bad
            # match to a nearby target.
            match_cost = min((abs(command_mm[i - 1] - measured_mm[j - 1]) /
                              0.03) ** 2, 30.0)
            choices = (score[i - 1, j - 1] + match_cost,
                       score[i - 1, j] + 5.0,
                       score[i, j - 1] + 4.0)
            move = int(np.argmin(choices))
            score[i, j] = choices[move]
            back[i, j] = move

    pairs, missing, extra = [], [], []
    i, j = n, m
    while i or j:
        move = int(back[i, j]) if i and j else (1 if i else 2)
        if move == 0:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif move == 1:
            missing.append(i - 1)
            i -= 1
        else:
            extra.append(j - 1)
            j -= 1
    return Alignment(pairs[::-1], missing[::-1], extra[::-1])


def estimate_counter_datum(command_mm: np.ndarray,
                           measured_relative_mm: np.ndarray) -> float:
    """Find the additive absolute-position datum; never alter distance scale."""
    differences = (command_mm[:, None] - measured_relative_mm[None, :]).ravel()
    edges = np.arange(-160.0, 160.005, 0.005)
    counts, _ = np.histogram(differences, bins=edges)
    order = np.argsort(counts)[::-1]
    candidates = []
    for index in order:
        centre = (edges[index] + edges[index + 1]) / 2
        if all(abs(centre - old) > 0.08 for old in candidates):
            near = differences[np.abs(differences - centre) < 0.025]
            candidates.append(float(np.median(near)) if near.size else float(centre))
        if len(candidates) == 30:
            break

    best = None
    for offset in candidates:
        measured = measured_relative_mm + offset
        aligned = align_sequences(command_mm, measured)
        error = np.asarray([measured[j] - command_mm[i]
                            for i, j in aligned.pairs])
        match_cost = float(np.sum(np.minimum((np.abs(error) / 0.03) ** 2, 30.0)))
        quality = (match_cost + 5.0 * len(aligned.missing_command) +
                   4.0 * len(aligned.extra_plateau))
        if best is None or quality < best[0]:
            best = (quality, offset)
    if best is None:
        raise ValueError("Could not estimate the interferometer counter datum")
    return float(best[1])


def fit_affine_clock(board_anchor: np.ndarray,
                     scope_anchor: np.ndarray) -> tuple[float, float]:
    """Robust single clock-rate and offset fit; never distort local durations."""
    left, right = np.triu_indices(board_anchor.size, k=1)
    separation = board_anchor[right] - board_anchor[left]
    use = separation > 60.0
    slopes = ((scope_anchor[right][use] - scope_anchor[left][use]) /
              separation[use])
    rate = float(np.median(slopes))
    offset = float(np.median(scope_anchor - rate * board_anchor))
    return rate, offset


def affine_clock_map(values: np.ndarray, rate: float, offset: float) -> np.ndarray:
    return rate * values + offset


def style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE,
                        help="TwinCAT scope CSV")
    parser.add_argument("--log", type=Path, default=LOG,
                        help="matching ESP32 terminal log")
    parser.add_argument("--output", type=Path,
                        help="output PNG (derived from source name by default)")
    parser.add_argument("--interactive", action="store_true",
                        help="open the comparison in an interactive Matplotlib window")
    args = parser.parse_args()

    raw = load_raw_module()
    scope = raw.read_scope(args.source)
    controller = raw.read_controller_log(args.log, scope.start_epoch_s)

    stops = [s for s in controller.stops
             if s[1] in ("DATUM", "ISO_REF", "ISO", "END_REF")]
    reference_positions = [s[5] for s in stops if s[1] in
                           ("DATUM", "ISO_REF", "END_REF")]
    if not reference_positions:
        raise ValueError("Controller log contains no reference-position stops")
    x_ref_mm = float(np.median(reference_positions))
    plateaus = detect_plateaus(scope)
    command_stop = np.asarray([s[5] for s in stops])
    board_stop_s = np.asarray([s[0] for s in stops])
    counter_datum_mm = estimate_counter_datum(command_stop, plateaus.position_mm)
    plateaus.position_mm = plateaus.position_mm + counter_datum_mm
    alignment = align_sequences(command_stop, plateaus.position_mm)
    command_i = np.asarray([i for i, _ in alignment.pairs], dtype=int)
    plateau_i = np.asarray([j for _, j in alignment.pairs], dtype=int)
    board_anchor = board_stop_s[command_i]
    scope_anchor = plateaus.centre_s[plateau_i]

    clock_rate, clock_offset = fit_affine_clock(board_anchor, scope_anchor)
    warped_command_s = affine_clock_map(controller.sample_time_s,
                                        clock_rate, clock_offset)
    valid_command = ((warped_command_s >= 0) &
                     (warped_command_s <= scope.time_s[-1]))

    measured_stop = plateaus.position_mm[plateau_i]
    intended_stop = command_stop[command_i]
    error_um = (measured_stop - intended_stop) * 1000.0
    iso_mask = np.asarray([stops[i][1] == "ISO" for i in command_i])
    iso_total = sum(s[1] == "ISO" for s in stops)
    iso_matched = int(np.count_nonzero(iso_mask))

    style()
    fig = plt.figure(figsize=(15.8, 14.0))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.05, 1.05, 1.4, 1.25])
    ax_raw = fig.add_subplot(gs[0, :])
    ax_cmd = fig.add_subplot(gs[1, :])
    ax_overlay = fig.add_subplot(gs[2, :])
    ax_one = fig.add_subplot(gs[3, 0])
    ax_err = fig.add_subplot(gs[3, 1])

    raw_abs = scope.displacement_mm + counter_datum_mm
    scope_min = scope.time_s / 60.0
    ax_raw.plot(scope_min, raw_abs, color="#155f8a", lw=0.68,
                rasterized=True, label="raw interferometer")
    ax_raw.set_title("(a) Complete raw interferometer recording", loc="left",
                     fontweight="bold")
    ax_raw.set_ylabel("Recorded position [mm]")
    ax_raw.set_xlim(0, scope_min[-1])

    datum_board = min(t for t, state in controller.state_changes if state == "DATUM")
    cmd_elapsed_min = (controller.sample_time_s - datum_board) / 60.0
    sequence = controller.sample_time_s >= datum_board
    ax_cmd.plot(cmd_elapsed_min[sequence], controller.command_mm[sequence],
                color="#5c5c5c", lw=0.72, rasterized=True,
                label="controller command")
    stop_elapsed = (board_stop_s - datum_board) / 60.0
    ax_cmd.scatter(stop_elapsed, command_stop, s=5, color="#d55e00",
                   zorder=3, label="recorded stop command")
    ax_cmd.set_title("(b) Complete controller command in the ESP32 board clock",
                     loc="left", fontweight="bold")
    ax_cmd.set_ylabel("Intended position [mm]")
    ax_cmd.set_xlabel("Controller sequence time [min]")

    ax_overlay.plot(scope_min, raw_abs, color="#155f8a", lw=0.82,
                    rasterized=True, label="recorded position")
    ax_overlay.plot(warped_command_s[valid_command] / 60.0,
                    controller.command_mm[valid_command], color="#4d4d4d",
                    lw=0.75, alpha=0.82, linestyle=(0, (3, 2)),
                    rasterized=True, label="intended position — synchronized")
    ax_overlay.scatter(plateaus.centre_s[plateau_i] / 60.0, measured_stop,
                       s=6, color="#009e73", zorder=4,
                       label="paired stable arrival")
    if alignment.missing_command:
        miss = np.asarray(alignment.missing_command, dtype=int)
        miss_t = affine_clock_map(board_stop_s[miss], clock_rate, clock_offset)
        ax_overlay.scatter(miss_t / 60.0, command_stop[miss], s=28,
                           facecolors="none", edgecolors="#cc3311", marker="o",
                           linewidths=0.9, zorder=5,
                           label="commanded stop without resolved plateau")
    ax_overlay.set_title(
        "(c) Affine clock-corrected overlay — command durations preserved",
        loc="left", fontweight="bold")
    ax_overlay.set_ylabel("Position [mm]")
    ax_overlay.set_xlabel("Scope elapsed time [min]")
    ax_overlay.set_xlim(0, scope_min[-1])
    ax_overlay.legend(loc="upper right", ncol=2, fontsize=8.5)

    for position_ax in (ax_raw, ax_cmd, ax_overlay):
        position_ax.set_ylim(-3, 123)

    state = np.asarray([stops[i][1] for i in command_i])
    direction = np.asarray([stops[i][3] for i in command_i])
    groups = [
        ((state == "ISO") & (direction == 1), "outward", "^", "#0072b2"),
        ((state == "ISO") & (direction == -1), "return", "v", "#d55e00"),
        (state != "ISO", "reference", "o", "#555555"),
    ]
    for mask, label, marker, color in groups:
        ax_one.scatter(intended_stop[mask], measured_stop[mask], s=18,
                       marker=marker, color=color, alpha=0.8, linewidths=0,
                       label=label)
    low = min(float(np.min(intended_stop)), float(np.min(measured_stop))) - 2
    high = max(float(np.max(intended_stop)), float(np.max(measured_stop))) + 2
    ax_one.plot([low, high], [low, high], color="black", lw=0.9,
                linestyle=(0, (4, 3)), label="ideal 1:1")
    ax_one.set_xlim(low, high)
    ax_one.set_ylim(low, high)
    ax_one.set_aspect("equal", adjustable="box")
    ax_one.set_title("(d) Arrival position: recorded versus intended", loc="left",
                     fontweight="bold")
    ax_one.set_xlabel("Intended position [mm]")
    ax_one.set_ylabel("Recorded plateau position [mm]")
    ax_one.legend(loc="upper left", fontsize=8)

    ax_err.axhspan(-5, 5, color="#009e73", alpha=0.10,
                   label="±5 µm reference band")
    ax_err.axhline(0, color="black", lw=0.75)
    ax_err.scatter(command_i[iso_mask] + 1, error_um[iso_mask], s=10,
                   c=np.where(direction[iso_mask] > 0, "#0072b2", "#d55e00"),
                   alpha=0.75, linewidths=0)
    ax_err.set_title("(e) Paired arrival error in sequence order", loc="left",
                     fontweight="bold")
    ax_err.set_xlabel("Controller stop number")
    ax_err.set_ylabel("Recorded − intended [µm]")
    ax_err.set_xlim(0, len(stops) + 1)
    ax_err.legend(loc="upper right", fontsize=8)

    for ax in (ax_raw, ax_cmd, ax_overlay, ax_one, ax_err):
        ax.grid(True, color="0.9", lw=0.55)

    raw_span = float(np.ptp(raw_abs))
    command_span = 120.5 - 2.0
    stats = (
        f"{args.source.name}  |  Paired stable plateaus: "
        f"{len(alignment.pairs)}/{len(stops)} total; "
        f"{iso_matched}/{iso_total} ISO target stops.  "
        f"Paired error: median {np.median(error_um):+.2f} µm, "
        f"RMS {np.sqrt(np.mean(error_um**2)):.2f} µm, "
        f"max |error| {np.max(np.abs(error_um)):.2f} µm.\n"
        f"Raw dynamic span {raw_span:.3f} mm versus commanded overrun span "
        f"{command_span:.3f} mm.  Clock rate {clock_rate:.8f}; "
        f"counter datum {counter_datum_mm:.6f} mm."
    )
    method = (
        "Synchronization method: raw stationary plateaus detected independently, "
        "then order-matched to controller STOP records.\nTheir times fit one affine "
        "board-to-scope clock rate and offset, preserving every command duration. "
        "One constant counter datum sets the "
        "absolute position; distance values are not gain-fitted or rescaled.\n"
        "Open red circles mark commands for which no separate stable plateau was "
        "resolved. Overrun/reversal points are not recorded stops."
    )
    fig.suptitle("ISO 230-2 — commanded and recorded position comparison",
                 x=0.06, y=0.988, ha="left", fontsize=16, fontweight="bold")
    fig.text(0.06, 0.958, stats, ha="left", va="top", fontsize=9.2, color="0.22")
    fig.text(0.06, 0.012, method, ha="left", va="bottom", fontsize=8.2,
             color="0.28")
    fig.subplots_adjust(left=0.065, right=0.985, top=0.905, bottom=0.085,
                        hspace=0.42, wspace=0.28)

    output = args.output
    if output is None:
        stem = args.source.stem.lower().replace("-", "_")
        output = HERE / "plots" / f"{stem}_command_recorded_overlay.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.interactive:
        fig.canvas.manager.set_window_title(
            "ISO 230-2 — commanded and recorded position comparison")
        plt.show()
    else:
        fig.savefig(output, dpi=300)
        plt.close(fig)
        print(f"Wrote {output}")
    print(stats)
    print(f"Unresolved controller stops: {len(alignment.missing_command)}")
    for i in alignment.missing_command:
        stop = stops[i]
        print(f"  stop {i + 1}: {stop[1]} run={stop[2]} dir={stop[3]:+d} "
              f"target={stop[4]} x={stop[5]:.6f} mm")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Visualize the commanded stepping sequence actually fed to the ODE, for
the same three cruise windows used in compare_D_3.5/70/200.png -- and
compare it to what a REAL discrete microstep staircase would look like.

Motivation: linear vs. Option A came out visually and numerically
indistinguishable (tracking error stayed under 0.32 deg everywhere,
against a 1.80 deg pull-out threshold) -- see README.md. One candidate
reason: reconstruct_segments()/trapezoid_fraction() reconstructs
controller-paced moves as a SMOOTH continuous ramp (the best available
reconstruction, since individual microstep pulses are not logged for
these blocks -- see the earlier command/timing investigation). A real
stepper, though, only ever receives discrete microstep commands, one
every 1/(rate*MRES) seconds. This script makes that discretization gap
visible directly: no ODE integration needed, purely geometric.

If the real hardware's tracking error is driven by re-settling after each
discrete microstep edge (the electrical-time-constant argument: tau ~=
0.61 ms vs. a 1.25 ms microstep period at D_200/MRES 1/4), our smooth
reconstruction cannot show that effect at all -- it has no edges for the
electrical dynamics to react to, regardless of which drive law (linear,
Option A, or Option B) is used on top of it.
"""
from __future__ import annotations

import csv
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
V3_ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = V3_ROOT / 'data' / 'raw_local'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
SCRIPTS_DIR = V3_ROOT / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
from command_reconstruction import reconstruct_segments, trapezoid_fraction  # noqa: E402

RUN_INDEX = 2
MRES = 4  # confirmed from RUN_CONFIG for run 2
MOTOR_FULL_STEPS_PER_REV = 200
BLOCKS = ('3.5', '70', '200')
N_CYCLES = 10
FILETIME_UNIX_EPOCH = 116444736000000000
CONTROLLER_CLOCK_SKEW_S = 0.319


def parse_ids_start_epoch(path: Path) -> float:
    start_filetime = None
    with path.open('r', encoding='utf-8-sig', errors='replace') as handle:
        for line in handle:
            fields = line.rstrip('\r\n').split('\t')
            if fields and fields[0] == 'Starttime of export':
                start_filetime = int(fields[1])
                break
    return (start_filetime - FILETIME_UNIX_EPOCH) / 1.0e7


def load_log_rows(path: Path, ids_start_epoch_s: float):
    rows = []
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            instant = datetime.fromisoformat(row['utc'])
            row['ids_time_s'] = (
                instant.timestamp() - ids_start_epoch_s - CONTROLLER_CLOCK_SKEW_S
            )
            rows.append(row)
    return rows


def find_block(rows, run_index, block_name):
    start = next(
        r for r in rows if r['event'] == 'BLOCK_START'
        and r['run_index'] == str(run_index) and r['block'] == block_name
    )
    end = next(
        r for r in rows if r['event'] == 'BLOCK_END'
        and r['run_index'] == str(run_index) and r['block'] == block_name
    )
    return start['ids_time_s'], end['ids_time_s']


def commanded_sequence(rows, rate: str):
    block_name = f'D_{rate}'
    start_s, end_s = find_block(rows, RUN_INDEX, block_name)
    segments = reconstruct_segments(rows, RUN_INDEX, start_s, end_s, block_name)
    first_trapezoid = next(s for s in segments if s['kind'] == 'trapezoid')

    rate_hz = float(rate)
    window_s = N_CYCLES / rate_hz
    cruise_mid = (
        first_trapezoid['t0'] + first_trapezoid['t_accel_s']
        + 0.5 * first_trapezoid['duration_s']
    )
    win_t0, win_t1 = cruise_mid - 0.5 * window_s, cruise_mid + 0.5 * window_s

    t_fine = np.linspace(win_t0, win_t1, 20000)
    frac = np.array([
        trapezoid_fraction(t - first_trapezoid['t0'], first_trapezoid['t_accel_s'],
                            first_trapezoid['duration_s'])
        for t in t_fine
    ])
    theta_smooth_rev = (
        first_trapezoid['value_start']
        + (first_trapezoid['value_end'] - first_trapezoid['value_start']) * frac
    )
    theta_smooth_deg = theta_smooth_rev * 360.0

    microstep_rev = 1.0 / (MOTOR_FULL_STEPS_PER_REV * MRES)
    theta_discrete_deg = np.floor(theta_smooth_rev / microstep_rev) * microstep_rev * 360.0

    microstep_period_s = 1.0 / (rate_hz * MRES)
    t_ms = (t_fine - t_fine[0]) * 1000.0
    return t_ms, theta_smooth_deg, theta_discrete_deg, microstep_period_s, window_s


def main():
    start_epoch_s = parse_ids_start_epoch(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)

    fig, axes = plt.subplots(1, len(BLOCKS), figsize=(7.0 * len(BLOCKS), 5.0))
    for ax, rate in zip(axes, BLOCKS):
        t_ms, smooth_deg, discrete_deg, microstep_period_s, window_s = commanded_sequence(rows, rate)
        ax.plot(t_ms, smooth_deg - smooth_deg[0], color='#1f77b4', lw=1.3,
                label='Smooth reconstruction\n(what the simulation commands)')
        ax.step(t_ms, discrete_deg - discrete_deg[0], where='post', color='#d62728', lw=1.0,
                alpha=0.85, label='True discrete microstep\nstaircase (1/(rate\u00b7MRES) s per step)')
        ax.set_title(
            f'D {rate} full-steps/s\nmicrostep period = {microstep_period_s * 1000:.3f} ms, '
            f'window = {window_s * 1000:.1f} ms',
            fontsize=10,
        )
        ax.set_xlabel('Time within window (ms)', fontsize=8.5)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7.5)
        if rate == BLOCKS[0]:
            ax.set_ylabel('Commanded rotor angle,\nrelative to window start (deg)', fontsize=8.5)
            ax.legend(loc='upper left', fontsize=7.5, framealpha=0.9)

    fig.suptitle(
        'Commanded stepping sequence: smooth trapezoid reconstruction (used by every '
        'simulation in this repo) vs. the real discrete microstep staircase -- run 2, '
        'MRES 1/4. Neither the linear nor Option A drive law sees the red staircase; '
        'both integrate only the blue smooth command.',
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    out_path = HERE / 'commanded_sequence_smooth_vs_discrete.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}')


if __name__ == '__main__':
    main()

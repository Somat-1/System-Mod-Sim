#!/usr/bin/env python3
"""Per-configuration 10-panel block montage: independent-y-scale D
plateaus, one C panel, and BLOCK_0 start/end overlaid.

Markers are excluded from every analysis panel -- each block's own
BLOCK_START/BLOCK_END timestamps bound the plotted window exactly, with
no marker lead-in/lead-out included.
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / 'data' / 'raw_local'
SUBPLOT_DIR = ROOT / 'rendered_assets' / 'individual_subplots'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'

FILETIME_UNIX_EPOCH = 116444736000000000
LEAD_PITCH_M = 2.0e-3
CURRENT_SHORT = {'I_50pct': '50% I', 'I_100pct': '100% I'}
D_RATES = ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')
MEASURED_COLOR = '#136f63'
COMMAND_COLOR = '#d1495b'


def parse_ids(path: Path):
    start_filetime = None
    sample_period_ms = None
    data_line = None
    with path.open('r', encoding='utf-8-sig', errors='replace') as handle:
        for line_number, line in enumerate(handle):
            fields = line.rstrip('\r\n').split('\t')
            if fields and fields[0] == 'Starttime of export':
                start_filetime = int(fields[1])
            elif fields and fields[0] == 'SampleTime[ms]':
                sample_period_ms = float(fields[1])
            if re.match(r'^\d+\t\d+\s*$', line):
                data_line = line_number
                break
    if start_filetime is None or sample_period_ms is None or data_line is None:
        raise RuntimeError('IDS export metadata or numeric data were not found')
    numeric = np.loadtxt(
        path, delimiter='\t', skiprows=data_line, usecols=(0, 1),
        dtype=np.uint64, comments='EOF',
    )
    raw = numeric[:, 1].astype(np.uint32)
    delta = np.diff(raw.astype(np.int64))
    delta[delta > 2**31] -= 2**32
    delta[delta < -(2**31)] += 2**32
    position_nm = np.empty(raw.size, dtype=np.float64)
    position_nm[0] = 0.0
    position_nm[1:] = np.cumsum(delta)
    time_s = np.arange(raw.size, dtype=np.float64) * sample_period_ms * 1.0e-3
    start_epoch_s = (start_filetime - FILETIME_UNIX_EPOCH) / 1.0e7
    return time_s, position_nm, start_epoch_s, sample_period_ms * 1.0e-3


def load_log_rows(path: Path, ids_start_epoch_s: float):
    rows = []
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            instant = datetime.fromisoformat(row['utc'])
            row['ids_time_s'] = instant.timestamp() - ids_start_epoch_s
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


def sample_bounds(start_s, end_s, sample_period_s, sample_count):
    start = max(0, int(np.ceil(start_s / sample_period_s)))
    end = min(sample_count, int(np.floor(end_s / sample_period_s)) + 1)
    return start, end


def measured_window(time_s, position_nm, sample_period_s, start_s, end_s):
    first, last = sample_bounds(start_s, end_s, sample_period_s, position_nm.size)
    baseline = float(np.median(position_nm[first:first + 20]))
    t = time_s[first:last] - start_s
    y = (position_nm[first:last] - baseline) / 1000.0
    return t, y


def command_series(rows, run_index, block_start_s, block_end_s):
    moves = [
        r for r in rows if r['event'] == 'MOVE_ACK'
        and r['run_index'] == str(run_index) and r['ideal_position_rev']
        and block_start_s <= r['ids_time_s'] <= block_end_s
    ]
    times = [0.0]
    positions_um = [0.0]
    for row in moves:
        times.append(row['ids_time_s'] - block_start_s)
        positions_um.append(float(row['ideal_position_rev']) * LEAD_PITCH_M * 1.0e6)
    times.append(block_end_s - block_start_s)
    positions_um.append(positions_um[-1])
    return np.asarray(times), np.asarray(positions_um)


def style_axis(ax):
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Time since block start (s)', fontsize=8)
    ax.set_ylabel('Position (µm)', fontsize=8)
    ax.tick_params(labelsize=7)


def build_montage(run_index, mres, current, rows, time_s, position_nm, sample_period_s):
    fig, axes = plt.subplots(2, 5, figsize=(22.0, 8.0))
    panels = list(axes.flat)

    # 8 D-plateau panels, each with its own independent y-scale.
    for ax, rate in zip(panels[:8], D_RATES):
        block_name = f'D_{rate}'
        start_s, end_s = find_block(rows, run_index, block_name)
        t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
        cmd_t, cmd_y = command_series(rows, run_index, start_s, end_s)
        ax.plot(t, y, color=MEASURED_COLOR, lw=0.8, label='Measured')
        ax.step(cmd_t, cmd_y, where='post', color=COMMAND_COLOR, lw=0.9,
                alpha=0.85, label='Commanded')
        ax.set_title(f'D {rate} full-steps/s', fontsize=9)
        style_axis(ax)

    # 1 panel for C.
    ax_c = panels[8]
    start_s, end_s = find_block(rows, run_index, 'C')
    t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
    cmd_t, cmd_y = command_series(rows, run_index, start_s, end_s)
    ax_c.plot(t, y, color=MEASURED_COLOR, lw=0.9, label='Measured')
    ax_c.step(cmd_t, cmd_y, where='post', color=COMMAND_COLOR, lw=1.0,
              alpha=0.85, label='Commanded')
    ax_c.set_title('C (creep/settling)', fontsize=9)
    style_axis(ax_c)

    # 1 panel: BLOCK_0_START and BLOCK_0_END overlaid.
    ax_ref = panels[9]
    for block_name, label, color in (
            ('BLOCK_0_START', 'Block start', MEASURED_COLOR),
            ('BLOCK_0_END', 'Block end', COMMAND_COLOR)):
        start_s, end_s = find_block(rows, run_index, block_name)
        t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
        ax_ref.plot(t, y, color=color, lw=0.9, label=label)
    ax_ref.set_title('BLOCK_0: start vs. end', fontsize=9)
    style_axis(ax_ref)

    for ax in (panels[0], ax_c, ax_ref):
        ax.legend(loc='best', fontsize=7, framealpha=0.9)

    fig.suptitle(
        f"Run {run_index} — MRES 1/{mres}, "
        f"{CURRENT_SHORT.get(current, current)}",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    folder_name = f'run_{run_index:02d}_mres_{mres}_{current.lower()}'
    out_dir = SUBPLOT_DIR / folder_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'block_montage.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main() -> None:
    time_s, position_nm, start_epoch_s, sample_period_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)

    for run_index in range(1, 7):
        start_row = next(
            r for r in rows if r['event'] == 'RUN_CONFIG'
            and r['run_index'] == str(run_index)
        )
        mres, current = start_row['mres'], start_row['current']
        out_path = build_montage(
            run_index, mres, current, rows, time_s, position_nm,
            sample_period_s,
        )
        print(f'Run {run_index}: {out_path}')


if __name__ == '__main__':
    main()

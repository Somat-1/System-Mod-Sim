#!/usr/bin/env python3
"""Plot the complete raw v3 IDS encoder record, with each of the six
run configurations shaded and labeled directly on the plot (no legend)."""

from __future__ import annotations

import csv
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / 'data' / 'raw_local'
ASSET_DIR = ROOT / 'rendered_assets'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
OUT_PATH = ASSET_DIR / 'full_raw_sequence.png'

FILETIME_UNIX_EPOCH = 116444736000000000


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
    # Encoder is nm-native (1 count = 1 nm). Relative to the first sample
    # (which is defined as zero), correcting for 32-bit counter wraparound
    # via signed deltas rather than a naive subtraction.
    delta = np.diff(raw.astype(np.int64))
    delta[delta > 2**31] -= 2**32
    delta[delta < -(2**31)] += 2**32
    position_nm = np.empty(raw.size, dtype=np.float64)
    position_nm[0] = 0.0
    position_nm[1:] = np.cumsum(delta)
    time_s = np.arange(raw.size, dtype=np.float64) * sample_period_ms * 1.0e-3
    start_epoch_s = (start_filetime - FILETIME_UNIX_EPOCH) / 1.0e7
    return time_s, position_nm, start_epoch_s


def parse_runs(path: Path, ids_start_epoch_s: float):
    rows = []
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            instant = datetime.fromisoformat(row['utc'])
            row['ids_time_s'] = instant.timestamp() - ids_start_epoch_s
            rows.append(row)
    runs = []
    for run_index in range(1, 7):
        start_row = next(
            r for r in rows
            if r['event'] == 'RUN_CONFIG' and r['run_index'] == str(run_index)
        )
        end_row = next(
            r for r in rows
            if r['event'] == 'RUN_COMPLETE' and r['run_index'] == str(run_index)
        )
        runs.append({
            'run_index': run_index,
            'mres': start_row['mres'],
            'current': start_row['current'],
            'start_s': start_row['ids_time_s'],
            'end_s': end_row['ids_time_s'],
        })
    return runs


CURRENT_SHORT = {'I_50pct': '50% I', 'I_100pct': '100% I'}


def main() -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    time_s, position_nm, start_epoch_s = parse_ids(IDS_PATH)
    runs = parse_runs(LOG_PATH, start_epoch_s)

    bin_samples = 50
    usable = position_nm.size - position_nm.size % bin_samples
    t_min = time_s[:usable].reshape(-1, bin_samples).mean(axis=1) / 60.0
    y_um = (
        position_nm[:usable].reshape(-1, bin_samples).mean(axis=1) / 1000.0
    )

    fig, ax = plt.subplots(figsize=(16.0, 6.0), constrained_layout=True)
    ax.plot(t_min, y_um, color='#136f63', lw=0.7)

    colors = ('#e8f1f2', '#f5e6cc')
    for run in runs:
        left = run['start_s'] / 60.0
        right = run['end_s'] / 60.0
        ax.axvspan(
            left, right, color=colors[(run['run_index'] - 1) % 2],
            alpha=0.6, zorder=-1,
        )
        ax.text(
            (left + right) / 2.0, 0.97,
            f"R{run['run_index']}  1/{run['mres']}, "
            f"{CURRENT_SHORT.get(run['current'], run['current'])}",
            transform=ax.get_xaxis_transform(), ha='center', va='top',
            fontsize=9,
        )

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Position (µm)')
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(t_min[0], t_min[-1])

    fig.savefig(OUT_PATH, dpi=160)
    plt.close(fig)
    print(f'Samples: {position_nm.size:,}')
    print(f'Saved: {OUT_PATH}')


if __name__ == '__main__':
    main()

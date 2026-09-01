#!/usr/bin/env python3
"""Plot the complete raw v3 IDS encoder record, with each of the six
run configurations shaded/labeled directly on the plot, the
run-transition marker moves highlighted, and per-run subplots + data
written to rendered_assets/individual_subplots/.

Each of the ten analysis-panel blocks (matching plot_block_montage.py)
is additionally shaded with its own opaque color on both the overview
and the per-run plots, so the raw trace shows exactly which stretch of
data feeds which montage panel. A shared color key applies across runs.

The run-transition marker is the "CONFIG_0N_..." data-visible separator
signature commanded by
../../v2/scripts/run_identification_dedicated_controller.py
(64 + 4*index full steps, negative-then-positive leap at 150 full
steps/s, 1.0 s reverse dwell, 0.5 s settle) -- it is the physical
"indicative move" that marks where one run ends and the next begins.
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / 'data' / 'raw_local'
ASSET_DIR = ROOT / 'rendered_assets'
SUBPLOT_DIR = ASSET_DIR / 'individual_subplots'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
OUT_PATH = ASSET_DIR / 'full_raw_sequence.png'

FILETIME_UNIX_EPOCH = 116444736000000000
CURRENT_SHORT = {'I_50pct': '50% I', 'I_100pct': '100% I'}
HIGHLIGHT_ALPHA = 0.15
# Conditioning uses the given yellow at the same low alpha as the rest, it
# is nearly invisible -- boost it specifically.
TYPE_ALPHA = {
    'reference': HIGHLIGHT_ALPHA,
    'creep': HIGHLIGHT_ALPHA,
    'conditioning': 0.6,
    'plateau': HIGHLIGHT_ALPHA,
}

# Four experiment types, one color each (reused across every block of that
# type -- e.g. all eight D plateaus share one color). Creep sits between
# the reference grey and the conditioning yellow so it doesn't read as the
# same color as reference at low alpha.
TYPE_COLORS = {
    'reference': '#454040',      # BLOCK_0_START / BLOCK_0_END
    'creep': '#9C975B',          # C -- olive, between grey and yellow
    'conditioning': '#D8D365',   # COND_BEFORE_C / COND_BEFORE_D
    'plateau': '#E6F082',        # D_0.125 ... D_200
}

# Chronological order of blocks shaded on the overview/per-run plots.
# BLOCK_0 stands in for both BLOCK_0_START and BLOCK_0_END.
BLOCK_ORDER = (
    'BLOCK_0', 'COND_BEFORE_C', 'C', 'COND_BEFORE_D', 'D_0.125', 'D_0.375',
    'D_1.25', 'D_3.5', 'D_9.5', 'D_27.5', 'D_70', 'D_200',
)
BLOCK_TYPE = {
    'BLOCK_0': 'reference',
    'COND_BEFORE_C': 'conditioning',
    'C': 'creep',
    'COND_BEFORE_D': 'conditioning',
    **{f'D_{rate}': 'plateau'
       for rate in ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')},
}
# Legend shows one swatch per type, not per block.
TYPE_LABELS = (
    ('reference', 'Reference (BLOCK_0)'),
    ('creep', 'Creep/settling (C)'),
    ('conditioning', 'Conditioning'),
    ('plateau', 'Velocity plateaus (D)'),
)


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
    # (defined as zero), correcting for 32-bit counter wraparound via
    # signed deltas rather than a naive subtraction.
    delta = np.diff(raw.astype(np.int64))
    delta[delta > 2**31] -= 2**32
    delta[delta < -(2**31)] += 2**32
    position_nm = np.empty(raw.size, dtype=np.float64)
    position_nm[0] = 0.0
    position_nm[1:] = np.cumsum(delta)
    time_s = np.arange(raw.size, dtype=np.float64) * sample_period_ms * 1.0e-3
    start_epoch_s = (start_filetime - FILETIME_UNIX_EPOCH) / 1.0e7
    return time_s, position_nm, start_epoch_s, sample_period_ms


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


def block_spans(rows, run_index):
    """(label, start_s, end_s) for the analysis blocks, in chronological order."""
    spans = []
    for label in BLOCK_ORDER:
        if label == 'BLOCK_0':
            s0, e0 = find_block(rows, run_index, 'BLOCK_0_START')
            s1, e1 = find_block(rows, run_index, 'BLOCK_0_END')
            spans.append(('BLOCK_0_START', s0, e0))
            spans.append(('BLOCK_0_END', s1, e1))
        else:
            s, e = find_block(rows, run_index, label)
            spans.append((label, s, e))
    return spans


def parse_runs(rows):
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
        marker_start = next(
            r for r in rows
            if r['event'] == 'BLOCK_START' and r['run_index'] == str(run_index)
            and r['block'].startswith('MARKER_CONFIG')
        )
        marker_end = next(
            r for r in rows
            if r['event'] == 'BLOCK_END' and r['run_index'] == str(run_index)
            and r['block'].startswith('MARKER_CONFIG')
        )
        runs.append({
            'run_index': run_index,
            'mres': start_row['mres'],
            'current': start_row['current'],
            'start_s': start_row['ids_time_s'],
            'end_s': end_row['ids_time_s'],
            'marker_start_s': marker_start['ids_time_s'],
            'marker_end_s': marker_end['ids_time_s'],
            'block_spans': block_spans(rows, run_index),
        })
    return runs


def sample_bounds(start_s, end_s, sample_period_s, sample_count):
    start = max(0, int(np.ceil(start_s / sample_period_s)))
    end = min(sample_count, int(np.floor(end_s / sample_period_s)) + 1)
    return start, end


def block_type(label):
    key = 'BLOCK_0' if label.startswith('BLOCK_0') else label
    return BLOCK_TYPE[key]


def block_color(label):
    return TYPE_COLORS[block_type(label)]


def add_block_key(ax):
    handles = [
        Patch(color=TYPE_COLORS[key], alpha=0.6, label=label)
        for key, label in TYPE_LABELS
    ]
    ax.legend(
        handles=handles, loc='upper center', ncol=4, fontsize=8,
        frameon=False, bbox_to_anchor=(0.5, -0.14),
    )


def plot_overview(t_s, position_nm, runs, sample_period_ms):
    bin_samples = 50
    usable = position_nm.size - position_nm.size % bin_samples
    t_min = t_s[:usable].reshape(-1, bin_samples).mean(axis=1) / 60.0
    y_um = position_nm[:usable].reshape(-1, bin_samples).mean(axis=1) / 1000.0

    fig, ax = plt.subplots(figsize=(16.0, 6.5), constrained_layout=True)
    ax.plot(t_min, y_um, color='#136f63', lw=0.7, zorder=2)

    zebra = ('#f4f4f4', '#ffffff')
    for run in runs:
        left = run['start_s'] / 60.0
        right = run['end_s'] / 60.0
        ax.axvspan(left, right, color=zebra[(run['run_index'] - 1) % 2],
                   alpha=0.5, zorder=-3)
        ax.text(
            (left + right) / 2.0, 1.02,
            f"R{run['run_index']}  1/{run['mres']}, "
            f"{CURRENT_SHORT.get(run['current'], run['current'])}",
            transform=ax.get_xaxis_transform(), ha='center', va='bottom',
            fontsize=9,
        )
        for label, start_s, end_s in run['block_spans']:
            ax.axvspan(
                start_s / 60.0, end_s / 60.0, color=block_color(label),
                alpha=TYPE_ALPHA[block_type(label)], zorder=-2, linewidth=0,
            )
        marker_left = run['marker_start_s'] / 60.0
        marker_right = run['marker_end_s'] / 60.0
        ax.axvspan(
            marker_left, marker_right, color='#c0392b',
            alpha=HIGHLIGHT_ALPHA, zorder=-1,
        )
        ax.axvline(
            marker_left, color='#c0392b', lw=1.1, linestyle='--', zorder=3,
        )

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Position (µm)')
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(t_min[0], t_min[-1])
    add_block_key(ax)

    fig.savefig(OUT_PATH, dpi=160)
    plt.close(fig)


def plot_and_save_run(run, t_s, position_nm, sample_period_ms):
    sample_period_s = sample_period_ms * 1.0e-3
    first, last = sample_bounds(
        run['start_s'], run['end_s'], sample_period_s, position_nm.size,
    )
    local_t_s = t_s[first:last] - run['start_s']
    local_y_um = (position_nm[first:last] - position_nm[first]) / 1000.0

    folder_name = (
        f"run_{run['run_index']:02d}_mres_{run['mres']}_"
        f"{run['current'].lower()}"
    )
    run_dir = SUBPLOT_DIR / folder_name
    run_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12.0, 5.5), constrained_layout=True)
    ax.plot(local_t_s, local_y_um, color='#136f63', lw=0.6, zorder=2)

    for label, start_s, end_s in run['block_spans']:
        ax.axvspan(
            start_s - run['start_s'], end_s - run['start_s'],
            color=block_color(label), alpha=TYPE_ALPHA[block_type(label)],
            zorder=-2, linewidth=0,
        )
    marker_span = (
        run['marker_start_s'] - run['start_s'],
        run['marker_end_s'] - run['start_s'],
    )
    ax.axvspan(*marker_span, color='#c0392b', alpha=HIGHLIGHT_ALPHA, zorder=-1)

    ax.set_xlabel('Time since run start (s)')
    ax.set_ylabel('Position (µm)')
    ax.set_title(
        f"Run {run['run_index']} — MRES 1/{run['mres']}, "
        f"{CURRENT_SHORT.get(run['current'], run['current'])}"
    )
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(local_t_s[0], local_t_s[-1])
    add_block_key(ax)
    fig.savefig(run_dir / 'sequence_plot.png', dpi=160)
    plt.close(fig)

    np.savez_compressed(
        run_dir / 'sequence_data.npz',
        time_s=local_t_s,
        position_um=local_y_um,
        run_index=np.asarray(run['run_index']),
        mres=np.asarray(run['mres']),
        current=np.asarray(run['current']),
        marker_span_s=np.asarray(marker_span),
    )


def main() -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    t_s, position_nm, start_epoch_s, sample_period_ms = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    runs = parse_runs(rows)

    plot_overview(t_s, position_nm, runs, sample_period_ms)
    for run in runs:
        plot_and_save_run(run, t_s, position_nm, sample_period_ms)

    print(f'Samples: {position_nm.size:,}')
    print(f'Overview: {OUT_PATH}')
    print(f'Per-run subfolders: {SUBPLOT_DIR}')


if __name__ == '__main__':
    main()

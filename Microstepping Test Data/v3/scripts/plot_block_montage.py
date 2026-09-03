#!/usr/bin/env python3
"""Per-configuration 10-panel block montage, in chronological order:
BLOCK_0 (start/end overlaid), C, then the 8 D plateaus ascending. Each
D panel keeps its own independent y-scale. Panel title colors match the
block-color key used in plot_full_raw_sequence.py's overview/per-run
highlight overlays, so a panel can be matched back to its source window
on the raw trace.

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

from command_reconstruction import reconstruct_segments, sample_segments

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / 'data' / 'raw_local'
SUBPLOT_DIR = ROOT / 'rendered_assets' / 'individual_subplots'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'

FILETIME_UNIX_EPOCH = 116444736000000000
# See plot_full_raw_sequence.py / README.md "Controller/IDS clock skew".
CONTROLLER_CLOCK_SKEW_S = 0.319
LEAD_PITCH_M = 2.0e-3  # stage lead screw pitch: 2 mm/rev
CURRENT_SHORT = {'I_50pct': '50% I', 'I_100pct': '100% I'}
D_RATES = ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')
SLOW_PLATEAU_D_RATES = ('0.125', '0.375', '1.25')
MEASURED_COLOR = '#136f63'
COMMAND_COLOR = '#9a9a9a'
ANNOTATION_BBOX = dict(facecolor='white', edgecolor='none', alpha=0.78, pad=1.5)

# D-rate/MRES combinations confirmed (README.md, "D_3.5 and D_9.5 panels
# look wrong at coarser MRES") to complete far faster than commanded --
# an execution-timing artifact in the controller, not a plotting bug.
# Excluded from analysis; marked with a cross rather than the C-panel's
# star.
ANOMALOUS_D_RATES_BY_MRES = {
    '4': (),
    '2': ('3.5', '9.5'),
    '1': ('3.5', '9.5', '27.5'),
}

# Same four experiment-type colors as plot_full_raw_sequence.py's key.
TYPE_COLORS = {
    'reference': '#454040',
    'creep': '#9C975B',
    'plateau': '#E6F082',
}
AXIS_COLOR = '#333333'


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


def command_series(rows, run_index, block_start_s, block_end_s, truncate=False,
                    block_name=None):
    """truncate=True stops the series right after the last real MOVE_ACK
    instead of extending it flat to the block end -- use where the stage
    is known not to hold the commanded value for the remainder of the
    block (see the C-block investigation in README.md).

    For controller-paced D-rate moves (3.5 full-steps/s and above), this
    traces the real ACCEL/CONST/DECEL ramp the controller executed instead
    of a single zero-order-hold step to the final target -- see
    command_reconstruction.py and README.md, "Controller-paced D-rate
    command reconstruction". block_name enables the same reconstruction
    for fixed-rate burst moves (e.g. the C block's approach/return) via
    BURST_RATE_BY_BLOCK_FS_S."""
    segments = reconstruct_segments(
        rows, run_index, block_start_s, block_end_s, block_name
    )
    return sample_segments(
        segments, scale=LEAD_PITCH_M * 1.0e6, truncate=truncate
    )


def style_axis(ax):
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    ax.set_xlabel('Time since block start (s)', fontsize=8)
    ax.set_ylabel('Position (µm)', fontsize=8)
    ax.tick_params(labelsize=7, colors=AXIS_COLOR)


def build_montage(run_index, mres, current, rows, time_s, position_nm, sample_period_s):
    fig, axes = plt.subplots(2, 5, figsize=(22.0, 8.0))
    panels = list(axes.flat)

    # Chronological order: BLOCK_0, C, then the 8 D plateaus ascending.
    ax_ref = panels[0]
    for block_name, label, color in (
            ('BLOCK_0_START', 'Block start', MEASURED_COLOR),
            ('BLOCK_0_END', 'Block end', COMMAND_COLOR)):
        start_s, end_s = find_block(rows, run_index, block_name)
        t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
        ax_ref.plot(t, y, color=color, lw=0.9, label=label)
    ax_ref.set_title('BLOCK_0: start vs. end', fontsize=9,
                      color=TYPE_COLORS['reference'])
    style_axis(ax_ref)
    ax_ref.legend(loc='best', fontsize=7, framealpha=0.9)

    ax_c = panels[1]
    start_s, end_s = find_block(rows, run_index, 'C')
    t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
    cmd_t, cmd_y = command_series(rows, run_index, start_s, end_s, block_name='C')
    ax_c.step(cmd_t, cmd_y, where='post', color=COMMAND_COLOR, lw=0.8,
              linestyle='--', alpha=0.6, zorder=1,
              label='Commanded (not sustained; see note)')
    ax_c.plot(t, y, color=MEASURED_COLOR, lw=1.1, zorder=3, label='Measured')
    ax_c.plot(
        0.06, 0.92, marker='*', markersize=16, color='#b30000',
        markeredgecolor='black', markeredgewidth=0.6, transform=ax_c.transAxes,
        zorder=5, clip_on=False,
    )
    ax_c.text(
        0.10, 0.92, 'did not track command', transform=ax_c.transAxes,
        fontsize=6.5, color='#b30000', ha='left', va='center', zorder=5,
        bbox=ANNOTATION_BBOX,
    )
    ax_c.set_title('C (creep/settling)', fontsize=9, color=TYPE_COLORS['creep'])
    style_axis(ax_c)
    ax_c.legend(loc='center right', fontsize=6.5, framealpha=0.9)

    anomalous_rates = ANOMALOUS_D_RATES_BY_MRES.get(str(mres), ())
    any_anomalous = False
    any_reconstructed = False
    for ax, rate in zip(panels[2:], D_RATES):
        label = f'D_{rate}'
        block_name = label
        start_s, end_s = find_block(rows, run_index, block_name)
        t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
        cmd_t, cmd_y = command_series(rows, run_index, start_s, end_s, block_name=block_name)
        cmd_label = (
            'Commanded (ramp, reconstructed)'
            if rate not in SLOW_PLATEAU_D_RATES else 'Commanded'
        )
        ax.step(cmd_t, cmd_y, where='post', color=COMMAND_COLOR, lw=0.7,
                linestyle='--', alpha=0.6, zorder=1, label=cmd_label)
        ax.plot(t, y, color=MEASURED_COLOR, lw=1.0, zorder=3, label='Measured')
        ax.set_title(f'D {rate} full-steps/s', fontsize=9, color=AXIS_COLOR)
        style_axis(ax)
        if rate not in SLOW_PLATEAU_D_RATES:
            any_reconstructed = True
            ax.text(
                0.97, 0.06, '† ramp reconstructed', transform=ax.transAxes,
                fontsize=6.0, color=AXIS_COLOR, ha='right', va='bottom',
                style='italic', zorder=5, bbox=ANNOTATION_BBOX,
            )
        if rate in anomalous_rates:
            any_anomalous = True
            ax.plot(
                0.06, 0.90, marker='x', markersize=13, markeredgewidth=3,
                color='#b30000', transform=ax.transAxes, zorder=5,
                clip_on=False,
            )
            ax.text(
                0.13, 0.90, 'excluded', transform=ax.transAxes,
                fontsize=6.5, color='#b30000', ha='left', va='center',
                zorder=5, bbox=ANNOTATION_BBOX,
            )

    panels[2].legend(loc='best', fontsize=7, framealpha=0.9)

    fig.suptitle(
        f"Run {run_index} — MRES 1/{mres}, "
        f"{CURRENT_SHORT.get(current, current)} "
        "(panels in chronological order: BLOCK_0 → C → D_0.125 → ... → D_200)",
        fontsize=12.5,
    )
    footnote_lines = [
        "★ C panel: the stage did not track the commanded ±40 µm "
        "approach/return as intended -- see the C-block investigation in "
        "README.md",
    ]
    if any_anomalous:
        footnote_lines.append(
            "✗ marked D panels completed far faster than commanded "
            "(controller execution artifact, not a plotting issue) and are "
            "excluded from analysis -- see README.md"
        )
    if any_reconstructed:
        footnote_lines.append(
            "† D_3.5 and above: commanded line is a reconstructed ACCEL/"
            "CONST/DECEL ramp (the controller executes these as one command "
            "and only the start/end are logged) -- see README.md, "
            "\"Controller-paced D-rate command reconstruction\""
        )
    fig.text(
        0.5, 0.012, "\n".join(footnote_lines),
        ha='center', va='bottom', fontsize=8.5, color='#b30000',
    )
    bottom_margin = 0.03 + 0.03 * any_anomalous + 0.03 * any_reconstructed
    fig.tight_layout(rect=(0, bottom_margin, 1, 0.94))

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

#!/usr/bin/env python3
"""Tracking-error companion to plot_block_montage_with_simulation.py: same
10-panel chronological layout (BLOCK_0, C, then the 8 D plateaus
ascending), but each panel plots (trace - reconstructed commanded) instead
of raw position -- how far the measured IDS data and the simulated model
each deviate from what was actually commanded (see command_reconstruction.py
for how "commanded" is reconstructed, and README.md, "Controller-paced
D-rate command reconstruction" and "Controller/IDS clock skew").

A second set of 10 panels below repeats the same chronological layout for
(measured - simulated) directly -- model-vs-hardware fidelity, independent
of the commanded trajectory -- with the simulated trace linearly
interpolated onto the measured trace's own (denser, uniform) time grid.

Carries the same C-panel/anomalous-D-rate annotations as the position
montages, since a block that didn't track its command, or whose command
reconstruction is itself compressed by the execution-timing anomaly, needs
the same caveat here.
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

from command_reconstruction import commanded_value_at, reconstruct_segments

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / 'data' / 'raw_local'
SUBPLOT_DIR = ROOT / 'rendered_assets' / 'individual_subplots'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'

FILETIME_UNIX_EPOCH = 116444736000000000
# See plot_full_raw_sequence.py / README.md "Controller/IDS clock skew".
CONTROLLER_CLOCK_SKEW_S = 0.319
LEAD_PITCH_M = 2.0e-3
CURRENT_SHORT = {'I_50pct': '50% I', 'I_100pct': '100% I'}
D_RATES = ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')
SLOW_PLATEAU_D_RATES = ('0.125', '0.375', '1.25')
MEASURED_COLOR = '#136f63'
SIMULATED_COLOR = '#e67e22'
MISMATCH_COLOR = '#7570b3'
ZERO_COLOR = '#9a9a9a'
MODEL_LABEL = 'Rev 4.2 parallel LuGre + nonlinear detent torque'
ANNOTATION_BBOX = dict(facecolor='white', edgecolor='none', alpha=0.78, pad=1.5)

ANOMALOUS_D_RATES_BY_MRES = {
    '4': (),
    '2': ('3.5', '9.5'),
    '1': ('3.5', '9.5', '27.5'),
}

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


def sim_trace(sim_npz, block_name):
    t_key, y_key = f'{block_name}_time_s', f'{block_name}_position_um'
    if sim_npz is None or t_key not in sim_npz:
        return None, None
    return sim_npz[t_key], sim_npz[y_key]


def style_axis(ax, ylabel='Tracking error (µm)'):
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    ax.set_xlabel('Time since block start (s)', fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7, colors=AXIS_COLOR)


def error_traces(rows, run_index, start_s, end_s, time_s, position_nm,
                  sample_period_s, sim_npz, block_name):
    segments = reconstruct_segments(rows, run_index, start_s, end_s, block_name)
    t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
    cmd_at_t = commanded_value_at(segments, t, scale=LEAD_PITCH_M * 1.0e6)
    err_measured = (t, y - cmd_at_t)

    sim_t, sim_y = sim_trace(sim_npz, block_name)
    err_simulated = None
    if sim_t is not None:
        cmd_at_sim_t = commanded_value_at(segments, sim_t, scale=LEAD_PITCH_M * 1.0e6)
        err_simulated = (sim_t, sim_y - cmd_at_sim_t)
    return err_measured, err_simulated


def measured_minus_sim_trace(time_s, position_nm, sample_period_s, start_s,
                              end_s, sim_npz, block_name):
    """(t, measured - simulated) on measured's own time grid, with the
    (sparser, non-uniform) simulated trace linearly interpolated onto it.
    None if no simulated trace is available for this block."""
    t, y = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
    sim_t, sim_y = sim_trace(sim_npz, block_name)
    if sim_t is None:
        return None
    sim_on_measured_grid = np.interp(t, sim_t, sim_y)
    return t, y - sim_on_measured_grid


def annotate_excluded(ax, bbox=ANNOTATION_BBOX):
    ax.plot(
        0.06, 0.90, marker='x', markersize=13, markeredgewidth=3,
        color='#b30000', transform=ax.transAxes, zorder=5, clip_on=False,
    )
    ax.text(
        0.13, 0.90, 'excluded', transform=ax.transAxes,
        fontsize=6.5, color='#b30000', ha='left', va='center',
        zorder=5, bbox=bbox,
    )


def annotate_reconstructed(ax, bbox=ANNOTATION_BBOX):
    ax.text(
        0.97, 0.06, '† ramp reconstructed', transform=ax.transAxes,
        fontsize=6.0, color=AXIS_COLOR, ha='right', va='bottom',
        style='italic', zorder=5, bbox=bbox,
    )


def annotate_did_not_track(ax, bbox=ANNOTATION_BBOX):
    ax.plot(
        0.06, 0.92, marker='*', markersize=16, color='#b30000',
        markeredgecolor='black', markeredgewidth=0.6, transform=ax.transAxes,
        zorder=5, clip_on=False,
    )
    ax.text(
        0.10, 0.92, 'did not track command', transform=ax.transAxes,
        fontsize=6.5, color='#b30000', ha='left', va='center', zorder=5,
        bbox=bbox,
    )


def build_vs_commanded_row(panels, run_index, mres, rows, time_s, position_nm,
                            sample_period_s, sim_npz):
    ax_ref, ax_c = panels[0], panels[1]

    for block_name, label, sim_style in (
            ('BLOCK_0_START', 'start', '-'), ('BLOCK_0_END', 'end', '--')):
        start_s, end_s = find_block(rows, run_index, block_name)
        err_measured, err_simulated = error_traces(
            rows, run_index, start_s, end_s, time_s, position_nm,
            sample_period_s, sim_npz, block_name,
        )
        ax_ref.axhline(0.0, color=ZERO_COLOR, lw=0.7, linestyle='--', zorder=1)
        if err_simulated is not None:
            ax_ref.plot(*err_simulated, color=SIMULATED_COLOR, lw=1.0,
                        linestyle=sim_style, alpha=0.9, zorder=2,
                        label=f'Sim error ({label})')
        ax_ref.plot(*err_measured, color=MEASURED_COLOR, lw=1.0,
                    linestyle=sim_style, zorder=3, label=f'Measured error ({label})')
    ax_ref.set_title('BLOCK_0: start vs. end', fontsize=9,
                      color=TYPE_COLORS['reference'])
    style_axis(ax_ref)
    ax_ref.legend(loc='best', fontsize=6.0, framealpha=0.9)

    start_s, end_s = find_block(rows, run_index, 'C')
    err_measured, err_simulated = error_traces(
        rows, run_index, start_s, end_s, time_s, position_nm,
        sample_period_s, sim_npz, 'C',
    )
    ax_c.axhline(0.0, color=ZERO_COLOR, lw=0.7, linestyle='--', zorder=1)
    if err_simulated is not None:
        ax_c.plot(*err_simulated, color=SIMULATED_COLOR, lw=1.1, alpha=0.9,
                  zorder=2, label='Simulated error')
    ax_c.plot(*err_measured, color=MEASURED_COLOR, lw=1.1, zorder=3,
              label='Measured error')
    annotate_did_not_track(ax_c)
    ax_c.set_title('C (creep/settling)', fontsize=9, color=TYPE_COLORS['creep'])
    style_axis(ax_c)
    ax_c.legend(loc='center right', fontsize=6.5, framealpha=0.9)

    anomalous_rates = ANOMALOUS_D_RATES_BY_MRES.get(str(mres), ())
    any_anomalous = False
    any_reconstructed = False
    for panel_index, (ax, rate) in enumerate(zip(panels[2:], D_RATES)):
        block_name = f'D_{rate}'
        start_s, end_s = find_block(rows, run_index, block_name)
        err_measured, err_simulated = error_traces(
            rows, run_index, start_s, end_s, time_s, position_nm,
            sample_period_s, sim_npz, block_name,
        )
        ax.axhline(0.0, color=ZERO_COLOR, lw=0.7, linestyle='--', zorder=1)
        if err_simulated is not None:
            ax.plot(*err_simulated, color=SIMULATED_COLOR, lw=0.9, alpha=0.9,
                    zorder=2, label='Simulated error')
        ax.plot(*err_measured, color=MEASURED_COLOR, lw=1.0, zorder=3,
                label='Measured error')
        ax.set_title(f'D {rate} full-steps/s', fontsize=9, color=AXIS_COLOR)
        style_axis(ax)
        if panel_index == 0:
            ax.legend(loc='best', fontsize=7, framealpha=0.9)
        if rate not in SLOW_PLATEAU_D_RATES:
            any_reconstructed = True
            annotate_reconstructed(ax)
        if rate in anomalous_rates:
            any_anomalous = True
            annotate_excluded(ax)

    return any_anomalous, any_reconstructed


def build_measured_minus_sim_row(panels, run_index, mres, rows, time_s,
                                  position_nm, sample_period_s, sim_npz):
    ax_ref, ax_c = panels[0], panels[1]
    ylabel = 'Measured − simulated (µm)'

    for block_name, label, sim_style in (
            ('BLOCK_0_START', 'start', '-'), ('BLOCK_0_END', 'end', '--')):
        start_s, end_s = find_block(rows, run_index, block_name)
        diff = measured_minus_sim_trace(
            time_s, position_nm, sample_period_s, start_s, end_s,
            sim_npz, block_name,
        )
        ax_ref.axhline(0.0, color=ZERO_COLOR, lw=0.7, linestyle='--', zorder=1)
        if diff is not None:
            ax_ref.plot(*diff, color=MISMATCH_COLOR, lw=1.0,
                        linestyle=sim_style, zorder=3, label=f'Measured − sim ({label})')
    ax_ref.set_title('BLOCK_0: start vs. end', fontsize=9,
                      color=TYPE_COLORS['reference'])
    style_axis(ax_ref, ylabel)
    ax_ref.legend(loc='best', fontsize=6.0, framealpha=0.9)

    start_s, end_s = find_block(rows, run_index, 'C')
    diff = measured_minus_sim_trace(
        time_s, position_nm, sample_period_s, start_s, end_s, sim_npz, 'C',
    )
    ax_c.axhline(0.0, color=ZERO_COLOR, lw=0.7, linestyle='--', zorder=1)
    if diff is not None:
        ax_c.plot(*diff, color=MISMATCH_COLOR, lw=1.1, zorder=3,
                  label='Measured − sim')
    annotate_did_not_track(ax_c)
    ax_c.set_title('C (creep/settling)', fontsize=9, color=TYPE_COLORS['creep'])
    style_axis(ax_c, ylabel)
    ax_c.legend(loc='center right', fontsize=6.5, framealpha=0.9)

    anomalous_rates = ANOMALOUS_D_RATES_BY_MRES.get(str(mres), ())
    any_anomalous = False
    any_reconstructed = False
    for panel_index, (ax, rate) in enumerate(zip(panels[2:], D_RATES)):
        block_name = f'D_{rate}'
        start_s, end_s = find_block(rows, run_index, block_name)
        diff = measured_minus_sim_trace(
            time_s, position_nm, sample_period_s, start_s, end_s,
            sim_npz, block_name,
        )
        ax.axhline(0.0, color=ZERO_COLOR, lw=0.7, linestyle='--', zorder=1)
        if diff is not None:
            ax.plot(*diff, color=MISMATCH_COLOR, lw=1.0, zorder=3,
                    label='Measured − sim')
        ax.set_title(f'D {rate} full-steps/s', fontsize=9, color=AXIS_COLOR)
        style_axis(ax, ylabel)
        if panel_index == 0:
            ax.legend(loc='best', fontsize=7, framealpha=0.9)
        if rate not in SLOW_PLATEAU_D_RATES:
            any_reconstructed = True
            annotate_reconstructed(ax)
        if rate in anomalous_rates:
            any_anomalous = True
            annotate_excluded(ax)

    return any_anomalous, any_reconstructed


def build_montage(run_index, mres, current, rows, time_s, position_nm,
                   sample_period_s, sim_npz):
    fig, axes = plt.subplots(4, 5, figsize=(22.0, 16.0))
    panels = list(axes.flat)

    any_anomalous_top, any_reconstructed_top = build_vs_commanded_row(
        panels[0:10], run_index, mres, rows, time_s, position_nm,
        sample_period_s, sim_npz,
    )
    any_anomalous_bottom, any_reconstructed_bottom = build_measured_minus_sim_row(
        panels[10:20], run_index, mres, rows, time_s, position_nm,
        sample_period_s, sim_npz,
    )
    any_anomalous = any_anomalous_top or any_anomalous_bottom
    any_reconstructed = any_reconstructed_top or any_reconstructed_bottom

    fig.suptitle(
        f"Run {run_index} — MRES 1/{mres}, "
        f"{CURRENT_SHORT.get(current, current)} "
        "(panels in chronological order: BLOCK_0 → C → D_0.125 → ... → D_200) "
        f"| Top: tracking error vs. reconstructed commanded | Bottom: "
        f"measured − simulated | Model: {MODEL_LABEL}",
        fontsize=11.0,
    )
    footnote_lines = [
        "Top 10 panels: trace − reconstructed commanded trajectory "
        "(command_reconstruction.py). Bottom 10 panels: measured − "
        "simulated directly (model-vs-hardware fidelity, independent of "
        "the commanded trajectory). Zero line shown dashed grey throughout.",
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
            "† D_3.5 and above: commanded line (top row) is a reconstructed "
            "ACCEL/CONST/DECEL ramp -- see README.md, \"Controller-paced "
            "D-rate command reconstruction\""
        )
    fig.text(
        0.5, 0.008, "\n".join(footnote_lines),
        ha='center', va='bottom', fontsize=8.0, color='#b30000',
    )
    bottom_margin = 0.035 + 0.012 * any_anomalous + 0.012 * any_reconstructed
    fig.tight_layout(rect=(0, bottom_margin, 1, 0.96))

    folder_name = f'run_{run_index:02d}_mres_{mres}_{current.lower()}'
    out_dir = SUBPLOT_DIR / folder_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'tracking_error_montage.png'
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
        folder_name = f'run_{run_index:02d}_mres_{mres}_{current.lower()}'
        sim_path = SUBPLOT_DIR / folder_name / 'lugre_simulation.npz'
        sim_npz = np.load(sim_path) if sim_path.exists() else None
        if sim_npz is None:
            print(f'Run {run_index}: WARNING no lugre_simulation.npz found, '
                  f'simulated error omitted')
        out_path = build_montage(
            run_index, mres, current, rows, time_s, position_nm,
            sample_period_s, sim_npz,
        )
        print(f'Run {run_index}: {out_path}')


if __name__ == '__main__':
    main()

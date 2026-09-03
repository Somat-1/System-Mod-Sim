#!/usr/bin/env python3
"""Per-configuration 10-panel torque/friction montage, in the same
chronological layout as plot_block_montage.py (BLOCK_0, C, then the 8 D
plateaus ascending): motor torque, detent torque, and the three LuGre
friction ports (way/nut/sb), decomposed, from the Rev 4.2 model's
simulated response to each block's reconstructed commanded trajectory
(see simulate_block_responses.py / command_reconstruction.py).

These are simulated quantities only -- the hardware has no torque or
friction-force sensor, so there is nothing measured to compare against.
Motor torque, detent torque, and the bearing friction port ("sb") are all
N*m and share the left axis; the two linear friction ports ("way", the
guideway, and "nut", the leadscrew-nut interface) are N and share the
right axis, since the two groups can differ by orders of magnitude.
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / 'data' / 'raw_local'
SUBPLOT_DIR = ROOT / 'rendered_assets' / 'individual_subplots'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'

CURRENT_SHORT = {'I_50pct': '50% I', 'I_100pct': '100% I'}
D_RATES = ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')
SLOW_PLATEAU_D_RATES = ('0.125', '0.375', '1.25')
MODEL_LABEL = 'Rev 4.2 parallel LuGre + nonlinear detent torque'

MOTOR_COLOR = '#1f77b4'
DETENT_COLOR = '#9467bd'
SB_COLOR = '#8c564b'
WAY_COLOR = '#2ca02c'
NUT_COLOR = '#d62728'
AXIS_COLOR = '#333333'
ANNOTATION_BBOX = dict(facecolor='white', edgecolor='none', alpha=0.78, pad=1.5)

ANOMALOUS_D_RATES_BY_MRES = {
    '4': (),
    '2': ('3.5', '9.5'),
    '1': ('3.5', '9.5', '27.5'),
}


def load_log_rows(path: Path):
    rows = []
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            rows.append(row)
    return rows


def style_axis(ax_left, ax_right, left_linthresh, right_linthresh):
    ax_left.grid(True, alpha=0.3, linewidth=0.6)
    ax_left.set_axisbelow(True)
    ax_left.spines['top'].set_visible(False)
    ax_right.spines['top'].set_visible(False)
    for spine in ax_left.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    for spine in ax_right.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    ax_left.set_xlabel('Time since block start (s)', fontsize=8)
    ax_left.set_ylabel('Motor / detent / bearing torque (N·m)', fontsize=7.5)
    ax_right.set_ylabel('Guideway / nut friction force (N)', fontsize=7.5)
    ax_left.tick_params(labelsize=7, colors=AXIS_COLOR)
    ax_right.tick_params(labelsize=7, colors=AXIS_COLOR)
    # Symmetric log scale: each axis overlays traces spanning orders of
    # magnitude (e.g. bearing friction vs. motor torque), so a linear scale
    # makes the smaller trace invisible. linthresh is set per panel from the
    # smallest trace's own peak, so it stays resolved in the linear region
    # while the largest trace's peak is log-compressed rather than blown out.
    ax_left.set_yscale('symlog', linthresh=left_linthresh, linscale=1.0)
    ax_right.set_yscale('symlog', linthresh=right_linthresh, linscale=1.0)


def _peak(sim_npz, block_name, key):
    return float(np.max(np.abs(sim_npz[f'{block_name}_{key}'])))


def _linthresh(peaks, floor=1.0e-9):
    nonzero = [p for p in peaks if p > floor]
    return max(min(nonzero), floor) if nonzero else floor


def plot_block(ax_left, sim_npz, block_name):
    ax_right = ax_left.twinx()
    t = sim_npz[f'{block_name}_time_s']
    left_linthresh = _linthresh([
        _peak(sim_npz, block_name, 'motor_torque_Nm'),
        _peak(sim_npz, block_name, 'detent_torque_Nm'),
        _peak(sim_npz, block_name, 'friction_sb_Nm'),
    ])
    right_linthresh = _linthresh([
        _peak(sim_npz, block_name, 'friction_way_N'),
        _peak(sim_npz, block_name, 'friction_nut_N'),
    ])
    handles = []
    handles += ax_left.plot(
        t, sim_npz[f'{block_name}_motor_torque_Nm'], color=MOTOR_COLOR,
        lw=0.7, alpha=0.55, label='Motor torque',
    )
    handles += ax_left.plot(
        t, sim_npz[f'{block_name}_detent_torque_Nm'], color=DETENT_COLOR,
        lw=0.7, alpha=0.55, label='Detent torque',
    )
    handles += ax_left.plot(
        t, sim_npz[f'{block_name}_friction_sb_Nm'], color=SB_COLOR,
        lw=0.7, alpha=0.8, label='Bearing friction (sb)',
    )
    handles += ax_right.plot(
        t, sim_npz[f'{block_name}_friction_way_N'], color=WAY_COLOR,
        lw=0.7, alpha=0.7, label='Guideway friction (way)',
    )
    handles += ax_right.plot(
        t, sim_npz[f'{block_name}_friction_nut_N'], color=NUT_COLOR,
        lw=0.7, alpha=0.8, label='Leadscrew-nut friction (nut)',
    )
    style_axis(ax_left, ax_right, left_linthresh, right_linthresh)
    return handles


def build_montage(run_index, mres, current, sim_npz):
    fig, axes = plt.subplots(2, 5, figsize=(22.0, 8.0))
    panels = list(axes.flat)

    ax_ref = panels[0]
    plot_block(ax_ref, sim_npz, 'BLOCK_0_START')
    ax_ref.set_title('BLOCK_0 (start)', fontsize=9)

    ax_c = panels[1]
    handles = plot_block(ax_c, sim_npz, 'C')
    ax_c.set_title('C (creep/settling)', fontsize=9)
    ax_c.legend(handles=handles, loc='best', fontsize=6.0, framealpha=0.9)

    anomalous_rates = ANOMALOUS_D_RATES_BY_MRES.get(str(mres), ())
    any_anomalous = False
    any_reconstructed = False
    for panel_index, (ax, rate) in enumerate(zip(panels[2:], D_RATES)):
        block_name = f'D_{rate}'
        panel_handles = plot_block(ax, sim_npz, block_name)
        ax.set_title(f'D {rate} full-steps/s', fontsize=9, color=AXIS_COLOR)
        if panel_index == 0:
            ax.legend(
                handles=panel_handles, loc='best', fontsize=6.0,
                framealpha=0.9,
            )
        if rate not in SLOW_PLATEAU_D_RATES:
            any_reconstructed = True
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

    fig.suptitle(
        f"Run {run_index} — MRES 1/{mres}, "
        f"{CURRENT_SHORT.get(current, current)} "
        "(panels in chronological order: BLOCK_0 → C → D_0.125 → ... → D_200) "
        f"| Simulated torque/friction: {MODEL_LABEL}",
        fontsize=11.5,
    )
    footnote_lines = [
        "All traces are simulated (no hardware torque/friction sensor); "
        "each block is driven independently from rest by its own "
        "reconstructed commanded trajectory -- see README.md.",
    ]
    if any_anomalous:
        footnote_lines.append(
            "✗ marked D panels: the controller completed this plateau "
            "far faster than commanded (execution-timing artifact); the "
            "underlying command reconstruction is compressed accordingly "
            "and these panels are excluded from analysis -- see README.md"
        )
    if any_reconstructed:
        footnote_lines.append(
            "D_3.5 and above are driven by the reconstructed ACCEL/CONST/"
            "DECEL ramp, not a zero-order-hold step -- see README.md, "
            "\"Controller-paced D-rate command reconstruction\""
        )
    fig.text(
        0.5, 0.012, "\n".join(footnote_lines),
        ha='center', va='bottom', fontsize=8.0, color='#b30000',
    )
    bottom_margin = 0.05 + 0.02 * any_anomalous + 0.02 * any_reconstructed
    fig.tight_layout(rect=(0, bottom_margin, 1, 0.93))

    folder_name = f'run_{run_index:02d}_mres_{mres}_{current.lower()}'
    out_dir = SUBPLOT_DIR / folder_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'torque_montage.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main() -> None:
    rows = load_log_rows(LOG_PATH)

    for run_index in range(1, 7):
        start_row = next(
            r for r in rows if r['event'] == 'RUN_CONFIG'
            and r['run_index'] == str(run_index)
        )
        mres, current = start_row['mres'], start_row['current']
        folder_name = f'run_{run_index:02d}_mres_{mres}_{current.lower()}'
        sim_path = SUBPLOT_DIR / folder_name / 'lugre_simulation.npz'
        if not sim_path.exists():
            print(f'Run {run_index}: WARNING no lugre_simulation.npz found, skipping')
            continue
        sim_npz = np.load(sim_path)
        if f'BLOCK_0_START_motor_torque_Nm' not in sim_npz:
            print(
                f'Run {run_index}: WARNING lugre_simulation.npz has no torque '
                'data -- rerun simulate_block_responses.py first, skipping'
            )
            continue
        out_path = build_montage(run_index, mres, current, sim_npz)
        print(f'Run {run_index}: {out_path}')


if __name__ == '__main__':
    main()

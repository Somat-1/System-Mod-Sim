#!/usr/bin/env python3
"""One-off diagnostic (not part of the regular pipeline, not written to
lugre_simulation.npz): investigate a reported guideway ("way") friction
polarity dip during otherwise one-directional stepping in
torque_montage.png, for run 2 (MRES 1/4, 100% I) only.

Re-simulates run 2's 11 blocks exactly as simulate_block_responses.py
does (same command reconstruction, same per-run current-dependent T_hold),
but additionally captures the "way" port's actual bristle state z(t) and
port velocity v(t), and compares z(t) against the LuGre quasi-steady
prediction z_ss(t) = sign(v) * g(v) / sigma0 (the value z would settle to
if v were held constant) -- a real, large, sustained disagreement between
z and z_ss while |v| is not near zero would indicate a genuine bug (e.g.
a state persistence/sign error); z lagging z_ss only transiently during
accel/decel or briefly overshooting through zero right at a direction
reversal is expected LuGre presliding behavior, not a bug.
"""

from __future__ import annotations

import csv
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = ROOT / 'data' / 'raw_local'
SUBPLOT_DIR = ROOT / 'rendered_assets' / 'individual_subplots'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
SCRIPTS_DIR = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3\scripts")

RUN_INDEX = 2
FILETIME_UNIX_EPOCH = 116444736000000000
CONTROLLER_CLOCK_SKEW_S = 0.319
BLOCK_NAMES = (
    'BLOCK_0_START', 'C',
    'D_0.125', 'D_0.375', 'D_1.25', 'D_3.5', 'D_9.5', 'D_27.5', 'D_70', 'D_200',
)
D_RATES = ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')
RTOL = 1.0e-6
ATOL = np.array([
    1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-12, 1.0e-12,
    1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-9, 1.0e-9,
    1.0e-11, 1.0e-11, 1.0e-9,
])
T_HOLD_BASE_NM = 0.060
I_RATED_BASE_A = 0.400
K_T_NM_PER_A = T_HOLD_BASE_NM / (np.sqrt(2.0) * I_RATED_BASE_A)
PEAK_MA_BY_CURRENT_LEVEL = {'I_50pct': 200, 'I_100pct': 400}

Z_COLOR = '#1f77b4'
ZSS_COLOR = '#d62728'
V_COLOR = '#9a9a9a'
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


def simulate_one_block_diag(run_index: int, block_name: str) -> dict:
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import (
        N_Q, N_STATES, PORTS, LuGreModelRev42, _port_values, load_parameters,
        lugre_terms,
    )
    sys.path.insert(0, str(SCRIPTS_DIR))
    from command_reconstruction import (
        reconstruct_segments, trapezoid_fraction, trapezoid_fraction_array,
    )

    start_epoch_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    start_s, end_s = find_block(rows, run_index, block_name)
    segments = reconstruct_segments(rows, run_index, start_s, end_s)

    run_config = next(
        r for r in rows if r['event'] == 'RUN_CONFIG'
        and r['run_index'] == str(run_index)
    )
    peak_ma = PEAK_MA_BY_CURRENT_LEVEL[run_config['current']]
    t_hold_run = K_T_NM_PER_A * np.sqrt(2.0) * (peak_ma / 1000.0)
    params = load_parameters()
    params['T_hold'] = t_hold_run
    model = LuGreModelRev42(parameters=params, enforce_interface_power=False)
    way_index = PORTS.index('way')
    way_values = _port_values(model.p, 'way')
    sigma0_way = way_values[0]

    state = np.zeros(N_STATES)
    t_parts, v_parts, z_parts, zss_parts, f_parts = [], [], [], [], []
    for seg in segments:
        t0, t1 = seg['t0'], seg['t1']
        if t1 <= t0:
            continue
        if seg['kind'] == 'hold':
            theta_command = seg['value_start'] * 2.0 * np.pi

            def rhs(t, y, cmd=theta_command):
                return model.rhs(t, y, cmd)
        else:
            value_start, value_end = seg['value_start'], seg['value_end']
            t_accel_s, duration_s = seg['t_accel_s'], seg['duration_s']

            def rhs(t, y, t0=t0, value_start=value_start, value_end=value_end,
                     t_accel_s=t_accel_s, duration_s=duration_s):
                frac = trapezoid_fraction(t - t0, t_accel_s, duration_s)
                cmd = (value_start + (value_end - value_start) * frac) * 2.0 * np.pi
                return model.rhs(t, y, cmd)
        sol = solve_ivp(
            rhs, (t0, t1), state, method='Radau',
            jac=lambda t, y: model.analytical_linearization(y)[0],
            t_eval=np.linspace(t0, t1, max(2, int((t1 - t0) * 200))),
            rtol=RTOL, atol=ATOL,
        )
        if not sol.success:
            raise RuntimeError(f'{block_name} segment failed: {sol.message}')
        state = sol.y[:, -1]

        velocity_arr = sol.y[N_Q:2 * N_Q, :]
        v_way = model.jacobians['way'] @ velocity_arr
        z_way = sol.y[2 * N_Q + way_index, :]
        sigma0, sigma1, sigma2, fc, fs, vs, eps = way_values
        g = fc + (fs - fc) * np.exp(-(v_way / vs) ** 2)
        z_ss_way = np.sign(v_way) * g / sigma0
        friction_way = lugre_terms(v_way, z_way, *way_values)[0]

        t_parts.append(sol.t)
        v_parts.append(v_way)
        z_parts.append(z_way)
        zss_parts.append(z_ss_way)
        f_parts.append(friction_way)

    return {
        'block_name': block_name,
        't': np.concatenate(t_parts),
        'v_way': np.concatenate(v_parts),
        'z_way': np.concatenate(z_parts),
        'z_ss_way': np.concatenate(zss_parts),
        'friction_way': np.concatenate(f_parts),
        'sigma0_way': sigma0_way,
    }


def style_axis(ax_left, ax_right):
    ax_left.grid(True, alpha=0.3, linewidth=0.6)
    ax_left.set_axisbelow(True)
    ax_left.spines['top'].set_visible(False)
    ax_right.spines['top'].set_visible(False)
    for spine in ax_left.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    for spine in ax_right.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    ax_left.axhline(0.0, color=AXIS_COLOR, lw=0.5, alpha=0.4)
    ax_left.set_xlabel('Time since block start (s)', fontsize=8)
    ax_left.set_ylabel('Bristle deflection z (µm)', fontsize=7.5)
    ax_right.set_ylabel('Port velocity v_way (mm/s)', fontsize=7.5)
    ax_left.tick_params(labelsize=7, colors=AXIS_COLOR)
    ax_right.tick_params(labelsize=7, colors=AXIS_COLOR)


def plot_block(ax_left, result):
    ax_right = ax_left.twinx()
    t = result['t']
    handles = []
    handles += ax_right.plot(
        t, result['v_way'] * 1.0e3, color=V_COLOR, lw=0.6, alpha=0.6,
        linestyle=':', label='v_way (mm/s)', zorder=1,
    )
    handles += ax_left.plot(
        t, result['z_ss_way'] * 1.0e6, color=ZSS_COLOR, lw=1.1, alpha=0.85,
        linestyle='--', label='z_ss = sign(v)·g(v)/σ0', zorder=2,
    )
    handles += ax_left.plot(
        t, result['z_way'] * 1.0e6, color=Z_COLOR, lw=1.0, zorder=3,
        label='z (actual bristle state)',
    )
    style_axis(ax_left, ax_right)
    return handles


def main() -> None:
    print(f'Diagnosing run {RUN_INDEX}, way port, {len(BLOCK_NAMES)} blocks', flush=True)
    started = time.perf_counter()
    results = {}
    workers = max(1, min(11, (os.cpu_count() or 4) - 2))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(simulate_one_block_diag, RUN_INDEX, block_name): block_name
            for block_name in BLOCK_NAMES
        }
        for future in as_completed(futures):
            block_name = futures[future]
            results[block_name] = future.result()
            print(f'  done {block_name} ({time.perf_counter()-started:.1f}s)', flush=True)

    sigma0_way = results['BLOCK_0_START']['sigma0_way']

    fig, axes = plt.subplots(2, 5, figsize=(22.0, 8.0))
    panels = list(axes.flat)

    handles = plot_block(panels[0], results['BLOCK_0_START'])
    panels[0].set_title('BLOCK_0 (start)', fontsize=9)
    panels[0].legend(handles=handles, loc='best', fontsize=6.5, framealpha=0.9)

    plot_block(panels[1], results['C'])
    panels[1].set_title('C (creep/settling)', fontsize=9)

    for ax, rate in zip(panels[2:], D_RATES):
        plot_block(ax, results[f'D_{rate}'])
        ax.set_title(f'D {rate} full-steps/s', fontsize=9)

    fig.suptitle(
        f'Run {RUN_INDEX} — MRES 1/4, 100% I | Guideway ("way") LuGre bristle '
        f'state diagnostic: actual z vs. quasi-steady z_ss = sign(v)·g(v)/σ0 '
        f'(σ0_way = {sigma0_way:.3g} N/m)',
        fontsize=11.5,
    )
    fig.text(
        0.5, 0.012,
        'z tracking z_ss during cruise, with a brief lag/overshoot through '
        'zero at accel/decel and direction reversals, is expected LuGre '
        'presliding behavior -- a sustained disagreement in sign while '
        '|v_way| is not near zero would indicate a state-persistence bug.',
        ha='center', va='bottom', fontsize=8.0, color='#b30000',
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))

    out_path = SUBPLOT_DIR / 'run_02_mres_4_i_100pct' / 'way_bristle_diagnostic.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}', flush=True)

    # Numeric check: sign(friction) vs sign(v) disagreement while |v| is
    # not near zero (an actual bug would show large, sustained fractions).
    print("\nSign-disagreement check (|v_way| > 1% of that block's own peak):")
    for block_name, result in results.items():
        v, f = result['v_way'], result['friction_way']
        thresh = 0.01 * np.max(np.abs(v)) if np.max(np.abs(v)) > 0 else 0.0
        moving = np.abs(v) > thresh
        if not np.any(moving):
            continue
        disagree = np.sign(f[moving]) != np.sign(v[moving])
        # kinetic friction opposes velocity, so disagreement means friction
        # points the SAME way as velocity (i.e. sign(f) == sign(v))
        same_sign_frac = np.mean(np.sign(f[moving]) == np.sign(v[moving]))
        print(f'  {block_name:12s} fraction of moving samples with friction '
              f'pointing WITH velocity (unexpected): {same_sign_frac:.4f}')

    print(f'Total wall time: {time.perf_counter()-started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

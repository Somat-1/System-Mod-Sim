#!/usr/bin/env python3
"""Step 2: run the C block (run 2) at epsilon=1e-9 and epsilon=1e-12,
compare simulated position during the 60 s dwells. If they agree, epsilon
is inert at the scale that matters. If they differ by tens of nm or more,
epsilon's regularization is spuriously relaxing the bristle state at rest
and the smaller value should be used when fitting from C.
"""
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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

V3_ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = V3_ROOT / 'data' / 'raw_local'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
SCRIPTS_DIR = V3_ROOT / 'scripts'
OUT_DIR = Path(__file__).resolve().parent

RUN_INDEX = 2
FILETIME_UNIX_EPOCH = 116444736000000000
CONTROLLER_CLOCK_SKEW_S = 0.319
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


def simulate_c_block(epsilon: float):
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import N_STATES, LuGreModelRev42, load_parameters
    sys.path.insert(0, str(SCRIPTS_DIR))
    from command_reconstruction import reconstruct_segments, trapezoid_fraction

    start_epoch_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    start_s, end_s = find_block(rows, RUN_INDEX, 'C')
    segments = reconstruct_segments(rows, RUN_INDEX, start_s, end_s, 'C')

    run_config = next(
        r for r in rows if r['event'] == 'RUN_CONFIG'
        and r['run_index'] == str(RUN_INDEX)
    )
    peak_ma = PEAK_MA_BY_CURRENT_LEVEL[run_config['current']]
    t_hold_run = K_T_NM_PER_A * np.sqrt(2.0) * (peak_ma / 1000.0)
    params = load_parameters()
    params['T_hold'] = t_hold_run
    params['smooth_velocity_epsilon'] = epsilon
    model = LuGreModelRev42(parameters=params, enforce_interface_power=False)

    state = np.zeros(N_STATES)
    t_parts, x_parts = [], []
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
            raise RuntimeError(f'C block segment failed at eps={epsilon}: {sol.message}')
        state = sol.y[:, -1]
        t_parts.append(sol.t)
        x_parts.append(sol.y[5, :])

    return np.concatenate(t_parts), np.concatenate(x_parts) * 1.0e6  # um


def main():
    started = time.perf_counter()
    results = {}
    for eps in (1.0e-9, 1.0e-12):
        t, x_um = simulate_c_block(eps)
        results[eps] = (t, x_um)
        print(f'eps={eps:.0e} done ({time.perf_counter()-started:.1f}s)', flush=True)

    t9, x9 = results[1.0e-9]
    t12, x12 = results[1.0e-12]
    # common grid for comparison
    t_common = t9 if len(t9) <= len(t12) else t12
    x9_i = np.interp(t_common, t9, x9)
    x12_i = np.interp(t_common, t12, x12)
    diff_um = x9_i - x12_i

    print(f'\nMax |difference| over whole C block: {np.max(np.abs(diff_um)):.4f} um '
          f'({np.max(np.abs(diff_um))*1000:.2f} nm)')
    print(f'Difference at end of block: {diff_um[-1]:.5f} um ({diff_um[-1]*1000:.2f} nm)')

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.0), sharex=True)
    axes[0].plot(t9, x9, label='eps=1e-9', color='#1f77b4', lw=1.0)
    axes[0].plot(t12, x12, label='eps=1e-12', color='#d62728', lw=1.0, linestyle='--')
    axes[0].set_ylabel('Simulated position (µm)')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title(f'Run {RUN_INDEX}, C block: position at two epsilon values')

    axes[1].plot(t_common, diff_um * 1000.0, color='#2ca02c', lw=1.0)
    axes[1].axhline(0.0, color='#9a9a9a', lw=0.7)
    axes[1].set_ylabel('Difference (nm)\n(eps=1e-9 minus eps=1e-12)')
    axes[1].set_xlabel('Time since block start (s)')
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = OUT_DIR / 'step2_epsilon_dwell_comparison.png'
    fig.savefig(out_path, dpi=150)
    print(f'Saved {out_path}')
    print(f'Total wall time: {time.perf_counter()-started:.1f}s')


if __name__ == '__main__':
    main()

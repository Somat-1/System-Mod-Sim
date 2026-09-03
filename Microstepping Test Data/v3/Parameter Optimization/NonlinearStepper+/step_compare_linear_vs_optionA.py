#!/usr/bin/env python3
"""Compare the existing linear-drive Rev 4.2 model against Option A (the
nonlinear current-projection drive law, see lugre_model_rev42_optionA.py)
on run 2, across the same 4 blocks used in the earlier detent-vs-friction
ablation (D_0.125, D_3.5, D_70, D_200): D_0.125 as the small-lag sanity
check (Option A should reduce to ~linear behavior there), D_70/D_200 as
the blocks with the large, parameter-insensitive model-vs-measurement
mismatch under investigation.

For each (block, model) pair: full-block simulation -> cost_vs_measured
(RMS vs. hardware, same metric as step6/step7), plus max|theta_err| over
the whole block (in degrees, relative to the pi/(2*N_r) = 1.8-degree
pull-out threshold for N_r=50). For the three controller-paced blocks
(D_3.5, D_70, D_200) also extracts a 10-detent-cycle window centered
mid-cruise (same windowing as torque_diagnostics/plot_cruise_zoom.py) for
a time-series comparison plot of tracking error and motor torque, both
models overlaid.
"""
from __future__ import annotations

import csv
import json
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

V3_ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = V3_ROOT / 'data' / 'raw_local'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
SCRIPTS_DIR = V3_ROOT / 'scripts'
OUT_DIR = HERE

RUN_INDEX = 2
BLOCKS = ('0.125', '3.5', '70', '200')
TRAPEZOID_BLOCKS = ('3.5', '70', '200')  # blocks with a real cruise phase to window into
MODEL_KINDS = ('linear', 'optionA')
N_CYCLES = 10
MOTOR_FULL_STEPS_PER_REV = 200

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

_CACHED = None


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


def _worker_init(time_s, position_nm, start_epoch_s, sample_period_s, rows):
    global _CACHED
    _CACHED = {
        'time_s': time_s, 'position_nm': position_nm,
        'sample_period_s': sample_period_s, 'rows': rows,
    }


def simulate_variant(run_index, block_name, model_kind):
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import N_STATES, LuGreModelRev42, load_parameters
    sys.path.insert(0, str(SCRIPTS_DIR))
    from command_reconstruction import reconstruct_segments, trapezoid_fraction
    from lugre_model_rev42_optionA import LuGreModelRev42OptionA

    rows = _CACHED['rows']
    time_s = _CACHED['time_s']
    position_nm = _CACHED['position_nm']
    sample_period_s = _CACHED['sample_period_s']

    start_s, end_s = find_block(rows, run_index, block_name)
    segments = reconstruct_segments(rows, run_index, start_s, end_s, block_name)

    run_config = next(
        r for r in rows if r['event'] == 'RUN_CONFIG' and r['run_index'] == str(run_index)
    )
    peak_ma = PEAK_MA_BY_CURRENT_LEVEL[run_config['current']]
    t_hold_run = K_T_NM_PER_A * np.sqrt(2.0) * (peak_ma / 1000.0)
    params = load_parameters()
    params['T_hold'] = t_hold_run
    k_em = params['N_r'] * params['T_hold']

    if model_kind == 'linear':
        model = LuGreModelRev42(parameters=params, enforce_interface_power=False)
    else:
        model = LuGreModelRev42OptionA(parameters=params, enforce_interface_power=False)

    first_trapezoid = next((s for s in segments if s['kind'] == 'trapezoid'), None)
    window_s = None
    if first_trapezoid is not None:
        rate_hz = float(block_name.split('_')[1])
        window_s = N_CYCLES / rate_hz
        cruise_mid = (
            first_trapezoid['t0'] + first_trapezoid['t_accel_s']
            + 0.5 * first_trapezoid['duration_s']
        )
        win_t0, win_t1 = cruise_mid - 0.5 * window_s, cruise_mid + 0.5 * window_s

    state = np.zeros(15)
    t_parts, x_parts = [], []
    max_abs_err_rad = 0.0
    win_t, win_theta_err, win_motor_torque, win_theta_m = [], [], [], []

    for seg in segments:
        t0, t1 = seg['t0'], seg['t1']
        if t1 <= t0:
            continue

        if seg['kind'] == 'hold':
            def cmd_at(t, value=seg['value_start']):
                return value * 2.0 * np.pi
        else:
            def cmd_at(t, t0=t0, value_start=seg['value_start'], value_end=seg['value_end'],
                       t_accel_s=seg['t_accel_s'], duration_s=seg['duration_s']):
                frac = trapezoid_fraction(t - t0, t_accel_s, duration_s)
                return (value_start + (value_end - value_start) * frac) * 2.0 * np.pi

        def rhs(t, y, cmd_at=cmd_at):
            return model.rhs(t, y, cmd_at(t))

        if model_kind == 'linear':
            jac = lambda t, y: model.analytical_linearization(y)[0]
        else:
            jac = lambda t, y, cmd_at=cmd_at: model.analytical_linearization(y, cmd_at(t))[0]

        if first_trapezoid is not None and seg is first_trapezoid:
            # t_eval only controls dense-output interpolation density (cheap --
            # solve_ivp's adaptive Radau stepping is unaffected), but the default
            # 200 samples/s is only ~2.9 samples/cycle at D_70's 70 Hz detent
            # forcing and exactly Nyquist at D_200's 200 Hz -- badly under-resolved
            # for a faithful windowed waveform plot. Boost density for this
            # segment (the only one that gets windowed/plotted) to >=40
            # samples/detent-cycle.
            rate_hz_for_seg = float(block_name.split('_')[1])
            samples_per_sec = max(200.0, rate_hz_for_seg * 40.0)
            n_eval = max(2, int((t1 - t0) * samples_per_sec))
        else:
            n_eval = max(2, int((t1 - t0) * 200))
        sol = solve_ivp(
            rhs, (t0, t1), state, method='Radau', jac=jac,
            t_eval=np.linspace(t0, t1, n_eval), rtol=RTOL, atol=ATOL,
        )
        if not sol.success:
            raise RuntimeError(f'{block_name} [{model_kind}] segment failed: {sol.message}')
        state = sol.y[:, -1]
        t_parts.append(sol.t)
        x_parts.append(sol.y[5, :])

        theta_cmd_arr = np.array([cmd_at(tt) for tt in sol.t])
        theta_err_arr = theta_cmd_arr - sol.y[0, :]
        max_abs_err_rad = max(max_abs_err_rad, float(np.max(np.abs(theta_err_arr))))

        if first_trapezoid is not None and seg is first_trapezoid:
            mask = (sol.t >= win_t0) & (sol.t <= win_t1)
            win_t.append(sol.t[mask])
            win_theta_err.append(theta_err_arr[mask])
            win_theta_m.append(sol.y[0, mask])
            if model_kind == 'linear':
                win_motor_torque.append(k_em * theta_err_arr[mask])
            else:
                win_motor_torque.append(model.motor_torque(theta_err_arr[mask]))

    sim_t = np.concatenate(t_parts)
    sim_y_um = np.concatenate(x_parts) * 1.0e6

    t_meas, y_meas = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
    sim_on_meas_grid = np.interp(t_meas, sim_t, sim_y_um, left=sim_y_um[0], right=sim_y_um[-1])
    cost = float(np.sqrt(np.mean((y_meas - sim_on_meas_grid) ** 2)))

    result = {
        'cost_um': cost,
        'max_abs_err_deg': float(np.degrees(max_abs_err_rad)),
        'k_em': k_em, 't_hold_run': t_hold_run,
    }
    if win_t:
        result['window_t'] = np.concatenate(win_t)
        result['window_theta_err_deg'] = np.degrees(np.concatenate(win_theta_err))
        result['window_motor_torque_Nm'] = np.concatenate(win_motor_torque)
        window_theta_m = np.concatenate(win_theta_m)
        result['window_detent_Nm'] = params['T_d'] * np.sin(4.0 * params['N_r'] * window_theta_m)
        result['window_net_drive_Nm'] = result['window_motor_torque_Nm'] - result['window_detent_Nm']
        result['window_s'] = window_s
        result['detent_hz'] = float(block_name.split('_')[1])
    return result


def plot_static_torque_law(params):
    n_r, t_hold = params['N_r'], params['T_hold']
    k_em = n_r * t_hold
    pullout_deg = np.degrees(np.pi / (2.0 * n_r))
    err_deg = np.linspace(-3.0 * pullout_deg, 3.0 * pullout_deg, 2000)
    err_rad = np.radians(err_deg)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.plot(err_deg, k_em * err_rad, color='#d62728', lw=1.6, label='Linear (existing): k_em*theta_err')
    ax.plot(err_deg, t_hold * np.sin(n_r * err_rad), color='#1f77b4', lw=1.8,
             label='Option A: T_hold*sin(N_r*theta_err)')
    ax.axhline(t_hold, color='#1f77b4', lw=0.8, ls=':')
    ax.axhline(-t_hold, color='#1f77b4', lw=0.8, ls=':')
    for sign in (-1, 1):
        ax.axvline(sign * pullout_deg, color='#555555', lw=0.8, ls='--')
    ax.text(pullout_deg, -t_hold * 1.15, f'  1 full step\n  ({pullout_deg:.2f}\u00b0)',
            fontsize=8, color='#555555', va='top')
    ax.set_xlabel('Tracking error theta_err = theta_cmd - theta_m (degrees)')
    ax.set_ylabel('Motor drive torque (N\u00b7m)')
    ax.set_title(f'Motor drive torque law: linear vs. Option A (T_hold={t_hold:.4f} N\u00b7m, N_r={n_r:.0f})')
    ax.axhline(0.0, color='#9a9a9a', lw=0.8)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=9)
    fig.tight_layout()
    out_path = OUT_DIR / 'torque_law_linear_vs_optionA.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_block_comparison(block_name, linear_result, optionA_result):
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True)
    t_ms = (linear_result['window_t'] - linear_result['window_t'][0]) * 1000.0
    pullout_deg = optionA_result.get('pullout_deg')

    axes[0].plot(t_ms, linear_result['window_theta_err_deg'], color='#d62728', lw=1.0, label='Linear')
    axes[0].plot(t_ms, optionA_result['window_theta_err_deg'], color='#1f77b4', lw=1.0, label='Option A')
    if pullout_deg:
        axes[0].axhline(pullout_deg, color='#555555', lw=0.8, ls='--')
        axes[0].axhline(-pullout_deg, color='#555555', lw=0.8, ls='--')
    axes[0].axhline(0.0, color='#9a9a9a', lw=0.7)
    axes[0].set_ylabel('Tracking error\ntheta_err (degrees)')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc='best', fontsize=8)

    axes[1].plot(t_ms, linear_result['window_motor_torque_Nm'], color='#d62728', lw=1.0, label='Linear')
    axes[1].plot(t_ms, optionA_result['window_motor_torque_Nm'], color='#1f77b4', lw=1.0, label='Option A')
    axes[1].axhline(0.0, color='#9a9a9a', lw=0.7)
    axes[1].set_ylabel('Motor drive torque (N\u00b7m)')
    axes[1].set_xlabel('Time within window (ms)')
    axes[1].grid(True, alpha=0.3)

    fig.suptitle(
        f'D_{block_name} full-steps/s -- run 2, cruise window '
        f'({optionA_result["detent_hz"]:.1f} Hz forcing, {N_CYCLES} cycles): '
        'linear vs. Option A drive law',
        fontsize=11.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out_path = OUT_DIR / f'compare_D_{block_name}.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    print('Parsing IDS + log once...', flush=True)
    t0 = time.perf_counter()
    time_s, position_nm, start_epoch_s, sample_period_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    print(f'  done in {time.perf_counter() - t0:.1f}s', flush=True)

    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import load_parameters
    base_params = load_parameters()
    static_plot_path = plot_static_torque_law(
        {'N_r': base_params['N_r'], 'T_hold': base_params['T_hold']}
    )
    print(f'Saved {static_plot_path}', flush=True)
    pullout_deg = float(np.degrees(np.pi / (2.0 * base_params['N_r'])))

    jobs = {(block, kind): (RUN_INDEX, f'D_{block}', kind)
            for block in BLOCKS for kind in MODEL_KINDS}
    workers = min(len(jobs), 12)
    print(f'{len(jobs)} jobs, {workers} workers', flush=True)

    results = {}
    started = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_worker_init,
        initargs=(time_s, position_nm, start_epoch_s, sample_period_s, rows),
    ) as pool:
        futures = {pool.submit(simulate_variant, *v): k for k, v in jobs.items()}
        for fut in as_completed(futures):
            key = futures[fut]
            res = fut.result()
            res['pullout_deg'] = pullout_deg
            results[key] = res
            print(f'{key}: cost={res["cost_um"]:.4f} um, max|err|={res["max_abs_err_deg"]:.4f} deg '
                  f'(t={time.perf_counter() - started:.1f}s)', flush=True)

    print('\n=== Summary: cost vs. measurement (um) and peak tracking error (deg) ===')
    print(f'{"block":>8s}  {"linear cost":>12s}  {"optA cost":>12s}  '
          f'{"linear max|err|":>16s}  {"optA max|err|":>14s}  {"pull-out @":>10s}')
    for block in BLOCKS:
        lin, opt = results[(block, 'linear')], results[(block, 'optionA')]
        print(f'D_{block:>6s}  {lin["cost_um"]:>12.4f}  {opt["cost_um"]:>12.4f}  '
              f'{lin["max_abs_err_deg"]:>16.4f}  {opt["max_abs_err_deg"]:>14.4f}  {pullout_deg:>9.2f}\u00b0')

    for block in TRAPEZOID_BLOCKS:
        out_path = plot_block_comparison(block, results[(block, 'linear')], results[(block, 'optionA')])
        print(f'Saved {out_path}', flush=True)

    serializable = {
        f'{k[0]}|{k[1]}': {kk: vv for kk, vv in v.items() if not isinstance(vv, np.ndarray)}
        for k, v in results.items()
    }
    out_json = OUT_DIR / 'compare_linear_vs_optionA.json'
    out_json.write_text(json.dumps(serializable, indent=2))
    print(f'Saved {out_json}', flush=True)
    print(f'Total wall time: {time.perf_counter() - started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

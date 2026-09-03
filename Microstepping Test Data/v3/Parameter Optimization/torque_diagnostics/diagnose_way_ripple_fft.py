#!/usr/bin/env python3
"""One-off diagnostic (not part of the regular pipeline): FFT the
simulated guideway ("way") port velocity during the cruise (constant
commanded-speed) phase of D_27.5, D_70, and D_200 for run 2 (MRES 1/4,
100% I) only, to distinguish detent-torque forcing from a freely-rung
resonance.

If the dominant ripple frequency tracks the commanded rate (27.5 Hz,
70 Hz, 200 Hz respectively), the ripple is locked to detent forcing:
the detent torque T_d*sin(4*N_r*theta_m) repeats 4*N_r=200 times per
rotor revolution, and a commanded rate of R full-steps/s corresponds to
R/200 rev/s of rotor speed, so its temporal ripple frequency is
(R/200 rev/s) * 200 cycles/rev = R Hz -- i.e. the forcing frequency in Hz
equals the commanded full-steps/s rate exactly.

If instead the same frequency shows up regardless of commanded rate,
that is a freely-rung structural/friction-coupled mode being excited,
not detent forcing -- and T_d has no business being tuned against it.
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
RATES = ('27.5', '70', '200')
FFT_SAMPLE_HZ = 5000.0
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


def cruise_v_way(run_index: int, rate: str) -> dict:
    """Simulate D_<rate> at high output sample rate and return v_way(t)
    restricted to the first (positive-direction) trapezoid segment's
    cruise sub-window (excludes accel/decel ramps)."""
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import (
        N_Q, N_STATES, PORTS, LuGreModelRev42, load_parameters,
    )
    sys.path.insert(0, str(SCRIPTS_DIR))
    from command_reconstruction import reconstruct_segments, trapezoid_fraction

    block_name = f'D_{rate}'
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
    way_jacobian = model.jacobians['way']

    first_trapezoid = next(s for s in segments if s['kind'] == 'trapezoid')
    cruise_t0 = first_trapezoid['t0'] + first_trapezoid['t_accel_s']
    cruise_t1 = cruise_t0 + first_trapezoid['duration_s']

    state = np.zeros(N_STATES)
    cruise_t, cruise_v = [], []
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
        n_eval = max(2, int(round((t1 - t0) * FFT_SAMPLE_HZ)))
        sol = solve_ivp(
            rhs, (t0, t1), state, method='Radau',
            jac=lambda t, y: model.analytical_linearization(y)[0],
            t_eval=np.linspace(t0, t1, n_eval),
            rtol=RTOL, atol=ATOL,
        )
        if not sol.success:
            raise RuntimeError(f'{block_name} segment failed: {sol.message}')
        state = sol.y[:, -1]

        if seg is first_trapezoid:
            velocity_arr = sol.y[N_Q:2 * N_Q, :]
            v_way = way_jacobian @ velocity_arr
            mask = (sol.t >= cruise_t0) & (sol.t <= cruise_t1)
            cruise_t.append(sol.t[mask])
            cruise_v.append(v_way[mask])

    return {
        'rate': rate,
        't': np.concatenate(cruise_t),
        'v_way': np.concatenate(cruise_v),
    }


def main() -> None:
    print(f'FFT-ing v_way cruise phase, run {RUN_INDEX}, rates {RATES}', flush=True)
    started = time.perf_counter()
    results = {}
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures = {
            pool.submit(cruise_v_way, RUN_INDEX, rate): rate for rate in RATES
        }
        for future in as_completed(futures):
            rate = futures[future]
            results[rate] = future.result()
            print(f'  done D_{rate} ({time.perf_counter()-started:.1f}s)', flush=True)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.5))
    print('\nDominant cruise-phase ripple frequency in v_way:')
    for ax, rate in zip(axes, RATES):
        result = results[rate]
        t, v = result['t'], result['v_way']
        dt = np.median(np.diff(t))
        n = len(v)
        window = np.hanning(n)
        v_ac = (v - np.mean(v)) * window
        spectrum = np.abs(np.fft.rfft(v_ac))
        freqs = np.fft.rfftfreq(n, d=dt)
        # ignore near-DC bins (below 1 Hz) so a slow drift doesn't win
        valid = freqs > 1.0
        peak_idx = np.argmax(spectrum[valid])
        peak_freq = freqs[valid][peak_idx]
        expected_hz = float(rate)
        print(f'  D_{rate:>5s}: dominant peak = {peak_freq:7.2f} Hz   '
              f'(expected if detent-locked: {expected_hz:.1f} Hz)   '
              f'window={t[-1]-t[0]:.2f}s, dt={dt*1000:.3f}ms')

        ax.plot(freqs, spectrum, color='#1f77b4', lw=0.9)
        ax.axvline(expected_hz, color='#d62728', linestyle='--', lw=1.2,
                   label=f'commanded rate = {expected_hz:g} Hz')
        ax.axvline(peak_freq, color='#2ca02c', linestyle=':', lw=1.5,
                   label=f'dominant peak = {peak_freq:.1f} Hz')
        ax.set_xlim(0, max(expected_hz * 2.5, 50))
        ax.set_yscale('log')
        ax.set_title(f'D {rate} full-steps/s', fontsize=10)
        ax.set_xlabel('Frequency (Hz)', fontsize=8)
        ax.set_ylabel('|FFT(v_way - mean)|', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        'Run 2 — MRES 1/4, 100% I | Cruise-phase v_way FFT: detent forcing '
        '(peak tracks commanded rate) vs. a freely-rung mode (peak fixed)',
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out_path = SUBPLOT_DIR / 'run_02_mres_4_i_100pct' / 'way_ripple_fft_diagnostic.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'\nSaved {out_path}', flush=True)
    print(f'Total wall time: {time.perf_counter()-started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

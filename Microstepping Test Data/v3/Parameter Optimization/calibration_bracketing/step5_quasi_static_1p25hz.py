#!/usr/bin/env python3
"""Step 5: extract the position ripple at 1.25 Hz from D_1.25 (run 2),
measured vs simulated, during the cruise phase only. 1.25 Hz is well
below any of the drivetrain resonances seen so far (~140 Hz), so this is
a quasi-static-transmission check: does the DC/low-frequency gain from
detent forcing to stage position match, independent of damping and
independent of the friction parameters (which mostly reshape dynamics
at/near resonance and at velocity reversals, not the quasi-static gain)?

Both measured and simulated data already exist (no new simulation) --
this is pure signal analysis on top of the existing pipeline outputs.
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

import numpy as np
from scipy.integrate import solve_ivp

V3_ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = V3_ROOT / 'data' / 'raw_local'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
SCRIPTS_DIR = V3_ROOT / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
from command_reconstruction import reconstruct_segments, trapezoid_fraction

RUN_INDEX = 2
FORCING_HZ = 1.25
OUT_DIR = Path(__file__).resolve().parent

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


def simulate_d1p25(run_index, rows):
    """Fresh D_1.25 simulation with current parameters (T_d=5 mN*m,
    c_nut=0, current-dependent T_hold) -- the existing lugre_simulation.npz
    predates all three of those and would bias this exact check, since
    T_d directly scales the detent-driven ripple being measured here."""
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import N_STATES, LuGreModelRev42, load_parameters

    start_s, end_s = find_block(rows, run_index, 'D_1.25')
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
            raise RuntimeError(f'D_1.25 segment failed: {sol.message}')
        state = sol.y[:, -1]
        t_parts.append(sol.t)
        x_parts.append(sol.y[5, :])
    return np.concatenate(t_parts), np.concatenate(x_parts) * 1.0e6, segments


def amplitude_at(t, y, freq_hz):
    """Amplitude of y's component at freq_hz via a least-squares
    sin/cos fit (robust to a short/non-power-of-two window, unlike a raw
    FFT bin read)."""
    omega = 2.0 * np.pi * freq_hz
    A = np.column_stack([np.cos(omega * t), np.sin(omega * t), np.ones_like(t)])
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b, _ = coeffs
    return np.hypot(a, b)


def main():
    time_s, position_nm, start_epoch_s, sample_period_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    start_s, end_s = find_block(rows, RUN_INDEX, 'D_1.25')
    segments = reconstruct_segments(rows, RUN_INDEX, start_s, end_s)
    # D_1.25 is software-paced (SLOW_PLATEAU_RATES_FS_S): every segment is a
    # 'hold' (one per pulse), there is no trapezoid segment. The "positive
    # direction" plateau runs from the first pulse to the peak value (the
    # point just before position starts decreasing again); trim 10% off
    # each end of that span to avoid the first/last few pulses' edge effects.
    jump_segments = [
        s for s in segments
        if s['kind'] == 'hold' and s['t0'] == s['t1'] and s['value_start'] != 0.0
    ]
    peak_idx = int(np.argmax([abs(s['value_start']) for s in jump_segments]))
    plateau_t0 = jump_segments[0]['t0']
    plateau_t1 = jump_segments[peak_idx]['t0']
    span = plateau_t1 - plateau_t0
    cruise_t0 = plateau_t0 + 0.10 * span
    cruise_t1 = plateau_t1 - 0.10 * span
    print(f'D_1.25 positive-direction plateau: [{plateau_t0:.3f}, {plateau_t1:.3f}] s, '
          f'{len(jump_segments[:peak_idx+1])} pulses')
    print(f'D_1.25 cruise window (middle 80%): [{cruise_t0:.3f}, {cruise_t1:.3f}] s '
          f'({cruise_t1-cruise_t0:.2f} s, {(cruise_t1-cruise_t0)*FORCING_HZ:.1f} cycles)')

    # --- measured ---
    first_idx = int(round(start_s / sample_period_s))
    idx0 = int(round((start_s + cruise_t0) / sample_period_s))
    idx1 = int(round((start_s + cruise_t1) / sample_period_s))
    y_meas_um = position_nm[idx0:idx1] / 1000.0
    t_meas = np.arange(idx0, idx1) * sample_period_s - (start_s + cruise_t0)
    y_meas_detrended = y_meas_um - np.polyval(np.polyfit(t_meas, y_meas_um, 1), t_meas)
    amp_meas = amplitude_at(t_meas, y_meas_detrended, FORCING_HZ)

    # --- simulated: fresh run with current parameters (see docstring) ---
    started = time.perf_counter()
    sim_t_full, sim_y_full, _ = simulate_d1p25(RUN_INDEX, rows)
    print(f'Simulated D_1.25 in {time.perf_counter()-started:.1f}s')
    mask = (sim_t_full >= cruise_t0) & (sim_t_full <= cruise_t1)
    t_sim = sim_t_full[mask] - cruise_t0
    y_sim = sim_y_full[mask]
    y_sim_detrended = y_sim - np.polyval(np.polyfit(t_sim, y_sim, 1), t_sim)
    amp_sim = amplitude_at(t_sim, y_sim_detrended, FORCING_HZ)

    print(f'\nMeasured amplitude at {FORCING_HZ} Hz:  {amp_meas:.4f} um')
    print(f'Simulated amplitude at {FORCING_HZ} Hz: {amp_sim:.4f} um')
    print(f'Ratio (simulated/measured): {amp_sim/amp_meas:.3f}')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.0))
    axes[0].plot(t_meas, y_meas_detrended, color='#136f63', lw=0.8, label='Measured (detrended)')
    axes[0].plot(t_sim, y_sim_detrended, color='#e67e22', lw=0.9, label='Simulated (detrended)', alpha=0.85)
    axes[0].set_xlabel('Time within cruise window (s)')
    axes[0].set_ylabel('Position ripple (µm)')
    axes[0].set_title(f'D_1.25 cruise-phase ripple (linear trend removed), run {RUN_INDEX}')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    labels = ['Measured', 'Simulated']
    values = [amp_meas, amp_sim]
    axes[1].bar(labels, values, color=['#136f63', '#e67e22'])
    axes[1].set_ylabel(f'Amplitude at {FORCING_HZ} Hz (µm)')
    axes[1].set_title('Quasi-static transmission check: 1.25 Hz ripple amplitude')
    for i, v in enumerate(values):
        axes[1].text(i, v, f'{v:.3f} µm', ha='center', va='bottom', fontsize=9)
    axes[1].grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    out_path = OUT_DIR / 'step5_quasi_static_1p25hz.png'
    fig.savefig(out_path, dpi=150)
    print(f'Saved {out_path}')


if __name__ == '__main__':
    main()

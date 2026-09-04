#!/usr/bin/env python3
"""Figure B (torque-diagnostic methodology): cruise-window zoom for run 2
(MRES 1/4, 100% I), D_9.5 and D_200 -- "is detent a sinusoid at the right
frequency with sane amplitude, and does v_way really cross zero."

Rules applied (see chat discussion that produced this folder):
  1. Never plot a whole block.  Both windows begin 5 s after their respective
     block starts and show 10 detent cycles.  Their durations are calculated
     from the reconstructed (real-time-scaled) cruise speed, rather than the
     nominal D-rate label, so the two columns contain comparable cycle counts.
  2. One signal per row, stacked, sharex, each independently autoscaled --
     never overlaid, never sharing a y-axis across rows of different units.
  3. Linear y-axis with a zero line for every row -- no symlog, since a
     sign-changing quantity's magnitude barely moves while symlog draws it
     leaping across the linear-threshold band.
  4. One unit throughout: motor-side torques (motor drive, detent, bearing
     friction) are converted to an equivalent load-side force via the
     lead-screw transmission (F = T / lead, lead = L / (2*pi)) so the
     motor-drive term, the three friction terms, and the detent term all
     sit in newtons, directly comparable to each other -- explicitly
     labelled motor-side vs. load-side throughout.

Motor drive torque row: T_motor = k_em*(theta_cmd - theta_m), k_em =
N_r*T_hold -- the linearized electromagnetic spring actually integrated by
lugre_model_rev42 (see build_structural_matrices/rhs; there is no sin()
pull-out law for the drive in this model, only for the detent). theta_cmd
isn't part of the ODE state, so it's recomputed here from the same
trapezoid_fraction the solver used internally.
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

V3_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = V3_ROOT / 'data' / 'raw_local'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = V3_ROOT.parents[1] / 'Rev 4' / 'lugre_friction' / 'Rev 4.2' / 'scripts'
SCRIPTS_DIR = V3_ROOT / 'scripts'
OUT_DIR = Path(__file__).resolve().parent

RUN_INDEX = 2
RATES = ('9.5',)
N_CYCLES = 10  # detent cycles shown per cruise window
WINDOW_START_FROM_BLOCK_S = 5.0
FFT_SAMPLE_HZ = 5000.0
FILETIME_UNIX_EPOCH = 116444736000000000
CONTROLLER_CLOCK_SKEW_S = 0.319
MOTOR_FULL_STEPS_PER_REV = 200
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

AXIS_COLOR = '#333333'
ROW_COLOR = '#1f77b4'
DRIVE_COLOR = '#d62728'
ZERO_COLOR = '#9a9a9a'


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


def cruise_zoom_data(run_index: int, rate: str) -> dict:
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import (
        N_Q, N_STATES, PORTS, LuGreModelRev42, _port_values, load_parameters,
        lugre_terms,
    )
    sys.path.insert(0, str(SCRIPTS_DIR))
    from command_reconstruction import (
        reconstruct_segments, trapezoid_fraction, trapezoid_fraction_array,
    )

    block_name = f'D_{rate}'
    # Only controller-event time differences are used below.  Keeping Unix
    # timestamps directly is equivalent to subtracting the IDS epoch/skew,
    # because that common offset cancels in every block-relative time.
    rows = load_log_rows(LOG_PATH, 0.0)
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
    lead = model.p['L'] / (2.0 * np.pi)
    t_d, n_r = model.p['T_d'], model.p['N_r']
    k_em = model.p['N_r'] * model.p['T_hold']
    way_values = _port_values(model.p, 'way')
    nut_values = _port_values(model.p, 'nut')
    sb_values = _port_values(model.p, 'sb')

    first_trapezoid = next(s for s in segments if s['kind'] == 'trapezoid')
    # reconstruct_segments rescales the theoretical trapezoid to the actual
    # logged move interval.  Use that reconstructed cruise speed: the nominal
    # D-rate label is not necessarily the speed represented by this timeline.
    cruise_rotor_rev_s = abs(
        first_trapezoid['value_end'] - first_trapezoid['value_start']
    ) / (first_trapezoid['duration_s'] + first_trapezoid['t_accel_s'])
    detent_hz = cruise_rotor_rev_s * MOTOR_FULL_STEPS_PER_REV
    window_s = N_CYCLES / detent_hz
    win_t0 = WINDOW_START_FROM_BLOCK_S
    win_t1 = win_t0 + window_s
    cruise_t0 = first_trapezoid['t0'] + first_trapezoid['t_accel_s']
    cruise_t1 = cruise_t0 + first_trapezoid['duration_s']
    if win_t0 < cruise_t0 or win_t1 > cruise_t1:
        raise RuntimeError(
            f'{block_name}: requested {win_t0:.3f}..{win_t1:.3f}s window '
            f'is outside cruise {cruise_t0:.3f}..{cruise_t1:.3f}s'
        )

    v_cmd = float(rate) * model.p['L'] / MOTOR_FULL_STEPS_PER_REV  # m/s, exact cruise value

    state = np.zeros(N_STATES)
    win_t, win_theta_m, win_vel, win_z = [], [], [], []
    for seg in segments:
        t0, t1 = seg['t0'], seg['t1']
        if t1 <= t0:
            continue
        # Figure B only uses the first trapezoid through the end of the
        # displayed window.  Integration after that instant cannot affect
        # earlier states, so do not spend several minutes solving the unseen
        # remainder of the outbound move, return move, and final hold.
        if seg is first_trapezoid:
            t1 = min(t1, win_t1)
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
            mask = (sol.t >= win_t0) & (sol.t <= win_t1)
            win_t.append(sol.t[mask])
            win_theta_m.append(sol.y[0, mask])
            win_vel.append(sol.y[N_Q:2 * N_Q, mask])
            win_z.append(sol.y[2 * N_Q:2 * N_Q + len(PORTS), mask])
            break

    t = np.concatenate(win_t)
    theta_m = np.concatenate(win_theta_m)
    vel = np.concatenate(win_vel, axis=1)
    z = np.concatenate(win_z, axis=1)  # rows: z_way, z_nut, z_sb (PORTS order)
    v_sim = vel[5, :]  # x_n_dot

    # theta_cmd isn't part of the ODE state -- it only exists transiently
    # inside the rhs() closure above -- so it's recomputed here from the
    # same trapezoid_fraction the solver used, restricted to the window
    # (which stays entirely within first_trapezoid by construction: the
    # window is a few cycles centered mid-cruise).
    frac = trapezoid_fraction_array(
        t - first_trapezoid['t0'], first_trapezoid['t_accel_s'], first_trapezoid['duration_s'],
    )
    theta_cmd = (
        first_trapezoid['value_start']
        + (first_trapezoid['value_end'] - first_trapezoid['value_start']) * frac
    ) * 2.0 * np.pi
    theta_err_deg = np.degrees(theta_cmd - theta_m)
    motor_torque = k_em * (theta_cmd - theta_m)  # T_motor = k_em*(theta_cmd - theta_m), k_em = N_r*T_hold
    motor_force = motor_torque / lead

    detent_torque = t_d * np.sin(4.0 * n_r * theta_m)
    # The rotor equation contains ``... + T_motor - T_detent``.  Plot the
    # signed contribution actually applied to that equation, rather than the
    # positive constitutive detent law, so this figure reads as a force
    # balance without requiring the viewer to infer the leading minus sign.
    detent_force = -detent_torque / lead

    z_way, z_nut, z_sb = z[PORTS.index('way')], z[PORTS.index('nut')], z[PORTS.index('sb')]
    v_way = model.jacobians['way'] @ vel
    v_nut = model.jacobians['nut'] @ vel
    v_sb = model.jacobians['sb'] @ vel
    friction_way = lugre_terms(v_way, z_way, *way_values)[0]
    friction_nut = lugre_terms(v_nut, z_nut, *nut_values)[0]
    friction_sb_torque_arr = lugre_terms(v_sb, z_sb, *sb_values)[0]
    friction_sb_force = friction_sb_torque_arr / lead

    return {
        'rate': rate, 't': t, 'v_cmd': v_cmd, 'v_sim': v_sim,
        'motor_force': motor_force, 'theta_err_deg': theta_err_deg,
        'detent_force': detent_force, 'detent_hz': detent_hz,
        'friction_way': friction_way, 'friction_nut': friction_nut,
        'friction_sb_force': friction_sb_force,
        'window_s': window_s, 'window_start_s': win_t0, 'lead': lead,
    }


def style_row(ax, ylabel):
    ax.grid(True, alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.axhline(0.0, color=ZERO_COLOR, lw=0.8, zorder=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    ax.set_ylabel(ylabel, fontsize=7.5)
    ax.tick_params(labelsize=7, colors=AXIS_COLOR)


def plot_column(fig, axes_col, result):
    t_ms = (result['t'] - result['t'][0]) * 1000.0

    axes_col[0].plot(t_ms, result['motor_force'], color=DRIVE_COLOR, lw=1.0, zorder=2)
    style_row(axes_col[0], 'Motor drive torque,\nload-side equiv. (N)')

    axes_col[1].plot(t_ms, result['v_sim'] * 1000.0, color=ROW_COLOR, lw=0.9, zorder=2)
    style_row(axes_col[1], 'Simulated stage\nvelocity v_sim (mm/s)')

    axes_col[2].plot(t_ms, result['theta_err_deg'], color=ROW_COLOR, lw=1.0, zorder=2)
    style_row(axes_col[2], 'Tracking error\ntheta_cmd - theta_m (deg)')

    axes_col[3].plot(t_ms, result['detent_force'], color=ROW_COLOR, lw=1.0, zorder=2)
    style_row(axes_col[3], 'Detent contribution\n(-T_detent / lead) (N)')

    axes_col[4].plot(t_ms, result['friction_way'], color=ROW_COLOR, lw=1.0, zorder=2)
    style_row(axes_col[4], 'Guideway friction\n(way), load-side (N)')

    axes_col[5].plot(t_ms, result['friction_nut'], color=ROW_COLOR, lw=1.0, zorder=2)
    style_row(axes_col[5], 'Leadscrew-nut friction\n(nut), load-side (N)')

    axes_col[6].plot(t_ms, result['friction_sb_force'], color=ROW_COLOR, lw=1.0, zorder=2)
    style_row(axes_col[6], 'Bearing friction (sb),\nload-side equiv. (N)')

    axes_col[0].set_title(
        f"D {result['rate']} full-steps/s\n"
        f"window = {result['window_s']*1000:.1f} ms ({N_CYCLES} cycles), "
        f"starts at {result['window_start_s']:.1f} s",
        fontsize=9.5,
    )
    axes_col[-1].set_xlabel('Time within window (ms)', fontsize=8)


def main() -> None:
    print(f'Figure B: cruise-window zoom, run {RUN_INDEX}, rates {RATES}', flush=True)
    started = time.perf_counter()
    results = {}
    with ProcessPoolExecutor(max_workers=len(RATES)) as pool:
        futures = {
            pool.submit(cruise_zoom_data, RUN_INDEX, rate): rate for rate in RATES
        }
        for future in as_completed(futures):
            rate = futures[future]
            results[rate] = future.result()
            print(f'  done D_{rate} ({time.perf_counter()-started:.1f}s)', flush=True)

    n_rows = 7
    fig, axes = plt.subplots(n_rows, len(RATES), figsize=(7.5 * len(RATES), 14.5), squeeze=False)
    for col, rate in enumerate(RATES):
        plot_column(fig, [axes[row][col] for row in range(n_rows)], results[rate])

    fig.suptitle(
        'Figure B -- Run 2 (MRES 1/4, 100% I): cruise-window force balance\n'
        'Motor torque converted with F = T / lead; detent shown as its '
        'signed opposing contribution.',
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))

    out_path = OUT_DIR / 'figureB_cruise_zoom_D9.5.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}', flush=True)
    print(f'Total wall time: {time.perf_counter()-started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

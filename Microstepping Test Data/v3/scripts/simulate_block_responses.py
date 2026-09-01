#!/usr/bin/env python3
"""Simulate the Rev 4.2 parallel-LuGre + nonlinear-detent model's response
to the ACTUAL recorded commanded trajectory for every (run, block) pair,
in parallel across CPU cores. Each block is simulated independently from
rest (matching how the measured/commanded montage panels are already
baselined per-block).

Saves one npz per run into rendered_assets/individual_subplots/run_0N_.../
containing every block's simulated (time_s, position_um) arrays.
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

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = ROOT / 'data' / 'raw_local'
SUBPLOT_DIR = ROOT / 'rendered_assets' / 'individual_subplots'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")

FILETIME_UNIX_EPOCH = 116444736000000000
BLOCK_NAMES = (
    'BLOCK_0_START', 'BLOCK_0_END', 'C',
    'D_0.125', 'D_0.375', 'D_1.25', 'D_3.5', 'D_9.5', 'D_27.5', 'D_70', 'D_200',
)
RTOL = 1.0e-6
ATOL = np.array([
    1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-12, 1.0e-12,
    1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-9, 1.0e-9,
    1.0e-11, 1.0e-11, 1.0e-9,
])
MODEL_LABEL = 'Rev 4.2 parallel LuGre + nonlinear detent torque'


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
    start_epoch_s = (start_filetime - FILETIME_UNIX_EPOCH) / 1.0e7
    return start_epoch_s


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


def real_command_segments(rows, run_index, block_start_s, block_end_s):
    """theta_command in RADIANS, piecewise-constant between real MOVE_ACK
    events (ideal_position_rev is in revolutions -> x2*pi for radians,
    matching command_series()'s ideal_position_rev*LEAD_PITCH_M for um)."""
    moves = [
        r for r in rows if r['event'] == 'MOVE_ACK'
        and r['run_index'] == str(run_index) and r['ideal_position_rev']
        and block_start_s <= r['ids_time_s'] <= block_end_s
    ]
    times = [0.0] + [m['ids_time_s'] - block_start_s for m in moves]
    theta = [0.0] + [float(m['ideal_position_rev']) * 2.0 * np.pi for m in moves]
    times.append(block_end_s - block_start_s)
    theta.append(theta[-1])
    return times, theta


def simulate_one_block(run_index: int, block_name: str) -> dict:
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import N_STATES, LuGreModelRev42

    start_epoch_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    start_s, end_s = find_block(rows, run_index, block_name)
    times, theta = real_command_segments(rows, run_index, start_s, end_s)

    model = LuGreModelRev42(enforce_interface_power=False)
    lead = model.p['L'] / (2.0 * np.pi)

    state = np.zeros(N_STATES)
    started = time.perf_counter()
    t_parts, x_parts = [], []
    for i in range(len(times) - 1):
        t0, t1 = times[i], times[i + 1]
        if t1 <= t0:
            continue
        theta_command = theta[i]
        sol = solve_ivp(
            lambda t, y, cmd=theta_command: model.rhs(t, y, cmd),
            (t0, t1), state, method='Radau',
            jac=lambda t, y: model.analytical_linearization(y)[0],
            t_eval=np.linspace(t0, t1, max(2, int((t1 - t0) * 200))),
            rtol=RTOL, atol=ATOL,
        )
        if not sol.success:
            raise RuntimeError(
                f'run={run_index} block={block_name} segment {i} failed: '
                f'{sol.message}'
            )
        state = sol.y[:, -1]
        t_parts.append(sol.t)
        x_parts.append(sol.y[5, :])  # x_n (nut/stage position, meters)
    elapsed = time.perf_counter() - started
    time_s = np.concatenate(t_parts)
    position_um = np.concatenate(x_parts) * 1.0e6
    return {
        'run_index': run_index, 'block_name': block_name,
        'time_s': time_s, 'position_um': position_um,
        'elapsed_s': elapsed, 'real_duration_s': end_s - start_s,
    }


def main() -> None:
    jobs = [
        (run_index, block_name)
        for run_index in range(1, 7)
        for block_name in BLOCK_NAMES
    ]
    print(f'{len(jobs)} jobs to run, model={MODEL_LABEL}', flush=True)
    workers = max(1, min(18, (os.cpu_count() or 4) - 2))
    print(f'Using {workers} worker processes', flush=True)

    results_by_run: dict[int, dict[str, dict]] = {i: {} for i in range(1, 7)}
    started_all = time.perf_counter()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(simulate_one_block, run_index, block_name): (run_index, block_name)
            for run_index, block_name in jobs
        }
        for future in as_completed(futures):
            run_index, block_name = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                print(f'FAILED run={run_index} block={block_name}: {exc}', flush=True)
                continue
            results_by_run[run_index][block_name] = result
            done += 1
            print(
                f'[{done}/{len(jobs)}] run={run_index} block={block_name} '
                f'elapsed={result["elapsed_s"]:.1f}s '
                f'real_duration={result["real_duration_s"]:.1f}s '
                f'(total wall {time.perf_counter()-started_all:.1f}s)',
                flush=True,
            )

    for run_index, blocks in results_by_run.items():
        if not blocks:
            continue
        mres = None
        current = None
        # recover mres/current for filename consistency with plot_block_montage.py
        start_epoch_s = parse_ids(IDS_PATH)
        rows = load_log_rows(LOG_PATH, start_epoch_s)
        start_row = next(
            r for r in rows if r['event'] == 'RUN_CONFIG'
            and r['run_index'] == str(run_index)
        )
        mres, current = start_row['mres'], start_row['current']
        folder_name = f'run_{run_index:02d}_mres_{mres}_{current.lower()}'
        out_dir = SUBPLOT_DIR / folder_name
        out_dir.mkdir(parents=True, exist_ok=True)
        payload = {'model_label': np.asarray(MODEL_LABEL)}
        for block_name, result in blocks.items():
            payload[f'{block_name}_time_s'] = result['time_s']
            payload[f'{block_name}_position_um'] = result['position_um']
        np.savez_compressed(out_dir / 'lugre_simulation.npz', **payload)
        print(f'Saved {out_dir / "lugre_simulation.npz"}', flush=True)

    print(f'Total wall time: {time.perf_counter()-started_all:.1f}s', flush=True)


if __name__ == '__main__':
    main()

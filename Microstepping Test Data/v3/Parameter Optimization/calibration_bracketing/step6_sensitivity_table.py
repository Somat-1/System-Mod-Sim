#!/usr/bin/env python3
"""Step 6: sensitivity table for the guideway ("way") LuGre port, run 2
only. Six parameters (sigma0, sigma1, sigma2, Fc, Fs, vs) x eight D-rate
blocks (D_0.125 ... D_200), each perturbed +-10% independently, one at a
time, with cost = RMS(measured - simulated) over the whole block (µm).

Sensitivity per cell = 0.5*(|cost(+10%)-cost(base)| + |cost(-10%)-cost(base)|)
  / cost(base)   -- a normalized, symmetric relative sensitivity.

104 total simulations (8 baseline + 6*2*8 perturbed), parallelized across
worker processes. Some individual D-rate blocks (D_0.375, D_1.25, D_27.5)
take several minutes each on their own -- this is a genuinely long run.
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.integrate import solve_ivp

STATUS_PATH = Path(__file__).resolve().parent / 'sensitivity_status.json'

V3_ROOT = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3")
RAW_DIR = V3_ROOT / 'data' / 'raw_local'
IDS_PATH = RAW_DIR / 'SteppingSequenceID.csv'
LOG_PATH = RAW_DIR / 'identification_controller_log.csv'
REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
SCRIPTS_DIR = V3_ROOT / 'scripts'
OUT_DIR = Path(__file__).resolve().parent

RUN_INDEX = 2
D_RATES = ('0.125', '0.375', '1.25', '3.5', '9.5', '27.5', '70', '200')
PARAMS = ('sigma0_way', 'sigma1_way', 'sigma2_way', 'Fc_way', 'Fs_way', 'vs_way')
PERTURBATION = 0.10
NOISE_FLOOR = 0.01  # 1% normalized sensitivity -- below this, treat as "no effect"

FILETIME_UNIX_EPOCH = 116444736000000000
CONTROLLER_CLOCK_SKEW_S = 0.319
LEAD_PITCH_M = 2.0e-3
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


def simulate_block(run_index: int, block_name: str, param_overrides: dict, rows) -> tuple:
    """Simulate one block from rest with param_overrides applied on top of
    the current best parameters (plus the usual per-run current-dependent
    T_hold). Returns (time_s, position_um). `rows` is the already-parsed
    log (shared with run_job's measured-side lookup, to avoid re-parsing
    the log CSV twice per job)."""
    sys.path.insert(0, str(REV42_SCRIPTS))
    from lugre_model_rev42 import N_STATES, LuGreModelRev42, load_parameters
    sys.path.insert(0, str(SCRIPTS_DIR))
    from command_reconstruction import reconstruct_segments, trapezoid_fraction

    start_s, end_s = find_block(rows, run_index, block_name)
    segments = reconstruct_segments(rows, run_index, start_s, end_s, block_name)

    run_config = next(
        r for r in rows if r['event'] == 'RUN_CONFIG'
        and r['run_index'] == str(run_index)
    )
    peak_ma = PEAK_MA_BY_CURRENT_LEVEL[run_config['current']]
    t_hold_run = K_T_NM_PER_A * np.sqrt(2.0) * (peak_ma / 1000.0)
    params = load_parameters()
    params['T_hold'] = t_hold_run
    params.update(param_overrides)
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
            raise RuntimeError(
                f'run={run_index} block={block_name} overrides={param_overrides} '
                f'segment failed: {sol.message}'
            )
        state = sol.y[:, -1]
        t_parts.append(sol.t)
        x_parts.append(sol.y[5, :])
    return np.concatenate(t_parts), np.concatenate(x_parts) * 1.0e6


_CACHED = None


def _worker_init(time_s, position_nm, start_epoch_s, sample_period_s, rows):
    """Runs once per worker process at pool startup. parse_ids() on the
    2.78M-row IDS CSV costs ~50s+ of network I/O over the UNC share;
    calling it fresh inside every one of the 104 jobs (as before) meant up
    to 18 worker processes hammering the same network file concurrently,
    which is what produced the wildly inconsistent per-job timings observed
    earlier (23s vs 260s+ for the same block) -- not parameter-dependent
    slowdown or network drift. Parsing once in the main process and handing
    the arrays to each worker via initargs removes that redundant I/O
    entirely."""
    global _CACHED
    _CACHED = {
        'time_s': time_s, 'position_nm': position_nm,
        'sample_period_s': sample_period_s, 'rows': rows,
    }


def run_job(run_index, block_name, param_overrides):
    """Full job: simulate + compute cost = RMS(measured - simulated) over
    the block, in um. Uses the per-worker cached IDS/log data set up once
    by _worker_init instead of re-parsing per job."""
    time_s = _CACHED['time_s']
    position_nm = _CACHED['position_nm']
    sample_period_s = _CACHED['sample_period_s']
    rows = _CACHED['rows']

    sim_t, sim_y = simulate_block(run_index, block_name, param_overrides, rows)

    start_s, end_s = find_block(rows, run_index, block_name)
    t_meas, y_meas = measured_window(time_s, position_nm, sample_period_s, start_s, end_s)

    sim_on_meas_grid = np.interp(t_meas, sim_t, sim_y, left=sim_y[0], right=sim_y[-1])
    cost = float(np.sqrt(np.mean((y_meas - sim_on_meas_grid) ** 2)))
    return cost


def _key_str(key):
    tag, param, block = key
    return f'{tag}|{param}|{block}'


def write_status(jobs, job_state, started, workers, note=None):
    """Single-writer status file (only the main process calls this, never
    the worker processes) -- safe to poll from outside while the run is in
    progress. Per-job timing on this network share has been observed to
    vary by an order of magnitude run-to-run, so estimated_remaining_s is
    a rough, low-confidence projection, not a real ETA."""
    now = datetime.now(timezone.utc)
    done_states = [s for s in job_state.values() if s['status'] == 'done']
    failed_states = [s for s in job_state.values() if s['status'] == 'failed']
    timed_done_states = [s for s in done_states if 'job_wall_s' in s]
    elapsed_s = time.perf_counter() - started
    avg_job_s = (
        sum(s['job_wall_s'] for s in timed_done_states) / len(timed_done_states)
        if timed_done_states else None
    )
    finished = len(done_states) + len(failed_states)
    remaining = len(jobs) - finished
    # Throughput-based estimate (elapsed / jobs-finished-so-far), not
    # avg_job_s * remaining/workers: with 104 jobs sharing ~18 workers,
    # individual job_wall_s includes queueing delay and is a poor basis
    # for an ETA on its own; wall-clock-per-completed-slot already nets
    # out the parallelism actually being achieved.
    est_remaining_s = (
        (elapsed_s / finished) * remaining if finished > 0 else None
    )
    payload = {
        'schema': 'step6-sensitivity-status-v1',
        'note': note or (
            'estimated_remaining_s is a rough projection from the average '
            'completed-job time; per-job time on this network share has '
            'varied by 10x+ run-to-run this session, so treat it as a very '
            'loose bound, not a real ETA.'
        ),
        'started_at': datetime.fromtimestamp(
            now.timestamp() - elapsed_s, tz=timezone.utc
        ).isoformat(),
        'last_updated': now.isoformat(),
        'workers': workers,
        'total_jobs': len(jobs),
        'completed_jobs': len(done_states),
        'failed_jobs': len(failed_states),
        'pending_jobs': remaining,
        'elapsed_s': round(elapsed_s, 1),
        'avg_completed_job_wall_s': round(avg_job_s, 1) if avg_job_s else None,
        'estimated_remaining_s': round(est_remaining_s, 1) if est_remaining_s else None,
        'jobs': {_key_str(k): v for k, v in job_state.items()},
    }
    tmp_path = STATUS_PATH.with_suffix('.json.tmp')
    tmp_path.write_text(json.dumps(payload, indent=2))
    tmp_path.replace(STATUS_PATH)  # atomic-ish swap so readers never see a half-written file


def main():
    jobs = {}
    for block_name in D_RATES:
        jobs[('baseline', None, block_name)] = (RUN_INDEX, f'D_{block_name}', {})
    for param in PARAMS:
        sys.path.insert(0, str(REV42_SCRIPTS))
        from lugre_model_rev42 import load_parameters
        base_value = load_parameters()[param]
        for sign, tag in ((+1, 'plus'), (-1, 'minus')):
            new_value = base_value * (1.0 + sign * PERTURBATION)
            for block_name in D_RATES:
                jobs[(tag, param, block_name)] = (
                    RUN_INDEX, f'D_{block_name}', {param: new_value}
                )

    print(f'{len(jobs)} total jobs (8 baseline + {len(PARAMS)}*2*8 perturbed)', flush=True)
    # Two prior launches at 18 workers both died with an unexplained
    # mid-run process-group crash (exit 127, no Python traceback, jobs
    # were completing with valid results right up to the moment it died)
    # at inconsistent points (~200s, ~430s). The machine was already at
    # 85% RAM usage / ~9.6GB free before launch. Kept deliberately
    # conservative rather than re-diagnosing further -- trades wall time
    # for not tipping an already memory-pressured machine over again.
    workers = max(1, min(6, (os.cpu_count() or 4) - 2))
    print(f'Using {workers} worker processes', flush=True)
    print(f'Live status: {STATUS_PATH}', flush=True)

    results = {}
    if STATUS_PATH.exists():
        try:
            prior = json.loads(STATUS_PATH.read_text())
            for key_str, state in prior.get('jobs', {}).items():
                if state.get('status') == 'done' and 'cost_um' in state:
                    tag, param, block = key_str.split('|')
                    param = None if param == 'None' else param
                    key = (tag, param, block)
                    if key in jobs:
                        results[key] = state['cost_um']
            if results:
                print(f'Resuming: {len(results)}/{len(jobs)} jobs already done '
                      f'per prior {STATUS_PATH.name}, skipping those.', flush=True)
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            print(f'Could not parse prior status file for resume ({exc}); '
                  f'starting fresh.', flush=True)

    print('Parsing IDS + log once in the main process (shared with all '
          'workers via pool initializer, avoids re-parsing per job)...', flush=True)
    parse_t0 = time.perf_counter()
    time_s, position_nm, start_epoch_s, sample_period_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    print(f'  done in {time.perf_counter() - parse_t0:.1f}s '
          f'({len(rows)} log rows, {position_nm.size} IDS samples)', flush=True)

    started = time.perf_counter()
    job_state = {
        key: ({'status': 'done', 'cost_um': results[key]} if key in results
              else {'status': 'pending'})
        for key in jobs
    }
    write_status(jobs, job_state, started, workers)
    done = len(results)
    remaining_jobs = {k: v for k, v in jobs.items() if k not in results}
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_worker_init,
        initargs=(time_s, position_nm, start_epoch_s, sample_period_s, rows),
    ) as pool:
        submit_times = {}
        futures = {}
        for key, (run_index, block_name, overrides) in remaining_jobs.items():
            futures[pool.submit(run_job, run_index, block_name, overrides)] = key
            submit_times[key] = time.perf_counter()
            job_state[key]['status'] = 'running'
        write_status(jobs, job_state, started, workers)
        for future in as_completed(futures):
            key = futures[future]
            job_wall_s = time.perf_counter() - submit_times[key]
            try:
                cost = future.result()
                job_state[key] = {
                    'status': 'done', 'cost_um': round(cost, 4),
                    'job_wall_s': round(job_wall_s, 1),
                }
            except Exception as exc:
                print(f'FAILED {key}: {exc}', flush=True)
                cost = float('nan')
                job_state[key] = {
                    'status': 'failed', 'error': str(exc),
                    'job_wall_s': round(job_wall_s, 1),
                }
            results[key] = cost
            done += 1
            print(f'[{done}/{len(jobs)}] {key} cost={cost:.4f} um '
                  f'(job {job_wall_s:.1f}s, wall {time.perf_counter()-started:.1f}s)', flush=True)
            write_status(jobs, job_state, started, workers)

    # --- build the table ---
    baseline_cost = {block: results[('baseline', None, block)] for block in D_RATES}
    print('\nBaseline cost per block (RMS measured-vs-simulated, um):')
    for block in D_RATES:
        print(f'  D_{block:>6s}: {baseline_cost[block]:.4f}')

    sensitivity = {}
    for param in PARAMS:
        row = []
        for block in D_RATES:
            base = baseline_cost[block]
            c_plus = results[('plus', param, block)]
            c_minus = results[('minus', param, block)]
            if base > 0 and np.isfinite(c_plus) and np.isfinite(c_minus):
                s = 0.5 * (abs(c_plus - base) + abs(c_minus - base)) / base
            else:
                s = float('nan')
            row.append(s)
        sensitivity[param] = row

    print(f'\nNormalized sensitivity table (0.5*(|dcost+|+|dcost-|)/cost_base), '
          f'noise floor = {NOISE_FLOOR:.2f}:')
    header = 'parameter'.ljust(12) + ''.join(f'D_{b:>8s}' for b in D_RATES)
    print(header)
    for param in PARAMS:
        row_str = param.ljust(12) + ''.join(
            f'{v:>10.4f}' if np.isfinite(v) else f'{"nan":>10s}' for v in sensitivity[param]
        )
        above = [v > NOISE_FLOOR for v in sensitivity[param] if np.isfinite(v)]
        verdict = 'FIT' if any(above) else 'FIX (no column above noise floor)'
        print(f'{row_str}   -> {verdict}')

    # collinearity check: Fs_way vs sigma0_way column patterns
    fs_row = np.array(sensitivity['Fs_way'])
    s0_row = np.array(sensitivity['sigma0_way'])
    valid = np.isfinite(fs_row) & np.isfinite(s0_row)
    if np.any(valid) and np.std(fs_row[valid]) > 0 and np.std(s0_row[valid]) > 0:
        corr = np.corrcoef(fs_row[valid], s0_row[valid])[0, 1]
        print(f'\nFs_way vs sigma0_way column-pattern correlation: {corr:+.4f} '
              f'({"collinear -- ratio constraint needed" if corr > 0.9 else "not obviously collinear"})')

    # save raw results
    out_json = OUT_DIR / 'step6_sensitivity_table.json'
    serializable = {f'{k[0]}|{k[1]}|{k[2]}': v for k, v in results.items()}
    out_json.write_text(json.dumps(serializable, indent=2))
    print(f'\nSaved raw results to {out_json}')
    print(f'Total wall time: {time.perf_counter()-started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

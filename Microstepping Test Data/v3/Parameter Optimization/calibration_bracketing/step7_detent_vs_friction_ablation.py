#!/usr/bin/env python3
"""Step 7: does the way-port friction actually shape the simulated
trajectory, or is it a small perturbation riding on a detent-torque-
dominated response?

Directly tests this by ablating each term from the Rev 4.2 model and
comparing the resulting trajectory to the *full* model's trajectory (not
to measurement -- that comparison is what step6's sensitivity table did,
and it's confounded by the large baseline model-vs-measurement mismatch at
D_70/D_200). Three variants per block, run 2:

  full             -- unmodified model
  no_detent        -- T_d = 0 (detent torque removed)
  no_way_friction  -- sigma0_way = sigma1_way = sigma2_way = 0 (way LuGre
                      force forced to exactly zero; safe because it removes
                      sigma0 from the z_dot decay term too, no risk of
                      dividing by a zero Fc/Fs floor)

For each block: RMS(full - no_detent) and RMS(full - no_way_friction) are
the trajectory-shape contributions of each term, directly comparable to
each other. Also reports each variant's own RMS-vs-measurement cost for
context against step6's baseline numbers.
"""
from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step6_sensitivity_table as s6

RUN_INDEX = s6.RUN_INDEX
BLOCKS = ('0.125', '3.5', '70', '200')  # slow anchor, peak Fc/Fs/vs-sensitivity block, and the two blocks with the large baseline model-vs-measurement mismatch
VARIANTS = {
    'full': {},
    'no_detent': {'T_d': 0.0},
    'no_way_friction': {'sigma0_way': 0.0, 'sigma1_way': 0.0, 'sigma2_way': 0.0},
}


def run_variant(run_index, block_name, overrides):
    time_s = s6._CACHED['time_s']
    position_nm = s6._CACHED['position_nm']
    sample_period_s = s6._CACHED['sample_period_s']
    rows = s6._CACHED['rows']

    sim_t, sim_y = s6.simulate_block(run_index, block_name, overrides, rows)
    start_s, end_s = s6.find_block(rows, run_index, block_name)
    t_meas, y_meas = s6.measured_window(time_s, position_nm, sample_period_s, start_s, end_s)
    sim_on_meas_grid = np.interp(t_meas, sim_t, sim_y, left=sim_y[0], right=sim_y[-1])
    cost_vs_measured = float(np.sqrt(np.mean((y_meas - sim_on_meas_grid) ** 2)))
    return sim_t, sim_y, cost_vs_measured


def main():
    print('Parsing IDS + log once...', flush=True)
    t0 = time.perf_counter()
    time_s, position_nm, start_epoch_s, sample_period_s = s6.parse_ids(s6.IDS_PATH)
    rows = s6.load_log_rows(s6.LOG_PATH, start_epoch_s)
    print(f'  done in {time.perf_counter() - t0:.1f}s', flush=True)

    jobs = {
        (block, variant): (RUN_INDEX, f'D_{block}', overrides)
        for block in BLOCKS for variant, overrides in VARIANTS.items()
    }
    workers = min(len(jobs), 12)
    print(f'{len(jobs)} jobs, {workers} workers', flush=True)

    results = {}
    started = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=workers, initializer=s6._worker_init,
        initargs=(time_s, position_nm, start_epoch_s, sample_period_s, rows),
    ) as pool:
        futures = {pool.submit(run_variant, *v): k for k, v in jobs.items()}
        for fut in as_completed(futures):
            key = futures[fut]
            sim_t, sim_y, cost = fut.result()
            results[key] = (sim_t, sim_y, cost)
            print(f'{key}: cost_vs_measured={cost:.4f} um '
                  f'(t={time.perf_counter() - started:.1f}s)', flush=True)

    print('\n=== Trajectory-shape contribution: RMS(full - ablated), um ===')
    print(f'{"block":>8s}  {"no_detent":>12s}  {"no_way_friction":>16s}  {"dominant term":>14s}')
    for block in BLOCKS:
        t_full, y_full, cost_full = results[(block, 'full')]
        _, y_nodetent, cost_nodetent = results[(block, 'no_detent')]
        _, y_nofric, cost_nofric = results[(block, 'no_way_friction')]
        rms_detent = float(np.sqrt(np.mean((y_full - y_nodetent) ** 2)))
        rms_fric = float(np.sqrt(np.mean((y_full - y_nofric) ** 2)))
        dominant = 'detent' if rms_detent > rms_fric else 'friction'
        print(f'D_{block:>6s}  {rms_detent:>12.4f}  {rms_fric:>16.4f}  {dominant:>14s}')

    print('\n=== cost vs measurement per variant, um (context vs step6 baseline) ===')
    print(f'{"block":>8s}  {"full":>10s}  {"no_detent":>12s}  {"no_way_friction":>16s}')
    for block in BLOCKS:
        _, _, cost_full = results[(block, 'full')]
        _, _, cost_nodetent = results[(block, 'no_detent')]
        _, _, cost_nofric = results[(block, 'no_way_friction')]
        print(f'D_{block:>6s}  {cost_full:>10.4f}  {cost_nodetent:>12.4f}  {cost_nofric:>16.4f}')

    print(f'\nTotal wall time: {time.perf_counter() - started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

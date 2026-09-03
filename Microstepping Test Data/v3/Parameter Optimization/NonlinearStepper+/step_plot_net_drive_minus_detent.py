#!/usr/bin/env python3
"""Informational plot: motor drive torque minus the detent torque acting
against it -- theta_m's local torque budget before it has to fight the
rest of the drivetrain (k_c*(theta_m-theta_c), c_em*theta_m_dot, etc.):

  I_m*theta_m_ddot = (T_motor - T_detent) - k_c*(theta_m-theta_c)
                      - c_c*(theta_m_dot-theta_c_dot) - c_EM*theta_m_dot

Net = T_motor - T_detent is literally that grouped term -- "what remains
after the motor and the detent have fought it out, and gets fed into the
coupling/shaft/screw chain."

Uses Option A (the nonlinear drive law) since linear and Option A were
already shown to be visually and numerically indistinguishable at these
blocks (see compare_linear_vs_optionA.json / README.md) -- no need to
duplicate both here.
"""
from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import step_compare_linear_vs_optionA as base

BLOCKS = base.TRAPEZOID_BLOCKS  # ('3.5', '70', '200')
MODEL_KIND = 'optionA'


def main():
    print('Parsing IDS + log once...', flush=True)
    t0 = time.perf_counter()
    time_s, position_nm, start_epoch_s, sample_period_s = base.parse_ids(base.IDS_PATH)
    rows = base.load_log_rows(base.LOG_PATH, start_epoch_s)
    print(f'  done in {time.perf_counter() - t0:.1f}s', flush=True)

    started = time.perf_counter()
    results = {}
    with ProcessPoolExecutor(
        max_workers=len(BLOCKS), initializer=base._worker_init,
        initargs=(time_s, position_nm, start_epoch_s, sample_period_s, rows),
    ) as pool:
        futures = {
            pool.submit(base.simulate_variant, base.RUN_INDEX, f'D_{block}', MODEL_KIND): block
            for block in BLOCKS
        }
        for fut in as_completed(futures):
            block = futures[fut]
            results[block] = fut.result()
            print(f'D_{block} done (t={time.perf_counter() - started:.1f}s)', flush=True)

    fig, axes = plt.subplots(1, len(BLOCKS), figsize=(7.5 * len(BLOCKS), 5.5))
    for ax, block in zip(axes, BLOCKS):
        res = results[block]
        t_ms = (res['window_t'] - res['window_t'][0]) * 1000.0
        ax.plot(t_ms, res['window_motor_torque_Nm'], color='#1f77b4', lw=1.0,
                alpha=0.75, label='Motor drive torque')
        ax.plot(t_ms, -res['window_detent_Nm'], color='#9467bd', lw=1.0,
                alpha=0.75, label='-Detent (opposing)')
        ax.plot(t_ms, res['window_net_drive_Nm'], color='#d62728', lw=1.6,
                label='Net = motor - detent\n(fed into the drivetrain)')
        ax.axhline(0.0, color='#9a9a9a', lw=0.7)
        ax.set_title(f'D {block} full-steps/s\n({res["detent_hz"]:.1f} Hz detent forcing)', fontsize=10)
        ax.set_xlabel('Time within window (ms)', fontsize=8.5)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7.5)
        if block == BLOCKS[0]:
            ax.set_ylabel('Torque (N·m)', fontsize=9)
            ax.legend(loc='best', fontsize=8, framealpha=0.9)

    fig.suptitle(
        'Motor drive torque net of detent opposition -- run 2, Option A drive law '
        '(linear is visually indistinguishable here, see README.md)',
        fontsize=11.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_path = HERE / 'net_drive_minus_detent.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}', flush=True)
    print(f'Total wall time: {time.perf_counter() - started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

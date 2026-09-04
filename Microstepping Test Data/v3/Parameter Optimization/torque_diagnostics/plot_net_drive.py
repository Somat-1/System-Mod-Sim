#!/usr/bin/env python3
"""Figure C: net drive -- motor torque with the detent contribution
already subtracted, i.e. what is actually left over to drive the stage
through the rest of the drivetrain (k_c/k_s1/k_s2/k_nut/k_brg chain).

Reuses cruise_zoom_data() from plot_cruise_zoom.py unmodified, for the
same run/rates/windowing (5 s into the block, 10 detent cycles, actual
reconstructed cruise speed -- see that module's docstring). detent_force
returned from there is already the signed -T_detent/lead contribution to
the rotor force balance, so:

    net_force = motor_force + detent_force  (both load-side, N)

is directly "T_motor - T_detent", converted to load-side force -- the
same grouped term that appears in the rotor's equation of motion,
I_m*theta_m_ddot = (T_motor - T_detent) - k_c*(theta_m - theta_c) - ...
"""
from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from plot_cruise_zoom import (  # noqa: E402
    RUN_INDEX, RATES, N_CYCLES, AXIS_COLOR, ZERO_COLOR, DRIVE_COLOR,
    cruise_zoom_data, style_row,
)

MOTOR_COLOR = DRIVE_COLOR
DETENT_COLOR = '#9467bd'
NET_COLOR = '#2ca02c'


def plot_column(axes_col, result):
    t_ms = (result['t'] - result['t'][0]) * 1000.0
    net_force = result['motor_force'] + result['detent_force']

    axes_col[0].plot(t_ms, result['motor_force'], color=MOTOR_COLOR, lw=1.0,
                      alpha=0.75, label='Motor drive, load-side (N)')
    axes_col[0].plot(t_ms, result['detent_force'], color=DETENT_COLOR, lw=1.0,
                      alpha=0.75, label='Detent contribution\n(-T_detent/lead), load-side (N)')
    axes_col[0].plot(t_ms, net_force, color=NET_COLOR, lw=1.6,
                      label='Net = motor + detent\n(fed into the drivetrain)')
    style_row(axes_col[0], 'Torque/force,\nload-side equiv. (N)')
    axes_col[0].legend(loc='best', fontsize=6.5, framealpha=0.9)

    axes_col[1].plot(t_ms, net_force, color=NET_COLOR, lw=1.4)
    style_row(axes_col[1], 'Net drive,\nload-side equiv. (N)')

    axes_col[0].set_title(
        f"D {result['rate']} full-steps/s\n"
        f"window = {result['window_s']*1000:.1f} ms ({N_CYCLES} cycles), "
        f"starts at {result['window_start_s']:.1f} s",
        fontsize=9.5,
    )
    axes_col[-1].set_xlabel('Time within window (ms)', fontsize=8)


def main() -> None:
    print(f'Figure C: net drive, run {RUN_INDEX}, rates {RATES}', flush=True)
    started = time.perf_counter()
    results = {}
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = {
            pool.submit(cruise_zoom_data, RUN_INDEX, rate): rate for rate in RATES
        }
        for future in as_completed(futures):
            rate = futures[future]
            results[rate] = future.result()
            print(f'  done D_{rate} ({time.perf_counter() - started:.1f}s)', flush=True)

    n_rows = 2
    fig, axes = plt.subplots(n_rows, len(RATES), figsize=(7.5 * len(RATES), 6.5), squeeze=False)
    for col, rate in enumerate(RATES):
        plot_column([axes[row][col] for row in range(n_rows)], results[rate])

    fig.suptitle(
        'Figure C -- Run 2 (MRES 1/4, 100% I): net drive torque, motor minus detent\n'
        '(motor_force + signed detent_force, both load-side N -- what is left over '
        'to actually drive the stage through the drivetrain)',
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.90))

    out_path = HERE / 'figureC_net_drive_D9.5.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}', flush=True)
    print(f'Total wall time: {time.perf_counter() - started:.1f}s', flush=True)


if __name__ == '__main__':
    main()

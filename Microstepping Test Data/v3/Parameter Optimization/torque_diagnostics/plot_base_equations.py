#!/usr/bin/env python3
"""Figure D: the two base torque equations actually used in
lugre_model_rev42.py (Rev 4.2, linear drive), plotted as pure analytical
curves -- no simulation, no time axis, just the constitutive laws vs.
their own natural variable.

  T_detent(theta_m) = T_d * sin(4*N_r*theta_m)              (lugre_model_rev42.py:216)
  T_motor(theta_err) = k_em * theta_err,  k_em = N_r*T_hold  (build_structural_matrices/rhs)
    theta_err = theta_cmd - theta_m

These are DIFFERENT functions of DIFFERENT variables (rotor position vs.
tracking error) -- they only interact through the dynamics of a specific
trajectory (see Figure B/C), not through their equations alone. This
figure is deliberately just the two laws side by side, for reference.

Detent's spatial period in theta_m is exactly one full step:
4*N_r*Delta_theta = 2*pi -> Delta_theta = pi/(2*N_r) = 1.8 deg (N_r=50) --
the same 1.8 deg that is the drive's own pull-out threshold in the
Option A nonlinear law (see ../NonlinearStepper+/README.md), plotted here
for comparison against the linear law actually in use.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]  # torque_diagnostics -> Parameter Optimization -> v3 -> Microstepping Test Data -> repo root
PARAMETER_FILE = REPO_ROOT / 'Rev 4' / 'lugre_friction' / 'Rev 4.2' / 'model_parameters.json'

T_HOLD_BASE_NM = 0.060
I_RATED_BASE_A = 0.400
K_T_NM_PER_A = T_HOLD_BASE_NM / (np.sqrt(2.0) * I_RATED_BASE_A)
CURRENTS_MA = {'50% I (200 mA)': 200.0, '100% I (400 mA)': 400.0}

DETENT_COLOR = '#9467bd'
MOTOR_COLORS = {'50% I (200 mA)': '#ff7f0e', '100% I (400 mA)': '#d62728'}
AXIS_COLOR = '#333333'
ZERO_COLOR = '#9a9a9a'


def main() -> None:
    params = json.loads(PARAMETER_FILE.read_text(encoding='utf-8'))['parameters']
    t_d, n_r = params['T_d'], params['N_r']
    pullout_deg = np.degrees(np.pi / (2.0 * n_r))  # one full step

    fig, (ax_detent, ax_motor) = plt.subplots(1, 2, figsize=(13.0, 5.5))

    # --- Panel 1: detent, vs. rotor position theta_m ---
    n_steps_shown = 2
    theta_deg = np.linspace(0.0, n_steps_shown * pullout_deg, 2000)
    theta_rad = np.radians(theta_deg)
    detent_torque = t_d * np.sin(4.0 * n_r * theta_rad)
    ax_detent.plot(theta_deg, detent_torque, color=DETENT_COLOR, lw=1.8)
    for k in range(n_steps_shown + 1):
        ax_detent.axvline(k * pullout_deg, color=AXIS_COLOR, lw=0.7, ls=':')
    ax_detent.axhline(0.0, color=ZERO_COLOR, lw=0.8)
    ax_detent.set_xlabel('Rotor position theta_m (degrees)')
    ax_detent.set_ylabel('Detent torque T_detent (N·m)')
    ax_detent.set_title(
        f'T_detent(theta_m) = T_d·sin(4·N_r·theta_m)\n'
        f'T_d={t_d*1000:.1f} mN·m, N_r={n_r:.0f} -- '
        f'period = {pullout_deg:.2f}° (one full step)',
        fontsize=10,
    )
    ax_detent.grid(True, alpha=0.3)

    # --- Panel 2: motor drive, vs. tracking error theta_err ---
    err_deg = np.linspace(-n_steps_shown * pullout_deg, n_steps_shown * pullout_deg, 2000)
    err_rad = np.radians(err_deg)
    for label, current_ma in CURRENTS_MA.items():
        t_hold = K_T_NM_PER_A * np.sqrt(2.0) * (current_ma / 1000.0)
        k_em = n_r * t_hold
        ax_motor.plot(err_deg, k_em * err_rad, color=MOTOR_COLORS[label], lw=1.8,
                      label=f'{label}: k_em={k_em:.3f} N·m/rad')
    for sign in (-1, 1):
        ax_motor.axvline(sign * pullout_deg, color=AXIS_COLOR, lw=0.7, ls=':')
    ax_motor.axhline(0.0, color=ZERO_COLOR, lw=0.8)
    ax_motor.set_xlabel('Tracking error theta_err = theta_cmd - theta_m (degrees)')
    ax_motor.set_ylabel('Motor drive torque T_motor (N·m)')
    ax_motor.set_title(
        'T_motor(theta_err) = k_em·theta_err,  k_em = N_r·T_hold\n'
        '(linear drive in lugre_model_rev42; see ../NonlinearStepper+/)',
        fontsize=10,
    )
    ax_motor.grid(True, alpha=0.3)
    ax_motor.legend(loc='upper left', fontsize=8.5)

    fig.suptitle(
        'Figure D -- the two base torque equations, plotted vs. their own '
        'natural variable (not vs. time -- see Figure B/C for time-domain behavior)',
        fontsize=11.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_path = HERE / 'figureD_base_torque_equations.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Step 1: confirm the T_hold anchor (k_em = N_r*T_hold for both drive
currents) and plot the combined rotor tangent stiffness

    k_eff(theta_m) = k_em + 4*N_r*T_d*cos(4*N_r*theta_m)

across one full mechanical step, for both current levels. k_eff going
negative anywhere would mean the detent's own restoring tendency locally
overpowers the electromagnetic spring at that current -- an instability
(cogging jump) the linear small-signal picture doesn't otherwise show.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
sys.path.insert(0, str(REV42_SCRIPTS))
from lugre_model_rev42 import load_parameters

p = load_parameters()
N_r, T_d = p['N_r'], p['T_d']
MOTOR_FULL_STEPS_PER_REV = 200

T_HOLD_BASE_NM = 0.060
I_RATED_BASE_A = 0.400
K_T_NM_PER_A = T_HOLD_BASE_NM / (np.sqrt(2.0) * I_RATED_BASE_A)
CURRENTS_MA = {'I_50pct (200 mA)': 200, 'I_100pct (400 mA)': 400}

theta_m = np.linspace(0.0, 2.0 * np.pi / MOTOR_FULL_STEPS_PER_REV, 2000)  # one full step
detent_tangent = 4.0 * N_r * T_d * np.cos(4.0 * N_r * theta_m)
k_d = 4.0 * N_r * T_d

print(f'N_r={N_r}, T_d={T_d*1000:.3f} mN*m, k_d=4*N_r*T_d={k_d:.4f} N*m/rad')
print(f'K_t={K_T_NM_PER_A:.6f} N*m/A\n')

fig, ax = plt.subplots(figsize=(8.5, 5.5))
theta_steps = theta_m * MOTOR_FULL_STEPS_PER_REV / (2.0 * np.pi)  # fraction of one full step, 0..1
colors = {'I_50pct (200 mA)': '#1f77b4', 'I_100pct (400 mA)': '#d62728'}

for label, peak_ma in CURRENTS_MA.items():
    t_hold = K_T_NM_PER_A * np.sqrt(2.0) * (peak_ma / 1000.0)
    k_em = N_r * t_hold
    k_eff = k_em + detent_tangent
    min_k, max_k = k_eff.min(), k_eff.max()
    stable = 'stable everywhere (k_eff > 0)' if min_k > 0 else f'GOES NEGATIVE (min={min_k:.3f})'
    print(f'{label}: T_hold={t_hold*1000:.1f} mN*m, k_em=N_r*T_hold={k_em:.3f} N*m/rad, '
          f'k_eff range=[{min_k:.3f}, {max_k:.3f}] N*m/rad -- {stable}')

    ax.plot(theta_steps, k_eff, color=colors[label], lw=1.6, label=f'{label}: k_em={k_em:.2f}')
    ax.axhline(k_em, color=colors[label], lw=0.8, linestyle=':', alpha=0.6)

ax.axhline(0.0, color='#333333', lw=1.0)
ax.set_xlabel('Position within one full step (fraction, 0-1)')
ax.set_ylabel('Combined rotor tangent stiffness k_eff (N·m/rad)')
ax.set_title(
    'Rotor tangent stiffness k_eff = k_em + 4·N_r·T_d·cos(4·N_r·θ_m)\n'
    f'N_r={N_r}, T_d={T_d*1000:.1f} mN·m (k_d={k_d:.2f} N·m/rad) -- dotted lines mark k_em anchors'
)
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, alpha=0.3)
fig.tight_layout()

out_path = Path(__file__).resolve().parent / 'step1_rotor_tangent_stiffness.png'
fig.savefig(out_path, dpi=150)
print(f'\nSaved {out_path}')

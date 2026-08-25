#!/usr/bin/env python3
'''Render a clearly labelled preliminary Bode from available checkpoints.'''

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from multisine_rev42 import AMPLITUDE_PERCENT, F_HI_HZ, F_LO_HZ, N_REALIZATIONS, load_case


ROOT = Path(__file__).resolve().parent.parent
CHECKPOINT_DIR = ROOT / 'rendered_assets' / 'npz' / 'multisine_checkpoints'
FIGURE = ROOT / 'rendered_assets' / 'multisine_bode_preliminary.png'
SUMMARY = ROOT / 'rendered_assets' / 'multisine_preliminary_summary.json'
SOLVERS = ('RK45', 'Radau')


def magnitude_db(response):
    return 20.0 * np.log10(np.maximum(np.abs(response), 1.0e-300))


def wrapped_phase(response):
    return (np.angle(response, deg=True) + 180.0) % 360.0 - 180.0


def blend(color, target, fraction):
    return tuple(
        (1.0 - fraction) * np.asarray(color[:3])
        + fraction * np.asarray(target)
    )


def load_groups():
    cases = [load_case(path) for path in sorted(CHECKPOINT_DIR.glob('*.npz'))]
    grouped = {}
    for solver in SOLVERS:
        for amplitude in AMPLITUDE_PERCENT:
            selected = [
                case for case in cases
                if case['solver'] == solver
                and np.isclose(case['amplitude_percent'], amplitude)
            ]
            if not selected:
                continue
            U = np.stack([case['U'] for case in selected])
            Y = np.stack([case['Y'] for case in selected])
            G = np.sum(np.conj(U) * Y, axis=0) / np.sum(
                np.abs(U) ** 2, axis=0
            )
            grouped[(solver, float(amplitude))] = {
                'frequency_hz': selected[0]['frequency_hz'],
                'G': G,
                'count': len(selected),
            }
    return grouped, cases


def main() -> None:
    grouped, cases = load_groups()
    if not grouped:
        raise RuntimeError('No multisine checkpoints are available')
    amplitudes = sorted({key[1] for key in grouped})
    cmap = plt.get_cmap('viridis')
    colors = {
        amplitude: cmap(index / max(len(amplitudes) - 1, 1))
        for index, amplitude in enumerate(amplitudes)
    }
    fig, (ax_mag, ax_phase) = plt.subplots(
        2, 1, figsize=(12.0, 9.0), sharex=True
    )
    for amplitude in amplitudes:
        for solver in SOLVERS:
            key = (solver, amplitude)
            if key not in grouped:
                continue
            data = grouped[key]
            complete = data['count'] == N_REALIZATIONS
            if solver == 'RK45':
                color = blend(colors[amplitude], (1.0, 1.0, 1.0), 0.25)
                linestyle = '--'
            else:
                color = blend(colors[amplitude], (0.0, 0.0, 0.0), 0.10)
                linestyle = '-'
            alpha = 1.0 if complete else 0.45
            linewidth = 1.15 if complete else 0.8
            frequency, response = data['frequency_hz'], data['G']
            count = data['count']
            label = (
                f'{amplitude:g}% {solver} '
                f'(M={count}/{N_REALIZATIONS})'
            )
            ax_mag.plot(
                frequency, magnitude_db(response), color=color,
                linestyle=linestyle, linewidth=linewidth, alpha=alpha,
                label=label,
            )
            ax_phase.plot(
                frequency, wrapped_phase(response), color=color,
                linestyle='none', marker='.', markersize=1.8, alpha=alpha,
            )
    for axis in (ax_mag, ax_phase):
        axis.set_xscale('log')
        axis.set_xlim(F_LO_HZ, F_HI_HZ)
        axis.grid(True, which='both', color='#cccccc', linewidth=0.45)
    ax_phase.set_ylim(-180.0, 180.0)
    ax_phase.set_yticks([-180, -90, 0, 90, 180])
    ax_mag.set_ylabel('Magnitude (dB re 1 m/rad)')
    ax_phase.set_ylabel('Phase (deg, wrapped)')
    ax_phase.set_xlabel('Frequency (Hz)')
    ax_mag.set_title(
        'PRELIMINARY Rev 4.2 multisine BLA from available checkpoints'
    )
    ax_phase.set_title(
        'Partial sets are faded; final curves require M=7 realizations'
    )
    ax_mag.legend(loc='best', fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)

    group_counts = {
        f'{solver}_{amplitude:g}pct': data['count']
        for (solver, amplitude), data in grouped.items()
    }
    summary = {
        'status': 'preliminary',
        'checkpoint_cases_used': len(cases),
        'required_realizations_per_group': N_REALIZATIONS,
        'group_realization_counts': group_counts,
        'figure': FIGURE.name,
        'warning': (
            'Groups with fewer than seven realizations are provisional and '
            'will change as the background sweep completes.'
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

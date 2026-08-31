#!/usr/bin/env python3
"""Local Bode overlays: frictionless, detent, and LuGre; L=1 and 2 mm."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from lugre_model_rev42 import (
    LuGreModelRev42, N_Q, build_structural_matrices, load_parameters,
)

ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = ROOT / 'rendered_assets'
NPZ_DIR = ASSET_DIR / 'npz'
V_STAGE = 5.0e-3
FREQ_HZ = np.geomspace(0.1, 8000.0, 8000)
PITCHES = (1.0e-3, 2.0e-3)
CASES = ('frictionless', 'frictionless_detent', 'lugre_detent')
LABELS = {
    'frictionless': 'Frictionless baseline (detent disabled)',
    'frictionless_detent': 'Frictionless + nonlinear detent tangent',
    'lugre_detent': 'LuGre friction + nonlinear detent tangent',
}
STYLES = {
    'frictionless': ('#2468a2', '-'),
    'frictionless_detent': ('#e08b22', '--'),
    'lugre_detent': ('#b83a3a', '-'),
}
PITCH_STYLES = {1.0e-3: '-', 2.0e-3: (0, (4, 2))}


def frictionless_state_space(p, include_detent):
    mass, damping, stiffness, command = build_structural_matrices(p)
    if include_detent:
        stiffness = stiffness.copy()
        stiffness[0, 0] += 4.0 * p['N_r'] * p['T_d']
    mass_inverse = np.diag(1.0 / np.diag(mass))
    system = np.zeros((2 * N_Q, 2 * N_Q))
    system[:N_Q, N_Q:] = np.eye(N_Q)
    system[N_Q:, :N_Q] = -mass_inverse @ stiffness
    system[N_Q:, N_Q:] = -mass_inverse @ damping
    input_vector = np.zeros(2 * N_Q)
    input_vector[N_Q:] = mass_inverse @ command
    output = np.zeros((1, 2 * N_Q))
    output[0, 5] = 1.0
    return system, input_vector, output


def frequency_response(system, input_vector, output):
    identity = np.eye(system.shape[0])
    values = np.empty(FREQ_HZ.size, dtype=complex)
    for index, frequency_hz in enumerate(FREQ_HZ):
        dynamic = 1j * 2.0 * np.pi * frequency_hz * identity - system
        values[index] = (output @ np.linalg.solve(dynamic, input_vector))[0]
    return values


def evaluate(system, input_vector, output):
    eigenvalues = np.linalg.eigvals(system)
    maximum_real = float(np.max(np.real(eigenvalues)))
    if maximum_real >= 0.0:
        raise RuntimeError(
            f'Unstable tangent: max real eigenvalue {maximum_real:g}'
        )
    response = frequency_response(system, input_vector, output)
    magnitude = 20.0 * np.log10(
        np.maximum(np.abs(response), 1.0e-300)
    )
    phase = np.unwrap(np.angle(response)) * 180.0 / np.pi
    dc = float(np.real((output @ np.linalg.solve(-system, input_vector))[0]))
    return {
        'response': response, 'magnitude_db': magnitude,
        'phase_deg': phase, 'eigenvalues': eigenvalues,
        'dc_gain_m_per_rad': dc,
        'maximum_real_eigenvalue_per_s': maximum_real,
    }

def build_pitch(pitch):
    parameters = load_parameters()
    parameters['L'] = pitch
    systems = {
        'frictionless': frictionless_state_space(parameters, False),
        'frictionless_detent': frictionless_state_space(parameters, True),
    }
    model = LuGreModelRev42(
        parameters, enforce_interface_power=False
    )
    operating_state = model.cruise_state(V_STAGE)
    systems['lugre_detent'] = model.analytical_linearization(
        operating_state
    )
    result = {case: evaluate(*systems[case]) for case in CASES}
    result['operating_state'] = operating_state
    return result


def padded_limits(arrays):
    minimum = min(float(np.min(values)) for values in arrays)
    maximum = max(float(np.max(values)) for values in arrays)
    padding = max(0.04 * (maximum - minimum), 1.0)
    return minimum - padding, maximum + padding


def plot_pitch(pitch, results, magnitude_limits, phase_limits):
    fig, (magnitude_axis, phase_axis) = plt.subplots(
        2, 1, figsize=(11.5, 8.2), sharex=True
    )
    for case in CASES:
        color, linestyle = STYLES[case]
        magnitude_axis.semilogx(
            FREQ_HZ, results[case]['magnitude_db'], color=color,
            linestyle=linestyle, linewidth=1.35, label=LABELS[case],
        )
        phase_axis.semilogx(
            FREQ_HZ, results[case]['phase_deg'], color=color,
            linestyle=linestyle, linewidth=1.35,
        )
    pitch_mm = pitch * 1.0e3
    magnitude_axis.set_ylabel('Magnitude (dB re 1 m/rad)')
    magnitude_axis.set_ylim(*magnitude_limits)
    magnitude_axis.set_title(
        'Command-to-stage Bode overlay: x_n(s) / theta_cmd(s) '
        f'- lead-screw pitch {pitch_mm:.0f} mm'
    )
    magnitude_axis.grid(
        True, which='both', linewidth=0.4, color='#cccccc'
    )
    magnitude_axis.legend(loc='best', framealpha=0.95)
    phase_axis.set_xlabel('Frequency (Hz)')
    phase_axis.set_ylabel('Unwrapped phase (deg)')
    phase_axis.set_xlim(FREQ_HZ[0], FREQ_HZ[-1])
    phase_axis.set_ylim(*phase_limits)
    phase_axis.grid(True, which='both', linewidth=0.4, color='#cccccc')
    phase_axis.text(
        0.01, 0.025,
        'Detent tangent at theta_m=0; LuGre analytical tangent at '
        f'V_stage={V_STAGE * 1e3:.1f} mm/s',
        transform=phase_axis.transAxes, fontsize=8.5, color='#444444',
    )
    fig.tight_layout()
    output = ASSET_DIR / (
        f'bode_three_case_overlay_pitch_{pitch_mm:.0f}mm.png'
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_combined(all_results, magnitude_limits, phase_limits):
    fig, (magnitude_axis, phase_axis) = plt.subplots(
        2, 1, figsize=(12.2, 8.8), sharex=True
    )
    for pitch in PITCHES:
        linestyle = PITCH_STYLES[pitch]
        for case in CASES:
            color = STYLES[case][0]
            magnitude_axis.semilogx(
                FREQ_HZ, all_results[pitch][case]['magnitude_db'],
                color=color, linestyle=linestyle, linewidth=1.4,
            )
            phase_axis.semilogx(
                FREQ_HZ, all_results[pitch][case]['phase_deg'],
                color=color, linestyle=linestyle, linewidth=1.4,
            )

    magnitude_axis.set_ylabel('Magnitude (dB re 1 m/rad)')
    magnitude_axis.set_ylim(*magnitude_limits)
    magnitude_axis.set_title(
        'Command-to-stage Bode: 1 mm versus 2 mm lead-screw pitch'
    )
    magnitude_axis.grid(
        True, which='both', linewidth=0.4, color='#cccccc'
    )
    phase_axis.set_xlabel('Frequency (Hz)')
    phase_axis.set_ylabel('Unwrapped phase (deg)')
    phase_axis.set_xlim(FREQ_HZ[0], FREQ_HZ[-1])
    phase_axis.set_ylim(*phase_limits)
    phase_axis.grid(True, which='both', linewidth=0.4, color='#cccccc')

    model_handles = [
        Line2D([0], [0], color=STYLES[case][0], linewidth=2.0,
               label=LABELS[case])
        for case in CASES
    ]
    pitch_handles = [
        Line2D([0], [0], color='#222222', linewidth=2.0,
               linestyle=PITCH_STYLES[pitch],
               label=f'{pitch * 1e3:.0f} mm pitch')
        for pitch in PITCHES
    ]
    fig.legend(
        handles=model_handles + pitch_handles, loc='lower center',
        ncol=3, framealpha=0.95, bbox_to_anchor=(0.5, 0.005),
    )
    fig.tight_layout(rect=(0.0, 0.09, 1.0, 1.0))
    output = ASSET_DIR / 'bode_three_case_overlay_pitch_1mm_vs_2mm.png'
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def main():
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    all_results = {pitch: build_pitch(pitch) for pitch in PITCHES}
    magnitude_limits = padded_limits([
        all_results[pitch][case]['magnitude_db']
        for pitch in PITCHES for case in CASES
    ])
    phase_limits = padded_limits([
        all_results[pitch][case]['phase_deg']
        for pitch in PITCHES for case in CASES
    ])
    figures = [
        plot_pitch(
            pitch, all_results[pitch], magnitude_limits, phase_limits
        )
        for pitch in PITCHES
    ]
    figures.append(
        plot_combined(all_results, magnitude_limits, phase_limits)
    )

    arrays = {'frequency_hz': FREQ_HZ}
    summary = {
        'model': 'Rev 4.2 locally linearized three-case comparison',
        'stage_operating_velocity_m_per_s': V_STAGE,
        'frequency_range_hz': [float(FREQ_HZ[0]), float(FREQ_HZ[-1])],
        'detent_linearization': 'theta_m=0; k_d=4*N_r*T_d',
        'case_labels': LABELS,
        'figures': [path.name for path in figures],
        'pitches': {},
    }
    for pitch in PITCHES:
        pitch_key = f'pitch_{pitch * 1e3:.0f}mm'
        summary['pitches'][pitch_key] = {}
        arrays[f'{pitch_key}_operating_state'] = (
            all_results[pitch]['operating_state']
        )
        for case in CASES:
            result = all_results[pitch][case]
            for field in (
                'response', 'magnitude_db', 'phase_deg', 'eigenvalues'
            ):
                arrays[f'{pitch_key}_{case}_{field}'] = result[field]
            summary['pitches'][pitch_key][case] = {
                'dc_gain_m_per_rad': result['dc_gain_m_per_rad'],
                'maximum_real_eigenvalue_per_s': (
                    result['maximum_real_eigenvalue_per_s']
                ),
            }

    data_path = NPZ_DIR / 'bode_three_case_pitch_comparison.npz'
    np.savez(data_path, **arrays)
    summary['data'] = str(data_path.relative_to(ROOT))
    summary_path = ASSET_DIR / 'bode_three_case_pitch_comparison.json'
    summary_path.write_text(
        json.dumps(summary, indent=2) + '\n', encoding='utf-8'
    )
    print(json.dumps(summary, indent=2))
    for path in figures:
        print(f'Wrote {path}')
    print(f'Wrote {data_path}')
    print(f'Wrote {summary_path}')


if __name__ == '__main__':
    main()

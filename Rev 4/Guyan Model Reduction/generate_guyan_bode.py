#!/usr/bin/env python3
'''Validate the supplied classical Guyan reduction and render its Bode plot.'''

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from guyan_model import (
    closed_form_mass,
    closed_form_ratios,
    closed_form_stiffness,
    closed_form_transformation,
    fixed_interface_frequencies_hz,
    frequency_response,
    load_parameters,
    modal_frequencies_hz,
    numerical_transformation,
    reduce_model,
)


ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / 'rendered_assets'
NPZ_DIR = ASSET_DIR / 'npz'
FIGURE = ASSET_DIR / 'guyan_bode.png'
OVERLAY_FIGURE = ASSET_DIR / 'guyan_vs_frictionless_bode.png'
SUMMARY = ASSET_DIR / 'guyan_model_summary.json'
DATA = NPZ_DIR / 'guyan_bode_data.npz'
OVERLAY_DATA = NPZ_DIR / 'guyan_vs_frictionless_bode_data.npz'


def resonance_peak_indices(
    frequencies_hz: np.ndarray,
    magnitude_db: np.ndarray,
    display_limit_hz: float,
) -> np.ndarray:
    '''Return the visible resonance maxima, excluding insignificant ripples.'''
    indices = np.flatnonzero(
        (magnitude_db[1:-1] > magnitude_db[:-2])
        & (magnitude_db[1:-1] >= magnitude_db[2:])
    ) + 1
    return indices[
        (frequencies_hz[indices] <= display_limit_hz)
        & (magnitude_db[indices] > -100.0)
    ]


def descending_threshold_crossing(
    frequencies_hz: np.ndarray,
    magnitude_db: np.ndarray,
    final_peak_index: int,
    threshold_db: float,
) -> float:
    '''Interpolate the first threshold crossing after the final resonance.'''
    candidates = np.flatnonzero(
        (np.arange(frequencies_hz.size) > final_peak_index)
        & (magnitude_db <= threshold_db)
    )
    if candidates.size == 0:
        raise AssertionError('No post-resonance threshold crossing found')
    upper = int(candidates[0])
    lower = upper - 1
    fraction = (
        threshold_db - magnitude_db[lower]
    ) / (
        magnitude_db[upper] - magnitude_db[lower]
    )
    return float(
        frequencies_hz[lower]
        + fraction * (frequencies_hz[upper] - frequencies_hz[lower])
    )


def main() -> None:
    p = load_parameters()
    reduced = reduce_model(p)
    ratios = closed_form_ratios(p)
    transformation_numeric = numerical_transformation(reduced['K_full'])
    mass_closed = closed_form_mass(p)
    stiffness_closed = closed_form_stiffness(p)
    frequencies_hz = np.linspace(0.0, 8000.0, 80001)
    response = frequency_response(
        reduced['M'], reduced['C'], reduced['K'], reduced['b'],
        np.array([0.0, 1.0]), frequencies_hz,
    )
    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(response), 1.0e-300))
    phase_deg = np.unwrap(np.angle(response)) * 180.0 / np.pi
    baseline_response = frequency_response(
        reduced['M_full'], reduced['C_full'], reduced['K_full'],
        reduced['b_full'], np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        frequencies_hz,
    )
    baseline_magnitude_db = 20.0 * np.log10(
        np.maximum(np.abs(baseline_response), 1.0e-300)
    )
    baseline_phase_deg = (
        np.unwrap(np.angle(baseline_response)) * 180.0 / np.pi
    )

    requested_display_limit_hz = 1300.0
    significance_threshold_db = -110.0
    peak_indices = resonance_peak_indices(
        frequencies_hz, magnitude_db, requested_display_limit_hz
    )
    baseline_peak_indices = resonance_peak_indices(
        frequencies_hz, baseline_magnitude_db, requested_display_limit_hz
    )
    if peak_indices.size == 0 or baseline_peak_indices.size == 0:
        raise AssertionError('No resonance peaks found in the display band')
    significance_crossing_hz = descending_threshold_crossing(
        frequencies_hz, magnitude_db, int(peak_indices[-1]),
        significance_threshold_db,
    )
    baseline_significance_crossing_hz = descending_threshold_crossing(
        frequencies_hz, baseline_magnitude_db,
        int(baseline_peak_indices[-1]), significance_threshold_db,
    )
    overlay_significance_crossing_hz = max(
        significance_crossing_hz, baseline_significance_crossing_hz
    )
    # The exact crossing is 1328.6 Hz, just outside the requested 1300 Hz.
    # Round up only enough to keep both the crossing and shaded tail visible.
    display_limit_hz = max(
        requested_display_limit_hz,
        float(np.ceil(overlay_significance_crossing_hz / 50.0) * 50.0),
    )
    display_mask = frequencies_hz <= display_limit_hz

    transformation_error = float(np.max(np.abs(
        reduced['T'] - transformation_numeric
    )))
    stiffness_error = float(np.max(np.abs(reduced['K'] - stiffness_closed)))
    mass_error = float(np.max(np.abs(reduced['M'] - mass_closed)))
    port_dependency_error = float(np.max(np.abs(
        reduced['J_nut']
        - (1.0 - ratios['beta'])
        * (ratios['ell'] * reduced['J_sb'] - reduced['J_way'])
    )))
    dc_projected = float(np.linalg.solve(reduced['K'], reduced['b'])[1])
    dc_full = float(np.linalg.solve(
        reduced['K_full'], reduced['b_full']
    )[5])
    if transformation_error > 1.0e-10:
        raise AssertionError('Closed-form and numerical transformations differ')
    if stiffness_error > 1.0e-6:
        raise AssertionError('Closed-form and projected stiffness differ')
    if mass_error > 1.0e-12:
        raise AssertionError('Closed-form and projected mass differ')
    if port_dependency_error > 1.0e-12:
        raise AssertionError('Reduced port dependency does not match')
    if not np.isclose(dc_projected, dc_full, rtol=1.0e-12, atol=1.0e-15):
        raise AssertionError('Guyan projection did not preserve the DC gain')

    summary = {
        'model': 'classical two-master Guyan reduction',
        'masters': ['theta_m', 'x_n'],
        'slaves': ['theta_c', 'theta_s', 'theta_sb', 'x_s'],
        'mass_tuning_applied': False,
        'irs_applied': False,
        'command_vector_convention': 'b=[k_EM,0]^T exactly as supplied',
        'ratios': ratios,
        'transformation': reduced['T'].tolist(),
        'reduced_mass': reduced['M'].tolist(),
        'reduced_damping': reduced['C'].tolist(),
        'reduced_stiffness': reduced['K'].tolist(),
        'reduced_command': reduced['b'].tolist(),
        'reduced_ports': {
            'way': reduced['J_way'].tolist(),
            'nut': reduced['J_nut'].tolist(),
            'sb': reduced['J_sb'].tolist(),
        },
        'reduced_modes_hz': modal_frequencies_hz(
            reduced['M'], reduced['K']
        ).tolist(),
        'full_modes_hz': modal_frequencies_hz(
            reduced['M_full'], reduced['K_full']
        ).tolist(),
        'fixed_interface_modes_hz': fixed_interface_frequencies_hz(
            reduced['M_full'], reduced['K_full']
        ).tolist(),
        'dc_gain_reduced_m_per_rad': dc_projected,
        'dc_gain_full_m_per_rad': dc_full,
        'closed_form_transformation_max_error': transformation_error,
        'closed_form_stiffness_max_error': stiffness_error,
        'closed_form_mass_max_error': mass_error,
        'port_dependency_max_error': port_dependency_error,
        'requested_display_limit_hz': requested_display_limit_hz,
        'rendered_display_limit_hz': display_limit_hz,
        'significance_threshold_db': significance_threshold_db,
        'post_resonance_threshold_crossing_hz': significance_crossing_hz,
        'baseline_post_resonance_threshold_crossing_hz': (
            baseline_significance_crossing_hz
        ),
        'overlay_shading_start_hz': overlay_significance_crossing_hz,
        'labelled_resonance_peak_frequencies_hz': (
            frequencies_hz[peak_indices].tolist()
        ),
        'baseline_labelled_resonance_peak_frequencies_hz': (
            frequencies_hz[baseline_peak_indices].tolist()
        ),
        'figure': FIGURE.name,
        'data': str(DATA.relative_to(ROOT)),
        'overlay_figure': OVERLAY_FIGURE.name,
        'overlay_data': str(OVERLAY_DATA.relative_to(ROOT)),
        'overlay_status': 'generated with documented b=[k_EM,0,...]^T',
    }
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    np.savez_compressed(
        DATA,
        frequencies_hz=frequencies_hz,
        response=response,
        magnitude_db=magnitude_db,
        phase_deg=phase_deg,
        transformation=reduced['T'],
        reduced_mass=reduced['M'],
        reduced_damping=reduced['C'],
        reduced_stiffness=reduced['K'],
        reduced_command=reduced['b'],
    )
    np.savez_compressed(
        OVERLAY_DATA,
        frequencies_hz=frequencies_hz,
        guyan_response=response,
        guyan_magnitude_db=magnitude_db,
        guyan_phase_deg=phase_deg,
        frictionless_response=baseline_response,
        frictionless_magnitude_db=baseline_magnitude_db,
        frictionless_phase_deg=baseline_phase_deg,
        display_limit_hz=display_limit_hz,
        shading_start_hz=overlay_significance_crossing_hz,
    )

    fig, (ax_mag, ax_phase) = plt.subplots(
        2, 1, figsize=(10.5, 8.0), sharex=True
    )
    ax_mag.plot(
        frequencies_hz[display_mask], magnitude_db[display_mask],
        color='#6a3d9a', linewidth=1.25,
    )
    ax_phase.plot(
        frequencies_hz[display_mask], phase_deg[display_mask],
        color='#6a3d9a', linewidth=1.25,
    )
    for axis in (ax_mag, ax_phase):
        axis.axvspan(
            significance_crossing_hz, display_limit_hz,
            color='#ef4444', alpha=0.12, linewidth=0.0, zorder=0,
        )
    ax_mag.scatter(
        frequencies_hz[peak_indices], magnitude_db[peak_indices],
        color='#6a3d9a', edgecolor='white', linewidth=0.7,
        s=30.0, zorder=3,
    )
    for peak_index in peak_indices:
        ax_mag.annotate(
            f'{frequencies_hz[peak_index]:.1f} Hz',
            xy=(frequencies_hz[peak_index], magnitude_db[peak_index]),
            xytext=(0, 9), textcoords='offset points',
            ha='center', va='bottom', fontsize=8.5, color='#4c1d95',
        )
    ax_mag.set_ylabel('Magnitude (dB re 1 m/rad)')
    ax_phase.set_ylabel('Phase (deg)')
    ax_phase.set_xlabel('Frequency (Hz)')
    ax_phase.set_xlim(0.0, display_limit_hz)
    ax_mag.set_title(
        'Classical Guyan reduction: command-to-stage Bode, '
        r'$x_n(s)/\theta_{cmd}(s)$'
    )
    ax_phase.set_title(
        'Masters [theta_m, x_n]; no mass tuning or IRS correction'
    )
    for axis in (ax_mag, ax_phase):
        axis.grid(True, color='#cccccc', linewidth=0.45)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=150)
    plt.close(fig)

    fig, (ax_mag, ax_phase) = plt.subplots(
        2, 1, figsize=(10.5, 8.0), sharex=True
    )
    ax_mag.plot(
        frequencies_hz[display_mask],
        baseline_magnitude_db[display_mask],
        color='#374151', linewidth=1.3,
        label='Full frictionless baseline (6 DOF)',
    )
    ax_mag.plot(
        frequencies_hz[display_mask], magnitude_db[display_mask],
        color='#6a3d9a', linewidth=1.3, linestyle='--',
        label='Classical Guyan reduction (2 DOF)',
    )
    ax_phase.plot(
        frequencies_hz[display_mask], baseline_phase_deg[display_mask],
        color='#374151', linewidth=1.3,
    )
    ax_phase.plot(
        frequencies_hz[display_mask], phase_deg[display_mask],
        color='#6a3d9a', linewidth=1.3, linestyle='--',
    )
    for axis in (ax_mag, ax_phase):
        axis.axvspan(
            overlay_significance_crossing_hz, display_limit_hz,
            color='#ef4444', alpha=0.12, linewidth=0.0, zorder=0,
        )
        axis.grid(True, color='#cccccc', linewidth=0.45)

    ax_mag.scatter(
        frequencies_hz[baseline_peak_indices],
        baseline_magnitude_db[baseline_peak_indices],
        color='#374151', edgecolor='white', linewidth=0.7,
        s=30.0, zorder=3,
    )
    ax_mag.scatter(
        frequencies_hz[peak_indices], magnitude_db[peak_indices],
        color='#6a3d9a', edgecolor='white', linewidth=0.7,
        s=30.0, zorder=3,
    )
    baseline_offsets = [(-8, 10), (-8, -12)]
    guyan_offsets = [(8, -16), (8, 10)]
    for peak_index, offset in zip(
        baseline_peak_indices, baseline_offsets, strict=True
    ):
        ax_mag.annotate(
            f'{frequencies_hz[peak_index]:.1f} Hz',
            xy=(
                frequencies_hz[peak_index],
                baseline_magnitude_db[peak_index],
            ),
            xytext=offset, textcoords='offset points',
            ha='right', va='bottom' if offset[1] > 0 else 'top',
            fontsize=8.5, color='#374151',
        )
    for peak_index, offset in zip(
        peak_indices, guyan_offsets, strict=True
    ):
        ax_mag.annotate(
            f'{frequencies_hz[peak_index]:.1f} Hz',
            xy=(frequencies_hz[peak_index], magnitude_db[peak_index]),
            xytext=offset, textcoords='offset points',
            ha='left', va='bottom' if offset[1] > 0 else 'top',
            fontsize=8.5, color='#4c1d95',
        )

    ax_mag.set_ylabel('Magnitude (dB re 1 m/rad)')
    ax_phase.set_ylabel('Phase (deg)')
    ax_phase.set_xlabel('Frequency (Hz)')
    ax_phase.set_xlim(0.0, display_limit_hz)
    ax_mag.set_title(
        'Frictionless baseline vs classical Guyan reduction, '
        r'$x_n(s)/\theta_{cmd}(s)$'
    )
    ax_phase.set_title(
        r'Same documented command vector $b=[k_{EM},0,\ldots]^T$'
    )
    ax_mag.legend(loc='upper right', frameon=True, framealpha=0.92)
    fig.tight_layout()
    fig.savefig(OVERLAY_FIGURE, dpi=150)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

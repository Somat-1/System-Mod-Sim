#!/usr/bin/env python3
"""Simplified single-edge detent/friction comparison, for presentation use.

Derived from generate_stepping_0p3um_settled.py: reuses the same LuGre
model, frictionless baselines, and full-periodic-detent frictionless
overlay, but only renders the "first edge from rest" panels (position and
tracking error), merges the Newton and Lagrange frictionless baselines
into a single curve (they are numerically identical -- see
maximum_newton_vs_lagrange_stage_difference_um in the 0.3 um summary),
relabels the LuGre curve, and drops the figure title/configuration-note
annotation. Runs at a configurable per-edge step size so the same
comparison can be regenerated at, e.g., 1-10 um instead of the original
0.3 um probe.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from generate_stepping_0p3um_settled import (
    FIRING_INTERVAL_S,
    MOVES,
    OUTPUT_DT_S,
    SUBSTEPS_PER_SEQUENCE_UNIT,
    LuGreModelRev42,
    baseline_systems,
    metrics,
    simulate_full_detent_frictionless,
    simulate_linear,
    simulate_lugre,
)

ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = ROOT / 'rendered_assets'
NPZ_DIR = ASSET_DIR / 'npz'


def command_travel_edges(step_travel_m: float) -> np.ndarray:
    travel = 0.0
    edges: list[float] = []
    for move in MOVES:
        for _ in range(abs(move) * SUBSTEPS_PER_SEQUENCE_UNIT):
            travel += np.sign(move) * step_travel_m
            edges.append(travel)
    result = np.asarray(edges)
    if not np.isclose(result[-1], 0.0, atol=1.0e-12):
        raise AssertionError('The stepping sequence must return to zero')
    return result


LEGACY_0P3UM_DATA = NPZ_DIR / 'stepping_0p3um_settled.npz'


def run_case(step_travel_um: float):
    step_travel_m = step_travel_um * 1.0e-6
    tag = f'{step_travel_um:g}um'
    data_path = NPZ_DIR / f'stepping_{tag}_settled.npz'
    summary_path = ASSET_DIR / f'stepping_{tag}_settled_summary.json'
    figure_path = ASSET_DIR / f'stepping_{tag}_settled_simplified.png'

    model = LuGreModelRev42(enforce_interface_power=False)
    lead = model.p['L'] / (2.0 * np.pi)
    travel_commands = command_travel_edges(step_travel_m)
    theta_commands = travel_commands / lead
    newton_system, lagrange_system, _conventions = baseline_systems()

    # The 0.3 um case was already computed by generate_stepping_0p3um_settled.py
    # under a different filename/key scheme; reuse it directly instead of
    # re-running ~256 Radau edges.
    legacy_key_map = {
        'lugre': 'lugre', 'newton': 'newton', 'lagrange': 'lagrange',
        'full_detent': 'full_detent_frictionless',
    }
    cached = None
    if step_travel_um == 0.3 and LEGACY_0P3UM_DATA.exists():
        try:
            with np.load(LEGACY_0P3UM_DATA) as payload:
                if np.array_equal(payload['command_edges_m'], travel_commands):
                    cached = {
                        'progress': payload['progress'].copy(),
                        'time_s': payload['time_s'].copy(),
                        'command_travel_m': payload['command_travel_m'].copy(),
                    }
                    for new_prefix, legacy_prefix in legacy_key_map.items():
                        cached[f'{new_prefix}_x_n_m'] = payload[
                            f'{legacy_prefix}_x_n_m'
                        ].copy()
                        cached[f'{new_prefix}_error_m'] = payload[
                            f'{legacy_prefix}_error_m'
                        ].copy()
        except (OSError, ValueError, EOFError, KeyError):
            cached = None
    if cached is None and data_path.exists():
        try:
            with np.load(data_path) as payload:
                if np.array_equal(payload['command_edges_m'], travel_commands):
                    cached = {
                        key: payload[key].copy() for key in payload.files
                    }
        except (OSError, ValueError, EOFError):
            cached = None

    def cached_or_run(key_prefix, runner):
        if cached is not None and f'{key_prefix}_x_n_m' in cached:
            print(f'Reusing cached {key_prefix} trajectory ({tag}).', flush=True)
            return {
                'progress': cached['progress'],
                'time_s': cached['time_s'],
                'command_travel_m': cached['command_travel_m'],
                'x_n_m': cached[f'{key_prefix}_x_n_m'],
                'error_m': cached[f'{key_prefix}_error_m'],
            }
        print(f'Running {key_prefix} ({tag})...', flush=True)
        return runner()

    lugre = cached_or_run(
        'lugre', lambda: simulate_lugre(model, theta_commands, lead)
    )
    newton = cached_or_run(
        'newton', lambda: simulate_linear(newton_system, theta_commands, lead)
    )
    lagrange = cached_or_run(
        'lagrange',
        lambda: simulate_linear(lagrange_system, theta_commands, lead),
    )
    full_detent = cached_or_run(
        'full_detent',
        lambda: simulate_full_detent_frictionless(model, theta_commands, lead),
    )

    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        data_path,
        progress=lugre['progress'],
        time_s=lugre['time_s'],
        command_travel_m=lugre['command_travel_m'],
        lugre_x_n_m=lugre['x_n_m'], lugre_error_m=lugre['error_m'],
        newton_x_n_m=newton['x_n_m'], newton_error_m=newton['error_m'],
        lagrange_x_n_m=lagrange['x_n_m'], lagrange_error_m=lagrange['error_m'],
        full_detent_x_n_m=full_detent['x_n_m'],
        full_detent_error_m=full_detent['error_m'],
        command_edges_m=travel_commands,
    )

    max_newton_lagrange_diff_um = float(
        np.max(np.abs(newton['x_n_m'] - lagrange['x_n_m'])) * 1.0e6
    )
    summary = {
        'step_travel_um': step_travel_um,
        'command_edges': int(theta_commands.size),
        'results': {
            'lugre_friction_and_detent': metrics(lugre),
            'frictionless_baseline': metrics(newton),
            'frictionless_full_periodic_detent': metrics(full_detent),
        },
        'maximum_newton_vs_lagrange_stage_difference_um': (
            max_newton_lagrange_diff_um
        ),
        'note': (
            'Newton and Lagrange frictionless baselines differ by at most '
            f'{max_newton_lagrange_diff_um:.3g} um and are merged into a '
            'single curve in the simplified figure.'
        ),
    }
    summary_path.write_text(
        json.dumps(summary, indent=2) + '\n', encoding='utf-8'
    )
    print(json.dumps(summary, indent=2))

    render_simplified(
        lugre, newton, full_detent, step_travel_um, figure_path,
    )
    return summary


def render_simplified(
    lugre, newton_or_lagrange, full_detent, step_travel_um: float,
    figure_path: Path,
) -> None:
    fig, (ax_position, ax_error) = plt.subplots(
        2, 1, figsize=(9.5, 8.0), sharex=True,
    )

    progress = lugre['progress']
    styles = [
        (lugre, '#c0392b', '-', 1.15, 'LuGre Friction & Nonlinear Detent'),
        (
            newton_or_lagrange, '#2b6cb0', '-', 1.25,
            'Frictionless baseline (Newton = Lagrange)',
        ),
        (
            full_detent, '#d97706', '-.', 1.25,
            'Frictionless, full periodic detent',
        ),
    ]

    ax_position.plot(
        progress, lugre['command_travel_m'] * 1.0e6,
        color='#777777', linestyle='--', linewidth=0.9, label='commanded',
    )
    for result, color, linestyle, linewidth, label in styles:
        ax_position.plot(
            result['progress'], result['x_n_m'] * 1.0e6,
            color=color, linestyle=linestyle, linewidth=linewidth,
            label=label,
        )
        ax_error.plot(
            result['progress'], result['error_m'] * 1.0e6,
            color=color, linestyle=linestyle, linewidth=linewidth,
        )
    ax_error.axhline(0.0, color='#333333', linestyle=':', linewidth=0.7)

    for axis in (ax_position, ax_error):
        axis.grid(True, color='#cccccc', linewidth=0.45)
        axis.set_xlim(0.0, 16.0)
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

    ax_position.set_ylabel(r'Stage position ($\mu$m)')
    ax_position.set_title(
        f'Complete 256-edge sequence ({step_travel_um:g} $\\mu$m increments)'
    )
    ax_position.legend(loc='best', framealpha=0.92, fontsize=9)
    ax_error.set_ylabel(r'Tracking error ($\mu$m)')
    ax_error.set_xlabel('Sequence progress (nominal full-step units)')
    ax_error.set_title('Complete-sequence tracking error')

    fig.tight_layout()
    temporary_figure = figure_path.with_name(
        f'{figure_path.stem}.tmp{figure_path.suffix}'
    )
    fig.savefig(temporary_figure, dpi=160)
    plt.close(fig)
    temporary_figure.replace(figure_path)


def main() -> None:
    run_case(0.3)
    run_case(10.0)


if __name__ == '__main__':
    main()

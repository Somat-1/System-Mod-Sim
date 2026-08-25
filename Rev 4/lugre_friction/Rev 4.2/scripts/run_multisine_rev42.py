#!/usr/bin/env python3
'''Run, resume, aggregate, and plot the Rev 4.2 multisine experiment.'''

from __future__ import annotations

import os

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lugre_model_rev42 import N_STATES, LuGreModelRev42
from multisine_rev42 import (
    AMPLITUDE_PERCENT,
    BASE_PERIOD_S,
    CONVERGENCE_TOL,
    DESIGN_SEED,
    F_HI_HZ,
    F_LO_HZ,
    FS_HZ,
    MAX_CONVERGENCE_PERIODS,
    MAX_STEP_S,
    N_REALIZATIONS,
    N_RETAINED_PERIODS,
    N_SAMPLES,
    RTOL,
    SOLVERS,
    breakaway_command_rms,
    design_random_odd_bins,
    load_case,
    run_case,
    save_case,
)


ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = ROOT / 'rendered_assets'
NPZ_DIR = ASSET_DIR / 'npz'
CHECKPOINT_DIR = NPZ_DIR / 'multisine_checkpoints'
RESULT_NPZ = NPZ_DIR / 'multisine_rev42_results.npz'
SUMMARY_JSON = ASSET_DIR / 'multisine_rev42_summary.json'
STATUS_JSON = ASSET_DIR / 'multisine_rev42_status.json'
SOLVER_FIGURE = ASSET_DIR / 'multisine_bode_solvers.png'
AMPLITUDE_FIGURE = ASSET_DIR / 'multisine_bode_amplitude_sweep.png'
STATISTICS_FIGURE = ASSET_DIR / 'multisine_run_statistics.png'


def case_path(solver: str, amplitude: float, realization: int) -> Path:
    amplitude_token = f'{amplitude:06.2f}'.replace('.', 'p')
    return CHECKPOINT_DIR / (
        f'multisine_{solver.lower()}_amp_{amplitude_token}_r{realization:02d}.npz'
    )


def work_items(
    solvers: list[str],
    amplitudes: list[float],
    realizations: int,
    max_convergence_periods: int,
) -> list[dict]:
    return [
        {
            'solver': solver,
            'amplitude_percent': amplitude,
            'realization': realization,
            'max_convergence_periods': max_convergence_periods,
        }
        for amplitude in amplitudes
        for realization in range(realizations)
        for solver in solvers
    ]


def status_payload(
    requested: int,
    completed: int,
    cached: int,
    failed: list[dict],
    wall_s: float,
) -> dict:
    return {
        'requested_cases': requested,
        'completed_cases': completed,
        'cached_cases': cached,
        'failed_cases': failed,
        'orchestrator_wall_s': wall_s,
        'last_update_local': time.strftime('%Y-%m-%d %H:%M:%S'),
    }


def write_status(payload: dict) -> None:
    STATUS_JSON.parent.mkdir(parents=True, exist_ok=True)
    STATUS_JSON.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def execute(items: list[dict], workers: int, rerun: bool):
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    results, pending = [], []
    for item in items:
        path = case_path(
            item['solver'], item['amplitude_percent'], item['realization']
        )
        if path.exists() and not rerun:
            results.append(load_case(path))
        else:
            pending.append(item)
    cached, failures = len(results), []
    start = time.perf_counter()
    print(
        f'{len(items)} cases; {cached} cached; {len(pending)} pending; '
        f'{workers} workers', flush=True,
    )
    write_status(status_payload(len(items), cached, cached, failures, 0.0))
    if pending:
        results, failures = execute_pending(
            pending, results, failures, items, cached, workers, start
        )
    wall_s = time.perf_counter() - start
    if failures:
        raise RuntimeError(f'{len(failures)} cases failed; see {STATUS_JSON}')
    return results, wall_s, cached


def execute_pending(
    pending, results, failures, items, cached, workers, start
):
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_case, item): item for item in pending}
        for future in as_completed(futures):
            item = futures[future]
            try:
                result = future.result()
                save_case(case_path(
                    result['solver'], result['amplitude_percent'],
                    result['realization'],
                ), result)
                results.append(result)
                report_completion(result, len(results), len(items))
            except Exception as exc:
                failure = dict(item)
                failure['error'] = repr(exc)
                failures.append(failure)
                print(f'FAILED {item}: {exc!r}', flush=True)
            write_status(status_payload(
                len(items), len(results), cached, failures,
                time.perf_counter() - start,
            ))
    return results, failures


def report_completion(result: dict, completed: int, requested: int) -> None:
    if completed % 35 != 0 and completed != requested:
        return
    stats = result['stats']
    solver = result['solver']
    amplitude = result['amplitude_percent']
    realization = result['realization']
    wall = stats['wall_s']
    periods = stats['total_periods']
    metric = stats['final_periodicity_metric']
    print(
        f'{completed:3d}/{requested} {solver:5s} A={amplitude:6.2f}% Fs '
        f'r={realization} {wall:8.1f}s periods={periods} metric={metric}',
        flush=True,
    )


def group_results(results, solvers, amplitudes, realizations):
    grouped = {}
    for solver in solvers:
        for amplitude in amplitudes:
            cases = sorted(
                (r for r in results if r['solver'] == solver
                 and np.isclose(r['amplitude_percent'], amplitude)),
                key=lambda r: r['realization'],
            )
            if len(cases) != realizations:
                raise RuntimeError(
                    f'Expected {realizations} {solver}/{amplitude}% cases; '
                    f'found {len(cases)}'
                )
            U = np.stack([case['U'] for case in cases])
            Y = np.stack([case['Y'] for case in cases])
            G_each = Y / U
            G = np.sum(np.conj(U) * Y, axis=0) / np.sum(np.abs(U) ** 2, axis=0)
            magnitude_each = 20.0 * np.log10(
                np.maximum(np.abs(G_each), 1.0e-300)
            )
            grouped[(solver, amplitude)] = {
                'cases': cases,
                'frequency_hz': cases[0]['frequency_hz'],
                'G': G,
                'G_each': G_each,
                'magnitude_std_db': np.std(
                    magnitude_each, axis=0,
                    ddof=1 if realizations > 1 else 0,
                ),
            }
    return grouped


def local_rest_response(frequency_hz: np.ndarray) -> np.ndarray:
    model = LuGreModelRev42(enforce_interface_power=False)
    system, input_vector, output = model.analytical_linearization(
        np.zeros(N_STATES)
    )
    identity = np.eye(N_STATES)
    response = np.empty(len(frequency_hz), dtype=np.complex128)
    for index, frequency in enumerate(frequency_hz):
        response[index] = (
            output @ np.linalg.solve(
                1j * 2.0 * np.pi * frequency * identity - system,
                input_vector,
            )
        )[0]
    return response


def wrapped_phase(response: np.ndarray) -> np.ndarray:
    return (np.angle(response, deg=True) + 180.0) % 360.0 - 180.0


def magnitude_db(response: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.maximum(np.abs(response), 1.0e-300))


def blend(color, target, fraction: float):
    color = np.asarray(color[:3])
    target = np.asarray(target)
    return tuple((1.0 - fraction) * color + fraction * target)


def style_axes(axis) -> None:
    axis.set_xscale('log')
    axis.set_xlim(F_LO_HZ, F_HI_HZ)
    axis.grid(True, which='both', color='#cccccc', linewidth=0.45)


def plot_solver_comparison(grouped, amplitudes, solvers) -> None:
    amplitude = min(amplitudes)
    frequency = grouped[(solvers[0], amplitude)]['frequency_hz']
    local = local_rest_response(frequency)
    colors = {'RK45': '#377eb8', 'Radau': '#e6550d'}
    fig, (ax_mag, ax_phase) = plt.subplots(
        2, 1, figsize=(11.0, 8.5), sharex=True
    )
    ax_mag.plot(
        frequency, magnitude_db(local), color='#222222', linewidth=1.2,
        label='local rest tangent',
    )
    ax_phase.plot(
        frequency, wrapped_phase(local), color='#222222', linewidth=1.2,
    )
    for solver in solvers:
        data = grouped[(solver, amplitude)]
        color = colors.get(solver)
        mag, std = magnitude_db(data['G']), data['magnitude_std_db']
        ax_mag.plot(
            frequency, mag, color=color, linewidth=1.15,
            label=f'{solver}, {amplitude:g}% Fs RMS',
        )
        ax_mag.fill_between(
            frequency, mag - std, mag + std, color=color, alpha=0.12,
            linewidth=0.0,
        )
        ax_phase.plot(
            frequency, wrapped_phase(data['G']), color=color, linewidth=0.0,
            marker='.', markersize=2.2,
        )
    style_axes(ax_mag)
    style_axes(ax_phase)
    ax_phase.set_ylim(-180.0, 180.0)
    ax_phase.set_yticks([-180, -90, 0, 90, 180])
    ax_mag.set_ylabel('Magnitude (dB re 1 m/rad)')
    ax_phase.set_ylabel('Phase (deg, wrapped)')
    ax_phase.set_xlabel('Frequency (Hz)')
    ax_mag.set_title(
        'Rev 4.2 random-odd multisine: solver agreement at smallest amplitude'
    )
    ax_phase.set_title('Phase at the 370 excited lines')
    ax_mag.legend(loc='best')
    fig.tight_layout()
    fig.savefig(SOLVER_FIGURE, dpi=150)
    plt.close(fig)


def plot_amplitude_sweep(grouped, amplitudes, solvers) -> None:
    cmap = plt.get_cmap('viridis')
    base_colors = {
        amplitude: cmap(index / max(len(amplitudes) - 1, 1))
        for index, amplitude in enumerate(amplitudes)
    }
    fig, (ax_mag, ax_phase) = plt.subplots(
        2, 1, figsize=(12.5, 9.5), sharex=True
    )
    for amplitude in amplitudes:
        base = base_colors[amplitude]
        for solver in solvers:
            response = grouped[(solver, amplitude)]['G']
            if solver == 'RK45':
                color = blend(base, (1.0, 1.0, 1.0), 0.28)
                linestyle, linewidth = '--', 0.9
            else:
                color = blend(base, (0.0, 0.0, 0.0), 0.10)
                linestyle, linewidth = '-', 1.05
            frequency = grouped[(solver, amplitude)]['frequency_hz']
            ax_mag.plot(
                frequency, magnitude_db(response), color=color,
                linestyle=linestyle, linewidth=linewidth,
            )
            ax_phase.plot(
                frequency, wrapped_phase(response), color=color,
                linestyle='none', marker='.', markersize=1.45,
            )
    style_axes(ax_mag)
    style_axes(ax_phase)
    ax_phase.set_ylim(-180.0, 180.0)
    ax_phase.set_yticks([-180, -90, 0, 90, 180])
    ax_mag.set_ylabel('Magnitude (dB re 1 m/rad)')
    ax_phase.set_ylabel('Phase (deg, wrapped)')
    ax_phase.set_xlabel('Frequency (Hz)')
    ax_mag.set_title(
        'Rev 4.2 multisine amplitude-dependent BLA: pre-sliding to gross slip'
    )
    ax_phase.set_title(
        'Phase shift with amplitude (dots are the 370 excited odd lines)'
    )
    amplitude_handles = [
        Patch(color=base_colors[a], label=f'{a:g}% Fs') for a in amplitudes
    ]
    solver_handles = [
        Line2D(
            [0], [0], color='#555555', linestyle='--',
            label='RK45 (lighter)',
        ),
        Line2D(
            [0], [0], color='#222222', linestyle='-',
            label='Radau (darker)',
        ),
    ]
    first_legend = ax_mag.legend(
        handles=amplitude_handles, loc='upper right', ncol=2, fontsize=7,
        title='RMS amplitude',
    )
    ax_mag.add_artist(first_legend)
    ax_mag.legend(handles=solver_handles, loc='lower left', fontsize=8)
    fig.tight_layout()
    fig.savefig(AMPLITUDE_FIGURE, dpi=150)
    plt.close(fig)


def statistic_row(amplitude: float, cases: list[dict]) -> dict:
    stats = [case['stats'] for case in cases]
    values = lambda key: np.array([item[key] for item in stats], dtype=float)
    return {
        'amplitude_percent_Fs': amplitude,
        'median_wall_s': float(np.median(values('wall_s'))),
        'min_wall_s': float(np.min(values('wall_s'))),
        'max_wall_s': float(np.max(values('wall_s'))),
        'median_convergence_period': float(
            np.median(values('convergence_period'))
        ),
        'maximum_total_periods': int(np.max(values('total_periods'))),
        'converged_realizations': int(sum(item['converged'] for item in stats)),
        'median_nfev': float(np.median(values('nfev'))),
        'median_njev': float(np.median(values('njev'))),
        'median_nlu': float(np.median(values('nlu'))),
        'median_way_saturation': float(np.median([
            item['saturation_max']['way'] for item in stats
        ])),
        'median_way_steady_slip_fraction': float(np.median([
            item['steady_slip_fraction']['way'] for item in stats
        ])),
        'median_crest_factor': float(np.median(values('crest_factor'))),
    }


def aggregate_stats(grouped, amplitudes, solvers):
    summary, arrays = {}, {}
    fields = (
        'median_wall_s', 'median_convergence_period', 'median_nfev',
        'median_way_saturation', 'median_way_steady_slip_fraction',
    )
    for solver in solvers:
        rows = [
            statistic_row(a, grouped[(solver, a)]['cases']) for a in amplitudes
        ]
        summary[solver] = rows
        for field in fields:
            arrays[f'{solver}_{field}'] = np.array(
                [row[field] for row in rows]
            )
    return summary, arrays


def plot_statistics(arrays, amplitudes, solvers) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), sharex=True)
    colors = {'RK45': '#377eb8', 'Radau': '#e6550d'}
    keys = (
        'median_wall_s', 'median_convergence_period',
        'median_nfev', 'median_way_steady_slip_fraction',
    )
    for solver in solvers:
        for axis, key in zip(axes.flat, keys):
            axis.plot(
                amplitudes, arrays[f'{solver}_{key}'], marker='o',
                color=colors.get(solver), label=solver,
            )
    axes[0, 0].set_ylabel('Median wall time / case (s)')
    axes[0, 1].set_ylabel('Median convergence period')
    axes[1, 0].set_ylabel('Median RHS evaluations')
    axes[1, 1].set_ylabel('Guideway steady-slip fraction')
    for axis in axes[1]:
        axis.set_xlabel('Multisine RMS amplitude (% of guideway Fs equivalent)')
    for axis in axes.flat:
        axis.grid(True, color='#cccccc', linewidth=0.45)
    axes[0, 0].legend()
    fig.suptitle('Rev 4.2 multisine run and convergence statistics')
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    fig.savefig(STATISTICS_FIGURE, dpi=150)
    plt.close(fig)


def save_aggregate(grouped, amplitudes, solvers, arrays) -> None:
    first = grouped[(solvers[0], amplitudes[0])]
    payload = {
        'amplitude_percent_Fs': np.asarray(amplitudes),
        'frequency_hz': first['frequency_hz'],
        'excited_bins': first['cases'][0]['excited_bins'],
        'omitted_bins': first['cases'][0]['omitted_bins'],
        'omission_counts': first['cases'][0]['omission_counts'],
    }
    for solver in solvers:
        payload[f'{solver}_G_bla'] = np.stack([
            grouped[(solver, amplitude)]['G'] for amplitude in amplitudes
        ])
        payload[f'{solver}_G_realizations'] = np.stack([
            grouped[(solver, amplitude)]['G_each'] for amplitude in amplitudes
        ])
    payload.update(arrays)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(RESULT_NPZ, **payload)


def solver_disagreement_rows(grouped, amplitudes, solvers) -> list[dict]:
    if not set(('RK45', 'Radau')).issubset(solvers):
        return []
    rows = []
    for amplitude in amplitudes:
        rk = grouped[('RK45', amplitude)]['G']
        radau = grouped[('Radau', amplitude)]['G']
        mag_delta = np.abs(magnitude_db(rk) - magnitude_db(radau))
        phase_delta = np.abs(wrapped_phase(rk / radau))
        rows.append({
            'amplitude_percent_Fs': amplitude,
            'maximum_magnitude_delta_db': float(np.max(mag_delta)),
            'median_magnitude_delta_db': float(np.median(mag_delta)),
            'maximum_phase_delta_deg': float(np.max(phase_delta)),
            'median_phase_delta_deg': float(np.median(phase_delta)),
        })
    return rows


def build_summary(
    grouped, stats_summary, amplitudes, solvers, realizations,
    orchestrator_wall_s, cached, workers,
) -> dict:
    model = LuGreModelRev42(enforce_interface_power=False)
    excited, omitted, omission_counts = design_random_odd_bins()
    serial_wall_s = float(sum(
        case['stats']['wall_s']
        for data in grouped.values() for case in data['cases']
    ))
    return {
        'experiment': 'Rev 4.2 random-odd periodic multisine BLA',
        'base_period_s': BASE_PERIOD_S,
        'sampling_frequency_hz': FS_HZ,
        'samples_per_period': N_SAMPLES,
        'frequency_resolution_hz': 1.0 / BASE_PERIOD_S,
        'requested_band_hz': [F_LO_HZ, F_HI_HZ],
        'actual_excited_band_hz': [
            float(excited[0] / BASE_PERIOD_S),
            float(excited[-1] / BASE_PERIOD_S),
        ],
        'excited_lines': int(len(excited)),
        'odd_detection_lines_omitted': int(len(omitted)),
        'omission_count_range': [
            int(np.min(omission_counts)), int(np.max(omission_counts))
        ],
        'design_seed': DESIGN_SEED,
        'realizations': realizations,
        'retained_periods_after_convergence': N_RETAINED_PERIODS,
        'convergence_tolerance': CONVERGENCE_TOL,
        'maximum_convergence_periods': MAX_CONVERGENCE_PERIODS,
        'rtol': RTOL,
        'maximum_internal_step_s': MAX_STEP_S,
        'solvers': {
            'RK45': 'Dormand-Prince explicit Runge-Kutta order 5(4)',
            'Radau': 'implicit Radau IIA order 5 with analytical Jacobian',
        },
        'amplitude_reference': (
            'RMS command angle equivalent to guideway static bristle '
            'deflection Fs_way/sigma0_way'
        ),
        'breakaway_command_rms_rad': breakaway_command_rms(model),
        'amplitude_percent_Fs': amplitudes,
        'workers': workers,
        'orchestrator_wall_s': orchestrator_wall_s,
        'cached_cases': cached,
        'serial_worker_wall_s': serial_wall_s,
        'effective_parallel_speedup': (
            serial_wall_s / orchestrator_wall_s
            if orchestrator_wall_s > 0.0 else None
        ),
        'per_solver_amplitude_statistics': stats_summary,
        'solver_disagreement': solver_disagreement_rows(
            grouped, amplitudes, solvers
        ),
        'figures': [
            SOLVER_FIGURE.name,
            AMPLITUDE_FIGURE.name,
            STATISTICS_FIGURE.name,
        ],
        'data': str(RESULT_NPZ.relative_to(ROOT)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--workers', type=int,
        default=max(1, min(12, (os.cpu_count() or 4) - 2)),
    )
    parser.add_argument('--realizations', type=int, default=N_REALIZATIONS)
    parser.add_argument(
        '--amplitudes', type=float, nargs='+',
        default=AMPLITUDE_PERCENT.tolist(),
    )
    parser.add_argument(
        '--solvers', nargs='+', choices=SOLVERS, default=list(SOLVERS),
    )
    parser.add_argument(
        '--max-convergence-periods', type=int,
        default=MAX_CONVERGENCE_PERIODS,
    )
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--plot-only', action='store_true')
    return parser.parse_args()


def collect_results(args, items):
    if not args.plot_only:
        return execute(items, args.workers, args.rerun)
    results = [
        load_case(case_path(
            item['solver'], item['amplitude_percent'], item['realization']
        ))
        for item in items
    ]
    return results, 0.0, len(results)


def main() -> None:
    args = parse_args()
    items = work_items(
        args.solvers, args.amplitudes, args.realizations,
        args.max_convergence_periods,
    )
    results, orchestrator_wall_s, cached = collect_results(args, items)
    grouped = group_results(
        results, args.solvers, args.amplitudes, args.realizations
    )
    stats_summary, stats_arrays = aggregate_stats(
        grouped, args.amplitudes, args.solvers
    )
    plot_solver_comparison(grouped, args.amplitudes, args.solvers)
    plot_amplitude_sweep(grouped, args.amplitudes, args.solvers)
    plot_statistics(stats_arrays, args.amplitudes, args.solvers)
    save_aggregate(grouped, args.amplitudes, args.solvers, stats_arrays)
    summary = build_summary(
        grouped, stats_summary, args.amplitudes, args.solvers,
        args.realizations, orchestrator_wall_s, cached, args.workers,
    )
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({
        'summary': str(SUMMARY_JSON),
        'solver_figure': str(SOLVER_FIGURE),
        'amplitude_figure': str(AMPLITUDE_FIGURE),
        'statistics_figure': str(STATISTICS_FIGURE),
        'data': str(RESULT_NPZ),
    }, indent=2))


if __name__ == '__main__':
    main()

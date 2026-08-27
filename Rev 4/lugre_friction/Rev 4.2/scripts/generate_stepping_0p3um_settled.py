#!/usr/bin/env python3
'''Generate the settled 0.300 um Rev 4.2 stepping comparison.

The nonlinear trajectory uses the smoothed parallel LuGre Rev 4.2 model.
Newton/free-body and Lagrange frictionless baselines are independently
assembled with the documented fixed-frame detent convention: k_d enters
the motor stiffness tangent and is absent from the command column.
'''

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path
from zipfile import BadZipFile

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import StateSpace, lsim

from lugre_model_rev42 import N_Q, N_STATES, LuGreModelRev42


ROOT = Path(__file__).resolve().parent.parent
REV4 = ROOT.parents[1]
ASSET_DIR = ROOT / 'rendered_assets'
NPZ_DIR = ASSET_DIR / 'npz'
FIGURE = ASSET_DIR / 'stepping_montage_0p3um_settled.png'
KD0_FIGURE = ASSET_DIR / 'stepping_montage_0p3um_settled_kd0.png'
SUMMARY = ASSET_DIR / 'stepping_0p3um_settled_summary.json'
DATA = NPZ_DIR / 'stepping_0p3um_settled.npz'

MOVES = (2, -1, 1, -1, 1, -4, 1, -1, 1, -1, 2)
assert sum(MOVES) == 0
SUBSTEPS_PER_SEQUENCE_UNIT = 16
STEP_TRAVEL_M = 0.300e-6
FIRING_INTERVAL_S = 250.0e-3
OUTPUT_DT_S = 100.0e-6
RTOL = 1.0e-6
ATOL = np.array([
    1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-12, 1.0e-12,
    1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-9, 1.0e-9,
    1.0e-11, 1.0e-11, 1.0e-9,
])


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Cannot import {path}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def command_travel_edges() -> np.ndarray:
    travel = 0.0
    edges: list[float] = []
    for move in MOVES:
        for _ in range(abs(move) * SUBSTEPS_PER_SEQUENCE_UNIT):
            travel += np.sign(move) * STEP_TRAVEL_M
            edges.append(travel)
    result = np.asarray(edges)
    if not np.isclose(result[-1], 0.0, atol=1.0e-18):
        raise AssertionError('The 0.300 um stepping sequence must return to zero')
    return result


def segment_grid() -> np.ndarray:
    count = int(round(FIRING_INTERVAL_S / OUTPUT_DT_S))
    if not np.isclose(count * OUTPUT_DT_S, FIRING_INTERVAL_S):
        raise AssertionError('The firing interval must divide the output grid')
    return np.linspace(0.0, FIRING_INTERVAL_S, count + 1)


def assemble_result(
    progress_parts: list[np.ndarray],
    time_parts: list[np.ndarray],
    command_parts: list[np.ndarray],
    position_parts: list[np.ndarray],
) -> dict[str, np.ndarray]:
    command_travel = np.concatenate(command_parts)
    position = np.concatenate(position_parts)
    return {
        'progress': np.concatenate(progress_parts),
        'time_s': np.concatenate(time_parts),
        'command_travel_m': command_travel,
        'x_n_m': position,
        'error_m': command_travel - position,
    }


def simulate_lugre(
    model: LuGreModelRev42,
    theta_commands: np.ndarray,
    lead: float,
) -> dict[str, np.ndarray | float | int]:
    local_time = segment_grid()
    kept_time = local_time[:-1]
    state = np.zeros(N_STATES)
    progress_parts: list[np.ndarray] = []
    time_parts: list[np.ndarray] = []
    command_parts: list[np.ndarray] = []
    position_parts: list[np.ndarray] = []
    nfev = njev = nlu = 0
    started = time.perf_counter()

    for segment, theta_command in enumerate(theta_commands):
        solution = solve_ivp(
            lambda t, y, command=theta_command: model.rhs(t, y, command),
            (0.0, FIRING_INTERVAL_S),
            state,
            method='Radau',
            jac=lambda _t, y: model.analytical_linearization(y)[0],
            t_eval=local_time,
            rtol=RTOL,
            atol=ATOL,
        )
        if not solution.success:
            raise RuntimeError(
                f'LuGre integration failed at edge {segment}: '
                f'{solution.message}'
            )
        state = solution.y[:, -1]
        nfev += solution.nfev
        njev += solution.njev
        nlu += solution.nlu
        progress_parts.append(
            (segment + kept_time / FIRING_INTERVAL_S)
            / SUBSTEPS_PER_SEQUENCE_UNIT
        )
        time_parts.append(segment * FIRING_INTERVAL_S + kept_time)
        command_parts.append(
            np.full(kept_time.size, lead * theta_command)
        )
        position_parts.append(solution.y[5, :-1])
        if (segment + 1) % 32 == 0:
            print(
                f'LuGre edges {segment + 1}/{theta_commands.size}; '
                f'elapsed {time.perf_counter() - started:.1f} s',
                flush=True,
            )

    result = assemble_result(
        progress_parts, time_parts, command_parts, position_parts
    )
    result.update({
        'elapsed_s': time.perf_counter() - started,
        'nfev': nfev,
        'njev': njev,
        'nlu': nlu,
        'final_state': state,
    })
    return result


def second_order_state_space(
    mass: np.ndarray,
    damping: np.ndarray,
    stiffness: np.ndarray,
    command: np.ndarray,
) -> StateSpace:
    n = mass.shape[0]
    inverse_mass = np.linalg.inv(mass)
    system = np.block([
        [np.zeros((n, n)), np.eye(n)],
        [-inverse_mass @ stiffness, -inverse_mass @ damping],
    ])
    input_vector = np.concatenate([
        np.zeros(n), inverse_mass @ command,
    ]).reshape(-1, 1)
    output = np.zeros((1, 2 * n))
    output[0, 5] = 1.0
    return StateSpace(system, input_vector, output, np.zeros((1, 1)))


def simulate_linear(
    system: StateSpace,
    theta_commands: np.ndarray,
    lead: float,
) -> dict[str, np.ndarray]:
    local_time = segment_grid()
    kept_time = local_time[:-1]
    state = np.zeros(system.A.shape[0])
    progress_parts: list[np.ndarray] = []
    time_parts: list[np.ndarray] = []
    command_parts: list[np.ndarray] = []
    position_parts: list[np.ndarray] = []

    for segment, theta_command in enumerate(theta_commands):
        values = np.full(local_time.size, theta_command)
        _, output, states = lsim(
            system, U=values, T=local_time, X0=state
        )
        state = states[-1]
        progress_parts.append(
            (segment + kept_time / FIRING_INTERVAL_S)
            / SUBSTEPS_PER_SEQUENCE_UNIT
        )
        time_parts.append(segment * FIRING_INTERVAL_S + kept_time)
        command_parts.append(
            np.full(kept_time.size, lead * theta_command)
        )
        position_parts.append(np.asarray(output[:-1]).reshape(-1))

    return assemble_result(
        progress_parts, time_parts, command_parts, position_parts
    )


def simulate_full_detent_frictionless(
    model: LuGreModelRev42,
    theta_commands: np.ndarray,
    lead: float,
) -> dict[str, np.ndarray | float | int]:
    '''Frictionless dynamics with the exact periodic detent restoring torque.'''
    local_time = segment_grid()
    kept_time = local_time[:-1]
    state = np.zeros(2 * N_Q)
    progress_parts: list[np.ndarray] = []
    time_parts: list[np.ndarray] = []
    command_parts: list[np.ndarray] = []
    position_parts: list[np.ndarray] = []
    nfev = njev = nlu = 0
    started = time.perf_counter()

    def rhs(y: np.ndarray, theta_command: float) -> np.ndarray:
        q = y[:N_Q]
        velocity = y[N_Q:]
        detent = np.zeros(N_Q)
        detent[0] = model.p['T_d'] * np.sin(
            4.0 * model.p['N_r'] * q[0]
        )
        acceleration = model.mass_inverse @ (
            model.command * theta_command
            - model.damping @ velocity
            - model.stiffness @ q
            - detent
        )
        return np.concatenate([velocity, acceleration])

    def jacobian(y: np.ndarray) -> np.ndarray:
        position_block = -model.stiffness.copy()
        position_block[0, 0] -= (
            4.0 * model.p['N_r'] * model.p['T_d']
            * np.cos(4.0 * model.p['N_r'] * y[0])
        )
        return np.block([
            [np.zeros((N_Q, N_Q)), np.eye(N_Q)],
            [
                model.mass_inverse @ position_block,
                -model.mass_inverse @ model.damping,
            ],
        ])

    for segment, theta_command in enumerate(theta_commands):
        solution = solve_ivp(
            lambda _t, y, command=theta_command: rhs(y, command),
            (0.0, FIRING_INTERVAL_S),
            state,
            method='Radau',
            jac=lambda _t, y: jacobian(y),
            t_eval=local_time,
            rtol=RTOL,
            atol=ATOL[:2 * N_Q],
        )
        if not solution.success:
            raise RuntimeError(
                f'Full-detent frictionless integration failed at edge '
                f'{segment}: {solution.message}'
            )
        state = solution.y[:, -1]
        nfev += solution.nfev
        njev += solution.njev
        nlu += solution.nlu
        progress_parts.append(
            (segment + kept_time / FIRING_INTERVAL_S)
            / SUBSTEPS_PER_SEQUENCE_UNIT
        )
        time_parts.append(segment * FIRING_INTERVAL_S + kept_time)
        command_parts.append(
            np.full(kept_time.size, lead * theta_command)
        )
        position_parts.append(solution.y[5, :-1])
        if (segment + 1) % 32 == 0:
            print(
                f'Full-detent frictionless edges '
                f'{segment + 1}/{theta_commands.size}; '
                f'elapsed {time.perf_counter() - started:.1f} s',
                flush=True,
            )

    result = assemble_result(
        progress_parts, time_parts, command_parts, position_parts
    )
    result.update({
        'elapsed_s': time.perf_counter() - started,
        'nfev': nfev,
        'njev': njev,
        'nlu': nlu,
    })
    return result


def baseline_systems(
    include_detent_tangent: bool = True,
) -> tuple[StateSpace, StateSpace, dict[str, str]]:
    newton = load_module(
        'newton_rev4',
        REV4 / 'scripts' / 'build_bode_rev4.py',
    )
    newton_parameters = newton.load_parameters()
    mass_n, damping_n, stiffness_n, inputs_n = newton.build_matrices(
        newton_parameters
    )
    k_d_newton = (
        4.0 * newton_parameters['N_r'] * newton_parameters['T_d']
    )
    if not include_detent_tangent:
        stiffness_n[0, 0] -= k_d_newton
    command_column = newton.INPUT_LABELS.index('theta_cmd')
    command_n = inputs_n[:, command_column].copy()
    command_n[0] = (
        newton_parameters['N_r'] * newton_parameters['T_hold']
    )
    newton_system = second_order_state_space(
        mass_n, damping_n, stiffness_n, command_n
    )

    lagrange = load_module(
        'lagrange_rev4',
        REV4 / 'Lagrange Derivation' / 'scripts'
        / 'build_bode_lagrange_frictionless.py',
    )
    lagrange_parameters = lagrange.load_parameters()
    mass_l, damping_l, stiffness_l, command_l = (
        lagrange.build_lagrange_matrices(lagrange_parameters)
    )
    k_d_lagrange = (
        4.0 * lagrange_parameters['N_r'] * lagrange_parameters['T_d']
    )
    if not include_detent_tangent:
        stiffness_l[0, 0] -= k_d_lagrange
    lagrange_system = second_order_state_space(
        mass_l, damping_l, stiffness_l, command_l
    )
    if include_detent_tangent:
        conventions = {
            'newton': (
                'documented fixed-frame detent: k_d in K[0,0], '
                'G[0]=k_EM'
            ),
            'lagrange': 'G[0]=k_EM with grounded k_d in K[0,0]',
        }
    else:
        conventions = {
            'newton': 'k_d=0 in K, G[0]=k_EM',
            'lagrange': 'k_d=0 in K, G[0]=k_EM',
        }
    return newton_system, lagrange_system, conventions


def kd_in_command_system() -> tuple[StateSpace, str]:
    '''Diagnostic model with the detent tangent tied to theta_cmd.'''
    newton = load_module(
        'newton_rev4_kd_in_command',
        REV4 / 'scripts' / 'build_bode_rev4.py',
    )
    parameters = newton.load_parameters()
    mass, damping, stiffness, inputs = newton.build_matrices(parameters)
    k_em = parameters['N_r'] * parameters['T_hold']
    k_d = 4.0 * parameters['N_r'] * parameters['T_d']
    stiffness[0, 0] -= k_d
    command_column = newton.INPUT_LABELS.index('theta_cmd')
    command = inputs[:, command_column].copy()
    command[0] = k_em + k_d
    system = second_order_state_space(
        mass, damping, stiffness, command
    )
    convention = (
        'diagnostic moving-detent well: k_d removed from K[0,0] and '
        'added to G[0], so G[0]=k_EM+k_d'
    )
    return system, convention


def metrics(result: dict[str, np.ndarray]) -> dict[str, float]:
    return {
        'peak_absolute_error_um': float(
            np.max(np.abs(result['error_m'])) * 1.0e6
        ),
        'final_error_nm': float(result['error_m'][-1] * 1.0e9),
        'maximum_stage_position_um': float(
            np.max(result['x_n_m']) * 1.0e6
        ),
        'minimum_stage_position_um': float(
            np.min(result['x_n_m']) * 1.0e6
        ),
    }


def load_cached_lugre(
    travel_commands: np.ndarray,
) -> dict[str, np.ndarray | float | int] | None:
    '''Reuse the already integrated nonlinear trajectory when it is exact.'''
    if not DATA.exists():
        return None
    try:
        with np.load(DATA) as payload:
            if 'command_edges_m' not in payload.files:
                return None
            if not np.array_equal(
                payload['command_edges_m'], travel_commands
            ):
                return None
            result: dict[str, np.ndarray | float | int] = {
                'progress': payload['progress'].copy(),
                'time_s': payload['time_s'].copy(),
                'command_travel_m': payload['command_travel_m'].copy(),
                'x_n_m': payload['lugre_x_n_m'].copy(),
                'error_m': payload['lugre_error_m'].copy(),
                'elapsed_s': 0.0,
                'nfev': 0,
                'njev': 0,
                'nlu': 0,
            }
    except (BadZipFile, EOFError, OSError, ValueError):
        return None
    if SUMMARY.exists():
        previous = json.loads(SUMMARY.read_text(encoding='utf-8'))
        result['elapsed_s'] = float(previous.get('solver_elapsed_s', 0.0))
        result['nfev'] = int(previous.get('solver_nfev', 0))
        result['njev'] = int(previous.get('solver_njev', 0))
        result['nlu'] = int(previous.get('solver_nlu', 0))
    return result


def load_cached_linear(
    prefix: str,
    travel_commands: np.ndarray,
) -> dict[str, np.ndarray] | None:
    '''Reuse a validated linear trajectory on the identical command grid.'''
    if not DATA.exists():
        return None
    position_key = f'{prefix}_x_n_m'
    error_key = f'{prefix}_error_m'
    try:
        with np.load(DATA) as payload:
            required = {
                'command_edges_m', 'progress', 'time_s',
                'command_travel_m', position_key, error_key,
            }
            if not required.issubset(payload.files):
                return None
            if not np.array_equal(
                payload['command_edges_m'], travel_commands
            ):
                return None
            return {
                'progress': payload['progress'].copy(),
                'time_s': payload['time_s'].copy(),
                'command_travel_m': payload['command_travel_m'].copy(),
                'x_n_m': payload[position_key].copy(),
                'error_m': payload[error_key].copy(),
            }
    except (BadZipFile, EOFError, OSError, ValueError):
        return None


def render(
    lugre: dict[str, np.ndarray],
    newton: dict[str, np.ndarray],
    lagrange: dict[str, np.ndarray],
    *,
    kd_in_command: dict[str, np.ndarray] | None = None,
    full_detent: dict[str, np.ndarray] | None = None,
    figure_path: Path = FIGURE,
    baseline_suffix: str = '',
    title: str = (
        r'Settled stepping montage: 0.300 $\mu$m increments, '
        '250 ms per edge'
    ),
    configuration_note: str = (
        r'LuGre: exact periodic $T_d\sin(4N_r\theta_m)$; '
        r'$k_d=4N_rT_d$ is its local tangent, not an added spring'
        '\n'
        r'Frictionless Newton/Lagrange: $k_d$ grounded in $K_{11}$; '
        r'command column $B_{u,1}$ uses $k_{EM}$ only'
    ),
) -> None:
    fig, axes = plt.subplots(
        2, 2, figsize=(14.0, 8.2),
        gridspec_kw={'width_ratios': [1.65, 1.0]},
    )
    ax_position, ax_first_position = axes[0]
    ax_error, ax_first_error = axes[1]
    progress = lugre['progress']
    styles = [
        (lugre, '#c0392b', '-', 1.15, 'Rev 4.2 parallel LuGre'),
        (
            newton, '#2b6cb0', '-', 1.25,
            f'Newton frictionless baseline{baseline_suffix}',
        ),
        (
            lagrange, '#2f855a', '--', 1.05,
            f'Lagrange frictionless baseline{baseline_suffix}',
        ),
    ]
    if kd_in_command is not None:
        styles.append(
            (
                kd_in_command, '#7b4ab5', '-.', 1.2,
                r'Frictionless diagnostic ($k_d$ in $G$)',
            )
        )
    if full_detent is not None:
        styles.append(
            (
                full_detent, '#d97706', '-.', 1.25,
                'Frictionless full periodic detent (Newton/Lagrange)',
            )
        )
    ax_position.plot(
        progress, lugre['command_travel_m'] * 1.0e6,
        color='#777777', linestyle='--', linewidth=0.9,
        label='commanded travel',
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
    ax_position.set_ylabel(r'Stage position ($\mu$m)')
    ax_error.set_ylabel(r'Tracking error ($\mu$m)')
    ax_error.set_xlabel('Sequence progress (nominal full-step units)')
    ax_position.set_title('Complete 256-edge sequence')
    ax_error.set_title('Complete-sequence tracking error')
    ax_position.legend(loc='best', framealpha=0.92)

    first_count = int(round(FIRING_INTERVAL_S / OUTPUT_DT_S))
    first = slice(0, first_count)
    first_time_ms = lugre['time_s'][first] * 1.0e3
    ax_first_position.plot(
        first_time_ms, lugre['command_travel_m'][first] * 1.0e6,
        color='#777777', linestyle='--', linewidth=0.9,
    )
    for result, color, linestyle, linewidth, _label in styles:
        ax_first_position.plot(
            first_time_ms, result['x_n_m'][first] * 1.0e6,
            color=color, linestyle=linestyle, linewidth=linewidth,
        )
        ax_first_error.plot(
            first_time_ms, result['error_m'][first] * 1.0e6,
            color=color, linestyle=linestyle, linewidth=linewidth,
        )
    ax_first_error.axhline(
        0.0, color='#333333', linestyle=':', linewidth=0.7
    )
    for axis in (ax_first_position, ax_first_error):
        axis.grid(True, color='#cccccc', linewidth=0.45)
        axis.set_xlim(0.0, FIRING_INTERVAL_S * 1.0e3)
    ax_first_position.set_ylabel(r'Stage position ($\mu$m)')
    ax_first_error.set_ylabel(r'Tracking error ($\mu$m)')
    ax_first_error.set_xlabel('Time after first edge (ms)')
    ax_first_position.set_title(r'First 0.300 $\mu$m edge from rest')
    ax_first_error.set_title('First-edge settling error')

    fig.suptitle(title, y=0.985)
    fig.text(
        0.5, 0.943, configuration_note,
        ha='center', va='top', fontsize=9.0, color='#263238',
        bbox={
            'boxstyle': 'round,pad=0.35',
            'facecolor': '#f4f7f8',
            'edgecolor': '#9aa7ad',
            'linewidth': 0.7,
        },
    )
    plot_top = 0.865 if configuration_note.count('\n') >= 2 else 0.895
    fig.tight_layout(rect=[0.0, 0.0, 1.0, plot_top])
    temporary_figure = figure_path.with_name(
        f'{figure_path.stem}.tmp{figure_path.suffix}'
    )
    fig.savefig(temporary_figure, dpi=150)
    plt.close(fig)
    temporary_figure.replace(figure_path)


def main() -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    model = LuGreModelRev42(enforce_interface_power=False)
    lead = model.p['L'] / (2.0 * np.pi)
    travel_commands = command_travel_edges()
    theta_commands = travel_commands / lead
    newton_system, lagrange_system, conventions = baseline_systems()
    newton_kd0_system, lagrange_kd0_system, kd0_conventions = (
        baseline_systems(include_detent_tangent=False)
    )
    kd_in_command_model, kd_in_command_convention = (
        kd_in_command_system()
    )

    lugre = load_cached_lugre(travel_commands)
    if lugre is None:
        print(
            f'Running {theta_commands.size} settled 0.300 um LuGre edges...',
            flush=True,
        )
        lugre = simulate_lugre(model, theta_commands, lead)
    else:
        print(
            'Reusing validated 0.300 um LuGre trajectory; '
            'rebuilding corrected derivation baselines.',
            flush=True,
        )
    baseline_jobs = (
        ('newton', 'Newton frictionless baseline', newton_system),
        ('lagrange', 'Lagrange frictionless baseline', lagrange_system),
        (
            'newton_kd0', 'Newton frictionless baseline with k_d=0',
            newton_kd0_system,
        ),
        (
            'lagrange_kd0', 'Lagrange frictionless baseline with k_d=0',
            lagrange_kd0_system,
        ),
        (
            'kd_in_command', 'frictionless diagnostic with k_d in G',
            kd_in_command_model,
        ),
    )
    baseline_results: dict[str, dict[str, np.ndarray]] = {}
    for prefix, label, system in baseline_jobs:
        result = load_cached_linear(prefix, travel_commands)
        if result is None:
            print(f'Running {label}...', flush=True)
            result = simulate_linear(system, theta_commands, lead)
        else:
            print(f'Reusing validated {label}.', flush=True)
        baseline_results[prefix] = result
    newton = baseline_results['newton']
    lagrange = baseline_results['lagrange']
    newton_kd0 = baseline_results['newton_kd0']
    lagrange_kd0 = baseline_results['lagrange_kd0']
    kd_in_command = baseline_results['kd_in_command']

    full_detent = load_cached_linear(
        'full_detent_frictionless', travel_commands
    )
    if full_detent is None:
        print(
            'Running frictionless baseline with full periodic detent...',
            flush=True,
        )
        full_detent = simulate_full_detent_frictionless(
            model, theta_commands, lead
        )
    else:
        print(
            'Reusing validated frictionless full-periodic-detent baseline.',
            flush=True,
        )

    render(
        lugre, newton, lagrange,
        full_detent=full_detent,
        configuration_note=(
            r'LuGre: exact $T_d\sin(4N_r\theta_m)$ plus LuGre friction; '
            r'no separate $k_d$ spring'
            '\n'
            r'Linear frictionless Newton/Lagrange: '
            r'$k_d=4N_rT_d$ grounded in $K_{11}$'
            '\n'
            r'Full-detent frictionless overlay: replace $k_d\theta_m$ by '
            r'$T_d\sin(4N_r\theta_m)$; no LuGre friction'
        ),
    )
    render(
        lugre, newton_kd0, lagrange_kd0,
        figure_path=KD0_FIGURE,
        baseline_suffix=r' ($k_d=0$)',
        title=(
            r'Settled stepping montage: frictionless $k_d=0$, '
            r'0.300 $\mu$m increments'
        ),
        configuration_note=(
            r'LuGre unchanged: exact periodic $T_d\sin(4N_r\theta_m)$; '
            r'$k_d=4N_rT_d$ is its local tangent, not an added spring'
            '\n'
            r'Frictionless Newton/Lagrange diagnostic: $k_d=0$ in '
            r'$K_{11}$; command column $B_{u,1}$ uses $k_{EM}$ only'
        ),
    )
    temporary_data = DATA.with_name(f'{DATA.stem}.tmp{DATA.suffix}')
    np.savez_compressed(
        temporary_data,
        progress=lugre['progress'],
        time_s=lugre['time_s'],
        command_travel_m=lugre['command_travel_m'],
        lugre_x_n_m=lugre['x_n_m'],
        lugre_error_m=lugre['error_m'],
        newton_x_n_m=newton['x_n_m'],
        newton_error_m=newton['error_m'],
        lagrange_x_n_m=lagrange['x_n_m'],
        lagrange_error_m=lagrange['error_m'],
        newton_kd0_x_n_m=newton_kd0['x_n_m'],
        newton_kd0_error_m=newton_kd0['error_m'],
        lagrange_kd0_x_n_m=lagrange_kd0['x_n_m'],
        lagrange_kd0_error_m=lagrange_kd0['error_m'],
        kd_in_command_x_n_m=kd_in_command['x_n_m'],
        kd_in_command_error_m=kd_in_command['error_m'],
        full_detent_frictionless_x_n_m=full_detent['x_n_m'],
        full_detent_frictionless_error_m=full_detent['error_m'],
        command_edges_m=travel_commands,
    )
    temporary_data.replace(DATA)

    summary = {
        'model': 'Rev 4.2 smoothed parallel LuGre settled stepping',
        'move_sequence': list(MOVES),
        'step_travel_um': STEP_TRAVEL_M * 1.0e6,
        'substeps_per_sequence_unit': SUBSTEPS_PER_SEQUENCE_UNIT,
        'command_edges': int(theta_commands.size),
        'nominal_sequence_unit_travel_um': (
            SUBSTEPS_PER_SEQUENCE_UNIT * STEP_TRAVEL_M * 1.0e6
        ),
        'firing_interval_s': FIRING_INTERVAL_S,
        'total_duration_s': float(theta_commands.size * FIRING_INTERVAL_S),
        'output_dt_s': OUTPUT_DT_S,
        'solver': 'piecewise Radau with analytical Rev 4.2 Jacobian',
        'solver_rtol': RTOL,
        'solver_elapsed_s': float(lugre['elapsed_s']),
        'solver_nfev': int(lugre['nfev']),
        'solver_njev': int(lugre['njev']),
        'solver_nlu': int(lugre['nlu']),
        'baseline_conventions': conventions,
        'kd0_baseline_conventions': kd0_conventions,
        'kd_in_command_convention': kd_in_command_convention,
        'full_detent_frictionless_convention': (
            'otherwise frictionless Newton/Lagrange dynamics with k_d '
            'removed from K[0,0] and replaced by the exact fixed-frame '
            'T_d*sin(4*N_r*theta_m) restoring torque'
        ),
        'detent_convention': (
            'fixed rotor-stator tooth well: k_d is grounded in K[0,0] '
            'and absent from the theta_cmd input column'
        ),
        'results': {
            'lugre': metrics(lugre),
            'newton_frictionless': metrics(newton),
            'lagrange_frictionless': metrics(lagrange),
            'newton_frictionless_kd0': metrics(newton_kd0),
            'lagrange_frictionless_kd0': metrics(lagrange_kd0),
            'frictionless_kd_in_command': metrics(kd_in_command),
            'frictionless_full_periodic_detent': metrics(full_detent),
        },
        'maximum_newton_vs_lagrange_stage_difference_um': float(
            np.max(np.abs(newton['x_n_m'] - lagrange['x_n_m'])) * 1.0e6
        ),
        'maximum_newton_vs_lagrange_kd0_stage_difference_um': float(
            np.max(
                np.abs(newton_kd0['x_n_m'] - lagrange_kd0['x_n_m'])
            ) * 1.0e6
        ),
        'maximum_fixed_kd_vs_kd0_stage_difference_um': float(
            np.max(np.abs(newton['x_n_m'] - newton_kd0['x_n_m']))
            * 1.0e6
        ),
        'maximum_fixed_kd_vs_kd_in_command_stage_difference_um': float(
            np.max(np.abs(newton['x_n_m'] - kd_in_command['x_n_m']))
            * 1.0e6
        ),
        'maximum_linear_tangent_vs_full_detent_stage_difference_um': float(
            np.max(np.abs(newton['x_n_m'] - full_detent['x_n_m']))
            * 1.0e6
        ),
        'figure': FIGURE.name,
        'kd0_figure': KD0_FIGURE.name,
        'data': str(DATA.relative_to(ROOT)),
    }
    SUMMARY.write_text(
        json.dumps(summary, indent=2) + '\n', encoding='utf-8'
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()

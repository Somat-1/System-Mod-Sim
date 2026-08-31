#!/usr/bin/env python3
'''Random-odd periodic multisine experiment for the Rev 4.2 LuGre model.

The excitation is defined on the exact 4 s DFT grid.  Only odd DFT bins are
excited; reproducible random phases provide M independent realizations and
nearby odd lines are intentionally left empty as nonlinear-distortion
detection lines.  Each realization is integrated one period at a time until
its full 15-state waveform is periodic, followed by three retained periods.
'''

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

from lugre_model_rev42 import N_Q, N_STATES, PORTS, LuGreModelRev42


BASE_PERIOD_S = 4.0
FS_HZ = 65536.0
N_SAMPLES = 262144
F_LO_HZ = 0.25
F_HI_HZ = 8000.0
N_EXCITED_LINES = 370
N_REALIZATIONS = 7
N_RETAINED_PERIODS = 3
AMPLITUDE_PERCENT = np.array(
    [5.0, 10.0, 20.0, 35.0, 50.0, 65.0, 80.0, 90.0,
     100.0, 110.0, 125.0, 150.0,
     175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,
     550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
)
SOLVERS = ('RK45', 'Radau')
RTOL = 1.0e-6
CONVERGENCE_TOL = 2.0e-3
MAX_CONVERGENCE_PERIODS = 12
MAX_STEP_S = 2.0 / FS_HZ
DESIGN_SEED = 42042
PHASE_SEED_BASE = 42100

assert N_SAMPLES == int(BASE_PERIOD_S * FS_HZ)


def breakaway_command_rms(model: LuGreModelRev42) -> float:
    '''Equivalent command angle for guideway Fs/sigma0 bristle deflection.'''
    p = model.p
    lead = p['L'] / (2.0 * np.pi)
    return (p['Fs_way'] / p['sigma0_way']) / lead


def _nearest_free_odd(target: float, used: set[int], k_max: int) -> int:
    center = int(np.clip(round(target), 1, k_max))
    if center % 2 == 0:
        center += 1 if center < k_max else -1
    for radius in range(0, k_max + 1, 2):
        for candidate in (center - radius, center + radius):
            if 1 <= candidate <= k_max and candidate % 2 == 1 and candidate not in used:
                return candidate
    raise RuntimeError('Unable to allocate requested odd multisine bins')


def design_random_odd_bins(
    seed: int = DESIGN_SEED,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Return excited bins, omitted odd detection bins, and omission counts.

    Excited lines are quasi-logarithmic to resolve the low-frequency range
    while reaching the highest admissible odd bin (7999.75 Hz).  For each
    excited line, one to three otherwise-unused odd bins are randomly reserved
    as detection lines.  This implements the requested random-odd omission
    without contaminating the 370 driven lines.
    '''
    k_min = int(round(F_LO_HZ * BASE_PERIOD_S))
    k_max = int(np.floor(F_HI_HZ * BASE_PERIOD_S))
    if k_min % 2 == 0:
        k_min += 1
    if k_max % 2 == 0:
        k_max -= 1

    targets = np.geomspace(k_min, k_max, N_EXCITED_LINES)
    used: set[int] = set()
    for target in targets:
        used.add(_nearest_free_odd(float(target), used, k_max))
    excited = np.array(sorted(used), dtype=np.int64)
    if len(excited) != N_EXCITED_LINES:
        raise AssertionError('Excited-line allocation did not produce exactly 370 lines')

    rng = np.random.default_rng(seed)
    omission_counts = rng.integers(1, 4, size=N_EXCITED_LINES)
    candidates = np.setdiff1d(
        np.arange(k_min, k_max + 1, 2, dtype=np.int64), excited,
        assume_unique=True,
    )
    omitted_count = int(np.sum(omission_counts))
    omitted = np.sort(rng.choice(candidates, size=omitted_count, replace=False))
    return excited, omitted, omission_counts


@dataclass
class PeriodicMultisine:
    amplitude_rms: float
    realization: int
    excited_bins: np.ndarray

    def __post_init__(self) -> None:
        rng = np.random.default_rng(PHASE_SEED_BASE + self.realization)
        self.phases = rng.uniform(0.0, 2.0 * np.pi, len(self.excited_bins))
        spectrum = np.zeros(N_SAMPLES // 2 + 1, dtype=np.complex128)
        line_peak = self.amplitude_rms * np.sqrt(2.0 / len(self.excited_bins))
        spectrum[self.excited_bins] = (
            0.5 * N_SAMPLES * line_peak * np.exp(1j * self.phases)
        )
        wave = np.fft.irfft(spectrum, n=N_SAMPLES)
        wave *= self.amplitude_rms / np.sqrt(np.mean(wave**2))
        self.wave = wave
        self.crest_factor = float(np.max(np.abs(wave)) / self.amplitude_rms)
        nodes_t = np.arange(N_SAMPLES + 1, dtype=float) / FS_HZ
        nodes_y = np.concatenate((wave, wave[:1]))
        self._spline = CubicSpline(nodes_t, nodes_y, bc_type='periodic')
        self.sampled_period = np.roll(wave, -1)
        self.input_fft = np.fft.rfft(self.sampled_period)

    def __call__(self, t: float) -> float:
        return float(self._spline(np.remainder(t, BASE_PERIOD_S)))


class FastRev42Dynamics:
    '''Allocation-light equivalent of LuGreModelRev42.rhs/Jacobian.'''

    def __init__(self, model: LuGreModelRev42, signal: PeriodicMultisine):
        self.model, self.signal, self.p = model, signal, model.p
        self.mass_inverse = np.diag(model.mass_inverse)
        self.damping = model.damping
        self.stiffness = model.stiffness
        self.command = model.command
        self.jacobians = np.stack([model.jacobians[port] for port in PORTS])
        self.outer_jacobians = np.stack([
            np.outer(row, row) for row in self.jacobians
        ])
        prefixes = ['F', 'F', 'T']
        self.sigma0 = np.array([self.p[f'sigma0_{port}'] for port in PORTS])
        self.sigma1 = np.array([self.p[f'sigma1_{port}'] for port in PORTS])
        self.sigma2 = np.array([self.p[f'sigma2_{port}'] for port in PORTS])
        self.fc = np.array([
            self.p[f'{prefix}c_{port}']
            for prefix, port in zip(prefixes, PORTS)
        ])
        self.fs = np.array([
            self.p[f'{prefix}s_{port}']
            for prefix, port in zip(prefixes, PORTS)
        ])
        self.vs = np.array([self.p[f'vs_{port}'] for port in PORTS])
        self.epsilon = self.p['smooth_velocity_epsilon']

    def port_terms(self, velocity: np.ndarray, z: np.ndarray):
        v = self.jacobians @ velocity
        exponent = np.exp(-(v / self.vs) ** 2)
        g = self.fc + (self.fs - self.fc) * exponent
        g_prime = (self.fs - self.fc) * exponent * (-2.0 * v / self.vs**2)
        speed = np.sqrt(v**2 + self.epsilon**2)
        speed_prime = v / speed
        decay = self.sigma0 * speed / g
        decay_prime = self.sigma0 * (
            speed_prime / g - speed * g_prime / g**2
        )
        z_dot = v - decay * z
        dzdv = 1.0 - decay_prime * z
        dzdz = -decay
        force = self.sigma0 * z + self.sigma1 * z_dot + self.sigma2 * v
        dfdv = self.sigma1 * dzdv + self.sigma2
        dfdz = self.sigma0 + self.sigma1 * dzdz
        return force, z_dot, dfdv, dfdz, dzdv, dzdz

    def rhs(self, t: float, state: np.ndarray) -> np.ndarray:
        q = state[:N_Q]
        velocity = state[N_Q:2 * N_Q]
        z = state[2 * N_Q:]
        force, z_dot, *_ = self.port_terms(velocity, z)
        reaction = -(force @ self.jacobians)
        detent = np.zeros(N_Q)
        detent[0] = (
            self.p['T_d']
            * np.sin(4.0 * self.p['N_r'] * q[0])
        )
        acceleration = self.mass_inverse * (
            self.command * self.signal(t)
            + reaction
            - self.damping @ velocity
            - self.stiffness @ q
            - detent
        )
        derivative = np.empty(N_STATES)
        derivative[:N_Q] = velocity
        derivative[N_Q:2 * N_Q] = acceleration
        derivative[2 * N_Q:] = z_dot
        return derivative

    def jacobian(self, _t: float, state: np.ndarray) -> np.ndarray:
        velocity = state[N_Q:2 * N_Q]
        z = state[2 * N_Q:]
        _, _, dfdv, dfdz, dzdv, dzdz = self.port_terms(velocity, z)
        system = np.zeros((N_STATES, N_STATES))
        system[:N_Q, N_Q:2 * N_Q] = np.eye(N_Q)
        position_block = -self.stiffness.copy()
        position_block[0, 0] -= (
            4.0 * self.p['N_r'] * self.p['T_d']
            * np.cos(4.0 * self.p['N_r'] * state[0])
        )
        velocity_block = -self.damping.copy()
        for index, jacobian in enumerate(self.jacobians):
            velocity_block -= dfdv[index] * self.outer_jacobians[index]
            system[N_Q:2 * N_Q, 2 * N_Q + index] = (
                self.mass_inverse * (-jacobian * dfdz[index])
            )
            system[2 * N_Q + index, N_Q:2 * N_Q] = (
                dzdv[index] * jacobian
            )
            system[2 * N_Q + index, 2 * N_Q + index] = dzdz[index]
        system[N_Q:2 * N_Q, :N_Q] = self.mass_inverse[:, None] * position_block
        system[N_Q:2 * N_Q, N_Q:2 * N_Q] = (
            self.mass_inverse[:, None] * velocity_block
        )
        return system


def state_scales(amplitude_rms: float, p: dict[str, float]) -> np.ndarray:
    lead = p['L'] / (2.0 * np.pi)
    omega = 2.0 * np.pi * F_HI_HZ
    angular = amplitude_rms
    axial = amplitude_rms * lead
    scales = np.array([
        angular, angular, angular, angular, axial, axial,
        angular * omega, angular * omega, angular * omega, angular * omega,
        axial * omega, axial * omega,
        p['Fs_way'] / p['sigma0_way'],
        p['Fs_nut'] / p['sigma0_nut'],
        p['Ts_sb'] / p['sigma0_sb'],
    ])
    return np.maximum(scales, 1.0e-15)


def solver_atol(amplitude_rms: float, p: dict[str, float]) -> np.ndarray:
    return np.maximum(1.0e-7 * state_scales(amplitude_rms, p), 1.0e-14)


def periodicity_metric(
    current: np.ndarray,
    previous: np.ndarray,
    amplitude_rms: float,
    p: dict[str, float],
) -> tuple[float, np.ndarray]:
    rms_current = np.sqrt(np.mean(current**2, axis=1))
    rms_difference = np.sqrt(np.mean((current - previous) ** 2, axis=1))
    floor = 1.0e-4 * state_scales(amplitude_rms, p)
    relative = rms_difference / np.maximum(rms_current, floor)
    return float(np.max(relative)), relative


def _port_metrics(
    model: LuGreModelRev42,
    states: np.ndarray,
) -> tuple[dict[str, float], dict[str, float]]:
    p = model.p
    velocity = states[N_Q:2 * N_Q]
    saturation: dict[str, float] = {}
    slip_fraction: dict[str, float] = {}
    for index, port in enumerate(PORTS):
        v = model.jacobians[port] @ velocity
        z = states[2 * N_Q + index]
        prefix = 'T' if port == 'sb' else 'F'
        static_z = p[f'{prefix}s_{port}'] / p[f'sigma0_{port}']
        saturation[port] = float(np.max(np.abs(z)) / static_z)
        fc = p[f'{prefix}c_{port}']
        fs = p[f'{prefix}s_{port}']
        vs = p[f'vs_{port}']
        g = fc + (fs - fc) * np.exp(-(v / vs) ** 2)
        speed = np.sqrt(v**2 + p['smooth_velocity_epsilon'] ** 2)
        z_dot = v - p[f'sigma0_{port}'] * speed / g * z
        moving = np.abs(v) >= vs
        steady_slip = moving & (np.abs(z_dot) <= 0.1 * np.maximum(np.abs(v), 1.0e-30))
        slip_fraction[port] = float(np.mean(steady_slip))
    return saturation, slip_fraction


def run_case(payload: dict) -> dict:
    solver = str(payload['solver'])
    amplitude_percent = float(payload['amplitude_percent'])
    realization = int(payload['realization'])
    max_convergence_periods = int(
        payload.get('max_convergence_periods', MAX_CONVERGENCE_PERIODS)
    )
    model = LuGreModelRev42(enforce_interface_power=False)
    amplitude_breakaway = breakaway_command_rms(model)
    amplitude_rms = amplitude_percent * 0.01 * amplitude_breakaway
    excited, omitted, omission_counts = design_random_odd_bins()
    signal = PeriodicMultisine(amplitude_rms, realization, excited)
    dynamics = FastRev42Dynamics(model, signal)
    atol = solver_atol(amplitude_rms, model.p)
    initial = np.zeros(N_STATES)
    previous: np.ndarray | None = None
    retained_fft: list[np.ndarray] = []
    omitted_fft: list[np.ndarray] = []
    convergence_history: list[float] = []
    component_history: list[list[float]] = []
    saturation_max = {port: 0.0 for port in PORTS}
    slip_sum = {port: 0.0 for port in PORTS}
    slip_periods = 0
    converged = False
    forced_after_limit = False
    convergence_period: int | None = None
    nfev = njev = nlu = 0
    wall_start = time.perf_counter()

    rhs = dynamics.rhs
    jac = dynamics.jacobian

    max_total = max_convergence_periods + N_RETAINED_PERIODS
    for period in range(1, max_total + 1):
        t_start = (period - 1) * BASE_PERIOD_S
        t_eval = t_start + np.arange(1, N_SAMPLES + 1, dtype=float) / FS_HZ
        kwargs = dict(
            method=solver,
            t_eval=t_eval,
            rtol=RTOL,
            atol=atol,
            max_step=MAX_STEP_S,
        )
        if solver == 'Radau':
            kwargs['jac'] = jac
        solution = solve_ivp(
            rhs, (t_start, t_start + BASE_PERIOD_S), initial, **kwargs
        )
        if not solution.success:
            raise RuntimeError(
                f'{solver} failed at {amplitude_percent:g}% Fs, '
                f'realization {realization}: {solution.message}'
            )
        current = solution.y
        initial = current[:, -1].copy()
        nfev += int(solution.nfev)
        njev += int(solution.njev)
        nlu += int(solution.nlu)

        saturation, slip_fraction = _port_metrics(model, current)
        for port in PORTS:
            saturation_max[port] = max(saturation_max[port], saturation[port])

        if previous is not None and convergence_period is None:
            metric, components = periodicity_metric(
                current, previous, amplitude_rms, model.p
            )
            convergence_history.append(metric)
            component_history.append(components.tolist())
            if metric <= CONVERGENCE_TOL:
                converged = True
                convergence_period = period

        if convergence_period is None and period >= max_convergence_periods:
            convergence_period = period
            forced_after_limit = True

        if convergence_period is not None and period > convergence_period:
            output_fft = np.fft.rfft(current[5])
            retained_fft.append(output_fft[excited])
            omitted_fft.append(output_fft[omitted])
            for port in PORTS:
                slip_sum[port] += slip_fraction[port]
            slip_periods += 1
            if len(retained_fft) == N_RETAINED_PERIODS:
                break
        previous = current

    wall_s = time.perf_counter() - wall_start
    y_excited = np.mean(np.stack(retained_fft), axis=0)
    y_omitted = np.mean(np.stack(omitted_fft), axis=0)
    u_excited = signal.input_fft[excited]
    transfer = y_excited / u_excited
    stats = {
        'solver': solver,
        'amplitude_percent_Fs': amplitude_percent,
        'amplitude_rms_rad': amplitude_rms,
        'amplitude_breakaway_rms_rad': amplitude_breakaway,
        'realization': realization,
        'wall_s': wall_s,
        'converged': converged,
        'forced_after_limit': forced_after_limit,
        'convergence_period': convergence_period,
        'total_periods': period,
        'retained_periods': len(retained_fft),
        'convergence_tolerance': CONVERGENCE_TOL,
        'final_periodicity_metric': (
            convergence_history[-1] if convergence_history else None
        ),
        'nfev': nfev,
        'njev': njev,
        'nlu': nlu,
        'rtol': RTOL,
        'atol_min': float(np.min(atol)),
        'atol_max': float(np.max(atol)),
        'max_step_s': MAX_STEP_S,
        'crest_factor': signal.crest_factor,
        'saturation_max': saturation_max,
        'steady_slip_fraction': {
            port: slip_sum[port] / max(slip_periods, 1) for port in PORTS
        },
        'periodicity_history': convergence_history,
        'periodicity_components': component_history,
    }
    return {
        'solver': solver,
        'amplitude_percent': amplitude_percent,
        'realization': realization,
        'excited_bins': excited,
        'omitted_bins': omitted,
        'omission_counts': omission_counts,
        'frequency_hz': excited / BASE_PERIOD_S,
        'omitted_frequency_hz': omitted / BASE_PERIOD_S,
        'U': u_excited,
        'Y': y_excited,
        'G': transfer,
        'Y_omitted': y_omitted,
        'stats': stats,
    }


def save_case(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        solver=np.array(result['solver']),
        amplitude_percent=np.array(result['amplitude_percent']),
        realization=np.array(result['realization']),
        excited_bins=result['excited_bins'],
        omitted_bins=result['omitted_bins'],
        omission_counts=result['omission_counts'],
        frequency_hz=result['frequency_hz'],
        omitted_frequency_hz=result['omitted_frequency_hz'],
        U=result['U'],
        Y=result['Y'],
        G=result['G'],
        Y_omitted=result['Y_omitted'],
        stats_json=np.array(json.dumps(result['stats'])),
    )


def load_case(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {
            'solver': str(data['solver']),
            'amplitude_percent': float(data['amplitude_percent']),
            'realization': int(data['realization']),
            'excited_bins': data['excited_bins'],
            'omitted_bins': data['omitted_bins'],
            'omission_counts': data['omission_counts'],
            'frequency_hz': data['frequency_hz'],
            'omitted_frequency_hz': data['omitted_frequency_hz'],
            'U': data['U'],
            'Y': data['Y'],
            'G': data['G'],
            'Y_omitted': data['Y_omitted'],
            'stats': json.loads(str(data['stats_json'])),
        }

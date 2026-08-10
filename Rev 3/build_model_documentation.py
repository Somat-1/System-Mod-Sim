#!/usr/bin/env python3
"""Build the model response figures and render both Markdown documents to HTML.

The script deliberately depends only on NumPy and Matplotlib.  It implements a
fixed-step RK4 integrator for nonlinear LuGre and GMS simulations and a compact
Markdown renderer tailored to the two project documents.
"""

from __future__ import annotations

import argparse
import html
import re
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"
DESCRIPTION_MD = ROOT / "ball_screw_stage_dynamic_derivation_v3.md"
DERIVATION_MD = ROOT / "Analytical_derivation_and_responses_v3.md"


# Executable defaults for the Revision 3 two-DOF reduction.
MODEL = {
    "lead": 1.0e-3,
    "rotor_teeth": 50,
    "m_d": 59.0,
    "m_s": 0.60,
    "k_ax": 1.14e7,
    "K_m": 1.20e8,
    # Provisional: retained structural damping, not identified in the source.
    "c_ax": 55.0,
    # Provisional electromagnetic modal damping ratio.  This represents the
    # current-regulator/back-EMF damping missing from the original equations.
    "zeta_m": 0.50,
}


# Revision 3 full-model values.  These are deliberately separate from MODEL:
# MODEL is the validated two-DOF reduction, whereas FULL retains all ten
# coordinates named in the source document.  Values not measured in the source
# are surfaced as highlighted assumptions in the Markdown documents.
FULL = {
    "J_m": 1.20e-6,
    # The source-table placeholder 1.2e-6 kg m^2 is inconsistent with the
    # stated 226 Hz mode.  5e-8 closes the 59 kg reflected-inertia budget.
    "J_c": 5.00e-8,
    "J_s1": 8.15e-8,
    "J_s2": 8.15e-8,
    "J_s3": 8.15e-8,
    "m_b": 0.015,
    "m_e": 0.015,
    "m_f": 0.010,
    "m_n": 0.050,
    "m_stage": 0.550,
    "k_c1": 100.0,
    "k_c2": 100.0,
    "k_theta_a": 211.0,
    "k_theta_b": 211.0,
    # 25 N/um is the closure-consistent bearing assumption discussed in Rev 3.
    "k_brg": 25.0e6,
    "k_sha": 67.0e6,
    "k_shb": 200.0e6,
    "k_mnt": 100.0e6,
    "zeta_internal": 0.010,
}

FULL_DOF_LABELS = (
    r"$\theta_m$", r"$\theta_c$", r"$\theta_{s1}$", r"$\theta_{s2}$",
    r"$\theta_{s3}$", r"$u_b$", r"$u_e$", r"$u_f$", r"$u_n$", r"$x_s$",
)


# Highlighted friction-port values used to make both law comparisons executable.
# sigma0_g is the estimate already quoted in the description; all other values
# below need experimental identification before quantitative use.
FRICTION = {
    "g": {"sigma0": 7.60e5, "sigma1": 3.0, "sigma2": 0.40,
          "F_s": 3.0, "F_c": 2.4, "v_s": 2.5e-4, "delta": 1.0, "C_gms": 5.0e3},
    "n": {"sigma0": 2.00e6, "sigma1": 5.0, "sigma2": 0.25,
          "F_s": 5.0, "F_c": 4.0, "v_s": 2.0e-4, "delta": 1.0, "C_gms": 5.0e3},
    "d": {"sigma0": 1.00e6, "sigma1": 4.0, "sigma2": 0.20,
          "F_s": 2.0, "F_c": 1.5, "v_s": 3.0e-4, "delta": 1.0, "C_gms": 5.0e3},
}


# Four GMS stop elements share each site's aggregate sigma0 and Stribeck
# force.  Opposing stiffness/force fractions create distinct yield distances
# and therefore non-local reversal memory while retaining the LuGre aggregate.
GMS_WEIGHTS = np.array([0.10, 0.20, 0.30, 0.40])
GMS_STIFFNESS_FRACTIONS = np.array([0.40, 0.30, 0.20, 0.10])
GMS_N = GMS_WEIGHTS.size
GMS_CONVERGENCE_DTS = (1.0e-5, 5.0e-6, 2.5e-6)
BODE_FOCUS_MIN_HZ = 100.0
BODE_FOCUS_MAX_HZ = 3000.0


# Nested reversals for the dedicated presliding-memory experiment.  One
# microstep is 1/8 of the already bounded quarter-step command (1/32 full
# step).  Repeated levels create return points without reaching gross sliding.
PRESLIDING_LEVELS = np.array([0, 7, 2, 6, 2, 7, 0, -6, -2, -5, -2, -6, 0], dtype=float)
PRESLIDING_START = 0.005
PRESLIDING_HOLD = 0.010
PRESLIDING_RETURN_PAIRS = ((1, 5), (2, 4), (7, 11), (8, 10))


CASES = OrderedDict([
    ("0", {"label": "Case 0 — frictionless", "sites": (), "friction": "none", "color": "#252525", "ls": "--"}),
    ("A", {"label": "Case A — guideway / LuGre", "sites": ("g",), "friction": "lugre", "color": "#277da1", "ls": "-"}),
    ("A2", {"label": "Case A2 — guideway / GMS", "sites": ("g",), "friction": "gms", "color": "#70b7cf", "ls": "--"}),
    ("B", {"label": "Case B — nut / LuGre", "sites": ("n",), "friction": "lugre", "color": "#e07a15", "ls": "-"}),
    ("B2", {"label": "Case B2 — nut / GMS", "sites": ("n",), "friction": "gms", "color": "#f5b35f", "ls": "--"}),
    ("C", {"label": "Case C — guideway + nut / LuGre", "sites": ("g", "n"), "friction": "lugre", "color": "#218c74", "ls": "-"}),
    ("C2", {"label": "Case C2 — guideway + nut / GMS", "sites": ("g", "n"), "friction": "gms", "color": "#72c9ad", "ls": "--"}),
])

PAIRS = (("A", "A2"), ("B", "B2"), ("C", "C2"))


H = {
    "g": np.array([0.0, 1.0]),
    "n": np.array([1.0, -1.0]),
    "d": np.array([1.0, 0.0]),
}

LUGRE_INDEX = {"g": 4, "n": 5, "d": 6}
GMS_START = {"g": 7, "n": 7 + GMS_N, "d": 7 + 2 * GMS_N}
STATE_SIZE = 7 + 3 * GMS_N


def validate_gms_partition() -> dict[str, object]:
    """Fail the build unless every executed GMS partition closes exactly."""
    if GMS_WEIGHTS.size != GMS_STIFFNESS_FRACTIONS.size:
        raise ValueError("GMS force weights and stiffness fractions must have equal length")
    if np.any(GMS_WEIGHTS <= 0.0) or np.any(GMS_STIFFNESS_FRACTIONS <= 0.0):
        raise ValueError("Every GMS force weight and stiffness fraction must be positive")
    weight_sum = float(np.sum(GMS_WEIGHTS))
    if not np.isclose(weight_sum, 1.0, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"GMS force weights must satisfy sum(nu_i)=1; got {weight_sum:.16g}")

    stiffness_sums: dict[str, float] = {}
    for site, parameters in FRICTION.items():
        element_stiffness = GMS_STIFFNESS_FRACTIONS * parameters["sigma0"]
        stiffness_sum = float(np.sum(element_stiffness))
        if not np.isclose(stiffness_sum, parameters["sigma0"],
                          rtol=1.0e-12, atol=1.0e-9):
            raise ValueError(
                f"GMS site {site} must satisfy sum(k_i)=sigma0; "
                f"got {stiffness_sum:.16g} versus {parameters['sigma0']:.16g}"
            )
        stiffness_sums[site] = stiffness_sum
    return {"weight_sum": weight_sum, "stiffness_sums": stiffness_sums}


def physical_constants() -> dict[str, float]:
    lead = MODEL["lead"]
    teeth = MODEL["rotor_teeth"]
    r = lead / (2.0 * np.pi)
    kappa = 2.0 * np.pi * teeth / lead
    t_max = MODEL["K_m"] * r * r / teeth
    f_max = t_max / r
    full_step = lead / (4.0 * teeth)
    c_m = 2.0 * MODEL["zeta_m"] * np.sqrt(MODEL["K_m"] * MODEL["m_d"])
    return {
        "r": r,
        "kappa": kappa,
        "T_max": t_max,
        "F_max": f_max,
        "c_m": c_m,
        "full_step": full_step,
        "quarter_step": full_step / 4.0,
    }


def full_parameters() -> dict[str, float]:
    """Return the ten-DOF parameters and close the measured axial compliance."""
    p = dict(FULL)
    remaining_compliance = (
        1.0 / MODEL["k_ax"]
        - 1.0 / p["k_brg"]
        - 1.0 / p["k_sha"]
        - 1.0 / p["k_mnt"]
    )
    if remaining_compliance <= 0.0:
        raise ValueError("Full-model axial compliance budget cannot be closed")
    p["k_ball"] = 1.0 / remaining_compliance
    r = physical_constants()["r"]
    p["J_total"] = p["J_m"] + p["J_c"] + p["J_s1"] + p["J_s2"] + p["J_s3"]
    p["m_d_reflected"] = p["J_total"] / r**2
    p["literal_table_m_d"] = (1.2e-6 + 1.2e-6 + p["J_s1"] + p["J_s2"] + p["J_s3"]) / r**2
    return p


def _add_pair(matrix: np.ndarray, i: int, j: int, value: float) -> None:
    matrix[i, i] += value
    matrix[j, j] += value
    matrix[i, j] -= value
    matrix[j, i] -= value


def _pair_damping(stiffness: float, mass_i: float, mass_j: float, zeta: float) -> float:
    reduced_mass = mass_i * mass_j / (mass_i + mass_j)
    return 2.0 * zeta * np.sqrt(stiffness * reduced_mass)


def full_linear_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Assemble the Revision 3 ten-DOF linear model from physical elements."""
    p = full_parameters()
    r = physical_constants()["r"]
    native_masses = np.array([
        p["J_m"], p["J_c"], p["J_s1"], p["J_s2"], p["J_s3"],
        p["m_b"], p["m_e"], p["m_f"], p["m_n"], p["m_stage"],
    ])
    mass = np.diag(native_masses)
    damping = np.zeros((10, 10), dtype=float)
    stiffness = np.zeros((10, 10), dtype=float)
    zeta = p["zeta_internal"]

    k_m_rot = MODEL["K_m"] * r**2
    stiffness[0, 0] += k_m_rot
    # The same damping repair validated in Rev 2, expressed in rotational units.
    damping[0, 0] += 2.0 * MODEL["zeta_m"] * np.sqrt(k_m_rot * p["J_total"])

    rotational_pairs = (
        (0, 1, p["k_c1"]),
        (1, 2, p["k_c2"]),
        (2, 3, p["k_theta_a"]),
        (3, 4, p["k_theta_b"]),
    )
    for i, j, k_value in rotational_pairs:
        _add_pair(stiffness, i, j, k_value)
        _add_pair(damping, i, j, _pair_damping(k_value, native_masses[i], native_masses[j], zeta))

    stiffness[5, 5] += p["k_brg"]
    damping[5, 5] += 2.0 * zeta * np.sqrt(p["k_brg"] * p["m_b"])
    for i, j, k_value in ((5, 6, p["k_sha"]), (6, 7, p["k_shb"]), (8, 9, p["k_mnt"])):
        _add_pair(stiffness, i, j, k_value)
        _add_pair(damping, i, j, _pair_damping(k_value, native_masses[i], native_masses[j], zeta))

    # delta_n = u_n - u_e - r theta_s2; outer products are the virtual-work
    # contribution of one conservative ball-contact element to K and C.
    h_nut = np.zeros(10)
    h_nut[3], h_nut[6], h_nut[8] = -r, -1.0, 1.0
    stiffness += p["k_ball"] * np.outer(h_nut, h_nut)
    relative_mass = 1.0 / (r**2 / p["J_s2"] + 1.0 / p["m_e"] + 1.0 / p["m_n"])
    c_ball = 2.0 * zeta * np.sqrt(p["k_ball"] * relative_mass)
    damping += c_ball * np.outer(h_nut, h_nut)
    p["c_ball"] = c_ball

    input_vector = np.zeros(10)
    input_vector[0] = k_m_rot / r
    return mass, damping, stiffness, input_vector, p


def _matrix_frequency_response(frequencies: np.ndarray, mass: np.ndarray, damping: np.ndarray,
                               stiffness: np.ndarray, input_vector: np.ndarray,
                               output_index: int) -> np.ndarray:
    response = np.empty(frequencies.size, dtype=complex)
    for i, frequency in enumerate(frequencies):
        omega = 2.0 * np.pi * frequency
        dynamic_stiffness = stiffness - omega**2 * mass + 1j * omega * damping
        response[i] = np.linalg.solve(dynamic_stiffness, input_vector)[output_index]
    return response


def _linear_modes(mass: np.ndarray, stiffness: np.ndarray) -> np.ndarray:
    omega_squared = np.linalg.eigvals(np.linalg.solve(mass, stiffness))
    positive = np.real(omega_squared[np.real(omega_squared) > 1e-5])
    return np.sort(np.sqrt(positive) / (2.0 * np.pi))


def _rk4_linear(mass: np.ndarray, damping: np.ndarray, stiffness: np.ndarray,
                input_vector: np.ndarray, constants: dict[str, float], dt: float = 2.5e-6,
                duration: float = 0.085) -> tuple[np.ndarray, np.ndarray]:
    """Integrate an arbitrary second-order linear model with a true ZOH input."""
    count = mass.shape[0]
    inverse_mass = np.diag(1.0 / np.diag(mass))
    times = np.arange(0.0, duration + 0.5 * dt, dt)
    states = np.zeros((times.size, 2 * count), dtype=float)

    def rhs(state: np.ndarray, held_command: float) -> np.ndarray:
        position = state[:count]
        velocity = state[count:]
        acceleration = inverse_mass @ (input_vector * held_command - damping @ velocity - stiffness @ position)
        return np.concatenate((velocity, acceleration))

    for i in range(times.size - 1):
        held_command = command_position(times[i] + 0.5 * dt, constants["quarter_step"])
        y = states[i]
        k1 = rhs(y, held_command)
        k2 = rhs(y + 0.5 * dt * k1, held_command)
        k3 = rhs(y + 0.5 * dt * k2, held_command)
        k4 = rhs(y + dt * k3, held_command)
        states[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return times, states


def full_reduced_verification(frequencies: np.ndarray, constants: dict[str, float]) -> dict[str, object]:
    full_m, full_c, full_k, full_b, p = full_linear_matrices()
    reduced_m, reduced_c, reduced_k, reduced_b = linear_matrices((), "none")
    full_response = _matrix_frequency_response(frequencies, full_m, full_c, full_k, full_b, 9)
    reduced_response = _matrix_frequency_response(frequencies, reduced_m, reduced_c, reduced_k, reduced_b, 1)
    times, full_states = _rk4_linear(full_m, full_c, full_k, full_b, constants)
    reduced_times, reduced_states = _rk4_linear(reduced_m, reduced_c, reduced_k, reduced_b, constants)
    if not np.array_equal(times, reduced_times):
        raise RuntimeError("Full and reduced verification time grids differ")
    command = np.array([command_position(t, constants["quarter_step"]) for t in times])
    full_stage = full_states[:, 9]
    reduced_stage = reduced_states[:, 1]
    residual = full_stage - reduced_stage
    return {
        "parameters": p,
        "full_modes": _linear_modes(full_m, full_k),
        "reduced_modes": _linear_modes(reduced_m, reduced_k),
        "full_response": full_response,
        "reduced_response": reduced_response,
        "times": times,
        "command": command,
        "full_stage": full_stage,
        "reduced_stage": reduced_stage,
        "residual": residual,
        "rms_residual_nm": float(np.sqrt(np.mean(residual**2)) * 1e9),
        "peak_residual_nm": float(np.max(np.abs(residual)) * 1e9),
    }


def linear_matrices(sites: tuple[str, ...], friction_model: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return M, C, K, B for the presliding linearization of one case."""
    m_d, m_s = MODEL["m_d"], MODEL["m_s"]
    k_ax, k_m, c_ax = MODEL["k_ax"], MODEL["K_m"], MODEL["c_ax"]
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    mass = np.diag([m_d, m_s])
    constants = physical_constants()
    damping = c_ax * coupling + constants["c_m"] * np.outer(H["d"], H["d"])
    stiffness = np.array([[k_m + k_ax, -k_ax], [-k_ax, k_ax]], dtype=float)
    for site in sites:
        outer = np.outer(H[site], H[site])
        p = FRICTION[site]
        stiffness += p["sigma0"] * outer
        site_damping = p["sigma2"] if friction_model == "gms" else p["sigma1"] + p["sigma2"]
        damping += site_damping * outer
    input_vector = np.array([k_m, 0.0])
    return mass, damping, stiffness, input_vector


def frequency_responses() -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, dict[str, float | np.ndarray]]]:
    frequencies = np.logspace(np.log10(5.0), np.log10(3000.0), 3200)
    s = 1j * 2.0 * np.pi * frequencies
    responses: dict[str, np.ndarray] = {}
    metrics: dict[str, dict[str, float | np.ndarray]] = {}
    for key, case in CASES.items():
        mass, damping, stiffness, input_vector = linear_matrices(case["sites"], case["friction"])
        z11 = mass[0, 0] * s**2 + damping[0, 0] * s + stiffness[0, 0]
        z12 = damping[0, 1] * s + stiffness[0, 1]
        z21 = damping[1, 0] * s + stiffness[1, 0]
        z22 = mass[1, 1] * s**2 + damping[1, 1] * s + stiffness[1, 1]
        determinant = z11 * z22 - z12 * z21
        response = (-z21 * input_vector[0] + z11 * input_vector[1]) / determinant
        responses[key] = response
        omega2 = np.linalg.eigvals(np.linalg.solve(mass, stiffness))
        modes = np.sort(np.sqrt(np.maximum(np.real(omega2), 0.0)) / (2.0 * np.pi))
        dc_gain = float((np.linalg.solve(stiffness, input_vector))[1])
        metrics[key] = {"modes": modes, "dc_gain": dc_gain}
    return frequencies, responses, metrics


def command_position(t: float, quarter_step: float) -> float:
    """Closed four-increment sequence; every move is at most one quarter step."""
    if t < 0.005:
        return 0.0
    if t < 0.025:
        return quarter_step
    if t < 0.045:
        return 0.0
    if t < 0.065:
        return -quarter_step
    return 0.0


def presliding_command_position(t: float, microstep: float) -> float:
    """Nested back-and-forth reversals quantized to 1/32 of a full step."""
    if t < PRESLIDING_START:
        return 0.0
    index = min(int((t - PRESLIDING_START) // PRESLIDING_HOLD),
                PRESLIDING_LEVELS.size - 1)
    return float(PRESLIDING_LEVELS[index] * microstep)


def stribeck(velocity: float, p: dict[str, float]) -> float:
    ratio = abs(velocity) / p["v_s"]
    return p["F_c"] + (p["F_s"] - p["F_c"]) * np.exp(-(ratio ** p["delta"]))


def lugre_site(velocity: float, z: float, p: dict[str, float]) -> tuple[float, float]:
    attraction = max(stribeck(velocity, p), 1e-12)
    z_dot = velocity - p["sigma0"] * abs(velocity) * z / attraction
    force = p["sigma0"] * z + p["sigma1"] * z_dot + p["sigma2"] * velocity
    return z_dot, force


def gms_site(velocity: float, element_forces: np.ndarray,
             p: dict[str, float]) -> tuple[np.ndarray, float]:
    """Return GMS derivatives after branch tests on the current RK trial state.

    Ordering is intentional: zero velocity holds the state; otherwise the
    reversal/re-stick test is evaluated first, then the yield test, and only
    then is either the stuck or slip derivative assigned.  No derivative is
    used to choose its own branch within the same RHS evaluation.
    """
    threshold = np.maximum(GMS_WEIGHTS * stribeck(velocity, p), 1e-12)
    stiffness = GMS_STIFFNESS_FRACTIONS * p["sigma0"]
    derivatives = np.zeros(GMS_N)
    if abs(velocity) > 1e-14:
        direction = np.sign(velocity)
        for i in range(GMS_N):
            re_stick = velocity * element_forces[i] <= 0.0
            below_yield = abs(element_forces[i]) < threshold[i]
            if re_stick:
                derivatives[i] = stiffness[i] * velocity
            elif below_yield:
                derivatives[i] = stiffness[i] * velocity
            else:
                # Stable slip attraction to F_i = sign(v) nu_i s(v).
                derivatives[i] = p["C_gms"] * (direction - element_forces[i] / threshold[i])
    total_force = float(np.sum(element_forces) + p["sigma2"] * velocity)
    return derivatives, total_force


def nonlinear_rhs(t: float, state: np.ndarray, case: dict[str, object], constants: dict[str, float],
                  held_command: float | None = None) -> np.ndarray:
    x_d, x_s, v_d, v_s = state[:4]
    command = command_position(t, constants["quarter_step"])
    if held_command is not None:
        command = held_command
    lag = constants["kappa"] * (command - x_d)
    magnetic_force = constants["F_max"] * np.sin(lag)
    electromagnetic_damping = constants["c_m"] * v_d
    axial_force = MODEL["k_ax"] * (x_d - x_s) + MODEL["c_ax"] * (v_d - v_s)

    velocities = {"g": v_s, "n": v_d - v_s, "d": v_d}
    forces = {"g": 0.0, "n": 0.0, "d": 0.0}
    derivative = np.zeros_like(state)
    if case["friction"] == "lugre":
        for site in case["sites"]:
            index = LUGRE_INDEX[site]
            derivative[index], forces[site] = lugre_site(
                velocities[site], state[index], FRICTION[site])
    elif case["friction"] == "gms":
        for site in case["sites"]:
            start = GMS_START[site]
            stop = start + GMS_N
            derivative[start:stop], forces[site] = gms_site(
                velocities[site], state[start:stop], FRICTION[site])

    a_d = (magnetic_force - electromagnetic_damping - axial_force
           - forces["n"] - forces["d"]) / MODEL["m_d"]
    a_s = (axial_force + forces["n"] - forces["g"]) / MODEL["m_s"]
    derivative[:4] = (v_d, v_s, a_d, a_s)
    return derivative


def rk4_case_with_command(case: dict[str, object], constants: dict[str, float],
                          command_function, duration: float,
                          dt: float = 5.0e-6) -> tuple[np.ndarray, np.ndarray]:
    """Integrate a case with an arbitrary zero-order-held position command."""
    times = np.arange(0.0, duration + 0.5 * dt, dt)
    states = np.zeros((times.size, STATE_SIZE), dtype=float)
    for i in range(times.size - 1):
        t = times[i]
        y = states[i]
        # Treat the discrete command as a true zero-order hold.  One midpoint
        # sample is fixed across all RK4 stages, so an endpoint discontinuity
        # cannot leak backward into the preceding integration interval.
        held_command = float(command_function(t + 0.5 * dt))
        k1 = nonlinear_rhs(t, y, case, constants, held_command)
        k2 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k1, case, constants, held_command)
        k3 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k2, case, constants, held_command)
        k4 = nonlinear_rhs(t + dt, y + dt * k3, case, constants, held_command)
        states[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return times, states


def rk4_case(case: dict[str, object], constants: dict[str, float], dt: float = 5.0e-6,
             duration: float = 0.085) -> tuple[np.ndarray, np.ndarray]:
    return rk4_case_with_command(
        case, constants,
        lambda t: command_position(t, constants["quarter_step"]),
        duration=duration, dt=dt,
    )


def friction_force_history(case: dict[str, object], states: np.ndarray,
                           site: str) -> np.ndarray:
    """Recover an executed site's constitutive force from integrated states."""
    if site not in case["sites"]:
        return np.zeros(states.shape[0])
    velocity = states[:, 3] if site == "g" else states[:, 2] - states[:, 3]
    p = FRICTION[site]
    if case["friction"] == "lugre":
        state = states[:, LUGRE_INDEX[site]]
        return np.array([
            lugre_site(float(v), float(z), p)[1] for v, z in zip(velocity, state)
        ])
    start = GMS_START[site]
    stop = start + GMS_N
    return np.sum(states[:, start:stop], axis=1) + p["sigma2"] * velocity


def presliding_responses(constants: dict[str, float]) -> dict[str, object]:
    """Run the matched guideway cases through a nested reversal sequence."""
    microstep = constants["quarter_step"] / 8.0
    duration = PRESLIDING_START + PRESLIDING_HOLD * PRESLIDING_LEVELS.size
    command_function = lambda t: presliding_command_position(t, microstep)
    results: dict[str, np.ndarray] = {}
    forces: dict[str, np.ndarray] = {}
    metrics: dict[str, dict[str, object]] = {}
    times: np.ndarray | None = None

    for key in ("A", "A2"):
        times, states = rk4_case_with_command(
            CASES[key], constants, command_function, duration=duration)
        results[key] = states
        forces[key] = friction_force_history(CASES[key], states, "g")

    assert times is not None
    command = np.array([command_function(t) for t in times])
    active = times >= PRESLIDING_START
    for key in ("A", "A2"):
        error = command - results[key][:, 1]
        endpoint_error = []
        endpoint_force = []
        endpoint_stage = []
        for level_index in range(PRESLIDING_LEVELS.size):
            plateau_end = PRESLIDING_START + (level_index + 1) * PRESLIDING_HOLD
            window = (times >= plateau_end - 0.002) & (times < plateau_end - 0.5e-9)
            if level_index == PRESLIDING_LEVELS.size - 1:
                window = times >= plateau_end - 0.002
            endpoint_error.append(float(np.mean(error[window])))
            endpoint_force.append(float(np.mean(forces[key][window])))
            endpoint_stage.append(float(np.mean(results[key][window, 1])))
        endpoint_error_array = np.asarray(endpoint_error)
        endpoint_force_array = np.asarray(endpoint_force)
        error_mismatch = np.array([
            abs(endpoint_error_array[first] - endpoint_error_array[second])
            for first, second in PRESLIDING_RETURN_PAIRS
        ])
        force_mismatch = np.array([
            abs(endpoint_force_array[first] - endpoint_force_array[second])
            for first, second in PRESLIDING_RETURN_PAIRS
        ])
        metrics[key] = {
            "whole_rms_nm": float(np.sqrt(np.mean(error[active] ** 2)) * 1e9),
            "max_abs_deviation_nm": float(np.max(np.abs(error[active])) * 1e9),
            "final_mean_nm": float(endpoint_error_array[-1] * 1e9),
            "return_error_mismatch_nm": float(np.mean(error_mismatch) * 1e9),
            "return_force_mismatch_N": float(np.mean(force_mismatch)),
            "max_force_N": float(np.max(np.abs(forces[key]))),
            "endpoint_error_nm": endpoint_error_array * 1e9,
            "endpoint_force_N": endpoint_force_array,
            "endpoint_stage_um": np.asarray(endpoint_stage) * 1e6,
            "pair_error_mismatch_nm": error_mismatch * 1e9,
            "pair_force_mismatch_N": force_mismatch,
        }
    return {
        "times": times,
        "command": command,
        "results": results,
        "forces": forces,
        "metrics": metrics,
        "microstep": microstep,
        "duration": duration,
    }


def final_window_rms_error_nm(times: np.ndarray, states: np.ndarray,
                              constants: dict[str, float]) -> float:
    """Return RMS(command-stage) over the final 2 ms on the given time grid."""
    command = np.array([command_position(t, constants["quarter_step"]) for t in times])
    final_window = times >= (times[-1] - 0.002)
    error = command - states[:, 1]
    return float(np.sqrt(np.mean(error[final_window] ** 2)) * 1e9)


def time_responses(constants: dict[str, float]) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, dict[str, float]]]:
    results: dict[str, np.ndarray] = {}
    metrics: dict[str, dict[str, float]] = {}
    times: np.ndarray | None = None
    for key, case in CASES.items():
        times, states = rk4_case(case, constants)
        results[key] = states
    assert times is not None
    command = np.array([command_position(t, constants["quarter_step"]) for t in times])
    final_window = times >= (times[-1] - 0.002)
    first_plateau = (times >= 0.005) & (times < 0.025)
    for key, states in results.items():
        error = command - states[:, 1]
        first_peak = float(np.max(states[first_plateau, 1]))
        metrics[key] = {
            "mean_final_error_nm": float(np.mean(error[final_window]) * 1e9),
            "rms_final_error_nm": final_window_rms_error_nm(times, states, constants),
            "rms_sequence_deviation_nm": float(np.sqrt(np.mean(error ** 2)) * 1e9),
            "max_abs_deviation_nm": float(np.max(np.abs(error)) * 1e9),
            "max_stage_um": float(np.max(np.abs(states[:, 1])) * 1e6),
            "first_peak_um": first_peak * 1e6,
            "first_overshoot_pct": max(0.0, (first_peak / constants["quarter_step"] - 1.0) * 100.0),
        }
    return times, command, results, metrics


def gms_step_halving_convergence(constants: dict[str, float], base_times: np.ndarray,
                                 base_results: dict[str, np.ndarray]) -> dict[str, dict[str, object]]:
    """Compare final-window RMS under h, h/2, and h/4 for all GMS cases."""
    study: dict[str, dict[str, object]] = {}
    for key in ("A2", "B2", "C2"):
        rms_values: list[float] = []
        for dt in GMS_CONVERGENCE_DTS:
            if np.isclose(dt, 5.0e-6, rtol=0.0, atol=1.0e-15):
                times, states = base_times, base_results[key]
            else:
                times, states = rk4_case(CASES[key], constants, dt=dt)
            rms_values.append(final_window_rms_error_nm(times, states, constants))
        coarse_difference = abs(rms_values[0] - rms_values[1])
        fine_difference = abs(rms_values[1] - rms_values[2])
        study[key] = {
            "dt_s": GMS_CONVERGENCE_DTS,
            "rms_nm": tuple(rms_values),
            "coarse_difference_nm": coarse_difference,
            "fine_difference_nm": fine_difference,
            "fine_relative_pct": 100.0 * fine_difference / max(abs(rms_values[2]), 1.0e-15),
            "difference_ratio": coarse_difference / max(fine_difference, 1.0e-15),
        }
    return study


def plot_case_responses(frequencies: np.ndarray, responses: dict[str, np.ndarray],
                        times: np.ndarray, command: np.ndarray,
                        results: dict[str, np.ndarray], constants: dict[str, float],
                        time_metrics: dict[str, dict[str, float]]) -> list[Path]:
    """Create one self-contained Bode/step/deviation figure beside each case."""
    outputs: list[Path] = []
    time_ms = times * 1e3
    for key, case in CASES.items():
        fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.4))
        ax_mag, ax_pos = axes[0]
        ax_phase, ax_err = axes[1]
        response = responses[key]
        magnitude = 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))
        phase = np.unwrap(np.angle(response)) * 180.0 / np.pi
        color, line_style = case["color"], case["ls"]

        ax_mag.semilogx(frequencies, magnitude, color=color, linestyle=line_style, linewidth=1.8)
        ax_mag.axhline(0.0, color="#888888", linewidth=0.8)
        ax_mag.set_ylabel("Magnitude (dB)")
        ax_mag.set_title("Command-to-stage Bode magnitude")
        ax_mag.set_ylim(-90.0, 35.0)
        ax_mag.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)

        ax_phase.semilogx(frequencies, phase, color=color, linestyle=line_style, linewidth=1.8)
        ax_phase.set_xlabel("Frequency (Hz)")
        ax_phase.set_ylabel("Phase (deg)")
        ax_phase.set_ylim(-380.0, 30.0)
        ax_phase.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)

        stage = results[key][:, 1]
        error = command - stage
        ax_pos.step(time_ms, command * 1e6, where="post", color="#111111", linewidth=1.9,
                    label="Command")
        ax_pos.plot(time_ms, stage * 1e6, color=color, linestyle=line_style, linewidth=1.6,
                    label="Actual stage")
        ax_pos.set_ylabel("Position (µm)")
        ax_pos.set_title("Bounded commanded / actual motion")
        ax_pos.set_ylabel("Position (um)")
        ax_pos.legend(loc="upper right", frameon=True)

        ax_err.plot(time_ms, error * 1e9, color=color, linestyle=line_style, linewidth=1.5)
        ax_err.axhline(0.0, color="#888888", linewidth=0.8)
        ax_err.set_xlabel("Time (ms)")
        ax_err.set_ylabel(r"Modeled deviation $x_{cmd}-x_s$ (nm)")
        ax_err.set_title("Open-loop command-stage deviation")

        for axis in axes.flat:
            axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
            axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
        metric = time_metrics[key]
        fig.suptitle(case["label"], fontsize=14, fontweight="bold")
        fig.text(
            0.5, 0.042,
            f"Open-loop modeled command-stage deviation (not a servo tracking specification): "
            f"RMS={metric['rms_sequence_deviation_nm']:.1f} nm; "
            f"peak |deviation|={metric['max_abs_deviation_nm']:.1f} nm; "
            f"final-2-ms RMS={metric['rms_final_error_nm']:.1f} nm.",
            ha="center", fontsize=8.2, color="#555555",
        )
        fig.text(0.5, 0.012,
                 f"Nonlinear magnetic force; zeta_m={MODEL['zeta_m']:.2f}; each command increment <= {constants['quarter_step'] * 1e6:.2f} µm.",
                 ha="center", fontsize=8.5, color="#555555")
        fig.tight_layout(rect=(0.02, 0.075, 0.99, 0.95))
        output = ASSET_DIR / f"response_case_{key}.svg"
        fig.savefig(output, format="svg", bbox_inches="tight")
        plt.close(fig)
        outputs.append(output)
    return outputs


def plot_pairwise_comparison(frequencies: np.ndarray, responses: dict[str, np.ndarray],
                             times: np.ndarray, command: np.ndarray,
                             results: dict[str, np.ndarray],
                             time_metrics: dict[str, dict[str, float]]) -> Path:
    """Compare each LuGre case only with its topology-matched GMS case."""
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 10.5))
    time_ms = times * 1e3
    for row, (lugre_key, gms_key) in enumerate(PAIRS):
        ax_bode, ax_error = axes[row]
        for key in (lugre_key, gms_key):
            case = CASES[key]
            metric = time_metrics[key]
            metric_label = (f"{case['label']} | RMS {metric['rms_sequence_deviation_nm']:.1f} nm; "
                            f"max {metric['max_abs_deviation_nm']:.1f} nm")
            magnitude = 20.0 * np.log10(np.maximum(np.abs(responses[key]), 1e-15))
            ax_bode.semilogx(frequencies, magnitude, color=case["color"],
                             linestyle=case["ls"], linewidth=1.7, label=metric_label)
            error = command - results[key][:, 1]
            ax_error.plot(time_ms, error * 1e9, color=case["color"],
                          linestyle=case["ls"], linewidth=1.45, label=metric_label)
        ax_bode.axhline(0.0, color="#888888", linewidth=0.7)
        ax_bode.set_ylabel("Magnitude (dB)")
        ax_bode.set_ylim(-90.0, 30.0)
        ax_bode.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
        ax_bode.set_title(f"{lugre_key}/{gms_key}: Bode magnitude")
        ax_bode.legend(loc="lower left", fontsize=8)
        ax_error.axhline(0.0, color="#888888", linewidth=0.7)
        ax_error.set_ylabel("Modeled command-stage deviation (nm)")
        ax_error.set_title(f"{lugre_key}/{gms_key}: nonlinear sequence")
        ax_error.legend(loc="upper right", fontsize=8)
        for axis in (ax_bode, ax_error):
            axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
            axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
    axes[-1, 0].set_xlabel("Frequency (Hz)")
    axes[-1, 1].set_xlabel("Time (ms)")
    fig.suptitle("Topology-matched friction-model response comparison", fontsize=15, fontweight="bold")
    fig.text(0.5, 0.008,
             "RMS and max labels summarize open-loop command-stage deviation, not closed-loop tracking performance.",
             ha="center", fontsize=8.3, color="#555555")
    fig.tight_layout(rect=(0.02, 0.03, 0.99, 0.96))
    output = ASSET_DIR / "lugre_gms_pairwise_comparison.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_presliding_memory(experiment: dict[str, object]) -> Path:
    """Visualize nested-reversal command following and friction return-point memory."""
    times = experiment["times"]
    command = experiment["command"]
    results = experiment["results"]
    forces = experiment["forces"]
    metrics = experiment["metrics"]
    time_ms = times * 1e3

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.8))
    ax_motion, ax_error = axes[0]
    ax_memory, ax_metrics = axes[1]

    ax_motion.step(time_ms, command * 1e6, where="post", color="#111111",
                   linewidth=2.0, label="Command")
    for key in ("A", "A2"):
        case = CASES[key]
        stage = results[key][:, 1]
        ax_motion.plot(time_ms, stage * 1e6, color=case["color"],
                       linestyle=case["ls"], linewidth=1.5, label=case["label"])
        ax_error.plot(time_ms, (command - stage) * 1e9, color=case["color"],
                      linestyle=case["ls"], linewidth=1.35, label=case["label"])
        ax_memory.plot(stage * 1e6, forces[key], color=case["color"],
                       linestyle=case["ls"], linewidth=1.0, alpha=0.65,
                       label=case["label"])
        ax_memory.plot(metrics[key]["endpoint_stage_um"],
                       metrics[key]["endpoint_force_N"], color=case["color"],
                       linestyle="none", marker="o" if key == "A" else "s",
                       markersize=4.0, markerfacecolor="white")

    ax_motion.set_title("Nested microstep command and actual stage motion")
    ax_motion.set_ylabel("Position (um)")
    ax_motion.legend(loc="upper right", fontsize=8)
    ax_error.set_title("Modeled command-stage deviation over the same reversal history")
    ax_error.set_ylabel(r"Modeled deviation $x_{cmd}-x_s$ (nm)")
    ax_error.axhline(0.0, color="#777777", linewidth=0.8)
    ax_error.legend(loc="upper right", fontsize=8)
    ax_memory.set_title("Guideway friction memory loops")
    ax_memory.set_xlabel("Actual stage position (um)")
    ax_memory.set_ylabel("Guideway friction force (N)")
    ax_memory.axhline(0.0, color="#888888", linewidth=0.7)
    ax_memory.axvline(0.0, color="#888888", linewidth=0.7)
    ax_memory.legend(loc="best", fontsize=8)

    categories = ("Whole-sequence\nRMS", "Peak absolute\ndeviation",
                  "Return-point\nmismatch", "Final-origin\nabsolute deviation")
    x_positions = np.arange(len(categories), dtype=float)
    width = 0.34
    lugre_values = np.array([
        metrics["A"]["whole_rms_nm"],
        metrics["A"]["max_abs_deviation_nm"],
        metrics["A"]["return_error_mismatch_nm"],
        abs(metrics["A"]["final_mean_nm"]),
    ])
    gms_values = np.array([
        metrics["A2"]["whole_rms_nm"],
        metrics["A2"]["max_abs_deviation_nm"],
        metrics["A2"]["return_error_mismatch_nm"],
        abs(metrics["A2"]["final_mean_nm"]),
    ])
    bars_a = ax_metrics.bar(x_positions - width / 2.0, lugre_values, width,
                            color=CASES["A"]["color"], label="LuGre A")
    bars_a2 = ax_metrics.bar(x_positions + width / 2.0, gms_values, width,
                             color=CASES["A2"]["color"], label="GMS A2")
    ax_metrics.set_yscale("log")
    ax_metrics.set_xticks(x_positions, categories)
    ax_metrics.set_ylabel("Command-stage deviation metric (nm, log scale)")
    ax_metrics.set_title("Open-loop response: global and memory-sensitive metrics")
    ax_metrics.legend(loc="upper right", fontsize=8)
    for bars in (bars_a, bars_a2):
        for bar in bars:
            value = bar.get_height()
            ax_metrics.text(bar.get_x() + bar.get_width() / 2.0, value * 1.08,
                            f"{value:.2f}", ha="center", va="bottom", fontsize=7.4,
                            rotation=90)

    for axis in (ax_motion, ax_error):
        axis.set_xlabel("Time (ms)")
    for axis in axes.flat:
        axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)

    max_force = max(metrics[key]["max_force_N"] for key in ("A", "A2"))
    macro_fraction = 100.0 * max_force / FRICTION["g"]["F_s"]
    fig.suptitle("Presliding nested-reversal experiment: LuGre versus GMS",
                 fontsize=15, fontweight="bold")
    fig.text(
        0.5, 0.012,
        f"1 microstep = {experiment['microstep'] * 1e9:.2f} nm = 1/32 full step; "
        f"peak friction = {max_force:.3f} N ({macro_fraction:.1f}% of macro breakaway). "
        "Markers are 2 ms plateau-end means.",
        ha="center", fontsize=8.4, color="#555555",
    )
    fig.tight_layout(rect=(0.02, 0.045, 0.99, 0.95), h_pad=2.0, w_pad=1.5)
    output = ASSET_DIR / "presliding_memory_comparison.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_kinematic_diagram() -> Path:
    """Render both the ten-DOF topology and its two-DOF reduction."""
    fig, (full_ax, reduced_ax) = plt.subplots(2, 1, figsize=(13.2, 8.3),
                                              gridspec_kw={"height_ratios": [1.55, 1.0]})
    for ax in (full_ax, reduced_ax):
        ax.axis("off")

    full_ax.set_xlim(0.0, 13.0)
    full_ax.set_ylim(0.0, 5.4)
    full_ax.text(6.5, 5.12, "Revision 3 full topology — ten independent generalized coordinates",
                 ha="center", fontsize=15, fontweight="bold")

    def node(ax: plt.Axes, x: float, y: float, label: str, index: int, color: str) -> None:
        box = FancyBboxPatch((x - 0.46, y - 0.34), 0.92, 0.68, boxstyle="round,pad=0.04",
                             facecolor=color, edgecolor="#39434d", linewidth=1.25)
        ax.add_patch(box)
        ax.text(x, y + 0.04, label, ha="center", va="center", fontsize=10)
        ax.text(x, y - 0.23, f"q{index}", ha="center", va="center", fontsize=7.5, color="#59636d")

    torsion_x = np.linspace(1.1, 8.7, 5)
    for i, (x, label) in enumerate(zip(torsion_x, FULL_DOF_LABELS[:5]), start=1):
        node(full_ax, x, 3.95, label, i, "#dceef6")
    for x1, x2, spring_label in zip(torsion_x[:-1], torsion_x[1:],
                                    (r"$k_{c1}$", r"$k_{c2}$", r"$k_{\theta a}$", r"$k_{\theta b}$")):
        full_ax.plot([x1 + 0.47, x2 - 0.47], [3.95, 3.95], color="#277da1", linewidth=2)
        full_ax.text((x1 + x2) / 2.0, 4.18, spring_label, ha="center", fontsize=9)
    full_ax.add_patch(FancyArrowPatch((0.05, 3.95), (0.62, 3.95), arrowstyle="-|>",
                                      mutation_scale=14, color="#c08a00", linewidth=2))
    full_ax.text(0.34, 4.25, "$T_{mag}$", ha="center", fontsize=9, color="#8a6200")

    axial_x = np.linspace(1.1, 10.6, 5)
    for index, x, label in zip(range(6, 11), axial_x, FULL_DOF_LABELS[5:]):
        node(full_ax, x, 1.45, label, index, "#dff2ea")
    for x1, x2, spring_label in zip(axial_x[:-1], axial_x[1:],
                                    (r"$k_{sha}$", r"$k_{shb}$", "nut contact", r"$k_{mnt}$")):
        full_ax.plot([x1 + 0.47, x2 - 0.47], [1.45, 1.45], color="#218c74", linewidth=2)
        full_ax.text((x1 + x2) / 2.0, 1.71, spring_label, ha="center", fontsize=9)
    full_ax.plot([0.28, 0.28], [0.84, 2.06], color="#59636d", linewidth=3)
    full_ax.plot([0.28, 0.64], [1.45, 1.45], color="#59636d", linewidth=2)
    full_ax.text(0.33, 0.61, "$k_{brg}$ to ground", ha="left", fontsize=8.5)

    # Highlight the only cross-branch coupling: delta_n = u_n-u_e-r theta_s2.
    theta_nut_x = torsion_x[3]
    ue_x, un_x = axial_x[1], axial_x[3]
    full_ax.plot([theta_nut_x, theta_nut_x], [3.59, 2.72], color="#d97800", linewidth=2)
    full_ax.plot([theta_nut_x, ue_x, un_x], [2.72, 2.25, 2.25], color="#d97800", linewidth=2)
    full_ax.plot([ue_x, ue_x], [2.25, 1.80], color="#d97800", linewidth=2)
    full_ax.plot([un_x, un_x], [2.25, 1.80], color="#d97800", linewidth=2)
    full_ax.text(7.45, 2.52, r"ball-contact port: $\delta_n=u_n-u_e-r\theta_{s2}$",
                 ha="center", fontsize=9.5, color="#a45600", fontweight="bold")
    full_ax.text(11.35, 3.98, "Rotational branch", color="#1f5d73", fontsize=10, fontweight="bold")
    full_ax.text(11.35, 1.48, "Axial branch", color="#176a55", fontsize=10, fontweight="bold")
    full_ax.text(6.5, 0.15, "The command is an input, not a DOF. Friction states add internal states but not mechanical DOFs.",
                 ha="center", fontsize=9, color="#555555")

    reduced_ax.set_xlim(0.0, 13.0)
    reduced_ax.set_ylim(0.0, 3.2)
    reduced_ax.text(6.5, 2.94, "Band-limited reduction — two retained mechanical DOFs",
                    ha="center", fontsize=14, fontweight="bold")
    node(reduced_ax, 3.15, 1.55, r"$x_d,\;m_d$", 1, "#dceef6")
    node(reduced_ax, 9.55, 1.55, r"$x_s,\;m_s$", 2, "#dff2ea")
    reduced_ax.add_patch(FancyArrowPatch((0.55, 1.55), (2.65, 1.55), arrowstyle="-|>",
                                         mutation_scale=14, color="#c08a00", linewidth=2))
    reduced_ax.text(1.55, 1.82, "$F_{mag}(x_{cmd}-x_d)$", ha="center", fontsize=9)
    reduced_ax.plot([3.62, 9.08], [1.72, 1.72], color="#555555", linewidth=2)
    reduced_ax.plot([3.62, 9.08], [1.33, 1.33], color="#777777", linewidth=1.6)
    reduced_ax.text(6.35, 2.02, "$k_{ax}$: retained series compliance", ha="center", fontsize=9.5)
    reduced_ax.text(6.35, 1.04, "$c_{ax}$ and internal nut-friction port", ha="center", fontsize=9)
    reduced_ax.plot([3.15, 3.15], [1.18, 0.55], color="#6a4c93", linewidth=2)
    reduced_ax.plot([9.55, 9.55], [1.18, 0.55], color="#b23a48", linewidth=2)
    reduced_ax.text(2.15, 0.42, "$c_m$ and $F_{f,d}$", fontsize=9, color="#6a4c93")
    reduced_ax.text(9.75, 0.42, "$F_{f,g}$", fontsize=9, color="#9b2f3d")
    reduced_ax.text(6.5, 0.12, "Internal masses are collapsed only after their compliance, frequency separation, and identifiability are audited.",
                    ha="center", fontsize=8.8, color="#555555")

    fig.tight_layout(h_pad=1.2)
    output = ASSET_DIR / "kinematic_diagram.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_full_reduced_verification(frequencies: np.ndarray, verification: dict[str, object]) -> Path:
    """Plot linear Bode and bounded stepping for the full/reduced audit."""
    full_response = verification["full_response"]
    reduced_response = verification["reduced_response"]
    times = verification["times"]
    command = verification["command"]
    full_stage = verification["full_stage"]
    reduced_stage = verification["reduced_stage"]
    residual = verification["residual"]
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.0))
    ax_mag, ax_step, ax_phase, ax_residual = axes.flat
    styles = ((full_response, "Ten-DOF full", "#d97800", "-"),
              (reduced_response, "Two-DOF reduced", "#277da1", "--"))
    for response, label, color, line_style in styles:
        ax_mag.semilogx(frequencies, 20.0 * np.log10(np.maximum(np.abs(response), 1e-15)),
                        label=label, color=color, linestyle=line_style, linewidth=1.7)
        ax_phase.semilogx(frequencies, np.unwrap(np.angle(response)) * 180.0 / np.pi,
                          label=label, color=color, linestyle=line_style, linewidth=1.7)
    ax_mag.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
    ax_phase.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
    ax_mag.axhline(0.0, color="#888888", linewidth=0.7)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title("Command-to-stage Bode magnitude")
    ax_mag.legend(fontsize=8.5)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_title("Phase and discarded-mode onset")

    time_ms = times * 1e3
    ax_step.step(time_ms, command * 1e6, where="post", color="#111111", linewidth=1.8, label="Command")
    ax_step.plot(time_ms, full_stage * 1e6, color="#d97800", linewidth=1.5, label="Ten-DOF full")
    ax_step.plot(time_ms, reduced_stage * 1e6, color="#277da1", linestyle="--", linewidth=1.4,
                 label="Two-DOF reduced")
    ax_step.set_ylabel("Position (µm)")
    ax_step.set_title("Linear quarter-step back-and-forth verification")
    ax_step.legend(fontsize=8.2, loc="upper right")
    ax_residual.plot(time_ms, residual * 1e9, color="#9b2f3d", linewidth=1.35)
    ax_residual.axhline(0.0, color="#888888", linewidth=0.7)
    ax_residual.set_xlabel("Time (ms)")
    ax_residual.set_ylabel("Full − reduced (nm)")
    ax_residual.set_title(f"Reduction residual: RMS {verification['rms_residual_nm']:.2f} nm; peak {verification['peak_residual_nm']:.2f} nm")
    for ax in axes.flat:
        ax.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        ax.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
    fig.suptitle("Revision 3 full-versus-reduced executable verification", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.95))
    output = ASSET_DIR / "full_vs_reduced_verification.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_position_dependence() -> Path:
    positions = np.array([50.0, 150.0, 250.0])
    k_sha = np.array([2.0e8, 6.7e7, 4.0e7])
    k_ax = np.array([1.29e7, 1.14e7, 1.02e7])
    mode = np.array([736.0, 694.0, 657.0])
    fig, left = plt.subplots(figsize=(10.2, 4.6))
    right = left.twinx()
    left.plot(positions, k_ax / 1e6, marker="o", color="#277da1", linewidth=2,
              label="$k_{ax}$")
    right.plot(positions, mode, marker="s", color="#d97800", linewidth=2,
               label="predicted stage mode")
    for x, sha in zip(positions, k_sha):
        left.annotate(f"$k_{{sha}}$={sha / 1e6:.0f} MN/m", (x, np.interp(x, positions, k_ax / 1e6)),
                      xytext=(0, 10), textcoords="offset points", ha="center", fontsize=8)
    left.set_xlabel("Nut position from support bearing (mm)")
    left.set_ylabel("Reduced axial stiffness (MN/m)", color="#277da1")
    right.set_ylabel("Predicted stage mode (Hz)", color="#d97800")
    left.grid(True, color="#dedede", linewidth=0.7)
    left.set_title("Falsifiable position-dependence prediction")
    lines = left.lines + right.lines
    left.legend(lines, [line.get_label() for line in lines], loc="center right")
    fig.tight_layout()
    output = ASSET_DIR / "position_dependence.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_stepper_resonance_visibility() -> Path:
    """Show how damping and output selection hide or expose the low motor mode."""
    frequencies = np.logspace(np.log10(80.0), np.log10(400.0), 1400)
    mass, _baseline_damping, stiffness, input_vector = linear_matrices((), "none")
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    low_mode = _linear_modes(mass, stiffness)[0]
    damping_ratios = (0.02, 0.10, MODEL["zeta_m"])
    colors = ("#b23a48", "#d97800", "#277da1")
    fig, (stage_ax, drive_ax) = plt.subplots(1, 2, figsize=(11.4, 4.7), sharex=True)
    for zeta, color in zip(damping_ratios, colors):
        c_m = 2.0 * zeta * np.sqrt(MODEL["K_m"] * MODEL["m_d"])
        damping = MODEL["c_ax"] * coupling + c_m * np.outer(H["d"], H["d"])
        stage_response = _matrix_frequency_response(
            frequencies, mass, damping, stiffness, input_vector, 1)
        drive_response = _matrix_frequency_response(
            frequencies, mass, damping, stiffness, input_vector, 0)
        label = rf"$\zeta_m={zeta:.2f}$" + (" (executed baseline)" if zeta == MODEL["zeta_m"] else "")
        stage_ax.semilogx(frequencies, 20.0 * np.log10(np.maximum(np.abs(stage_response), 1e-15)),
                          color=color, linewidth=1.8, label=label)
        drive_ax.semilogx(frequencies, 20.0 * np.log10(np.maximum(np.abs(drive_response), 1e-15)),
                          color=color, linewidth=1.8, label=label)
    for axis in (stage_ax, drive_ax):
        axis.axvspan(155.0, 190.0, color="#8f6bb3", alpha=0.13,
                     label="155–190 Hz broad test feature")
        axis.axvline(low_mode, color="#252525", linestyle="--", linewidth=1.1,
                    label=f"model pole {low_mode:.1f} Hz")
        axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
        axis.set_xlabel("Frequency (Hz)")
        axis.set_ylabel("Magnitude (dB)")
        axis.set_xlim(80.0, 400.0)
    stage_ax.set_title(r"Measured output: $X_s/X_{cmd}$")
    drive_ax.set_title(r"Internal drive output: $X_d/X_{cmd}$")
    stage_ax.legend(fontsize=7.8, loc="best")
    drive_ax.legend(fontsize=7.8, loc="best")
    fig.suptitle("Low-frequency stepper-mode visibility: damping and output selection",
                 fontsize=14, fontweight="bold")
    fig.text(0.5, 0.012,
             "Detent torque is not included in this sensitivity plot because its amplitude and equilibrium phase are not identified.",
             ha="center", fontsize=8.5, color="#555555")
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.93))
    output = ASSET_DIR / "stepper_resonance_visibility.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_rotor_stage_transfer_functions(frequencies: np.ndarray) -> Path:
    """Plot command/drive/stage transfers and the mechanical drive-to-stage ratio."""
    mass, damping, stiffness, input_vector = linear_matrices((), "none")
    drive_response = _matrix_frequency_response(
        frequencies, mass, damping, stiffness, input_vector, 0)
    stage_response = _matrix_frequency_response(
        frequencies, mass, damping, stiffness, input_vector, 1)
    rotor_to_stage = stage_response / drive_response
    transfers = (
        (drive_response, r"$X_d/X_{cmd}$ — command to rotor-equivalent drive", "#b23a48", "-"),
        (stage_response, r"$X_s/X_{cmd}$ — command to stage", "#277da1", "-"),
        (rotor_to_stage, r"$X_s/X_d$ — rotor-equivalent drive to stage", "#218c74", "--"),
    )
    modes = _linear_modes(mass, stiffness)
    fig, (magnitude_ax, phase_ax) = plt.subplots(1, 2, figsize=(12.0, 4.9), sharex=True)
    for response, label, color, line_style in transfers:
        magnitude_ax.semilogx(
            frequencies, 20.0 * np.log10(np.maximum(np.abs(response), 1e-15)),
            color=color, linestyle=line_style, linewidth=1.8, label=label)
        phase_ax.semilogx(
            frequencies, np.unwrap(np.angle(response)) * 180.0 / np.pi,
            color=color, linestyle=line_style, linewidth=1.8, label=label)
    for axis in (magnitude_ax, phase_ax):
        for index, mode in enumerate(modes):
            axis.axvline(mode, color="#777777", linestyle=":" if index else "--", linewidth=1.0)
        axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
        axis.set_xlabel("Frequency (Hz)")
        axis.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
        axis.legend(fontsize=8.0, loc="best")
    magnitude_ax.set_ylabel("Magnitude (dB)")
    magnitude_ax.set_title("Transfer-function magnitude")
    magnitude_ax.set_ylim(-100.0, 35.0)
    phase_ax.set_ylabel("Phase (deg)")
    phase_ax.set_title("Transfer-function phase")
    phase_ax.text(modes[0] * 1.03, -345.0, f"{modes[0]:.1f} Hz low pole", fontsize=8, color="#555555")
    phase_ax.text(modes[1] * 1.03, -345.0, f"{modes[1]:.1f} Hz axial pole", fontsize=8, color="#555555")
    fig.suptitle("Rotor-equivalent drive and stage transfer functions",
                 fontsize=14.5, fontweight="bold")
    fig.text(0.5, 0.012,
             r"$X_s/X_d$ treats drive motion as prescribed and therefore does not contain the common motor/drive pole by itself.",
             ha="center", fontsize=8.5, color="#555555")
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.93))
    output = ASSET_DIR / "rotor_stage_transfer_functions.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def generated_summary(linear_metrics: dict[str, dict[str, float | np.ndarray]],
                      time_metrics: dict[str, dict[str, float]],
                      verification: dict[str, object]) -> str:
    lines = [
        "<!-- BEGIN GENERATED RESPONSE SUMMARY -->",
        "| Case | Friction law | Presliding modes (Hz) | DC gain $X_s/X_{cmd}$ | First-step overshoot | Full-sequence RMS deviation | Peak absolute deviation | Final-window RMS deviation |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, case in CASES.items():
        modes = linear_metrics[key]["modes"]
        mode_text = f"{modes[0]:.1f}, {modes[1]:.1f}"
        friction_label = {"none": "none", "lugre": "LuGre", "gms": "GMS"}[case["friction"]]
        lines.append(
            f"| {key} | {friction_label} | {mode_text} | {linear_metrics[key]['dc_gain']:.5f} | "
            f"{time_metrics[key]['first_overshoot_pct']:.1f}% | "
            f"{time_metrics[key]['rms_sequence_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['max_abs_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['rms_final_error_nm']:.1f} nm |"
        )
    lines.extend([
        "",
        "The three deviation columns use $d(t)=x_{cmd}(t)-x_s(t)$. They describe the open-loop modeled plant response under each friction law, not closed-loop servo tracking performance. The final column summarizes the last 2 ms of the nonlinear run and is not an identified settling specification. "
        "All cases include the separately highlighted electromagnetic damping assumption; Case 0 remains frictionless.",
        "",
        "### Generated reduction audit",
        "",
        "| Quantity | Executed value |",
        "|---|---:|",
        f"| Closure-derived $k_{{ball}}$ | {verification['parameters']['k_ball'] / 1e6:.3f} MN/m |",
        f"| Full-model reflected drivetrain mass | {verification['parameters']['m_d_reflected']:.3f} kg |",
        f"| Literal source-table reflected mass | {verification['parameters']['literal_table_m_d']:.3f} kg |",
        f"| Full/reduced sequence RMS residual | {verification['rms_residual_nm']:.3f} nm |",
        f"| Full/reduced sequence peak residual | {verification['peak_residual_nm']:.3f} nm |",
        "",
        "The literal table value is reported as an audit only; it is not silently used. The executable default uses the highlighted coupling-inertia assumption that closes the stated 59 kg reduction.",
        "<!-- END GENERATED RESPONSE SUMMARY -->",
    ])
    return "\n".join(lines)


def update_generated_summary(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED RESPONSE SUMMARY -->.*?<!-- END GENERATED RESPONSE SUMMARY -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated response summary markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_presliding_summary(experiment: dict[str, object]) -> str:
    metrics = experiment["metrics"]
    lines = [
        "<!-- BEGIN GENERATED PRESLIDING SUMMARY -->",
        "| Executed metric | LuGre A | GMS A2 | GMS change relative to LuGre |",
        "|---|---:|---:|---:|",
    ]
    rows = (
        ("Whole-sequence RMS command-stage deviation", "whole_rms_nm", "nm", False),
        ("Peak absolute command-stage deviation", "max_abs_deviation_nm", "nm", False),
        ("Mean repeated-return deviation mismatch", "return_error_mismatch_nm", "nm", False),
        ("Mean repeated-return friction-force mismatch", "return_force_mismatch_N", "N", False),
        ("Absolute mean error after final return to zero", "final_mean_nm", "nm", True),
    )
    for label, key, unit, use_absolute in rows:
        lugre = float(metrics["A"][key])
        gms = float(metrics["A2"][key])
        if use_absolute:
            lugre, gms = abs(lugre), abs(gms)
        reduction = 100.0 * (1.0 - gms / lugre) if lugre > 0.0 else 0.0
        precision = 4 if unit == "N" else 2
        lines.append(
            f"| {label} | {lugre:.{precision}f} {unit} | {gms:.{precision}f} {unit} | "
            f"{reduction:.1f}% lower |"
        )
    max_force = max(float(metrics[key]["max_force_N"]) for key in ("A", "A2"))
    lines.extend([
        "",
        f"The maximum executed guideway friction magnitude is **{max_force:.3f} N**, "
        f"or **{100.0 * max_force / FRICTION['g']['F_s']:.1f}%** of the provisional "
        f"{FRICTION['g']['F_s']:.1f} N macro breakaway level. The sequence therefore "
        "probes partial slip rather than gross sliding.",
        "",
        "The whole-sequence RMS includes the unavoidable error at every instantaneous command edge. "
        "The repeated-return and final-origin measures isolate the history dependence that this "
        "experiment is intended to distinguish.",
        "<!-- END GENERATED PRESLIDING SUMMARY -->",
    ])
    return "\n".join(lines)


def update_generated_presliding_summary(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED PRESLIDING SUMMARY -->.*?<!-- END GENERATED PRESLIDING SUMMARY -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated presliding summary markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_convergence_summary(study: dict[str, dict[str, object]]) -> str:
    lines = [
        "<!-- BEGIN GENERATED STEP HALVING SUMMARY -->",
        "| Case | 10.0 us | 5.0 us | 2.5 us | $\\Delta R_{10\\to5}$ | $\\Delta R_{5\\to2.5}$ | Difference ratio |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key in ("A2", "B2", "C2"):
        result = study[key]
        rms = result["rms_nm"]
        lines.append(
            f"| {key} | {rms[0]:.5f} nm | {rms[1]:.5f} nm | {rms[2]:.5f} nm | "
            f"{result['coarse_difference_nm']:.5f} nm | {result['fine_difference_nm']:.5f} nm | "
            f"{result['difference_ratio']:.2f} |"
        )
    monotone = all(
        float(study[key]["fine_difference_nm"]) < float(study[key]["coarse_difference_nm"])
        for key in ("A2", "B2", "C2")
    )
    max_relative = max(float(study[key]["fine_relative_pct"]) for key in ("A2", "B2", "C2"))
    if monotone:
        interpretation = (
            "The successive change decreases for all three GMS cases, which is consistent with "
            "time-step convergence for this reported metric."
        )
    else:
        interpretation = (
            "At least one case does not show a smaller second difference, so this metric is not yet "
            "demonstrably converged and a finer run is required."
        )
    lines.extend([
        "",
        interpretation + f" The largest 5.0-to-2.5 us relative change is **{max_relative:.4f}%**.",
        "",
        "These values use the identical 85 ms zero-order-held command and the identical final 2 ms "
        "RMS definition. Since GMS branch switching is evaluated at RK trial states without event "
        "localization, the difference ratio is a sensitivity indicator, not a claimed fourth-order "
        "convergence rate for the hybrid trajectory.",
        "<!-- END GENERATED STEP HALVING SUMMARY -->",
    ])
    return "\n".join(lines)


def update_generated_convergence_summary(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED STEP HALVING SUMMARY -->.*?<!-- END GENERATED STEP HALVING SUMMARY -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated step-halving summary markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def slugify(text: str) -> str:
    text = re.sub(r"\$.*?\$", "", text)
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"[^a-zA-Z0-9]+", "-", text).strip("-").lower()
    return text or "section"


def render_inline(text: str) -> str:
    stored: list[str] = []

    def keep(value: str) -> str:
        token = f"@@KEEP{len(stored)}@@"
        stored.append(value)
        return token

    def parameter_input(match: re.Match[str]) -> str:
        kind, key, value = match.group(1), match.group(2), match.group(3)
        css_class = "parameter-input assumed-input" if kind == "assumed" else "parameter-input"
        escaped_key = html.escape(key, quote=True)
        escaped_value = html.escape(value.strip(), quote=True)
        return keep(
            f'<input class="{css_class}" data-param="{escaped_key}" '
            f'data-default="{escaped_value}" value="{escaped_value}" '
            f'aria-label="Editable parameter {escaped_key}" spellcheck="false">'
        )

    text = re.sub(r"\[\[(input|assumed):([A-Za-z0-9_]+)=([^\]]+)\]\]", parameter_input, text)
    text = re.sub(r"`([^`]+)`", lambda m: keep(f"<code>{html.escape(m.group(1))}</code>"), text)
    text = re.sub(r"\$([^$\n]+)\$", lambda m: keep(f"\\({m.group(1)}\\)"), text)
    text = re.sub(
        r"!\[([^\]]*)\]\(([^)]+)\)",
        lambda m: keep(f'<img src="{html.escape(m.group(2), quote=True)}" alt="{html.escape(m.group(1), quote=True)}">'),
        text,
    )
    text = re.sub(
        r"\[([^\]]+)\]\(([^)]+)\)",
        lambda m: keep(f'<a href="{html.escape(m.group(2), quote=True)}">{html.escape(m.group(1))}</a>'),
        text,
    )
    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text)
    for index, value in enumerate(stored):
        text = text.replace(f"@@KEEP{index}@@", value)
    return text


def split_table_row(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def markdown_to_html(source: str) -> tuple[str, list[tuple[int, str, str]]]:
    lines = source.splitlines()
    output: list[str] = []
    toc: list[tuple[int, str, str]] = []
    used_ids: dict[str, int] = {}
    i = 0
    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()
        if not stripped:
            i += 1
            continue
        if stripped.startswith("```"):
            language = stripped[3:].strip()
            code_lines: list[str] = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i])
                i += 1
            i += 1
            output.append(f'<pre><code class="language-{html.escape(language)}">{html.escape(chr(10).join(code_lines))}</code></pre>')
            continue
        if stripped.startswith("$$"):
            # Accept all common display forms: a delimiter-only line, one line
            # with both delimiters, or multiline content sharing either line.
            first = stripped[2:]
            if first.endswith("$$"):
                equation = [first[:-2]]
                i += 1
            else:
                equation = [first] if first else []
                i += 1
                while i < len(lines):
                    candidate = lines[i]
                    if candidate.strip().endswith("$$"):
                        end_position = candidate.rfind("$$")
                        if candidate[:end_position]:
                            equation.append(candidate[:end_position])
                        i += 1
                        break
                    equation.append(candidate)
                    i += 1
            output.append('<div class="display-math">\\[' + "\n".join(equation) + '\\]</div>')
            continue
        heading = re.match(r"^(#{1,6})\s+(.+)$", raw)
        if heading:
            level = len(heading.group(1))
            title = heading.group(2).strip()
            base = slugify(title)
            count = used_ids.get(base, 0)
            used_ids[base] = count + 1
            section_id = base if count == 0 else f"{base}-{count + 1}"
            output.append(f'<h{level} id="{section_id}">{render_inline(title)}</h{level}>')
            if level in (2, 3):
                toc.append((level, re.sub(r"[*`$]", "", title), section_id))
            i += 1
            continue
        if stripped in ("---", "***", "___"):
            output.append("<hr>")
            i += 1
            continue
        if stripped.startswith("<!--"):
            comment = [raw]
            while "-->" not in comment[-1] and i + 1 < len(lines):
                i += 1
                comment.append(lines[i])
            output.append("\n".join(comment))
            i += 1
            continue
        if re.match(r"^</?(details|summary)(\s|>|$)", stripped):
            if stripped.startswith("<summary>") and stripped.endswith("</summary>"):
                middle = stripped[len("<summary>"):-len("</summary>")]
                output.append(f"<summary>{render_inline(middle)}</summary>")
            else:
                output.append(stripped)
            i += 1
            continue
        if stripped.startswith("|") and i + 1 < len(lines) and re.match(r"^\s*\|?\s*:?-{3,}", lines[i + 1]):
            header_cells = split_table_row(raw)
            align_cells = split_table_row(lines[i + 1])
            aligns = []
            for cell in align_cells:
                left, right = cell.startswith(":"), cell.endswith(":")
                aligns.append("center" if left and right else "right" if right else "left")
            table = ["<div class=\"table-wrap\"><table><thead><tr>"]
            for j, cell in enumerate(header_cells):
                align = aligns[j] if j < len(aligns) else "left"
                table.append(f'<th style="text-align:{align}">{render_inline(cell)}</th>')
            table.append("</tr></thead><tbody>")
            i += 2
            while i < len(lines) and lines[i].strip().startswith("|"):
                table.append("<tr>")
                for j, cell in enumerate(split_table_row(lines[i])):
                    align = aligns[j] if j < len(aligns) else "left"
                    table.append(f'<td style="text-align:{align}">{render_inline(cell)}</td>')
                table.append("</tr>")
                i += 1
            table.append("</tbody></table></div>")
            output.append("".join(table))
            continue
        if stripped.startswith(">"):
            quote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                quote_lines.append(lines[i].strip()[1:].strip())
                i += 1
            output.append(f"<blockquote>{render_inline(' '.join(quote_lines))}</blockquote>")
            continue
        list_match = re.match(r"^\s*([-+*]|\d+\.)\s+(.+)$", raw)
        if list_match:
            ordered = list_match.group(1)[0].isdigit()
            tag = "ol" if ordered else "ul"
            items: list[str] = []
            while i < len(lines):
                item = re.match(r"^\s*([-+*]|\d+\.)\s+(.+)$", lines[i])
                if not item or item.group(1)[0].isdigit() != ordered:
                    break
                items.append(f"<li>{render_inline(item.group(2).strip())}</li>")
                i += 1
            output.append(f"<{tag}>" + "".join(items) + f"</{tag}>")
            continue
        if stripped.startswith("<") and stripped.endswith(">"):
            output.append(raw)
            i += 1
            continue

        paragraph = [stripped]
        i += 1
        while i < len(lines):
            candidate = lines[i].strip()
            if not candidate:
                break
            if (candidate.startswith(("#", "```", "$$", ">", "|", "<", "---")) or
                    re.match(r"^\s*([-+*]|\d+\.)\s+", lines[i])):
                break
            paragraph.append(candidate)
            i += 1
        output.append(f"<p>{render_inline(' '.join(paragraph))}</p>")
    return "\n".join(output), toc


def html_page(markdown_path: Path) -> str:
    source = markdown_path.read_text(encoding="utf-8")
    body, toc = markdown_to_html(source)
    title_match = re.search(r"^#\s+(.+)$", source, flags=re.MULTILINE)
    title = title_match.group(1) if title_match else markdown_path.stem
    toc_html = []
    for level, label, section_id in toc:
        toc_html.append(f'<a class="toc-level-{level}" href="#{section_id}">{html.escape(label)}</a>')
    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<script>
window.MathJax = {{tex: {{inlineMath: [['\\\\(','\\\\)']], displayMath: [['\\\\[','\\\\]']]}}, svg: {{fontCache: 'global'}}}};
</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
<style>
:root {{ --bg:#eef2f5; --card:#fff; --text:#20262d; --muted:#65717d; --line:#d9e0e6; --accent:#1f6f8b; --soft:#f4f8fa; --code:#17212b; --assumed:#fff0b8; --assumed-line:#d49b00; }}
html[data-theme="dark"] {{ --bg:#11171d; --card:#182129; --text:#e9eef2; --muted:#aab5be; --line:#33414c; --accent:#73c2df; --soft:#202c35; --code:#0d1217; --assumed:#5b4810; --assumed-line:#f1bf36; }}
* {{ box-sizing:border-box; }} html {{ scroll-behavior:smooth; }}
body {{ margin:0; background:var(--bg); color:var(--text); font:16px/1.67 Inter,Segoe UI,system-ui,sans-serif; }}
.topbar {{ position:sticky; top:0; z-index:20; display:flex; gap:.75rem; align-items:center; padding:.65rem 1rem; background:color-mix(in srgb,var(--card) 94%,transparent); border-bottom:1px solid var(--line); backdrop-filter:blur(9px); }}
.topbar .name {{ font-weight:700; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; margin-right:auto; }}
button {{ color:var(--text); background:var(--soft); border:1px solid var(--line); border-radius:7px; padding:.42rem .7rem; cursor:pointer; }}
.layout {{ display:grid; grid-template-columns:280px minmax(0,1fr); gap:1.5rem; max-width:1510px; margin:0 auto; padding:1.5rem; }}
nav {{ position:sticky; top:4.4rem; align-self:start; max-height:calc(100vh - 5.5rem); overflow:auto; background:var(--card); border:1px solid var(--line); border-radius:12px; padding:1rem; }}
nav .caption {{ font-size:.78rem; color:var(--muted); text-transform:uppercase; letter-spacing:.08em; margin-bottom:.6rem; }}
nav a {{ display:block; color:var(--muted); text-decoration:none; border-left:2px solid transparent; padding:.24rem .35rem; font-size:.88rem; }} nav a:hover {{ color:var(--accent); border-color:var(--accent); }} nav .toc-level-3 {{ padding-left:1.15rem; font-size:.81rem; }}
article {{ width:100%; max-width:1100px; background:var(--card); border:1px solid var(--line); border-radius:12px; padding:clamp(1.25rem,4vw,3.7rem); box-shadow:0 12px 32px rgba(22,36,46,.07); }}
h1,h2,h3,h4 {{ line-height:1.24; scroll-margin-top:5rem; }} h1 {{ font-size:clamp(2rem,4vw,3rem); margin-top:0; }} h2 {{ margin-top:2.8rem; padding-bottom:.38rem; border-bottom:1px solid var(--line); }} h3 {{ margin-top:2rem; color:var(--accent); }}
p {{ margin:.85rem 0; }} a {{ color:var(--accent); }} strong {{ color:var(--text); }} hr {{ border:0; border-top:1px solid var(--line); margin:2rem 0; }}
blockquote {{ margin:1.2rem 0; padding:.75rem 1rem; background:var(--soft); border-left:4px solid var(--accent); border-radius:0 8px 8px 0; color:var(--muted); }}
.table-wrap {{ overflow-x:auto; margin:1.2rem 0; }} table {{ width:100%; border-collapse:collapse; font-size:.92rem; }} th,td {{ border:1px solid var(--line); padding:.55rem .65rem; vertical-align:top; }} th {{ background:var(--soft); }}
.parameter-input {{ width:100%; min-width:7rem; padding:.38rem .48rem; color:var(--text); background:var(--card); border:1px solid var(--line); border-radius:5px; font:inherit; font-variant-numeric:tabular-nums; }}
.parameter-input:focus {{ outline:2px solid var(--accent); outline-offset:1px; }}
.assumed-input {{ background:var(--assumed); border-color:var(--assumed-line); font-weight:700; }}
.edit-note {{ margin:0 0 1.2rem; padding:.7rem .9rem; border:1px solid var(--line); border-radius:8px; background:var(--soft); color:var(--muted); font-size:.86rem; }}
.save-status {{ display:inline-block; margin-left:.55rem; padding:.14rem .45rem; border-radius:999px; border:1px solid var(--line); background:var(--card); color:var(--muted); font-weight:650; }}
.save-status.ok {{ color:#176a55; border-color:#4fa88e; }} .save-status.warn {{ color:#9b5b00; border-color:#d49b00; }}
.assumed-swatch {{ display:inline-block; width:1.1rem; height:.8rem; margin:0 .25rem; vertical-align:middle; background:var(--assumed); border:1px solid var(--assumed-line); border-radius:3px; }}
.live-transfer-panel {{ margin:1.2rem 0 1.8rem; padding:1rem; border:2px solid var(--accent); border-radius:10px; background:var(--soft); }}
.live-transfer-panel .live-summary {{ margin:0 0 .7rem; color:var(--muted); font-size:.9rem; font-variant-numeric:tabular-nums; }}
.live-plot-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:.85rem; }}
.live-plot-card {{ min-width:0; padding:.6rem; border:1px solid var(--line); border-radius:8px; background:var(--card); }}
.live-plot-card h4 {{ margin:.1rem 0 .35rem; color:var(--accent); }}
.live-plot-card svg {{ display:block; width:100%; height:auto; color:var(--text); }}
details {{ margin:1rem 0; border:1px solid var(--line); border-radius:9px; background:color-mix(in srgb,var(--soft) 45%,var(--card)); padding:.2rem .9rem .8rem; }} details details {{ margin-left:.45rem; }} summary {{ cursor:pointer; font-weight:700; padding:.75rem .1rem; color:var(--accent); }}
pre {{ overflow:auto; background:var(--code); color:#e8edf2; border-radius:9px; padding:1rem; font-size:.87rem; }} code {{ font-family:Cascadia Code,Consolas,monospace; }} p code,li code,td code {{ background:var(--soft); border:1px solid var(--line); border-radius:4px; padding:.1rem .28rem; }}
.display-math {{ overflow-x:auto; padding:.5rem 0; }} img {{ display:block; max-width:100%; height:auto; margin:1.3rem auto; border-radius:6px; }}
.footer {{ color:var(--muted); font-size:.78rem; margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line); }}
@media (max-width:920px) {{ .layout {{ grid-template-columns:1fr; padding:.7rem; }} nav {{ position:relative; top:auto; max-height:18rem; }} article {{ padding:1.2rem; }} .hide-small {{ display:none; }} .live-plot-grid {{ grid-template-columns:1fr; }} }}
@media print {{ .topbar,nav {{ display:none; }} body {{ background:white; }} .layout {{ display:block; padding:0; }} article {{ max-width:none; border:0; box-shadow:none; }} details {{ break-inside:avoid; }} }}
</style>
</head>
<body>
<div class="topbar"><span class="name">{html.escape(title)}</span><button onclick="setDetails(true)">Expand derivations</button><button onclick="setDetails(false)">Collapse</button><button onclick="saveParameterInputs()">Save variables</button><button onclick="saveEditedHtml()">Save HTML copy</button><button onclick="resetParameterInputs()">Reset inputs</button><button class="hide-small" onclick="toggleTheme()">Theme</button><button class="hide-small" onclick="window.print()">Print</button></div>
<div class="layout"><nav><div class="caption">On this page</div>{''.join(toc_html)}</nav><article><div class="edit-note"><span class="assumed-swatch"></span>Amber inputs are unidentified assumptions. Values auto-save to browser storage and the page URL; “Save HTML copy” embeds them in a chosen HTML file. The live transfer panel recalculates immediately. Publication SVGs are regenerated only by the Python build.<span id="parameter-save-status" class="save-status">Loading values…</span></div>{body}<div class="footer">Rendered from {html.escape(markdown_path.name)} · {generated}</div></article></div>
<script>
function setDetails(open) {{ document.querySelectorAll('details').forEach(d => d.open=open); }}
function toggleTheme() {{ const root=document.documentElement; root.dataset.theme=root.dataset.theme==='dark'?'light':'dark'; }}
const parameterStorageKey = 'model-parameters:' + document.title + ':' + location.pathname;
let parameterSaveTimer = null;
function currentParameterValues() {{
  const values = {{}};
  document.querySelectorAll('.parameter-input').forEach(input => values[input.dataset.param] = input.value);
  return values;
}}
function setParameterStatus(message, kind='') {{
  const status = document.getElementById('parameter-save-status');
  if (!status) return;
  status.textContent = message;
  status.className = 'save-status' + (kind ? ' ' + kind : '');
}}
function readHashParameterValues() {{
  if (!location.hash.startsWith('#params=')) return {{}};
  try {{ return JSON.parse(decodeURIComponent(location.hash.slice(8))); }} catch (_) {{ return {{}}; }}
}}
function writeHashParameterValues(values) {{
  try {{
    const url = new URL(location.href);
    url.hash = 'params=' + encodeURIComponent(JSON.stringify(values));
    history.replaceState(null, '', url.href);
    return true;
  }} catch (_) {{ return false; }}
}}
function persistParameterInputs(showStatus=true) {{
  const values = currentParameterValues();
  let browserStorageSaved = true;
  try {{ localStorage.setItem(parameterStorageKey, JSON.stringify(values)); }}
  catch (_) {{ browserStorageSaved = false; }}
  const urlSaved = writeHashParameterValues(values);
  if (showStatus) {{
    const time = new Date().toLocaleTimeString();
    if (browserStorageSaved || urlSaved) setParameterStatus('Variables saved · ' + time, 'ok');
    else setParameterStatus('Browser storage unavailable — use Save HTML copy', 'warn');
  }}
  return browserStorageSaved || urlSaved;
}}
function saveParameterInputs() {{
  persistParameterInputs(true);
  refreshInteractivePlots();
}}
function scheduleParameterUpdate() {{
  setParameterStatus('Saving and updating live plots…', 'warn');
  if (parameterSaveTimer) clearTimeout(parameterSaveTimer);
  parameterSaveTimer = setTimeout(() => {{
    persistParameterInputs(true);
    refreshInteractivePlots();
  }}, 160);
}}
function resetParameterInputs() {{
  document.querySelectorAll('.parameter-input').forEach(input => {{
    input.value = input.dataset.default;
    input.setAttribute('value', input.value);
  }});
  try {{ localStorage.removeItem(parameterStorageKey); }} catch (_) {{}}
  try {{ const url = new URL(location.href); url.hash = ''; history.replaceState(null, '', url.href); }} catch (_) {{}}
  refreshInteractivePlots();
  setParameterStatus('Defaults restored; saved overrides cleared', 'ok');
}}
async function saveEditedHtml() {{
  persistParameterInputs(false);
  const originalInputs = Array.from(document.querySelectorAll('.parameter-input'));
  originalInputs.forEach(input => input.setAttribute('value', input.value));
  const clonedRoot = document.documentElement.cloneNode(true);
  const clonedInputs = Array.from(clonedRoot.querySelectorAll('.parameter-input'));
  clonedInputs.forEach((input, index) => {{
    const value = originalInputs[index].value;
    input.setAttribute('value', value);
    input.setAttribute('data-default', value);
  }});
  const source = '<!doctype html>\\n' + clonedRoot.outerHTML;
  const suggestedName = location.pathname.split('/').pop() || 'model-report.html';
  if ('showSaveFilePicker' in window) {{
    try {{
      const handle = await window.showSaveFilePicker({{suggestedName, types:[{{description:'HTML document', accept:{{'text/html':['.html']}}}}]}});
      const writable = await handle.createWritable();
      await writable.write(source);
      await writable.close();
      setParameterStatus('HTML copy saved with current values', 'ok');
      return;
    }} catch (error) {{ if (error.name === 'AbortError') return; }}
  }}
  const blob = new Blob([source], {{type:'text/html;charset=utf-8'}});
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url; link.download = suggestedName; document.body.appendChild(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  setParameterStatus('HTML copy downloaded with current values', 'ok');
}}

const SVG_NS = 'http://www.w3.org/2000/svg';
function svgNode(name, attributes={{}}, text='') {{
  const node = document.createElementNS(SVG_NS, name);
  Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, String(value)));
  if (text) node.textContent = text;
  return node;
}}
function parameterNumber(name, fallback) {{
  const input = document.querySelector('.parameter-input[data-param="' + name + '"]');
  const value = input ? Number(input.value) : fallback;
  return Number.isFinite(value) ? value : fallback;
}}
function cMul(a, b) {{ return {{re:a.re*b.re-a.im*b.im, im:a.re*b.im+a.im*b.re}}; }}
function cSub(a, b) {{ return {{re:a.re-b.re, im:a.im-b.im}}; }}
function cDiv(a, b) {{
  const denominator = b.re*b.re + b.im*b.im;
  return {{re:(a.re*b.re+a.im*b.im)/denominator, im:(a.im*b.re-a.re*b.im)/denominator}};
}}
function liveTransferData() {{
  const md = parameterNumber('reduced_drive_mass', 59.0);
  const ms = parameterNumber('reduced_stage_mass', 0.60);
  const km = parameterNumber('magnetic_stiffness', 1.20e8);
  const kax = parameterNumber('reduced_axial_stiffness', 1.14e7);
  const cax = parameterNumber('axial_damping', 55.0);
  const zeta = parameterNumber('electromagnetic_zeta', 0.50);
  if (!(md>0 && ms>0 && km>0 && kax>0 && cax>=0 && zeta>=0))
    throw new Error('Masses and stiffnesses must be positive; damping values must be non-negative.');
  const cm = 2*zeta*Math.sqrt(km*md);
  const frequencies=[], drive=[], stage=[], rotorStage=[];
  const count=560, logMin=Math.log10(100), logMax=Math.log10(3000);
  for (let i=0; i<count; i++) {{
    const frequency = Math.pow(10, logMin + (logMax-logMin)*i/(count-1));
    const omega = 2*Math.PI*frequency;
    const a = {{re:km+kax-md*omega*omega, im:omega*(cm+cax)}};
    const b = {{re:-kax, im:-omega*cax}};
    const d = {{re:kax-ms*omega*omega, im:omega*cax}};
    const determinant = cSub(cMul(a,d), cMul(b,b));
    const gd = cDiv(cMul({{re:km,im:0}},d), determinant);
    const gs = cDiv(cMul({{re:-km,im:0}},b), determinant);
    frequencies.push(frequency); drive.push(gd); stage.push(gs); rotorStage.push(cDiv(gs,gd));
  }}
  const qa=md*ms, qb=md*kax+ms*(km+kax), qc=km*kax;
  const discriminant=Math.max(qb*qb-4*qa*qc,0);
  const roots=[(qb-Math.sqrt(discriminant))/(2*qa),(qb+Math.sqrt(discriminant))/(2*qa)];
  const modes=roots.map(value => Math.sqrt(Math.max(value,0))/(2*Math.PI));
  return {{frequencies, drive, stage, rotorStage, modes, md, ms, km, kax, cax, zeta}};
}}
function unwrapPhases(values) {{
  const phases=[]; let previous=null, offset=0;
  values.forEach(value => {{
    let phase=Math.atan2(value.im,value.re)*180/Math.PI;
    if (previous!==null) {{
      while (phase+offset-previous>180) offset-=360;
      while (phase+offset-previous<-180) offset+=360;
    }}
    phase+=offset; phases.push(phase); previous=phase;
  }});
  return phases;
}}
function drawLiveBode(svgId, data, phasePlot=false) {{
  const svg=document.getElementById(svgId); if (!svg) return;
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const width=760, height=360, left=64, right=18, top=22, bottom=48;
  const plotWidth=width-left-right, plotHeight=height-top-bottom;
  const xMin=Math.log10(100), xMax=Math.log10(3000);
  const yMin=phasePlot?-370:-100, yMax=phasePlot?30:40;
  const mapX=f => left+(Math.log10(f)-xMin)/(xMax-xMin)*plotWidth;
  const mapY=y => top+(yMax-y)/(yMax-yMin)*plotHeight;
  const xTicks=[100,200,500,1000,2000,3000];
  const yTicks=phasePlot?[-360,-270,-180,-90,0]:[-100,-80,-60,-40,-20,0,20,40];
  xTicks.forEach(tick => {{
    const x=mapX(tick); svg.appendChild(svgNode('line',{{x1:x,y1:top,x2:x,y2:top+plotHeight,stroke:'#cbd3da','stroke-width':0.7}}));
    svg.appendChild(svgNode('text',{{x:x,y:height-22,'text-anchor':'middle','font-size':11,fill:'currentColor'}},tick>=1000?(tick/1000)+'k':String(tick)));
  }});
  yTicks.forEach(tick => {{
    const y=mapY(tick); svg.appendChild(svgNode('line',{{x1:left,y1:y,x2:left+plotWidth,y2:y,stroke:'#cbd3da','stroke-width':0.7}}));
    svg.appendChild(svgNode('text',{{x:left-8,y:y+4,'text-anchor':'end','font-size':11,fill:'currentColor'}},String(tick)));
  }});
  svg.appendChild(svgNode('rect',{{x:left,y:top,width:plotWidth,height:plotHeight,fill:'none',stroke:'#697680','stroke-width':1}}));
  const series=[
    {{name:'Xd / Xcmd',color:'#b23a48',values:data.drive}},
    {{name:'Xs / Xcmd',color:'#277da1',values:data.stage}},
    {{name:'Xs / Xd',color:'#218c74',values:data.rotorStage}}
  ];
  series.forEach((item,index) => {{
    const values=phasePlot?unwrapPhases(item.values):item.values.map(value => 20*Math.log10(Math.max(Math.hypot(value.re,value.im),1e-15)));
    let path='';
    values.forEach((value,i) => {{ const x=mapX(data.frequencies[i]), y=mapY(Math.max(yMin,Math.min(yMax,value))); path+=(i?'L':'M')+x.toFixed(2)+' '+y.toFixed(2)+' '; }});
    svg.appendChild(svgNode('path',{{d:path,fill:'none',stroke:item.color,'stroke-width':2}}));
    const legendX=left+12, legendY=top+16+index*17;
    svg.appendChild(svgNode('line',{{x1:legendX,y1:legendY-4,x2:legendX+25,y2:legendY-4,stroke:item.color,'stroke-width':2.4}}));
    svg.appendChild(svgNode('text',{{x:legendX+32,y:legendY,'font-size':11,fill:'currentColor'}},item.name));
  }});
  data.modes.forEach((mode,index) => {{
    const x=mapX(mode); svg.appendChild(svgNode('line',{{x1:x,y1:top,x2:x,y2:top+plotHeight,stroke:'#555','stroke-width':1,'stroke-dasharray':index?'2 3':'6 4'}}));
    svg.appendChild(svgNode('text',{{x:x+4,y:top+plotHeight-7,'font-size':10,fill:'currentColor'}},mode.toFixed(1)+' Hz'));
  }});
  svg.appendChild(svgNode('text',{{x:left+plotWidth/2,y:height-4,'text-anchor':'middle','font-size':12,fill:'currentColor'}},'Frequency (Hz)'));
  svg.appendChild(svgNode('text',{{x:15,y:top+plotHeight/2,transform:'rotate(-90 15 '+(top+plotHeight/2)+')','text-anchor':'middle','font-size':12,fill:'currentColor'}},phasePlot?'Phase (deg)':'Magnitude (dB)'));
}}
function refreshInteractivePlots() {{
  const panel=document.querySelector('[data-live-transfer-plots]'); if (!panel) return;
  if (!panel.dataset.initialized) {{
    panel.innerHTML='<div id="live-model-summary" class="live-summary"></div><div class="live-plot-grid"><div class="live-plot-card"><h4>Live magnitude</h4><svg id="live-bode-magnitude" viewBox="0 0 760 360" role="img" aria-label="Live transfer-function magnitude"></svg></div><div class="live-plot-card"><h4>Live phase</h4><svg id="live-bode-phase" viewBox="0 0 760 360" role="img" aria-label="Live transfer-function phase"></svg></div></div>';
    panel.dataset.initialized='true';
  }}
  try {{
    const data=liveTransferData();
    document.getElementById('live-model-summary').textContent='Live values: md='+data.md.toPrecision(5)+' kg, ms='+data.ms.toPrecision(5)+' kg, Km='+data.km.toPrecision(5)+' N/m, kax='+data.kax.toPrecision(5)+' N/m, zeta_m='+data.zeta.toPrecision(4)+' · modes '+data.modes.map(value=>value.toFixed(2)+' Hz').join(', ');
    drawLiveBode('live-bode-magnitude',data,false); drawLiveBode('live-bode-phase',data,true);
  }} catch (error) {{ document.getElementById('live-model-summary').textContent='Live plot error: '+error.message; }}
}}
document.addEventListener('DOMContentLoaded', () => {{
  let saved = {{}};
  let browserStorageLoaded = true;
  try {{ saved = JSON.parse(localStorage.getItem(parameterStorageKey) || '{{}}'); }} catch (_) {{ saved = {{}}; browserStorageLoaded = false; }}
  saved = Object.assign(saved, readHashParameterValues());
  document.querySelectorAll('.parameter-input').forEach(input => {{
    if (Object.prototype.hasOwnProperty.call(saved, input.dataset.param)) input.value = saved[input.dataset.param];
    input.setAttribute('value', input.value);
    input.addEventListener('input', () => {{ input.setAttribute('value', input.value); scheduleParameterUpdate(); }});
  }});
  refreshInteractivePlots();
  if (Object.keys(saved).length) setParameterStatus('Saved variables restored', 'ok');
  else if (browserStorageLoaded) setParameterStatus('Defaults loaded · auto-save active', 'ok');
  else setParameterStatus('URL saving active; browser storage unavailable', 'warn');
}});
</script>
</body></html>"""


def render_document(markdown_path: Path) -> Path:
    output = markdown_path.with_suffix(".html")
    output.write_text(html_page(markdown_path), encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-summary-update", action="store_true",
                        help="Render without refreshing the generated metrics table")
    args = parser.parse_args()
    if not DESCRIPTION_MD.exists() or not DERIVATION_MD.exists():
        raise FileNotFoundError("Both Markdown source documents must exist before building")
    ASSET_DIR.mkdir(exist_ok=True)
    gms_audit = validate_gms_partition()
    constants = physical_constants()
    frequencies, bode, linear_metrics = frequency_responses()
    times, command, time_data, time_metrics = time_responses(constants)
    convergence = gms_step_halving_convergence(constants, times, time_data)
    presliding = presliding_responses(constants)
    verification = full_reduced_verification(frequencies, constants)
    case_paths = plot_case_responses(
        frequencies, bode, times, command, time_data, constants, time_metrics)
    comparison_path = plot_pairwise_comparison(
        frequencies, bode, times, command, time_data, time_metrics)
    presliding_path = plot_presliding_memory(presliding)
    diagram_path = plot_kinematic_diagram()
    verification_path = plot_full_reduced_verification(frequencies, verification)
    position_path = plot_position_dependence()
    resonance_path = plot_stepper_resonance_visibility()
    rotor_stage_path = plot_rotor_stage_transfer_functions(frequencies)
    for obsolete_name in ("bode_all_cases.svg", "step_tracking_all_cases.svg"):
        obsolete_path = ASSET_DIR / obsolete_name
        if obsolete_path.exists():
            obsolete_path.unlink()
    if not args.skip_summary_update:
        update_generated_summary(generated_summary(linear_metrics, time_metrics, verification))
        update_generated_presliding_summary(generated_presliding_summary(presliding))
        update_generated_convergence_summary(generated_convergence_summary(convergence))
    description_html = render_document(DESCRIPTION_MD)
    derivation_html = render_document(DERIVATION_MD)
    for case_path in case_paths:
        print(f"Built {case_path.relative_to(ROOT)}")
    print(f"Built {comparison_path.relative_to(ROOT)}")
    print(f"Built {presliding_path.relative_to(ROOT)}")
    print(f"Built {diagram_path.relative_to(ROOT)}")
    print(f"Built {verification_path.relative_to(ROOT)}")
    print(f"Built {position_path.relative_to(ROOT)}")
    print(f"Built {resonance_path.relative_to(ROOT)}")
    print(f"Built {rotor_stage_path.relative_to(ROOT)}")
    print(f"Built {description_html.name}")
    print(f"Built {derivation_html.name}")
    print(f"GMS partition: sum(nu_i)={gms_audit['weight_sum']:.12f}; "
          + "; ".join(
              f"sum(k_i)_{site}={value:.6g}=sigma0_{site}"
              for site, value in gms_audit["stiffness_sums"].items()
          ))
    for key in CASES:
        modes = linear_metrics[key]["modes"]
        print(f"Case {key}: modes={modes[0]:.2f}, {modes[1]:.2f} Hz; "
              f"DC gain={linear_metrics[key]['dc_gain']:.6f}; "
              f"overshoot={time_metrics[key]['first_overshoot_pct']:.2f}%; "
              f"sequence RMS deviation={time_metrics[key]['rms_sequence_deviation_nm']:.2f} nm; "
              f"peak deviation={time_metrics[key]['max_abs_deviation_nm']:.2f} nm; "
              f"final-window RMS deviation={time_metrics[key]['rms_final_error_nm']:.2f} nm")
    print(f"Full/reduced residual: RMS={verification['rms_residual_nm']:.3f} nm; "
          f"peak={verification['peak_residual_nm']:.3f} nm")
    for key in ("A", "A2"):
        metric = presliding["metrics"][key]
        print(f"Presliding {key}: RMS={metric['whole_rms_nm']:.3f} nm; "
              f"return mismatch={metric['return_error_mismatch_nm']:.3f} nm; "
              f"force closure={metric['return_force_mismatch_N']:.6f} N; "
              f"final mean={metric['final_mean_nm']:.3f} nm")
    for key in ("A2", "B2", "C2"):
        result = convergence[key]
        rms = result["rms_nm"]
        print(f"Step halving {key}: RMS(10/5/2.5 us)="
              f"{rms[0]:.6f}/{rms[1]:.6f}/{rms[2]:.6f} nm; "
              f"fine relative change={result['fine_relative_pct']:.6f}%")
    print("Full-model modes below 3 kHz: " + ", ".join(
        f"{mode:.2f}" for mode in verification["full_modes"] if mode < 3000.0))


if __name__ == "__main__":
    main()

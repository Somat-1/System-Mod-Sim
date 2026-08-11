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
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"
DESCRIPTION_MD = ROOT / "ball_screw_stage_dynamic_derivation_v3.md"
DERIVATION_MD = ROOT / "Analytical_derivation_and_responses_v3.md"


# Executable defaults for the Revision 3 two-DOF reduction.
MODEL = {
    "lead": 1.0e-3,
    "rotor_teeth": 50,
    # Rated-current holding torque for the lower-current motor variant.
    "T_max": 0.060,
    # Published detent torque.  Phase zero places the report origin at a
    # stable detent equilibrium and enables the detent term in every model.
    "T_det": 0.005,
    "detent_phase": 0.0,
    "m_s": 0.60,
    "k_ax": 1.14e7,
    # Provisional: retained structural damping, not identified in the source.
    "c_ax": 55.0,
    # Provisional open-loop drive damping ratio.  Driver mode and tuning are
    # not recorded, so the report executes a low-damping baseline and shows a
    # sensitivity sweep rather than presenting one value as identified.
    "zeta_m": 0.05,
    # Conservative external STEP/DIR resolution.  The TMC2209 accepts up to
    # 64 microsteps per full step and may interpolate internally to 256.
    "microstep_divisor": 64,
}


# Revision 3 full-model values.  These are deliberately separate from MODEL:
# MODEL is the validated two-DOF reduction, whereas FULL retains all ten
# coordinates named in the source document.  Values not measured in the source
# are surfaced as highlighted assumptions in the Markdown documents.
FULL = {
    # Component values.  Screw inertia and axial masses are derived below.
    "J_m": 9.00e-7,
    "J_c": 1.18e-6,
    "screw_length": 0.320,
    "screw_diameter": 8.00e-3,
    "screw_density": 7850.0,
    "m_n": 0.050,
    "m_stage": 0.550,
    # Datasheet series stiffness is 1.2 N m/deg = 68.7549 N m/rad.
    # Two equal half-springs must each be twice the series value.
    "k_c_series": 1.2 * 180.0 / np.pi,
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
    # Differential nut-contact microslip.  Its first GMS element yields at
    # 0.25*F_s/sigma0 = 0.20 um, so this port can express actual partial slip.
    "n": {"sigma0": 2.00e6, "sigma1": 5.0, "sigma2": 0.25,
          "F_s": 1.6, "F_c": 1.2, "v_s": 2.0e-4, "delta": 1.0, "C_gms": 5.0e3},
    # Gross ball-nut rolling drag acts on common drivetrain motion, not on the
    # differential elastic-deformation rate represented by site n.
    "r": {"sigma0": 2.00e6, "sigma1": 5.0, "sigma2": 0.25,
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


# Nested reversals for the dedicated memory experiment.  Counts use the
# conservative 64-microstep STEP/DIR quantum.  The outer excursion is large
# enough to yield two guideway GMS elements, while the force signal remains the
# primary discriminator.
PRESLIDING_LEVELS = np.array([0, 48, 12, 42, 12, 48, 0, -46, -12, -40, -12, -46, 0], dtype=float)
PRESLIDING_START = 0.005
PRESLIDING_HOLD = 0.010
PRESLIDING_RETURN_PAIRS = ((1, 5), (2, 4), (7, 11), (8, 10))


CASES = OrderedDict([
    ("0", {"label": "Case 0: frictionless", "sites": (), "friction": "none", "color": "#252525", "ls": "--"}),
    ("A", {"label": "Case A: drivetrain + guideway / LuGre", "sites": ("d", "g"), "friction": "lugre", "color": "#277da1", "ls": "-"}),
    ("A2", {"label": "Case A2: drivetrain + guideway / GMS", "sites": ("d", "g"), "friction": "gms", "color": "#70b7cf", "ls": "--"}),
    ("B", {"label": "Case B: drivetrain + nut rolling/microslip / LuGre", "sites": ("d", "r", "n"), "friction": "lugre", "color": "#e07a15", "ls": "-"}),
    ("B2", {"label": "Case B2: drivetrain + nut rolling/microslip / GMS", "sites": ("d", "r", "n"), "friction": "gms", "color": "#f5b35f", "ls": "--"}),
    ("C", {"label": "Case C: all friction ports / LuGre", "sites": ("d", "g", "r", "n"), "friction": "lugre", "color": "#218c74", "ls": "-"}),
    ("C2", {"label": "Case C2: all friction ports / GMS", "sites": ("d", "g", "r", "n"), "friction": "gms", "color": "#72c9ad", "ls": "--"}),
])

PAIRS = (("A", "A2"), ("B", "B2"), ("C", "C2"))


H = {
    "g": np.array([0.0, 1.0]),
    "n": np.array([1.0, -1.0]),
    "r": np.array([1.0, 0.0]),
    "d": np.array([1.0, 0.0]),
}

SITE_KEYS = tuple(FRICTION)
LUGRE_INDEX = {site: 4 + index for index, site in enumerate(SITE_KEYS)}
GMS_BASE = 4 + len(SITE_KEYS)
GMS_START = {site: GMS_BASE + index * GMS_N for index, site in enumerate(SITE_KEYS)}
STATE_SIZE = GMS_BASE + len(SITE_KEYS) * GMS_N


def save_svg(fig: plt.Figure, output: Path) -> None:
    """Write deterministic SVG text without Matplotlib's trailing spaces."""
    fig.savefig(output, format="svg", bbox_inches="tight")
    svg = output.read_text(encoding="utf-8")
    output.write_text(
        "\n".join(line.rstrip() for line in svg.splitlines()) + "\n",
        encoding="utf-8",
    )


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


def validate_case_topology() -> None:
    """Reject regressions in the friction-port allocation."""
    expected_sites = {
        "0": set(),
        "A": {"d", "g"}, "A2": {"d", "g"},
        "B": {"d", "r", "n"}, "B2": {"d", "r", "n"},
        "C": {"d", "g", "r", "n"}, "C2": {"d", "g", "r", "n"},
    }
    for key, expected in expected_sites.items():
        actual = set(CASES[key]["sites"])
        if actual != expected:
            raise ValueError(f"Case {key} friction sites are {sorted(actual)}, expected {sorted(expected)}")
    nut_first_yield = float(np.min(
        GMS_WEIGHTS * FRICTION["n"]["F_s"] /
        (GMS_STIFFNESS_FRACTIONS * FRICTION["n"]["sigma0"])
    ))
    if not np.isclose(nut_first_yield, 0.20e-6, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"Nut microslip first yield is {nut_first_yield:.6g} m, expected 0.20 um")


def physical_constants() -> dict[str, float]:
    lead = MODEL["lead"]
    teeth = MODEL["rotor_teeth"]
    r = lead / (2.0 * np.pi)
    kappa = 2.0 * np.pi * teeth / lead
    t_max = MODEL["T_max"]
    t_det = MODEL["T_det"]
    k_m = teeth * t_max / r**2
    k_det = 4.0 * teeth * t_det * np.cos(MODEL["detent_phase"]) / r**2
    k_drive = k_m + k_det
    f_max = t_max / r
    full_step = lead / (4.0 * teeth)
    component = component_parameters()
    m_d = component["J_total"] / r**2
    c_m = 2.0 * MODEL["zeta_m"] * np.sqrt(k_drive * m_d)
    return {
        "r": r,
        "kappa": kappa,
        "T_max": t_max,
        "T_det": t_det,
        "F_max": f_max,
        "K_m": k_m,
        "K_det": k_det,
        "K_drive": k_drive,
        "m_d": m_d,
        "c_m": c_m,
        "full_step": full_step,
        "quarter_step": full_step / 4.0,
        "command_step": full_step / MODEL["microstep_divisor"],
        "interpolated_step": full_step / 256.0,
    }


def component_parameters() -> dict[str, float]:
    """Derive screw inertia and lumped masses from the 0.320 m component."""
    p = dict(FULL)
    radius = 0.5 * p["screw_diameter"]
    area = np.pi * radius**2
    screw_mass = p["screw_density"] * area * p["screw_length"]
    screw_inertia = 0.5 * screw_mass * radius**2
    p["screw_mass"] = screw_mass
    p["screw_inertia"] = screw_inertia
    for key in ("J_s1", "J_s2", "J_s3"):
        p[key] = screw_inertia / 3.0
    for key in ("m_b", "m_e", "m_f"):
        p[key] = screw_mass / 3.0
    p["k_c1"] = 2.0 * p["k_c_series"]
    p["k_c2"] = 2.0 * p["k_c_series"]
    p["J_total"] = p["J_m"] + p["J_c"] + screw_inertia
    return p


def full_parameters() -> dict[str, float]:
    """Return the ten-DOF parameters and close the measured axial compliance."""
    p = component_parameters()
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
    p["m_d_reflected"] = p["J_total"] / r**2
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

    constants = physical_constants()
    k_m_rot = constants["K_m"] * r**2
    k_det_rot = constants["K_det"] * r**2
    stiffness[0, 0] += k_m_rot + k_det_rot
    # The same damping repair validated in Rev 2, expressed in rotational units.
    damping[0, 0] += 2.0 * MODEL["zeta_m"] * np.sqrt(
        (k_m_rot + k_det_rot) * p["J_total"]
    )

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
        held_command = command_position(times[i] + 0.5 * dt, constants["command_step"])
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
    command = np.array([command_position(t, constants["command_step"]) for t in times])
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
    constants = physical_constants()
    m_d, m_s = constants["m_d"], MODEL["m_s"]
    k_ax, k_m, k_det, c_ax = (
        MODEL["k_ax"], constants["K_m"], constants["K_det"], MODEL["c_ax"]
    )
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    mass = np.diag([m_d, m_s])
    damping = c_ax * coupling + constants["c_m"] * np.outer(H["d"], H["d"])
    stiffness = np.array([[k_m + k_det + k_ax, -k_ax], [-k_ax, k_ax]], dtype=float)
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
        tangent_dc_gain = float((np.linalg.solve(stiffness, input_vector))[1])
        first_yield = min(
            (float(np.min(GMS_WEIGHTS * FRICTION[site]["F_s"] /
                          (GMS_STIFFNESS_FRACTIONS * FRICTION[site]["sigma0"])))
             for site in case["sites"]),
            default=np.inf,
        )
        metrics[key] = {
            "modes": modes,
            "tangent_dc_gain": tangent_dc_gain,
            "first_yield_m": first_yield,
        }
    return frequencies, responses, metrics


def command_position(t: float, command_step: float) -> float:
    """Closed four-increment sequence at the configured STEP/DIR quantum."""
    if t < 0.005:
        return 0.0
    if t < 0.025:
        return command_step
    if t < 0.045:
        return 0.0
    if t < 0.065:
        return -command_step
    return 0.0


def presliding_command_position(t: float, microstep: float) -> float:
    """Nested back-and-forth reversals quantized to the STEP/DIR input."""
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
    command = command_position(t, constants["command_step"])
    if held_command is not None:
        command = held_command
    lag = constants["kappa"] * (command - x_d)
    magnetic_force = constants["F_max"] * np.sin(lag)
    detent_force = -(constants["T_det"] / constants["r"]) * np.sin(
        4.0 * constants["kappa"] * x_d + MODEL["detent_phase"]
    )
    electromagnetic_damping = constants["c_m"] * v_d
    axial_force = MODEL["k_ax"] * (x_d - x_s) + MODEL["c_ax"] * (v_d - v_s)

    velocities = {site: float(H[site] @ np.array([v_d, v_s])) for site in SITE_KEYS}
    forces = {site: 0.0 for site in SITE_KEYS}
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

    friction_generalized = sum(H[site] * forces[site] for site in SITE_KEYS)
    a_d = (magnetic_force + detent_force - electromagnetic_damping - axial_force
           - friction_generalized[0]) / constants["m_d"]
    a_s = (axial_force - friction_generalized[1]) / MODEL["m_s"]
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
        lambda t: command_position(t, constants["command_step"]),
        duration=duration, dt=dt,
    )


def friction_force_history(case: dict[str, object], states: np.ndarray,
                           site: str) -> np.ndarray:
    """Recover an executed site's constitutive force from integrated states."""
    if site not in case["sites"]:
        return np.zeros(states.shape[0])
    velocity = H[site][0] * states[:, 2] + H[site][1] * states[:, 3]
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
    microstep = constants["command_step"]
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
    command = np.array([command_position(t, constants["command_step"]) for t in times])
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
    command = np.array([command_position(t, constants["command_step"]) for t in times])
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
            "first_overshoot_pct": max(0.0, (first_peak / constants["command_step"] - 1.0) * 100.0),
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
                 f"Nonlinear magnetic and detent force; zeta_m={MODEL['zeta_m']:.2f}; "
                 f"each command increment = {constants['command_step'] * 1e9:.2f} nm at 1/{MODEL['microstep_divisor']} full step.",
                 ha="center", fontsize=8.5, color="#555555")
        fig.tight_layout(rect=(0.02, 0.075, 0.99, 0.95))
        output = ASSET_DIR / f"response_case_{key}.svg"
        save_svg(fig, output)
        plt.close(fig)
        outputs.append(output)
    return outputs


def plot_case_response_overlay(frequencies: np.ndarray,
                               responses: dict[str, np.ndarray]) -> Path:
    """Overlay every case and quantify the small differences near resonance."""
    fig = plt.figure(figsize=(12.2, 9.0))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.08, 1.0), hspace=0.30, wspace=0.25)
    ax_full = fig.add_subplot(grid[0, :])
    ax_zoom = fig.add_subplot(grid[1, 0])
    ax_delta = fig.add_subplot(grid[1, 1])

    magnitudes = {
        key: 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))
        for key, response in responses.items()
    }
    short_labels = {
        "0": "0: frictionless",
        "A": "A: guideway LuGre",
        "A2": "A2: guideway GMS",
        "B": "B: nut LuGre",
        "B2": "B2: nut GMS",
        "C": "C: both LuGre",
        "C2": "C2: both GMS",
    }
    for key, case in CASES.items():
        for axis in (ax_full, ax_zoom):
            axis.semilogx(
                frequencies, magnitudes[key], color=case["color"],
                linestyle=case["ls"], linewidth=1.8,
                label=short_labels[key] if axis is ax_full else None,
            )

    ax_full.axhline(0.0, color="#888888", linewidth=0.7)
    ax_full.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
    ax_full.set_ylim(-90.0, 30.0)
    ax_full.set_title("All command-to-stage Bode responses")
    ax_full.set_ylabel("Magnitude (dB)")
    ax_full.legend(loc="lower left", ncol=2, fontsize=8)

    zoom_mask = (frequencies >= 620.0) & (frequencies <= 830.0)
    zoom_indices = np.flatnonzero(zoom_mask)
    baseline_index = zoom_indices[np.argmax(magnitudes["0"][zoom_mask])]
    baseline_frequency = float(frequencies[baseline_index])
    annotations = (
        ("0", "0", (682.0, 14.2)),
        ("A2", "A/A2", (720.0, 15.2)),
        ("B2", "B/B2", (750.0, 12.8)),
        ("C2", "C/C2", (796.0, 14.8)),
    )
    for key, label, text_position in annotations:
        peak_index = zoom_indices[np.argmax(magnitudes[key][zoom_mask])]
        peak_frequency = float(frequencies[peak_index])
        peak_magnitude = float(magnitudes[key][peak_index])
        if key == "0":
            note = f"{label}: {peak_frequency:.0f} Hz"
        else:
            shift_hz = peak_frequency - baseline_frequency
            shift_pct = 100.0 * shift_hz / baseline_frequency
            note = f"{label}: {peak_frequency:.0f} Hz\n+{shift_hz:.1f} Hz ({shift_pct:.1f}%)"
        ax_zoom.annotate(
            note, xy=(peak_frequency, peak_magnitude),
            xytext=text_position, textcoords="data", ha="center", fontsize=7.8,
            color=CASES[key]["color"],
            arrowprops={"arrowstyle": "->", "color": CASES[key]["color"], "lw": 0.8},
        )
    ax_zoom.set_xlim(620.0, 830.0)
    ax_zoom.set_ylim(2.0, 16.2)
    ax_zoom.set_title("Higher resonance: topology shifts the peak")
    ax_zoom.set_xlabel("Frequency (Hz)")
    ax_zoom.set_ylabel("Magnitude (dB)")

    for lugre_key, gms_key in PAIRS:
        difference = magnitudes[gms_key] - magnitudes[lugre_key]
        ax_delta.semilogx(
            frequencies, difference, color=CASES[gms_key]["color"],
            linewidth=1.9, label=f"{lugre_key}/{gms_key}",
        )
        maximum_index = int(np.argmax(np.abs(difference)))
        maximum = float(difference[maximum_index])
        maximum_frequency = float(frequencies[maximum_index])
        ax_delta.plot(maximum_frequency, maximum, "o", color=CASES[gms_key]["color"], ms=4)
        ax_delta.annotate(
            f"{abs(maximum):.2f} dB at {maximum_frequency:.0f} Hz",
            xy=(maximum_frequency, maximum), xytext=(7, 7), textcoords="offset points",
            fontsize=7.6, color=CASES[gms_key]["color"],
        )
    ax_delta.axhline(0.0, color="#888888", linewidth=0.7)
    ax_delta.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
    ax_delta.set_ylim(-0.1, 1.2)
    ax_delta.set_title("GMS minus LuGre magnitude")
    ax_delta.set_xlabel("Frequency (Hz)")
    ax_delta.set_ylabel("Difference (dB)")
    ax_delta.legend(loc="upper left", fontsize=8)

    for axis in (ax_full, ax_zoom, ax_delta):
        axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)

    fig.suptitle("Revision 3 response comparison", fontsize=15, fontweight="bold")
    fig.text(
        0.5, 0.012,
        "Matched LuGre and GMS cases share presliding stiffness. Their visible Bode gap comes from tangent damping.",
        ha="center", fontsize=8.4, color="#555555",
    )
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.08, top=0.91,
                        hspace=0.34, wspace=0.26)
    output = ASSET_DIR / "lugre_gms_pairwise_comparison.svg"
    save_svg(fig, output)
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
    ax_metrics.text(
        0.02, 0.03,
        f"Primary force metric: return mismatch "
        f"{metrics['A']['return_force_mismatch_N']:.4f} N LuGre, "
        f"{metrics['A2']['return_force_mismatch_N']:.4f} N GMS",
        transform=ax_metrics.transAxes, fontsize=7.6, color="#555555",
        bbox={"boxstyle": "round,pad=0.24", "facecolor": "white", "edgecolor": "#c9cfd4"},
    )
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
    fig.suptitle("Force-instrumented partial-slip experiment: LuGre versus GMS",
                 fontsize=15, fontweight="bold")
    fig.text(
        0.5, 0.012,
        f"1 STEP/DIR quantum = {experiment['microstep'] * 1e9:.2f} nm = 1/{MODEL['microstep_divisor']} full step; "
        f"peak friction = {max_force:.3f} N ({macro_fraction:.1f}% of macro breakaway). "
        "Markers are 2 ms plateau-end means.",
        ha="center", fontsize=8.4, color="#555555",
    )
    fig.tight_layout(rect=(0.02, 0.045, 0.99, 0.95), h_pad=2.0, w_pad=1.5)
    output = ASSET_DIR / "presliding_memory_comparison.svg"
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_kinematic_diagram() -> Path:
    """Render the full topology, retained reduction, and rejected one-DOF collapse."""
    constants = physical_constants()
    parameters = full_parameters()
    r = constants["r"]
    drive_color = "#dceef6"
    stage_color = "#dff2ea"
    dropped_color = "#eeeeee"
    spring_color = "#d97800"
    discarded_spring_color = "#8a8a8a"
    friction_color = "#b23a48"
    damping_color = "#6a4c93"
    rigid_color = "#39434d"

    fig = plt.figure(figsize=(16.2, 11.2))
    grid = fig.add_gridspec(2, 3, height_ratios=(1.78, 1.0), width_ratios=(1.0, 1.0, 0.84),
                            hspace=0.20, wspace=0.14)
    full_ax = fig.add_subplot(grid[0, :])
    reduced_ax = fig.add_subplot(grid[1, :2])
    rejected_ax = fig.add_subplot(grid[1, 2])
    for axis in (full_ax, reduced_ax, rejected_ax):
        axis.axis("off")

    def node(ax: plt.Axes, x: float, y: float, label: str, index: int | None,
             color: str, subtitle: str = "", width: float = 0.92) -> None:
        box = FancyBboxPatch((x - width / 2.0, y - 0.31), width, 0.62,
                             boxstyle="round,pad=0.04", facecolor=color,
                             edgecolor=rigid_color, linewidth=1.2, zorder=4)
        ax.add_patch(box)
        ax.text(x, y + 0.06, label, ha="center", va="center", fontsize=9.4, zorder=5)
        lower = subtitle if subtitle else (f"q{index}" if index is not None else "")
        ax.text(x, y - 0.19, lower, ha="center", va="center", fontsize=7.2,
                color="#59636d", zorder=5)

    def spring(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
               color: str = spring_color, linewidth: float = 1.7, amplitude: float = 0.075,
               linestyle: str = "-") -> None:
        p0, p1 = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
        delta = p1 - p0
        length = np.linalg.norm(delta)
        unit = delta / length
        normal = np.array([-unit[1], unit[0]])
        points = [p0, p0 + 0.15 * delta]
        for i in range(9):
            fraction = 0.18 + i * 0.08
            points.append(p0 + fraction * delta + (amplitude if i % 2 == 0 else -amplitude) * normal)
        points.extend([p0 + 0.85 * delta, p1])
        points = np.asarray(points)
        ax.plot(points[:, 0], points[:, 1], color=color, linewidth=linewidth,
                linestyle=linestyle, solid_capstyle="round", zorder=2)

    def dashpot(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
                color: str = damping_color, linewidth: float = 1.5, width: float = 0.10) -> None:
        p0, p1 = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
        delta = p1 - p0
        length = np.linalg.norm(delta)
        unit = delta / length
        normal = np.array([-unit[1], unit[0]])
        a, b, c, d = (p0 + f * delta for f in (0.28, 0.40, 0.65, 0.78))
        ax.plot(*zip(p0, a), color=color, linewidth=linewidth)
        ax.plot(*zip(a - width * normal, a + width * normal), color=color, linewidth=linewidth)
        ax.plot(*zip(a, c), color=color, linewidth=linewidth)
        ax.plot(*zip(b - width * normal, d - width * normal), color=color, linewidth=linewidth)
        ax.plot(*zip(b + width * normal, d + width * normal), color=color, linewidth=linewidth)
        ax.plot(*zip(d - width * normal, d + width * normal), color=color, linewidth=linewidth)
        ax.plot(*zip(d, p1), color=color, linewidth=linewidth)

    def friction(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
                 label: str, fontsize: float = 7.4, tag: str = "",
                 defined_only: bool = False) -> None:
        x0, y0 = start
        x1, y1 = end
        horizontal = abs(x1 - x0) >= abs(y1 - y0)
        cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        if horizontal:
            block_w, block_h = min(0.62, abs(x1 - x0) * 0.43), 0.23
            ax.plot([x0, cx - block_w / 2.0], [y0, cy], color=friction_color, linewidth=1.5)
            ax.plot([cx + block_w / 2.0, x1], [cy, y1], color=friction_color, linewidth=1.5)
        else:
            block_w, block_h = 0.34, min(0.48, abs(y1 - y0) * 0.43)
            ax.plot([x0, cx], [y0, cy - block_h / 2.0], color=friction_color, linewidth=1.5)
            ax.plot([cx, x1], [cy + block_h / 2.0, y1], color=friction_color, linewidth=1.5)
        block = FancyBboxPatch((cx - block_w / 2.0, cy - block_h / 2.0), block_w, block_h,
                               boxstyle="round,pad=0.025", facecolor="#f8dce1",
                               edgecolor=friction_color, linewidth=1.1, zorder=3,
                               linestyle="--" if defined_only else "-")
        ax.add_patch(block)
        ax.text(cx, cy, label, ha="center", va="center", fontsize=fontsize,
                color="#8d2936", zorder=4)
        if tag:
            ax.text(cx, cy - block_h / 2.0 - 0.10, tag, ha="center", va="top",
                    fontsize=5.8, color="#8d2936", zorder=4)

    def ground(ax: plt.Axes, x: float, y: float, width: float = 0.48) -> None:
        ax.plot([x - width / 2.0, x + width / 2.0], [y, y], color=rigid_color, linewidth=2.0)
        for offset in np.linspace(-width / 2.0, width / 2.0, 5):
            ax.plot([x + offset, x + offset - 0.08], [y, y - 0.10], color=rigid_color, linewidth=0.8)

    def wall(ax: plt.Axes, x: float, y: float, height: float = 0.72) -> None:
        ax.plot([x, x], [y - height / 2.0, y + height / 2.0], color=rigid_color, linewidth=2.2)
        for offset in np.linspace(-height / 2.0, height / 2.0, 5):
            ax.plot([x, x - 0.10], [y + offset, y + offset - 0.08], color=rigid_color, linewidth=0.8)

    components = component_parameters()
    coupling_reflected = parameters["k_c_series"] / r**2
    torsion_reflected = parameters["k_theta_a"] / r**2
    shares = {
        "brg": 100.0 * MODEL["k_ax"] / parameters["k_brg"],
        "sha": 100.0 * MODEL["k_ax"] / parameters["k_sha"],
        "ball": 100.0 * MODEL["k_ax"] / parameters["k_ball"],
        "mnt": 100.0 * MODEL["k_ax"] / parameters["k_mnt"],
    }
    reduced_mass, _, reduced_stiffness, _ = linear_matrices((), "none")
    reduced_modes = _linear_modes(reduced_mass, reduced_stiffness)
    full_mass, _, full_stiffness, _, _ = full_linear_matrices()
    eigenvalues, eigenvectors = np.linalg.eig(np.linalg.solve(full_mass, full_stiffness))
    order = np.argsort(np.real(eigenvalues))
    modal_frequencies = np.sqrt(np.maximum(np.real(eigenvalues[order]), 0.0)) / (2.0 * np.pi)
    modal_vectors = np.real(eigenvectors[:, order])
    participation = np.diag(full_mass)[:, None] * modal_vectors**2
    participation /= np.maximum(np.sum(participation, axis=0, keepdims=True), 1e-30)

    def first_stub_mode(dof: int) -> float:
        candidates = np.where((modal_frequencies > 900.0) & (participation[dof] > 0.15))[0]
        return float(modal_frequencies[candidates[0]]) if candidates.size else float("nan")

    torsional_stub_mode = first_stub_mode(4)
    axial_stub_mode = first_stub_mode(7)

    full_ax.set_xlim(0.0, 15.6)
    full_ax.set_ylim(0.0, 7.75)
    full_ax.text(7.8, 7.50, "Ten-DOF physical topology and reduction map",
                 ha="center", fontsize=15, fontweight="bold")
    stations = ((2.0, "motor"), (3.55, "coupling"), (5.15, "bearing"),
                (6.55, "screw"), (8.35, "nut"), (10.55, "stage"),
                (11.55, "guideways"))
    for x_pos, label in stations:
        full_ax.plot([x_pos, x_pos], [0.24, 6.86], color="#dfe3e6", linewidth=0.7,
                     linestyle=(0, (2, 4)), zorder=0)
        full_ax.text(x_pos, 7.06, label, ha="center", fontsize=8.2, color="#59636d",
                     bbox={"boxstyle": "round,pad=0.20", "facecolor": "#f5f6f7",
                           "edgecolor": "#d4d8dc"})

    # Commanded magnetic field and motor-side grounding.
    wall(full_ax, 0.48, 5.30, 0.95)
    full_ax.text(0.48, 5.94, r"moving datum $x_{cmd}$", ha="center", fontsize=8.3)
    node(full_ax, 2.0, 5.30, r"$\theta_m$", 1, drive_color,
         rf"q1; $J={parameters['J_m']:.2e}$")
    spring(full_ax, (0.50, 5.30), (1.52, 5.30), color="#c08a00")
    full_ax.text(1.02, 5.57, r"$K_m$", ha="center", fontsize=8.2, color="#8a6200")
    spring(full_ax, (1.72, 4.97), (1.72, 4.05), color="#c08a00", amplitude=0.06)
    ground(full_ax, 1.72, 3.91, 0.42)
    full_ax.text(1.08, 4.43, r"$K_{det}(x_d)$", fontsize=7.4, color="#8a6200")
    dashpot(full_ax, (2.00, 4.97), (2.00, 4.05), width=0.08)
    ground(full_ax, 2.00, 3.91, 0.42)
    friction(full_ax, (2.33, 4.97), (2.33, 4.05), r"$T_{mb}$", 6.7)
    ground(full_ax, 2.33, 3.91, 0.42)

    # Rotational masses migrate to m_d; their internal spring compliance is discarded.
    node(full_ax, 3.55, 5.30, r"$\theta_c$", 2, drive_color,
         rf"q2; $J={parameters['J_c']:.2e}$")
    node(full_ax, 5.15, 5.30, r"$\theta_{{s1}}$", 3, drive_color,
         rf"q3; $J={parameters['J_s1']:.2e}$")
    node(full_ax, 8.35, 5.30, r"$\theta_{{s2}}$", 4, drive_color,
         rf"q4; $J={parameters['J_s2']:.2e}$")
    spring(full_ax, (2.48, 5.45), (3.07, 5.45), color=discarded_spring_color,
           amplitude=0.055, linestyle="--")
    friction(full_ax, (2.48, 5.08), (3.07, 5.08), r"$T_{h1}$", 6.5)
    spring(full_ax, (4.03, 5.45), (4.67, 5.45), color=discarded_spring_color,
           amplitude=0.055, linestyle="--")
    friction(full_ax, (4.03, 5.08), (4.67, 5.08), r"$T_{h2}$", 6.5)
    spring(full_ax, (5.63, 5.30), (7.87, 5.30), color=discarded_spring_color,
           amplitude=0.065, linestyle="--")
    full_ax.text(6.75, 5.56, r"$k_{\theta a}$", ha="center", fontsize=8.0, color="#666666")
    friction(full_ax, (5.15, 4.97), (5.15, 4.05), r"$T_{brg}$", 6.6)
    ground(full_ax, 5.15, 3.91, 0.44)

    node(full_ax, 8.35, 6.42, r"$\theta_{{s3}}$", 5, drive_color,
         rf"q5; $J={parameters['J_s3']:.2e}$")
    spring(full_ax, (8.35, 5.63), (8.35, 6.09), color=discarded_spring_color,
           amplitude=0.045, linestyle="--")
    full_ax.text(8.70, 5.92, r"$k_{\theta b}$", fontsize=7.4, color="#666666")

    # The four full drivetrain friction torques aggregate into the reduced F_f,d port.
    aggregate_sources = (2.33, 2.78, 4.35, 5.15)
    for x_pos in aggregate_sources:
        full_ax.plot([x_pos, x_pos], [4.02 if x_pos in (2.33, 5.15) else 4.91, 3.67],
                     color="#9a6aae", linewidth=0.8, linestyle="--", zorder=1)
    full_ax.plot([aggregate_sources[0], aggregate_sources[-1]], [3.67, 3.67],
                 color="#9a6aae", linewidth=1.1, linestyle="--")
    full_ax.add_patch(FancyArrowPatch((5.15, 3.67), (5.84, 3.67), arrowstyle="-|>",
                                      mutation_scale=9, color="#9a6aae", linewidth=1.1))
    friction(full_ax, (5.84, 3.67), (6.72, 3.67), r"$F_{f,d}$", 6.5,
             tag="A/A2, B/B2, C/C2")
    full_ax.text(4.02, 3.45, r"$T_{mb},T_{h1},T_{h2},T_{brg}$ aggregate; executed in every friction case",
                 ha="center", fontsize=6.7, color="#6a4c93")

    # Axial masses u_b, u_e, and u_f are dropped; their retained load-path springs stay orange.
    wall(full_ax, 4.03, 2.05, 0.86)
    node(full_ax, 5.15, 2.05, r"$u_b$", 6, dropped_color,
         rf"q6; $m={parameters['m_b']:.3f}$ kg")
    node(full_ax, 6.55, 2.05, r"$u_e$", 7, dropped_color,
         rf"q7; $m={parameters['m_e']:.3f}$ kg")
    node(full_ax, 6.55, 0.63, r"$u_f$", 8, dropped_color,
         rf"q8; $m={parameters['m_f']:.3f}$ kg")
    node(full_ax, 8.35, 2.05, r"$u_n$", 9, stage_color,
         rf"q9; $m={parameters['m_n']:.3f}$ kg")
    node(full_ax, 10.55, 2.05, r"$x_s$", 10, stage_color,
         rf"q10; $m={parameters['m_stage']:.3f}$ kg", width=1.04)
    spring(full_ax, (4.05, 2.05), (4.67, 2.05), amplitude=0.055)
    spring(full_ax, (5.63, 2.05), (6.07, 2.05), amplitude=0.052)
    spring(full_ax, (6.55, 1.72), (6.55, 0.96), color=discarded_spring_color,
           amplitude=0.055, linestyle="--")
    full_ax.text(6.82, 1.30, r"$k_{shb}$", fontsize=7.2, color="#666666")
    full_ax.text(4.36, 2.38, rf"$k_{{brg}}$ {shares['brg']:.0f}%", ha="center",
                 fontsize=7.1, color="#9a5600")
    full_ax.text(5.85, 2.38, rf"$k_{{sha}}(x_s)$ {shares['sha']:.0f}%", ha="center",
                 fontsize=7.1, color="#9a5600")

    # Nut-column transformer. Rotation drops vertically; ball stiffness and friction share two nodes.
    transformer = FancyBboxPatch((7.97, 2.66), 0.76, 0.66, boxstyle="round,pad=0.04",
                                 facecolor="#fff2dc", edgecolor=spring_color,
                                 linewidth=1.4, zorder=4)
    full_ax.add_patch(transformer)
    full_ax.text(8.35, 3.08, "TF", ha="center", va="center", fontsize=8.8,
                 fontweight="bold", color="#9a5600", zorder=5)
    full_ax.text(8.35, 2.87, r"$r$: $u_t=u_e+r\theta_{s2}$", ha="center", va="center",
                 fontsize=6.1, color="#9a5600", zorder=5)
    full_ax.add_patch(FancyArrowPatch((8.35, 4.97), (8.35, 3.34), arrowstyle="-|>",
                                      mutation_scale=9, color=rigid_color, linewidth=1.4))
    full_ax.plot([7.03, 7.45], [2.05, 2.05], color=rigid_color, linewidth=1.4)
    full_ax.plot([7.97, 7.45], [2.78, 2.18], color=rigid_color, linewidth=1.3)
    spring(full_ax, (7.45, 2.18), (7.87, 2.18), amplitude=0.043)
    friction(full_ax, (7.45, 1.73), (7.87, 1.73), r"$F_{f,n}$", 6.1,
             tag="B/B2, C/C2")
    full_ax.plot([7.45, 7.45], [1.73, 2.18], color=rigid_color, linewidth=1.0)
    full_ax.plot([7.87, 7.87], [1.73, 2.18], color=rigid_color, linewidth=1.0)
    full_ax.text(7.66, 2.42, rf"$k_{{ball}}$ {shares['ball']:.0f}%", ha="center",
                 fontsize=7.0, color="#9a5600")
    friction(full_ax, (8.92, 3.18), (9.78, 3.18), r"$F_{f,r}$", 6.1,
             tag="B/B2, C/C2")
    full_ax.text(9.35, 3.49, r"gross nut rolling drag, $v_r=\dot x_d$",
                 ha="center", fontsize=6.5, color="#8d2936")
    spring(full_ax, (8.83, 2.05), (10.03, 2.05), amplitude=0.062)
    full_ax.text(9.43, 2.38, rf"$k_{{mnt}}$ {shares['mnt']:.0f}%", ha="center",
                 fontsize=7.1, color="#9a5600")

    # Guideway port remains stage referenced and carries case tags.
    full_ax.plot([11.06, 11.55], [2.05, 2.05], color=rigid_color, linewidth=1.5)
    full_ax.plot([11.15, 11.95], [1.53, 1.53], color=rigid_color, linewidth=1.8)
    for pad_x in (11.27, 11.47, 11.67, 11.87):
        full_ax.add_patch(Rectangle((pad_x - 0.065, 1.34), 0.13, 0.13,
                                    facecolor="#f8dce1", edgecolor=friction_color, linewidth=0.8))
    friction(full_ax, (11.55, 1.34), (11.55, 0.54), r"$F_{f,g}$", 6.4,
             tag="A/A2, C/C2")
    ground(full_ax, 11.55, 0.40, 0.58)

    # Right-hand evidence column removes numeric labels from the mechanism drawing.
    evidence = FancyBboxPatch((12.20, 4.22), 3.08, 2.68, boxstyle="round,pad=0.08",
                              facecolor="#f7fafc", edgecolor="#8da2b2", linewidth=1.0)
    full_ax.add_patch(evidence)
    full_ax.text(13.74, 6.66, "reduction evidence", ha="center", fontsize=8.6,
                 fontweight="bold", color="#425b6b")
    evidence_lines = (
        rf"$K_m={constants['K_m']:.2e}$ N/m",
        rf"drive pole {reduced_modes[0]:.0f} Hz; 140 to 190 Hz step-dependent band",
        rf"$k_c/r^2={coupling_reflected:.2e}$ N/m, discarded compliance",
        rf"$k_\theta/r^2={torsion_reflected:.2e}$ N/m per link, discarded",
        rf"$k_{{ax}}={MODEL['k_ax']:.2e}$ N/m, relative mode {reduced_modes[1]:.0f} Hz",
        rf"stub-rich modes {torsional_stub_mode/1000:.2f}/{axial_stub_mode/1000:.2f} kHz, above band",
    )
    full_ax.text(12.40, 6.35, "\n".join(evidence_lines), ha="left", va="top",
                 fontsize=7.1, linespacing=1.42, color="#425b6b")

    full_ax.text(12.22, 3.95, r"sign key: $+x$ right; $+\theta$ by RH rule about $+x$; $r=L/(2\pi)$",
                 fontsize=6.9, color="#176a55")
    full_ax.text(12.22, 3.68, r"$\delta_n=u_n-u_e-r\theta_{s2}$; $v_n=-\dot\delta_n$",
                 fontsize=7.6, color="#9a5600", fontweight="bold")

    # Separate legends for mass destination and spring retention.
    full_ax.text(12.22, 3.36, "mass fill", fontsize=7.2, fontweight="bold", color="#4d555c")
    for y_pos, color, label in ((3.10, drive_color, "migrates to $m_d$"),
                                (2.84, stage_color, "migrates to $m_s$"),
                                (2.58, dropped_color, "mass dropped")):
        full_ax.add_patch(Rectangle((12.22, y_pos - 0.09), 0.22, 0.18,
                                    facecolor=color, edgecolor="#777777", linewidth=0.7))
        full_ax.text(12.53, y_pos, label, va="center", fontsize=6.7, color="#4d555c")
    full_ax.text(14.05, 3.36, "spring stroke", fontsize=7.2, fontweight="bold", color="#4d555c")
    spring(full_ax, (14.05, 3.08), (14.55, 3.08), color=spring_color, amplitude=0.04)
    full_ax.text(14.63, 3.08, "retained", va="center", fontsize=6.7, color="#4d555c")
    spring(full_ax, (14.05, 2.75), (14.55, 2.75), color=discarded_spring_color,
           amplitude=0.04, linestyle="--")
    full_ax.text(14.63, 2.75, "discarded", va="center", fontsize=6.7, color="#4d555c")

    exclusion = FancyBboxPatch((12.20, 0.28), 3.08, 1.98, boxstyle="round,pad=0.08",
                               facecolor="#fbfbfb", edgecolor="#8a8a8a", linewidth=1.0,
                               linestyle="--")
    full_ax.add_patch(exclusion)
    full_ax.text(13.74, 2.02, "deliberately absent", ha="center", fontsize=8.0,
                 fontweight="bold", color="#555555")
    full_ax.text(12.43, 1.69, "base compliance\npayload bracket\npitch, yaw, roll\nscrew bending",
                 ha="left", va="top", fontsize=7.0, linespacing=1.35, color="#555555")

    # Two-DOF panel with distinct structural, damping, and friction elements.
    reduced_ax.set_xlim(0.0, 10.0)
    reduced_ax.set_ylim(0.0, 4.35)
    reduced_ax.text(5.0, 4.06, "Retained two-DOF model: executable through the axial mode",
                    ha="center", fontsize=12.7, fontweight="bold")
    wall(reduced_ax, 0.38, 2.30, 0.82)
    reduced_ax.text(0.38, 2.87, r"$x_{cmd}$", ha="center", fontsize=8.2)
    node(reduced_ax, 2.35, 2.30, r"$x_d$", 1, drive_color,
         rf"$m_d={constants['m_d']:.2f}$ kg", width=1.18)
    node(reduced_ax, 7.70, 2.30, r"$x_s$", 2, stage_color,
         rf"$m_s={MODEL['m_s']:.2f}$ kg", width=1.12)
    spring(reduced_ax, (0.40, 2.30), (1.76, 2.30), color="#c08a00")
    reduced_ax.text(1.15, 2.66, r"$K_m(x_{cmd}-x_d)$", ha="center", fontsize=7.7)
    reduced_ax.text(1.15, 2.47, rf"drive pole {reduced_modes[0]:.0f} Hz", ha="center",
                    fontsize=6.5, color="#8a6200")
    spring(reduced_ax, (2.95, 2.60), (7.14, 2.60), amplitude=0.08)
    dashpot(reduced_ax, (2.95, 2.05), (7.14, 2.05), width=0.10)
    friction(reduced_ax, (2.95, 1.60), (7.14, 1.60), r"$F_{f,n}$", 7.0,
             tag="B/B2, C/C2")
    reduced_ax.text(5.03, 2.87, rf"$k_{{ax}}$ structural compliance; relative mode {reduced_modes[1]:.0f} Hz",
                    ha="center", fontsize=7.7,
                    color="#9a5600")
    reduced_ax.text(5.03, 1.83, r"$c_{ax}$", ha="center", fontsize=7.4, color=damping_color)
    reduced_ax.text(5.03, 1.18, r"internal equal-opposite port, $v_n=\dot x_d-\dot x_s$",
                    ha="center", fontsize=7.4, color="#8d2936")
    spring(reduced_ax, (2.05, 1.97), (2.05, 0.66), color="#c08a00", amplitude=0.06)
    ground(reduced_ax, 2.05, 0.52, 0.44)
    reduced_ax.text(1.15, 1.20, r"$K_{det}(x_d)$", fontsize=7.0, color="#8a6200")
    dashpot(reduced_ax, (2.35, 1.97), (2.35, 0.66), width=0.08)
    ground(reduced_ax, 2.35, 0.52, 0.44)
    friction(reduced_ax, (2.68, 1.97), (2.68, 0.66), r"$F_{f,d}$", 6.2,
             tag="A/B/C")
    ground(reduced_ax, 2.68, 0.52, 0.44)
    friction(reduced_ax, (3.18, 1.40), (3.18, 0.66), r"$F_{f,r}$", 6.2,
             tag="B/B2, C/C2")
    ground(reduced_ax, 3.18, 0.52, 0.44)
    friction(reduced_ax, (7.70, 1.97), (7.70, 0.66), r"$F_{f,g}$", 6.2,
             tag="A/A2, C/C2")
    ground(reduced_ax, 7.70, 0.52, 0.50)
    reduced_ax.text(5.0, 0.08, "F_f,d is active in every friction case. F_f,r carries gross nut rolling drag. F_f,n remains internal microslip.",
                    ha="center", fontsize=7.6, color="#555555")

    # One-DOF panel makes the identifiability loss explicit.
    rejected_ax.set_xlim(0.0, 4.2)
    rejected_ax.set_ylim(0.0, 4.35)
    rejected_ax.text(2.1, 4.06, "One DOF: rejected",
                     ha="center", fontsize=12.7, fontweight="bold", color="#9b2f3d")
    rejected_ax.add_patch(FancyBboxPatch((0.12, 0.22), 3.96, 3.45, boxstyle="round,pad=0.08",
                                        facecolor="#fffafa", edgecolor="#b23a48",
                                        linewidth=1.3, linestyle="--"))
    wall(rejected_ax, 0.42, 2.25, 0.76)
    node(rejected_ax, 2.10, 2.25, r"$x,\ m_d+m_s$", 1, "#eef1f3", width=1.20)
    spring(rejected_ax, (0.44, 2.25), (1.48, 2.25), color="#c08a00", amplitude=0.06)
    friction(rejected_ax, (2.10, 1.92), (2.10, 1.25), r"merged $F_f$", 6.0)
    ground(rejected_ax, 2.10, 1.10, 0.50)
    rejected_ax.text(2.10, 0.88, r"$v_n=\dot x_d-\dot x_s=0$", ha="center",
                     fontsize=8.4, color="#9b2f3d", fontweight="bold")
    rejected_ax.plot([0.96, 3.24], [0.38, 0.68], color="#b23a48", linewidth=2.0)
    rejected_ax.plot([0.96, 3.24], [0.68, 0.38], color="#b23a48", linewidth=2.0)
    rejected_ax.text(2.10, 3.42, "loses the 696 Hz relative mode\nmerges nut and guideway sites\nremoves the nut-friction port",
                     ha="center", va="top", fontsize=7.7, linespacing=1.35, color="#6a3038")

    fig.subplots_adjust(left=0.025, right=0.985, top=0.975, bottom=0.035)
    output = ASSET_DIR / "kinematic_diagram.svg"
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_reduced_bond_graph() -> Path:
    """Render the reduced-model bond graph and friction incidence audit."""
    fig, ax = plt.subplots(figsize=(13.2, 6.3))
    ax.axis("off")
    ax.set_xlim(0.0, 13.2)
    ax.set_ylim(0.0, 6.3)
    ax.text(6.6, 6.02, "Reduced-model bond graph and power-port audit",
            ha="center", fontsize=15, fontweight="bold")

    def junction(x: float, y: float, label: str, color: str) -> None:
        circle = Circle((x, y), 0.38, facecolor=color, edgecolor="#39434d", linewidth=1.4, zorder=4)
        ax.add_patch(circle)
        ax.text(x, y + 0.07, "1", ha="center", va="center", fontsize=11,
                fontweight="bold", zorder=5)
        ax.text(x, y - 0.16, label, ha="center", va="center", fontsize=7.2, zorder=5)

    def element(x: float, y: float, label: str, color: str = "#f5f6f7", width: float = 1.28) -> None:
        box = FancyBboxPatch((x - width / 2.0, y - 0.27), width, 0.54,
                             boxstyle="round,pad=0.04", facecolor=color,
                             edgecolor="#59636d", linewidth=1.0, zorder=3)
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center", fontsize=8.2, zorder=4)

    def bond(start: tuple[float, float], end: tuple[float, float], color: str = "#59636d",
             arrow_fraction: float = 0.67) -> None:
        ax.plot([start[0], end[0]], [start[1], end[1]], color=color, linewidth=1.8, zorder=1)
        x = start[0] + arrow_fraction * (end[0] - start[0])
        y = start[1] + arrow_fraction * (end[1] - start[1])
        dx, dy = end[0] - start[0], end[1] - start[1]
        scale = max(np.hypot(dx, dy), 1e-9)
        ax.add_patch(FancyArrowPatch((x - 0.12 * dx / scale, y - 0.12 * dy / scale),
                                     (x + 0.12 * dx / scale, y + 0.12 * dy / scale),
                                     arrowstyle="-|>", mutation_scale=10, color=color,
                                     linewidth=1.2, zorder=2))

    drive_x, center_x, stage_x, graph_y = 3.25, 6.60, 9.95, 3.65
    junction(drive_x, graph_y, r"$\dot x_d$", "#dceef6")
    junction(stage_x, graph_y, r"$\dot x_s$", "#dff2ea")
    zero = Circle((center_x, graph_y), 0.36, facecolor="#fde8ca", edgecolor="#d97800",
                  linewidth=1.4, zorder=4)
    ax.add_patch(zero)
    ax.text(center_x, graph_y + 0.07, "0", ha="center", va="center", fontsize=11,
            fontweight="bold", zorder=5)
    ax.text(center_x, graph_y - 0.16, r"$F_{ax}$", ha="center", va="center", fontsize=7.2,
            zorder=5)
    bond((3.63, graph_y), (6.24, graph_y), "#d97800")
    bond((9.57, graph_y), (6.96, graph_y), "#d97800")

    element(1.05, graph_y, r"Se: $F_{mag}+F_{det}$", "#fff2c7", 1.75)
    bond((1.93, graph_y), (2.87, graph_y), "#c08a00")
    element(drive_x, 5.05, r"I: $m_d$", "#dceef6")
    bond((drive_x, 4.78), (drive_x, 4.03), "#277da1")
    element(stage_x, 5.05, r"I: $m_s$", "#dff2ea")
    bond((stage_x, 4.78), (stage_x, 4.03), "#218c74")

    element(5.02, 2.25, r"C: $1/k_{ax}$", "#fde8ca")
    element(6.60, 2.25, r"R: $c_{ax}$", "#e7e0ef")
    element(8.18, 2.25, r"R: $F_{f,n}$", "#f8dce1")
    bond((center_x, 3.29), (5.02, 2.52), "#d97800")
    bond((center_x, 3.29), (6.60, 2.52), "#6a4c93")
    bond((center_x, 3.29), (8.18, 2.52), "#b23a48")

    element(2.00, 2.25, r"R: $c_m,F_{f,d},F_{f,r}$", "#eee8f3", 2.08)
    bond((3.03, 3.38), (2.36, 2.52), "#6a4c93")
    element(11.18, 2.25, r"R: $F_{f,g}$", "#f8dce1", 1.45)
    bond((10.17, 3.38), (10.82, 2.52), "#b23a48")

    ax.text(6.60, 4.30, r"internal port: $v_n=\dot x_d-\dot x_s$",
            ha="center", fontsize=8.5, color="#8d2936")
    ax.text(4.70, 3.92, r"$-F_{f,n}$", fontsize=8.0, color="#8d2936")
    ax.text(8.05, 3.92, r"$+F_{f,n}$", fontsize=8.0, color="#8d2936")
    ax.text(6.60, 1.48,
            r"$\mathbf{H}_g=[0,1],\quad \mathbf{H}_n=[1,-1],\quad \mathbf{H}_r=\mathbf{H}_d=[1,0]$",
            ha="center", fontsize=10.0, color="#39434d")
    ax.text(6.60, 0.92,
            r"$\mathbf{Q}_f=-\mathbf{H}^T F_f,\qquad P_f=\dot{\mathbf{x}}^{T}\mathbf{Q}_f=-v_fF_f\leq0$",
            ha="center", fontsize=10.0, color="#39434d")
    ax.text(6.60, 0.35,
            "Each bond carries force and velocity. The paired nut bonds show equal-opposite generalized force and preserve power.",
            ha="center", fontsize=8.4, color="#555555")
    fig.tight_layout()
    output = ASSET_DIR / "reduced_bond_graph.svg"
    save_svg(fig, output)
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
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_position_dependence() -> Path:
    positions = np.array([50.0, 150.0, 250.0])
    k_sha = np.array([2.0e8, 6.7e7, 4.0e7])
    k_ax = np.array([1.29e7, 1.14e7, 1.02e7])
    mode = np.array([739.9, 695.6, 657.9])
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
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_stepper_resonance_visibility() -> Path:
    """Show how damping and output selection hide or expose the low motor mode."""
    constants = physical_constants()
    frequencies = np.logspace(np.log10(80.0), np.log10(400.0), 1400)
    mass, _baseline_damping, stiffness, input_vector = linear_matrices((), "none")
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    low_mode = _linear_modes(mass, stiffness)[0]
    damping_ratios = (0.02, MODEL["zeta_m"], 0.10, 0.50)
    colors = ("#b23a48", "#277da1", "#d97800", "#7a7a7a")
    fig, (stage_ax, drive_ax) = plt.subplots(1, 2, figsize=(11.4, 4.7), sharex=True)
    for zeta, color in zip(damping_ratios, colors):
        c_m = 2.0 * zeta * np.sqrt(constants["K_drive"] * constants["m_d"])
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
             f"Published detent torque {constants['T_det']:.3f} N m is enabled at the stable zero-phase equilibrium.",
             ha="center", fontsize=8.5, color="#555555")
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.93))
    output = ASSET_DIR / "stepper_resonance_visibility.svg"
    save_svg(fig, output)
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
        (drive_response, r"$X_d/X_{cmd}$: command to rotor-equivalent drive", "#b23a48", "-"),
        (stage_response, r"$X_s/X_{cmd}$: command to stage", "#277da1", "-"),
        (rotor_to_stage, r"$X_s/X_d$: rotor-equivalent drive to stage", "#218c74", "--"),
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
    save_svg(fig, output)
    plt.close(fig)
    return output


def generated_summary(linear_metrics: dict[str, dict[str, float | np.ndarray]],
                      time_metrics: dict[str, dict[str, float]],
                      verification: dict[str, object]) -> str:
    lines = [
        "<!-- BEGIN GENERATED RESPONSE SUMMARY -->",
        "| Case | Friction law | Presliding modes (Hz) | Presliding tangent gain $X_s/X_{cmd}$ | Smallest first-yield travel | First-step overshoot | Full-sequence RMS deviation | Peak absolute deviation | Final-window RMS deviation |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for key, case in CASES.items():
        modes = linear_metrics[key]["modes"]
        mode_text = f"{modes[0]:.1f}, {modes[1]:.1f}"
        friction_label = {"none": "none", "lugre": "LuGre", "gms": "GMS"}[case["friction"]]
        first_yield = float(linear_metrics[key]["first_yield_m"])
        yield_text = "not applicable" if not np.isfinite(first_yield) else f"{first_yield * 1e6:.3f} µm"
        lines.append(
            f"| {key} | {friction_label} | {mode_text} | {linear_metrics[key]['tangent_dc_gain']:.5f} | "
            f"{yield_text} | "
            f"{time_metrics[key]['first_overshoot_pct']:.1f}% | "
            f"{time_metrics[key]['rms_sequence_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['max_abs_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['rms_final_error_nm']:.1f} nm |"
        )
    lines.extend([
        "",
        "The tangent gain is a local presliding linearization. It is valid only below the listed first-yield travel and is not a full-range tracking gain. Sustained travel produces bounded friction offsets in the nonlinear model. "
        "The three deviation columns use $d(t)=x_{cmd}(t)-x_s(t)$. They describe the open-loop modeled plant response under each friction law, not closed-loop servo tracking performance. The final column summarizes the last 2 ms of the nonlinear run and is not an identified settling specification. "
        "All cases include rated-current commutation, enabled detent torque, and the highlighted electromagnetic damping assumption. Case 0 remains frictionless.",
        "",
        "### Generated reduction audit",
        "",
        "| Quantity | Executed value |",
        "|---|---:|",
        f"| Closure-derived $k_{{ball}}$ | {verification['parameters']['k_ball'] / 1e6:.3f} MN/m |",
        f"| Motor rotor inertia | {verification['parameters']['J_m']:.3e} kg m² |",
        f"| Coupling inertia | {verification['parameters']['J_c']:.3e} kg m² |",
        f"| 0.320 m screw inertia | {verification['parameters']['screw_inertia']:.3e} kg m² |",
        f"| 0.320 m screw mass | {verification['parameters']['screw_mass']:.4f} kg |",
        f"| Full-model reflected drivetrain mass | {verification['parameters']['m_d_reflected']:.3f} kg |",
        f"| Rated-current holding torque | {MODEL['T_max']:.3f} N m |",
        f"| Enabled detent torque | {MODEL['T_det']:.3f} N m |",
        f"| Full/reduced sequence RMS residual | {verification['rms_residual_nm']:.3f} nm |",
        f"| Full/reduced sequence peak residual | {verification['peak_residual_nm']:.3f} nm |",
        "",
        "The reduced drive mass is derived from the listed component inertias and the current lead. It is not an independent input.",
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


def generated_bode_comparison(frequencies: np.ndarray,
                              responses: dict[str, np.ndarray]) -> str:
    magnitudes = {
        key: 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))
        for key, response in responses.items()
    }
    peak_mask = (frequencies >= 620.0) & (frequencies <= 830.0)
    peak_indices = np.flatnonzero(peak_mask)

    def peak(key: str) -> tuple[float, float]:
        index = peak_indices[np.argmax(magnitudes[key][peak_mask])]
        return float(frequencies[index]), float(magnitudes[key][index])

    baseline_frequency, _ = peak("0")
    lines = [
        "<!-- BEGIN GENERATED BODE COMPARISON -->",
        "| Topology | Local peak | Shift from Case 0 | Largest GMS/LuGre gap | Cause |",
        "|---|---:|---:|---:|---|",
        f"| Case 0 | {baseline_frequency:.1f} Hz | reference | not applicable | No friction tangent |",
    ]
    rows = (
        ("A/A2", "A", "A2", "Guideway presliding stiffness acts against ground"),
        ("B/B2", "B", "B2", "Nut microslip shifts the relative mode; rolling and drivetrain tangents act on the drive"),
        ("C/C2", "C", "C2", "All four friction tangents are active"),
    )
    for label, lugre_key, gms_key, cause in rows:
        peak_frequency, _ = peak(gms_key)
        shift = peak_frequency - baseline_frequency
        shift_percent = 100.0 * shift / baseline_frequency
        difference = magnitudes[gms_key] - magnitudes[lugre_key]
        maximum_index = int(np.argmax(np.abs(difference)))
        maximum = abs(float(difference[maximum_index]))
        maximum_frequency = float(frequencies[maximum_index])
        lines.append(
            f"| {label} | {peak_frequency:.1f} Hz | +{shift:.1f} Hz, +{shift_percent:.1f}% | "
            f"{maximum:.2f} dB at {maximum_frequency:.0f} Hz | {cause} |"
        )
    lines.append("<!-- END GENERATED BODE COMPARISON -->")
    return "\n".join(lines)


def update_generated_bode_comparison(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED BODE COMPARISON -->.*?<!-- END GENERATED BODE COMPARISON -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated Bode comparison markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_presliding_summary(experiment: dict[str, object]) -> str:
    metrics = experiment["metrics"]
    lines = [
        "<!-- BEGIN GENERATED PRESLIDING SUMMARY -->",
        "| Executed metric | LuGre A | GMS A2 | GMS minus LuGre |",
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
        precision = 4 if unit == "N" else 2
        lines.append(
            f"| {label} | {lugre:.{precision}f} {unit} | {gms:.{precision}f} {unit} | "
            f"{gms - lugre:+.{precision}f} {unit} |"
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
        "experiment is intended to distinguish. The provisional parameters do not guarantee that "
        "GMS closes more tightly than LuGre; measured force loops must decide that question.",
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

    def derived_output(match: re.Match[str]) -> str:
        key, value = match.group(1), match.group(2)
        escaped_key = html.escape(key, quote=True)
        escaped_value = html.escape(value.strip())
        return keep(
            f'<output class="derived-output" data-derived="{escaped_key}" '
            f'aria-label="Derived value {escaped_key}">{escaped_value}</output>'
        )

    text = re.sub(r"\[\[derived:([A-Za-z0-9_]+)=([^\]]+)\]\]", derived_output, text)
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
            heading_class = "appendix-heading" if level == 2 and title.startswith("Appendix") else ""
            class_attribute = f' class="{heading_class}"' if heading_class else ""
            output.append(
                f'<h{level} id="{section_id}"{class_attribute}>{render_inline(title)}</h{level}>'
            )
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
h1,h2,h3,h4 {{ line-height:1.24; scroll-margin-top:5rem; }} h1 {{ font-size:clamp(2rem,4vw,3rem); margin-top:0; }}
h2 {{ margin:4.5rem -.8rem 1.5rem; padding:1rem 1.1rem; border:1px solid var(--line); border-left:6px solid var(--accent); border-radius:9px; background:var(--soft); box-shadow:0 8px 22px rgba(22,36,46,.06); }}
h2.appendix-heading {{ border-left-color:var(--assumed-line); background:color-mix(in srgb,var(--assumed) 28%,var(--card)); }}
h3 {{ margin-top:2.2rem; color:var(--accent); }}
p {{ margin:.85rem 0; }} a {{ color:var(--accent); }} strong {{ color:var(--text); }} hr {{ border:0; border-top:1px solid var(--line); margin:2rem 0; }}
blockquote {{ margin:1.2rem 0; padding:.75rem 1rem; background:var(--soft); border-left:4px solid var(--accent); border-radius:0 8px 8px 0; color:var(--muted); }}
.table-wrap {{ overflow-x:auto; margin:1.2rem 0; }} table {{ width:100%; border-collapse:collapse; font-size:.92rem; }} th,td {{ border:1px solid var(--line); padding:.55rem .65rem; vertical-align:top; }} th {{ background:var(--soft); }}
.parameter-input {{ width:100%; min-width:7rem; padding:.38rem .48rem; color:var(--text); background:var(--card); border:1px solid var(--line); border-radius:5px; font:inherit; font-variant-numeric:tabular-nums; }}
.parameter-input:focus {{ outline:2px solid var(--accent); outline-offset:1px; }}
.assumed-input {{ background:var(--assumed); border-color:var(--assumed-line); font-weight:700; }}
.derived-output {{ display:inline-block; width:100%; min-width:7rem; padding:.38rem .48rem; color:var(--accent); background:var(--soft); border:1px dashed var(--accent); border-radius:5px; font-variant-numeric:tabular-nums; font-weight:700; }}
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
.parameter-group {{ margin:1.15rem 0; border-left:5px solid var(--accent); background:var(--soft); }} .parameter-group summary {{ font-size:1.04rem; }}
pre {{ overflow:auto; background:var(--code); color:#e8edf2; border-radius:9px; padding:1rem; font-size:.87rem; }} code {{ font-family:Cascadia Code,Consolas,monospace; }} p code,li code,td code {{ background:var(--soft); border:1px solid var(--line); border-radius:4px; padding:.1rem .28rem; }}
.display-math {{ overflow-x:auto; padding:.5rem 0; }} img {{ display:block; max-width:100%; height:auto; margin:1.3rem auto; border-radius:6px; }}
.footer {{ color:var(--muted); font-size:.78rem; margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line); }}
@media (max-width:920px) {{ .layout {{ grid-template-columns:1fr; padding:.7rem; }} nav {{ position:relative; top:auto; max-height:18rem; }} article {{ padding:1.2rem; }} .hide-small {{ display:none; }} .live-plot-grid {{ grid-template-columns:1fr; }} }}
@media print {{ .topbar,nav {{ display:none; }} body {{ background:white; }} .layout {{ display:block; padding:0; }} article {{ max-width:none; border:0; box-shadow:none; }} details {{ break-inside:avoid; }} }}
</style>
</head>
<body>
<div class="topbar"><span class="name">{html.escape(title)}</span><button onclick="setDetails(true)">Expand derivations</button><button onclick="setDetails(false)">Collapse</button><button onclick="saveParameterInputs()">Save variables</button><button onclick="saveEditedHtml()">Save HTML copy</button><button onclick="resetParameterInputs()">Reset inputs</button><button class="hide-small" onclick="toggleTheme()">Theme</button><button class="hide-small" onclick="window.print()">Print</button></div>
<div class="layout"><nav><div class="caption">On this page</div>{''.join(toc_html)}</nav><article><div class="edit-note"><span class="assumed-swatch"></span>Amber inputs are unidentified assumptions. Values auto-save to browser storage and the page URL; “Save HTML copy” embeds them in a chosen HTML file. Derived outputs and the live transfer panel recalculate immediately. Publication SVGs require a Python build.<span id="parameter-save-status" class="save-status">Loading values…</span></div>{body}<div class="footer">Rendered from {html.escape(markdown_path.name)} · {generated}</div></article></div>
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
    else setParameterStatus('Browser storage unavailable; use Save HTML copy', 'warn');
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
  const lead = parameterNumber('lead', 1.0e-3);
  const teeth = parameterNumber('rotor_teeth', 50);
  const jm = parameterNumber('J_m', 9.0e-7);
  const jc = parameterNumber('J_c', 1.18e-6);
  const screwLength = parameterNumber('screw_length', 0.320);
  const screwDiameter = parameterNumber('screw_diameter', 8.0e-3);
  const screwDensity = parameterNumber('screw_density', 7850.0);
  const tmax = parameterNumber('holding_torque', 0.060);
  const tdet = parameterNumber('detent_torque', 0.005);
  const detentPhase = parameterNumber('detent_phase', 0.0);
  const couplingSeries = parameterNumber('k_c_series', 68.7549);
  const ms = parameterNumber('reduced_stage_mass', 0.60);
  const kax = parameterNumber('reduced_axial_stiffness', 1.14e7);
  const cax = parameterNumber('axial_damping', 55.0);
  const zeta = parameterNumber('electromagnetic_zeta', 0.05);
  const microstepDivisor = parameterNumber('microstep_divisor', 64);
  if (!(lead>0 && teeth>0 && jm>0 && jc>=0 && screwLength>0 && screwDiameter>0 &&
        screwDensity>0 && tmax>0 && tdet>=0 && ms>0 && kax>0 && cax>=0 && zeta>=0 && microstepDivisor>=1))
    throw new Error('Geometry, masses, torque, and stiffness must be positive; damping and detent torque cannot be negative.');
  const r = lead/(2*Math.PI);
  const screwRadius = screwDiameter/2;
  const screwMass = screwDensity*Math.PI*screwRadius*screwRadius*screwLength;
  const screwInertia = 0.5*screwMass*screwRadius*screwRadius;
  const jTotal = jm+jc+screwInertia;
  const md = jTotal/(r*r);
  const km = teeth*tmax/(r*r);
  const kdet = 4*teeth*tdet*Math.cos(detentPhase)/(r*r);
  const kdrive = km+kdet;
  if (!(kdrive>0)) throw new Error('The selected detent phase makes the net drive tangent non-positive.');
  const cm = 2*zeta*Math.sqrt(kdrive*md);
  const frequencies=[], drive=[], stage=[], rotorStage=[];
  const count=560, logMin=Math.log10(100), logMax=Math.log10(3000);
  for (let i=0; i<count; i++) {{
    const frequency = Math.pow(10, logMin + (logMax-logMin)*i/(count-1));
    const omega = 2*Math.PI*frequency;
    const a = {{re:kdrive+kax-md*omega*omega, im:omega*(cm+cax)}};
    const b = {{re:-kax, im:-omega*cax}};
    const d = {{re:kax-ms*omega*omega, im:omega*cax}};
    const determinant = cSub(cMul(a,d), cMul(b,b));
    const gd = cDiv(cMul({{re:km,im:0}},d), determinant);
    const gs = cDiv(cMul({{re:-km,im:0}},b), determinant);
    frequencies.push(frequency); drive.push(gd); stage.push(gs); rotorStage.push(cDiv(gs,gd));
  }}
  const qa=md*ms, qb=md*kax+ms*(kdrive+kax), qc=kdrive*kax;
  const discriminant=Math.max(qb*qb-4*qa*qc,0);
  const roots=[(qb-Math.sqrt(discriminant))/(2*qa),(qb+Math.sqrt(discriminant))/(2*qa)];
  const modes=roots.map(value => Math.sqrt(Math.max(value,0))/(2*Math.PI));
  return {{
    frequencies, drive, stage, rotorStage, modes, md, ms, km, kdet, kdrive,
    kax, cax, zeta, lead, teeth, r, jm, jc, screwLength, screwDiameter,
    screwDensity, screwMass, screwInertia, jTotal, tmax, tdet, detentPhase,
    couplingSeries, couplingHalf:2*couplingSeries, kappa:teeth/r,
    fullStep:lead/(4*teeth), quarterStep:lead/(16*teeth),
    commandStep:lead/(4*teeth*microstepDivisor), interpolatedStep:lead/(4*teeth*256),
    microstepDivisor, fmax:tmax/r, cm
  }};
}}
function formatDerivedValue(key, value) {{
  const scientific = new Set([
    'transmission_ratio','magnetic_stiffness','detent_stiffness','screw_inertia',
    'screw_segment_inertia','full_step_pitch','quarter_step_bound',
    'command_step','interpolated_step'
  ]);
  if (scientific.has(key)) return value.toExponential(5);
  if (key==='reduced_drive_mass') return value.toFixed(3);
  if (key==='screw_mass' || key==='screw_segment_mass') return value.toFixed(6);
  if (key==='k_c_half') return value.toFixed(3);
  if (key.endsWith('_hz')) return value.toFixed(2);
  return Number(value).toPrecision(6);
}}
function refreshDerivedOutputs(data) {{
  const values = {{
    transmission_ratio:data.r,
    reduced_drive_mass:data.md,
    magnetic_stiffness:data.km,
    detent_stiffness:data.kdet,
    full_step_pitch:data.fullStep,
    quarter_step_bound:data.quarterStep,
    command_step:data.commandStep,
    interpolated_step:data.interpolatedStep,
    screw_inertia:data.screwInertia,
    screw_segment_inertia:data.screwInertia/3,
    screw_mass:data.screwMass,
    screw_segment_mass:data.screwMass/3,
    k_c_half:data.couplingHalf,
    mode_1_hz:data.modes[0],
    mode_2_hz:data.modes[1],
    drive_stiffness:data.kdrive,
    force_limit:data.fmax,
    spatial_wavenumber:data.kappa
  }};
  document.querySelectorAll('[data-derived]').forEach(output => {{
    const key=output.dataset.derived;
    if (Object.prototype.hasOwnProperty.call(values,key)) output.textContent=formatDerivedValue(key,values[key]);
  }});
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
  let data;
  try {{ data=liveTransferData(); refreshDerivedOutputs(data); }}
  catch (error) {{
    const summary=document.getElementById('live-model-summary');
    if (summary) summary.textContent='Live calculation error: '+error.message;
    return;
  }}
  const panel=document.querySelector('[data-live-transfer-plots]'); if (!panel) return;
  if (!panel.dataset.initialized) {{
    panel.innerHTML='<div id="live-model-summary" class="live-summary"></div><div class="live-plot-grid"><div class="live-plot-card"><h4>Live magnitude</h4><svg id="live-bode-magnitude" viewBox="0 0 760 360" role="img" aria-label="Live transfer-function magnitude"></svg></div><div class="live-plot-card"><h4>Live phase</h4><svg id="live-bode-phase" viewBox="0 0 760 360" role="img" aria-label="Live transfer-function phase"></svg></div></div>';
    panel.dataset.initialized='true';
  }}
  document.getElementById('live-model-summary').textContent='Live derived values: r='+data.r.toExponential(5)+' m/rad, md='+data.md.toFixed(3)+' kg, Km='+data.km.toExponential(5)+' N/m, Kdet='+data.kdet.toExponential(5)+' N/m · modes '+data.modes.map(value=>value.toFixed(2)+' Hz').join(', ');
  drawLiveBode('live-bode-magnitude',data,false); drawLiveBode('live-bode-phase',data,true);
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
    validate_case_topology()
    constants = physical_constants()
    frequencies, bode, linear_metrics = frequency_responses()
    times, command, time_data, time_metrics = time_responses(constants)
    convergence = gms_step_halving_convergence(constants, times, time_data)
    presliding = presliding_responses(constants)
    verification = full_reduced_verification(frequencies, constants)
    comparison_path = plot_case_response_overlay(frequencies, bode)
    presliding_path = plot_presliding_memory(presliding)
    diagram_path = plot_kinematic_diagram()
    bond_graph_path = plot_reduced_bond_graph()
    verification_path = plot_full_reduced_verification(frequencies, verification)
    position_path = plot_position_dependence()
    resonance_path = plot_stepper_resonance_visibility()
    rotor_stage_path = plot_rotor_stage_transfer_functions(frequencies)
    for obsolete_name in ("bode_all_cases.svg", "step_tracking_all_cases.svg"):
        obsolete_path = ASSET_DIR / obsolete_name
        if obsolete_path.exists():
            obsolete_path.unlink()
    if not args.skip_summary_update:
        update_generated_bode_comparison(generated_bode_comparison(frequencies, bode))
        update_generated_summary(generated_summary(linear_metrics, time_metrics, verification))
        update_generated_presliding_summary(generated_presliding_summary(presliding))
        update_generated_convergence_summary(generated_convergence_summary(convergence))
    description_html = render_document(DESCRIPTION_MD)
    derivation_html = render_document(DERIVATION_MD)
    print(f"Built {comparison_path.relative_to(ROOT)}")
    print(f"Built {presliding_path.relative_to(ROOT)}")
    print(f"Built {diagram_path.relative_to(ROOT)}")
    print(f"Built {bond_graph_path.relative_to(ROOT)}")
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
              f"presliding tangent gain={linear_metrics[key]['tangent_dc_gain']:.6f}; "
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

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
    # Published detent torque.  It is executed as a periodic nonlinear torque;
    # its tangent is reported separately and is never used as a global spring.
    "T_det": 0.005,
    "detent_phase": 0.0,
    # Measured upper axial mode used to calibrate the reduced axial chain.
    # Stage and nut masses live in FULL; m_s, k_ax, and k_ball are derived.
    "axial_mode_target_hz": 695.82,
    # Provisional: retained structural damping, not identified in the source.
    "c_ax": 55.0,
    # Provisional open-loop drive damping ratio.  Driver mode and tuning are
    # not recorded.  The requested baseline is 10% of critical damping and a
    # sensitivity sweep is retained rather than presenting it as identified.
    "zeta_m": 0.10,
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
    "screw_length": 0.192,
    "usable_screw_travel": 0.170,
    "stage_travel": 0.150,
    "screw_diameter": 8.00e-3,
    "screw_density": 7850.0,
    "m_n": 0.050,
    # Measured stage body mass.  The retained stage-side mass also includes
    # the nut body after the internal nut coordinate is collapsed.
    "m_stage": 0.355,
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
    # Identifiable drive-side lump.  Motor/support-bearing drag and gross nut
    # rolling were formerly two laws on the same H=[1,0] port.  Their force,
    # tangent, and damping budgets are combined here; the aggregate still
    # requires identification from a drive-side measurement.
    "d": {"sigma0": 3.00e6, "sigma1": 9.0, "sigma2": 0.45,
          "F_s": 7.0, "F_c": 5.5, "v_s": 2.3e-4, "delta": 1.0, "C_gms": 5.0e3},
}


# Four GMS stop elements share each site's aggregate sigma0 and Stribeck
# force.  Opposing stiffness/force fractions create distinct yield distances
# and therefore non-local reversal memory while retaining the LuGre aggregate.
GMS_WEIGHTS = np.array([0.10, 0.20, 0.30, 0.40])
GMS_STIFFNESS_FRACTIONS = np.array([0.40, 0.30, 0.20, 0.10])
GMS_N = GMS_WEIGHTS.size
PRODUCTION_DT = 2.5e-5
GMS_CONVERGENCE_DTS = (5.0e-5, 2.5e-5, 1.25e-5)
BODE_FOCUS_MIN_HZ = 100.0
BODE_FOCUS_MAX_HZ = 3000.0
# The main response is quantized at 78.125 nm and deliberately spans the
# provisional first-yield distances.  Adjacent increments remain below one
# quarter of a full step.  The final move is positive and returns to zero.
MAIN_LEVELS = np.array([3, -3, 6, -6, 0, 13, 0, -13, -26, -13,
                        0, 13, 26, 13, 0], dtype=float)
MAIN_START = 0.010


# Nested reversals for the dedicated memory experiment.  Counts use the
# conservative 64-microstep STEP/DIR quantum.  The outer excursion is large
# enough to yield two guideway GMS elements, while the force signal remains the
# primary discriminator.
GUIDEWAY_PRESLIDING_LEVELS = np.array(
    [0, 48, 12, 42, 12, 48, 0, -46, -12, -40, -12, -46, 0],
    dtype=float,
)
# With the stage blocked, the drive coordinate is the nut-port deflection.
# The +/-10-count outer loop crosses the provisional 0.20 and 0.533 um stop
# thresholds but remains below the third threshold at 1.20 um.
NUT_PRESLIDING_LEVELS = np.array(
    [0, 10, 3, 9, 3, 10, 0, -10, -3, -9, -3, -10, 0],
    dtype=float,
)
PRESLIDING_START = 0.005
PRESLIDING_RETURN_PAIRS = ((1, 5), (2, 4), (7, 11), (8, 10))


CASES = OrderedDict([
    ("0", {"label": "Case 0: frictionless", "sites": (), "friction": "none", "color": "#252525", "ls": "--"}),
    ("A", {"label": "Case A: lumped drive drag + guideway / LuGre", "sites": ("d", "g"), "friction": "lugre", "color": "#277da1", "ls": "-"}),
    ("A2", {"label": "Case A2: lumped drive drag + guideway / GMS", "sites": ("d", "g"), "friction": "gms", "color": "#70b7cf", "ls": "--"}),
    ("B", {"label": "Case B: lumped drive drag + nut microslip / LuGre", "sites": ("d", "n"), "friction": "lugre", "color": "#e07a15", "ls": "-"}),
    ("B2", {"label": "Case B2: lumped drive drag + nut microslip / GMS", "sites": ("d", "n"), "friction": "gms", "color": "#f5b35f", "ls": "--"}),
    ("C", {"label": "Case C: all identifiable ports / LuGre", "sites": ("d", "g", "n"), "friction": "lugre", "color": "#218c74", "ls": "-"}),
    ("C2", {"label": "Case C2: all identifiable ports / GMS", "sites": ("d", "g", "n"), "friction": "gms", "color": "#72c9ad", "ls": "--"}),
])

PAIRS = (("A", "A2"), ("B", "B2"), ("C", "C2"))
# Parameter provenance used by the browser registry and the standalone
# dependency flowcharts.  Every derived token emitted into either Markdown
# document has an explicit dependency list.  Modal-calibrated k_ax and
# closure-derived k_ball are derived outputs, not independent inputs.
PARAMETER_REGISTRY: dict[str, dict[str, object]] = {
    "lead": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "rotor_teeth": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "holding_torque": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "detent_torque": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "detent_phase": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "J_m": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "J_c": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "screw_length": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "screw_diameter": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "screw_density": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "stage_mass": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "nut_mass": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "reduced_stage_mass": {
        "category": "derived", "dependencies": ("stage_mass", "nut_mass"),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
    "axial_mode_target_hz": {"category": "input", "dependencies": (), "section": "6-reduction-from-ten-dofs-to-two"},
    "reduced_axial_stiffness": {
        "category": "derived",
        "dependencies": ("reduced_drive_mass", "reduced_stage_mass",
                         "magnetic_stiffness", "axial_mode_target_hz"),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
    "k_c_series": {"category": "input", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_theta_a": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_theta_b": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_brg": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_sha": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_shb": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_ball": {
        "category": "derived",
        "dependencies": ("reduced_axial_stiffness", "k_brg", "k_sha", "k_mnt"),
        "section": "4-full-ten-dof-derivation",
    },
    "k_mnt": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "microstep_divisor": {"category": "assumed", "dependencies": (), "section": "5-stepper-input-nonlinear-law-linearization-and-bound"},
    "transmission_ratio": {
        "category": "derived", "dependencies": ("lead",),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "screw_mass": {
        "category": "derived",
        "dependencies": ("screw_length", "screw_diameter", "screw_density"),
        "section": "2-entry-parameters",
    },
    "screw_inertia": {
        "category": "derived",
        "dependencies": ("screw_mass", "screw_diameter"),
        "section": "2-entry-parameters",
    },
    "screw_segment_inertia": {
        "category": "derived", "dependencies": ("screw_inertia",),
        "section": "2-entry-parameters",
    },
    "screw_segment_mass": {
        "category": "derived", "dependencies": ("screw_mass",),
        "section": "2-entry-parameters",
    },
    "total_rotational_inertia": {
        "category": "derived",
        "dependencies": ("J_m", "J_c", "screw_segment_inertia"),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
    "reduced_drive_mass": {
        "category": "derived",
        "dependencies": ("total_rotational_inertia", "transmission_ratio"),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
    "magnetic_stiffness": {
        "category": "derived",
        "dependencies": ("rotor_teeth", "holding_torque", "transmission_ratio"),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "detent_stiffness": {
        "category": "derived",
        "dependencies": ("rotor_teeth", "detent_torque", "detent_phase", "transmission_ratio"),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "k_c_half": {
        "category": "derived", "dependencies": ("k_c_series",),
        "section": "4-full-ten-dof-derivation",
    },
    "full_step_pitch": {
        "category": "derived", "dependencies": ("lead", "rotor_teeth"),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "quarter_step_bound": {
        "category": "derived", "dependencies": ("full_step_pitch",),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "command_step": {
        "category": "derived",
        "dependencies": ("full_step_pitch", "microstep_divisor"),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "interpolated_step": {
        "category": "derived", "dependencies": ("full_step_pitch",),
        "section": "5-stepper-input-nonlinear-law-linearization-and-bound",
    },
    "mode_1_hz": {
        "category": "output",
        "dependencies": ("reduced_drive_mass", "reduced_stage_mass",
                         "magnetic_stiffness", "reduced_axial_stiffness"),
        "section": "5-1-why-the-150-250-hz-stepper-feature-is-difficult-to-see",
    },
    "mode_2_hz": {
        "category": "output",
        "dependencies": ("reduced_drive_mass", "reduced_stage_mass",
                         "magnetic_stiffness", "reduced_axial_stiffness"),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
}

PARAMETER_CATEGORY_STYLE = {
    "input": {"face": "#dceef6", "edge": "#52778b"},
    "assumed": {"face": "#fff0b8", "edge": "#c28a00"},
    "derived": {"face": "#dff2ea", "edge": "#4b806d"},
    "output": {"face": "#eee7f8", "edge": "#76549b"},
}
H = {
    "g": np.array([0.0, 1.0]),
    "n": np.array([1.0, -1.0]),
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
def validate_parameter_registry() -> None:
    """Fail if a rendered derived token lacks explicit dependency metadata."""
    allowed_categories = set(PARAMETER_CATEGORY_STYLE)
    for key, entry in PARAMETER_REGISTRY.items():
        category = entry.get("category")
        dependencies = entry.get("dependencies")
        section = entry.get("section")
        if category not in allowed_categories:
            raise ValueError(f"Parameter {key} has unsupported category {category!r}")
        if not isinstance(dependencies, tuple):
            raise TypeError(f"Parameter {key} dependencies must be an ordered tuple")
        if not isinstance(section, str) or not section:
            raise ValueError(f"Parameter {key} needs a section anchor")
        missing_dependencies = [name for name in dependencies
                                if name not in PARAMETER_REGISTRY]
        if missing_dependencies:
            raise KeyError(f"Parameter {key} has unregistered dependencies: "
                           + ", ".join(missing_dependencies))

    rendered_derived: set[str] = set()
    pattern = re.compile(r"\[\[derived:([A-Za-z0-9_]+)=")
    for document in (DESCRIPTION_MD, DERIVATION_MD):
        rendered_derived.update(pattern.findall(document.read_text(encoding="utf-8")))
    missing = sorted(rendered_derived - PARAMETER_REGISTRY.keys())
    if missing:
        raise KeyError("Derived Markdown tokens lack dependency metadata: "
                       + ", ".join(missing))
    empty = sorted(
        key for key in rendered_derived
        if not PARAMETER_REGISTRY[key]["dependencies"]
    )
    if empty:
        raise ValueError("Derived Markdown tokens need non-empty dependency lists: "
                         + ", ".join(empty))
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
        "B": {"d", "n"}, "B2": {"d", "n"},
        "C": {"d", "g", "n"}, "C2": {"d", "g", "n"},
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


def modal_calibrated_axial_stiffness(
        m_d: float, m_s: float, k_m: float, target_hz: float) -> float:
    """Invert the two-DOF characteristic equation for the axial stiffness.

    The measured upper mode is the calibration datum.  Consequently k_ax is
    not independent of the retained masses or commutation stiffness.
    """
    eigenvalue = (2.0 * np.pi * target_hz) ** 2
    numerator = eigenvalue * m_s * (k_m - eigenvalue * m_d)
    denominator = k_m - eigenvalue * (m_d + m_s)
    if abs(denominator) <= np.finfo(float).eps * max(abs(k_m), 1.0):
        raise ValueError("Axial-mode calibration is singular for the current inputs")
    stiffness = numerator / denominator
    if not np.isfinite(stiffness) or stiffness <= 0.0:
        raise ValueError("Axial-mode calibration does not yield a positive k_ax")
    return float(stiffness)


def closure_ball_stiffness(k_ax: float, k_brg: float,
                           k_sha: float, k_mnt: float) -> float:
    """Close the four-element axial compliance chain for k_ball."""
    remaining_compliance = (
        1.0 / k_ax - 1.0 / k_brg - 1.0 / k_sha - 1.0 / k_mnt
    )
    if remaining_compliance <= 0.0:
        raise ValueError("Axial compliance budget cannot yield a positive k_ball")
    return float(1.0 / remaining_compliance)


def physical_constants() -> dict[str, float]:
    lead = MODEL["lead"]
    teeth = MODEL["rotor_teeth"]
    r = lead / (2.0 * np.pi)
    kappa = 2.0 * np.pi * teeth / lead
    t_max = MODEL["T_max"]
    t_det = MODEL["T_det"]
    k_m = teeth * t_max / r**2
    k_det_amplitude = 4.0 * teeth * t_det / r**2
    k_det = k_det_amplitude * np.cos(MODEL["detent_phase"])
    f_max = t_max / r
    full_step = lead / (4.0 * teeth)
    component = component_parameters()
    m_d = component["J_total"] / r**2
    m_stage = component["m_stage"]
    m_n = component["m_n"]
    m_s = m_stage + m_n
    axial_mode_target_hz = MODEL["axial_mode_target_hz"]
    k_ax = modal_calibrated_axial_stiffness(
        m_d, m_s, k_m, axial_mode_target_hz)
    k_ball = closure_ball_stiffness(
        k_ax, component["k_brg"], component["k_sha"], component["k_mnt"])
    c_m = 2.0 * MODEL["zeta_m"] * np.sqrt(k_m * m_d)
    minimum_local_tangent = k_m - k_det_amplitude
    if minimum_local_tangent <= 0.0:
        raise ValueError("Detent tangent can exceed the commutation tangent")
    minimum_local_omega = np.sqrt(minimum_local_tangent / m_d)
    settling_time_2pct = 4.0 / (MODEL["zeta_m"] * minimum_local_omega)
    plateau_dwell = max(0.100, settling_time_2pct)
    return {
        "r": r,
        "kappa": kappa,
        "T_max": t_max,
        "T_det": t_det,
        "F_max": f_max,
        "K_m": k_m,
        "K_det": k_det,
        "K_det_amplitude": k_det_amplitude,
        "m_d": m_d,
        "m_stage": m_stage,
        "m_n": m_n,
        "m_s": m_s,
        "axial_mode_target_hz": axial_mode_target_hz,
        "k_ax": k_ax,
        "k_ball": k_ball,
        "c_m": c_m,
        "full_step": full_step,
        "quarter_step": full_step / 4.0,
        "command_step": full_step / MODEL["microstep_divisor"],
        "interpolated_step": full_step / 256.0,
        "detent_period": full_step,
        "settling_time_2pct": settling_time_2pct,
        "plateau_dwell": plateau_dwell,
        "metric_window": min(0.020, 0.20 * plateau_dwell),
    }


def component_parameters() -> dict[str, float]:
    """Derive screw inertia and lumped masses from the 0.192 m component."""
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
    """Return the ten-DOF parameters and close the calibrated axial compliance."""
    p = component_parameters()
    constants = physical_constants()
    p["k_ball"] = constants["k_ball"]
    r = constants["r"]
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
    # Global command-to-position model: commutation is the only origin spring.
    # Detent is periodic and belongs in the nonlinear law or a declared local
    # tangent about one selected equilibrium, never in this global K matrix.
    stiffness[0, 0] += k_m_rot
    # The same damping repair validated in Rev 2, expressed in rotational units.
    damping[0, 0] += 2.0 * MODEL["zeta_m"] * np.sqrt(
        k_m_rot * p["J_total"]
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
        held_command = verification_command_position(
            times[i] + 0.5 * dt, constants["command_step"])
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
    command = np.array([
        verification_command_position(t, constants["command_step"]) for t in times
    ])
    full_stage = full_states[:, 9]
    reduced_stage = reduced_states[:, 1]
    residual = full_stage - reduced_stage
    command_amplitude = float(np.max(np.abs(command)))
    rms_residual = float(np.sqrt(np.mean(residual**2)))
    peak_residual = float(np.max(np.abs(residual)))
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
        "rms_residual_nm": rms_residual * 1e9,
        "peak_residual_nm": peak_residual * 1e9,
        "command_amplitude_nm": command_amplitude * 1e9,
        "rms_residual_pct_command": 100.0 * rms_residual / command_amplitude,
        "peak_residual_pct_command": 100.0 * peak_residual / command_amplitude,
    }


def linear_matrices(sites: tuple[str, ...], friction_model: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return M, C, K, B for the presliding linearization of one case."""
    constants = physical_constants()
    m_d, m_s = constants["m_d"], constants["m_s"]
    k_ax, k_m, c_ax = constants["k_ax"], constants["K_m"], MODEL["c_ax"]
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    mass = np.diag([m_d, m_s])
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
def detent_local_mode_band() -> tuple[float, float]:
    """Return the low-mode limits from the two extreme local detent tangents."""
    constants = physical_constants()
    mass, _damping, global_stiffness, _input = linear_matrices((), "none")
    local_modes = []
    for tangent in (-constants["K_det_amplitude"], constants["K_det_amplitude"]):
        stiffness = global_stiffness.copy()
        stiffness[0, 0] += tangent
        local_modes.append(float(_linear_modes(mass, stiffness)[0]))
    return min(local_modes), max(local_modes)


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


def verification_command_position(t: float, command_step: float) -> float:
    """Short linear sequence used only to normalize the reduction residual."""
    if t < 0.005:
        return 0.0
    if t < 0.025:
        return command_step
    if t < 0.045:
        return 0.0
    if t < 0.065:
        return -command_step
    return 0.0
def command_position(t: float, command_step: float, plateau_dwell: float) -> float:
    """Finite-amplitude, yield-spanning, quarter-step-bounded main sequence."""
    if t < MAIN_START:
        return 0.0
    index = min(int((t - MAIN_START) // plateau_dwell), MAIN_LEVELS.size - 1)
    return float(MAIN_LEVELS[index] * command_step)
def main_duration(constants: dict[str, float]) -> float:
    return MAIN_START + constants["plateau_dwell"] * MAIN_LEVELS.size
def presliding_command_position(t: float, microstep: float,
                                plateau_dwell: float,
                                levels: np.ndarray = GUIDEWAY_PRESLIDING_LEVELS) -> float:
    """Nested back-and-forth reversals quantized to the STEP/DIR input."""
    if t < PRESLIDING_START:
        return 0.0
    index = min(int((t - PRESLIDING_START) // plateau_dwell),
                levels.size - 1)
    return float(levels[index] * microstep)


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
                  held_command: float | None = None,
                  blocked_stage: bool = False) -> np.ndarray:
    x_d, x_s, v_d, v_s = state[:4]
    if blocked_stage:
        x_s, v_s = 0.0, 0.0
    command = command_position(t, constants["command_step"], constants["plateau_dwell"])
    if held_command is not None:
        command = held_command
    lag = constants["kappa"] * (command - x_d)
    magnetic_force = constants["F_max"] * np.sin(lag)
    detent_force = -(constants["T_det"] / constants["r"]) * np.sin(
        4.0 * constants["kappa"] * x_d + MODEL["detent_phase"]
    )
    electromagnetic_damping = constants["c_m"] * v_d
    axial_force = constants["k_ax"] * (x_d - x_s) + MODEL["c_ax"] * (v_d - v_s)

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
    a_s = (axial_force - friction_generalized[1]) / constants["m_s"]
    derivative[:4] = (v_d, 0.0 if blocked_stage else v_s,
                      a_d, 0.0 if blocked_stage else a_s)
    return derivative


def rk4_case_with_command(case: dict[str, object], constants: dict[str, float],
                          command_function, duration: float,
                          dt: float = PRODUCTION_DT,
                          blocked_stage: bool = False) -> tuple[np.ndarray, np.ndarray]:
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
        k1 = nonlinear_rhs(t, y, case, constants, held_command, blocked_stage)
        k2 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k1, case, constants,
                           held_command, blocked_stage)
        k3 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k2, case, constants,
                           held_command, blocked_stage)
        k4 = nonlinear_rhs(t + dt, y + dt * k3, case, constants,
                           held_command, blocked_stage)
        states[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return times, states


def rk4_case(case: dict[str, object], constants: dict[str, float], dt: float = PRODUCTION_DT,
             duration: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    if duration is None:
        duration = main_duration(constants)
    return rk4_case_with_command(
        case, constants,
        lambda t: command_position(t, constants["command_step"], constants["plateau_dwell"]),
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


def presliding_responses(constants: dict[str, float], keys: tuple[str, str],
                         site: str) -> dict[str, object]:
    """Run one matched LuGre/GMS pair through a settled reversal sequence.

    The guideway experiment uses the normal free-stage plant.  The nut
    identification experiment imposes x_s=0 so that the commanded drive
    coordinate prescribes a measurable nut-port deflection.  This boundary
    condition belongs only to the dedicated B/B2 material test; it does not
    alter the normal free-stage B/B2 simulations.
    """
    microstep = constants["command_step"]
    plateau_dwell = constants["plateau_dwell"]
    blocked_stage = site == "n"
    levels = NUT_PRESLIDING_LEVELS if blocked_stage else GUIDEWAY_PRESLIDING_LEVELS
    duration = PRESLIDING_START + plateau_dwell * levels.size
    command_function = lambda t: presliding_command_position(
        t, microstep, plateau_dwell, levels)
    results: dict[str, np.ndarray] = {}
    forces: dict[str, np.ndarray] = {}
    metrics: dict[str, dict[str, object]] = {}
    times: np.ndarray | None = None

    for key in keys:
        times, states = rk4_case_with_command(
            CASES[key], constants, command_function, duration=duration,
            dt=PRODUCTION_DT, blocked_stage=blocked_stage)
        results[key] = states
        forces[key] = friction_force_history(CASES[key], states, site)

    assert times is not None
    command = np.array([command_function(t) for t in times])
    active = times >= PRESLIDING_START
    for key in keys:
        observed_coordinate = results[key][:, 0] if blocked_stage else results[key][:, 1]
        error = command - observed_coordinate
        site_coordinate = (results[key][:, 1] if site == "g"
                           else results[key][:, 0] - results[key][:, 1])
        endpoint_error = []
        endpoint_force = []
        endpoint_coordinate = []
        for level_index in range(levels.size):
            plateau_end = PRESLIDING_START + (level_index + 1) * plateau_dwell
            window = ((times >= plateau_end - constants["metric_window"])
                      & (times < plateau_end - 0.5e-9))
            if level_index == levels.size - 1:
                window = times >= plateau_end - constants["metric_window"]
            endpoint_error.append(float(np.mean(error[window])))
            endpoint_force.append(float(np.mean(forces[key][window])))
            endpoint_coordinate.append(float(np.mean(site_coordinate[window])))
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
            "endpoint_coordinate_um": np.asarray(endpoint_coordinate) * 1e6,
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
        "plateau_dwell": plateau_dwell,
        "keys": keys,
        "site": site,
        "blocked_stage": blocked_stage,
        "levels": levels,
    }


def final_window_rms_error_nm(times: np.ndarray, states: np.ndarray,
                              constants: dict[str, float]) -> float:
    """Return RMS(command-stage) over the final settled-window samples."""
    command = np.array([
        command_position(t, constants["command_step"], constants["plateau_dwell"])
        for t in times
    ])
    final_window = times >= (times[-1] - constants["metric_window"])
    error = command - states[:, 1]
    return float(np.sqrt(np.mean(error[final_window] ** 2)) * 1e9)


def time_responses(constants: dict[str, float]) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, dict[str, float]]]:
    results: dict[str, np.ndarray] = {}
    metrics: dict[str, dict[str, float]] = {}
    times: np.ndarray | None = None
    for key, case in CASES.items():
        times, states = rk4_case(case, constants, dt=PRODUCTION_DT)
        results[key] = states
    assert times is not None
    command = np.array([
        command_position(t, constants["command_step"], constants["plateau_dwell"])
        for t in times
    ])
    final_window = times >= (times[-1] - constants["metric_window"])
    first_plateau = ((times >= MAIN_START)
                     & (times < MAIN_START + constants["plateau_dwell"]))
    settled_mask = np.zeros(times.size, dtype=bool)
    for level_index in range(MAIN_LEVELS.size):
        plateau_end = MAIN_START + (level_index + 1) * constants["plateau_dwell"]
        settled_mask |= ((times >= plateau_end - constants["metric_window"])
                         & (times < plateau_end - 0.5e-9))
    settled_mask |= final_window
    first_target = MAIN_LEVELS[0] * constants["command_step"]
    for key, states in results.items():
        error = command - states[:, 1]
        first_peak = float(np.max(states[first_plateau, 1]))
        metrics[key] = {
            "mean_final_error_nm": float(np.mean(error[final_window]) * 1e9),
            "rms_final_error_nm": final_window_rms_error_nm(times, states, constants),
            "rms_sequence_deviation_nm": float(np.sqrt(np.mean(error ** 2)) * 1e9),
            "max_abs_deviation_nm": float(np.max(np.abs(error)) * 1e9),
            "rms_settled_deviation_nm": float(np.sqrt(np.mean(error[settled_mask] ** 2)) * 1e9),
            "max_settled_deviation_nm": float(np.max(np.abs(error[settled_mask])) * 1e9),
            "max_stage_um": float(np.max(np.abs(states[:, 1])) * 1e6),
            "first_peak_um": first_peak * 1e6,
            "first_overshoot_pct": max(0.0, (first_peak - first_target)
                                        / abs(first_target) * 100.0),
        }
    return times, command, results, metrics


def gms_step_halving_convergence(constants: dict[str, float], base_times: np.ndarray,
                                 base_results: dict[str, np.ndarray]) -> dict[str, dict[str, object]]:
    """Compare final-window RMS under h, h/2, and h/4 for all GMS cases."""
    study: dict[str, dict[str, object]] = {}
    for key in ("A2", "B2", "C2"):
        rms_values: list[float] = []
        for dt in GMS_CONVERGENCE_DTS:
            if np.isclose(dt, PRODUCTION_DT, rtol=0.0, atol=1.0e-15):
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
             f"Settled-window modeled command-stage deviation (not a servo tracking specification): "
             f"RMS={metric['rms_settled_deviation_nm']:.1f} nm; "
             f"max={metric['max_settled_deviation_nm']:.1f} nm; "
             f"final-window RMS={metric['rms_final_error_nm']:.1f} nm.",
            ha="center", fontsize=8.2, color="#555555",
        )
        fig.text(0.5, 0.012,
                 f"Nonlinear magnetic and periodic detent force; zeta_m={MODEL['zeta_m']:.2f}; "
                 f"{constants['plateau_dwell'] * 1e3:.0f} ms derived dwell; "
                 f"one command quantum = {constants['command_step'] * 1e9:.2f} nm.",
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


def plot_presliding_memory(experiment: dict[str, object], output_name: str) -> Path:
    """Visualize nested-reversal command following and friction return-point memory."""
    times = experiment["times"]
    command = experiment["command"]
    results = experiment["results"]
    forces = experiment["forces"]
    metrics = experiment["metrics"]
    keys = experiment["keys"]
    site = experiment["site"]
    site_title = "Guideway" if site == "g" else "Nut microslip"
    blocked_stage = experiment["blocked_stage"]
    time_ms = times * 1e3

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.8))
    ax_motion, ax_error = axes[0]
    ax_memory, ax_metrics = axes[1]

    ax_motion.step(time_ms, command * 1e6, where="post", color="#111111",
                   linewidth=2.0, label="Command")
    for key in keys:
        case = CASES[key]
        stage = results[key][:, 1]
        observed_coordinate = results[key][:, 0] if blocked_stage else stage
        site_coordinate = stage if site == "g" else results[key][:, 0] - stage
        ax_motion.plot(time_ms, observed_coordinate * 1e6, color=case["color"],
                       linestyle=case["ls"], linewidth=1.5, label=case["label"])
        ax_error.plot(time_ms, (command - observed_coordinate) * 1e9, color=case["color"],
                      linestyle=case["ls"], linewidth=1.35, label=case["label"])
        ax_memory.plot(site_coordinate * 1e6, forces[key], color=case["color"],
                       linestyle=case["ls"], linewidth=0.7, alpha=0.20,
                       label=case["label"])
        ax_memory.plot(metrics[key]["endpoint_coordinate_um"],
                       metrics[key]["endpoint_force_N"], color=case["color"],
                       linestyle=case["ls"], linewidth=1.5,
                       marker="o" if key == keys[0] else "s",
                       markersize=4.0, markerfacecolor="white")

    ax_motion.set_title("Nested command and drive coordinate; stage blocked"
                        if blocked_stage else "Nested microstep command and actual stage motion")
    ax_motion.set_ylabel("Position (um)")
    ax_motion.legend(loc="upper right", fontsize=8)
    ax_error.set_title("Modeled command-drive deviation over the same reversal history"
                       if blocked_stage else "Modeled command-stage deviation over the same reversal history")
    ax_error.set_ylabel(r"Modeled deviation $x_{cmd}-x_d$ (nm)" if blocked_stage
                        else r"Modeled deviation $x_{cmd}-x_s$ (nm)")
    ax_error.axhline(0.0, color="#777777", linewidth=0.8)
    ax_error.legend(loc="upper right", fontsize=8)
    ax_memory.set_title(f"{site_title} friction memory loops")
    ax_memory.set_xlabel("Stage position (um)" if site == "g" else r"Nut-port deflection $x_d-x_s$ (um)")
    ax_memory.set_ylabel(f"{site_title} friction force (N)")
    ax_memory.axhline(0.0, color="#888888", linewidth=0.7)
    ax_memory.axvline(0.0, color="#888888", linewidth=0.7)
    ax_memory.legend(loc="best", fontsize=8)

    categories = ("Whole-sequence\nRMS", "Peak absolute\ndeviation",
                  "Return-point\nmismatch", "Final-origin\nabsolute deviation")
    x_positions = np.arange(len(categories), dtype=float)
    width = 0.34
    lugre_values = np.array([
        metrics[keys[0]]["whole_rms_nm"],
        metrics[keys[0]]["max_abs_deviation_nm"],
        metrics[keys[0]]["return_error_mismatch_nm"],
        abs(metrics[keys[0]]["final_mean_nm"]),
    ])
    gms_values = np.array([
        metrics[keys[1]]["whole_rms_nm"],
        metrics[keys[1]]["max_abs_deviation_nm"],
        metrics[keys[1]]["return_error_mismatch_nm"],
        abs(metrics[keys[1]]["final_mean_nm"]),
    ])
    bars_a = ax_metrics.bar(x_positions - width / 2.0, lugre_values, width,
                            color=CASES[keys[0]]["color"], label=f"LuGre {keys[0]}")
    bars_a2 = ax_metrics.bar(x_positions + width / 2.0, gms_values, width,
                             color=CASES[keys[1]]["color"], label=f"GMS {keys[1]}")
    ax_metrics.set_yscale("log")
    ax_metrics.set_xticks(x_positions, categories)
    ax_metrics.set_ylabel("Command-output deviation metric (nm, log scale)")
    ax_metrics.set_title("Open-loop response: global and memory-sensitive metrics")
    ax_metrics.legend(loc="upper right", fontsize=8)
    ax_metrics.text(
        0.02, 0.03,
        f"Primary force metric: return mismatch "
        f"{metrics[keys[0]]['return_force_mismatch_N']:.4f} N LuGre, "
        f"{metrics[keys[1]]['return_force_mismatch_N']:.4f} N GMS",
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

    max_force = max(metrics[key]["max_force_N"] for key in keys)
    macro_fraction = 100.0 * max_force / FRICTION[site]["F_s"]
    fig.suptitle(f"{site_title} partial-slip experiment: LuGre versus GMS",
                 fontsize=15, fontweight="bold")
    fig.text(
        0.5, 0.012,
        f"1 STEP/DIR quantum = {experiment['microstep'] * 1e9:.2f} nm = 1/{MODEL['microstep_divisor']} full step; "
        f"peak friction = {max_force:.3f} N ({macro_fraction:.1f}% of macro breakaway). "
        f"Each plateau is {experiment['plateau_dwell'] * 1e3:.0f} ms; markers are "
        f"{min(20.0, 0.2 * experiment['plateau_dwell'] * 1e3):.0f} ms settled-window means.",
        ha="center", fontsize=8.4, color="#555555",
    )
    fig.tight_layout(rect=(0.02, 0.045, 0.99, 0.95), h_pad=2.0, w_pad=1.5)
    output = ASSET_DIR / output_name
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_kinematic_diagram() -> tuple[Path, Path, Path]:
    """Render separate ten-DOF topology and retained two-DOF figures.

    Geometry is intentionally symbol-only.  The topology, compliance shares,
    and collapse checks are all sourced from MODEL/FULL through the same
    derivation functions used by the Markdown registry.
    """
    constants = physical_constants()
    parameters = full_parameters()
    required_registry_keys = {
        "lead", "screw_length", "screw_diameter", "screw_density",
        "reduced_stage_mass", "reduced_axial_stiffness", "k_c_series",
        "k_theta_a", "k_theta_b", "k_brg", "k_sha", "k_shb",
        "k_ball", "k_mnt", "screw_inertia", "total_rotational_inertia",
        "reduced_drive_mass",
    }
    missing_registry_keys = required_registry_keys.difference(PARAMETER_REGISTRY)
    if missing_registry_keys:
        raise ValueError(
            "Kinematic figure has unregistered parameters: "
            + ", ".join(sorted(missing_registry_keys))
        )
    if not np.isclose(parameters["m_d_reflected"], constants["m_d"]):
        raise ValueError("Full and reduced reflected drive masses have drifted")
    axial_elements = ("k_brg", "k_sha", "k_ball", "k_mnt")
    compliance = np.array([1.0 / parameters[key] for key in axial_elements])
    compliance_shares = compliance / np.sum(compliance)
    if not np.isclose(np.sum(compliance_shares), 1.0):
        raise ValueError("Axial compliance bar does not close")

    drive_color = "#dceef6"
    stage_color = "#dff2ea"
    dropped_color = "#eeeeee"
    spring_color = "#d97800"
    friction_color = "#b23a48"
    damping_color = "#6a4c93"
    rigid_color = "#39434d"
    ground_color = "#7b858c"
    annotation_color = "#59636d"
    detent_color = "#c08a00"
    detent_fill = "#fff2c7"
    TITLE_FS = 15.5
    BOX_FS = 10.5
    ANNO_FS = 8.2

    def setup(width: float, height: float, xlim: float, ylim: float) -> tuple[plt.Figure, plt.Axes]:
        figure, axis = plt.subplots(figsize=(width, height))
        axis.axis("off")
        axis.set_xlim(0.0, xlim)
        axis.set_ylim(0.0, ylim)
        return figure, axis

    def node(ax: plt.Axes, x: float, y: float, label: str, color: str,
             subtitle: str = "", width: float = 1.02, height: float = 0.66,
             alpha: float = 1.0) -> None:
        patch = FancyBboxPatch(
            (x - width / 2.0, y - height / 2.0), width, height,
            boxstyle="round,pad=0.04", facecolor=color,
            edgecolor=rigid_color, linewidth=1.2, zorder=5, alpha=alpha)
        ax.add_patch(patch)
        ax.text(x, y + (0.08 if subtitle else 0.0), label,
                ha="center", va="center", fontsize=BOX_FS, zorder=6,
                alpha=alpha)
        if subtitle:
            ax.text(x, y - 0.20, subtitle, ha="center", va="center",
                    fontsize=ANNO_FS, color=annotation_color, zorder=6,
                    alpha=alpha)

    def spring(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
               color: str = spring_color, amplitude: float = 0.075,
               linewidth: float = 1.80, alpha: float = 1.0) -> None:
        p0, p1 = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
        delta = p1 - p0
        length = np.linalg.norm(delta)
        if length <= 0.0:
            raise ValueError("Spring endpoints must differ")
        unit = delta / length
        normal = np.array([-unit[1], unit[0]])
        terminal = min(0.18, 0.22 * length)
        active_start = p0 + terminal * unit
        active_end = p1 - terminal * unit
        active_delta = active_end - active_start
        points = [p0, active_start]
        for index in range(10):
            fraction = (index + 1.0) / 11.0
            sign = 1.0 if index % 2 == 0 else -1.0
            points.append(active_start + fraction * active_delta
                          + sign * amplitude * normal)
        points.extend([active_end, p1])
        points = np.asarray(points)
        ax.plot(points[:, 0], points[:, 1], color=color,
                linewidth=linewidth, solid_capstyle="round", zorder=3,
                alpha=alpha)

    def dashpot(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
                width: float = 0.11) -> None:
        p0, p1 = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
        delta = p1 - p0
        unit = delta / np.linalg.norm(delta)
        normal = np.array([-unit[1], unit[0]])
        a, b, c, d = (p0 + fraction * delta for fraction in (0.25, 0.39, 0.64, 0.78))
        ax.plot(*zip(p0, a), color=damping_color, linewidth=1.5, zorder=3)
        ax.plot(*zip(a - width * normal, a + width * normal),
                color=damping_color, linewidth=1.5, zorder=3)
        ax.plot(*zip(a, c), color=damping_color, linewidth=1.5, zorder=3)
        ax.plot(*zip(b - width * normal, d - width * normal),
                color=damping_color, linewidth=1.5, zorder=3)
        ax.plot(*zip(b + width * normal, d + width * normal),
                color=damping_color, linewidth=1.5, zorder=3)
        ax.plot(*zip(d - width * normal, d + width * normal),
                color=damping_color, linewidth=1.5, zorder=3)
        ax.plot(*zip(d, p1), color=damping_color, linewidth=1.5, zorder=3)

    def friction(ax: plt.Axes, start: tuple[float, float],
                 end: tuple[float, float], label: str,
                 block_width: float = 0.70) -> None:
        x0, y0 = start
        x1, y1 = end
        horizontal = abs(x1 - x0) >= abs(y1 - y0)
        cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        if horizontal:
            box_w = min(block_width, 0.46 * abs(x1 - x0))
            box_h = 0.26
            ax.plot([x0, cx - box_w / 2.0], [y0, cy],
                    color=friction_color, linewidth=1.45, zorder=3)
            ax.plot([cx + box_w / 2.0, x1], [cy, y1],
                    color=friction_color, linewidth=1.45, zorder=3)
        else:
            box_w = 0.48
            box_h = min(0.44, 0.42 * abs(y1 - y0))
            ax.plot([x0, cx], [y0, cy + box_h / 2.0],
                    color=friction_color, linewidth=1.45, zorder=3)
            ax.plot([cx, x1], [cy - box_h / 2.0, y1],
                    color=friction_color, linewidth=1.45, zorder=3)
        patch = FancyBboxPatch(
            (cx - box_w / 2.0, cy - box_h / 2.0), box_w, box_h,
            boxstyle="round,pad=0.025", facecolor="#f8dce1",
            edgecolor=friction_color, linewidth=1.1, zorder=4)
        ax.add_patch(patch)
        ax.text(cx, cy, label, ha="center", va="center",
                fontsize=ANNO_FS, color="#8d2936", zorder=5)

    def ground_baseline(ax: plt.Axes, x0: float, x1: float, y: float,
                        label: str) -> None:
        ax.plot([x0, x1], [y, y], color=ground_color,
                linewidth=1.8, zorder=1)
        for x_pos in np.arange(x0 + 0.10, x1 + 1e-9, 0.26):
            ax.plot([x_pos, x_pos - 0.10], [y, y - 0.11],
                    color=ground_color, linewidth=0.7, zorder=1)
        ax.text(x0, y - 0.23, label, ha="left", va="top",
                fontsize=ANNO_FS, color=annotation_color)

    def grounded_friction(ax: plt.Axes, x: float, source_y: float,
                          baseline_y: float, block_y: float, label: str) -> None:
        box_w, box_h = 0.50, 0.34
        ax.plot([x, x], [source_y, block_y + box_h / 2.0],
                color=ground_color, linewidth=1.15, zorder=2)
        ax.plot([x, x], [block_y - box_h / 2.0, baseline_y],
                color=ground_color, linewidth=1.15, zorder=2)
        patch = FancyBboxPatch(
            (x - box_w / 2.0, block_y - box_h / 2.0), box_w, box_h,
            boxstyle="round,pad=0.025", facecolor="#f8dce1",
            edgecolor=friction_color, linewidth=1.1, zorder=4)
        ax.add_patch(patch)
        ax.text(x, block_y, label, ha="center", va="center",
                fontsize=ANNO_FS, color="#8d2936", zorder=5)

    def moving_wall(ax: plt.Axes, x: float, y: float, label: str,
                    height: float = 0.82) -> None:
        """Draw the imposed command as a moving hatched boundary."""
        y0, y1 = y - height / 2.0, y + height / 2.0
        ax.plot([x, x], [y0, y1], color="#72569a",
                linewidth=2.0, zorder=4)
        for ypos in np.linspace(y0 + 0.04, y1 - 0.04, 6):
            ax.plot([x, x - 0.12], [ypos, ypos - 0.09],
                    color="#72569a", linewidth=0.8, zorder=4)
        ax.add_patch(FancyArrowPatch(
            (x - 0.28, y1 + 0.12), (x + 0.28, y1 + 0.12),
            arrowstyle="<->", mutation_scale=8, color="#72569a",
            linewidth=1.0, zorder=4))
        ax.text(x, y1 + 0.29, label, ha="center", va="bottom",
                fontsize=ANNO_FS, color="#5b437c", fontweight="bold")

    def sensor(ax: plt.Axes, x: float, y: float) -> None:
        circle = Circle((x, y), 0.18, facecolor="white",
                        edgecolor="#277da1", linewidth=1.3, zorder=6)
        ax.add_patch(circle)
        ax.plot([x - 0.09, x + 0.08], [y - 0.08, y + 0.09],
                color="#277da1", linewidth=1.2, zorder=7)
        ax.add_patch(FancyArrowPatch(
            (x - 0.09, y - 0.08), (x + 0.08, y + 0.09),
            arrowstyle="-|>", mutation_scale=7, color="#277da1",
            linewidth=1.0, zorder=7))

    # ------------------------------------------------------------------
    # Figure 1: full ten-coordinate topology.
    # ------------------------------------------------------------------
    fig1, ax1 = setup(17.2, 7.9, 17.2, 10.0)
    ax1.set_ylim(2.15, 9.92)
    ax1.text(8.6, 9.72, "Figure 1 — Ten-DOF physical topology",
             ha="center", fontsize=TITLE_FS, fontweight="bold")

    columns = {
        "ground": 0.75,
        "motor": 2.25,
        "coupling": 4.05,
        "bearing + screw drive end": 5.85,
        "screw nut plane": 8.15,
        "beyond nut": 10.15,
        "nut body": 12.05,
        "stage": 14.00,
        "guideways": 15.85,
    }
    for label, x_pos in columns.items():
        ax1.plot([x_pos, x_pos], [2.65, 9.05], color="#dfe3e6",
                 linewidth=0.7, linestyle=(0, (2, 4)), zorder=0)
        ax1.text(
            x_pos, 9.08, label, ha="center", va="center",
            fontsize=ANNO_FS, color=annotation_color,
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "#f5f6f7",
                  "edgecolor": "#d4d8dc"})

    y_torsion = 8.10
    y_torsion_stub = 7.20
    y_torsion_ground = 6.42
    y_transformer = 5.82
    y_axial = 4.72
    y_axial_stub = 3.50
    y_nut_friction = 3.98
    y_axial_ground = 2.78
    for label, y_pos in (
        (r"Torsional coordinates $\theta$ [rad]", y_torsion),
        ("torsional constitutive stubs", y_torsion_stub),
        ("transformer bridge", y_transformer),
        (r"Axial coordinates $u,x$ [m]", y_axial),
        ("axial constitutive stubs", y_axial_stub),
    ):
        ax1.plot([0.45, 16.75], [y_pos, y_pos],
                 color="#edf0f2", linewidth=0.65, zorder=0)
        ax1.text(0.40, y_pos, label, ha="right", va="center",
                 fontsize=ANNO_FS, color="#7a838a")

    ground_baseline(
        ax1, columns["motor"] - 0.55, columns["screw nut plane"] + 0.55,
        y_torsion_ground, "torsional ground datum")
    ground_baseline(
        ax1, columns["bearing + screw drive end"] - 0.70,
        columns["guideways"] + 0.45, y_axial_ground,
        "axial ground datum")

    # Command input is an imposed moving datum, not another coordinate.
    moving_wall(ax1, columns["ground"], y_torsion, r"$x_{cmd}$")
    spring(
        ax1, (columns["ground"], y_torsion),
        (columns["motor"] - 0.56, y_torsion), color=detent_color)
    ax1.text(1.48, y_torsion + 0.24, r"$K_m$", ha="center",
             fontsize=ANNO_FS, color="#8a6200")
    ax1.plot(1.47, y_torsion - 0.23, marker="D", markersize=6.0,
             markerfacecolor=detent_fill, markeredgecolor=detent_color,
             zorder=5)
    ax1.text(1.47, y_torsion - 0.48,
             r"$K_{det}$: periodic conservative tangent",
             ha="center", fontsize=ANNO_FS, color="#8a6200")

    rotational_nodes = (
        ("motor", r"$\theta_m$", "q1"),
        ("coupling", r"$\theta_c$", "q2"),
        ("bearing + screw drive end", r"$\theta_{s1}$", "q3"),
        ("screw nut plane", r"$\theta_{s2}$", "q4"),
        ("beyond nut", r"$\theta_{s3}$", "q5"),
    )
    for column, label, index in rotational_nodes:
        node(ax1, columns[column], y_torsion, label, drive_color, index)
    rotational_springs = (
        ("motor", "coupling", r"$k_{c1}$"),
        ("coupling", "bearing + screw drive end", r"$k_{c2}$"),
        ("bearing + screw drive end", "screw nut plane", r"$k_{\theta a}$"),
        ("screw nut plane", "beyond nut", r"$k_{\theta b}$"),
    )
    for left, right, label in rotational_springs:
        x0 = columns[left] + 0.56
        x1 = columns[right] - 0.56
        spring(ax1, (x0, y_torsion), (x1, y_torsion))
        ax1.text((columns[left] + columns[right]) / 2.0,
                 y_torsion + 0.24, label, ha="center",
                 fontsize=ANNO_FS, color="#9a5600")

    # Hub losses use one fixed stub row; grounded losses use neutral drops.
    for left, right, label in (
        ("motor", "coupling", r"$T_{h1}$"),
        ("coupling", "bearing + screw drive end", r"$T_{h2}$"),
    ):
        x0 = columns[left] + 0.56
        x1 = columns[right] - 0.56
        ax1.plot([x0, x0], [y_torsion - 0.33, y_torsion_stub],
                 color=ground_color, linewidth=0.9)
        ax1.plot([x1, x1], [y_torsion - 0.33, y_torsion_stub],
                 color=ground_color, linewidth=0.9)
        friction(ax1, (x0, y_torsion_stub), (x1, y_torsion_stub), label)
    grounded_friction(
        ax1, columns["motor"], y_torsion - 0.33,
        y_torsion_ground, 6.83, r"$T_{mb}$")
    grounded_friction(
        ax1, columns["bearing + screw drive end"], y_torsion - 0.33,
        y_torsion_ground, 6.83, r"$T_{brg}$")
    grounded_friction(
        ax1, columns["screw nut plane"], y_torsion - 0.33,
        y_torsion_ground, 6.83, r"$T_{f,r}$")
    # Transformer and axial coordinate rail.
    tf = FancyBboxPatch(
        (columns["screw nut plane"] - 0.43, y_transformer - 0.29),
        0.86, 0.58, boxstyle="round,pad=0.04",
        facecolor="#fff2dc", edgecolor=spring_color,
        linewidth=1.25, zorder=4)
    ax1.add_patch(tf)
    ax1.text(columns["screw nut plane"], y_transformer + 0.07,
             "TF", ha="center", fontsize=BOX_FS,
             fontweight="bold", color="#9a5600", zorder=6)
    ax1.text(columns["screw nut plane"], y_transformer - 0.15,
             r"$r\theta_{s2}$", ha="center", fontsize=ANNO_FS,
             color="#9a5600", zorder=6)
    ax1.add_patch(FancyArrowPatch(
        (columns["screw nut plane"], y_torsion - 0.34),
        (columns["screw nut plane"], y_transformer + 0.31),
        arrowstyle="-|>", mutation_scale=8, color=ground_color,
        linewidth=1.1))

    node(ax1, columns["bearing + screw drive end"], y_axial,
         r"$u_b$", dropped_color, "q6", alpha=0.45)
    node(ax1, columns["screw nut plane"], y_axial,
         r"$u_e$", dropped_color, "q7", alpha=0.45)
    node(ax1, columns["beyond nut"], y_axial_stub,
         r"$u_f$", dropped_color, "q8", alpha=0.35)
    node(ax1, columns["nut body"], y_axial,
         r"$u_n$", stage_color, "q9")
    node(ax1, columns["stage"], y_axial,
         r"$x_s$", stage_color, "q10", width=1.12)

    spring(
        ax1, (columns["bearing + screw drive end"] + 0.55, y_axial),
        (columns["screw nut plane"] - 0.55, y_axial))
    ax1.text(
        (columns["bearing + screw drive end"] + columns["screw nut plane"]) / 2.0,
        y_axial + 0.24, r"$k_{sha}$", ha="center",
        fontsize=ANNO_FS, color="#9a5600")

    sum_x = columns["screw nut plane"] + 0.76
    ax1.plot(
        [columns["screw nut plane"] + 0.55, sum_x],
        [y_axial, y_axial], color=ground_color, linewidth=1.3)
    ax1.plot(
        [columns["screw nut plane"], sum_x],
        [y_transformer - 0.31, y_transformer - 0.31],
        color=ground_color, linewidth=1.1)
    ax1.plot(
        [sum_x, sum_x], [y_transformer - 0.31, y_axial],
        color=ground_color, linewidth=1.1)
    ax1.add_patch(Circle(
        (sum_x, y_axial), 0.065, facecolor=rigid_color,
        edgecolor="white", linewidth=0.6, zorder=7))
    ax1.text(sum_x + 0.12, y_axial + 0.20, r"$\Sigma$",
             fontsize=ANNO_FS, color=rigid_color, fontweight="bold")
    ax1.text(sum_x, y_transformer - 0.52,
             r"$u_t=u_e+r\theta_{s2}$", ha="center",
             fontsize=ANNO_FS, color="#9a5600")

    ball_end = columns["nut body"] - 0.55
    spring(ax1, (sum_x + 0.08, y_axial), (ball_end, y_axial))
    ax1.text((sum_x + ball_end) / 2.0, y_axial + 0.24,
             r"$k_{ball}$", ha="center", fontsize=ANNO_FS,
             color="#9a5600")
    spring(
        ax1, (columns["nut body"] + 0.55, y_axial),
        (columns["stage"] - 0.60, y_axial))
    ax1.text(
        (columns["nut body"] + columns["stage"]) / 2.0,
        y_axial + 0.24, r"$k_{mnt}$", ha="center",
        fontsize=ANNO_FS, color="#9a5600")

    # Beyond-nut overhang, internal nut port, bearing support, and guideways.
    ax1.plot(
        [columns["screw nut plane"], columns["screw nut plane"]],
        [y_axial - 0.33, y_axial_stub], color=ground_color,
        linewidth=0.9, alpha=0.45)
    spring(
        ax1, (columns["screw nut plane"] + 0.08, y_axial_stub),
        (columns["beyond nut"] - 0.55, y_axial_stub),
        color="#8a8a8a", amplitude=0.055, alpha=0.40)
    ax1.text(
        (columns["screw nut plane"] + columns["beyond nut"]) / 2.0,
        y_axial_stub + 0.22, r"$k_{shb}$", ha="center",
        fontsize=ANNO_FS, color="#777777", alpha=0.75)
    ax1.plot([sum_x, sum_x], [y_axial - 0.06, y_nut_friction],
             color=ground_color, linewidth=0.9)
    ax1.plot(
        [columns["nut body"], columns["nut body"]],
        [y_axial - 0.33, y_nut_friction], color=ground_color,
        linewidth=0.9)
    friction(
        ax1, (sum_x, y_nut_friction),
        (columns["nut body"], y_nut_friction), r"$F_{f,n}$")

    support_x = columns["bearing + screw drive end"]
    spring(
        ax1, (support_x, y_axial - 0.33),
        (support_x, y_axial_ground + 0.04), amplitude=0.055)
    ax1.text(support_x - 0.20, 3.48, r"$k_{brg}$",
             ha="right", fontsize=ANNO_FS, color="#9a5600")
    terminal_x = columns["guideways"]
    ax1.plot(
        [columns["stage"] + 0.60, terminal_x],
        [y_axial, y_axial], color=ground_color, linewidth=1.80,
        solid_capstyle="round", zorder=3)
    ax1.add_patch(Circle(
        (terminal_x, y_axial), 0.055, facecolor=ground_color,
        edgecolor="white", linewidth=0.5, zorder=6))
    for pad_x in np.linspace(columns["stage"] + 0.88,
                             columns["guideways"] - 0.18, 4):
        ax1.add_patch(Rectangle(
            (pad_x - 0.075, y_axial - 0.16), 0.15, 0.16,
            facecolor="#f8dce1", edgecolor=friction_color,
            linewidth=0.9, zorder=4))
    grounded_friction(
        ax1, terminal_x, y_axial,
        y_axial_ground, 3.38, r"$F_{f,g}$")
    sensor_x, sensor_y = columns["stage"] + 0.82, y_axial + 0.62
    ax1.plot(
        [columns["stage"] + 0.56, sensor_x - 0.18],
        [y_axial + 0.20, sensor_y], color="#277da1",
        linewidth=1.2, solid_capstyle="round", zorder=4)
    sensor(ax1, sensor_x, sensor_y)


    fig1.subplots_adjust(left=0.02, right=0.985, top=0.985, bottom=0.02)
    topology_output = ASSET_DIR / "kinematic_diagram.svg"
    save_svg(fig1, topology_output)
    plt.close(fig1)

    # ------------------------------------------------------------------
    # Figure 2: retained two-coordinate plant, drawn independently.
    # ------------------------------------------------------------------
    fig2, ax2 = setup(13.4, 5.8, 13.4, 7.4)
    ax2.set_ylim(1.02, 7.32)
    ax2.text(6.70, 7.08, "Figure 2 — Retained two-DOF model",
             ha="center", fontsize=TITLE_FS, fontweight="bold")

    x_input, x_drive, x_stage = 0.72, 3.85, 9.45
    y_mass = 4.18
    mass_width, mass_height = 2.05, 2.18
    node(ax2, x_drive, y_mass, r"$x_d$", drive_color,
         r"reflected drivetrain $m_d$", width=mass_width,
         height=mass_height)
    node(ax2, x_stage, y_mass, r"$x_s$", stage_color,
         r"stage assembly $m_s$", width=mass_width,
         height=mass_height)

    moving_wall(ax2, x_input, y_mass, r"$x_{cmd}$", height=1.02)
    spring(
        ax2, (x_input, y_mass),
        (x_drive - mass_width / 2.0, y_mass), color=detent_color,
        amplitude=0.08)
    ax2.text(1.84, y_mass + 0.27,
             r"$K_m(x_{cmd}-x_d)$", ha="center",
             fontsize=ANNO_FS, color="#8a6200")

    left_edge = x_drive + mass_width / 2.0
    right_edge = x_stage - mass_width / 2.0
    element_rows = {
        "spring": y_mass + 0.68,
        "damper": y_mass,
        "friction": y_mass - 0.68,
    }
    spring(
        ax2, (left_edge, element_rows["spring"]),
        (right_edge, element_rows["spring"]), amplitude=0.085)
    dashpot(
        ax2, (left_edge, element_rows["damper"]),
        (right_edge, element_rows["damper"]), width=0.12)
    friction(
        ax2, (left_edge, element_rows["friction"]),
        (right_edge, element_rows["friction"]), r"$F_{f,n}$",
        block_width=0.86)
    midpoint = (left_edge + right_edge) / 2.0
    ax2.text(midpoint, element_rows["spring"] + 0.24,
             r"$k_{ax}$", ha="center", fontsize=ANNO_FS,
             color="#9a5600")
    ax2.text(midpoint, element_rows["damper"] + 0.23,
             r"$c_{ax}$", ha="center", fontsize=ANNO_FS,
             color=damping_color)
    ax2.text(
        midpoint, element_rows["friction"] - 0.28,
        r"internal equal-and-opposite port: $v_n=\dot x_d-\dot x_s$",
        ha="center", fontsize=ANNO_FS, color="#8d2936")

    y_ground = 1.32
    ground_baseline(ax2, 2.35, 10.95, y_ground,
                    "retained-model ground datum")
    grounded_friction(
        ax2, x_drive - 0.32, y_mass - mass_height / 2.0,
        y_ground, 2.15, r"$F_{f,d}$")
    dashpot(
        ax2, (x_drive + 0.38, y_mass - mass_height / 2.0),
        (x_drive + 0.38, y_ground), width=0.07)
    ax2.text(x_drive + 0.62, 2.16, r"$c_m$", ha="left",
             va="center", fontsize=ANNO_FS, color=damping_color)
    grounded_friction(
        ax2, x_stage, y_mass - mass_height / 2.0,
        y_ground, 2.15, r"$F_{f,g}$")
    sensor_x, sensor_y = x_stage + 1.28, y_mass + 0.72
    ax2.plot(
        [x_stage + mass_width / 2.0, sensor_x - 0.18],
        [y_mass + 0.30, sensor_y], color="#277da1",
        linewidth=1.2, solid_capstyle="round", zorder=4)
    sensor(ax2, sensor_x, sensor_y)

    # The diamond is defined once in the shared legend.
    detent_x, detent_y = x_drive + 0.55, y_mass + mass_height / 2.0
    ax2.plot(detent_x, detent_y, marker="D", markersize=6.5,
             markerfacecolor=detent_fill, markeredgecolor=detent_color,
             zorder=7)

    fig2.subplots_adjust(left=0.025, right=0.985, top=0.985, bottom=0.025)
    reduced_output = ASSET_DIR / "kinematic_diagram_reduced.svg"
    save_svg(fig2, reduced_output)
    plt.close(fig2)

    # ------------------------------------------------------------------
    # Shared legend: one key for both independent diagrams.
    # ------------------------------------------------------------------
    fig3, ax3 = setup(15.2, 4.1, 15.2, 4.1)
    ax3.text(7.60, 3.86, "Shared kinematic-diagram legend",
             ha="center", fontsize=TITLE_FS, fontweight="bold")

    ax3.text(0.35, 3.45, "Coordinate aggregation",
             fontsize=BOX_FS, fontweight="bold", color="#425b6b")
    aggregation = (
        (0.35, drive_color,
         r"$\theta_m,\theta_c,\theta_{s1},\theta_{s2},\theta_{s3}"
         r"\ \rightarrow\ m_d,x_d$"),
        (3.05, dropped_color,
         r"$u_b,u_e,u_f\ \rightarrow\ k_{ax}$ path; inertia omitted"),
        (5.75, stage_color,
         r"$u_n,x_s\ \rightarrow\ m_s,x_s$"),
    )
    for xpos, color, label in aggregation:
        ax3.add_patch(FancyBboxPatch(
            (xpos, 2.86), 2.48, 0.42, boxstyle="round,pad=0.04",
            facecolor=color, edgecolor="#aab3b9", linewidth=0.9))
        ax3.text(xpos + 1.24, 3.07, label, ha="center", va="center",
                 fontsize=ANNO_FS, color="#3e474e")

    ax3.text(0.35, 2.48, "Elements and boundaries",
             fontsize=BOX_FS, fontweight="bold", color="#425b6b")
    spring(ax3, (0.45, 1.93), (1.65, 1.93), amplitude=0.065)
    ax3.text(1.05, 1.66, r"spring $k_j$", ha="center",
             fontsize=ANNO_FS, color="#9a5600")
    dashpot(ax3, (2.05, 1.93), (3.25, 1.93), width=0.09)
    ax3.text(2.65, 1.66, r"damper $c_j$", ha="center",
             fontsize=ANNO_FS, color=damping_color)
    friction(ax3, (3.65, 1.93), (4.85, 1.93), r"$F_f$",
             block_width=0.58)
    ax3.text(4.25, 1.66, "friction port", ha="center",
             fontsize=ANNO_FS, color="#8d2936")
    moving_wall(ax3, 5.45, 1.93, r"$x_{cmd}$", height=0.62)
    ground_baseline(ax3, 6.10, 7.25, 1.93, "fixed datum")
    sensor(ax3, 7.78, 1.93)
    ax3.text(7.78, 1.66, "sensor", ha="center",
             fontsize=ANNO_FS, color="#277da1")
    ax3.plot(8.42, 1.93, marker="D", markersize=6.5,
             markerfacecolor=detent_fill, markeredgecolor=detent_color)
    ax3.text(8.62, 1.93, r"periodic $F_{det}$",
             ha="left", va="center", fontsize=ANNO_FS, color="#8a6200")

    ax3.text(0.35, 0.96, "Drive-side reduction",
             fontsize=BOX_FS, fontweight="bold", color="#425b6b")
    ax3.text(
        0.35, 0.54,
        r"$\{T_{mb},T_{h1},T_{h2},T_{brg},T_{f,r}\}"
        r"\ \rightarrow\ F_{f,d}$",
        fontsize=ANNO_FS, color="#8d2936",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "#fff7f8",
              "edgecolor": "#d79aa3"})
    ax3.text(
        3.80, 0.54,
        r"$F_{f,n}$: internal equal–opposite port"
        r"        $F_{f,g}$: stage–ground port",
        fontsize=ANNO_FS, color="#8d2936")

    ax3.text(9.52, 3.45, "Series-compliance share",
             fontsize=BOX_FS, fontweight="bold", color="#425b6b")
    bar_left, bar_right = 9.52, 14.85
    bar_y, bar_h = 3.05, 0.32
    colors = ("#f4c27a", "#ed9f4a", "#d97800", "#b65f00")
    labels = (r"$k_{brg}$", r"$k_{sha}$", r"$k_{ball}$", r"$k_{mnt}$")
    cursor = bar_left
    for share, label, color in zip(compliance_shares, labels, colors):
        width = (bar_right - bar_left) * float(share)
        ax3.add_patch(Rectangle(
            (cursor, bar_y - bar_h / 2.0), width, bar_h,
            facecolor=color, edgecolor="white", linewidth=0.8, zorder=4))
        ax3.text(cursor + width / 2.0, bar_y, label,
                 ha="center", va="center", fontsize=ANNO_FS,
                 color="#332619", zorder=5)
        cursor += width
    ax3.add_patch(Rectangle(
        (bar_left, bar_y - bar_h / 2.0), bar_right - bar_left, bar_h,
        fill=False, edgecolor="#9a5600", linewidth=0.9, zorder=5))

    ax3.text(9.52, 2.48, "Friction port / case matrix",
             fontsize=BOX_FS, fontweight="bold", color="#425b6b")
    case_labels = ("0", "A/A2", "B/B2", "C/C2")
    port_labels = (r"$F_{f,d}$", r"$F_{f,n}$", r"$F_{f,g}$")
    active = (
        (0, 1, 1, 1),
        (0, 0, 1, 1),
        (0, 1, 0, 1),
    )
    for column_index, case_label in enumerate(case_labels):
        ax3.text(11.75 + 0.80 * column_index, 2.17,
                 case_label, ha="center", fontsize=ANNO_FS,
                 color="#4d555c")
    for row_index, (port_label, row_active) in enumerate(zip(port_labels, active)):
        ypos = 1.84 - 0.34 * row_index
        ax3.text(9.58, ypos, port_label, ha="left", va="center",
                 fontsize=ANNO_FS, color="#4d555c")
        for column_index, enabled in enumerate(row_active):
            xpos = 11.75 + 0.80 * column_index
            ax3.add_patch(Circle(
                (xpos, ypos), 0.065,
                facecolor=friction_color if enabled else "white",
                edgecolor=friction_color if enabled else "#aeb6bc",
                linewidth=0.85))
    ax3.plot(9.62, 0.56, marker="D", markersize=6.0,
             markerfacecolor=detent_fill, markeredgecolor=detent_color)
    ax3.text(9.84, 0.56,
             r"$F_{det}$ is conservative and excluded from the matrix",
             ha="left", va="center", fontsize=ANNO_FS, color="#8a6200")

    fig3.subplots_adjust(left=0.02, right=0.985, top=0.98, bottom=0.04)
    legend_output = ASSET_DIR / "kinematic_diagram_legend.svg"
    save_svg(fig3, legend_output)
    plt.close(fig3)
    return topology_output, reduced_output, legend_output
def _flow_section_url(section: str) -> str:
    # SVGs live in rendered_assets/, while the analytical HTML is one level up.
    return f"../Analytical_derivation_and_responses_v3.html#{section}"


def _flow_node(ax: plt.Axes, node_id: str, x: float, y: float,
               label: str, category: str, section: str,
               width: float = 2.45, height: float = 0.78,
               edge_override: str | None = None) -> dict[str, object]:
    """Draw one linked registry-category node and return its routing geometry."""
    style = PARAMETER_CATEGORY_STYLE[category]
    is_output = category == "output"
    boxstyle = ("round,pad=0.06,rounding_size=0.22"
                if is_output else "round,pad=0.04,rounding_size=0.045")
    patch = FancyBboxPatch(
        (x - width / 2.0, y - height / 2.0), width, height,
        boxstyle=boxstyle, facecolor=style["face"],
        edgecolor=edge_override or style["edge"],
        linewidth=1.45 if is_output else 1.15, zorder=4)
    url = _flow_section_url(section)
    patch.set_gid(f"flow-node-{node_id}")
    patch.set_url(url)
    ax.add_patch(patch)
    text_artist = ax.text(
        x, y, label, ha="center", va="center", fontsize=7.6,
        linespacing=1.22, color="#27323a", zorder=5)
    text_artist.set_url(url)
    return {
        "id": node_id, "x": x, "y": y, "width": width,
        "height": height, "category": category, "patch": patch,
    }


def _flow_edge(ax: plt.Axes, source: dict[str, object],
               target: dict[str, object], label: str = "",
               dashed: bool = False, rad: float = 0.0,
               label_offset: tuple[float, float] = (0.0, 0.0),
               color: str | None = None) -> None:
    """Route a dependency edge with the chart's two-style grammar."""
    edge_color = color or ("#9a6a00" if dashed else "#5d6a73")
    line_style = (0, (5, 4)) if dashed else "-"
    arrow = FancyArrowPatch(
        (float(source["x"]), float(source["y"])),
        (float(target["x"]), float(target["y"])),
        arrowstyle="-|>", mutation_scale=9.5, linewidth=1.15,
        linestyle=line_style, color=edge_color,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=31.0, shrinkB=31.0, zorder=2)
    ax.add_patch(arrow)
    if label:
        mid_x = (float(source["x"]) + float(target["x"])) / 2.0 + label_offset[0]
        mid_y = (float(source["y"]) + float(target["y"])) / 2.0 + label_offset[1]
        ax.text(
            mid_x, mid_y, label, ha="center", va="center",
            fontsize=6.4, color=edge_color, zorder=6,
            bbox={"boxstyle": "round,pad=0.14", "facecolor": "white",
                  "edgecolor": "none", "alpha": 0.93})


def _flow_bands(ax: plt.Axes,
                bands: tuple[tuple[str, float, float], ...]) -> None:
    """Draw subtle horizontal bands without introducing another encoding."""
    for index, (label, y0, y1) in enumerate(bands):
        ax.add_patch(Rectangle(
            (0.28, y0), 17.44, y1 - y0,
            facecolor="#f7f9fa" if index % 2 == 0 else "#ffffff",
            edgecolor="none", zorder=0))
        ax.plot([0.28, 17.72], [y0, y0],
                color="#dce2e6", linewidth=0.7, zorder=1)
        ax.text(0.42, y1 - 0.20, label, ha="left", va="top",
                fontsize=8.0, fontweight="bold", color="#69757e")


def _flow_category_key(ax: plt.Axes, y: float = 0.20) -> None:
    ax.text(0.58, y + 0.11, "registry category", ha="left", va="center",
            fontsize=7.2, fontweight="bold", color="#59636d")
    x_cursor = 2.10
    for category in ("input", "assumed", "derived", "output"):
        style = PARAMETER_CATEGORY_STYLE[category]
        rounding = 0.10 if category == "output" else 0.02
        patch = FancyBboxPatch(
            (x_cursor, y), 0.38, 0.22,
            boxstyle=f"round,pad=0.01,rounding_size={rounding}",
            facecolor=style["face"], edgecolor=style["edge"],
            linewidth=0.8)
        ax.add_patch(patch)
        ax.text(x_cursor + 0.46, y + 0.11, category,
                ha="left", va="center", fontsize=6.8, color="#59636d")
        x_cursor += 1.75
    ax.plot([9.45, 10.10], [y + 0.11, y + 0.11],
            color="#5d6a73", linewidth=1.2)
    ax.text(10.20, y + 0.11, "dependency", va="center",
            fontsize=6.8, color="#59636d")
    ax.plot([11.65, 12.30], [y + 0.11, y + 0.11],
            color="#9a6a00", linewidth=1.2, linestyle=(0, (5, 4)))
    ax.text(12.40, y + 0.11, "closure / calibration back edge",
            va="center", fontsize=6.8, color="#59636d")
    ax.text(17.42, y + 0.11, "click any node to open its derivation",
            ha="right", va="center", fontsize=6.8, color="#59636d")


def plot_flowchart_provenance_structure() -> Path:
    """Render standalone Chart A: provenance, scalar closure, and reduction."""
    validate_parameter_registry()
    fig, ax = plt.subplots(figsize=(17.2, 10.2))
    ax.axis("off")
    ax.set_xlim(0.0, 18.0)
    ax.set_ylim(0.0, 11.6)
    _flow_bands(ax, (
        ("BAND 1 — SOURCES", 9.25, 10.75),
        ("BAND 2 — DERIVED SCALARS AND CLOSURES", 5.95, 9.25),
        ("BAND 3 — REDUCTION TESTS AND ARTIFACTS", 0.62, 5.95),
    ))
    ax.text(9.0, 11.32, "Chart A — Parameter provenance and model structure",
            ha="center", fontsize=15, fontweight="bold", color="#24313a")
    ax.text(
        9.0, 11.02,
        "Forward dependencies are solid. The only dashed edges are the two declared closure/calibration loops.",
        ha="center", fontsize=8.0, color="#61707a")

    nodes: dict[str, dict[str, object]] = {}
    def add(node_id: str, *args: object, **kwargs: object) -> None:
        nodes[node_id] = _flow_node(ax, node_id, *args, **kwargs)

    # Band 1: sources.
    add("source_motor", 2.35, 9.90,
        "Motor + coupling datasheets\n$T_{max},\\ J_m,\\ m_c,\\ k_c$",
        "input", "2-entry-parameters", width=3.10)
    add("source_screw", 6.70, 9.90,
        "Screw geometry\n$L,\\ L_s,\\ d_s,\\rho_s$",
        "input", "2-entry-parameters", width=2.95)
    add("source_driver", 11.00, 9.90,
        "Driver configuration\ncurrent, mode, microstep divisor",
        "assumed", "5-stepper-input-nonlinear-law-linearization-and-bound",
        width=3.15)
    add("source_modal", 15.45, 9.90,
        "Modal measurement\n$155$–$190$ Hz and $\\approx690$ Hz",
        "input", "5-3-can-the-present-model-reproduce-every-measured-feature",
        width=3.20)

    # Band 2: the two explicit convergences and executable scalars.
    add("ratio", 1.65, 8.15, "$r=L/(2\\pi)$",
        "derived", "5-stepper-input-nonlinear-law-linearization-and-bound",
        width=1.85, height=0.66)
    add("inertia_bundle", 4.35, 8.15,
        "Component inertias\n$J_m$ • $J_c$ • $J_{s,drive}$ • $J_{s,tail}$",
        "assumed", "2-entry-parameters", width=3.30)
    add("j_sum", 7.42, 8.15, "$J_\\Sigma=\\sum J_i$\ninformation collapsed",
        "derived", "6-reduction-from-ten-dofs-to-two", width=2.35)
    add("m_drive", 10.15, 8.15, "$m_d=J_\\Sigma/r^2$",
        "derived", "6-reduction-from-ten-dofs-to-two", width=2.20)
    add("k_mag", 13.05, 8.15, "$K_m=N_rT_{max}/r^2$",
        "derived", "5-stepper-input-nonlinear-law-linearization-and-bound",
        width=2.55)
    add("k_det", 16.15, 8.15,
        "$K_{det}(x_0)$\nlocal tangent only",
        "derived", "5-stepper-input-nonlinear-law-linearization-and-bound",
        width=2.45)
    add("series_bundle", 3.00, 6.72,
        "Independent spring inputs\n$k_{brg}$ • $k_{sha}$ • $k_{mnt}$",
        "assumed", "4-full-ten-dof-derivation", width=3.45)
    add("k_ball", 6.35, 6.72,
        "$k_{ball}$ remainder\nclosure-derived, not measured",
        "derived", "4-full-ten-dof-derivation", width=2.75)
    add("k_ax", 9.40, 6.72,
        "$k_{ax}=\\left(\\sum 1/k_i\\right)^{-1}$\nmodal-calibrated target",
        "derived", "6-reduction-from-ten-dofs-to-two", width=2.85)
    add("full_plant", 13.05, 6.72,
        "Ten-DOF plant\n$\\mathbf{M},\\mathbf{C},\\mathbf{K},\\mathbf{b}$",
        "derived", "4-full-ten-dof-derivation", width=2.90)

    # Band 3: three sequential tests, retained plant, rejected dead end.
    add("test_ratio", 13.05, 5.20,
        "TEST 1 — reflected stiffness\n$173\\times$ and $731\\times$\n"
        "torsional coordinates rigid → collapse",
        "derived", "6-reduction-from-ten-dofs-to-two",
        width=3.60, height=0.94)
    add("test_load", 13.05, 3.92,
        "TEST 2 — load path\n$u_f$ and $\\theta_{s3}$ are stubs → drop",
        "derived", "6-reduction-from-ten-dofs-to-two",
        width=3.60, height=0.84)
    add("test_port", 13.05, 2.68,
        "TEST 3 — port survival\n1 DOF gives $v_n=0$ → reject",
        "derived", "6-reduction-from-ten-dofs-to-two",
        width=3.60, height=0.84)
    add("mode_compare", 4.30, 1.30,
        "Modal comparison\nmeasured $690$ Hz → modeled $696$ Hz",
        "output", "7-full-versus-reduced-verification", width=3.35)
    add("reduced_plant", 9.25, 1.30,
        "Retained two-DOF plant\n$x_d,\\ x_s$; relative port survives",
        "derived", "6-reduction-from-ten-dofs-to-two", width=3.25)
    add("one_dof", 15.80, 1.30,
        "TERMINATOR — one DOF rejected\nrelative mode and $F_{f,n}$ lost",
        "output", "6-reduction-from-ten-dofs-to-two",
        width=3.25, edge_override="#b23a48")

    if len(nodes) > 22:
        raise ValueError(f"Chart A exceeds the 22-node cap: {len(nodes)}")

    # Registry-derived dependencies.
    _flow_edge(ax, nodes["source_screw"], nodes["ratio"], "$L$")
    _flow_edge(ax, nodes["source_motor"], nodes["inertia_bundle"], "$J_m,J_c$")
    _flow_edge(ax, nodes["source_screw"], nodes["inertia_bundle"], "screw $J_i$")
    _flow_edge(ax, nodes["inertia_bundle"], nodes["j_sum"], "four components")
    _flow_edge(ax, nodes["j_sum"], nodes["m_drive"])
    _flow_edge(ax, nodes["ratio"], nodes["m_drive"], "$r^{-2}$", rad=-0.10)
    _flow_edge(ax, nodes["source_motor"], nodes["k_mag"], "$T_{max}$", rad=-0.12)
    _flow_edge(ax, nodes["source_driver"], nodes["k_mag"], "current", rad=0.10)
    _flow_edge(ax, nodes["ratio"], nodes["k_mag"], "$r^{-2}$", rad=-0.16)
    _flow_edge(ax, nodes["source_motor"], nodes["k_det"], "$N_r,\\hat T_{det}$", rad=-0.18)
    _flow_edge(ax, nodes["source_driver"], nodes["k_det"], "$\\phi_{det}$", rad=0.12)
    _flow_edge(ax, nodes["series_bundle"], nodes["k_ball"])
    _flow_edge(ax, nodes["k_ball"], nodes["k_ax"])

    # Required dashed closure/calibration back edges: exactly two.
    _flow_edge(
        ax, nodes["k_ax"], nodes["k_ball"],
        "closure, not independent validation", dashed=True,
        rad=0.34, label_offset=(0.0, 0.50))
    _flow_edge(
        ax, nodes["source_modal"], nodes["k_ax"],
        "closure, not independent validation", dashed=True,
        rad=-0.40, label_offset=(1.10, 0.15))

    # Plant assembly and sequential reduction decisions.
    for source_id, label, rad in (
        ("m_drive", "$m_d$", 0.04),
        ("k_mag", "$K_m$", 0.00),
        ("k_det", "local sensitivity", -0.08),
        ("k_ax", "$k_{ax}$", 0.05),
    ):
        _flow_edge(ax, nodes[source_id], nodes["full_plant"], label, rad=rad)
    _flow_edge(ax, nodes["full_plant"], nodes["test_ratio"])
    _flow_edge(ax, nodes["test_ratio"], nodes["test_load"])
    _flow_edge(ax, nodes["test_load"], nodes["test_port"])
    _flow_edge(
        ax, nodes["test_port"], nodes["reduced_plant"],
        "Masses → $m_d,m_s$", rad=0.16,
        label_offset=(-0.15, 0.22), color="#52778b")
    _flow_edge(
        ax, nodes["test_port"], nodes["reduced_plant"],
        "Springs → $k_{ax}$", rad=-0.16,
        label_offset=(-0.15, -0.22), color="#b06b00")
    _flow_edge(ax, nodes["test_port"], nodes["one_dof"],
               "dead-end branch", rad=-0.10, color="#b23a48")
    _flow_edge(ax, nodes["reduced_plant"], nodes["mode_compare"],
               "predicts $696$ Hz", rad=0.08)
    _flow_edge(ax, nodes["source_modal"], nodes["mode_compare"],
               "measured $690$ Hz", rad=0.28)

    _flow_category_key(ax)
    chart_link = ax.text(
        17.35, 0.57, "continues in Chart B →", ha="right",
        fontsize=7.4, fontweight="bold", color="#1f6f8b")
    chart_link.set_url("flowchart_B_friction_and_results.svg")
    ax.text(17.68, 10.57, f"{len(nodes)} nodes", ha="right",
            fontsize=6.8, color="#7a858d")
    fig.subplots_adjust(left=0.015, right=0.985, top=0.985, bottom=0.02)
    output = ASSET_DIR / "flowchart_A_provenance_and_structure.svg"
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_flowchart_friction_results() -> Path:
    """Render standalone Chart B: friction ports, paths, and valid outputs."""
    validate_parameter_registry()
    fig, ax = plt.subplots(figsize=(17.2, 10.8))
    ax.axis("off")
    ax.set_xlim(0.0, 18.0)
    ax.set_ylim(0.0, 11.8)
    _flow_bands(ax, (
        ("BAND 1 — PROVISIONAL FRICTION PARAMETERS", 8.55, 10.05),
        ("BAND 2 — PHYSICAL POWER PORTS", 6.60, 8.55),
        ("BAND 3 — LINEAR AND NONLINEAR ENTRY PATHS", 4.75, 6.60),
        ("BAND 4 — MECHANISMS THAT ALTER THE RESPONSE", 2.55, 4.75),
        ("BAND 5 — REPORTED OUTPUTS AND VALIDITY CONDITIONS", 0.62, 2.55),
    ))
    ax.text(9.0, 11.48, "Chart B — Friction insertion and result provenance",
            ha="center", fontsize=15, fontweight="bold", color="#24313a")
    ax.text(
        9.0, 11.16,
        "The same physical port feeds either a local tangent model or time-domain internal states; it is never added as an external correction.",
        ha="center", fontsize=8.0, color="#61707a")

    nodes: dict[str, dict[str, object]] = {}

    def add(node_id: str, *args: object, **kwargs: object) -> None:
        nodes[node_id] = _flow_node(ax, node_id, *args, **kwargs)

    add(
        "reduced_entry", 9.0, 10.55,
        "FROM CHART A — retained two-DOF plant\n"
        "$\\mathbf{M}_r,\\mathbf{C}_r,\\mathbf{K}_r,\\mathbf{b}_r$; $x_d,x_s$",
        "derived", "6-reduction-from-ten-dofs-to-two",
        width=4.35, height=0.76)

    # Band 1: parameter groups are amber because the executed values are
    # preemptive/provisional rather than independently identified measurements.
    add(
        "param_g", 2.35, 9.25,
        "Guideway set g\n$F_c,F_s,v_s,\\sigma_0,\\sigma_1,\\sigma_2$",
        "assumed", "8-1-executed-provisional-friction-values",
        width=3.25)
    add(
        "param_n", 6.75, 9.25,
        "Nut set n\nincludes $\\sigma_{0,n}$ and GMS $k_i,\\nu_i,C$",
        "assumed", "8-1-executed-provisional-friction-values",
        width=3.35)
    add(
        "param_r", 11.25, 9.25,
        "Rolling allocation r\nprovisional gross drive-side share",
        "assumed", "8-1-executed-provisional-friction-values",
        width=3.25)
    add(
        "param_d", 15.65, 9.25,
        "Drive allocation d\nprovisional gross drive-side share",
        "assumed", "8-1-executed-provisional-friction-values",
        width=3.25)

    # Band 2: physical ports. r and d intentionally have the same incidence row.
    add(
        "port_g", 2.35, 7.42,
        "PORT g — guideway to ground\n$v_g=\\dot x_s$\nA/A2 and C/C2",
        "derived", "8-2-force-locations", width=3.25, height=0.90)
    add(
        "port_n", 6.75, 7.42,
        "PORT n — nut internal pair\n$v_n=\\dot x_d-\\dot x_s$\nB/B2 and C/C2",
        "derived", "8-2-force-locations", width=3.35, height=0.90)
    add(
        "port_r", 11.25, 7.42,
        "PORT r — drive side to ground\n$v_r=\\dot x_d$; $H_r=[1,0]$\nB/B2 and C/C2",
        "derived", "8-2-force-locations", width=3.35, height=0.90)
    add(
        "port_d", 15.65, 7.42,
        "PORT d — drive side to ground\n$v_d=\\dot x_d$; $H_d=[1,0]$\nall friction cases",
        "derived", "8-2-force-locations", width=3.35, height=0.90)

    # Band 3: the two legitimate insertion paths.
    add(
        "linear_path", 6.10, 5.45,
        "LINEAR TANGENT PATH\n"
        "$\\Delta\\mathbf{K}=\\sigma_0\\mathbf{H}^T\\mathbf{H}$\n"
        "local presliding role only",
        "derived", "8-friction-constitutive-laws",
        width=3.75, height=0.98)
    add(
        "nonlinear_path", 11.90, 5.45,
        "NONLINEAR STATE PATH\n"
        "LuGre $z$ or GMS $F_i$ integrated\nwith the mechanical states",
        "derived", "8-friction-constitutive-laws",
        width=3.75, height=0.98)

    # Band 4: mechanism-level interpretation.
    add(
        "mech_tangent", 2.35, 3.58,
        "Differential tangent stiffening\n"
        "$k_{ax}+\\sigma_{0,n}H_n^TH_n$\nmode shift 696 → 775 Hz",
        "derived", "8-friction-constitutive-laws",
        width=3.35, height=1.00)
    add(
        "mech_preslide", 6.75, 3.58,
        "Presliding compliance\n"
        "grounded ports change tangent gain\nbefore first yield",
        "derived", "9-force-instrumented-partial-slip-memory-experiment",
        width=3.30, height=1.00)
    add(
        "mech_memory", 11.25, 3.58,
        "Yield and memory\n"
        "LuGre $z$ / GMS $F_i$\nreturn-point and final-origin closure",
        "derived", "9-force-instrumented-partial-slip-memory-experiment",
        width=3.35, height=1.00)
    add(
        "mech_drag", 15.65, 3.58,
        "Identifiable gross drive drag\n"
        "r and d are one lump at $\\dot x_d$\nsteady position deviation",
        "derived", "8-2-force-locations",
        width=3.35, height=1.00)

    # Band 5: outputs are pill-shaped and include the condition under which
    # the number is meaningful.
    add(
        "out_modes", 2.35, 1.48,
        "Mode table\nvalid at stated operating point\nand local detent tangent",
        "output", "11-generated-numerical-summary",
        width=3.25, height=0.86)
    add(
        "out_gain", 6.75, 1.48,
        "Tangent gain\nvalid before first yield\n—not finite travel tracking—",
        "output", "10-response-comparison-across-friction-cases",
        width=3.25, height=0.86)
    add(
        "out_time", 11.25, 1.48,
        "Time-domain deviation metrics\nvalid only after derived settling dwell\nand integration convergence",
        "output", "13-verification-checks-and-limitations",
        width=3.45, height=0.86)
    add(
        "out_residual", 15.65, 1.48,
        "Reduction residual\nnormalized by command amplitude\nwithin retained bandwidth",
        "output", "7-full-versus-reduced-verification",
        width=3.25, height=0.86)

    if len(nodes) > 22:
        raise ValueError(f"Chart B exceeds the 22-node cap: {len(nodes)}")

    # Parameter group to physical port.
    for suffix in ("g", "n", "r", "d"):
        _flow_edge(ax, nodes[f"param_{suffix}"], nodes[f"port_{suffix}"])

    # A single incidence bus keeps the fork readable: every physical port and
    # the retained mechanics enter the same junction. r and d therefore have
    # visibly identical downstream connectivity, not merely similar labels.
    bus_y = 6.28
    ax.plot([2.35, 15.65], [bus_y, bus_y], color="#5d6a73",
            linewidth=1.35, zorder=2)
    for suffix in ("g", "n", "r", "d"):
        x = float(nodes[f"port_{suffix}"]["x"])
        ax.add_patch(FancyArrowPatch(
            (x, 6.97), (x, bus_y),
            arrowstyle="-|>", mutation_scale=8.5, linewidth=1.1,
            color="#5d6a73", shrinkA=0.0, shrinkB=0.0, zorder=3))
    ax.add_patch(FancyArrowPatch(
        (9.0, 10.17), (9.0, bus_y),
        arrowstyle="-|>", mutation_scale=8.5, linewidth=1.15,
        color="#5d6a73", shrinkA=0.0, shrinkB=0.0, zorder=2))
    ax.text(9.13, 8.62, "retained mechanics: $k_{ax}$",
            ha="left", va="center", fontsize=6.5, color="#5d6a73",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white",
                  "edgecolor": "none", "alpha": 0.94}, zorder=6)
    ax.text(6.75, 8.78, "$\\sigma_{0,n}$",
            ha="center", va="center", fontsize=6.5, color="#9a6a00",
            bbox={"boxstyle": "round,pad=0.10", "facecolor": "white",
                  "edgecolor": "none", "alpha": 0.94}, zorder=6)
    ax.text(9.0, bus_y + 0.12, "port-incidence bus $\\mathbf{H}$",
            ha="center", va="bottom", fontsize=6.6, color="#5d6a73",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white",
                  "edgecolor": "none", "alpha": 0.95}, zorder=6)
    for target_id in ("linear_path", "nonlinear_path"):
        target = nodes[target_id]
        x = float(target["x"])
        target_top = float(target["y"]) + float(target["height"]) / 2.0
        ax.add_patch(FancyArrowPatch(
            (x, bus_y), (x, target_top),
            arrowstyle="-|>", mutation_scale=8.5, linewidth=1.15,
            color="#5d6a73", shrinkA=0.0, shrinkB=0.0, zorder=3))

    # The displayed tangent expression is the explicit sigma0,n versus k_ax
    # convergence. Finite-amplitude reversal data are needed to separate them.
    _flow_edge(
        ax, nodes["linear_path"], nodes["mech_tangent"],
        "$k_{ax}$ + $\\sigma_{0,n}H_n^TH_n$", rad=0.13,
        label_offset=(-0.25, 0.08))
    _flow_edge(ax, nodes["linear_path"], nodes["mech_preslide"], rad=-0.06)
    _flow_edge(ax, nodes["nonlinear_path"], nodes["mech_memory"], rad=0.06)
    _flow_edge(ax, nodes["nonlinear_path"], nodes["mech_drag"], rad=-0.16)

    _flow_edge(ax, nodes["mech_tangent"], nodes["out_modes"])
    _flow_edge(ax, nodes["mech_preslide"], nodes["out_gain"])
    _flow_edge(ax, nodes["mech_memory"], nodes["out_time"])
    _flow_edge(ax, nodes["mech_drag"], nodes["out_time"], rad=0.20)
    # Reduction residual bypasses the friction mechanisms. Route it around the
    # right margin so that independence is visible without crossing the graph.
    ax.plot([11.18, 17.62, 17.62], [10.55, 10.55, 1.48],
            color="#5d6a73", linewidth=1.1, zorder=2)
    ax.add_patch(FancyArrowPatch(
        (17.62, 1.48), (17.28, 1.48),
        arrowstyle="-|>", mutation_scale=8.5, linewidth=1.1,
        color="#5d6a73", shrinkA=0.0, shrinkB=0.0, zorder=3))
    ax.text(17.49, 5.25, "full/reduced comparison",
            ha="center", va="center", rotation=90,
            fontsize=6.3, color="#5d6a73",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white",
                  "edgecolor": "none", "alpha": 0.94}, zorder=6)

    # Tiny topology glyphs make grounding versus an internal equal/opposite port
    # visible without adding nodes or a second edge grammar.
    for suffix in ("g", "r", "d"):
        port = nodes[f"port_{suffix}"]
        x = float(port["x"]) + 1.16
        y_ground = 7.08
        ax.plot([x, x], [7.27, y_ground], color="#66727a",
                linewidth=0.95, zorder=5)
        ax.plot([x - 0.19, x + 0.19], [y_ground, y_ground],
                color="#66727a", linewidth=1.05, zorder=5)
        for hatch_x in np.linspace(x - 0.15, x + 0.15, 4):
            ax.plot([hatch_x - 0.05, hatch_x + 0.01],
                    [y_ground - 0.07, y_ground],
                    color="#8a949b", linewidth=0.70, zorder=5)
    internal_x = float(nodes["port_n"]["x"]) + 1.13
    internal_y = 7.10
    ax.add_patch(Rectangle(
        (internal_x - 0.34, internal_y - 0.08), 0.18, 0.16,
        facecolor="#dceef6", edgecolor="#55798a", linewidth=0.75, zorder=5))
    ax.add_patch(Rectangle(
        (internal_x + 0.16, internal_y - 0.08), 0.18, 0.16,
        facecolor="#dff2ea", edgecolor="#557f70", linewidth=0.75, zorder=5))
    ax.add_patch(FancyArrowPatch(
        (internal_x - 0.02, internal_y), (internal_x - 0.15, internal_y),
        arrowstyle="-|>", mutation_scale=7, linewidth=0.8,
        color="#9b4d52", zorder=6))
    ax.add_patch(FancyArrowPatch(
        (internal_x + 0.02, internal_y), (internal_x + 0.15, internal_y),
        arrowstyle="-|>", mutation_scale=7, linewidth=0.8,
        color="#9b4d52", zorder=6))

    _flow_category_key(ax)
    chart_link = ax.text(
        0.62, 0.57, "← back to Chart A", ha="left",
        fontsize=7.4, fontweight="bold", color="#1f6f8b")
    chart_link.set_url("flowchart_A_provenance_and_structure.svg")
    ax.text(17.68, 10.02, f"{len(nodes)} nodes", ha="right",
            fontsize=6.8, color="#7a858d")
    fig.subplots_adjust(left=0.015, right=0.985, top=0.985, bottom=0.02)
    output = ASSET_DIR / "flowchart_B_friction_and_results.svg"
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

    element(2.00, 2.25, r"R: $c_m,F_{f,d}$", "#eee8f3", 2.08)
    bond((3.03, 3.38), (2.36, 2.52), "#6a4c93")
    element(11.18, 2.25, r"R: $F_{f,g}$", "#f8dce1", 1.45)
    bond((10.17, 3.38), (10.82, 2.52), "#b23a48")

    ax.text(6.60, 4.30, r"internal port: $v_n=\dot x_d-\dot x_s$",
            ha="center", fontsize=8.5, color="#8d2936")
    ax.text(4.70, 3.92, r"$-F_{f,n}$", fontsize=8.0, color="#8d2936")
    ax.text(8.05, 3.92, r"$+F_{f,n}$", fontsize=8.0, color="#8d2936")
    ax.text(6.60, 1.48,
            r"$\mathbf{H}_g=[0,1],\quad \mathbf{H}_n=[1,-1],\quad \mathbf{H}_d=[1,0]$",
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
    ax_residual.set_title(
        f"Reduction residual: RMS {verification['rms_residual_pct_command']:.2f}% command; "
        f"peak {verification['peak_residual_pct_command']:.2f}% command")
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
    # A 150 mm stage stroke within about 170 mm usable screw travel implies a
    # 20 mm minimum support-to-nut free length in this illustrative datum.
    positions = np.array([0.0, 75.0, FULL["stage_travel"] * 1e3])
    free_lengths = (FULL["usable_screw_travel"] - FULL["stage_travel"]
                    + positions * 1e-3)
    axial_rigidity = FULL["k_sha"] * 0.150
    k_sha = axial_rigidity / free_lengths
    constants = physical_constants()
    fixed_compliance = 1.0 / constants["k_ax"] - 1.0 / FULL["k_sha"]
    k_ax = 1.0 / (fixed_compliance + 1.0 / k_sha)
    mass = np.diag([constants["m_d"], constants["m_s"]])
    mode = []
    for stiffness_value in k_ax:
        stiffness = np.array([
            [constants["K_m"] + stiffness_value, -stiffness_value],
            [-stiffness_value, stiffness_value],
        ])
        mode.append(_linear_modes(mass, stiffness)[1])
    mode = np.asarray(mode)
    fig, left = plt.subplots(figsize=(10.2, 4.6))
    right = left.twinx()
    left.plot(positions, k_ax / 1e6, marker="o", color="#277da1", linewidth=2,
              label="$k_{ax}$")
    right.plot(positions, mode, marker="s", color="#d97800", linewidth=2,
               label="predicted stage mode")
    for x, free_length, sha in zip(positions, free_lengths, k_sha):
        left.annotate(f"$L_{{free}}$={free_length * 1e3:.0f} mm\n$k_{{sha}}$={sha / 1e6:.0f} MN/m",
                      (x, np.interp(x, positions, k_ax / 1e6)),
                      xytext=(0, 10), textcoords="offset points", ha="center", fontsize=8)
    left.set_xlabel("Stage position across 150 mm travel (mm)")
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
    local_low, local_high = detent_local_mode_band()
    damping_ratios = (0.02, 0.05, MODEL["zeta_m"], 0.50)
    colors = ("#b23a48", "#277da1", "#d97800", "#7a7a7a")
    fig, (stage_ax, drive_ax) = plt.subplots(1, 2, figsize=(11.4, 4.7), sharex=True)
    for zeta, color in zip(damping_ratios, colors):
        c_m = 2.0 * zeta * np.sqrt(constants["K_m"] * constants["m_d"])
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
        axis.axvspan(local_low, local_high, color="#d97800", alpha=0.12,
                     label=f"local detent-tangent sweep {local_low:.0f}-{local_high:.0f} Hz")
        axis.axvline(low_mode, color="#252525", linestyle="--", linewidth=1.1,
                    label=f"global commutation pole {low_mode:.1f} Hz")
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
             f"Global Bode curves exclude detent as an origin spring; the shaded band uses its local tangent extremes. "
             f"The nonlinear model retains the full {constants['detent_period'] * 1e6:.1f} um-period torque.",
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
    constants = physical_constants()
    local_low, local_high = detent_local_mode_band()
    lines = [
        "<!-- BEGIN GENERATED RESPONSE SUMMARY -->",
        "| Case | Friction law | Global-linear modes (Hz) | Local friction-tangent gain $X_s/X_{cmd}$ | Smallest first-yield travel | First-step overshoot | Settled-window RMS deviation | Settled-window maximum | All-time peak deviation | Final-window RMS deviation |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
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
            f"{time_metrics[key]['rms_settled_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['max_settled_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['max_abs_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['rms_final_error_nm']:.1f} nm |"
        )
    lines.extend([
        "",
        "The displayed modes and gains are the global commutation linearization: periodic detent is deliberately excluded from the global stiffness matrix. The friction tangent is local and valid only below the listed first-yield travel. "
        f"The nonlinear cases include the periodic detent torque and use a {constants['plateau_dwell'] * 1e3:.0f} ms dwell derived from the 2% settling estimate ({constants['settling_time_2pct'] * 1e3:.1f} ms, with a 100 ms floor). "
        f"Settled values collect the last {constants['metric_window'] * 1e3:.0f} ms of every plateau. All deviation columns use $d(t)=x_{{cmd}}(t)-x_s(t)$ and describe open-loop modeled plant behavior, not servo tracking performance. Case 0 remains frictionless.",
        "",
        "### Generated reduction audit",
        "",
        "| Quantity | Executed value |",
        "|---|---:|",
        f"| Measured stage body mass | {constants['m_stage']:.3f} kg |",
        f"| Nut body mass retained at stage node | {constants['m_n']:.3f} kg |",
        f"| Derived retained stage-side mass | [[derived:reduced_stage_mass={constants['m_s']:.3f}]] kg |",
        f"| Upper-mode calibration target | {constants['axial_mode_target_hz']:.2f} Hz |",
        f"| Modal-calibrated $k_{{ax}}$ | [[derived:reduced_axial_stiffness={constants['k_ax']:.6e}]] N/m |",
        f"| Closure-derived $k_{{ball}}$ | [[derived:k_ball={verification['parameters']['k_ball']:.6e}]] N/m |",
        f"| Motor rotor inertia | {verification['parameters']['J_m']:.3e} kg m² |",
        f"| Coupling inertia | {verification['parameters']['J_c']:.3e} kg m² |",
        f"| 0.192 m screw inertia | {verification['parameters']['screw_inertia']:.3e} kg m² |",
        f"| 0.192 m screw mass | {verification['parameters']['screw_mass']:.4f} kg |",
        f"| Stage travel / usable screw distance | {FULL['stage_travel'] * 1e3:.0f} / {FULL['usable_screw_travel'] * 1e3:.0f} mm |",
        f"| Full-model reflected drivetrain mass | {verification['parameters']['m_d_reflected']:.3f} kg |",
        f"| Rated-current holding torque | {MODEL['T_max']:.3f} N m |",
        f"| Enabled detent torque | {MODEL['T_det']:.3f} N m |",
        f"| Global commutation low pole | {linear_metrics['0']['modes'][0]:.2f} Hz |",
        f"| Local detent-tangent low-pole band | {local_low:.2f} to {local_high:.2f} Hz |",
        f"| Full/reduced sequence RMS residual | {verification['rms_residual_nm']:.3f} nm |",
        f"| Full/reduced sequence peak residual | {verification['peak_residual_nm']:.3f} nm |",
        f"| RMS residual / command amplitude | {verification['rms_residual_pct_command']:.3f}% |",
        f"| Peak residual / command amplitude | {verification['peak_residual_pct_command']:.3f}% |",
        "",
        "The reduced drive mass is derived from the listed component inertias and the current lead. It is not an independent input. The normalized residual, unlike its nanometre value, is invariant to a simple rescaling of this linear verification command.",
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
        ("B/B2", "B", "B2", "Nut microslip shifts the relative mode; the same lumped drive tangent is shared by every friction case"),
        ("C/C2", "C", "C2", "All three identifiable friction tangents are active"),
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


def generated_presliding_summary(experiments: dict[str, dict[str, object]]) -> str:
    lines = [
        "<!-- BEGIN GENERATED PRESLIDING SUMMARY -->",
    ]
    rows = (
        ("Whole-sequence RMS command-output deviation", "whole_rms_nm", "nm", False),
        ("Peak absolute command-output deviation", "max_abs_deviation_nm", "nm", False),
        ("Mean repeated-return deviation mismatch", "return_error_mismatch_nm", "nm", False),
        ("Mean repeated-return friction-force mismatch", "return_force_mismatch_N", "N", False),
        ("Absolute mean error after final return to zero", "final_mean_nm", "nm", True),
    )
    for experiment_name, experiment in experiments.items():
        metrics = experiment["metrics"]
        lugre_key, gms_key = experiment["keys"]
        site = experiment["site"]
        site_title = "Guideway" if site == "g" else "Nut microslip"
        boundary_note = (
            "Normal free-stage plant; the observed output is the stage coordinate."
            if not experiment["blocked_stage"] else
            "Dedicated blocked-stage identification boundary, $x_s=0$; the observed output is the drive coordinate."
        )
        lines.extend([
            f"### {site_title}: {lugre_key}/{gms_key}",
            "",
            boundary_note,
            "",
            f"| Executed metric | LuGre {lugre_key} | GMS {gms_key} | GMS minus LuGre |",
            "|---|---:|---:|---:|",
        ])
        for label, key, unit, use_absolute in rows:
            lugre = float(metrics[lugre_key][key])
            gms = float(metrics[gms_key][key])
            if use_absolute:
                lugre, gms = abs(lugre), abs(gms)
            precision = 4 if unit == "N" else 2
            lines.append(
                f"| {label} | {lugre:.{precision}f} {unit} | {gms:.{precision}f} {unit} | "
                f"{gms - lugre:+.{precision}f} {unit} |"
            )
        max_force = max(float(metrics[key]["max_force_N"])
                        for key in (lugre_key, gms_key))
        lines.extend([
            "",
            f"Maximum executed {site_title.lower()} friction is **{max_force:.3f} N** "
            f"({100.0 * max_force / FRICTION[site]['F_s']:.1f}% of the provisional "
            f"{FRICTION[site]['F_s']:.1f} N macro breakaway level).",
            "",
        ])
    first_experiment = next(iter(experiments.values()))
    lines.extend([
        f"Every plateau is held for **{first_experiment['plateau_dwell'] * 1e3:.0f} ms**, so return-point means are settled samples rather than drive-mode ringing. "
        r"The dedicated blocked-stage B/B2 force-deflection loop is the finite-amplitude test of the exact small-signal correlation between $k_{ax}$ and $\sigma_{0,n}$: both multiply $[1,-1]^T[1,-1]$ before microslip yields. This identification-fixture boundary is not used for the normal plant-response plots.",
        "",
        "The whole-sequence RMS still includes instantaneous command edges. Repeated-return force and final-origin measures target constitutive history. The provisional parameters do not predetermine that GMS is better; measured loops must select and fit the law.",
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
    dt_us = tuple(value * 1e6 for value in GMS_CONVERGENCE_DTS)
    constants = physical_constants()
    lines = [
        "<!-- BEGIN GENERATED STEP HALVING SUMMARY -->",
        f"| Case | {dt_us[0]:.1f} us | {dt_us[1]:.1f} us | {dt_us[2]:.1f} us | "
        f"$\\Delta R_{{{dt_us[0]:g}\\to{dt_us[1]:g}}}$ | "
        f"$\\Delta R_{{{dt_us[1]:g}\\to{dt_us[2]:g}}}$ | Difference ratio |",
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
        interpretation + f" The largest {dt_us[1]:.1f}-to-{dt_us[2]:.1f} us relative change is **{max_relative:.4f}%**.",
        "",
        f"These values use the identical {main_duration(constants) * 1e3:.0f} ms zero-order-held, yield-spanning command and the identical final {constants['metric_window'] * 1e3:.0f} ms "
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


def _format_default_like(value: object, template: str) -> str:
    """Format an authoritative Python value with the Markdown token's style."""
    if isinstance(value, str):
        return value
    numeric = float(value)
    stripped = template.strip()
    if "e" in stripped.lower():
        mantissa = stripped.lower().split("e", 1)[0]
        digits = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
        return f"{numeric:.{digits}e}"
    if "." in stripped:
        digits = len(stripped.split(".", 1)[1])
        return f"{numeric:.{digits}f}"
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:.12g}"


def browser_parameter_defaults() -> dict[str, object]:
    """Single-source editable defaults emitted into both rendered documents."""
    p = component_parameters()
    defaults: dict[str, object] = {
        "lead": MODEL["lead"],
        "rotor_teeth": MODEL["rotor_teeth"],
        "holding_torque": MODEL["T_max"],
        "detent_torque": MODEL["T_det"],
        "detent_phase": MODEL["detent_phase"],
        "axial_mode_target_hz": MODEL["axial_mode_target_hz"],
        "axial_damping": MODEL["c_ax"],
        "electromagnetic_zeta": MODEL["zeta_m"],
        "microstep_divisor": MODEL["microstep_divisor"],
        "J_m": p["J_m"],
        "J_c": p["J_c"],
        "screw_length": p["screw_length"],
        "usable_screw_travel": p["usable_screw_travel"],
        "stage_travel": p["stage_travel"],
        "lead_accuracy_class": "IT1",
        "screw_diameter": p["screw_diameter"],
        "screw_density": p["screw_density"],
        "nut_mass": p["m_n"],
        "stage_mass": p["m_stage"],
        "k_c_series": p["k_c_series"],
        "k_theta_a": p["k_theta_a"],
        "k_theta_b": p["k_theta_b"],
        "k_brg": p["k_brg"],
        "k_sha": p["k_sha"],
        "k_shb": p["k_shb"],
        "k_mnt": p["k_mnt"],
        "zeta_internal": p["zeta_internal"],
    }
    for site, values in FRICTION.items():
        defaults.update({
            f"{site}_sigma0": values["sigma0"],
            f"{site}_sigma1": values["sigma1"],
            f"{site}_sigma2": values["sigma2"],
            f"{site}_Fs": values["F_s"],
            f"{site}_Fc": values["F_c"],
            f"{site}_vs": values["v_s"],
            f"{site}_C": values["C_gms"],
        })
        for index, fraction in enumerate(GMS_STIFFNESS_FRACTIONS, start=1):
            defaults[f"{site}_k{index}"] = fraction * values["sigma0"]
    for index, weight in enumerate(GMS_WEIGHTS, start=1):
        defaults[f"gms_nu{index}"] = weight
    return defaults


def browser_derived_defaults() -> dict[str, float]:
    """Derived defaults generated from the same equations as the simulations."""
    constants = physical_constants()
    component = component_parameters()
    mass, _damping, stiffness, _input = linear_matrices((), "none")
    modes = _linear_modes(mass, stiffness)
    return {
        "transmission_ratio": constants["r"],
        "total_rotational_inertia": component["J_total"],
        "reduced_drive_mass": constants["m_d"],
        "reduced_stage_mass": constants["m_s"],
        "magnetic_stiffness": constants["K_m"],
        "detent_stiffness": constants["K_det"],
        "reduced_axial_stiffness": constants["k_ax"],
        "k_ball": constants["k_ball"],
        "full_step_pitch": constants["full_step"],
        "quarter_step_bound": constants["quarter_step"],
        "command_step": constants["command_step"],
        "interpolated_step": constants["interpolated_step"],
        "screw_inertia": component["screw_inertia"],
        "screw_segment_inertia": component["screw_inertia"] / 3.0,
        "screw_mass": component["screw_mass"],
        "screw_segment_mass": component["screw_mass"] / 3.0,
        "k_c_half": component["k_c1"],
        "mode_1_hz": float(modes[0]),
        "mode_2_hz": float(modes[1]),
    }


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
        authoritative = browser_parameter_defaults().get(key, value.strip())
        formatted = _format_default_like(authoritative, value)
        escaped_value = html.escape(formatted, quote=True)
        return keep(
            f'<input class="{css_class}" data-param="{escaped_key}" '
            f'data-default="{escaped_value}" value="{escaped_value}" '
            f'aria-label="Editable parameter {escaped_key}" spellcheck="false">'
        )

    text = re.sub(r"\[\[(input|assumed):([A-Za-z0-9_]+)=([^\]]+)\]\]", parameter_input, text)

    def derived_output(match: re.Match[str]) -> str:
        key, value = match.group(1), match.group(2)
        escaped_key = html.escape(key, quote=True)
        authoritative = browser_derived_defaults().get(key, value.strip())
        formatted = _format_default_like(authoritative, value)
        escaped_value = html.escape(formatted)
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
p .derived-output,li .derived-output,blockquote .derived-output {{ width:auto; min-width:0; padding:.08rem .32rem; vertical-align:baseline; }}
.live-equation {{ margin:.7rem 0; padding:.7rem .85rem; overflow-x:auto; border:1px dashed var(--accent); border-radius:7px; background:var(--soft); color:var(--text); font:600 .92rem/1.5 "Cascadia Code",Consolas,monospace; font-variant-numeric:tabular-nums; }}
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
<div class="topbar"><span class="name">{html.escape(title)}</span><button onclick="setDetails(true)">Expand derivations</button><button onclick="setDetails(false)">Collapse</button><button onclick="saveParameterInputs()">Save in browser</button><button onclick="saveEditedHtml()">Save HTML copy</button><button onclick="resetParameterInputs()">Reset inputs</button><button class="hide-small" onclick="toggleTheme()">Theme</button><button class="hide-small" onclick="window.print()">Print</button></div>
<div class="layout"><nav><div class="caption">On this page</div>{''.join(toc_html)}</nav><article><div class="edit-note"><span class="assumed-swatch"></span>Amber inputs are unidentified assumptions. “Save in browser” stores overrides only in this browser and the page URL; it does not rewrite the workspace file. “Save HTML copy” embeds them in a chosen file. Dependent scalar values, live equations, and the live Bode panel recalculate immediately. Publication SVGs and nonlinear simulations require a Python rebuild.<span id="parameter-save-status" class="save-status">Loading values…</span></div>{body}<div class="footer">Rendered from {html.escape(markdown_path.name)} · {generated}</div></article></div>
<script>
function setDetails(open) {{ document.querySelectorAll('details').forEach(d => d.open=open); }}
function toggleTheme() {{ const root=document.documentElement; root.dataset.theme=root.dataset.theme==='dark'?'light':'dark'; }}
const parameterStorageKey = 'model-parameters:dependency-v2:' + document.title + ':' + location.pathname;
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
  setParameterStatus('Saved in this browser · live dependencies updated · rebuild required for publication simulations', 'warn');
}}
function scheduleParameterUpdate() {{
  setParameterStatus('Updating dependent values and live plots…', 'warn');
  if (parameterSaveTimer) clearTimeout(parameterSaveTimer);
  parameterSaveTimer = setTimeout(() => {{
    persistParameterInputs(false);
    refreshInteractivePlots();
    setParameterStatus('Auto-saved in this browser · live dependencies updated · static simulations are now stale', 'warn');
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
  const screwLength = parameterNumber('screw_length', 0.192);
  const usableScrewTravel = parameterNumber('usable_screw_travel', 0.170);
  const stageTravel = parameterNumber('stage_travel', 0.150);
  const screwDiameter = parameterNumber('screw_diameter', 8.0e-3);
  const screwDensity = parameterNumber('screw_density', 7850.0);
  const tmax = parameterNumber('holding_torque', 0.060);
  const tdet = parameterNumber('detent_torque', 0.005);
  const detentPhase = parameterNumber('detent_phase', 0.0);
  const couplingSeries = parameterNumber('k_c_series', 68.7549);
  const mStage = parameterNumber('stage_mass', 0.355);
  const mNut = parameterNumber('nut_mass', 0.050);
  const ms = mStage + mNut;
  const axialModeTarget = parameterNumber('axial_mode_target_hz', 695.82);
  const kbrg = parameterNumber('k_brg', 25.0e6);
  const ksha = parameterNumber('k_sha', 67.0e6);
  const kmnt = parameterNumber('k_mnt', 100.0e6);
  const cax = parameterNumber('axial_damping', 55.0);
  const zeta = parameterNumber('electromagnetic_zeta', 0.10);
  const microstepDivisor = parameterNumber('microstep_divisor', 64);
  if (!(lead>0 && teeth>0 && jm>0 && jc>=0 && screwLength>0 && usableScrewTravel>0 &&
        stageTravel>0 && stageTravel<=usableScrewTravel && usableScrewTravel<=screwLength && screwDiameter>0 &&
        screwDensity>0 && tmax>0 && tdet>=0 && mStage>0 && mNut>=0 &&
        axialModeTarget>0 && kbrg>0 && ksha>0 && kmnt>0 &&
        cax>=0 && zeta>=0 && microstepDivisor>=1))
    throw new Error('Geometry, masses, torque, and stiffness must be positive; damping and detent torque cannot be negative.');
  const r = lead/(2*Math.PI);
  const screwRadius = screwDiameter/2;
  const screwMass = screwDensity*Math.PI*screwRadius*screwRadius*screwLength;
  const screwInertia = 0.5*screwMass*screwRadius*screwRadius;
  const jTotal = jm+jc+screwInertia;
  const md = jTotal/(r*r);
  const km = teeth*tmax/(r*r);
  const modalLambda = Math.pow(2*Math.PI*axialModeTarget,2);
  const modalDenominator = km-modalLambda*(md+ms);
  if (Math.abs(modalDenominator) <= Number.EPSILON*Math.max(Math.abs(km),1))
    throw new Error('The axial-mode calibration is singular for the current inputs.');
  const kax = modalLambda*ms*(km-modalLambda*md)/modalDenominator;
  if (!(kax>0 && Number.isFinite(kax)))
    throw new Error('The selected mass and modal target do not yield a positive axial stiffness.');
  const remainingBallCompliance = 1/kax-1/kbrg-1/ksha-1/kmnt;
  if (!(remainingBallCompliance>0))
    throw new Error('The current axial inputs leave no positive compliance for k_ball.');
  const kBall = 1/remainingBallCompliance;
  const kdetAmplitude = 4*teeth*tdet/(r*r);
  const kdet = kdetAmplitude*Math.cos(detentPhase);
  if (!(km>kdetAmplitude)) throw new Error('The detent tangent amplitude must remain below commutation stiffness.');
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
  function modePair(driveTangent) {{
    const qa=md*ms, qb=md*kax+ms*(driveTangent+kax), qc=driveTangent*kax;
    const discriminant=Math.max(qb*qb-4*qa*qc,0);
    const roots=[(qb-Math.sqrt(discriminant))/(2*qa),(qb+Math.sqrt(discriminant))/(2*qa)];
    return roots.map(value => Math.sqrt(Math.max(value,0))/(2*Math.PI));
  }}
  const modes=modePair(km);
  const localSoftModes=modePair(km-kdetAmplitude), localHardModes=modePair(km+kdetAmplitude);
  const localLowBand=[Math.min(localSoftModes[0],localHardModes[0]),Math.max(localSoftModes[0],localHardModes[0])];
  return {{
    frequencies, drive, stage, rotorStage, modes, localLowBand, md, ms, mStage, mNut,
    axialModeTarget, km, kdet, kdetAmplitude, kax, kBall, kbrg, ksha, kmnt,
    cax, zeta, lead, teeth, r, jm, jc, screwLength, screwDiameter,
    screwDensity, screwMass, screwInertia, jTotal, tmax, tdet, detentPhase,
    usableScrewTravel, stageTravel,
    couplingSeries, couplingHalf:2*couplingSeries, kappa:teeth/r,
    fullStep:lead/(4*teeth), quarterStep:lead/(16*teeth),
    commandStep:lead/(4*teeth*microstepDivisor), interpolatedStep:lead/(4*teeth*256),
    microstepDivisor, fmax:tmax/r, cm, detentPeriod:lead/(4*teeth)
  }};
}}
function formatDerivedValue(key, value) {{
  const scientific = new Set([
    'transmission_ratio','magnetic_stiffness','detent_stiffness','screw_inertia',
    'screw_segment_inertia','full_step_pitch','quarter_step_bound',
    'command_step','interpolated_step','total_rotational_inertia',
    'reduced_axial_stiffness','k_ball'
  ]);
  if (scientific.has(key)) return value.toExponential(5);
  if (key==='reduced_drive_mass') return value.toFixed(3);
  if (key==='reduced_stage_mass') return value.toFixed(3);
  if (key==='screw_mass' || key==='screw_segment_mass') return value.toFixed(6);
  if (key==='k_c_half') return value.toFixed(3);
  if (key.endsWith('_hz')) return value.toFixed(2);
  return Number(value).toPrecision(6);
}}
function refreshDerivedOutputs(data) {{
  const values = {{
    transmission_ratio:data.r,
    total_rotational_inertia:data.jTotal,
    reduced_drive_mass:data.md,
    reduced_stage_mass:data.ms,
    magnetic_stiffness:data.km,
    detent_stiffness:data.kdet,
    reduced_axial_stiffness:data.kax,
    k_ball:data.kBall,
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
    drive_stiffness:data.km,
    force_limit:data.fmax,
    spatial_wavenumber:data.kappa
  }};
  document.querySelectorAll('[data-derived]').forEach(output => {{
    const key=output.dataset.derived;
    if (Object.prototype.hasOwnProperty.call(values,key)) output.textContent=formatDerivedValue(key,values[key]);
  }});
}}
function refreshLiveEquations(data) {{
  const mn = value => (value/1e6).toFixed(3);
  const reflectedKg = inertia => inertia/(data.r*data.r);
  const compliancePct = stiffness => 100*data.kax/stiffness;
  const equations = {{
    'reduced-mass':
      'm_s = m_stage + m_n = ' + data.mStage.toFixed(3) + ' + ' +
      data.mNut.toFixed(3) + ' = ' + data.ms.toFixed(3) + ' kg',
    'inertia-aggregation':
      'J_sum = J_m + J_c + J_s = ' + data.jTotal.toExponential(5) +
      ' kg m^2; reflected contributions: motor ' + reflectedKg(data.jm).toFixed(3) +
      ' kg + coupling ' + reflectedKg(data.jc).toFixed(3) +
      ' kg + screw ' + reflectedKg(data.screwInertia).toFixed(3) +
      ' kg = m_d ' + data.md.toFixed(3) + ' kg',
    'compliance-breakdown':
      'Current compliance shares: bearing ' + compliancePct(data.kbrg).toFixed(2) +
      '%, loaded screw ' + compliancePct(data.ksha).toFixed(2) +
      '%, ball contact ' + compliancePct(data.kBall).toFixed(2) +
      '%, nut mount ' + compliancePct(data.kmnt).toFixed(2) + '% (sum 100.00%).',
    'modal-stiffness':
      'f₂,target = ' + data.axialModeTarget.toFixed(2) + ' Hz  →  k_ax = ' +
      mn(data.kax) + ' MN/m',
    'axial-compliance':
      '1/' + mn(data.kbrg) + ' + 1/' + mn(data.ksha) + ' + 1/' +
      mn(data.kBall) + ' + 1/' + mn(data.kmnt) + ' = 1/' +
      mn(data.kax) + '  (MN/m)⁻¹'
  }};
  document.querySelectorAll('[data-live-equation]').forEach(element => {{
    const equation = equations[element.dataset.liveEquation];
    if (equation) element.textContent = equation;
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
  try {{ data=liveTransferData(); refreshDerivedOutputs(data); refreshLiveEquations(data); }}
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
  document.getElementById('live-model-summary').textContent='Live global model: mstage='+data.mStage.toFixed(3)+' kg, mn='+data.mNut.toFixed(3)+' kg, ms='+data.ms.toFixed(3)+' kg; target f2='+data.axialModeTarget.toFixed(2)+' Hz gives kax='+(data.kax/1e6).toFixed(3)+' MN/m and kball='+(data.kBall/1e6).toFixed(3)+' MN/m; global modes '+data.modes.map(value=>value.toFixed(2)+' Hz').join(', ')+'; local detent low-pole band '+data.localLowBand.map(value=>value.toFixed(2)+' Hz').join(' to ')+'.';
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
    validate_parameter_registry()
    validate_case_topology()
    constants = physical_constants()
    frequencies, bode, linear_metrics = frequency_responses()
    times, command, time_data, time_metrics = time_responses(constants)
    convergence = gms_step_halving_convergence(constants, times, time_data)
    memory_experiments = {
        "guideway": presliding_responses(constants, ("A", "A2"), "g"),
        "nut": presliding_responses(constants, ("B", "B2"), "n"),
    }
    verification = full_reduced_verification(frequencies, constants)
    case_response_paths = plot_case_responses(
        frequencies, bode, times, command, time_data, constants, time_metrics)
    comparison_path = plot_case_response_overlay(frequencies, bode)
    guide_memory_path = plot_presliding_memory(
        memory_experiments["guideway"], "presliding_memory_comparison.svg")
    nut_memory_path = plot_presliding_memory(
        memory_experiments["nut"], "nut_memory_comparison.svg")
    diagram_paths = plot_kinematic_diagram()
    flowchart_a_path = plot_flowchart_provenance_structure()
    flowchart_b_path = plot_flowchart_friction_results()
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
        update_generated_presliding_summary(generated_presliding_summary(memory_experiments))
        update_generated_convergence_summary(generated_convergence_summary(convergence))
    description_html = render_document(DESCRIPTION_MD)
    derivation_html = render_document(DERIVATION_MD)
    print(f"Built {comparison_path.relative_to(ROOT)}")
    print(f"Built {len(case_response_paths)} per-case response figures")
    print(f"Built {guide_memory_path.relative_to(ROOT)}")
    print(f"Built {nut_memory_path.relative_to(ROOT)}")
    for diagram_path in diagram_paths:
        print(f"Built {diagram_path.relative_to(ROOT)}")
    print(f"Built {flowchart_a_path.relative_to(ROOT)}")
    print(f"Built {flowchart_b_path.relative_to(ROOT)}")
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
          f"peak={verification['peak_residual_nm']:.3f} nm; "
          f"normalized={verification['rms_residual_pct_command']:.3f}%/"
          f"{verification['peak_residual_pct_command']:.3f}%")
    for experiment in memory_experiments.values():
        for key in experiment["keys"]:
            metric = experiment["metrics"][key]
            print(f"Memory {key}: RMS={metric['whole_rms_nm']:.3f} nm; "
                  f"return mismatch={metric['return_error_mismatch_nm']:.3f} nm; "
                  f"force closure={metric['return_force_mismatch_N']:.6f} N; "
                  f"final mean={metric['final_mean_nm']:.3f} nm")
    for key in ("A2", "B2", "C2"):
        result = convergence[key]
        rms = result["rms_nm"]
        dt_text = "/".join(f"{dt * 1e6:g}" for dt in GMS_CONVERGENCE_DTS)
        print(f"Step halving {key}: RMS({dt_text} us)="
              f"{rms[0]:.6f}/{rms[1]:.6f}/{rms[2]:.6f} nm; "
              f"fine relative change={result['fine_relative_pct']:.6f}%")
    print("Full-model modes below 3 kHz: " + ", ".join(
        f"{mode:.2f}" for mode in verification["full_modes"] if mode < 3000.0))


if __name__ == "__main__":
    main()

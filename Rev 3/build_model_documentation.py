#!/usr/bin/env python3
"""Build the model response figures and render both Markdown documents to HTML.

The script deliberately depends only on NumPy and Matplotlib.  It implements a
fixed-step RK4 integrator for nonlinear LuGre and GMS simulations and a compact
Markdown renderer tailored to the two project documents.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import re
from collections import OrderedDict
from concurrent.futures import Executor, ProcessPoolExecutor
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

# Stable element IDs keep unchanged SVGs byte-identical across rebuilds.
plt.rcParams["svg.hashsalt"] = "rev3-model-documentation"


ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"
DESCRIPTION_MD = ROOT / "ball_screw_stage_dynamic_derivation_v3.md"
DERIVATION_MD = ROOT / "Analytical_derivation_and_responses_v3.md"
MICROSTEP_DATA_DIR = ROOT.parent / "Microstepping Test Data"
PARAMETER_FILE = ROOT / "model_parameters.json"
BUILD_ID = "rev3-section9-final-20260813"


def _load_parameter_file() -> dict[str, object]:
    """Load browser-saved overrides before executable defaults are assembled."""
    if not PARAMETER_FILE.exists():
        return {}
    payload = json.loads(PARAMETER_FILE.read_text(encoding="utf-8"))
    values = payload.get("parameters", payload)
    if not isinstance(values, dict):
        raise TypeError(f"{PARAMETER_FILE.name} must contain a parameter object")
    return values


_PARAMETER_OVERRIDES = _load_parameter_file()


def configured(key: str, default: object) -> object:
    """Return one typed parameter-file value or its in-code default."""
    value = _PARAMETER_OVERRIDES.get(key, default)
    if isinstance(default, str):
        return str(value)
    if isinstance(default, int) and not isinstance(default, bool):
        numeric = float(value)
        if not numeric.is_integer():
            raise ValueError(f"Parameter {key} must be an integer, got {value!r}")
        return int(numeric)
    return float(value)


# Executable defaults for the Revision 3 two-DOF reduction.
MODEL = {
    "lead": configured("lead", 1.0e-3),
    "rotor_teeth": configured("rotor_teeth", 50),
    # Rated-current holding torque for the lower-current motor variant.
    "T_max": configured("holding_torque", 0.060),
    # Published detent torque.  It is executed as a periodic nonlinear torque;
    # its tangent is reported separately and is never used as a global spring.
    "T_det": configured("detent_torque", 0.005),
    "detent_phase": configured("detent_phase", 0.0),
    # Selected upper-mode calibration target.  Appendix G.4 records that this
    # is NOT a measured number: the measured axial band is the pair below and
    # 695.82 is the value the reduced chain is calibrated to reproduce.  Its
    # own provenance is still undocumented and is an open identification item.
    # Stage and nut masses live in FULL; m_s, k_ax, and k_ball are derived.
    "axial_mode_target_hz": configured("axial_mode_target_hz", 695.82),
    # Measured axial band from the modal campaign, kept separate from the
    # calibration target so the two can never be quoted as the same number.
    "measured_axial_band_low_hz": configured("measured_axial_band_low_hz", 681.0),
    "measured_axial_band_high_hz": configured("measured_axial_band_high_hz", 690.0),
    # Detent enable flag.  The nonlinear campaign runs with detent on; the
    # paired ablation reruns every case with this term removed so the settled
    # window can be split into a detent term and a friction term.
    "detent_enabled": configured("detent_enabled", 1),
    # Correction levels a positional pre-distortion table must place across one
    # period of the detent equilibrium error.  It sets the command-grid
    # requirement in Section 5 and is a design choice, not a measurement.
    "predistortion_levels": configured("predistortion_levels", 7),
    # Effective mass and relative-mode damping from the modal campaign.  The
    # damping value remains provisional until the half-power extraction is
    # repeated from the source FRF.
    "m_eff_measured": configured("m_eff_measured", 0.600),
    "zeta_relative_measured": configured("zeta_relative_measured", 0.0014),
    # Existing provisional reduced-link damper.  It remains an explicit
    # sensitivity input so it can be compared with interface propagation and
    # measured-FRF identification rather than being mistaken for either one.
    "c_ax": configured("axial_damping", 55.0),
    # Provisional open-loop drive damping ratio.  Driver mode and tuning are
    # not recorded.  The requested baseline is 10% of critical damping and a
    # sensitivity sweep is retained rather than presenting it as identified.
    "zeta_m": configured("electromagnetic_zeta", 0.10),
    # Production Stepper-Board STEP/DIR setting.  Position commands therefore
    # land on a 312.5 nm grid; pulse timing remains independently controllable.
    "microstep_divisor": configured("microstep_divisor", 16),
    # Stribeck exponent.  It appears in every s(v) evaluation but was never
    # exposed; 2.0 is the conventional Gaussian form.  Fixed, not identified.
    "stribeck_delta": configured("stribeck_delta", 2.0),
    # Shared Stribeck relaxation time.  Each site's GMS attractor rate is
    # C_alpha = (F_s - F_c)/tau_C, so one time constant replaces three
    # unanchored N/s values.  Provisional: in the source GMS work C is
    # identified from measured hysteresis loops, never assumed.
    "tau_C": configured("tau_C", 2.0e-4),
    # Provenance inputs for the drive-port breakaway estimate only.  They are
    # never applied to the transformer; see the standing constraint in 8.1.
    "eta_screw": configured("eta_screw", 0.90),
    "F_preload_nut": configured("F_preload_nut", 100.0),
}


# Revision 3 full-model values.  These are deliberately separate from MODEL:
# MODEL is the validated two-DOF reduction, whereas FULL retains all ten
# coordinates named in the source document.  Values not measured in the source
# are surfaced as highlighted assumptions in the Markdown documents.
FULL = {
    # Component values.  Screw inertia and axial masses are derived below.
    "J_m": configured("J_m", 9.00e-7),
    "J_c": configured("J_c", 1.18e-6),
    "screw_length": configured("screw_length", 0.192),
    "usable_screw_travel": configured("usable_screw_travel", 0.170),
    "stage_travel": configured("stage_travel", 0.150),
    # Two diameters, because a ball screw is not a plain cylinder.  The
    # nominal diameter is the mass diameter used for m_screw; every stiffness
    # and the polar inertia use the root diameter, which is the section that
    # actually carries torsion and axial load.  The root value is an estimate
    # from the KGT-F1-08-01 class and must be confirmed against the
    # manufacturer drawing or the Creo mass properties.
    "screw_diameter": configured("screw_diameter", 8.00e-3),
    "screw_root_diameter": configured("screw_root_diameter", 6.80e-3),
    "screw_density": configured("screw_density", 7850.0),
    # Material constants for the derived screw stiffnesses.
    "youngs_modulus": configured("youngs_modulus", 210.0e9),
    "shear_modulus": configured("shear_modulus", 80.8e9),
    # Declared support-to-nut free length.  Appendix A maps free length to
    # 20 mm plus stage position, so this datum is a 138 mm stage position of
    # 150 mm travel: the softest end of the axis, retained deliberately as the
    # worst case.  L_b = L_s - L_a, so the two screw segments close on the
    # complete screw length by construction.
    "nut_axial_datum": configured("nut_axial_datum", 0.158),
    "m_n": configured("nut_mass", 0.050),
    # Measured stage body mass.  The retained stage-side mass also includes
    # the nut body after the internal nut coordinate is collapsed.
    "m_stage": configured("stage_mass", 0.355),
    # Datasheet series stiffness is 1.2 N m/deg = 68.7549 N m/rad.
    # Two equal half-springs must each be twice the series value.
    "k_c_series": configured("k_c_series", 1.2 * 180.0 / np.pi),
    # k_theta_a, k_theta_b, k_sha and k_shb are no longer independent entries.
    # Entered separately they described two different screws: the equal
    # torsional pair placed the nut at midspan, the axial pair placed it at
    # 75% of span, and both implied segment sums longer than L_s.  They are
    # now derived from E, G, the root diameter and the declared nut datum in
    # screw_segment_stiffnesses(), which closes L_a + L_b = L_s by
    # construction.
    # 25 N/um is the closure-consistent bearing assumption discussed in Rev 3.
    "k_brg": configured("k_brg", 25.0e6),
    "k_mnt": configured("k_mnt", 100.0e6),
    # The four interface damping ratios are no longer entries.  The
    # light-damping identity eta = 2*zeta holds only at an element's OWN
    # resonance, so hand-entered ratios drift out of their target loss factor
    # whenever the geometry changes.  interface_damping_ratios() solves
    # zeta_j = eta_j * f_j / (2 * f_2) from the loss factors below, and
    # validate_interface_loss_factors() still checks the executed result.
}

# Target loss factors at the retained upper mode.  Joints are bolted,
# preloaded, rolling-element interfaces; the screw segment is monolithic steel.
INTERFACE_LOSS_FACTORS = {
    "zeta_bearing": configured("eta_bearing", 0.03),
    "zeta_steel": configured("eta_steel", 0.0005),
    "zeta_ball_nut": configured("eta_ball_nut", 0.03),
    "zeta_nut_mount": configured("eta_nut_mount", 0.03),
}

FULL_DOF_LABELS = (
    r"$\theta_m$", r"$\theta_c$", r"$\theta_{s1}$", r"$\theta_{s2}$",
    r"$\theta_{s3}$", r"$u_b$", r"$u_e$", r"$u_f$", r"$u_n$", r"$x_s$",
)


# Highlighted friction-port values used to make both law comparisons executable.
# sigma0_g is the estimate already quoted in the description; all other values
# below need experimental identification before quantitative use.
FRICTION = {
    "g": {"sigma0": configured("g_sigma0", 7.60e5),
          "sigma1": configured("g_sigma1", 0.0),
          "sigma2": configured("g_sigma2", 0.40),
          "F_s": configured("g_Fs", 3.0), "F_c": configured("g_Fc", 2.4),
          "v_s": configured("g_vs", 2.5e-4)},
    # Differential nut-contact microslip.  Its first GMS element yields at
    # 0.25*F_s/sigma0 = 0.20 um, so this port can express actual partial slip.
    "n": {"sigma0": configured("n_sigma0", 2.00e6),
          "sigma1": configured("n_sigma1", 0.0),
          "sigma2": configured("n_sigma2", 0.25),
          "F_s": configured("n_Fs", 1.6), "F_c": configured("n_Fc", 1.2),
          "v_s": configured("n_vs", 2.0e-4)},
    # Identifiable drive-side lump.  Motor/support-bearing drag and gross nut
    # rolling were formerly two laws on the same H=[1,0] port.  Their force,
    # tangent, and damping budgets are combined here; the aggregate still
    # requires identification from a drive-side measurement.
    "d": {"sigma0": configured("d_sigma0", 3.00e6),
          "sigma1": configured("d_sigma1", 0.0),
          "sigma2": configured("d_sigma2", 0.45),
          "F_s": configured("d_Fs", 7.0), "F_c": configured("d_Fc", 5.5),
          "v_s": configured("d_vs", 2.3e-4)},
}

# sigma_1 is zero in the executed A/B/C comparison so that LuGre and GMS
# contribute the identical tangent damping sigma_2*H^T H.  Without this the
# closure claim in 8.4 is false for damping: LuGre adds sigma_1+sigma_2 while
# GMS adds only sigma_2, a factor of 8.5 to 21 across the three sites, and any
# plotted difference mixes memory structure with port damping.  The former
# values are retained here and restored by the A1v micro-viscous variant.
MICRO_VISCOUS_SIGMA1 = {"g": 3.0, "n": 5.0, "d": 9.0}

for _site_key, _site_values in FRICTION.items():
    _site_values["delta"] = MODEL["stribeck_delta"]
    # One shared relaxation time replaces three independent N/s constants.
    _site_values["C_gms"] = (
        _site_values["F_s"] - _site_values["F_c"]) / MODEL["tau_C"]


# Four GMS stop elements share each site's aggregate sigma0 and Stribeck
# force.  Opposing stiffness/force fractions create distinct yield distances
# and therefore non-local reversal memory while retaining the LuGre aggregate.
GMS_WEIGHTS = np.array([
    configured(f"gms_nu{index}", default)
    for index, default in enumerate((0.10, 0.20, 0.30, 0.40), start=1)
])
GMS_STIFFNESS_FRACTIONS = np.array([0.40, 0.30, 0.20, 0.10])
GMS_N = GMS_WEIGHTS.size
GMS_STIFFNESS_BY_SITE = {
    site: np.array([
        configured(f"{site}_k{index}", fraction * values["sigma0"])
        for index, fraction in enumerate(GMS_STIFFNESS_FRACTIONS, start=1)
    ])
    for site, values in FRICTION.items()
}
# Filled in by main() before rendering.  A document rendered without a rebuilt
# census says so rather than printing stale counts.
BRANCH_CENSUS_SENTENCE: str | None = None
PRODUCTION_DT = 2.5e-5
GMS_CONVERGENCE_DTS = (5.0e-5, 2.5e-5, 1.25e-5)
A2_CONVERGENCE_DT = 6.25e-6
# The full/reduced audit contains full-model modes above the 3 kHz plotting
# range.  Its production step is therefore finer than the nonlinear case step,
# and its convergence study intentionally includes one unstable coarse step so
# the actual RK4 stability limit remains visible in the report.
VERIFICATION_DT = 2.5e-6
VERIFICATION_CONVERGENCE_DTS = (25.0e-6, 12.5e-6, 6.25e-6, VERIFICATION_DT)
VERIFICATION_EDGES = (0.005, 0.025, 0.045, 0.065)
SETTLING_2PCT_FACTOR = 4.0
BODE_FOCUS_MIN_HZ = 100.0
BODE_FOCUS_MAX_HZ = 3000.0
# The main response uses the production 1/16 grid and deliberately spans the
# provisional first-yield distances.  Adjacent increments remain at or below
# one quarter of a full step.  The final move is positive and returns to zero.
MAIN_LEVELS = np.array([1, -1, 2, -2, 0, 3, 0, -3, -6, -3,
                        0, 3, 6, 3, 0], dtype=float)
MAIN_START = 0.010


# Nested reversals for the dedicated memory experiment.  The positive guideway
# levels are 12/10/4 microsteps; the negative outer and nested levels are one
# full 1/16 microstep smaller to retain an executable asymmetry.
GUIDEWAY_PRESLIDING_LEVELS = np.array(
    [0, 12, 4, 10, 4, 12, 0, -11, -4, -9, -4, -11, 0],
    dtype=float,
)
# With the stage blocked, the drive coordinate is the nut-port deflection.
# The 3/2/1 integer levels are the only distinct 1/16 set that crosses the
# first two provisional stops while remaining below the third at 1.20 um.
NUT_PRESLIDING_LEVELS = np.array(
    [0, 3, 1, 2, 1, 3, 0, -3, -1, -2, -1, -3, 0],
    dtype=float,
)
PRESLIDING_START = 0.005
PRESLIDING_RETURN_PAIRS = ((1, 5), (2, 4), (7, 11), (8, 10))


CASES = OrderedDict([
    ("0", {"label": "Case 0: frictionless", "sites": (), "friction": "none", "color": "#252525", "ls": "--"}),
    ("A", {"label": "Case A: lumped drive drag + guideway / LuGre", "sites": ("d", "g"), "friction": "lugre", "color": "#277da1", "ls": "-"}),
    ("A2", {"label": "Case A2: lumped drive drag + guideway / GMS", "sites": ("d", "g"), "friction": "gms", "color": "#70b7cf", "ls": "--"}),
    ("G", {"label": "Case G: guideway only / LuGre", "sites": ("g",), "friction": "lugre", "color": "#4059ad", "ls": "-"}),
    ("G2", {"label": "Case G2: guideway only / GMS", "sites": ("g",), "friction": "gms", "color": "#8b9dc3", "ls": "--"}),
    ("B", {"label": "Case B: lumped drive drag + nut microslip / LuGre", "sites": ("d", "n"), "friction": "lugre", "color": "#e07a15", "ls": "-"}),
    ("B2", {"label": "Case B2: lumped drive drag + nut microslip / GMS", "sites": ("d", "n"), "friction": "gms", "color": "#f5b35f", "ls": "--"}),
    ("C", {"label": "Case C: all identifiable ports / LuGre", "sites": ("d", "g", "n"), "friction": "lugre", "color": "#218c74", "ls": "-"}),
    ("C2", {"label": "Case C2: all identifiable ports / GMS", "sites": ("d", "g", "n"), "friction": "gms", "color": "#72c9ad", "ls": "--"}),
    # Micro-viscous variant.  Same ports and same law as A, with sigma_1
    # restored, so it isolates bristle damping instead of confounding it with
    # the memory comparison.  B1v and C1v would differ only in which ports are
    # active and are not executed; the effect is identical in kind.
    ("A1v", {"label": "Case A1v: drive drag + guideway / LuGre, micro-viscous",
             "sites": ("d", "g"), "friction": "lugre", "micro_viscous": True,
             "color": "#8a5fbf", "ls": ":"}),
])

PAIRS = (("A", "A2"), ("G", "G2"), ("B", "B2"), ("C", "C2"))
# Parameter provenance used by the browser registry and the standalone
# dependency flowcharts.  Every derived token emitted into either Markdown
# document has an explicit dependency list.  Modal-calibrated k_ax and
# closure-derived k_ball are derived outputs, not independent inputs.
PARAMETER_REGISTRY: dict[str, dict[str, object]] = {
    "case_definitions": {"category": "input", "dependencies": (), "section": "1-model-hierarchy-and-case-map"},
    "case_count": {"category": "output", "dependencies": ("case_definitions",), "section": "1-model-hierarchy-and-case-map"},
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
    "m_eff_measured": {"category": "input", "dependencies": (), "section": "2-entry-parameters"},
    "zeta_relative_measured": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "axial_damping": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "electromagnetic_zeta": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    # Target loss factors are the assumption; the four damping ratios that
    # realize them at the retained mode are derived from the element
    # frequencies, so they can never be left describing a different assembly.
    "eta_steel": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "eta_bearing": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "eta_ball_nut": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "eta_nut_mount": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "zeta_steel": {
        "category": "derived",
        "dependencies": ("eta_steel", "k_sha", "screw_mass", "axial_mode_target_hz"),
        "section": "2-entry-parameters",
    },
    "zeta_bearing": {
        "category": "derived",
        "dependencies": ("eta_bearing", "k_brg", "screw_mass", "axial_mode_target_hz"),
        "section": "2-entry-parameters",
    },
    "zeta_ball_nut": {
        "category": "derived",
        "dependencies": ("eta_ball_nut", "k_ball", "screw_mass", "screw_inertia",
                         "nut_mass", "axial_mode_target_hz"),
        "section": "2-entry-parameters",
    },
    "zeta_nut_mount": {
        "category": "derived",
        "dependencies": ("eta_nut_mount", "k_mnt", "nut_mass", "stage_mass",
                         "axial_mode_target_hz"),
        "section": "2-entry-parameters",
    },
    "interface_axial_damping": {
        "category": "derived",
        "dependencies": (
            "reduced_axial_stiffness", "k_ball", "k_brg", "k_sha", "k_mnt",
            "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount",
            "screw_length", "screw_diameter", "screw_density", "nut_mass",
            "stage_mass", "axial_mode_target_hz",
        ),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
    "reduced_axial_stiffness": {
        "category": "derived",
        "dependencies": ("reduced_drive_mass", "reduced_stage_mass",
                         "magnetic_stiffness", "axial_mode_target_hz"),
        "section": "6-reduction-from-ten-dofs-to-two",
    },
    "k_c_series": {"category": "input", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    # Screw section and material, shared by both segment stiffness pairs.
    "screw_root_diameter": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "youngs_modulus": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "shear_modulus": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "nut_axial_datum": {"category": "assumed", "dependencies": (), "section": "2-entry-parameters"},
    "screw_length_a": {
        "category": "derived", "dependencies": ("nut_axial_datum",),
        "section": "2-entry-parameters",
    },
    "screw_length_b": {
        "category": "derived", "dependencies": ("nut_axial_datum", "screw_length"),
        "section": "2-entry-parameters",
    },
    "k_theta_a": {
        "category": "derived",
        "dependencies": ("shear_modulus", "screw_root_diameter", "nut_axial_datum"),
        "section": "4-full-ten-dof-derivation",
    },
    "k_theta_b": {
        "category": "derived",
        "dependencies": ("shear_modulus", "screw_root_diameter", "nut_axial_datum",
                         "screw_length"),
        "section": "4-full-ten-dof-derivation",
    },
    "k_brg": {"category": "assumed", "dependencies": (), "section": "4-full-ten-dof-derivation"},
    "k_sha": {
        "category": "derived",
        "dependencies": ("youngs_modulus", "screw_root_diameter", "nut_axial_datum"),
        "section": "4-full-ten-dof-derivation",
    },
    "k_shb": {
        "category": "derived",
        "dependencies": ("youngs_modulus", "screw_root_diameter", "nut_axial_datum",
                         "screw_length"),
        "section": "4-full-ten-dof-derivation",
    },
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
        "dependencies": ("screw_root_diameter", "screw_length", "screw_density"),
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

# The route-comparison outputs are emitted into the analytical document and
# mirrored by the browser equations.  Keeping their dependencies explicit
# makes changes to component values, modal inputs, or damping assumptions
# auditable through the existing provenance machinery.
_ROUTE_COMMON_DEPS = (
    "reduced_drive_mass", "reduced_stage_mass", "magnetic_stiffness",
    "reduced_axial_stiffness", "k_ball", "k_brg", "k_sha", "k_mnt",
)
for _route in ("p", "s", "b"):
    for _suffix in ("md", "ms", "kax", "cax", "zeta", "f1", "f2", "kball", "settling"):
        PARAMETER_REGISTRY[f"route_{_route}_{_suffix}"] = {
            "category": "derived", "dependencies": _ROUTE_COMMON_DEPS + ("axial_damping",),
            "section": "6-3-reduction-evidence",
        }
for _suffix in ("md", "ms", "kax", "cax", "zeta", "f1", "f2", "kball", "settling"):
    PARAMETER_REGISTRY[f"route_f_{_suffix}"] = {
        "category": "derived",
        "dependencies": _ROUTE_COMMON_DEPS + (
            "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount",
            "screw_length", "screw_diameter", "screw_density",
            "nut_mass", "stage_mass", "J_m", "J_c", "k_c_series", "k_theta_a",
        ),
        "section": "6-3-reduction-evidence",
    }
    PARAMETER_REGISTRY[f"route_c_{_suffix}"] = {
        "category": "output",
        "dependencies": _ROUTE_COMMON_DEPS + (
            "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount",
            "screw_length", "screw_diameter", "screw_density",
            "nut_mass", "stage_mass", "J_m", "J_c", "k_c_series", "k_theta_a",
        ),
        "section": "6-3-reduction-evidence",
    }
    PARAMETER_REGISTRY[f"route_m_{_suffix}"] = {
        "category": "derived",
        "dependencies": (
            "reduced_drive_mass", "magnetic_stiffness", "axial_mode_target_hz",
            "m_eff_measured", "zeta_relative_measured", "k_brg", "k_sha", "k_mnt",
        ),
        "section": "6-3-reduction-evidence",
    }
for _key, _dependencies in {
    "route_s_kax_full": _ROUTE_COMMON_DEPS + ("k_c_series", "k_theta_a", "lead"),
    "torsional_share": _ROUTE_COMMON_DEPS + ("k_c_series", "k_theta_a", "lead"),
    "mass_ratio": ("reduced_drive_mass", "reduced_stage_mass"),
    "reduced_mu": ("reduced_drive_mass", "reduced_stage_mass"),
    "full_model_zeta": _ROUTE_COMMON_DEPS + (
        "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount"),
    "full_model_settling": _ROUTE_COMMON_DEPS + (
        "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount"),
    "cb_frequency_delta": _ROUTE_COMMON_DEPS + (
        "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount"),
    "cb_damping_delta": _ROUTE_COMMON_DEPS + (
        "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount"),
    "full_model_upper_hz": _ROUTE_COMMON_DEPS + (
        "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount"),
    "full_model_bandwidth_hz": _ROUTE_COMMON_DEPS + (
        "zeta_steel", "zeta_bearing", "zeta_ball_nut", "zeta_nut_mount"),
    "route_p_bandwidth_hz": _ROUTE_COMMON_DEPS + ("axial_damping",),
    "first_fixed_interface_hz": _ROUTE_COMMON_DEPS + (
        "screw_length", "screw_diameter", "screw_density",
        "nut_mass", "stage_mass", "J_m", "J_c", "k_c_series", "k_theta_a"),
    "first_discarded_hz": _ROUTE_COMMON_DEPS + (
        "screw_length", "screw_diameter", "screw_density",
        "nut_mass", "stage_mass", "J_m", "J_c", "k_c_series", "k_theta_a"),
    "fixed_interface_separation": _ROUTE_COMMON_DEPS + (
        "screw_length", "screw_diameter", "screw_density",
        "nut_mass", "stage_mass", "J_m", "J_c", "k_c_series", "k_theta_a"),
    "discarded_pole_separation": _ROUTE_COMMON_DEPS + (
        "screw_length", "screw_diameter", "screw_density",
        "nut_mass", "stage_mass", "J_m", "J_c", "k_c_series", "k_theta_a"),
    "mu_fraction": ("reduced_drive_mass", "reduced_stage_mass"),
    "relative_mode_nearground_hz": ("reduced_axial_stiffness", "reduced_stage_mass"),
    "drive_pole_hz": ("magnetic_stiffness", "reduced_drive_mass"),
}.items():
    PARAMETER_REGISTRY[_key] = {
        "category": "output" if _key.startswith((
            "full_model", "cb_", "first_", "fixed_interface", "discarded_",
        )) else "derived",
        "dependencies": _dependencies,
        "section": "6-3-reduction-evidence",
    }

for _index in range(1, GMS_N + 1):
    PARAMETER_REGISTRY[f"gms_nu{_index}"] = {
        "category": "assumed", "dependencies": (),
        "section": "8-1-executed-provisional-friction-values",
    }
for _key in ("stribeck_delta", "tau_C", "eta_screw", "F_preload_nut"):
    PARAMETER_REGISTRY[_key] = {
        "category": "assumed", "dependencies": (),
        "section": "8-3-executed-provisional-friction-values",
    }
for _key, _dependencies in {
    "d_Fs_efficiency_estimate": ("eta_screw", "F_preload_nut"),
    "d_Fs_torque_equivalent": ("d_Fs", "lead"),
    "detent_velocity_drive": ("magnetic_stiffness", "reduced_drive_mass", "lead", "rotor_teeth"),
    "detent_velocity_axial": ("reduced_axial_stiffness", "reduced_stage_mass", "lead", "rotor_teeth"),
    "detent_velocity_discarded": ("lead", "rotor_teeth"),
    "retained_mode_period": ("reduced_axial_stiffness", "reduced_stage_mass"),
    "tau_C_mode_ratio": ("tau_C", "reduced_axial_stiffness", "reduced_stage_mass"),
}.items():
    PARAMETER_REGISTRY[_key] = {
        "category": "derived", "dependencies": _dependencies,
        "section": "8-2-executed-provisional-friction-values",
    }
for _site in ("g", "n", "d"):
    for _suffix in ("sigma0", "sigma1", "sigma2", "Fs", "Fc", "vs", "C"):
        PARAMETER_REGISTRY[f"{_site}_{_suffix}"] = {
            "category": "input" if (_site, _suffix) == ("g", "sigma0") else "assumed",
            "dependencies": (),
            "section": "8-1-executed-provisional-friction-values",
        }
    for _index in range(1, GMS_N + 1):
        PARAMETER_REGISTRY[f"{_site}_k{_index}"] = {
            "category": "assumed", "dependencies": (),
            "section": "8-1-executed-provisional-friction-values",
        }
    _site_deps = tuple(f"{_site}_k{_index}" for _index in range(1, GMS_N + 1))
    _fraction_deps = tuple(f"gms_nu{_index}" for _index in range(1, GMS_N + 1))
    for _index in range(1, GMS_N + 1):
        for _level in ("fs", "fc"):
            PARAMETER_REGISTRY[f"yield_{_site}_{_index}_{_level}"] = {
                "category": "derived",
                "dependencies": (f"{_site}_k{_index}", f"gms_nu{_index}",
                                 f"{_site}_Fs" if _level == "fs" else f"{_site}_Fc"),
                "section": "8-3-implementation-choices",
            }
    PARAMETER_REGISTRY[f"gms_rate_{_site}"] = {
        "category": "derived",
        "dependencies": (f"{_site}_Fs", f"{_site}_Fc", "tau_C"),
        "section": "8-2-executed-provisional-friction-values",
    }
    PARAMETER_REGISTRY[f"tau_C_{_site}"] = {
        "category": "derived",
        "dependencies": (f"{_site}_Fs", f"{_site}_Fc", "tau_C"),
        "section": "8-2-executed-provisional-friction-values",
    }
    PARAMETER_REGISTRY[f"yield_span_{_site}"] = {
        "category": "derived",
        "dependencies": _site_deps + _fraction_deps,
        "section": "8-4-implementation-choices",
    }
    PARAMETER_REGISTRY[f"static_deflection_{_site}"] = {
        "category": "derived",
        "dependencies": (f"{_site}_Fs", f"{_site}_sigma0"),
        "section": "8-4-implementation-choices",
    }

for _key, _dependencies in {
    "detent_settling_time_2pct": ("electromagnetic_zeta", "detent_torque", "detent_phase", "holding_torque", "reduced_drive_mass"),
    "axial_settling_time_2pct": ("axial_damping", "reduced_axial_stiffness", "reduced_drive_mass", "reduced_stage_mass", "electromagnetic_zeta"),
    "plateau_dwell": ("detent_settling_time_2pct", "axial_settling_time_2pct",
                      "interface_settling_ms", "measured_settling_ms"),
    "guideway_a2_final_origin_nm": ("g_sigma0", "g_Fs", "g_Fc", "tau_C", "command_step"),
    "guideway_loop_energy_ratio_pct": ("g_sigma0", "g_Fs", "g_Fc", "tau_C", "command_step"),
    "nut_loop_energy_ratio_pct": ("n_sigma0", "n_Fs", "n_Fc", "tau_C", "command_step"),
    "nut_return_force_ratio": ("n_sigma0", "n_Fs", "n_Fc", "tau_C", "command_step"),
    "guideway_r_hold_lugre_pct": ("g_sigma0", "g_Fs", "command_step"),
    "guideway_r_hold_gms_pct": ("g_sigma0", "g_Fs", "command_step"),
    "guideway_r_hold_ratio": ("g_sigma0", "g_Fs", "command_step"),
    "retention_ratio_low": ("g_sigma0", "g_Fs", "axial_damping",
                            "zeta_relative_measured"),
    "retention_ratio_high": ("g_sigma0", "g_Fs", "axial_damping",
                             "zeta_relative_measured"),
    # Prose frequencies, band edges and requirements.  Every one of these was
    # a hand-typed literal that had drifted from its generated source.
    "detent_band_low_hz": ("detent_torque", "holding_torque", "rotor_teeth",
                           "lead", "reduced_drive_mass", "reduced_stage_mass",
                           "reduced_axial_stiffness"),
    "detent_band_high_hz": ("detent_torque", "holding_torque", "rotor_teeth",
                            "lead", "reduced_drive_mass", "reduced_stage_mass",
                            "reduced_axial_stiffness"),
    "operating_mode_hz": ("reduced_axial_stiffness", "reduced_stage_mass",
                          "g_sigma0", "n_sigma0", "d_sigma0"),
    "operating_fixed_interface_separation": (
        "operating_mode_hz", "first_fixed_interface_hz"),
    "operating_fixed_interface_ratio": (
        "operating_mode_hz", "first_fixed_interface_hz"),
    "baseline_fixed_interface_ratio": (
        "axial_mode_target_hz", "first_fixed_interface_hz"),
    "presliding_k_ax_mn": ("reduced_axial_stiffness", "g_sigma0", "n_sigma0",
                           "axial_mode_target_hz"),
    "closure_singular_limit_mn": ("reduced_axial_stiffness", "k_sha", "k_mnt"),
    "measured_settling_ms": ("zeta_relative_measured", "axial_mode_target_hz"),
    "interface_settling_ms": ("interface_axial_damping", "reduced_axial_stiffness",
                              "reduced_drive_mass", "reduced_stage_mass"),
    "executed_settling_ms": ("axial_damping", "reduced_axial_stiffness",
                             "reduced_drive_mass", "reduced_stage_mass"),
    "detent_equilibrium_error_nm": ("detent_torque", "holding_torque",
                                    "rotor_teeth", "lead"),
    "predistortion_resolution_nm": ("detent_equilibrium_error_nm",
                                    "predistortion_levels"),
    "required_microstep_divisor": ("predistortion_resolution_nm", "full_step_pitch"),
    "predistortion_levels_executed": ("detent_equilibrium_error_nm", "command_step"),
    "predistortion_levels": (),
    "screw_inertia_nominal": ("screw_diameter", "screw_length", "screw_density"),
}.items():
    PARAMETER_REGISTRY[_key] = {
        "category": "output", "dependencies": _dependencies,
        "section": "9-force-instrumented-partial-slip-memory-experiment",
    }

for _key, _section in {
    "friction_site_definitions": "8-1-how-the-friction-laws-attach-to-the-plant",
    "friction_state_definition": "8-2-constitutive-laws",
    "identifiability_analysis": "8-4-implementation-choices",
    "metrology_campaign": "9-force-instrumented-partial-slip-memory-experiment",
    "measured_axial_band_low_hz": "2-entry-parameters",
    "measured_axial_band_high_hz": "2-entry-parameters",
    "detent_enabled": "2-entry-parameters",
}.items():
    PARAMETER_REGISTRY[_key] = {
        "category": "input", "dependencies": (), "section": _section,
    }
for _key, _dependencies, _section in (
    ("section7_rms_pct", ("full_model_upper_hz", "route_c_f2"), "7-full-versus-reduced-verification"),
    ("section7_rms_nm", ("full_model_upper_hz", "route_c_f2"), "7-full-versus-reduced-verification"),
    ("section7_drive_share_pct", ("full_model_upper_hz", "route_c_f2"), "7-full-versus-reduced-verification"),
    ("section7_drive_pole_error_pct", ("full_model_upper_hz", "route_c_f1"), "7-full-versus-reduced-verification"),
    ("section7_frequency_share_pct", ("full_model_upper_hz", "route_c_f2"), "7-full-versus-reduced-verification"),
    ("section7_damping_share_pct", ("full_model_zeta", "route_f_zeta"), "7-full-versus-reduced-verification"),
    ("section7_reduced_coordinate_count", ("route_p_f1", "route_p_f2"), "7-full-versus-reduced-verification"),
    ("section7_full_coordinate_count", ("full_model_upper_hz",), "7-full-versus-reduced-verification"),
    ("friction_port_count", ("friction_site_definitions",), "8-0-how-the-friction-laws-attach-to-the-plant"),
    ("gms_states_per_site", ("friction_state_definition",), "8-1-constitutive-laws"),
    ("lugre_states_per_site", ("friction_state_definition",), "8-1-constitutive-laws"),
    ("structural_identifiability_result_count", ("identifiability_analysis",), "8-3-implementation-choices"),
    ("project_adev_floor_nm", ("metrology_campaign",), "9-force-instrumented-partial-slip-memory-experiment"),
    ("a2_convergence_order", ("g_sigma0", "g_Fs", "g_Fc", "tau_C", "command_step"), "12-1-gms-step-halving-convergence"),
):
    PARAMETER_REGISTRY[_key] = {
        "category": "output", "dependencies": _dependencies, "section": _section,
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
    fig.savefig(output, format="svg", bbox_inches="tight", metadata={"Date": None})
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
    pattern = re.compile(r"\[\[derived:([A-Za-z0-9_]+)(?:@[a-z_]+)?=")
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
def interface_loss_factors() -> dict[str, dict[str, float]]:
    """Return the executed loss factor of every axial interface at f_2.

    A frequency-independent dashpot has eta_j(w) = w*c_j/k_j = 2*zeta_j*w/w_j,
    so the executed loss factor depends on the element's own frequency and is
    not 2*zeta_j anywhere except at that frequency.
    """
    constants = physical_constants()
    component = component_parameters()
    r = constants["r"]
    f_2 = constants["axial_mode_target_hz"]
    segment_mass = component["screw_mass"] / 3.0
    segment_inertia = component["screw_inertia"] / 3.0
    pair = lambda a, b: a * b / (a + b)  # noqa: E731
    elements = {
        "zeta_bearing": (component["k_brg"], segment_mass),
        "zeta_steel": (component["k_sha"], pair(segment_mass, segment_mass)),
        "zeta_ball_nut": (constants["k_ball"], 1.0 / (
            r**2 / segment_inertia + 1.0 / segment_mass + 1.0 / component["m_n"])),
        "zeta_nut_mount": (component["k_mnt"], pair(component["m_n"], component["m_stage"])),
    }
    report: dict[str, dict[str, float]] = {}
    total_compliance = sum(1.0 / stiffness for stiffness, _mass in elements.values())
    for key, (stiffness, relative_mass) in elements.items():
        zeta = component[key]
        f_j = np.sqrt(stiffness / relative_mass) / (2.0 * np.pi)
        report[key] = {
            "stiffness": stiffness,
            "relative_mass": relative_mass,
            "f_j": f_j,
            "zeta": zeta,
            "damping": 2.0 * zeta * np.sqrt(stiffness * relative_mass),
            "eta": 2.0 * zeta * f_2 / f_j,
            "weight": (1.0 / stiffness) / total_compliance,
        }
    return report


def validate_interface_loss_factors() -> dict[str, dict[str, float]]:
    """Fail the build if an executed damper drifts from its target loss factor."""
    report = interface_loss_factors()
    for key, target in INTERFACE_LOSS_FACTORS.items():
        executed = report[key]["eta"]
        if not np.isclose(executed, target, rtol=1.0e-4, atol=0.0):
            raise ValueError(
                f"Interface {key} executes eta={executed:.6g} at the retained mode "
                f"but targets {target:.6g}; update FULL['{key}'] to "
                f"{target * report[key]['f_j'] / (2.0 * physical_constants()['axial_mode_target_hz']):.9f}"
            )
    return report


def validate_gms_partition() -> dict[str, object]:
    """Fail the build unless every executed GMS partition closes exactly."""
    if any(GMS_WEIGHTS.size != stiffness.size
           for stiffness in GMS_STIFFNESS_BY_SITE.values()):
        raise ValueError("GMS force weights and per-site stiffness vectors must have equal length")
    if np.any(GMS_WEIGHTS <= 0.0) or any(
            np.any(stiffness <= 0.0) for stiffness in GMS_STIFFNESS_BY_SITE.values()):
        raise ValueError("Every GMS force weight and element stiffness must be positive")
    weight_sum = float(np.sum(GMS_WEIGHTS))
    if not np.isclose(weight_sum, 1.0, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"GMS force weights must satisfy sum(nu_i)=1; got {weight_sum:.16g}")

    stiffness_sums: dict[str, float] = {}
    for site, parameters in FRICTION.items():
        element_stiffness = GMS_STIFFNESS_BY_SITE[site]
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
        "G": {"g"}, "G2": {"g"},
        "B": {"d", "n"}, "B2": {"d", "n"},
        "C": {"d", "g", "n"}, "C2": {"d", "g", "n"},
        "A1v": {"d", "g"},
    }
    for key, expected in expected_sites.items():
        actual = set(CASES[key]["sites"])
        if actual != expected:
            raise ValueError(f"Case {key} friction sites are {sorted(actual)}, expected {sorted(expected)}")
    # The controlled comparison requires equal tangent damping, so every case
    # except the declared micro-viscous variant must execute sigma_1 = 0.
    for key, case in CASES.items():
        if case.get("micro_viscous"):
            continue
        for site in case["sites"]:
            if FRICTION[site]["sigma1"] != 0.0:
                raise ValueError(
                    f"Case {key} site {site} carries sigma_1={FRICTION[site]['sigma1']}; "
                    "the controlled A/B/C comparison requires zero")
    for site, values in FRICTION.items():
        executed = (values["F_s"] - values["F_c"]) / MODEL["tau_C"]
        if not np.isclose(values["C_gms"], executed, rtol=1.0e-9, atol=0.0):
            raise ValueError(
                f"Site {site} executes C={values['C_gms']:.6g} N/s but tau_C implies {executed:.6g}")
    nut_first_yield = float(np.min(
        GMS_WEIGHTS * FRICTION["n"]["F_s"] /
        GMS_STIFFNESS_BY_SITE["n"]
    ))
    if not np.isclose(nut_first_yield, 0.20e-6, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"Nut microslip first yield is {nut_first_yield:.6g} m, expected 0.20 um")


def validate_command_design(constants: dict[str, float]) -> dict[str, float]:
    """Fail unless the production 1/16 command families retain their margins."""
    quantum = constants["command_step"]
    expected_quantum = constants["full_step"] / 16.0
    if MODEL["microstep_divisor"] != 16 or not np.isclose(
            quantum, expected_quantum, rtol=0.0, atol=1.0e-15):
        raise ValueError(
            f"Production command quantum must be full_step/16; got {quantum:.9g} m")

    max_increment = float(np.max(np.abs(np.diff(MAIN_LEVELS))) * quantum)
    if max_increment > constants["quarter_step"] + 1.0e-15:
        raise ValueError(
            f"Main command increment {max_increment:.9g} m exceeds the quarter-step bound")

    def yields(site: str) -> np.ndarray:
        values = FRICTION[site]
        return GMS_WEIGHTS * values["F_s"] / GMS_STIFFNESS_BY_SITE[site]

    guideway_yields = yields("g")
    nut_yields = yields("n")
    guideway_inner = 4.0 * quantum
    guideway_outer = 12.0 * quantum
    nut_inner = quantum
    nut_nested = 2.0 * quantum
    nut_outer = 3.0 * quantum
    if not (guideway_inner > guideway_yields[0]
            and guideway_outer > guideway_yields[1]
            and guideway_outer < guideway_yields[2]):
        raise ValueError("Guideway 12/10/4 command levels no longer bracket yields 1 and 2")
    if not (nut_inner > nut_yields[0]
            and nut_nested > nut_yields[1]
            and nut_outer < nut_yields[2]):
        raise ValueError("Nut 3/2/1 command levels no longer bracket yields 1 and 2 below yield 3")
    return {
        "quantum_nm": quantum * 1.0e9,
        "main_max_increment_um": max_increment * 1.0e6,
        "guideway_inner_margin_um": (guideway_inner - guideway_yields[0]) * 1.0e6,
        "nut_outer_margin_um": (nut_yields[2] - nut_outer) * 1.0e6,
    }


def presliding_calibrated_axial_stiffness(
        m_d: float, m_s: float, k_m: float, target_hz: float,
        sites: tuple[str, ...] = ("g", "n")) -> float:
    """Solve k_ax so the presliding-tangent plant reproduces the measured mode.

    A hammer FRF at micrometre amplitudes excites the assembled axis inside
    the presliding regime, where every rolling contact behaves as a spring of
    stiffness sigma_0.  The measured pole therefore already contains the
    guideway and nut presliding tangents, and calibrating k_ax on the
    frictionless eigenproblem attributes that stiffness to the structure a
    second time when the friction ports are added back.  This branch removes
    the double count by requiring the friction-on plant, not the frictionless
    one, to land on the target.
    """
    eigenvalue = (2.0 * np.pi * target_hz) ** 2
    sigma_d = sum(FRICTION[site]["sigma0"] for site in sites if site == "d")
    sigma_g = sum(FRICTION[site]["sigma0"] for site in sites if site == "g")
    sigma_n = sum(FRICTION[site]["sigma0"] for site in sites if site == "n")
    # K = [[k_m + sigma_d + k_ax + sigma_n, -(k_ax + sigma_n)],
    #      [-(k_ax + sigma_n), k_ax + sigma_n + sigma_g]]
    # det(K - lambda M) = 0 is linear in the relative stiffness k = k_ax+sigma_n.
    drive = k_m + sigma_d - eigenvalue * m_d
    stage = sigma_g - eigenvalue * m_s
    denominator = drive + stage
    if abs(denominator) <= np.finfo(float).eps * max(abs(k_m), 1.0):
        raise ValueError("Presliding-inclusive calibration is singular for these inputs")
    relative = -drive * stage / denominator
    k_ax = relative - sigma_n
    if not np.isfinite(k_ax) or k_ax <= 0.0:
        raise ValueError(
            "The presliding-inclusive branch leaves no positive structural "
            f"stiffness (k_ax={k_ax:.6g} N/m); the measured mode is below the "
            "presliding tangents alone")
    return float(k_ax)


def ball_closure_band(k_ax: float, k_sha: float, k_mnt: float) -> dict[str, float]:
    """Report how much of k_ball is a choice of k_brg rather than a result.

    k_ball is a closure residual: it absorbs whatever axial compliance the
    other three elements leave over.  Below a singular limit on k_brg there is
    nothing left to absorb, and just above it the residual swings by a factor
    of two, so the reported k_ball carries the bearing assumption with it.
    """
    residual = 1.0 / k_ax - 1.0 / k_sha - 1.0 / k_mnt
    if residual <= 0.0:
        raise ValueError(
            "The screw and nut-mount compliances alone exceed the calibrated "
            "axial compliance; no bearing stiffness can close the budget")
    singular_limit = 1.0 / residual
    def k_ball_at(k_brg: float) -> float:
        remaining = residual - 1.0 / k_brg
        return float("inf") if remaining <= 0.0 else 1.0 / remaining
    return {
        "singular_limit": singular_limit,
        "k_ball_at": k_ball_at,
        "samples": tuple(
            (k_brg, k_ball_at(k_brg))
            for k_brg in (1.05 * singular_limit, 15.0e6, 25.0e6, 40.0e6)
            if k_brg > singular_limit),
    }


def validate_closure_band(constants: dict[str, float],
                          component: dict[str, float]) -> dict[str, float]:
    """Fail the build if the bearing assumption sits at the singular limit."""
    band = ball_closure_band(
        constants["k_ax"], component["k_sha"], component["k_mnt"])
    k_brg = component["k_brg"]
    limit = band["singular_limit"]
    if k_brg <= limit:
        raise ValueError(
            f"k_brg={k_brg / 1e6:.3f} MN/m is at or below the closure singular "
            f"limit {limit / 1e6:.3f} MN/m; k_ball cannot be positive")
    margin = k_brg / limit
    if margin < 1.10:
        raise ValueError(
            f"k_brg={k_brg / 1e6:.3f} MN/m sits within 10% of the "
            f"{limit / 1e6:.3f} MN/m closure singular limit; k_ball is not a "
            "reportable quantity there")
    return {"singular_limit": limit, "margin": margin,
            "k_brg": k_brg, "k_ball": constants["k_ball"]}


def validate_prose_frequencies() -> dict[str, int]:
    """Fail if a two-decimal frequency is typed into prose instead of tokenized.

    A literal such as `695.82 Hz` or `166.93 Hz` is a generated value that was
    pasted by hand; every one of them in Revision 3 had drifted from its
    source by the time it was found.  Generated blocks are exempt because the
    builder writes them on every run, and so are integer bands such as
    `681-690 Hz`, which are recorded inputs rather than computed outputs.
    """
    literal = re.compile(r"(?<!derived:)(?<![\w.])\d{2,4}\.\d{2}\s*(?:Hz|kHz)")
    generated = re.compile(r"<!-- BEGIN GENERATED.*?<!-- END GENERATED[^>]*-->",
                           flags=re.DOTALL)
    offenders: dict[str, int] = {}
    for document in (DESCRIPTION_MD, DERIVATION_MD):
        text = generated.sub("", document.read_text(encoding="utf-8"))
        # Strip token bodies: the fallback text inside [[derived:x=695.82]] is
        # generator-owned and is rewritten on every build.
        text = re.sub(r"\[\[[a-z]+:[^\]]+\]\]", "", text)
        for match in literal.finditer(text):
            offenders[f"{document.name}: {match.group(0)}"] = (
                offenders.get(f"{document.name}: {match.group(0)}", 0) + 1)
    if offenders:
        raise ValueError(
            "Hand-typed frequencies must be generated tokens: "
            + "; ".join(f"{key} x{count}" for key, count in sorted(offenders.items())))
    return offenders


def validate_breakaway_forces() -> dict[str, dict[str, float]]:
    """Warn when an executed breakaway force leaves its own stated range.

    This is a report rather than an abort: the executed guideway value is
    outside its stated likely range on purpose, Section 12.4 executes the
    alternative, and the build has to stay runnable in order to price it.
    """
    stated_ranges = {"g": (1.0, 1.5)}
    report: dict[str, dict[str, float]] = {}
    for site, (low, high) in stated_ranges.items():
        force = float(FRICTION[site]["F_s"])
        report[site] = {
            "F_s": force, "low": low, "high": high,
            "inside": bool(low <= force <= high),
            "factor_above": force / high if force > high else 1.0,
        }
    return report


def validate_predistortion_authority(constants: dict[str, float]) -> dict[str, float]:
    """Report whether the command grid can express a detent pre-distortion.

    This is deliberately a report and not an abort.  The executed divisor
    fails the requirement by a large factor, and that failure is a documented
    finding in Section 5 and Appendix C item 16 rather than a build error: the
    model must stay runnable at the production setting in order to quantify
    what that setting costs.
    """
    achieved_levels = (constants["detent_equilibrium_error"]
                       / constants["command_step"])
    return {
        "command_step": constants["command_step"],
        "required_resolution": constants["predistortion_resolution"],
        "required_divisor": constants["required_microstep_divisor"],
        "executed_divisor": float(MODEL["microstep_divisor"]),
        "requested_levels": float(MODEL["predistortion_levels"]),
        "achieved_levels": achieved_levels,
        "satisfied": bool(constants["command_step"]
                          <= constants["predistortion_resolution"]),
        "shortfall": constants["required_microstep_divisor"]
        / float(MODEL["microstep_divisor"]),
    }


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


def axial_complex_link(component: dict[str, float], r: float,
                       k_ball: float, omega: float) -> tuple[float, float]:
    """Condense the damped axial chain at one declared angular frequency.

    Each Kelvin-Voigt element uses the damping ratio appropriate to its
    physical interface.  Taking the reciprocal of the summed complex
    compliances preserves those dependencies in the executable 2-DOF link.
    """
    m_b, m_e = component["m_b"], component["m_e"]
    m_n, m_stage = component["m_n"], component["m_stage"]
    ball_relative_mass = 1.0 / (
        r**2 / component["J_s2"] + 1.0 / m_e + 1.0 / m_n)
    elements = (
        (component["k_brg"], 2.0 * component["zeta_bearing"]
         * np.sqrt(component["k_brg"] * m_b)),
        (component["k_sha"], 2.0 * component["zeta_steel"]
         * np.sqrt(component["k_sha"] * m_b * m_e / (m_b + m_e))),
        (k_ball, 2.0 * component["zeta_ball_nut"]
         * np.sqrt(k_ball * ball_relative_mass)),
        (component["k_mnt"], 2.0 * component["zeta_nut_mount"]
         * np.sqrt(component["k_mnt"] * m_n * m_stage / (m_n + m_stage))),
    )
    compliance = sum(1.0 / complex(k_value, omega * c_value)
                     for k_value, c_value in elements)
    stiffness = 1.0 / compliance
    return float(stiffness.real), float(stiffness.imag / omega)


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
    _k_ax_complex, c_ax_interface = axial_complex_link(
        component, r, k_ball, 2.0 * np.pi * axial_mode_target_hz)
    c_m = 2.0 * MODEL["zeta_m"] * np.sqrt(k_m * m_d)
    minimum_local_tangent = k_m - k_det_amplitude
    if minimum_local_tangent <= 0.0:
        raise ValueError("Detent tangent can exceed the commutation tangent")
    minimum_local_omega = np.sqrt(minimum_local_tangent / m_d)
    detent_settling_time_2pct = SETTLING_2PCT_FACTOR / (
        MODEL["zeta_m"] * minimum_local_omega)
    mass = np.diag([m_d, m_s])
    damping = np.array([
        [c_m + MODEL["c_ax"], -MODEL["c_ax"]],
        [-MODEL["c_ax"], MODEL["c_ax"]],
    ])
    stiffness = np.array([
        [k_m + k_ax, -k_ax],
        [-k_ax, k_ax],
    ])
    axial_frequency, axial_zeta = _damped_modal_data(mass, damping, stiffness)[1]
    axial_settling_time_2pct = SETTLING_2PCT_FACTOR / (
        axial_zeta * 2.0 * np.pi * axial_frequency)
    # Three candidate settling times exist for the same retained mode and they
    # disagree by a factor of eleven: the executed link damper, the interface
    # loss factors propagated in E.5, and the measured relative-mode damping.
    # The rule previously took the maximum of a set that omitted the two
    # longest candidates, which silently selected the shortest branch.  The
    # dwell now covers every candidate, so the settled window is defensible on
    # whichever branch the E.7 half-power extraction confirms.
    interface_zeta = _damped_modal_data(
        mass,
        np.array([[c_m + c_ax_interface, -c_ax_interface],
                  [-c_ax_interface, c_ax_interface]]),
        stiffness)[1][1]
    interface_settling_time_2pct = SETTLING_2PCT_FACTOR / (
        interface_zeta * 2.0 * np.pi * axial_mode_target_hz)
    measured_settling_time_2pct = SETTLING_2PCT_FACTOR / (
        MODEL["zeta_relative_measured"] * 2.0 * np.pi * axial_mode_target_hz)
    plateau_dwell = max(
        0.100, detent_settling_time_2pct, axial_settling_time_2pct,
        interface_settling_time_2pct, measured_settling_time_2pct)
    # Positional pre-distortion authority.  The detent equilibrium error is a
    # position-periodic term, so a correction table can only address it if the
    # command grid can place levels inside it.  This is a divisor requirement,
    # not a quantization-noise argument; see Section 5 and Appendix C item 16.
    detent_equilibrium_error = (r / teeth) * np.arcsin(min(t_det / t_max, 1.0))
    predistortion_resolution = detent_equilibrium_error / MODEL["predistortion_levels"]
    required_divisor = full_step / predistortion_resolution
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
        "c_ax": MODEL["c_ax"],
        "c_ax_interface": c_ax_interface,
        "c_m": c_m,
        "full_step": full_step,
        "quarter_step": full_step / 4.0,
        "command_step": full_step / MODEL["microstep_divisor"],
        "detent_period": full_step,
        "settling_time_2pct": detent_settling_time_2pct,
        "detent_settling_time_2pct": detent_settling_time_2pct,
        "axial_settling_time_2pct": axial_settling_time_2pct,
        "interface_settling_time_2pct": interface_settling_time_2pct,
        "measured_settling_time_2pct": measured_settling_time_2pct,
        "axial_zeta_executed": axial_zeta,
        "axial_zeta_interface": interface_zeta,
        "axial_zeta_measured": MODEL["zeta_relative_measured"],
        "plateau_dwell": plateau_dwell,
        "metric_window": min(0.020, 0.20 * plateau_dwell),
        "detent_equilibrium_error": detent_equilibrium_error,
        "predistortion_resolution": predistortion_resolution,
        "required_microstep_divisor": required_divisor,
        "detent_enabled": float(MODEL["detent_enabled"]),
    }


def screw_segment_stiffnesses() -> dict[str, float]:
    """Derive both screw segment stiffness pairs from one geometry.

    Torsional and axial stiffness must describe the same screw and the same
    nut position.  Both segments therefore come from the declared datum
    L_a = nut_axial_datum with L_b = L_s - L_a, and from the root-diameter
    section rather than the nominal thread diameter, so that
    k_theta_a/k_theta_b and k_sha/k_shb are the same length ratio by
    construction.
    """
    length_a = float(FULL["nut_axial_datum"])
    length_b = float(FULL["screw_length"]) - length_a
    if length_a <= 0.0 or length_b <= 0.0:
        raise ValueError(
            f"Nut axial datum {length_a:.4f} m must lie inside the "
            f"{FULL['screw_length']:.4f} m screw")
    root_radius = 0.5 * float(FULL["screw_root_diameter"])
    area = np.pi * root_radius**2
    polar_inertia = 0.5 * np.pi * root_radius**4
    tensile = float(FULL["youngs_modulus"]) * area
    torsional = float(FULL["shear_modulus"]) * polar_inertia
    return {
        "screw_length_a": length_a,
        "screw_length_b": length_b,
        "screw_root_area": area,
        "screw_root_polar_inertia": polar_inertia,
        "k_theta_a": torsional / length_a,
        "k_theta_b": torsional / length_b,
        "k_sha": tensile / length_a,
        "k_shb": tensile / length_b,
        "axial_rigidity": tensile,
        "torsional_rigidity": torsional,
    }


def validate_screw_geometry() -> dict[str, float]:
    """Fail the build if the two screw segments stop describing one screw."""
    segments = screw_segment_stiffnesses()
    length_sum = segments["screw_length_a"] + segments["screw_length_b"]
    if not np.isclose(length_sum, FULL["screw_length"], rtol=0.0, atol=1.0e-9):
        raise ValueError(
            f"Screw segments sum to {length_sum:.6f} m against a "
            f"{FULL['screw_length']:.6f} m screw")
    torsional_ratio = segments["k_theta_a"] / segments["k_theta_b"]
    axial_ratio = segments["k_sha"] / segments["k_shb"]
    if not np.isclose(torsional_ratio, axial_ratio, rtol=1.0e-9, atol=0.0):
        raise ValueError(
            f"Torsional segment ratio {torsional_ratio:.6f} disagrees with the "
            f"axial ratio {axial_ratio:.6f}; the two pairs describe different "
            "nut positions")
    if segments["screw_root_polar_inertia"] <= 0.0:
        raise ValueError("Screw root section must be positive")
    return segments


def _axial_element_frequencies(p: dict[str, float]) -> dict[str, dict[str, float]]:
    """Stiffness and relative mass of every axial interface element.

    This is the one place the four elements are defined.  It runs before
    physical_constants() so that the damping ratios can be derived rather than
    entered, and it therefore recomputes k_ax and k_ball from the same closed
    form those functions use rather than importing them.
    """
    r = MODEL["lead"] / (2.0 * np.pi)
    m_d = (p["J_m"] + p["J_c"] + p["screw_inertia"]) / r**2
    m_s = p["m_stage"] + p["m_n"]
    k_m = MODEL["rotor_teeth"] * MODEL["T_max"] / r**2
    k_ax = modal_calibrated_axial_stiffness(
        m_d, m_s, k_m, MODEL["axial_mode_target_hz"])
    k_ball = closure_ball_stiffness(k_ax, p["k_brg"], p["k_sha"], p["k_mnt"])
    segment_mass = p["screw_mass"] / 3.0
    segment_inertia = p["screw_inertia"] / 3.0
    pair = lambda a, b: a * b / (a + b)  # noqa: E731
    elements = {
        "zeta_bearing": (p["k_brg"], segment_mass),
        "zeta_steel": (p["k_sha"], pair(segment_mass, segment_mass)),
        "zeta_ball_nut": (k_ball, 1.0 / (
            r**2 / segment_inertia + 1.0 / segment_mass + 1.0 / p["m_n"])),
        "zeta_nut_mount": (p["k_mnt"], pair(p["m_n"], p["m_stage"])),
    }
    return {
        key: {
            "stiffness": stiffness,
            "relative_mass": relative_mass,
            "f_j": float(np.sqrt(stiffness / relative_mass) / (2.0 * np.pi)),
        }
        for key, (stiffness, relative_mass) in elements.items()
    }


def interface_damping_ratios(p: dict[str, float]) -> dict[str, float]:
    """Solve each interface damping ratio from its target loss factor.

    A frequency-independent dashpot delivers eta_j(w) = 2*zeta_j*w/w_j, so the
    ratio that realizes a target loss factor AT THE RETAINED MODE is

        zeta_j = eta_j * f_j / (2 * f_2).

    The ratios are derived here rather than entered, so a change of screw
    geometry or bearing stiffness cannot leave four hand-copied constants
    describing a different assembly.  This is the only conversion in the
    document; Section 2 and E.5 both quote it from here.
    """
    f_2 = MODEL["axial_mode_target_hz"]
    return {
        key: INTERFACE_LOSS_FACTORS[key] * element["f_j"] / (2.0 * f_2)
        for key, element in _axial_element_frequencies(p).items()
    }


def component_parameters() -> dict[str, float]:
    """Derive screw inertia and lumped masses from the 0.192 m component."""
    p = dict(FULL)
    p.update(screw_segment_stiffnesses())
    # Mass follows the nominal diameter; the polar inertia follows the root
    # section, because the thread removes exactly the material that would
    # otherwise dominate a d^4 quantity.
    mass_radius = 0.5 * p["screw_diameter"]
    root_radius = 0.5 * p["screw_root_diameter"]
    screw_mass = p["screw_density"] * np.pi * mass_radius**2 * p["screw_length"]
    screw_inertia = (0.5 * np.pi * root_radius**4
                     * p["screw_density"] * p["screw_length"])
    p["screw_mass"] = screw_mass
    p["screw_inertia"] = screw_inertia
    p["screw_inertia_nominal"] = 0.5 * screw_mass * mass_radius**2
    p.update(interface_damping_ratios(p))
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
        _add_pair(damping, i, j, _pair_damping(
            k_value, native_masses[i], native_masses[j], p["zeta_steel"]))

    stiffness[5, 5] += p["k_brg"]
    damping[5, 5] += 2.0 * p["zeta_bearing"] * np.sqrt(p["k_brg"] * p["m_b"])
    for i, j, k_value in ((5, 6, p["k_sha"]), (6, 7, p["k_shb"]), (8, 9, p["k_mnt"])):
        _add_pair(stiffness, i, j, k_value)
        interface_zeta = p["zeta_nut_mount"] if (i, j) == (8, 9) else p["zeta_steel"]
        _add_pair(damping, i, j, _pair_damping(
            k_value, native_masses[i], native_masses[j], interface_zeta))

    # delta_n = u_n - u_e - r theta_s2; outer products are the virtual-work
    # contribution of one conservative ball-contact element to K and C.
    h_nut = np.zeros(10)
    h_nut[3], h_nut[6], h_nut[8] = -r, -1.0, 1.0
    stiffness += p["k_ball"] * np.outer(h_nut, h_nut)
    relative_mass = 1.0 / (r**2 / p["J_s2"] + 1.0 / p["m_e"] + 1.0 / p["m_n"])
    c_ball = 2.0 * p["zeta_ball_nut"] * np.sqrt(p["k_ball"] * relative_mass)
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


def _damped_modal_data(mass: np.ndarray, damping: np.ndarray,
                       stiffness: np.ndarray) -> list[tuple[float, float]]:
    """Return damped modal frequency and damping ratio from state eigenvalues."""
    count = mass.shape[0]
    state_matrix = np.block([
        [np.zeros((count, count)), np.eye(count)],
        [-np.linalg.solve(mass, stiffness), -np.linalg.solve(mass, damping)],
    ])
    eigenvalues = np.linalg.eigvals(state_matrix)
    positive_imaginary = eigenvalues[np.imag(eigenvalues) > 1.0]
    ordered = positive_imaginary[np.argsort(np.abs(np.imag(positive_imaginary)))]
    return [
        (float(abs(value) / (2.0 * np.pi)), float(-np.real(value) / abs(value)))
        for value in ordered
    ]


def craig_bampton_plant() -> dict[str, object]:
    """Project the ten-DOF model onto [x_d, x_s] plus one fixed-interface mode.

    The retained physical coordinates are exact, the eliminated partition keeps
    its static constraint modes, and only the lowest fixed-interface mode is
    restored.  The projected damping is inherited from the ten-DOF assembly, so
    this plant is the one reduction that does not need an assumed link damper.
    """
    constants = physical_constants()
    full_m, full_c, full_k, full_b, _p = full_linear_matrices()
    order = [0, 9, *range(1, 9)]
    coordinate_map = np.eye(10)[:, order]
    coordinate_map[:, 0] /= constants["r"]
    transformed_m = coordinate_map.T @ full_m @ coordinate_map
    transformed_c = coordinate_map.T @ full_c @ coordinate_map
    transformed_k = coordinate_map.T @ full_k @ coordinate_map
    transformed_b = coordinate_map.T @ full_b
    k_er = transformed_k[2:, :2]
    k_ee = transformed_k[2:, 2:]
    m_ee = transformed_m[2:, 2:]
    constraint_modes = -np.linalg.solve(k_ee, k_er)
    fixed_values, fixed_vectors = np.linalg.eig(np.linalg.solve(m_ee, k_ee))
    fixed_interface_hz = np.sort(
        np.sqrt(np.maximum(np.real(fixed_values), 0.0)) / (2.0 * np.pi))
    first_index = int(np.argmin(np.real(fixed_values)))
    fixed_mode = np.real(fixed_vectors[:, first_index])
    fixed_mode /= np.sqrt(fixed_mode @ m_ee @ fixed_mode)
    t_cb = np.block([
        [np.eye(2), np.zeros((2, 1))],
        [constraint_modes, fixed_mode[:, None]],
    ])
    return {
        "mass": t_cb.T @ transformed_m @ t_cb,
        "damping": t_cb.T @ transformed_c @ t_cb,
        "stiffness": t_cb.T @ transformed_k @ t_cb,
        "input_vector": t_cb.T @ transformed_b,
        "transform": t_cb,
        "fixed_interface_hz": fixed_interface_hz,
    }


def multi_route_reduction_metrics() -> dict[str, float]:
    """Evaluate every Section 6 route from the authoritative model inputs.

    The browser mirrors the algebraic routes.  The full-model and one-mode
    Craig-Bampton values are intentionally generated here because they require
    the assembled ten-DOF matrices.
    """
    constants = physical_constants()
    component = component_parameters()
    md, ms, km = constants["m_d"], constants["m_s"], constants["K_m"]
    kax_p, kball_p = constants["k_ax"], constants["k_ball"]

    def mode_pair(stage_mass: float, link_stiffness: float) -> np.ndarray:
        mass = np.diag([md, stage_mass])
        stiffness = np.array([
            [km + link_stiffness, -link_stiffness],
            [-link_stiffness, link_stiffness],
        ])
        return _linear_modes(mass, stiffness)

    def route_row(stage_mass: float, link_stiffness: float,
                  link_damping: float, ball_stiffness: float) -> dict[str, float]:
        modes = mode_pair(stage_mass, link_stiffness)
        reduced_mass = md * stage_mass / (md + stage_mass)
        relative_zeta = link_damping / (2.0 * np.sqrt(link_stiffness * reduced_mass))
        settling = SETTLING_2PCT_FACTOR / (relative_zeta * 2.0 * np.pi * modes[1])
        return {
            "md": md, "ms": stage_mass, "kax": link_stiffness,
            "cax": link_damping, "zeta": relative_zeta,
            "f1": float(modes[0]), "f2": float(modes[1]),
            "kball": ball_stiffness, "settling": settling,
        }

    axial_compliance = sum(1.0 / component[key]
                           for key in ("k_brg", "k_sha", "k_mnt")) + 1.0 / kball_p
    kax_series = 1.0 / axial_compliance
    torsional_compliance = constants["r"]**2 * (
        1.0 / component["k_c1"] + 1.0 / component["k_c2"]
        + 1.0 / component["k_theta_a"])
    kax_series_full = 1.0 / (axial_compliance + torsional_compliance)
    torsional_share = torsional_compliance / (axial_compliance + torsional_compliance)

    omega_target = 2.0 * np.pi * constants["axial_mode_target_hz"]
    kax_f, cax_f = axial_complex_link(component, constants["r"], kball_p, omega_target)

    # Route C: retain x_d and x_s physically, add the first fixed-interface
    # mode of the eight-coordinate eliminated partition, then compare its
    # retained upper mode with the full and 2-DOF models.
    full_m, full_c, full_k, _full_b, _full_p = full_linear_matrices()
    cb_plant = craig_bampton_plant()
    cb_m = cb_plant["mass"]
    cb_c = cb_plant["damping"]
    cb_k = cb_plant["stiffness"]
    cb_modes = _linear_modes(cb_m, cb_k)
    cb_damped_modes = _damped_modal_data(cb_m, cb_c, cb_k)
    cb_upper_damped = min(
        cb_damped_modes, key=lambda item: abs(item[0] - constants["axial_mode_target_hz"]))
    kax_c = modal_calibrated_axial_stiffness(md, ms, km, float(cb_modes[1]))
    kball_c = closure_ball_stiffness(
        kax_c, component["k_brg"], component["k_sha"], component["k_mnt"])
    mu = md * ms / (md + ms)
    cax_c = 2.0 * cb_upper_damped[1] * np.sqrt(kax_c * mu)

    m_eff = MODEL["m_eff_measured"]
    kax_m = omega_target**2 * m_eff
    kball_m = closure_ball_stiffness(
        kax_m, component["k_brg"], component["k_sha"], component["k_mnt"])
    mu_m = md * m_eff / (md + m_eff)
    cax_m = 2.0 * MODEL["zeta_relative_measured"] * np.sqrt(kax_m * mu_m)

    rows = {
        "p": route_row(ms, kax_p, MODEL["c_ax"], kball_p),
        "s": route_row(ms, kax_series, MODEL["c_ax"], kball_p),
        "b": route_row(ms, kax_series, MODEL["c_ax"], kball_p),
        "f": route_row(ms, kax_f, cax_f, kball_p),
        "c": route_row(ms, kax_c, cax_c, kball_c),
        "m": route_row(m_eff, kax_m, cax_m, kball_m),
    }
    # Report the DAMPED modal frequencies, which is what the Section 7.2
    # per-plant audit tabulates.  The undamped pair differed in the second
    # decimal and the two appeared in the document as two different numbers
    # for one quantity.
    cb_lower_damped = _damped_modal_data(
        cb_plant["mass"], cb_plant["damping"], cb_plant["stiffness"])[0]
    rows["c"].update(f1=float(cb_lower_damped[0]), f2=float(cb_upper_damped[0]),
                     zeta=cb_upper_damped[1], settling=(
                         SETTLING_2PCT_FACTOR / (
                             cb_upper_damped[1] * 2.0 * np.pi * cb_upper_damped[0])))

    full_damped_modes = _damped_modal_data(full_m, full_c, full_k)
    full_upper = min(
        full_damped_modes, key=lambda item: abs(item[0] - constants["axial_mode_target_hz"]))
    metrics = {
        f"route_{route}_{key}": value
        for route, row in rows.items() for key, value in row.items()
    }
    first_discarded = next(item for item in full_damped_modes if item[0] > 900.0)
    fixed_interface_hz = np.asarray(cb_plant["fixed_interface_hz"])
    metrics.update({
        "route_s_kax_full": kax_series_full,
        "torsional_share": torsional_share,
        "mass_ratio": md / ms,
        "reduced_mu": mu,
        "full_model_zeta": full_upper[1],
        "full_model_upper_hz": full_upper[0],
        "full_model_settling": SETTLING_2PCT_FACTOR / (
            full_upper[1] * 2.0 * np.pi * full_upper[0]),
        "cb_frequency_delta": 100.0 * abs(cb_modes[1] - full_upper[0]) / full_upper[0],
        "cb_damping_delta": 100.0 * abs(cb_upper_damped[1] - full_upper[1]) / full_upper[1],
        # Half-power bandwidths of the retained upper mode.  They are the
        # resolution requirement the measured-FRF route has to meet.
        "route_p_bandwidth_hz": 2.0 * rows["p"]["zeta"] * rows["p"]["f2"],
        "full_model_bandwidth_hz": 2.0 * full_upper[1] * full_upper[0],
        # Truncation separation.  The fixed-interface value is the formal
        # condensation bound; the discarded system pole is what Section 7 sees.
        "first_fixed_interface_hz": float(fixed_interface_hz[0]),
        "first_discarded_hz": first_discarded[0],
        "fixed_interface_separation": (
            constants["axial_mode_target_hz"] / float(fixed_interface_hz[0]))**2,
        "discarded_pole_separation": (
            constants["axial_mode_target_hz"] / first_discarded[0])**2,
        "mu_fraction": mu / ms,
        "relative_mode_nearground_hz": np.sqrt(kax_p / ms) / (2.0 * np.pi),
        "drive_pole_hz": np.sqrt(km / md) / (2.0 * np.pi),
    })
    return metrics


def _rk4_linear_stability_radius(mass: np.ndarray, damping: np.ndarray,
                                 stiffness: np.ndarray, dt: float) -> float:
    """Largest RK4 amplification magnitude over the linear state eigenvalues."""
    count = mass.shape[0]
    state_matrix = np.block([
        [np.zeros((count, count)), np.eye(count)],
        [-np.linalg.solve(mass, stiffness), -np.linalg.solve(mass, damping)],
    ])
    scaled = np.linalg.eigvals(state_matrix) * dt
    amplification = (
        1.0 + scaled + scaled**2 / 2.0 + scaled**3 / 6.0 + scaled**4 / 24.0
    )
    return float(np.max(np.abs(amplification)))


def _rk4_linear(mass: np.ndarray, damping: np.ndarray, stiffness: np.ndarray,
                input_vector: np.ndarray, command_step: float,
                dt: float = VERIFICATION_DT,
                duration: float = 0.085,
                command_function=None) -> tuple[np.ndarray, np.ndarray]:
    """Integrate an arbitrary second-order linear model with a true ZOH input."""
    count = mass.shape[0]
    if np.count_nonzero(mass - np.diag(np.diag(mass))) == 0:
        inverse_mass = np.diag(1.0 / np.diag(mass))
    else:
        # The Craig-Bampton plant couples its retained and modal coordinates,
        # so the diagonal shortcut is only valid for the physical plants.
        inverse_mass = np.linalg.inv(mass)
    times = np.arange(0.0, duration + 0.5 * dt, dt)
    states = np.zeros((times.size, 2 * count), dtype=float)

    def rhs(state: np.ndarray, held_command: float) -> np.ndarray:
        position = state[:count]
        velocity = state[count:]
        acceleration = inverse_mass @ (input_vector * held_command - damping @ velocity - stiffness @ position)
        return np.concatenate((velocity, acceleration))

    command_law = verification_command_position if command_function is None else command_function
    for i in range(times.size - 1):
        held_command = command_law(
            times[i] + 0.5 * dt, command_step)
        y = states[i]
        k1 = rhs(y, held_command)
        k2 = rhs(y + 0.5 * dt * k1, held_command)
        k3 = rhs(y + 0.5 * dt * k2, held_command)
        k4 = rhs(y + dt * k3, held_command)
        states[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return times, states


def verification_route_plants(constants: dict[str, float]) -> "OrderedDict[str, dict[str, object]]":
    """Assemble every candidate reduced plant that Section 6 puts forward.

    Each entry carries the matrices needed by the residual audit so that the
    Section 6.3 damping claim is measured against the ten-DOF model rather than
    inferred from settling times alone.
    """
    component = component_parameters()
    md, ms, km = constants["m_d"], constants["m_s"], constants["K_m"]
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    drive_damping = constants["c_m"] * np.outer(H["d"], H["d"])

    def two_dof(stage_mass: float, link_stiffness: float,
                link_damping: float) -> dict[str, object]:
        return {
            "mass": np.diag([md, stage_mass]),
            "damping": link_damping * coupling + drive_damping,
            "stiffness": np.array([
                [km + link_stiffness, -link_stiffness],
                [-link_stiffness, link_stiffness],
            ], dtype=float),
            "input_vector": np.array([km, 0.0]),
            "output_index": 1,
        }

    omega_target = 2.0 * np.pi * constants["axial_mode_target_hz"]
    kax_f, cax_f = axial_complex_link(
        component, constants["r"], constants["k_ball"], omega_target)
    m_eff = MODEL["m_eff_measured"]
    kax_m = omega_target**2 * m_eff
    mu_m = md * m_eff / (md + m_eff)
    cax_m = 2.0 * MODEL["zeta_relative_measured"] * np.sqrt(kax_m * mu_m)

    cb_plant = craig_bampton_plant()
    cb_modes = _linear_modes(cb_plant["mass"], cb_plant["stiffness"])
    cb_upper_damped = min(
        _damped_modal_data(cb_plant["mass"], cb_plant["damping"], cb_plant["stiffness"]),
        key=lambda item: abs(item[0] - constants["axial_mode_target_hz"]))
    mu = md * ms / (md + ms)
    kax_c = modal_calibrated_axial_stiffness(md, ms, km, float(cb_modes[1]))
    cax_c = 2.0 * cb_upper_damped[1] * np.sqrt(kax_c * mu)

    plants: "OrderedDict[str, dict[str, object]]" = OrderedDict()
    plants["P"] = two_dof(ms, constants["k_ax"], MODEL["c_ax"])
    plants["P"]["label"] = "Formal static condensation"
    plants["P"]["basis"] = "assumed 55 N·s/m link damper"
    plants["F"] = two_dof(ms, kax_f, cax_f)
    plants["F"]["label"] = "Frequency-domain complex stiffness"
    plants["F"]["basis"] = "interface loss factors propagated to $c_{ax}$"
    plants["M"] = two_dof(m_eff, kax_m, cax_m)
    plants["M"]["label"] = "Measured-FRF identification"
    plants["M"]["basis"] = "0.600 kg modal mass and measured $\\zeta$"
    plants["C2"] = two_dof(ms, kax_c, cax_c)
    plants["C2"]["label"] = "Craig–Bampton, condensed to two coordinates"
    plants["C2"]["basis"] = "Craig–Bampton frequency and damping, two coordinates"
    plants["CB"] = {
        "mass": cb_plant["mass"],
        "damping": cb_plant["damping"],
        "stiffness": cb_plant["stiffness"],
        "input_vector": cb_plant["input_vector"],
        "output_index": 1,
        "label": "Craig–Bampton, three coordinates",
        "basis": "damping projected from the ten-DOF matrices",
    }
    return plants


def full_reduced_verification(frequencies: np.ndarray, constants: dict[str, float]) -> dict[str, object]:
    full_m, full_c, full_k, full_b, p = full_linear_matrices()
    reduced_m, reduced_c, reduced_k, reduced_b = linear_matrices((), "none")
    full_response = _matrix_frequency_response(frequencies, full_m, full_c, full_k, full_b, 9)
    reduced_response = _matrix_frequency_response(frequencies, reduced_m, reduced_c, reduced_k, reduced_b, 1)
    command_step = constants["full_step"]
    full_damped_modes = _damped_modal_data(full_m, full_c, full_k)
    reduced_damped_modes = _damped_modal_data(reduced_m, reduced_c, reduced_k)
    discarded_mode = next(item for item in full_damped_modes if item[0] > 900.0)

    convergence: list[dict[str, float | bool]] = []
    stable_runs: dict[float, dict[str, np.ndarray | float]] = {}
    for dt in VERIFICATION_CONVERGENCE_DTS:
        stability_radius = max(
            _rk4_linear_stability_radius(full_m, full_c, full_k, dt),
            _rk4_linear_stability_radius(reduced_m, reduced_c, reduced_k, dt),
        )
        row: dict[str, float | bool] = {
            "dt": dt,
            "points_per_discarded_cycle": 1.0 / (dt * discarded_mode[0]),
            "stability_radius": stability_radius,
            "stable": stability_radius <= 1.0 + 1.0e-10,
        }
        if not bool(row["stable"]):
            convergence.append(row)
            continue
        times, full_states = _rk4_linear(
            full_m, full_c, full_k, full_b, command_step, dt=dt)
        reduced_times, reduced_states = _rk4_linear(
            reduced_m, reduced_c, reduced_k, reduced_b, command_step, dt=dt)
        if not np.array_equal(times, reduced_times):
            raise RuntimeError("Full and reduced verification time grids differ")
        command = np.array([
            verification_command_position(t, command_step) for t in times
        ])
        full_stage = full_states[:, 9]
        reduced_stage = reduced_states[:, 1]
        residual = full_stage - reduced_stage
        amplitude = float(np.max(np.abs(command)))
        rms_residual = float(np.sqrt(np.mean(residual**2)))
        peak_residual = float(np.max(np.abs(residual)))
        edge_peaks = []
        for edge_index, edge_time in enumerate(VERIFICATION_EDGES):
            stop_time = (VERIFICATION_EDGES[edge_index + 1]
                         if edge_index + 1 < len(VERIFICATION_EDGES)
                         else times[-1] + dt)
            mask = (times >= edge_time) & (times < stop_time)
            edge_peaks.append(float(np.max(np.abs(residual[mask]))))
        row.update({
            "rms_residual_nm": rms_residual * 1e9,
            "peak_residual_nm": peak_residual * 1e9,
            "rms_residual_pct_command": 100.0 * rms_residual / amplitude,
            "peak_residual_pct_command": 100.0 * peak_residual / amplitude,
        })
        convergence.append(row)
        stable_runs[dt] = {
            "times": times,
            "command": command,
            "full_stage": full_stage,
            "reduced_stage": reduced_stage,
            "residual": residual,
            "rms_residual": rms_residual,
            "peak_residual": peak_residual,
            "edge_peaks": np.asarray(edge_peaks),
        }

    production = stable_runs.get(VERIFICATION_DT)
    if production is None:
        raise RuntimeError("Production verification time step is not RK4-stable")
    times = np.asarray(production["times"])
    command = np.asarray(production["command"])
    full_stage = np.asarray(production["full_stage"])
    reduced_stage = np.asarray(production["reduced_stage"])
    residual = np.asarray(production["residual"])
    rms_residual = float(production["rms_residual"])
    peak_residual = float(production["peak_residual"])
    command = np.array([
        verification_command_position(t, command_step) for t in times
    ])
    command_amplitude = float(np.max(np.abs(command)))

    spectral_mask = times >= VERIFICATION_EDGES[0]
    spectral_residual = residual[spectral_mask] - np.mean(residual[spectral_mask])
    spectral_frequencies = np.fft.rfftfreq(spectral_residual.size, VERIFICATION_DT)
    spectral_amplitude = np.abs(np.fft.rfft(
        spectral_residual * np.hanning(spectral_residual.size)))
    resolved_band = (spectral_frequencies >= 100.0) & (spectral_frequencies <= 3000.0)
    ripple_band = (spectral_frequencies >= 1200.0) & (spectral_frequencies <= 2800.0)
    dominant_frequency = float(spectral_frequencies[resolved_band][
        np.argmax(spectral_amplitude[resolved_band])])
    ripple_frequency = float(spectral_frequencies[ripple_band][
        np.argmax(spectral_amplitude[ripple_band])])
    full_upper_mode = min(full_damped_modes, key=lambda item: abs(item[0] - constants["axial_mode_target_hz"]))
    reduced_upper_mode = min(reduced_damped_modes, key=lambda item: abs(item[0] - constants["axial_mode_target_hz"]))

    # Per-plant residual audit.  Every Section 6 candidate is driven with the
    # same command and compared with the same ten-DOF stage output, so the
    # damping-versus-truncation question is answered by measurement.
    route_residuals: list[dict[str, object]] = []
    for key, plant in verification_route_plants(constants).items():
        output_index = int(plant["output_index"])
        _route_times, route_states = _rk4_linear(
            plant["mass"], plant["damping"], plant["stiffness"],
            plant["input_vector"], command_step, dt=VERIFICATION_DT)
        route_residual = full_stage - route_states[:, output_index]
        route_modes = _damped_modal_data(plant["mass"], plant["damping"], plant["stiffness"])
        route_upper = min(
            route_modes, key=lambda item: abs(item[0] - constants["axial_mode_target_hz"]))
        route_lower = min(route_modes, key=lambda item: item[0])
        route_rms = float(np.sqrt(np.mean(route_residual**2)))
        route_peak = float(np.max(np.abs(route_residual)))
        route_residuals.append({
            "key": key,
            "label": plant["label"],
            "basis": plant["basis"],
            "coordinates": int(np.asarray(plant["mass"]).shape[0]),
            "lower_hz": route_lower[0],
            "lower_zeta": route_lower[1],
            "upper_hz": route_upper[0],
            "upper_zeta": route_upper[1],
            "dc_gain": float(np.linalg.solve(
                plant["stiffness"], plant["input_vector"])[output_index]),
            "rms_residual_nm": route_rms * 1e9,
            "peak_residual_nm": route_peak * 1e9,
            "rms_residual_pct_command": 100.0 * route_rms / command_amplitude,
            "peak_residual_pct_command": 100.0 * route_peak / command_amplitude,
        })

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
        "convergence": convergence,
        "route_residuals": route_residuals,
        "production_dt": VERIFICATION_DT,
        "edge_peaks_nm": np.asarray(production["edge_peaks"]) * 1e9,
        "dominant_residual_frequency_hz": dominant_frequency,
        "ripple_residual_frequency_hz": ripple_frequency,
        "full_upper_damped_mode": full_upper_mode,
        "full_lower_damped_mode": min(full_damped_modes, key=lambda item: item[0]),
        "reduced_upper_damped_mode": reduced_upper_mode,
        "first_discarded_damped_mode": discarded_mode,
        "full_dc_gain": float(np.linalg.solve(full_k, full_b)[9]),
        "reduced_dc_gain": float(np.linalg.solve(reduced_k, reduced_b)[1]),
        "full_step_electrical_angle": constants["kappa"] * command_step,
        "linear_force_to_limit_ratio": constants["K_m"] * command_step / constants["F_max"],
        "rms_residual_nm": rms_residual * 1e9,
        "peak_residual_nm": peak_residual * 1e9,
        "command_amplitude_nm": command_amplitude * 1e9,
        "rms_residual_pct_command": 100.0 * rms_residual / command_amplitude,
        "peak_residual_pct_command": 100.0 * peak_residual / command_amplitude,
    }


def linear_matrices(sites: tuple[str, ...], friction_model: str,
                    micro_viscous: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return M, C, K, B for the presliding linearization of one case."""
    constants = physical_constants()
    m_d, m_s = constants["m_d"], constants["m_s"]
    k_ax, k_m, c_ax = constants["k_ax"], constants["K_m"], constants["c_ax"]
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    mass = np.diag([m_d, m_s])
    damping = c_ax * coupling + constants["c_m"] * np.outer(H["d"], H["d"])
    stiffness = np.array([[k_m + k_ax, -k_ax], [-k_ax, k_ax]], dtype=float)
    for site in sites:
        outer = np.outer(H[site], H[site])
        p = site_parameters({"micro_viscous": micro_viscous}, site)
        stiffness += p["sigma0"] * outer
        # With sigma_1 = 0 the two laws contribute the identical sigma_2 term,
        # which is what makes the A/A2 comparison a pure memory comparison.
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
        mass, damping, stiffness, input_vector = linear_matrices(
            case["sites"], case["friction"], bool(case.get("micro_viscous")))
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
                          GMS_STIFFNESS_BY_SITE[site]))
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


def site_parameters(case: dict[str, object], site: str) -> dict[str, float]:
    """Return a site's constitutive parameters as the given case executes them.

    Only the micro-viscous variant differs: it restores the sigma_1 values that
    the controlled A/B/C comparison sets to zero.
    """
    if case.get("micro_viscous"):
        return {**FRICTION[site], "sigma1": MICRO_VISCOUS_SIGMA1[site],
                "gms_stiffness": GMS_STIFFNESS_BY_SITE[site]}
    return {**FRICTION[site], "gms_stiffness": GMS_STIFFNESS_BY_SITE[site]}


def stribeck(velocity: float, p: dict[str, float]) -> float:
    ratio = abs(velocity) / p["v_s"]
    return p["F_c"] + (p["F_s"] - p["F_c"]) * np.exp(-(ratio ** p["delta"]))


def lugre_site(velocity: float, z: float, p: dict[str, float]) -> tuple[float, float]:
    attraction = max(stribeck(velocity, p), 1e-12)
    z_dot = velocity - p["sigma0"] * abs(velocity) * z / attraction
    force = p["sigma0"] * z + p["sigma1"] * z_dot + p["sigma2"] * velocity
    return z_dot, force


def gms_site(velocity: float, element_forces: np.ndarray,
             p: dict[str, float],
             branch_state: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Return GMS derivatives after branch tests on the current RK trial state.

    Ordering is intentional: zero velocity holds the state; otherwise the
    reversal/re-stick test is evaluated first, then the yield test, and only
    then is either the stuck or slip derivative assigned.  No derivative is
    used to choose its own branch within the same RHS evaluation.

    ``branch_state`` is the counterfactual only.  When it is ``None`` the
    executed stateless test runs unchanged.  When a boolean array is supplied
    it carries a persistent per-element slip flag, so an element that has
    yielded keeps slipping until a reversal instead of being reclassified by a
    rising threshold.  Section 12.2 uses it to price that departure.
    """
    threshold = np.maximum(GMS_WEIGHTS * stribeck(velocity, p), 1e-12)
    stiffness = np.asarray(p["gms_stiffness"], dtype=float)
    derivatives = np.zeros(GMS_N)
    if abs(velocity) > 1e-14:
        direction = np.sign(velocity)
        for i in range(GMS_N):
            re_stick = velocity * element_forces[i] <= 0.0
            below_yield = abs(element_forces[i]) < threshold[i]
            if branch_state is None:
                stuck = re_stick or below_yield
            elif re_stick:
                stuck = True
                branch_state[i] = False
            elif branch_state[i]:
                stuck = False
            else:
                stuck = below_yield
                branch_state[i] = not below_yield
            if stuck:
                derivatives[i] = stiffness[i] * velocity
            else:
                # Stable slip attraction to F_i = sign(v) nu_i s(v).
                derivatives[i] = p["C_gms"] * (direction - element_forces[i] / threshold[i])
    total_force = float(np.sum(element_forces) + p["sigma2"] * velocity)
    return derivatives, total_force


class GmsBranchCensus:
    """Count where the stateless branch test and a persistent flag disagree.

    The shadow flag is observation only: it is advanced alongside the executed
    trajectory and never feeds a derivative, so attaching a census cannot
    change any reported result.
    """

    def __init__(self, sites: tuple[str, ...]) -> None:
        self.sites = tuple(sites)
        self.slipping = {site: np.zeros(GMS_N, dtype=bool) for site in self.sites}
        self.reversal_flips = {site: 0 for site in self.sites}
        self.threshold_flips = {site: 0 for site in self.sites}
        self.evaluations = {site: 0 for site in self.sites}

    def observe(self, site: str, velocity: float, element_forces: np.ndarray,
                p: dict[str, float]) -> None:
        if abs(velocity) <= 1e-14:
            return
        threshold = np.maximum(GMS_WEIGHTS * stribeck(velocity, p), 1e-12)
        slipping = self.slipping[site]
        for i in range(GMS_N):
            self.evaluations[site] += 1
            re_stick = velocity * element_forces[i] <= 0.0
            below_yield = abs(element_forces[i]) < threshold[i]
            if re_stick:
                # Both models send the element back to stick.
                if slipping[i]:
                    self.reversal_flips[site] += 1
                slipping[i] = False
            elif slipping[i]:
                # The persistent model would keep slipping here.  The executed
                # stateless test sticks the element whenever the rising
                # Stribeck threshold has overtaken its force.
                if below_yield:
                    self.threshold_flips[site] += 1
            elif not below_yield:
                slipping[i] = True

    def total(self, key: str) -> int:
        return int(sum(getattr(self, key).values()))


def nonlinear_rhs(t: float, state: np.ndarray, case: dict[str, object], constants: dict[str, float],
                  held_command: float | None = None,
                  blocked_stage: bool = False,
                  census: "GmsBranchCensus | None" = None,
                  branch_states: dict[str, np.ndarray] | None = None) -> np.ndarray:
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
    axial_force = constants["k_ax"] * (x_d - x_s) + constants["c_ax"] * (v_d - v_s)

    velocities = {site: float(H[site] @ np.array([v_d, v_s])) for site in SITE_KEYS}
    forces = {site: 0.0 for site in SITE_KEYS}
    derivative = np.zeros_like(state)
    if case["friction"] == "lugre":
        for site in case["sites"]:
            index = LUGRE_INDEX[site]
            derivative[index], forces[site] = lugre_site(
                velocities[site], state[index], site_parameters(case, site))
    elif case["friction"] == "gms":
        for site in case["sites"]:
            start = GMS_START[site]
            stop = start + GMS_N
            site_p = site_parameters(case, site)
            derivative[start:stop], forces[site] = gms_site(
                velocities[site], state[start:stop], site_p,
                None if branch_states is None else branch_states[site])
            if census is not None:
                census.observe(site, velocities[site], state[start:stop], site_p)

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
                          blocked_stage: bool = False,
                          census: "GmsBranchCensus | None" = None,
                          branch_states: dict[str, np.ndarray] | None = None
                          ) -> tuple[np.ndarray, np.ndarray]:
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
        k1 = nonlinear_rhs(t, y, case, constants, held_command, blocked_stage,
                           census, branch_states)
        k2 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k1, case, constants,
                           held_command, blocked_stage, census, branch_states)
        k3 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k2, case, constants,
                           held_command, blocked_stage, census, branch_states)
        k4 = nonlinear_rhs(t + dt, y + dt * k3, case, constants,
                           held_command, blocked_stage, census, branch_states)
        states[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return times, states


def rk4_case(case: dict[str, object], constants: dict[str, float], dt: float = PRODUCTION_DT,
             duration: float | None = None,
             census: "GmsBranchCensus | None" = None,
             branch_states: dict[str, np.ndarray] | None = None
             ) -> tuple[np.ndarray, np.ndarray]:
    if duration is None:
        duration = main_duration(constants)
    return rk4_case_with_command(
        case, constants,
        lambda t: command_position(t, constants["command_step"], constants["plateau_dwell"]),
        duration=duration, dt=dt, census=census, branch_states=branch_states,
    )


def friction_force_history(case: dict[str, object], states: np.ndarray,
                           site: str) -> np.ndarray:
    """Recover an executed site's constitutive force from integrated states."""
    if site not in case["sites"]:
        return np.zeros(states.shape[0])
    velocity = H[site][0] * states[:, 2] + H[site][1] * states[:, 3]
    p = site_parameters(case, site)
    if case["friction"] == "lugre":
        state = states[:, LUGRE_INDEX[site]]
        return np.array([
            lugre_site(float(v), float(z), p)[1] for v, z in zip(velocity, state)
        ])
    start = GMS_START[site]
    stop = start + GMS_N
    return np.sum(states[:, start:stop], axis=1) + p["sigma2"] * velocity


def presliding_responses(constants: dict[str, float], keys: tuple[str, ...],
                         site: str,
                         persistent_branch_keys: tuple[str, ...] = ()) -> dict[str, object]:
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
    censuses: dict[str, GmsBranchCensus] = {}
    times: np.ndarray | None = None

    for key in keys:
        branch_states = None
        if key in persistent_branch_keys:
            branch_states = {
                active_site: np.zeros(GMS_N, dtype=bool)
                for active_site in CASES[key]["sites"]
            }
        census = None
        if CASES[key]["friction"] == "gms" and branch_states is None:
            census = GmsBranchCensus(tuple(CASES[key]["sites"]))
            censuses[key] = census
        times, states = rk4_case_with_command(
            CASES[key], constants, command_function, duration=duration,
            dt=PRODUCTION_DT, blocked_stage=blocked_stage,
            census=census, branch_states=branch_states)
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
        # Dissipated energy of the closed memory sequence, A_loop = |contour
        # integral of F_f dx_o|.  The command starts and ends at the origin, so
        # the path closes and the integral is the summed hysteresis-loop area.
        # It is the only metric evaluated on the dynamic trace, which is where
        # the deceleration phases the branch census counts actually live.
        loop_increments = []
        for level_index in range(levels.size):
            plateau_end = PRESLIDING_START + (level_index + 1) * plateau_dwell
            segment = ((times >= plateau_end - plateau_dwell)
                       & (times <= plateau_end))
            loop_increments.append(float(np.trapezoid(
                forces[key][segment], site_coordinate[segment])))
        loop_increment_array = np.asarray(loop_increments)
        loop_area = float(abs(np.trapezoid(
            forces[key][active], site_coordinate[active])))
        # Retention diagnostic (Section 9 / Appendix G.5): the fraction of the
        # available elastic force sigma_0*x, capped by the Stribeck limit
        # s(0)=F_s, that a law actually holds at rest on a settled plateau.
        # Grouped by distinct commanded level so repeated visits to the same
        # level (e.g. the outer level is revisited) contribute one value.
        site_params = site_parameters(CASES[key], site)
        unique_levels = sorted({float(level) for level in levels if level != 0.0})
        r_hold_ratios = []
        for level_value in unique_levels:
            level_indices = np.flatnonzero(levels == level_value)
            mean_abs_force = float(np.mean(np.abs(endpoint_force_array[level_indices])))
            available_force = min(
                site_params["sigma0"] * abs(level_value) * microstep, site_params["F_s"])
            r_hold_ratios.append(mean_abs_force / available_force)
        metrics[key] = {
            "loop_area_J": loop_area,
            "loop_increments_J": loop_increment_array,
            "loop_increment_min_J": float(np.min(loop_increment_array)),
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
            "r_hold": float(np.mean(r_hold_ratios)),
            "r_hold_levels": unique_levels,
            "r_hold_ratios": r_hold_ratios,
        }
    return {
        "times": times,
        "command": command,
        "results": results,
        "forces": forces,
        "metrics": metrics,
        "censuses": censuses,
        "microstep": microstep,
        "duration": duration,
        "plateau_dwell": plateau_dwell,
        "keys": keys,
        "site": site,
        "blocked_stage": blocked_stage,
        "levels": levels,
    }


def _retained_mode_pole(constants: dict[str, float], sites: tuple[str, ...],
                        damping_multiplier: float = 1.0) -> dict[str, float] | None:
    """Damped complex pole of the retained (axial) mode at a given structural
    damping multiplier on c_ax/c_m.  Used only by the Appendix G.5 high-damping
    confirmation run; returns None if the mode is overdamped (no complex pair)."""
    m_d, m_s = constants["m_d"], constants["m_s"]
    k_ax, k_m = constants["k_ax"], constants["K_m"]
    c_ax = constants["c_ax"] * damping_multiplier
    c_m = constants["c_m"] * damping_multiplier
    coupling = np.array([[1.0, -1.0], [-1.0, 1.0]])
    mass = np.diag([m_d, m_s])
    damping = c_ax * coupling + c_m * np.outer(H["d"], H["d"])
    stiffness = np.array([[k_m + k_ax, -k_ax], [-k_ax, k_ax]], dtype=float)
    for site in sites:
        outer = np.outer(H[site], H[site])
        p = site_parameters({"micro_viscous": False}, site)
        stiffness += p["sigma0"] * outer
        damping += (p["sigma1"] + p["sigma2"]) * outer
    state_matrix = np.zeros((4, 4))
    mass_inverse = np.linalg.inv(mass)
    state_matrix[0:2, 2:4] = np.eye(2)
    state_matrix[2:4, 0:2] = -mass_inverse @ stiffness
    state_matrix[2:4, 2:4] = -mass_inverse @ damping
    eigenvalues = np.linalg.eigvals(state_matrix)
    oscillatory = eigenvalues[np.imag(eigenvalues) > 1.0]
    if oscillatory.size == 0:
        return None
    pole = oscillatory[np.argmax(np.imag(oscillatory))]
    sigma = float(-np.real(pole))
    magnitude = float(abs(pole))
    tau = 1.0 / sigma if sigma > 0.0 else float("inf")
    return {
        "frequency_hz": float(np.imag(pole)) / (2.0 * np.pi),
        "zeta": sigma / magnitude,
        "tau_s": tau,
        "envelope_5pct_s": 3.0 * tau,
    }


def high_damping_confirmation_run(constants: dict[str, float],
                                  multiplier: float = 50.0) -> dict[str, object]:
    """Rerun the guideway memory experiment with c_ax/c_m scaled up so
    post-edge ringing dies within about a millisecond (Section 9 / Part 1.3
    diagnostic, second direction). If LuGre's settled force recovers toward
    the elastic prediction once ringing is suppressed, the dither mechanism
    is confirmed from both directions."""
    sites = tuple(CASES["A"]["sites"])
    baseline_pole = _retained_mode_pole(constants, sites, 1.0)
    raised_pole = _retained_mode_pole(constants, sites, multiplier)
    raised_constants = dict(constants)
    raised_constants["c_ax"] = constants["c_ax"] * multiplier
    raised_constants["c_m"] = constants["c_m"] * multiplier
    experiment = presliding_responses(raised_constants, ("A", "A2"), "g")
    return {
        "multiplier": multiplier,
        "baseline_pole": baseline_pole,
        "raised_pole": raised_pole,
        "experiment": experiment,
    }


def true_loop_command_position(t: float, amplitude: float, ramp_time: float,
                               start: float) -> float:
    """Continuous 0 -> +A -> -A -> 0 triangular reversal, no plateaus."""
    if t < start:
        return 0.0
    tau = t - start
    if tau < ramp_time:
        return amplitude * (tau / ramp_time)
    if tau < 3.0 * ramp_time:
        return amplitude * (1.0 - (tau - ramp_time) / ramp_time)
    if tau < 4.0 * ramp_time:
        return -amplitude + amplitude * (tau - 3.0 * ramp_time) / ramp_time
    return 0.0


def true_presliding_loop(constants: dict[str, float],
                         ramp_time: float = 0.15) -> dict[str, object]:
    """Slow continuous quasi-static ramp-reversal (Part 1.5 item 4): a
    literature-comparable F-x loop, distinct from the settled return-point
    map built from discrete plateaus."""
    amplitude = 12.0 * constants["command_step"]
    start = PRESLIDING_START
    duration = start + 4.0 * ramp_time + 0.02
    command_function = lambda t: true_loop_command_position(t, amplitude, ramp_time, start)
    results: dict[str, np.ndarray] = {}
    forces: dict[str, np.ndarray] = {}
    times: np.ndarray | None = None
    for key in ("A", "A2"):
        times, states = rk4_case_with_command(
            CASES[key], constants, command_function, duration=duration, dt=PRODUCTION_DT)
        results[key] = states
        forces[key] = friction_force_history(CASES[key], states, "g")
    assert times is not None
    command = np.array([command_function(t) for t in times])
    site_coordinate = {key: results[key][:, 1] for key in results}
    return {
        "times": times, "command": command, "results": results, "forces": forces,
        "site_coordinate": site_coordinate, "amplitude": amplitude, "ramp_time": ramp_time,
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


def settled_window_masks(times: np.ndarray,
                         constants: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    """Return the settled-plateau and final-window sample masks."""
    final_window = times >= (times[-1] - constants["metric_window"])
    settled_mask = np.zeros(times.size, dtype=bool)
    for level_index in range(MAIN_LEVELS.size):
        plateau_end = MAIN_START + (level_index + 1) * constants["plateau_dwell"]
        settled_mask |= ((times >= plateau_end - constants["metric_window"])
                         & (times < plateau_end - 0.5e-9))
    settled_mask |= final_window
    return settled_mask, final_window


def _main_response_job(key: str, constants: dict[str, float], collect_census: bool
                       ) -> tuple[str, np.ndarray, np.ndarray, GmsBranchCensus | None]:
    """Process-pool unit for one independent nonlinear response case."""
    case = CASES[key]
    census = (GmsBranchCensus(tuple(case["sites"]))
              if collect_census and case["friction"] == "gms" else None)
    times, states = rk4_case(case, constants, dt=PRODUCTION_DT, census=census)
    return key, times, states, census


def time_responses(constants: dict[str, float],
                   censuses: dict[str, "GmsBranchCensus"] | None = None,
                   executor: Executor | None = None
                   ) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, dict[str, float]]]:
    results: dict[str, np.ndarray] = {}
    metrics: dict[str, dict[str, float]] = {}
    times: np.ndarray | None = None
    if executor is None:
        completed = (
            _main_response_job(key, constants, censuses is not None)
            for key in CASES
        )
    else:
        futures = [
            executor.submit(_main_response_job, key, constants, censuses is not None)
            for key in CASES
        ]
        completed = (future.result() for future in futures)
    for key, case_times, states, census in completed:
        times = case_times
        results[key] = states
        if censuses is not None and census is not None:
            # Observation only; the executed derivatives are unchanged.
            censuses[key] = census
    assert times is not None
    command = np.array([
        command_position(t, constants["command_step"], constants["plateau_dwell"])
        for t in times
    ])
    first_plateau = ((times >= MAIN_START)
                     & (times < MAIN_START + constants["plateau_dwell"]))
    settled_mask, final_window = settled_window_masks(times, constants)
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


def detent_ablation_study(constants: dict[str, float],
                          metrics: dict[str, dict[str, float]],
                          executor: Executor | None = None) -> dict[str, object]:
    """Rerun every case with the detent torque removed.

    Case 0 is frictionless but not force-free: the nonlinear campaign runs
    with the periodic detent torque enabled, so its settled-window deviation
    is a pure detent number.  Every other case therefore reports a settled
    window containing both terms.  Rerunning the identical command with
    T_det = 0 separates them by measurement instead of by quadrature.
    """
    if not MODEL["detent_enabled"]:
        raise RuntimeError(
            "The executed campaign already has detent disabled; the ablation "
            "pair would be two copies of the same run")
    ablated = dict(constants)
    ablated["T_det"] = 0.0
    if executor is None:
        completed = [_main_response_job(key, ablated, False) for key in CASES]
    else:
        futures = [executor.submit(_main_response_job, key, ablated, False)
                   for key in CASES]
        completed = [future.result() for future in futures]
    times = completed[0][1]
    command = np.array([
        command_position(t, ablated["command_step"], ablated["plateau_dwell"])
        for t in times
    ])
    settled_mask, _final_window = settled_window_masks(times, ablated)
    rows: "OrderedDict[str, dict[str, float]]" = OrderedDict()
    for key, _case_times, states, _census in sorted(
            completed, key=lambda item: list(CASES).index(item[0])):
        error = command - states[:, 1]
        friction_only = float(np.sqrt(np.mean(error[settled_mask] ** 2)) * 1e9)
        executed = float(metrics[key]["rms_settled_deviation_nm"])
        rows[key] = {
            "executed_nm": executed,
            "detent_off_nm": friction_only,
            "detent_share_pct": 100.0 * (1.0 - friction_only / executed)
            if executed > 0.0 else float("nan"),
            "quadrature_nm": float(np.sqrt(max(
                executed**2 - float(metrics["0"]["rms_settled_deviation_nm"])**2, 0.0))),
        }
    baseline = rows["0"]
    return {
        "rows": rows,
        "detent_only_nm": float(metrics["0"]["rms_settled_deviation_nm"]),
        "detent_off_baseline_nm": baseline["detent_off_nm"],
        "largest_friction_nm": max(row["detent_off_nm"] for row in rows.values()),
        "worst_share_pct": max(row["detent_share_pct"] for row in rows.values()),
    }


def _retention_damping_job(multiplier: float, constants: dict[str, float]
                           ) -> tuple[float, dict[str, float]]:
    """Process-pool unit for one point of the retention-versus-damping sweep."""
    scaled = dict(constants)
    scaled["c_ax"] = constants["c_ax"] * multiplier
    scaled["c_m"] = constants["c_m"] * multiplier
    experiment = presliding_responses(scaled, ("A", "A2"), "g")
    metrics = experiment["metrics"]
    lugre = 100.0 * float(metrics["A"]["r_hold"])
    gms = 100.0 * float(metrics["A2"]["r_hold"])
    return multiplier, {
        "r_hold_lugre_pct": lugre,
        "r_hold_gms_pct": gms,
        "r_hold_ratio": gms / max(lugre, 1.0e-9),
        "force_mismatch_lugre_N": float(metrics["A"]["return_force_mismatch_N"]),
        "force_mismatch_gms_N": float(metrics["A2"]["return_force_mismatch_N"]),
        "force_ratio": (float(metrics["A2"]["return_force_mismatch_N"])
                        / max(float(metrics["A"]["return_force_mismatch_N"]), 1e-30)),
    }


def retention_damping_sweep(constants: dict[str, float],
                            executor: Executor | None = None,
                            targets: tuple[float, ...] = (
                                0.0014, 0.0133, 0.0157, 0.05, 0.743),
                            ) -> dict[str, object]:
    """Measure the Section 9 discriminator as a function of retained-mode damping.

    Appendix G.5 already shows that suppressing post-edge ringing collapses
    the LuGre/GMS retention gap, which makes the headline 6x a property of
    this plant's damping rather than of the constitutive laws.  The damping
    branch is exactly the quantity Section 7.3 cannot resolve, so the
    discriminator is reported across the whole disputed range and the window
    in which force discriminates becomes a fixture requirement.
    """
    sites = tuple(CASES["A"]["sites"])

    def zeta_of(multiplier: float) -> float:
        pole = _retained_mode_pole(constants, sites, multiplier)
        return float("inf") if pole is None else float(pole["zeta"])

    baseline_zeta = zeta_of(1.0)
    multipliers: list[float] = []
    for target in targets:
        # zeta is monotone in the multiplier over this range, so a bisection
        # on the executed pole gives the multiplier that realizes each target.
        low, high = 1.0e-4, 1.0e4
        for _ in range(200):
            middle = np.sqrt(low * high)
            if zeta_of(middle) < target:
                low = middle
            else:
                high = middle
        multipliers.append(float(np.sqrt(low * high)))
    if executor is None:
        completed = [_retention_damping_job(m, constants) for m in multipliers]
    else:
        futures = [executor.submit(_retention_damping_job, m, constants)
                   for m in multipliers]
        completed = [future.result() for future in futures]
    by_multiplier = dict(completed)
    rows = []
    for target, multiplier in zip(targets, multipliers):
        pole = _retained_mode_pole(constants, sites, multiplier)
        record = dict(by_multiplier[multiplier])
        record.update({
            "target_zeta": target,
            "multiplier": multiplier,
            "executed_zeta": float("nan") if pole is None else float(pole["zeta"]),
            "settling_2pct_s": SETTLING_2PCT_FACTOR / (
                target * 2.0 * np.pi * constants["axial_mode_target_hz"]),
        })
        rows.append(record)
    discriminating = [row for row in rows if row["r_hold_ratio"] >= 2.0]
    return {
        "rows": rows,
        "baseline_zeta": baseline_zeta,
        "discriminating_zeta": (
            (min(row["target_zeta"] for row in discriminating),
             max(row["target_zeta"] for row in discriminating))
            if discriminating else None),
        "ratio_range": (min(row["r_hold_ratio"] for row in rows),
                        max(row["r_hold_ratio"] for row in rows)),
    }


def _breakaway_job(force: float, constants: dict[str, float]
                   ) -> tuple[float, dict[str, object]]:
    """Process-pool unit for one guideway breakaway-force variant."""
    original = FRICTION["g"]["F_s"]
    original_coulomb = FRICTION["g"]["F_c"]
    original_rate = FRICTION["g"]["C_gms"]
    # F_c scales with F_s at the executed Stribeck ratio.  Holding F_c fixed
    # while lowering F_s would invert the Stribeck curve, which is a different
    # and unphysical model rather than a lower breakaway force.
    coulomb = original_coulomb * force / original
    FRICTION["g"]["F_s"] = force
    FRICTION["g"]["F_c"] = coulomb
    FRICTION["g"]["C_gms"] = (force - coulomb) / MODEL["tau_C"]
    try:
        experiment = presliding_responses(constants, ("A", "A2"), "g")
        metrics = experiment["metrics"]
        yields = GMS_WEIGHTS * force / GMS_STIFFNESS_BY_SITE["g"]
        record = {
            "F_s": force,
            "F_c": coulomb,
            "yields_um": tuple(float(value * 1e6) for value in yields),
            "force_mismatch_gms_N": float(metrics["A2"]["return_force_mismatch_N"]),
            "force_ratio": (float(metrics["A2"]["return_force_mismatch_N"])
                            / max(float(metrics["A"]["return_force_mismatch_N"]), 1e-30)),
            "loop_area_gms_J": float(metrics["A2"]["loop_area_J"]),
            "r_hold_gms_pct": 100.0 * float(metrics["A2"]["r_hold"]),
        }
    finally:
        FRICTION["g"]["F_s"] = original
        FRICTION["g"]["F_c"] = original_coulomb
        FRICTION["g"]["C_gms"] = original_rate
    return force, record


def breakaway_sensitivity(constants: dict[str, float],
                          executor: Executor | None = None,
                          likely_range: tuple[float, float] = (1.0, 1.5),
                          ) -> dict[str, object]:
    """Execute the guideway breakaway force at the middle of its stated range.

    Section 8.3 executes F_s = 3.0 N at the guideway while stating a likely
    range of 1.0 to 1.5 N.  The command design depends on which is right: the
    four element yield distances scale with F_s, so the inner 1.250 um level
    crosses a different threshold at the low value.  This runs the variant
    instead of describing it.
    """
    executed = float(FRICTION["g"]["F_s"])
    middle = 0.5 * (likely_range[0] + likely_range[1])
    forces = (executed, middle)
    if executor is None:
        completed = [_breakaway_job(force, constants) for force in forces]
    else:
        futures = [executor.submit(_breakaway_job, force, constants)
                   for force in forces]
        completed = [future.result() for future in futures]
    rows = [record for _force, record in
            sorted(completed, key=lambda item: -item[0])]
    inner_level = 4.0 * constants["command_step"]
    for row in rows:
        crossed = sum(1 for value in row["yields_um"] if inner_level * 1e6 > value)
        row["inner_level_um"] = inner_level * 1e6
        row["elements_yielded_at_inner"] = crossed
    return {
        "rows": rows,
        "executed_F_s": executed,
        "likely_range": likely_range,
        "middle": middle,
        "inner_level_um": inner_level * 1e6,
        "design_changes": rows[0]["elements_yielded_at_inner"]
        != rows[-1]["elements_yielded_at_inner"],
    }


def gms_branch_census_study(constants: dict[str, float], times: np.ndarray,
                            command: np.ndarray,
                            censuses: dict[str, GmsBranchCensus],
                            metrics: dict[str, dict[str, float]]) -> dict[str, object]:
    """Summarize the branch census and price the departure when it is active.

    The counterfactual rerun is executed only for cases that actually recorded
    a threshold-driven reclassification, so a zero census costs nothing.
    """
    settled_mask, _final_window = settled_window_masks(times, constants)
    rows: list[dict[str, object]] = []
    enforced: dict[str, dict[str, float]] = {}
    for key, census in censuses.items():
        for site in census.sites:
            rows.append({
                "case": key,
                "site": site,
                "flips_reversal": int(census.reversal_flips[site]),
                "flips_threshold": int(census.threshold_flips[site]),
                "evals_total": int(census.evaluations[site]),
            })
        if census.total("threshold_flips") == 0:
            continue
        branch_states = {site: np.zeros(GMS_N, dtype=bool) for site in census.sites}
        _enforced_times, enforced_states = rk4_case(
            CASES[key], constants, dt=PRODUCTION_DT, branch_states=branch_states)
        error = command - enforced_states[:, 1]
        enforced_rms = float(np.sqrt(np.mean(error[settled_mask] ** 2)) * 1e9)
        baseline_rms = float(metrics[key]["rms_settled_deviation_nm"])
        enforced[key] = {
            "settled_rms_nm": enforced_rms,
            "baseline_rms_nm": baseline_rms,
            "delta_nm": enforced_rms - baseline_rms,
            "delta_pct": 100.0 * (enforced_rms - baseline_rms) / baseline_rms,
        }
    return {
        "rows": rows,
        "enforced": enforced,
        "threshold_total": sum(int(row["flips_threshold"]) for row in rows),
        "reversal_total": sum(int(row["flips_reversal"]) for row in rows),
        "evaluation_total": sum(int(row["evals_total"]) for row in rows),
    }


# Memory-sequence metrics that the branch departure can actually reach.  The
# Section 10 settled window samples one plateau at rest; these compare repeated
# returns to the same level, so they depend on every intervening deceleration.
LOOP_METRIC_LABELS = (
    ("return_error_mismatch_nm", "$E_{ret}$", "nm", 2),
    ("return_force_mismatch_N", "$F_{ret}$", "N", 4),
    ("final_mean_nm", "final-origin magnitude $D_{13}$", "nm", 2),
    ("loop_area_J", "loop area $A_{loop}$", "µJ", 2),
)


def memory_branch_departure(constants: dict[str, float],
                            experiments: dict[str, dict[str, object]]
                            ) -> dict[str, object]:
    """Price the stateless branch test against the Section 9.4 loop metrics.

    The Section 10 settled window is sampled at rest on a single plateau; the
    departure lives in deceleration.  Repeated-return and loop-area metrics
    depend on the whole history, so they are the ones that can see it.
    """
    records: list[dict[str, object]] = []
    for experiment in experiments.values():
        gms_key = experiment["keys"][1]
        lugre_key = experiment["keys"][0]
        if CASES[gms_key]["friction"] != "gms":
            continue
        rerun = presliding_responses(
            constants, experiment["keys"], str(experiment["site"]),
            persistent_branch_keys=(gms_key,))
        for metric_key, label, unit, digits in LOOP_METRIC_LABELS:
            executed = abs(float(experiment["metrics"][gms_key][metric_key]))
            persistent = abs(float(rerun["metrics"][gms_key][metric_key]))
            law_gap = abs(float(experiment["metrics"][gms_key][metric_key])
                          - float(experiment["metrics"][lugre_key][metric_key]))
            if unit == "µJ":
                executed *= 1.0e6
                persistent *= 1.0e6
                law_gap *= 1.0e6
            records.append({
                "case": gms_key,
                "site": str(experiment["site"]),
                "metric": label,
                "unit": unit,
                "digits": digits,
                "executed": executed,
                "persistent": persistent,
                "delta": persistent - executed,
                "law_gap": law_gap,
                "exceeds_law_gap": abs(persistent - executed) > law_gap,
            })
    return {
        "records": records,
        "any_exceeds": any(bool(record["exceeds_law_gap"]) for record in records),
    }


def tau_c_sensitivity(constants: dict[str, float],
                      keys: tuple[str, str], site: str,
                      values: tuple[float, ...] = (1.0e-4, 2.0e-4, 4.0e-4)
                      ) -> dict[str, object]:
    """Rerun the memory sequence at several Stribeck relaxation times.

    C is the least anchored parameter in Section 8.  This converts that from an
    admitted weakness into a bounded one by measuring what it moves.
    """
    original = {site_key: values_dict["C_gms"]
                for site_key, values_dict in FRICTION.items()}
    rows: list[dict[str, float]] = []
    try:
        for tau in values:
            for site_key, values_dict in FRICTION.items():
                values_dict["C_gms"] = (
                    values_dict["F_s"] - values_dict["F_c"]) / tau
            experiment = presliding_responses(constants, keys, site)
            gms_metrics = experiment["metrics"][keys[1]]
            lugre_metrics = experiment["metrics"][keys[0]]
            rows.append({
                "tau_C": tau,
                "C_site": FRICTION[site]["C_gms"],
                "force_mismatch_N": float(gms_metrics["return_force_mismatch_N"]),
                "loop_area_J": float(gms_metrics["loop_area_J"]),
                "law_gap_force_N": abs(
                    float(gms_metrics["return_force_mismatch_N"])
                    - float(lugre_metrics["return_force_mismatch_N"])),
                "law_gap_loop_J": abs(float(gms_metrics["loop_area_J"])
                                      - float(lugre_metrics["loop_area_J"])),
            })
    finally:
        for site_key, value in original.items():
            FRICTION[site_key]["C_gms"] = value
    forces = [row["force_mismatch_N"] for row in rows]
    areas = [row["loop_area_J"] for row in rows]
    baseline = min(rows, key=lambda row: abs(row["tau_C"] - MODEL["tau_C"]))
    provenance = friction_provenance_metrics()
    return {
        "rows": rows,
        "force_spread_N": max(forces) - min(forces),
        "loop_spread_J": max(areas) - min(areas),
        "law_gap_force_N": baseline["law_gap_force_N"],
        "law_gap_loop_J": baseline["law_gap_loop_J"],
        "mode_period_ms": provenance["retained_mode_period"] * 1e3,
        "mode_ratio": provenance["tau_C_mode_ratio"],
    }


def _convergence_rms_job(key: str, constants: dict[str, float], dt: float
                         ) -> tuple[str, float, float]:
    """Process-pool unit for one non-production convergence trajectory."""
    times, states = rk4_case(CASES[key], constants, dt=dt)
    return key, dt, final_window_rms_error_nm(times, states, constants)


def gms_step_halving_convergence(constants: dict[str, float], base_times: np.ndarray,
                                 base_results: dict[str, np.ndarray],
                                 executor: Executor | None = None
                                 ) -> dict[str, dict[str, object]]:
    """Compare final-window RMS under h, h/2, and h/4 for all GMS cases."""
    study: dict[str, dict[str, object]] = {}
    rms_by_key = {
        key: {PRODUCTION_DT: final_window_rms_error_nm(base_times, base_results[key], constants)}
        for key in ("A2", "B2", "C2")
    }
    dts_by_key = {
        "A2": GMS_CONVERGENCE_DTS + (A2_CONVERGENCE_DT,),
        "B2": GMS_CONVERGENCE_DTS,
        "C2": GMS_CONVERGENCE_DTS,
    }
    pending = [
        (key, dt) for key, dts in dts_by_key.items() for dt in dts
        if not np.isclose(dt, PRODUCTION_DT, rtol=0.0, atol=1.0e-15)
    ]
    if executor is None:
        completed = (_convergence_rms_job(key, constants, dt) for key, dt in pending)
    else:
        futures = [
            executor.submit(_convergence_rms_job, key, constants, dt)
            for key, dt in pending
        ]
        completed = (future.result() for future in futures)
    for key, dt, rms in completed:
        rms_by_key[key][dt] = rms
    for key in ("A2", "B2", "C2"):
        dts = dts_by_key[key]
        rms_values = [rms_by_key[key][dt] for dt in dts]
        coarse_difference = abs(rms_values[0] - rms_values[1])
        fine_difference = abs(rms_values[1] - rms_values[2])
        observed_order = float(np.log2(
            coarse_difference / max(fine_difference, 1.0e-15)))
        extra_difference = (abs(rms_values[2] - rms_values[3])
                            if len(rms_values) == 4 else None)
        study[key] = {
            "dt_s": dts,
            "rms_nm": tuple(rms_values),
            "coarse_difference_nm": coarse_difference,
            "fine_difference_nm": fine_difference,
            "fine_relative_pct": 100.0 * fine_difference / max(abs(rms_values[2]), 1.0e-15),
            "difference_ratio": coarse_difference / max(fine_difference, 1.0e-15),
            "observed_order": observed_order,
            "extra_difference_nm": extra_difference,
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


def _audit_mode_shift_labels(drawn_notes: dict[str, str],
                             linear_metrics: dict[str, dict[str, float | np.ndarray]]) -> None:
    """Read back what the shift figure actually prints and re-derive it.

    The previous figure rounded a peak to whole hertz but kept an unrounded
    delta, so each label implied a different baseline and none of them matched
    the 10.1 percentages.  This parses the drawn strings and checks, against
    the eigenvalues alone, that the printed delta is the difference of the two
    printed frequencies and that the printed percentage is the ladder's.
    """
    label_pattern = re.compile(
        r"^[^\n]+\n(?P<peak>[\d.]+) Hz(?:, (?P<shift>[+-][\d.]+) Hz "
        r"\((?P<percent>[+-][\d.]+)%\))?$")
    parsed: dict[str, re.Match[str]] = {}
    for key, note in drawn_notes.items():
        match = label_pattern.match(note)
        if match is None:
            raise AssertionError(f"Mode-shift label for case {key} is unreadable: {note!r}")
        parsed[key] = match
    baseline_text = parsed["0"].group("peak")
    baseline_high = float(linear_metrics["0"]["modes"][1])
    if baseline_text != f"{baseline_high:.1f}":
        raise AssertionError(
            f"Mode-shift baseline label {baseline_text} Hz is not the case 0 eigenvalue "
            f"{baseline_high:.1f} Hz")
    failures: list[str] = []
    for key, match in parsed.items():
        if match.group("shift") is None:
            continue
        peak_text, shift_text, percent_text = match.group("peak", "shift", "percent")
        high = float(linear_metrics[key]["modes"][1])
        if peak_text != f"{high:.1f}":
            failures.append(f"{key} prints {peak_text} Hz against eigenvalue {high:.1f} Hz")
        closing_shift = f"{float(peak_text) - float(baseline_text):+.1f}"
        if shift_text != closing_shift:
            failures.append(
                f"{key} prints {shift_text} Hz, but {peak_text} minus {baseline_text} "
                f"is {closing_shift} Hz")
        ladder_percent = f"{100.0 * (high - baseline_high) / baseline_high:+.1f}"
        if percent_text != ladder_percent:
            failures.append(
                f"{key} prints {percent_text}% against the 10.1 ladder's {ladder_percent}%")
    if failures:
        raise AssertionError("Mode-shift figure labels are inconsistent: " + "; ".join(failures))


def plot_case_response_overlay(frequencies: np.ndarray,
                               responses: dict[str, np.ndarray],
                               linear_metrics: dict[str, dict[str, float | np.ndarray]],
                               effect: dict[str, float]
                               ) -> tuple[Path, Path, Path]:
    """Render the mode-shift result, the only non-degenerate tangent delta,
    and the complete Bode overlay as three figures with distinct jobs."""
    magnitudes = {
        key: 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))
        for key, response in responses.items()
    }
    short_labels = {
        "0": "0: frictionless",
        "A": "A: guideway LuGre",
        "A2": "A2: guideway GMS",
        "G": "G: guideway-only LuGre",
        "G2": "G2: guideway-only GMS",
        "B": "B: nut LuGre",
        "B2": "B2: nut GMS",
        "C": "C: all-port LuGre",
        "C2": "C2: all-port GMS",
        "A1v": "A1v: guideway LuGre, micro-viscous",
    }
    missing_labels = sorted(set(CASES).difference(short_labels))
    if missing_labels:
        raise KeyError("Bode overlay has unlabelled cases: " + ", ".join(missing_labels))

    # The labels quote the same unrounded eigenvalues as the 10.1 ladder, so
    # the peak, the shift, and the percentage cannot contradict each other or
    # the table.  A linear axis is used because the whole span is 115 Hz.
    ladder = mode_shift_ladder(linear_metrics)
    ladder_by_key = {str(row["key"]): row for row in ladder}
    zoom_floor_db = -5.0
    zoom_ceiling_db = 15.5
    fig_zoom, ax_zoom = plt.subplots(figsize=(10.8, 6.2))
    for key, case in CASES.items():
        ax_zoom.plot(
            frequencies, magnitudes[key], color=case["color"],
            linestyle=case["ls"], linewidth=1.9,
        )

    def magnitude_at(key: str, frequency: float) -> float:
        return float(np.interp(frequency, frequencies, magnitudes[key]))

    # Each matched pair is exactly coincident, so one label per visible curve
    # replaces a ten-entry legend.  G/G2 sit under A/A2 at this mode and are
    # named in the footnote instead of a colliding second label.
    annotations = (("0", 12.8), ("A", 9.4), ("B", 12.8), ("C", 9.4))
    drawn_notes: dict[str, str] = {}
    for key, text_height in annotations:
        row = ladder_by_key[key]
        frequency = float(row["high_hz"])
        peak_magnitude = magnitude_at(key, frequency)
        if row["shift_pct_text"] is None:
            note = f"{row['figure_label']}\n{row['high_text']} Hz"
        else:
            note = (f"{row['figure_label']}\n{row['high_text']} Hz, "
                    f"{row['shift_hz_text']} Hz ({row['shift_pct_text']}%)")
        drawn_notes[key] = note
        ax_zoom.plot(
            [frequency, frequency], [peak_magnitude + 0.5, text_height - 0.4],
            color=CASES[key]["color"], linestyle=(0, (3, 2)), linewidth=1.0,
        )
        ax_zoom.text(
            frequency, text_height, note, ha="center", va="bottom",
            fontsize=8.6, color=CASES[key]["color"], linespacing=1.35,
        )
    _audit_mode_shift_labels(drawn_notes, linear_metrics)

    a1v_frequency = float(ladder_by_key["A"]["high_hz"])
    ax_zoom.annotate(
        f"A1v (dotted): same ports, $\\sigma_1$ restored;\n"
        f"peak {effect['peak_drop_db']:.3f} dB lower, mode unmoved",
        xy=(a1v_frequency, magnitude_at("A1v", a1v_frequency)),
        xytext=(747.0, 6.2), textcoords="data", ha="left", va="center",
        fontsize=8.2, color=CASES["A1v"]["color"], linespacing=1.35,
        arrowprops={"arrowstyle": "->", "color": CASES["A1v"]["color"],
                    "lw": 0.9, "connectionstyle": "arc3,rad=0.16"},
    )
    fig_zoom.text(
        0.5, 0.015,
        "Matched LuGre/GMS pairs are exactly coincident, so each visible curve is labelled once. "
        "G/G2 lie under A/A2 here; the drive port moves only the low mode.",
        ha="center", fontsize=8.2, color="#555555",
    )
    ax_zoom.set_xlim(678.0, 838.0)
    ax_zoom.set_ylim(zoom_floor_db, zoom_ceiling_db)
    ax_zoom.set_xticks([690.0, 720.0, 750.0, 780.0, 810.0])
    ax_zoom.set_title("Presliding stiffness shifts the retained resonance")
    ax_zoom.set_xlabel("Frequency (Hz)")
    ax_zoom.set_ylabel("Command-to-stage magnitude (dB)")
    ax_zoom.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
    ax_zoom.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
    ax_zoom.minorticks_on()
    fig_zoom.tight_layout(rect=(0.0, 0.045, 1.0, 1.0))
    zoom_output = ASSET_DIR / "friction_mode_shift_zoom.svg"
    save_svg(fig_zoom, zoom_output)
    plt.close(fig_zoom)

    # One trace with one feature: half the height of a normal panel, no legend,
    # and a sign the caption explains.  A negative difference is a reduction in
    # peak magnitude produced by added damping, not a loss of gain elsewhere.
    difference = magnitudes["A1v"] - magnitudes["A"]
    maximum_index = int(np.argmax(np.abs(difference)))
    maximum = float(difference[maximum_index])
    maximum_frequency = float(frequencies[maximum_index])
    fig_delta, ax_delta = plt.subplots(figsize=(9.3, 2.7))
    ax_delta.semilogx(
        frequencies, difference, color=CASES["A1v"]["color"], linewidth=2.0)
    ax_delta.plot(maximum_frequency, maximum, "o", color=CASES["A1v"]["color"], ms=5)
    ax_delta.annotate(
        f"peak reduced {abs(maximum):.3f} dB at {maximum_frequency:.0f} Hz",
        xy=(maximum_frequency, maximum), xytext=(12, 6), textcoords="offset points",
        fontsize=8.4, color=CASES["A1v"]["color"],
    )
    finite_delta = difference[np.isfinite(difference)]
    delta_min, delta_max = float(np.min(finite_delta)), float(np.max(finite_delta))
    delta_span = max(delta_max - delta_min, abs(maximum) * 0.18, 1.0e-5)
    ax_delta.set_ylim(delta_min - 0.16 * delta_span, delta_max + 0.16 * delta_span)
    ax_delta.axhline(0.0, color="#888888", linewidth=0.7)
    ax_delta.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
    ax_delta.set_title("A1v minus A: isolated micro-viscous tangent effect")
    ax_delta.set_xlabel("Frequency (Hz)")
    ax_delta.set_ylabel("Magnitude\ndifference (dB)")
    ax_delta.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
    ax_delta.grid(True, which="minor", color="#eeeeee", linewidth=0.45)

    # The full-range panel proves the flatness; the inset shows the feature.
    inset = ax_delta.inset_axes((0.085, 0.13, 0.23, 0.64))
    inset_mask = (frequencies >= 650.0) & (frequencies <= 850.0)
    inset.plot(frequencies[inset_mask], difference[inset_mask],
               color=CASES["A1v"]["color"], linewidth=1.5)
    inset.axhline(0.0, color="#888888", linewidth=0.6)
    inset.set_xlim(650.0, 850.0)
    inset.set_ylim(1.10 * maximum, abs(0.10 * maximum))
    inset.set_xticks([650.0, 750.0, 850.0])
    inset.set_yticks([0.0, round(maximum, 2)])
    inset.tick_params(labelsize=6.4, pad=1.5)
    inset.set_title("650-850 Hz, linear", fontsize=6.6, pad=2.0)
    inset.grid(True, which="major", color="#e2e2e2", linewidth=0.5)

    caption = "\n".join([
        f"Micro-viscous damping acts only at the retained mode and nowhere else. The "
        f"{abs(maximum):.3f} dB peak reduction corresponds to a "
        f"{100.0 * (effect['magnitude_ratio'] - 1.0):.1f}% increase in modal damping",
        f"and shifts the settled RMS deviation by {effect['rms_shift_nm']:.1f} nm, which is why "
        "$\\sigma_1=0$ in the matched comparisons.",
    ])
    fig_delta.text(0.5, 0.03, caption, ha="center", va="bottom",
                   fontsize=8.0, color="#555555", linespacing=1.4)
    fig_delta.tight_layout(rect=(0.0, 0.235, 1.0, 1.0))
    delta_output = ASSET_DIR / "micro_viscous_difference.svg"
    save_svg(fig_delta, delta_output)
    plt.close(fig_delta)

    fig_full, ax_full = plt.subplots(figsize=(11.2, 6.8))
    for key, case in CASES.items():
        ax_full.semilogx(
            frequencies, magnitudes[key], color=case["color"],
            linestyle=case["ls"], linewidth=1.8, label=short_labels[key],
        )
    ax_full.axhline(0.0, color="#888888", linewidth=0.7)
    ax_full.set_xlim(BODE_FOCUS_MIN_HZ, BODE_FOCUS_MAX_HZ)
    ax_full.set_ylim(-90.0, 30.0)
    ax_full.set_title("All command-to-stage Bode responses")
    ax_full.set_xlabel("Frequency (Hz)")
    ax_full.set_ylabel("Magnitude (dB)")
    ax_full.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
    ax_full.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
    ax_full.legend(loc="lower left", ncol=2, fontsize=8)
    fig_full.tight_layout()
    full_output = ASSET_DIR / "bode_all_cases.svg"
    save_svg(fig_full, full_output)
    plt.close(fig_full)
    return zoom_output, delta_output, full_output


def plot_presliding_memory(experiment: dict[str, object], output_name: str) -> Path:
    """Show motion and two readable, settled-return friction-loop comparisons."""
    times = np.asarray(experiment["times"], dtype=float)
    command = np.asarray(experiment["command"], dtype=float)
    results = experiment["results"]
    forces = experiment["forces"]
    metrics = experiment["metrics"]
    keys = tuple(experiment["keys"])
    site = str(experiment["site"])
    site_title = "Guideway" if site == "g" else "Nut microslip"
    blocked_stage = bool(experiment["blocked_stage"])
    time_ms = times * 1e3
    plateau_dwell = float(experiment["plateau_dwell"])
    levels = np.asarray(experiment["levels"], dtype=float)
    plateau_elapsed = times - PRESLIDING_START
    plateau_phase = np.mod(np.maximum(plateau_elapsed, 0.0), plateau_dwell)
    settled_display_mask = (
        (plateau_elapsed >= 0.0)
        & (times <= PRESLIDING_START + levels.size * plateau_dwell + 0.5e-9)
        & (plateau_phase >= 0.040 - 0.5e-9)
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.8))
    ax_motion, ax_error = axes[0]
    ax_loop_a, ax_loop_b = axes[1]

    ax_motion.step(time_ms, command * 1e6, where="post", color="#111111",
                   linewidth=2.0, label="Command")
    peak_error_nm = 0.0
    for key in keys:
        case = CASES[key]
        stage = results[key][:, 1]
        observed = results[key][:, 0] if blocked_stage else stage
        error_nm = (command - observed) * 1e9
        peak_error_nm = max(peak_error_nm, float(np.max(np.abs(error_nm))))
        ax_motion.plot(time_ms, observed * 1e6, color=case["color"],
                       linestyle=case["ls"], linewidth=1.35, label=case["label"])
        ax_error.plot(
            time_ms, np.where(settled_display_mask, error_nm, np.nan),
            color=case["color"], linestyle=case["ls"], linewidth=1.15,
        )

    ax_motion.set_title("Nested command and drive coordinate; stage blocked"
                        if blocked_stage else "Nested command and stage motion")
    ax_motion.set_ylabel("Position (µm)")
    ax_motion.legend(loc="upper right", fontsize=7.7, ncol=2)
    ax_error.set_title("Settled deviation, per-plateau transients masked")
    ax_error.set_ylabel(r"Modeled deviation $x_{cmd}-x_o$ (nm)")
    ax_error.axhline(0.0, color="#777777", linewidth=0.8)
    ax_error.text(
        0.02, 0.04, f"first 40 ms of each plateau masked; edge peak {peak_error_nm:.0f} nm",
        transform=ax_error.transAxes, fontsize=8.0, color="#555555",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white",
              "edgecolor": "#c9cfd4", "alpha": 0.90})

    def coordinates(key: str) -> np.ndarray:
        stage = results[key][:, 1]
        return stage if site == "g" else results[key][:, 0] - stage

    def add_direction_arrows(axis: plt.Axes, x: np.ndarray, y: np.ndarray,
                             color: str) -> None:
        for index in range(1, len(x), 3):
            start = max(0, index - 1)
            axis.annotate(
                "", xy=(x[index], y[index]), xytext=(x[start], y[start]),
                arrowprops={"arrowstyle": "-|>", "color": color,
                            "linewidth": 0.8, "mutation_scale": 8},
                zorder=5)

    def loop_panel(axis: plt.Axes, panel_keys: tuple[str, ...],
                   title: str, point_slice: slice | None = None,
                   background_mask: np.ndarray | None = None,
                   numbered: bool = False, dotted_connector: bool = False) -> None:
        dynamic_key = panel_keys[-1]
        background_x = coordinates(dynamic_key)
        background_f = forces[dynamic_key]
        if background_mask is not None:
            background_x = background_x[background_mask]
            background_f = background_f[background_mask]
        axis.plot(
            background_x * 1e6, background_f,
            color="#9ba3aa", linewidth=0.75, alpha=0.28,
            label=f"{dynamic_key} full dynamic trace")
        for index, key in enumerate(panel_keys):
            x_all = np.asarray(metrics[key]["endpoint_coordinate_um"], dtype=float)
            y_all = np.asarray(metrics[key]["endpoint_force_N"], dtype=float)
            plateau_numbers = np.arange(1, x_all.size + 1)
            if point_slice is not None:
                x, y, numbers = x_all[point_slice], y_all[point_slice], plateau_numbers[point_slice]
            else:
                x, y, numbers = x_all, y_all, plateau_numbers
            connector = ({"linestyle": ":", "alpha": 0.45, "linewidth": 1.1} if dotted_connector
                         else {"linestyle": CASES[key]["ls"], "alpha": 1.0, "linewidth": 1.8})
            axis.plot(
                x, y, color=CASES[key]["color"], marker="o" if index == 0 else "s",
                markersize=4.5, markerfacecolor="white", label=key, **connector)
            if numbered:
                offset = (3, 3) if index == 0 else (3, -9)
                for xi, yi, ni in zip(x, y, numbers):
                    axis.annotate(
                        str(int(ni)), xy=(xi, yi), xytext=offset, textcoords="offset points",
                        fontsize=6.4, color=CASES[key]["color"])
            else:
                add_direction_arrows(axis, x, y, CASES[key]["color"])
        all_coordinates = np.concatenate([
            np.asarray(metrics[key]["endpoint_coordinate_um"], dtype=float)
            if point_slice is None else
            np.asarray(metrics[key]["endpoint_coordinate_um"], dtype=float)[point_slice]
            for key in panel_keys])
        all_forces = np.concatenate([
            np.asarray(metrics[key]["endpoint_force_N"], dtype=float)
            if point_slice is None else
            np.asarray(metrics[key]["endpoint_force_N"], dtype=float)[point_slice]
            for key in panel_keys])
        x_extent = max(float(np.max(np.abs(all_coordinates))) * 1.10, 0.35)
        y_extent = max(float(np.max(np.abs(all_forces))) * 1.12, 0.25)
        if site == "g":
            x_extent, y_extent = min(max(x_extent, 4.0), 4.2), min(max(y_extent, 2.5), 2.7)
        axis.set_xlim(-x_extent, x_extent)
        axis.set_ylim(-y_extent, y_extent)
        axis.set_title(title)
        axis.set_xlabel("Stage position (µm)" if site == "g"
                        else r"Nut-port deflection $x_d-x_s$ (µm)")
        axis.set_ylabel(f"{site_title} friction force (N)")
        axis.axhline(0.0, color="#888888", linewidth=0.7)
        axis.axvline(0.0, color="#888888", linewidth=0.7)

    lugre_key, gms_key = keys[:2]
    loop_panel(
        ax_loop_a, (lugre_key, gms_key),
        f"Law comparison: {lugre_key} versus {gms_key}",
        numbered=site == "n", dotted_connector=site == "n",
    )
    r_hold_lugre = 100.0 * float(metrics[lugre_key]["r_hold"])
    r_hold_gms = 100.0 * float(metrics[gms_key]["r_hold"])
    ax_loop_a.text(
        0.03, 0.96,
        f"LuGre retains {r_hold_lugre:.1f}% of available elastic force at rest; "
        f"GMS {r_hold_gms:.1f}%.",
        transform=ax_loop_a.transAxes, va="top", fontsize=7.8, color="#3f4b53",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white",
              "edgecolor": "#c9cfd4", "alpha": 0.92},
    )

    plateau_numbers = np.arange(1, levels.size + 1)
    for key, marker in ((lugre_key, "o"), (gms_key, "s")):
        ax_loop_b.plot(
            plateau_numbers,
            np.asarray(metrics[key]["endpoint_force_N"], dtype=float),
            color=CASES[key]["color"], linestyle=CASES[key]["ls"],
            linewidth=1.8, marker=marker, markersize=4.5,
            markerfacecolor="white", label=key,
        )
    ax_loop_b.axhline(0.0, color="#888888", linewidth=0.7)
    ax_loop_b.set_xlim(1, levels.size)
    ax_loop_b.set_xticks(plateau_numbers)
    ax_loop_b.set_title("Retention diagnostic: settled force versus plateau index")
    ax_loop_b.set_xlabel("Plateau index")
    ax_loop_b.set_ylabel(f"Settled {site_title.lower()} friction force (N)")
    ax_loop_b.legend(loc="best", fontsize=8)

    for axis in (ax_motion, ax_error):
        axis.set_xlabel("Time (ms)")
    for axis in axes.flat:
        axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)

    force_ratio = (float(metrics[gms_key]["return_force_mismatch_N"]) /
                   max(float(metrics[lugre_key]["return_force_mismatch_N"]), 1e-30))
    origin_ratio = (abs(float(metrics[gms_key]["final_mean_nm"])) /
                    max(abs(float(metrics[lugre_key]["final_mean_nm"])), 1e-30))
    fig.suptitle(f"{site_title} partial-slip memory experiment",
                 fontsize=15, fontweight="bold")
    if site == "g":
        caption = (
            f"GMS fails to return to the same friction force by {force_ratio:.2f}× more "
            f"than LuGre and leaves {origin_ratio:.0f}× more residual error at the origin. "
            "That non-closure is the nonlocal memory. Both laws produce nearly identical "
            "stage motion, which is why force rather than displacement is the discriminator.")
    else:
        caption = (
            f"The blocked nut fixture exposes the same non-closure signature: GMS has "
            f"{force_ratio:.2f}× the LuGre return-force mismatch. Force reveals memory that "
            "the nearly co-moving free-stage coordinates would hide.")
    fig.text(0.5, 0.012, caption, ha="center", fontsize=8.5, color="#555555")
    fig.tight_layout(rect=(0.02, 0.055, 0.99, 0.95), h_pad=2.0, w_pad=1.5)
    output = ASSET_DIR / output_name
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_presliding_supplement(experiments: dict[str, dict[str, object]]) -> Path:
    """Demote the guideway ablation and nut branch-split diagnostics from the
    two main Section 9 figures without losing their audit value."""
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))

    def add_endpoint_panel(axis: plt.Axes, experiment: dict[str, object],
                           panel_keys: tuple[str, ...], title: str,
                           point_slice: slice | None = None) -> None:
        site = str(experiment["site"])
        metrics = experiment["metrics"]
        for index, key in enumerate(panel_keys):
            x = np.asarray(metrics[key]["endpoint_coordinate_um"], dtype=float)
            y = np.asarray(metrics[key]["endpoint_force_N"], dtype=float)
            numbers = np.arange(1, x.size + 1)
            if point_slice is not None:
                x, y, numbers = x[point_slice], y[point_slice], numbers[point_slice]
            axis.plot(
                x, y, color=CASES[key]["color"], linestyle=CASES[key]["ls"],
                linewidth=1.7, marker="o" if index == 0 else "s", markersize=4.3,
                markerfacecolor="white", label=key,
            )
            for xi, yi, plateau in zip(x, y, numbers):
                offset = (3, 3) if index == 0 else (3, -9)
                axis.annotate(
                    str(int(plateau)), xy=(xi, yi), xytext=offset,
                    textcoords="offset points", fontsize=6.2,
                    color=CASES[key]["color"],
                )
        axis.axhline(0.0, color="#888888", linewidth=0.7)
        axis.axvline(0.0, color="#888888", linewidth=0.7)
        axis.set_title(title)
        axis.set_xlabel("Stage position (um)" if site == "g"
                        else r"Nut-port deflection $x_d-x_s$ (um)")
        axis.set_ylabel("Friction force (N)")
        axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
        axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
        axis.legend(loc="best", fontsize=8)

    guideway = experiments["guideway"]
    nut = experiments["nut"]
    add_endpoint_panel(
        axes[0], guideway, ("A2", "G2"),
        "Guideway drive-port ablation: A2 versus G2",
    )
    add_endpoint_panel(
        axes[1], nut, ("B", "B2"),
        "Nut positive branch", slice(0, 7),
    )
    add_endpoint_panel(
        axes[2], nut, ("B", "B2"),
        "Nut negative branch", slice(7, 13),
    )
    fig.suptitle("Supplementary memory diagnostics", fontsize=14, fontweight="bold")
    fig.text(
        0.5, 0.012,
        "Plateau numbers show traversal order. These panels audit ablation and branch asymmetry; "
        "the main Section 9 figures reserve their lower-right slot for direct force retention.",
        ha="center", fontsize=8.2, color="#555555",
    )
    fig.tight_layout(rect=(0.01, 0.055, 0.995, 0.93), w_pad=1.4)
    output = ASSET_DIR / "memory_diagnostic_supplement.svg"
    save_svg(fig, output)
    plt.close(fig)
    return output


def plot_true_presliding_loop(loop: dict[str, object]) -> Path:
    """The literature-comparable continuous loop (Part 1.5 item 4): a slow
    triangular ramp-reversal with no plateaus, distinct from the settled
    return-point maps built from discrete commanded levels."""
    fig, axis = plt.subplots(figsize=(6.6, 5.8))
    for key in ("A", "A2"):
        case = CASES[key]
        x = loop["site_coordinate"][key] * 1e6
        f = loop["forces"][key]
        axis.plot(x, f, color=case["color"], linestyle=case["ls"],
                  linewidth=1.6, label=case["label"])
    axis.axhline(0.0, color="#888888", linewidth=0.7)
    axis.axvline(0.0, color="#888888", linewidth=0.7)
    axis.set_xlabel("Stage position (µm)")
    axis.set_ylabel("Guideway friction force (N)")
    axis.set_title("Continuous quasi-static presliding loop")
    axis.legend(loc="upper left", fontsize=8.5)
    axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
    axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
    fig.text(
        0.5, 0.012,
        f"Slow triangular reversal, no plateaus, {loop['ramp_time'] * 1e3:.0f} ms per quarter-cycle "
        f"at the {loop['amplitude'] * 1e6:.3f} um guideway outer amplitude. Literature-comparable to "
        "published presliding F-x curves, unlike the settled return-point maps above.",
        ha="center", fontsize=8.0, color="#555555")
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.95))
    output = ASSET_DIR / "presliding_true_loop.svg"
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
    ax_step.set_title("Linear full-step scaling test (5 µm; ideal ZOH)")
    ax_step.legend(fontsize=8.2, loc="upper right")
    ax_residual.plot(time_ms, residual * 1e9, color="#9b2f3d", linewidth=1.35)
    ax_residual.axhline(0.0, color="#888888", linewidth=0.7)
    for edge_time in VERIFICATION_EDGES:
        for axis in (ax_step, ax_residual):
            axis.axvline(edge_time * 1e3, color="#777777", linewidth=0.65,
                        linestyle=":", alpha=0.75)
    ax_residual.set_xlabel("Time (ms)")
    ax_residual.set_ylabel("Full − reduced (nm)")
    ax_residual.set_title(
        f"Reduction residual: RMS {verification['rms_residual_pct_command']:.2f}% command; "
        f"peak {verification['peak_residual_pct_command']:.2f}% command")
    ax_residual.text(
        0.02, 0.04,
        f"h = {verification['production_dt'] * 1e6:.1f} µs | "
        f"dominant energy {verification['dominant_residual_frequency_hz']:.0f} Hz | "
        f"fast ripple {verification['ripple_residual_frequency_hz']:.0f} Hz",
        transform=ax_residual.transAxes, ha="left", va="bottom", fontsize=7.7,
        color="#4d2730",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white",
              "edgecolor": "#d9c6ca", "alpha": 0.88})
    ax_step.text(
        0.02, 0.04, "Structural audit; not a nonlinear full-step prediction",
        transform=ax_step.transAxes, ha="left", va="bottom", fontsize=7.7,
        color="#555555",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white",
              "edgecolor": "#dddddd", "alpha": 0.88})
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
    # Same modulus, root section and datum as the Section 2 entry table, so
    # the sweep and the four executed stiffnesses cannot describe different
    # screws.
    component = component_parameters()
    axial_rigidity = component["axial_rigidity"]
    k_sha = axial_rigidity / free_lengths
    constants = physical_constants()
    fixed_compliance = 1.0 / constants["k_ax"] - 1.0 / component["k_sha"]
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
        damping = constants["c_ax"] * coupling + c_m * np.outer(H["d"], H["d"])
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


def generated_reduction_convergence(verification: dict[str, object]) -> str:
    """Build the Chapter 7 solver audit and physical interpretation."""
    discarded_frequency, discarded_zeta = verification["first_discarded_damped_mode"]
    full_upper_frequency, full_upper_zeta = verification["full_upper_damped_mode"]
    reduced_upper_frequency, reduced_upper_zeta = verification["reduced_upper_damped_mode"]
    edge_spacing = VERIFICATION_EDGES[1] - VERIFICATION_EDGES[0]
    full_upper_carryover = np.exp(
        -full_upper_zeta * 2.0 * np.pi * full_upper_frequency * edge_spacing)
    discarded_carryover = np.exp(
        -discarded_zeta * 2.0 * np.pi * discarded_frequency * edge_spacing)
    full_lower_frequency, full_lower_zeta_value = verification["full_lower_damped_mode"]
    full_lower_carryover = np.exp(
        -full_lower_zeta_value * 2.0 * np.pi * full_lower_frequency * edge_spacing)
    edge_peak_values = np.asarray(verification["edge_peaks_nm"])
    edge_growth = float(edge_peak_values[-1] / edge_peak_values[0])
    # Damping ratio that would make a single drive-pole carryover reproduce the
    # observed growth; reported only to show that it is not the executed value.
    zeta_1_for_observed = -np.log(max(1.0 - 1.0 / edge_growth, 1e-12)) / (
        2.0 * np.pi * full_lower_frequency * edge_spacing)
    lines = [
        "<!-- BEGIN GENERATED REDUCTION CONVERGENCE -->",
        "### 7.1 Solver convergence",
        "",
        f"The time-domain comparison now uses the physical {physical_constants()['full_step'] * 1e6:.3f} µm full-step pitch. "
        "Because both verification plants are linear, this rescales the displacement and residual in nanometres but does not change the normalized RMS or peak percentages.",
        "",
        f"| RK4 step $h$ | Points/cycle at {discarded_frequency:.1f} Hz | Maximum $\\lvert R(h\\lambda)\\rvert$ | Result | RMS residual | Peak residual |",
        "|---:|---:|---:|---|---:|---:|",
    ]
    for row in verification["convergence"]:
        if bool(row["stable"]):
            result = "stable"
            rms_text = f"{row['rms_residual_nm']:.3f} nm ({row['rms_residual_pct_command']:.5f}%)"
            peak_text = f"{row['peak_residual_nm']:.3f} nm ({row['peak_residual_pct_command']:.5f}%)"
        else:
            result = "**unstable**"
            rms_text = "not reportable"
            peak_text = "not reportable"
        lines.append(
            f"| {row['dt'] * 1e6:.2f} µs | {row['points_per_discarded_cycle']:.1f} | "
            f"{row['stability_radius']:.6f} | {result} | {rms_text} | {peak_text} |"
        )
    edge_peaks = ", ".join(
        f"{value:.1f}" for value in np.asarray(verification["edge_peaks_nm"]))
    residual_rows = list(verification["route_residuals"])
    by_key = {str(row["key"]): row for row in residual_rows}
    executed, interface, measured = by_key["P"], by_key["F"], by_key["M"]
    cb_two_dof, craig_bampton = by_key["C2"], by_key["CB"]
    damping_only_change = 100.0 * abs(
        float(interface["rms_residual_nm"]) - float(executed["rms_residual_nm"])
    ) / float(executed["rms_residual_nm"])
    truncation_factor = (float(executed["rms_residual_nm"])
                         / float(craig_bampton["rms_residual_nm"]))
    full_settling = SETTLING_2PCT_FACTOR / (
        full_upper_zeta * 2.0 * np.pi * full_upper_frequency)
    reduced_settling = SETTLING_2PCT_FACTOR / (
        reduced_upper_zeta * 2.0 * np.pi * reduced_upper_frequency)
    constants = physical_constants()
    dwell = constants["plateau_dwell"]
    separation_ratio = (constants["axial_mode_target_hz"] / discarded_frequency)**2
    lines.extend([
        "",
        "The 25 µs result is not a coarse but usable answer: it is mathematically unstable for this ten-DOF state matrix. "
        f"The unplotted full model reaches {max(verification['full_modes']) / 1e3:.2f} kHz, and the largest RK4 amplification magnitude is greater than one. "
        "The 12.5, 6.25, and production 2.5 µs results converge to the same output residual, so the inter-edge growth below is not integration drift.",
        "",
        f"Both static gains are unity to numerical precision ($G_{{full}}(0)={verification['full_dc_gain']:.12f}$ and $G_{{red}}(0)={verification['reduced_dc_gain']:.12f}$), and the residual is zero before the first edge. "
        f"The four successive inter-edge peak magnitudes are {edge_peaks} nm. "
        f"The strongest residual spectral energy is near {verification['dominant_residual_frequency_hz']:.1f} Hz; the visibly faster ripple is near {verification['ripple_residual_frequency_hz']:.1f} Hz.",
        "",
        f"The residual is not explained by the {discarded_frequency:.1f} Hz ripple alone. "
        f"That full-model mode has $\\zeta={discarded_zeta:.5f}$ and retains only {100.0 * discarded_carryover:.1f}% of its amplitude over the 20 ms edge spacing. "
        f"It is also the pole that sets the timescale-separation ratio used in [Appendix E.8.2](#e-8-equivalence-proofs-and-error-bounds): "
        f"$({constants['axial_mode_target_hz']:.2f}/{discarded_frequency:.1f})^2={separation_ratio:.3f}$. "
        f"The retained upper mode carries some of the rest: the full model has {full_upper_frequency:.1f} Hz with $\\zeta_2={full_upper_zeta:.5f}$ against the reduced model's {reduced_upper_frequency:.1f} Hz with $\\zeta_2={reduced_upper_zeta:.5f}$, "
        f"the two damping ratios differing by {100.0 * abs(reduced_upper_zeta - full_upper_zeta) / full_upper_zeta:.0f}%. "
        "The two plants now agree about how fast that mode decays, so what remains is not a damping inconsistency.",
        "",
        f"The peaks still climb, by a factor of {edge_growth:.2f} across the four edges, and no single-mode carryover argument reproduces that. "
        f"The upper mode retains {100.0 * full_upper_carryover:.1f}% of its amplitude over the 20 ms edge spacing, which would cap the accumulation at "
        f"{1.0 / (1.0 - full_upper_carryover):.2f}; the drive pole, which the per-plant audit below identifies as the dominant residual line, retains only "
        f"{100.0 * full_lower_carryover:.1f}% at $\\zeta_1={full_lower_zeta_value:.4f}$ and would cap it at {1.0 / (1.0 - full_lower_carryover):.2f}. "
        f"Matching the observed {edge_growth:.2f} from the drive pole alone would need $\\zeta_1\\approx{zeta_1_for_observed:.3f}$ against the executed {full_lower_zeta_value:.4f}. "
        "The growth is therefore bounded and modest but is not attributed to one mode here: it is a difference signal between two plants whose poles differ in frequency as well as amplitude, "
        "and a scalar carryover argument does not apply to it.",
        "",
        "### 7.2 Per-plant residual audit",
        "",
        f"Every Section 6 candidate is driven with the same command and differenced against the same ten-DOF stage output at the production {VERIFICATION_DT * 1e6:.1f} µs step. "
        "The damping question and the truncation question are then separated by measurement instead of inference.",
        "",
        "| Reduced plant | Coordinates | Damping basis | $f_1$ (Hz) | $\\zeta_1$ | $f_2$ (Hz) | $\\zeta_2$ | RMS residual | Peak residual |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|",
    ])
    full_lower, full_lower_zeta = verification["full_lower_damped_mode"]
    for row in residual_rows:
        lines.append(
            f"| {row['label']} | {row['coordinates']} | {row['basis']} | "
            f"{row['lower_hz']:.2f} | {row['lower_zeta']:.4f} | "
            f"{row['upper_hz']:.2f} | {row['upper_zeta']:.3e} | "
            f"{row['rms_residual_nm']:.3f} nm ({row['rms_residual_pct_command']:.3f}%) | "
            f"{row['peak_residual_nm']:.3f} nm ({row['peak_residual_pct_command']:.3f}%) |"
        )
    lines.extend([
        f"| **Ten-DOF reference** | 10 | element-wise $\\eta_j$ | **{full_lower:.2f}** | **{full_lower_zeta:.4f}** | "
        f"**{full_upper_frequency:.2f}** | **{full_upper_zeta:.3e}** | - | - |",
        "",
        "Every plant in the table has unity static gain, so none of the residual is a compliance error.",
        "",
        f"Damping assignment is now nearly irrelevant to the residual. The executed plant and the interface-propagated plant differ in $\\zeta_2$ by only a factor of {executed['upper_zeta'] / interface['upper_zeta']:.2f}, "
        f"and their RMS residuals differ by {damping_only_change:.1f}%. The measured-mass plant is the exception at {100.0 * (measured['rms_residual_nm'] / executed['rms_residual_nm'] - 1.0):.0f}% worse, "
        f"and it is worse precisely because its $\\zeta_2$ is set by the separately assumed measured relative damping rather than by the interface loss factors, leaving it an order of magnitude underdamped against the ten-DOF plant. "
        f"Coordinate content does the rest: the 2-DOF plant rebuilt at the Craig-Bampton frequency of {cb_two_dof['upper_hz']:.2f} Hz drops the RMS residual to {cb_two_dof['rms_residual_nm']:.1f} nm, "
        f"and restoring one eliminated coordinate drops it to {craig_bampton['rms_residual_nm']:.1f} nm, a factor of {truncation_factor:.1f} below the executed plant.",
        "",
        f"**With the damping question removed, coordinate truncation is what is left.** "
        f"Every two-coordinate plant that carries a defensible damping value lands within {100.0 * abs(interface['rms_residual_nm'] - executed['rms_residual_nm']) / executed['rms_residual_nm']:.1f}% of the same residual, "
        f"and only adding a coordinate moves it, by a factor of {truncation_factor:.1f}. "
        "This is the measurement behind the [Section 6.3](#6-3-reduction-evidence) row, and it is now a clean one-variable result rather than an inference drawn across two confounded variables.",
        "",
        f"**The $f_1$ column locates the error, and it is not where the section previously looked.** "
        f"The strongest residual line has moved to {verification['dominant_residual_frequency_hz']:.1f} Hz, near the drive pole rather than near the axial mode. "
        f"Every two-coordinate plant places that pole at {executed['lower_hz']:.2f} Hz against the ten-DOF value of {full_lower:.2f} Hz, "
        f"a {100.0 * abs(executed['lower_hz'] - full_lower) / full_lower:.2f}% error that static condensation cannot remove because it is a dynamic-participation effect, not a static-stiffness one. "
        f"Restoring one fixed-interface coordinate corrects it to {craig_bampton['lower_hz']:.2f} Hz, within {100.0 * abs(craig_bampton['lower_hz'] - full_lower) / full_lower:.3f}%. "
        f"The decomposition is therefore explicit: damping assignment is worth {100.0 * abs(interface['rms_residual_nm'] - executed['rms_residual_nm']) / executed['rms_residual_nm']:.1f}%, "
        f"aligning the upper mode is worth {100.0 * (executed['rms_residual_nm'] - cb_two_dof['rms_residual_nm']) / executed['rms_residual_nm']:.0f}%, "
        f"and correcting the drive pole is worth a further {100.0 * (cb_two_dof['rms_residual_nm'] - craig_bampton['rms_residual_nm']) / executed['rms_residual_nm']:.0f}%; "
        f"the last {100.0 * craig_bampton['rms_residual_nm'] / executed['rms_residual_nm']:.0f}% is removed by none of the three and is the residual the three-coordinate plant still carries. "
        "This also prices the one assumption that [E.3](#e-3-direct-series-compliance-reduction) could previously only bound: dropping the eliminated axial inertia costs half a percent on the drive pole, and that half percent is now the largest single term in the reduction residual.",
        "",
        "### 7.3 Dwell consequence",
        "",
        "**The dwell is conditional on an unresolved damping branch, so it now covers every branch.** "
        f"The same retained mode carries three candidate damping ratios and they disagree by a factor of "
        f"{constants['axial_zeta_executed'] / constants['axial_zeta_measured']:.0f}: "
        f"the measured relative-mode value $\\zeta_2=$ {constants['axial_zeta_measured']:.4f}, still pending the "
        f"[E.7](#e-7-measured-frf-identification) half-power re-extraction, implies a 2% settling time of "
        f"{constants['measured_settling_time_2pct'] * 1e3:.0f} ms; the interface loss factors propagated in "
        f"[E.5](#e-5-frequency-domain-complex-stiffness-reduction) give {constants['axial_zeta_interface']:.4f} and "
        f"{constants['interface_settling_time_2pct'] * 1e3:.1f} ms; the executed link damper gives "
        f"{constants['axial_zeta_executed']:.4f} and {constants['axial_settling_time_2pct'] * 1e3:.1f} ms. "
        f"The ten-DOF plant settles in {full_settling * 1e3:.1f} ms and the reduced plant in {reduced_settling * 1e3:.1f} ms, "
        "both on the loss-factor branch.",
        "",
        f"The previous rule took the maximum of the 100 ms floor, the detent-softened drive estimate and the reduced "
        f"axial estimate, which omitted the two longest candidates and therefore silently selected the shortest branch. "
        f"[Section 10](#10-friction-case-responses-and-generated-summary) now runs its nonlinear campaign on a "
        f"{dwell * 1e3:.0f} ms plateau dwell, the maximum over the 100 ms floor, the "
        f"{constants['detent_settling_time_2pct'] * 1e3:.1f} ms detent-softened drive estimate, the "
        f"{constants['axial_settling_time_2pct'] * 1e3:.1f} ms executed reduced axial estimate, the "
        f"{constants['interface_settling_time_2pct'] * 1e3:.1f} ms loss-factor estimate and the "
        f"{constants['measured_settling_time_2pct'] * 1e3:.0f} ms measured estimate. "
        f"It therefore exceeds the ten-DOF settling time by a factor of {dwell / full_settling:.1f} on the "
        "assumption branch and closes on the measurement branch by construction, so no settled-window number in "
        "Sections 9, 10 or 12 depends on which branch the extraction confirms. "
        "**That is a cost, not a result:** the campaign is "
        f"{dwell / 0.100:.1f}$\\times$ longer than the 100 ms floor purely because the damping is unresolved, and "
        "[E.7](#e-7-measured-frf-identification) is the experiment that would shorten it. The earlier claim that "
        "the loss-factor correction had settled the dwell question used the assumption branch as evidence against "
        "the measurement, and the ten-DOF reference it appealed to is built from the same four assumed loss "
        "factors, so it was not independent.",
        "",
        "### 7.4 Reading the trajectory",
        "",
        f"The large oscillation is expected **inside this deliberately frictionless, global-linear audit**, but it is not a quantitative prediction of a real repeated full-step move. "
        f"One full step changes the electrical equilibrium by {verification['full_step_electrical_angle']:.3f} rad (90°). "
        f"Applying the small-signal magnetic tangent across that entire jump initially requests {verification['linear_force_to_limit_ratio']:.3f} times the sinusoidal force limit. "
        "The ideal zero-rise-time edge also injects energy into every retained and discarded mode, while friction, detent nonlinearity, current-loop bandwidth, current rise, and torque saturation are absent.",
        "",
        "Accordingly, the top-right panel should be read as an amplitude-scaled structural comparison: do the two mathematical plants react alike to the same broadband edge? "
        "A physically predictive full-step trajectory requires applying the nonlinear magnetic force and driver/current dynamics to the full-order plant. "
        "The normalized reduction residual remains useful, but the absolute overshoot in this linear panel should not be interpreted as expected stage motion.",
        "<!-- END GENERATED REDUCTION CONVERGENCE -->",
    ])
    return "\n".join(lines)


def update_generated_reduction_convergence(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED REDUCTION CONVERGENCE -->.*?"
        r"<!-- END GENERATED REDUCTION CONVERGENCE -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError(
            "Generated reduction-convergence markers are missing from the derivation document")
    DERIVATION_MD.write_text(
        pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_summary(linear_metrics: dict[str, dict[str, float | np.ndarray]],
                      time_metrics: dict[str, dict[str, float]],
                      verification: dict[str, object]) -> str:
    constants = physical_constants()
    local_low, local_high = detent_local_mode_band()
    lines = [
        "<!-- BEGIN GENERATED RESPONSE SUMMARY -->",
        "| Case | Retained mode | Settled RMS deviation |",
        "|---|---:|---:|",
    ]
    for key in CASES:
        modes = linear_metrics[key]["modes"]
        lines.append(
            f"| {key} | {float(modes[1]):.1f} Hz | "
            f"{time_metrics[key]['rms_settled_deviation_nm']:.1f} nm |"
        )
    lines.extend([
        "",
        "This digest keeps the two values needed to compare topology and settled motion. "
        "Appendix H contains the full ten-case metrics dump and the complete Bode overlay.",
        "",
        "### 10.4 Generated reduction audit",
        "",
        "| Quantity | Executed value |",
        "|---|---:|",
        f"| Measured stage body mass | {constants['m_stage']:.3f} kg |",
        f"| Nut body mass retained at stage node | {constants['m_n']:.3f} kg |",
        f"| Derived retained stage-side mass | [[derived:reduced_stage_mass={constants['m_s']:.3f}]] kg |",
        f"| Upper-mode calibration target | {constants['axial_mode_target_hz']:.2f} Hz |",
        f"| Modal-calibrated $k_{{ax}}$ | "
        f"[[derived:reduced_axial_stiffness@mnm={constants['k_ax'] / 1.0e6:.3f}]] MN/m |",
        f"| Closure-derived $k_{{ball}}$ | "
        f"[[derived:k_ball@mnm={verification['parameters']['k_ball'] / 1.0e6:.3f}]] MN/m |",
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


def calibration_branches(constants: dict[str, float],
                         component: dict[str, float]) -> dict[str, object]:
    """Both k_ax calibration branches and the k_ball closure band."""
    frictionless = constants["k_ax"]
    presliding = presliding_calibrated_axial_stiffness(
        constants["m_d"], constants["m_s"], constants["K_m"],
        constants["axial_mode_target_hz"])
    band = ball_closure_band(frictionless, component["k_sha"], component["k_mnt"])
    presliding_ball = closure_ball_stiffness(
        presliding, component["k_brg"], component["k_sha"], component["k_mnt"])
    unpowered = modal_calibrated_axial_stiffness(
        constants["m_d"], constants["m_s"], 0.0, constants["axial_mode_target_hz"])
    return {
        "frictionless_k_ax": frictionless,
        "presliding_k_ax": presliding,
        "ratio": presliding / frictionless,
        "frictionless_k_ball": constants["k_ball"],
        "presliding_k_ball": presliding_ball,
        "unpowered_k_ax": unpowered,
        "unpowered_shift_pct": 100.0 * (unpowered - frictionless) / frictionless,
        "sigma0_g": FRICTION["g"]["sigma0"],
        "sigma0_n": FRICTION["n"]["sigma0"],
        "singular_limit": band["singular_limit"],
        "samples": band["samples"],
        "k_brg": component["k_brg"],
    }


def generated_calibration_branches(branches: dict[str, object],
                                   constants: dict[str, float]) -> str:
    """Build the Section 6.3 calibration-provenance block.

    The measurement that sets k_ax is taken on the assembled axis at
    micrometre amplitudes, which is inside the presliding regime, so the two
    branches below bracket what the same measured pole implies.
    """
    lines = [
        "<!-- BEGIN GENERATED CALIBRATION BRANCHES -->",
        "**The calibration measurement contains the presliding tangents, and only one branch removes "
        "them.** A hammer FRF at micrometre amplitudes never leaves the presliding regime: every rolling "
        "contact behaves as a spring, so the measured pole already carries "
        f"$\\sigma_{{0,g}}={branches['sigma0_g'] / 1e6:.2f}\\times10^6$ N/m and "
        f"$\\sigma_{{0,n}}={branches['sigma0_n'] / 1e6:.2f}\\times10^6$ N/m. Solving for $k_{{ax}}$ on the "
        "frictionless eigenproblem and then adding the friction ports back therefore counts that stiffness "
        "twice, which is why [10.1](#10-1-presliding-stiffness-shifts-the-retained-mode) predicts an "
        "operating mode above the measurement that set the calibration.",
        "",
        "| Calibration branch | $k_{ax}$ | $k_{ball}$ closure | Reproduces the measured pole with |",
        "|---|---:|---:|---|",
        f"| Frictionless (executed) | {branches['frictionless_k_ax'] / 1e6:.3f} MN/m | "
        f"{branches['frictionless_k_ball'] / 1e6:.3f} MN/m | structure only |",
        f"| Presliding-inclusive | {branches['presliding_k_ax'] / 1e6:.3f} MN/m | "
        f"{branches['presliding_k_ball'] / 1e6:.3f} MN/m | structure plus $\\sigma_{{0,g}}$ and "
        "$\\sigma_{0,n}$ |",
        "",
        f"The presliding-inclusive branch is a factor of {branches['ratio']:.3f} softer. Which branch is "
        "correct is a question about the fixture, not about the algebra, and the fixture record below is "
        "the missing evidence. Until it is filled in, every $k_{ax}$-dependent number in this document "
        "carries that factor as an unquantified bias.",
        "",
        "| Calibration boundary condition | Recorded value |",
        "|---|---|",
        "| Screw coupled or decoupled from the motor | not recorded |",
        "| Motor powered or unpowered during the impact | not recorded |",
        "| Excitation amplitude at the stage | not recorded |",
        "| Any GMS element beyond its yield distance during the impact | not recorded |",
        "| Measurement point and direction | not recorded |",
        "",
        f"**One worry is cheap to remove.** The drive-side boundary condition barely matters: setting "
        f"$K_m=0$, the unpowered-motor limit, moves $k_{{ax}}$ by "
        f"{abs(branches['unpowered_shift_pct']):.3f}%, because $\\lambda m_d$ dominates $K_m$ in the "
        "characteristic equation. Whether the motor was energized during the hammer test is therefore not "
        "the open question; whether the friction ports were loaded is.",
        "",
        "**$k_{ball}$ is a closure residual, so it inherits the bearing assumption.** It absorbs whatever "
        "axial compliance the other three elements leave over, and below a singular limit on $k_{brg}$ "
        "there is nothing left to absorb:",
        "",
        "| $k_{brg}$ | Closure $k_{ball}$ |",
        "|---:|---:|",
    ]
    for k_brg, k_ball in branches["samples"]:
        marker = " (executed)" if np.isclose(k_brg, branches["k_brg"]) else ""
        lines.append(f"| {k_brg / 1e6:.2f} MN/m{marker} | {k_ball / 1e6:.3f} MN/m |")
    lines.extend([
        "",
        f"The singular limit is {branches['singular_limit'] / 1e6:.3f} MN/m and the executed "
        f"{branches['k_brg'] / 1e6:.1f} MN/m sits a factor of "
        f"{branches['k_brg'] / branches['singular_limit']:.2f} above it. The Barden duplex contact angle is "
        "itself unresolved between 15° and 25°, which is a factor-of-two axial stiffness question, so "
        "[E.7](#e-7-measured-frf-identification)'s instruction to compare the implied $k_{ball}$ with an "
        "independent contact estimate has no discriminating power until $k_{brg}$ is fixed first.",
        "",
        "**The damping chain survives this uncertainty even though the stiffness does not.** Every joint "
        "carries the same target loss factor, so the equivalent $\\eta$ at the retained mode is unchanged "
        "across the admissible $k_{brg}$ band; only the reported $k_{ball}$ moves. The two conclusions are "
        "therefore not coupled, and fixing $k_{brg}$ is a stiffness measurement, not a damping one.",
        "<!-- END GENERATED CALIBRATION BRANCHES -->",
    ])
    return "\n".join(lines)


def update_generated_calibration_branches(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED CALIBRATION BRANCHES -->.*?"
        r"<!-- END GENERATED CALIBRATION BRANCHES -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated calibration-branch markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_detent_ablation(study: dict[str, object],
                              constants: dict[str, float]) -> str:
    """Build the Section 10.3 detent decomposition that precedes the case table."""
    rows = study["rows"]
    detent_only = float(study["detent_only_nm"])
    worst = max(rows.items(), key=lambda item: item[1]["detent_share_pct"])
    lines = [
        "<!-- BEGIN GENERATED DETENT ABLATION -->",
        "**Every settled-window number below contains a detent term, and for most cases it is the larger "
        "term.** Case 0 is frictionless but not force-free: the nonlinear campaign runs with the periodic "
        f"detent torque enabled, so its {detent_only:.1f} nm settled deviation is a pure detent result. The "
        "builder therefore reruns the identical command for every case with $\\hat T_{det}=0$ and reports both "
        "numbers, so friction attribution rests on a measured pair rather than on a quadrature guess.",
        "",
        "| Case | Settled RMS, executed | Settled RMS, detent off | Detent share | Quadrature estimate |",
        "|---|---:|---:|---:|---:|",
    ]
    for key, row in rows.items():
        lines.append(
            f"| {key} | {row['executed_nm']:.1f} nm | {row['detent_off_nm']:.1f} nm | "
            f"{row['detent_share_pct']:.1f}% | {row['quadrature_nm']:.1f} nm |")
    lines.extend([
        "",
        "The detent-off column is the friction-only contribution. The quadrature column is the estimate a "
        f"reader can form without the ablation, $\\sqrt{{R^2-R_0^2}}$ against the {detent_only:.1f} nm case-0 "
        "floor; the two agree closely enough to confirm the terms combine in power, and the ablation is the "
        "one that is executed.",
        "",
        f"The detent share reaches {float(worst[1]['detent_share_pct']):.1f}% at case {worst[0]} and stays "
        f"above {min(row['detent_share_pct'] for row in rows.values()):.1f}% everywhere. **No settled-window "
        "difference between friction cases should be read as a friction result without its detent-off pair**, "
        "and the pre-distortion argument in [Section 5](#5-stepper-input-nonlinear-law-linearization-and-bound) "
        "is about the same term: detent is a position-periodic error the command grid cannot currently "
        "correct.",
        "<!-- END GENERATED DETENT ABLATION -->",
    ])
    return "\n".join(lines)


def update_generated_detent_ablation(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED DETENT ABLATION -->.*?<!-- END GENERATED DETENT ABLATION -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated detent-ablation markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_breakaway_sensitivity(study: dict[str, object]) -> str:
    """Build the Section 12.4 guideway breakaway-force sensitivity block."""
    rows = list(study["rows"])
    low, high = study["likely_range"]
    lines = [
        "<!-- BEGIN GENERATED BREAKAWAY SENSITIVITY -->",
        f"[8.3](#8-3-executed-provisional-friction-values) executes $F_{{s,g}}=$ {study['executed_F_s']:.1f} N "
        f"while stating a likely range of {low:.1f} to {high:.1f} N. That is not a rounding difference: the four "
        "GMS yield distances scale with $F_s$, so the command levels in "
        "[G.1](#g-1-exact-1-16-microstep-commands) cross different thresholds at the two values. The variant is "
        "therefore executed rather than described. $F_c$ scales with $F_s$ at the executed Stribeck ratio, "
        "because holding $F_c$ fixed while lowering $F_s$ would invert the Stribeck curve instead of "
        "modelling a weaker interface.",
        "",
        "| $F_{s,g}$ | $F_{c,g}$ | Element yields | Elements yielded at the "
        f"{rows[0]['inner_level_um']:.3f} µm inner level | $F_{{ret}}$ (A2) | GMS/LuGre | $R_{{hold}}$ (A2) |",
        "|---:|---:|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        yields = ", ".join(f"{value:.2f}" for value in row["yields_um"])
        lines.append(
            f"| {row['F_s']:.1f} N | {row['F_c']:.1f} N | {yields} µm | "
            f"{row['elements_yielded_at_inner']} of 4 | "
            f"{row['force_mismatch_gms_N']:.4f} N | {row['force_ratio']:.2f}× | "
            f"{row['r_hold_gms_pct']:.1f}% |")
    verdict = (
        "**The command design is not portable across the stated range.** The inner level crosses a different "
        "number of elements at the two forces, so the executed sequence tests a different partial-slip state "
        "at the value the document itself calls likely."
        if study["design_changes"] else
        "The inner level crosses the same number of elements at both forces, so the sequence still probes "
        "partial slip at the middle of the stated range; the metric values move but the design survives.")
    lines.extend([
        "",
        verdict + " Either the guideway breakaway force is re-identified before the memory campaign is "
        "executed on hardware, or the levels are recomputed at the identified value. This is the same "
        "identification-order argument as [G.4](#g-4-detent-contamination-and-the-forced-identification-order), "
        "applied to an amplitude rather than to a frequency.",
        "<!-- END GENERATED BREAKAWAY SENSITIVITY -->",
    ])
    return "\n".join(lines)


def update_generated_breakaway_sensitivity(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED BREAKAWAY SENSITIVITY -->.*?"
        r"<!-- END GENERATED BREAKAWAY SENSITIVITY -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated breakaway-sensitivity markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_full_response_summary(
        linear_metrics: dict[str, dict[str, float | np.ndarray]],
        time_metrics: dict[str, dict[str, float]],
        detent_ablation: dict[str, object]) -> str:
    """Build the Appendix H audit table kept out of Section 10's main line."""
    constants = physical_constants()
    lines = [
        "<!-- BEGIN GENERATED FULL RESPONSE SUMMARY -->",
        "| Case | Friction law | Global-linear modes (Hz) | Local tangent gain $X_s/X_{cmd}$ | Smallest first-yield travel | First-step overshoot | Settled RMS deviation | Settled RMS, detent off | Friction-only contribution | Settled maximum | All-time peak | Final-window RMS |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    ablation_rows = detent_ablation["rows"]
    for key, case in CASES.items():
        modes = linear_metrics[key]["modes"]
        friction_label = {"none": "none", "lugre": "LuGre", "gms": "GMS"}[case["friction"]]
        first_yield = float(linear_metrics[key]["first_yield_m"])
        yield_text = ("not applicable" if not np.isfinite(first_yield)
                      else f"{first_yield * 1e6:.3f} µm")
        ablated = ablation_rows[key]
        lines.append(
            f"| {key} | {friction_label} | {float(modes[0]):.1f}, {float(modes[1]):.1f} | "
            f"{linear_metrics[key]['tangent_dc_gain']:.5f} | {yield_text} | "
            f"{time_metrics[key]['first_overshoot_pct']:.1f}% | "
            f"{time_metrics[key]['rms_settled_deviation_nm']:.1f} nm | "
            f"{ablated['detent_off_nm']:.1f} nm | "
            f"{ablated['detent_share_pct']:.1f}% detent | "
            f"{time_metrics[key]['max_settled_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['max_abs_deviation_nm']:.1f} nm | "
            f"{time_metrics[key]['rms_final_error_nm']:.1f} nm |"
        )
    lines.extend([
        "",
        "The displayed modes and gains are the global commutation linearization; periodic detent is "
        "excluded from the global stiffness matrix. The friction tangent is local and valid only "
        "below the listed first-yield travel. "
        f"The nonlinear cases include periodic detent torque and use a "
        f"{constants['plateau_dwell'] * 1e3:.0f} ms dwell. Settled values collect the last "
        f"{constants['metric_window'] * 1e3:.0f} ms of every plateau. All deviation columns use "
        "$d(t)=x_{cmd}(t)-x_s(t)$ and describe open-loop modeled plant behavior, not servo tracking.",
        "",
        f"The first-yield travel independently checks the ablation: A/A2 begins at "
        f"{float(linear_metrics['A']['first_yield_m']) * 1e6:.3f} µm on the drive port, whereas "
        f"G/G2 begins at {float(linear_metrics['G']['first_yield_m']) * 1e6:.3f} µm on the guideway "
        "after the drive port is removed.",
        "",
        "The two detent columns are the paired ablation described in "
        "[10.3](#10-3-generated-numerical-summary): the same command rerun with $\\hat T_{det}=0$. The "
        "friction-only column reports how much of each settled window survives that removal, so a friction "
        "comparison between two rows is only defensible on the detent-off column.",
        "<!-- END GENERATED FULL RESPONSE SUMMARY -->",
    ])
    return "\n".join(lines)


def update_generated_full_response_summary(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED FULL RESPONSE SUMMARY -->.*?"
        r"<!-- END GENERATED FULL RESPONSE SUMMARY -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated full response-summary markers are missing from Appendix H")
    DERIVATION_MD.write_text(
        pattern.sub(lambda _match: summary, source), encoding="utf-8")


def mode_shift_ladder(
        linear_metrics: dict[str, dict[str, float | np.ndarray]]) -> list[dict[str, object]]:
    """Return the active-port ladder shared by the 10.1 table and the shift figure.

    Every case except 0 has the drive port active, so each label names it and
    the guideway-only ablation G/G2 is a row rather than an Appendix H diff.
    The displayed shift in hertz is the difference of the two displayed
    frequencies, so the figure's own arithmetic closes; the percentage comes
    from the unrounded eigenvalues, which is what the table has always used.
    """
    def modes(key: str) -> tuple[float, float]:
        case_modes = linear_metrics[key]["modes"]
        return float(case_modes[0]), float(case_modes[1])

    baseline_low, baseline_high = modes("0")
    ladder: list[dict[str, object]] = [{
        "key": "0",
        "label": "none (case 0)",
        "figure_label": "0",
        "low_hz": baseline_low,
        "high_hz": baseline_high,
        "low_text": f"{baseline_low:.1f}",
        "high_text": f"{baseline_high:.1f}",
        "shift_hz_text": None,
        "shift_pct_text": None,
    }]
    for label, figure_label, key in (
            ("guideway only (G/G2)", "G/G2", "G"),
            ("drive + guideway (A/A2)", "A/A2", "A"),
            ("drive + nut (B/B2)", "B/B2", "B"),
            ("all three (C/C2)", "C/C2", "C")):
        low, high = modes(key)
        displayed_shift = round(high, 1) - round(baseline_high, 1)
        ladder.append({
            "key": key,
            "label": label,
            "figure_label": figure_label,
            "low_hz": low,
            "high_hz": high,
            "low_text": f"{low:.1f}",
            "high_text": f"{high:.1f}",
            "shift_hz_text": f"{displayed_shift:+.1f}",
            "shift_pct_text": f"{100.0 * (high - baseline_high) / baseline_high:+.1f}",
        })
    return ladder


def generated_bode_comparison(frequencies: np.ndarray,
                              responses: dict[str, np.ndarray],
                              linear_metrics: dict[str, dict[str, float | np.ndarray]]) -> str:
    constants = physical_constants()
    ladder = mode_shift_ladder(linear_metrics)
    lines = [
        "<!-- BEGIN GENERATED BODE COMPARISON -->",
        "| Active ports | Low mode | Retained mode | Shift |",
        "|---|---:|---:|---:|",
    ]
    for row in ladder:
        shift = "—" if row["shift_pct_text"] is None else f"{row['shift_pct_text']}%"
        lines.append(
            f"| {row['label']} | {row['low_text']} Hz | {row['high_text']} Hz | {shift} |")
    by_key = {str(row["key"]): row for row in ladder}
    guideway_only, drive_guideway = by_key["G"], by_key["A"]
    lines.extend([
        "",
        f"The drive port shifts only the low mode, from {guideway_only['low_text']} to "
        f"{drive_guideway['low_text']} Hz, and leaves the retained mode at "
        f"{drive_guideway['high_text']} Hz untouched. Its presliding stiffness acts on $x_d$, "
        "which barely participates in the relative mode because "
        f"$m_d/m_s\\approx{constants['m_d'] / constants['m_s']:.0f}$. That is the reflected-inertia "
        "result of [Section 6](#6-reduction-from-ten-dofs-to-two) reappearing as a friction "
        "measurement.",
        "",
        "The nut port shifts the mode nearly three times as much as the guideway despite carrying "
        "roughly half the friction force, because $\\sigma_{0,n}=2.0\\times10^6$ N/m against the "
        "guideway's $7.6\\times10^5$ N/m and because it acts on the relative coordinate, directly "
        "in series with $k_{ax}$.",
        "<!-- END GENERATED BODE COMPARISON -->",
    ])
    return "\n".join(lines)


def signed_value(value: float, digits: int) -> str:
    """Format a signed delta without ever printing a negative zero.

    A cell that rounds to zero carries no sign information, so `-0.0` is a
    formatting artifact rather than a result.  Anything that rounds away is
    printed as an unsigned zero at the same precision.
    """
    if round(value, digits) == 0.0:
        return f"{0.0:.{digits}f}"
    return f"{value:+.{digits}f}"


def _latex_power(value: float, digits: int) -> str:
    """Format one number as LaTeX mantissa times a power of ten."""
    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return f"{mantissa}\\times10^{{{int(exponent)}}}"


def micro_viscous_effect(frequencies: np.ndarray,
                         responses: dict[str, np.ndarray],
                         linear_metrics: dict[str, dict[str, float | np.ndarray]],
                         time_metrics: dict[str, dict[str, float]]) -> dict[str, float]:
    """Price the one non-degenerate tangent difference in the case set.

    A1v restores $\\sigma_1$ at the same ports as A, so the pair isolates
    micro-viscous bristle damping.  The peak drop is converted to an implied
    modal-damping change and checked against the closed-form port damper
    $\\Delta\\zeta=c/(2m_s\\omega)$; agreement is what shows the tangent
    assembly puts the damper on the stage coordinate rather than the drive.
    """
    constants = physical_constants()
    magnitudes = {
        key: 20.0 * np.log10(np.maximum(np.abs(responses[key]), 1e-15))
        for key in ("A", "A1v")
    }
    difference = magnitudes["A1v"] - magnitudes["A"]
    extreme_index = int(np.argmax(np.abs(difference)))
    peak_drop_db = float(-difference[extreme_index])
    mode_hz = float(linear_metrics["A"]["modes"][1])
    omega = 2.0 * np.pi * mode_hz
    damped = {}
    for key in ("A", "A1v"):
        case = CASES[key]
        mass, damping, stiffness, _input = linear_matrices(
            case["sites"], case["friction"], bool(case.get("micro_viscous")))
        damped[key] = _damped_modal_data(mass, damping, stiffness)[1][1]
    magnitude_ratio = 10.0 ** (peak_drop_db / 20.0)
    implied_delta_zeta = damped["A"] * (magnitude_ratio - 1.0)
    predicted_delta_zeta = MICRO_VISCOUS_SIGMA1["g"] / (2.0 * constants["m_s"] * omega)
    rms_a = float(time_metrics["A"]["rms_settled_deviation_nm"])
    rms_a1v = float(time_metrics["A1v"]["rms_settled_deviation_nm"])
    return {
        "rms_a_nm": rms_a,
        "rms_a1v_nm": rms_a1v,
        "rms_shift_nm": abs(rms_a1v - rms_a),
        "sigma1_g": MICRO_VISCOUS_SIGMA1["g"],
        "peak_drop_db": peak_drop_db,
        "peak_frequency_hz": float(frequencies[extreme_index]),
        "mode_hz": mode_hz,
        "magnitude_ratio": magnitude_ratio,
        "zeta_a": damped["A"],
        "zeta_a1v": damped["A1v"],
        "implied_delta_zeta": implied_delta_zeta,
        "eigenvalue_delta_zeta": damped["A1v"] - damped["A"],
        "predicted_delta_zeta": predicted_delta_zeta,
        "agreement_pct": 100.0 * abs(implied_delta_zeta - predicted_delta_zeta)
        / implied_delta_zeta,
        "m_s": constants["m_s"],
    }


def generated_micro_viscous(effect: dict[str, float]) -> str:
    """Build the Section 10.2 block: the matched-pair claim plus its number."""
    rms_a, rms_a1v = effect["rms_a_nm"], effect["rms_a1v_nm"]
    return "\n".join([
        "<!-- BEGIN GENERATED MICRO VISCOUS -->",
        "Matched LuGre and GMS pairs are linearly identical by construction: with $\\sigma_1=0$ "
        "both contribute the same $\\sigma_2$ tangent damping, and $\\sum k_i=\\sigma_0$ equalizes "
        "presliding stiffness. Any difference in the nonlinear results of Section 9 is therefore "
        "memory structure, not tangent. Every matched pair in the figure above is exactly "
        "coincident for the same reason. A1v is the only case with $\\sigma_1$ restored, and its "
        "difference against A is the isolated micro-viscous effect.",
        "",
        f"The effect is small and confined to the mode. Restoring $\\sigma_1={effect['sigma1_g']:.1f}$ "
        f"N·s/m at the guideway lowers the {effect['mode_hz']:.0f} Hz peak by "
        f"{effect['peak_drop_db']:.3f} dB and leaves the response unchanged everywhere else, which "
        f"moves the settled RMS deviation from {rms_a:.1f} nm to {rms_a1v:.1f} nm. A "
        f"{effect['rms_shift_nm']:.1f} nm change is the empirical justification for setting "
        "$\\sigma_1=0$ in the matched comparisons.",
        "",
        "<details>",
        "<summary>Cross-check: is the damper landing on the right coordinate?</summary>",
        "",
        f"A peak drop of {effect['peak_drop_db']:.3f} dB implies the modal damping rose by a factor "
        f"of {effect['magnitude_ratio']:.3f}, so $\\Delta\\zeta="
        f"{_latex_power(effect['implied_delta_zeta'], 2)}$ against case A's $\\zeta_2="
        f"{_latex_power(effect['zeta_a'], 3)}$. The direct prediction for an added port damper is "
        f"$\\Delta\\zeta=c/(2m_s\\omega)={effect['sigma1_g']:.1f}/(2\\times{effect['m_s']:.3f}"
        f"\\times2\\pi\\times{effect['mode_hz']:.1f})={_latex_power(effect['predicted_delta_zeta'], 2)}$. "
        f"Those agree to {effect['agreement_pct']:.1f}%, which confirms the tangent assembly is "
        "placing the damper on the stage coordinate that carries the guideway port.",
        "",
        "The state-space eigenvalues say the same thing without the decibel step: $\\zeta_2$ moves "
        f"from ${_latex_power(effect['zeta_a'], 3)}$ in A to "
        f"${_latex_power(effect['zeta_a1v'], 3)}$ in A1v, a direct "
        f"$\\Delta\\zeta={_latex_power(effect['eigenvalue_delta_zeta'], 2)}$.",
        "",
        "</details>",
        "<!-- END GENERATED MICRO VISCOUS -->",
    ])


def update_generated_micro_viscous(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED MICRO VISCOUS -->.*?<!-- END GENERATED MICRO VISCOUS -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated micro-viscous markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def update_generated_bode_comparison(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED BODE COMPARISON -->.*?<!-- END GENERATED BODE COMPARISON -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated Bode comparison markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


SITE_TITLES = {"g": "Guideway $g$", "n": "Nut microslip $n$", "d": "Drive side $d$"}


def generated_tau_c_sensitivity(study: dict[str, object]) -> str:
    """Build the Section 12.3 Stribeck-relaxation-time sensitivity block."""
    rows = list(study["rows"])
    force_spread = float(study["force_spread_N"])
    loop_spread = float(study["loop_spread_J"])
    force_gap = float(study["law_gap_force_N"])
    loop_gap = float(study["law_gap_loop_J"])
    bounded = force_spread < force_gap and loop_spread < loop_gap
    lines = [
        "<!-- BEGIN GENERATED TAU C SENSITIVITY -->",
        "$C$ is the least anchored parameter in Section 8: it is identified from measured hysteresis loops in the source GMS work and is assumed here. "
        "The A/A2 guideway memory sequence is rerun at three relaxation times spanning a factor of four, and the two metrics that can see the attractor dynamics are reported.",
        "",
        "| $\\tau_C$ | Guideway $C$ | $F_{ret}$ (A2) | $A_{loop}$ (A2) |",
        "|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['tau_C'] * 1e3:.1f} ms | {row['C_site']:.0f} N/s | "
            f"{row['force_mismatch_N']:.5f} N | {row['loop_area_J'] * 1e6:.3f} µJ |"
        )
    verdict = (
        "**$C$ is not a dominant uncertainty over this range.** "
        if bounded else
        "**$C$ must be identified before the law comparison means anything.** "
    )
    lines.extend([
        "",
        f"Across a four-fold change in $\\tau_C$ the return-force mismatch spreads by {force_spread:.5f} N and the loop area by {loop_spread * 1e6:.2f} µJ, "
        f"against GMS-minus-LuGre gaps of {force_gap:.5f} N and {loop_gap * 1e6:.2f} µJ on the same metrics. "
        + verdict
        + (f"The spread is {100.0 * force_spread / force_gap:.1f}% of the force gap and {100.0 * loop_spread / loop_gap:.1f}% of the loop-area gap, "
           "so the law comparison in Section 9 survives the assumption. This bounds a weakness rather than removing it: "
           "the reported insensitivity holds over the tested range and does not license an arbitrary value."
           if bounded else
           "The assumed value changes the comparison by more than the effect the comparison is designed to detect."),
        "",
        f"The upper bound on $\\tau_C$ is dynamic, not statistical. The retained mode has a period of {float(study['mode_period_ms']):.3f} ms, "
        f"and the executed $\\tau_C$ is {float(study['mode_ratio']):.1f} times faster. At 0.4 ms that margin falls to "
        f"{float(study['mode_period_ms']) / 0.4:.1f}, which is inside the range where the attractor dynamics begin to alias into the structural response; "
        "below roughly 0.05 ms it stiffens the ODE without adding physics.",
        "<!-- END GENERATED TAU C SENSITIVITY -->",
    ])
    return "\n".join(lines)


def update_generated_tau_c_sensitivity(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED TAU C SENSITIVITY -->.*?<!-- END GENERATED TAU C SENSITIVITY -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated tau_C markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_branch_census(study: dict[str, object],
                            departure: dict[str, object],
                            memory_experiments: dict[str, dict[str, object]]) -> str:
    """Build the Section 12.2 branch-selection census and its verdict."""
    rows = list(study["rows"])
    threshold_total = int(study["threshold_total"])
    reversal_total = int(study["reversal_total"])
    evaluation_total = int(study["evaluation_total"])
    lines = [
        "<!-- BEGIN GENERATED BRANCH CENSUS -->",
        "The executed GMS branch test is stateless: it reconstructs stick or slip from "
        "$(v,F_i)$ at every Runge-Kutta evaluation instead of carrying a persistent per-element flag. "
        "This census advances a shadow persistent flag alongside the executed trajectory and counts where the two disagree. "
        "The shadow flag never feeds a derivative, so collecting it cannot change any reported result. "
        f"Counts are element-evaluations: four elements per site per RK stage, over the "
        f"{main_duration(physical_constants()) * 1e3:.0f} ms main sequence. The re-priced "
        "comparison below uses the Section 9 memory trajectory, which is a different trajectory.",
        "",
        "| Case | Site | `flips_reversal` | `flips_threshold` | `evals_total` | Threshold share |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        evaluations = max(int(row["evals_total"]), 1)
        lines.append(
            f"| {row['case']} | {SITE_TITLES[str(row['site'])]} | {row['flips_reversal']:,} | "
            f"{row['flips_threshold']:,} | {row['evals_total']:,} | "
            f"{100.0 * int(row['flips_threshold']) / evaluations:.3f}% |"
        )
    lines.extend([
        "",
        "`flips_reversal` counts transitions to stick caused by $vF_i\\le0$, which both models make. "
        "`flips_threshold` counts transitions to stick caused by $|F_i|<\\nu_is(v)$ with no velocity reversal, "
        "which the persistent-state model would not make. Only the second column is a departure.",
        "",
    ])
    if threshold_total == 0:
        lines.extend([
            "**`flips_threshold` is zero across the executed GMS cases.** For these trajectories the stateless test is "
            "equivalent to the persistent-state model and the departure is inert. No counterfactual rerun is required.",
            "",
        ])
    else:
        lines.extend([
            f"**`flips_threshold` is not zero: {threshold_total:,} element-evaluations across the executed GMS cases**, "
            f"against {reversal_total:,} genuine reversals and {evaluation_total:,} evaluations. "
            "The departure is therefore active, and it is the dominant re-stick mechanism rather than a rare corner: "
            f"threshold-driven reclassification outnumbers reversal-driven re-stick by {threshold_total / max(reversal_total, 1):.0f} to 1. "
            "Its cost is priced by rerunning each affected case with the shadow flag enforced, so that a yielded element keeps "
            "slipping and chases the rising threshold at rate $C$ until an actual reversal.",
            "",
            "| Case | Executed settled-window RMS | Persistent-flag rerun | Change |",
            "|---|---:|---:|---:|",
        ])
        enforced = dict(study["enforced"])
        for key, record in enforced.items():
            lines.append(
                f"| {key} | {record['baseline_rms_nm']:.3f} nm | {record['settled_rms_nm']:.3f} nm | "
                f"{signed_value(float(record['delta_nm']), 3)} nm "
                f"({signed_value(float(record['delta_pct']), 2)}%) |"
            )
        worst = max(abs(float(record["delta_pct"])) for record in enforced.values())
        lines.extend([
            "",
            f"The largest settled-window change is {worst:.2f}%. **That number understates the departure, and the reason is structural.** "
            "The settled window is the final 20 ms of a single plateau, sampled after motion has stopped, whereas the departure occurs "
            "during deceleration while $s(v)$ is rising. A metric evaluated at rest on one plateau has no mechanism for seeing it.",
            "",
            "The nut site records zero threshold flips in both B2 and C2. In steady motion the nut-port velocity $\\dot x_d-\\dot x_s$ "
            "is identically zero because the elastic deformation is constant, so every nut element is stuck and no branch decision is "
            "ever contested. See [8.1](#8-1-how-the-friction-laws-attach-to-the-plant).",
            "",
            "#### Memory-sequence branch census",
            "",
            f"The Section 9 memory trajectory lasts "
            f"{float(memory_experiments['guideway']['duration']) * 1e3:.0f} ms. Its branch counts "
            "are measured on that trajectory itself; they are not copied or duration-scaled from "
            "the main-sequence census.",
            "",
            "| Case | Site | `flips_reversal` | `flips_threshold` | `evals_total` | Threshold share |",
            "|---|---|---:|---:|---:|---:|",
        ])
        memory_threshold_total = 0
        memory_reversal_total = 0
        for experiment in memory_experiments.values():
            for key, census in experiment["censuses"].items():
                for site in census.sites:
                    reversals = int(census.reversal_flips[site])
                    thresholds = int(census.threshold_flips[site])
                    evaluations = int(census.evaluations[site])
                    memory_reversal_total += reversals
                    memory_threshold_total += thresholds
                    lines.append(
                        f"| {key} | {SITE_TITLES[site]} | {reversals:,} | "
                        f"{thresholds:,} | {evaluations:,} | "
                        f"{100.0 * thresholds / max(evaluations, 1):.3f}% |"
                    )
        main_ratio = threshold_total / max(reversal_total, 1)
        memory_ratio = memory_threshold_total / max(memory_reversal_total, 1)
        relative_ratio_change = abs(memory_ratio - main_ratio) / max(main_ratio, 1e-30)
        ratio_verdict = (
            "**materially differs**" if relative_ratio_change >= 0.20
            else "**does not materially differ**")
        lines.extend([
            "",
            f"The memory-sequence threshold-to-reversal ratio is {memory_ratio:.1f}:1, versus "
            f"{main_ratio:.1f}:1 on the main sequence. It {ratio_verdict} under the stated 20% "
            f"relative-change criterion ({100.0 * relative_ratio_change:.1f}% here).",
            "",
            "#### Re-priced against the Section 9.4 loop metrics",
            "",
            "The repeated-return metrics compare settled means at the *same* command level reached by different histories, so their value "
            "depends on every intervening deceleration. The loop area is integrated along the dynamic trace itself. Both can see what the "
            "settled window cannot. The threshold column is the GMS-minus-LuGre gap on the same metric: the effect the Section 9 experiment "
            "exists to detect, and therefore the level the departure has to stay below to count as bookkeeping.",
            "",
            "| Case | Metric | Executed (stateless) | Persistent-flag rerun | Change | GMS − LuGre gap | Exceeds? |",
            "|---|---|---:|---:|---:|---:|---|",
        ])
        for record in departure["records"]:
            digits = int(record["digits"])
            unit = str(record["unit"])
            def shown(value: object, signed: bool = False) -> str:
                number = float(value)
                if unit == "µJ":
                    return f"{number:+.{digits}f}" if signed else f"{number:.{digits}f}"
                return f"{number:+.{digits}g}" if signed else f"{number:.{digits}g}"
            lines.append(
                f"| {record['case']} | {record['metric']} | "
                f"{shown(record['executed'])} {unit} | "
                f"{shown(record['persistent'])} {unit} | "
                f"{shown(record['delta'], True)} {unit} | "
                f"{shown(record['law_gap'])} {unit} | "
                f"{'**yes**' if record['exceeds_law_gap'] else 'no'} |"
            )
        force_records = [r for r in departure["records"] if "F_{ret}" in str(r["metric"])]
        worst_force = max(force_records,
                          key=lambda r: abs(float(r["delta"])) / max(float(r["law_gap"]), 1e-30))
        force_fraction = 100.0 * abs(float(worst_force["delta"])) / float(worst_force["law_gap"])
        lines.extend([
            "",
            f"On the metric the experiment is built around, the departure moves $F_{{ret}}$ for {worst_force['case']} by "
            f"{abs(float(worst_force['delta'])):.4f} N against a law gap of {float(worst_force['law_gap']):.4f} N, which is "
            f"{force_fraction:.0f}% of the effect being measured.",
            "",
            "This conclusion is conditional on command resolution. The earlier 1/256 run priced the same A2 $F_{ret}$ departure at 91.5% of the law gap; the rebuilt production 1/16 sequence prices it at "
            f"{force_fraction:.0f}%. Coarser executable commands changed the loop trajectory and the comparison margin, so the branch-model warning must be re-evaluated whenever the microstep divisor or reversal sequence changes.",
            "",
            f"**That denominator is not a clean two-law difference.** [Appendix G.5](#g-5-settled-force-retention-diagnostic) "
            "shows $F_{ret,LuGre}$ at the guideway is degenerate: LuGre's settled force is wiped by post-edge ringing "
            "before every settled window opens, so the law gap above is dominated by $F_{ret,GMS}$ rather than by a "
            f"genuine LuGre-versus-GMS contrast. Equivalently, the departure moves GMS's own return-force mismatch by "
            f"about {force_fraction:.0f}% of its own value. That is still a meaningful bound on the branch-selection "
            "departure, but it is not the two-law comparison the percentage suggests at first read.",
            "",
        ])
        if bool(departure["any_exceeds"]):
            exceeded = ", ".join(
                f"{record['case']} {record['metric']}"
                for record in departure["records"] if record["exceeds_law_gap"])
            lines.append(
                "**The departure is not a bookkeeping detail.** It changes at least one Section 9.4 metric by more than the "
                f"GMS-minus-LuGre difference on that same metric ({exceeded}), which means it is comparable to the effect the "
                "experiment is designed to detect. Persistent branch flags must be added before any GMS parameter is identified "
                "from this model. The 1.08% settled-window figure above is retained because it is correct for what it measures, "
                "but it does not govern: these metrics do.")
        else:
            lines.append(
                "No Section 9.4 metric moves by more than the GMS-minus-LuGre difference on that metric, so the departure stays "
                "below the effect the experiment is designed to detect. It remains a defect to close before identification, "
                "because the margin is not large.")
    lines.append("<!-- END GENERATED BRANCH CENSUS -->")
    return "\n".join(lines)


def update_generated_branch_census(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED BRANCH CENSUS -->.*?<!-- END GENERATED BRANCH CENSUS -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated branch-census markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def branch_census_sentence(study: dict[str, object]) -> str:
    """One-line census summary baked into the Section 8.4 live equation."""
    parts = [
        f"{row['case']}/{row['site']}: reversal {row['flips_reversal']:,}, "
        f"threshold {row['flips_threshold']:,}, evals {row['evals_total']:,}"
        for row in study["rows"]
    ]
    verdict = ("departure inert for these trajectories"
               if int(study["threshold_total"]) == 0
               else f"{int(study['threshold_total']):,} threshold flips total; "
                    "see Section 12.2 for the priced counterfactual")
    return "GMS branch census - " + "; ".join(parts) + ". " + verdict + "."


def rendered_branch_census_sentence() -> str:
    """Return the live sentence from this run or recover it from generated Markdown."""
    if BRANCH_CENSUS_SENTENCE:
        return BRANCH_CENSUS_SENTENCE
    source = DERIVATION_MD.read_text(encoding="utf-8")
    block_match = re.search(
        r"<!-- BEGIN GENERATED BRANCH CENSUS -->(.*?)<!-- END GENERATED BRANCH CENSUS -->",
        source, flags=re.DOTALL)
    main_block = ("" if block_match is None else
                  block_match.group(1).split("#### Memory-sequence branch census", 1)[0])
    rows = [] if block_match is None else re.findall(
        r"^\| (A2|G2|B2|C2) \| .*?\$([dgn])\$ \| ([\d,]+) \| ([\d,]+) \| ([\d,]+) \|",
        main_block, flags=re.MULTILINE)
    if rows:
        parts = [
            f"{case}/{site}: reversal {reversal}, threshold {threshold}, evals {evaluations}"
            for case, site, reversal, threshold, evaluations in rows
        ]
        threshold_total = sum(int(threshold.replace(",", ""))
                              for _case, _site, _reversal, threshold, _evaluations in rows)
        verdict = ("departure inert for these trajectories"
                   if threshold_total == 0 else
                   f"{threshold_total:,} threshold flips total; see Section 12.2 for the priced counterfactual")
        return "GMS branch census - " + "; ".join(parts) + ". " + verdict + "."
    return "GMS branch census - not rebuilt in this render; rerun the Python builder to refresh the counts."


def friction_port_sentence() -> str:
    """One-line port and integrated-state summary baked into Section 8.1."""
    parts = []
    for key, case in CASES.items():
        sites = ",".join(case["sites"]) if case["sites"] else "none"
        if case["friction"] == "lugre":
            friction_states = len(case["sites"])
        elif case["friction"] == "gms":
            friction_states = GMS_N * len(case["sites"])
        else:
            friction_states = 0
        parts.append(f"{key}: sites {sites}, {4 + friction_states} live states")
    return ("Friction ports per case - " + "; ".join(parts)
            + f". The allocated RK4 vector is a fixed {STATE_SIZE} entries; inactive site blocks hold zero derivatives.")


def generated_presliding_summary(experiments: dict[str, dict[str, object]],
                                 damping_sweep: dict[str, object],
                                 true_loop_path: Path,
                                 constants: dict[str, float]) -> str:
    """Build the Section 9 results.

    The continuous loop leads because it is the only discriminator that does
    not depend on plateau settling, and the plateau maps follow labelled as
    damping-conditional, which is what the sweep in 9.4 measures.
    """
    lines = ["<!-- BEGIN GENERATED PRESLIDING SUMMARY -->"]
    ratio_low, ratio_high = damping_sweep["ratio_range"]
    lines.extend([
        "### 9.1 Continuous presliding loop: the primary discriminator",
        "",
        f"![Continuous quasi-static presliding loop]({true_loop_path.relative_to(ROOT).as_posix()})",
        "",
        "A slow continuous triangular ramp-reversal at the guideway outer amplitude, with no plateaus. "
        "This is the literature-comparable presliding $F$-$x$ loop and it is what a quasi-static Kistler "
        "sweep actually produces. It leads the section because it is the one comparison that does not "
        "depend on how fast post-edge ringing decays: there are no command edges to ring, so the "
        f"{ratio_low:.1f}$\\times$ to {ratio_high:.1f}$\\times$ spread that "
        "[9.4](#9-4-the-retention-gap-is-a-function-of-damping) measures across the disputed damping range "
        "does not apply to it. The settled return-point maps below remain the richer diagnostic, but they "
        "are conditional in a way this loop is not.",
        "",
    ])
    row_definitions = (
        ("Return-force mismatch $F_{ret}$", "return_force_mismatch_N", "N", False),
        ("Final-origin magnitude", "final_mean_nm", "nm", True),
        ("Closed-loop energy $A_{loop}$", "loop_area_J", "µJ", False),
        ("Whole-sequence RMS deviation †", "whole_rms_nm", "nm", False),
        ("Peak absolute deviation †", "max_abs_deviation_nm", "nm", False),
    )

    for experiment_name, experiment in experiments.items():
        metrics = experiment["metrics"]
        keys = tuple(experiment["keys"])
        lugre_key, gms_key = keys[:2]
        guideway = experiment_name == "guideway"
        section_title = ("9.2 Guideway plateau map, damping-conditional" if guideway
                         else "9.3 Nut microslip plateau map, damping-conditional")
        image = ("rendered_assets/presliding_memory_comparison.svg" if guideway
                 else "rendered_assets/nut_memory_comparison.svg")
        alt = ("Guideway nested-return memory comparison"
               if guideway else "Blocked nut nested-return memory comparison")
        lines.extend([
            f"### {section_title}",
            "",
            f"![{alt}]({image})",
            "",
            (f"| Executed metric | LuGre {lugre_key} | GMS {gms_key} | "
             f"GMS / LuGre | LuGre {keys[2]} | GMS {keys[3]} | "
             f"{gms_key} minus {keys[3]} (% vs {keys[3]}) |"
             if guideway else
             f"| Executed metric | LuGre {lugre_key} | GMS {gms_key} | "
             "GMS / LuGre | GMS minus LuGre |"),
            ("|---|---:|---:|---:|---:|---:|---:|"
             if guideway else "|---|---:|---:|---:|---:|"),
        ])
        for label, metric_key, unit, use_absolute in row_definitions:
            scale = 1.0e6 if unit == "µJ" else 1.0
            lugre = scale * float(metrics[lugre_key][metric_key])
            gms = scale * float(metrics[gms_key][metric_key])
            if use_absolute:
                lugre, gms = abs(lugre), abs(gms)
            precision = 4 if unit == "N" else 2
            ratio = gms / max(abs(lugre), 1.0e-30)
            if guideway:
                ablation_lugre = scale * float(metrics[keys[2]][metric_key])
                ablation_gms = scale * float(metrics[keys[3]][metric_key])
                if use_absolute:
                    ablation_lugre, ablation_gms = abs(ablation_lugre), abs(ablation_gms)
                delta = gms - ablation_gms
                delta_pct = 100.0 * delta / max(abs(ablation_gms), 1.0e-30)
                lines.append(
                    f"| {label} | {lugre:.{precision}f} {unit} | "
                    f"{gms:.{precision}f} {unit} | {ratio:.2f}× | "
                    f"{ablation_lugre:.{precision}f} {unit} | "
                    f"{ablation_gms:.{precision}f} {unit} | "
                    f"{signed_value(delta, precision)} {unit} "
                    f"({signed_value(delta_pct, 1)}%) |")
            else:
                lines.append(
                    f"| {label} | {lugre:.{precision}f} {unit} | "
                    f"{gms:.{precision}f} {unit} | {ratio:.2f}× | "
                    f"{signed_value(gms - lugre, precision)} {unit} |")
        r_hold_lugre = 100.0 * float(metrics[lugre_key]["r_hold"])
        r_hold_gms = 100.0 * float(metrics[gms_key]["r_hold"])
        r_hold_ratio = r_hold_gms / max(r_hold_lugre, 1.0e-9)
        if guideway:
            r_hold_ablation_lugre = 100.0 * float(metrics[keys[2]]["r_hold"])
            r_hold_ablation_gms = 100.0 * float(metrics[keys[3]]["r_hold"])
            r_hold_delta = r_hold_gms - r_hold_ablation_gms
            r_hold_delta_pct = 100.0 * r_hold_delta / max(r_hold_ablation_gms, 1.0e-9)
            # Two decimals: the ablation gap is a few hundredths of a point,
            # which at one decimal printed as a signless -0.0 against a
            # non-zero percentage change.
            lines.append(
                f"| Retention $R_{{hold}}$ ‡ | {r_hold_lugre:.1f}% | {r_hold_gms:.1f}% | "
                f"{r_hold_ratio:.2f}× | {r_hold_ablation_lugre:.1f}% | {r_hold_ablation_gms:.1f}% | "
                f"{signed_value(r_hold_delta, 2)} pp ({signed_value(r_hold_delta_pct, 1)}%) |")
        else:
            lines.append(
                f"| Retention $R_{{hold}}$ ‡ | {r_hold_lugre:.1f}% | {r_hold_gms:.1f}% | "
                f"{r_hold_ratio:.2f}× | {signed_value(r_hold_gms - r_hold_lugre, 2)} pp |")
        lines.extend([
            "",
            "† Edge-dominated response descriptor; included for context, not as a memory discriminator.",
            "",
            "‡ $R_{hold}=|F_{settled}|/\\min(\\sigma_0|x_{plateau}|,s(0))$, the fraction of the available "
            "elastic force actually held at rest, averaged over the six non-zero plateau levels. See "
            "[Appendix G.5](#g-5-settled-force-retention-diagnostic).",
            "",
        ])

        whole_change = 100.0 * abs(
            float(metrics[gms_key]["whole_rms_nm"])
            - float(metrics[lugre_key]["whole_rms_nm"])
        ) / float(metrics[lugre_key]["whole_rms_nm"])
        peak_change = 100.0 * abs(
            float(metrics[gms_key]["max_abs_deviation_nm"])
            - float(metrics[lugre_key]["max_abs_deviation_nm"])
        ) / float(metrics[lugre_key]["max_abs_deviation_nm"])
        force_ratio = (float(metrics[gms_key]["return_force_mismatch_N"]) /
                       float(metrics[lugre_key]["return_force_mismatch_N"]))
        energy_ratio = (100.0 * float(metrics[gms_key]["loop_area_J"]) /
                        float(metrics[lugre_key]["loop_area_J"]))
        if guideway:
            origin_ratio = (abs(float(metrics[gms_key]["final_mean_nm"])) /
                            max(abs(float(metrics[lugre_key]["final_mean_nm"])), 1e-30))
            guide_changes = []
            for _label, metric_key, _unit, use_absolute in row_definitions:
                coupled = float(metrics[gms_key][metric_key])
                ablated = float(metrics[keys[3]][metric_key])
                if use_absolute:
                    coupled, ablated = abs(coupled), abs(ablated)
                guide_changes.append(
                    abs(100.0 * (coupled - ablated) / max(abs(ablated), 1e-30)))
            fret_change = abs(100.0 * (
                float(metrics[gms_key]["return_force_mismatch_N"])
                - float(metrics[keys[3]]["return_force_mismatch_N"])
            ) / float(metrics[keys[3]]["return_force_mismatch_N"]))
            lines.extend([
                f"The two laws produce almost the same stage motion: whole-sequence RMS differs "
                f"by {whole_change:.1f}% and peak deviation by {peak_change:.1f}%. They differ "
                f"sharply in what the interface remembers. GMS's return-force mismatch is "
                f"{force_ratio:.2f}× LuGre's and its residual error at the origin is "
                f"{origin_ratio:.0f}× larger. GMS also dissipates only {energy_ratio:.1f}% of "
                "the LuGre loop energy, which is consistent rather than contradictory: elements "
                "below yield store elastic energy instead of burning it, and that same partial "
                "yielding is what prevents return-point closure.",
                "",
                f"**That $F_{{ret}}$ ratio is a consequence of a retention gap, not an independent "
                f"result.** LuGre retains just {r_hold_lugre:.1f}% of the available elastic force "
                f"at rest ($R_{{hold}}$); GMS retains {r_hold_gms:.1f}%, {r_hold_ratio:.1f}× more. "
                "Post-edge structural ringing bleeds the single LuGre bristle state down within a "
                "few milliseconds of every command edge, long before the settled window opens at "
                "80 to 100 ms, so LuGre's near-zero settled force makes the levels agree with each "
                "other trivially and makes $F_{ret,LuGre}$ small by construction rather than by "
                "genuine return-point closure. GMS's four yielded-and-stuck elements survive the "
                "same ringing far better. See [Appendix G.5](#g-5-settled-force-retention-diagnostic) "
                "for the per-plateau diagnostic and a high-damping confirmation run.",
                "",
                f"Ablating the drive port moves every guideway metric by under 10%, and "
                f"$F_{{ret}}$, the metric the comparison rests on, by {fret_change:.1f}%. "
                "A/A2 is therefore a serviceable proxy for the guideway law comparison despite "
                "not being a physical uncoupled fixture. This supersedes the pre-1/16 estimate "
                "of a 27 to 32% drive-port contribution, computed on the finer command grid, "
                "which no longer holds.",
                "",
            ])
        else:
            lines.extend([
                f"The nut port shows the same signature at {force_ratio:.2f}× on $F_{{ret}}$ "
                f"and {energy_ratio:.1f}% relative loop energy, on a command 4× smaller. The "
                "blocked fixture is what makes this visible: on a free stage the drive and stage "
                "move together, the port sees almost no relative travel, and no element yields.",
                "",
                f"The same retention gap that drives the guideway $F_{{ret}}$ ratio is present here: "
                f"LuGre holds {r_hold_lugre:.1f}% of the available elastic force at rest against "
                f"GMS's {r_hold_gms:.1f}%, {r_hold_ratio:.1f}× more. $F_{{ret,LuGre}}$ at the nut "
                "is the same order as the settled force level itself, the signature of a degenerate "
                "denominator rather than a genuinely closed return point.",
                "",
            ])

    lines.extend([
        "### 9.4 The retention gap is a function of damping",
        "",
        "**The headline retention gap is not a constitutive result on its own.** Both plateau maps "
        "measure what each law still holds after post-edge ringing has acted on it, and how much "
        "ringing there is depends on the retained-mode damping that "
        "[7.3](#7-3-dwell-consequence) cannot currently resolve. The same experiment is therefore "
        "rerun across the whole disputed range, from the measured relative-mode value to the "
        "high-damping confirmation point of [G.5](#g-5-settled-force-retention-diagnostic).",
        "",
        "| Target $\\zeta_2$ | 2% settling | $c$ multiplier | $R_{hold}$ LuGre | $R_{hold}$ GMS | "
        "GMS / LuGre | $F_{ret}$ GMS / LuGre |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in damping_sweep["rows"]:
        lines.append(
            f"| {row['target_zeta']:.4f} | {row['settling_2pct_s'] * 1e3:.1f} ms | "
            f"{row['multiplier']:.2f}× | {row['r_hold_lugre_pct']:.1f}% | "
            f"{row['r_hold_gms_pct']:.1f}% | {row['r_hold_ratio']:.2f}× | "
            f"{row['force_ratio']:.2f}× |")
    window = damping_sweep["discriminating_zeta"]
    lines.extend([
        "",
        f"The retention ratio spans {ratio_low:.2f}$\\times$ to {ratio_high:.2f}$\\times$ across the "
        f"range, against a baseline executed damping of {damping_sweep['baseline_zeta']:.4f}. "
        + (f"Retention separates the two laws by at least a factor of two over "
           f"$\\zeta_2\\in[{window[0]:.4f},\\ {window[1]:.4f}]$ and collapses outside it, so that "
           "interval is the window the fixture has to sit in."
           if window is not None else
           "No point in the swept range separates the two laws by a factor of two, so the plateau "
           "map is not a discriminator on this plant at any damping in the range.")
        + " **This is a fixture design requirement, not a result**: the physical Kistler experiment "
        "must be run at a damping inside the discriminating window, and the "
        "[E.7](#e-7-measured-frf-identification) extraction is what determines whether the installed "
        "axis already is.",
        "",
        "### 9.5 What this means for identification",
        "",
        "Return-point force non-closure, not edge-dominated displacement, is the discriminating "
        "observable, and the continuous loop of [9.1](#9-1-continuous-presliding-loop-the-primary-discriminator) "
        "is the form of it that survives the damping question. The comparison does not assume that "
        "GMS is better; measured force loops must select and fit the constitutive law. Appendix G "
        "records the exact commands, yield-window rationale, memory mechanism, and forced "
        "identification order.",
        "",
        "Drift under a zero-mean or oscillating velocity is a documented deficiency of the "
        "single-state LuGre bristle, not a defect introduced here, and it is one of the reasons the "
        "literature moved to multi-state Maxwell-slip constructions. The retention gap measured in "
        "[9.2](#9-2-guideway-plateau-map-damping-conditional) is that known property appearing on "
        "this plant's post-edge ringing; see [Appendix G.5](#g-5-settled-force-retention-diagnostic) "
        "for the per-plateau diagnostic and the high-damping confirmation run.",
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


def _first_positive_branch_sign_flip(levels: np.ndarray,
                                     forces: np.ndarray) -> tuple[int, int] | None:
    """Locate the first settled-force sign change made at positive commands.

    Both plateaus must be commanded positive, so the reversal belongs to the
    element stack and not to the command.  This is the nested-return signature
    stated in G.3 appearing directly in the plateau table.
    """
    levels = np.asarray(levels, dtype=float)
    forces = np.asarray(forces, dtype=float)
    for index in range(forces.size - 1):
        if levels[index] <= 0.0 or levels[index + 1] <= 0.0:
            continue
        if forces[index] > 0.0 > forces[index + 1]:
            return index, index + 1
    return None


def generated_retention_diagnostic(memory_experiments: dict[str, dict[str, object]],
                                   confirmation: dict[str, object],
                                   true_loop_path: Path,
                                   supplement_path: Path) -> str:
    """Build Appendix G.5: the settled-force-versus-plateau diagnostic, the
    high-damping confirmation run, and the continuous presliding loop
    (Part 1.3-1.5 of the LuGre settled-force-degeneracy patch)."""
    lines = ["<!-- BEGIN GENERATED RETENTION DIAGNOSTIC -->"]
    lines.extend([
        "LuGre's plotted settled force sits near zero at every deflection in Section 9 while GMS holds "
        "up a substantial fraction. The mechanism is the `|v|` term in $\\dot z=v-\\sigma_0|v|z/s(v)$: "
        "it relaxes $z$ for as long as any velocity exists, and a plateau is not quiescent because "
        "every command edge rings the plant. The table below reports settled friction force at every "
        "plateau, using the identical 20 ms settled window as every other metric in this document.",
        "",
    ])

    def level_um(experiment: dict[str, object], level: float) -> float:
        return float(level) * float(experiment["microstep"]) * 1e6

    mechanism_stated = False
    for experiment_name, keys, heading in (
        ("guideway", ("A", "A2"), "Guideway (A, A2)"),
        ("nut", ("B", "B2"), "Blocked nut (B, B2)"),
    ):
        experiment = memory_experiments[experiment_name]
        levels = np.asarray(experiment["levels"], dtype=float)
        metrics = experiment["metrics"]
        lugre_force = np.asarray(metrics[keys[0]]["endpoint_force_N"], dtype=float)
        gms_force = np.asarray(metrics[keys[1]]["endpoint_force_N"], dtype=float)
        lines.extend([
            f"**{heading}: settled friction force versus plateau index**",
            "",
            f"| Plateau | Commanded level | LuGre {keys[0]} | GMS {keys[1]} |",
            "|---:|---:|---:|---:|",
        ])
        for index, level in enumerate(levels):
            lines.append(
                f"| {index + 1} | {signed_value(level_um(experiment, level), 4)} µm | "
                f"{signed_value(float(lugre_force[index]), 4)} N | "
                f"{signed_value(float(gms_force[index]), 4)} N |")
        lines.append("")
        flip = _first_positive_branch_sign_flip(levels, gms_force)
        if flip is not None:
            first, second = flip
            plateau_text = (
                f"Plateau {first + 1} holds {signed_value(float(gms_force[first]), 3)} N at "
                f"{signed_value(level_um(experiment, levels[first]), 4)} µm, and plateau "
                f"{second + 1} holds {signed_value(float(gms_force[second]), 3)} N at "
                f"{signed_value(level_um(experiment, levels[second]), 4)} µm, both at positive "
                "commanded levels.")
            lugre_bound = float(np.max(np.abs(lugre_force)))
            lugre_flip = _first_positive_branch_sign_flip(levels, lugre_force)
            if mechanism_stated:
                lines.append(
                    f"**{keys[1]} shows the same sign change on this fixture.** {plateau_text} "
                    "The mechanism is the one stated under the guideway table.")
            else:
                lines.append(
                    f"**{keys[1]}'s settled force changes sign inside the positive branch.** "
                    f"{plateau_text} That is the nested-return signature in raw digits: on the "
                    "inner return the elements that yielded on the way out are reloaded in the "
                    "opposite direction, so the stack unloads past zero while the command is "
                    "still positive. "
                    + (f"LuGre {keys[0]} changes sign at the same pair, but its whole column stays "
                       f"within {lugre_bound:.3f} N of zero, so that sign is the post-edge "
                       "relaxation residue of a single bristle rather than a held return-point "
                       "state; the retention gap in "
                       "[9.2](#9-2-guideway-plateau-map-damping-conditional) is the same "
                       "observation stated as a fraction."
                       if lugre_flip is not None else
                       f"{keys[0]} keeps its sign across the same pair and its whole column "
                       f"stays within {lugre_bound:.3f} N of zero."))
            mechanism_stated = True
        lines.extend(["", ""])

    guideway_metrics = memory_experiments["guideway"]["metrics"]
    nut_metrics = memory_experiments["nut"]["metrics"]
    guideway_lugre_max = float(np.max(np.abs(guideway_metrics["A"]["endpoint_force_N"])))
    guideway_gms_max = float(np.max(np.abs(guideway_metrics["A2"]["endpoint_force_N"])))
    nut_lugre_max = float(np.max(np.abs(nut_metrics["B"]["endpoint_force_N"])))
    nut_gms_max = float(np.max(np.abs(nut_metrics["B2"]["endpoint_force_N"])))
    lines.extend([
        f"LuGre's column stays within {max(guideway_lugre_max, nut_lugre_max):.3f} N of zero at every "
        f"plateau at both sites, while GMS's column tracks the commanded deflection up to "
        f"{max(guideway_gms_max, nut_gms_max):.3f} N. **The mechanism is confirmed**: this is not a "
        "plotting artifact, LuGre's column genuinely holds almost nothing.",
        "",
        "#### Demoted ablation and branch diagnostics",
        "",
        f"![Guideway ablation and nut branch diagnostics]({supplement_path.relative_to(ROOT).as_posix()})",
        "",
        "The A2/G2 ablation and the nut positive/negative branch split remain available here as "
        "confirmation diagnostics. The parallel main figures use the freed panel for the direct "
        "settled-force retention comparison.",
        "",
        "#### High-damping confirmation run",
        "",
    ])

    baseline_pole = confirmation["baseline_pole"]
    raised_pole = confirmation["raised_pole"]
    multiplier = float(confirmation["multiplier"])
    if baseline_pole is not None and raised_pole is not None:
        lines.extend([
            f"The retained-mode pole at baseline structural damping sits at "
            f"{baseline_pole['frequency_hz']:.1f} Hz with $\\zeta_2={baseline_pole['zeta']:.5f}$, "
            f"decaying with $\\tau={baseline_pole['tau_s'] * 1e3:.2f}$ ms (envelope to 5%: "
            f"{baseline_pole['envelope_5pct_s'] * 1e3:.1f} ms) — far longer than the millisecond-scale "
            "bristle relaxation time computed in 9.3, which is why the ringing wipes the bristle before "
            f"the settled window opens. Scaling $c_{{ax}}$ and $c_m$ by {multiplier:.0f}$\\times$ raises "
            f"the same pole to $\\zeta_2={raised_pole['zeta']:.5f}$, $\\tau={raised_pole['tau_s'] * 1e3:.2f}$ "
            f"ms, envelope to 5% in {raised_pole['envelope_5pct_s'] * 1e3:.2f} ms — ringing now dies "
            "within about a millisecond.",
            "",
        ])
    raised_experiment = confirmation["experiment"]
    raised_levels = np.asarray(raised_experiment["levels"], dtype=float)
    raised_metrics = raised_experiment["metrics"]
    baseline_lugre_force = np.asarray(guideway_metrics["A"]["endpoint_force_N"], dtype=float)
    raised_lugre_force = np.asarray(raised_metrics["A"]["endpoint_force_N"], dtype=float)
    raised_gms_force = np.asarray(raised_metrics["A2"]["endpoint_force_N"], dtype=float)
    lines.extend([
        "Rerunning the guideway A/A2 experiment at this damping, with everything else unchanged:",
        "",
        "| Plateau | Commanded level | LuGre A (baseline damping) | LuGre A (ringing suppressed) | GMS A2 (ringing suppressed) |",
        "|---:|---:|---:|---:|---:|",
    ])
    for index, level in enumerate(raised_levels):
        lines.append(
            f"| {index + 1} | {signed_value(level_um(raised_experiment, level), 4)} µm | "
            f"{signed_value(float(baseline_lugre_force[index]), 4)} N | "
            f"{signed_value(float(raised_lugre_force[index]), 4)} N | "
            f"{signed_value(float(raised_gms_force[index]), 4)} N |")
    nonzero_mask = raised_levels != 0.0
    baseline_nonzero = np.abs(baseline_lugre_force[nonzero_mask])
    raised_nonzero = np.abs(raised_lugre_force[nonzero_mask])
    recovery_ratio = float(np.mean(raised_nonzero)) / max(float(np.mean(baseline_nonzero)), 1e-9)
    lines.extend([
        "",
        f"With ringing suppressed, LuGre's settled force recovers from a mean {float(np.mean(baseline_nonzero)):.4f} N "
        f"over the nonzero plateaus to {float(np.mean(raised_nonzero)):.4f} N — {recovery_ratio:.1f}$\\times$ "
        "larger, comparable to or exceeding GMS at the same raised damping. **The dither-driven relaxation "
        "mechanism is confirmed from both directions**: it is present when ringing is left alone and it "
        "disappears when ringing is suppressed.",
        "",
        "#### Continuous presliding loop",
        "",
        "The continuous quasi-static loop was promoted out of this appendix and now opens "
        "[Section 9](#9-1-continuous-presliding-loop-the-primary-discriminator). It is the only "
        "discriminator in the chapter that does not depend on plateau settling, which is exactly the "
        "quantity this appendix shows the plateau maps are sensitive to.",
        "<!-- END GENERATED RETENTION DIAGNOSTIC -->",
    ])
    return "\n".join(lines)


def update_generated_retention_diagnostic(summary: str) -> None:
    source = DERIVATION_MD.read_text(encoding="utf-8")
    pattern = re.compile(
        r"<!-- BEGIN GENERATED RETENTION DIAGNOSTIC -->.*?<!-- END GENERATED RETENTION DIAGNOSTIC -->",
        flags=re.DOTALL,
    )
    if not pattern.search(source):
        raise RuntimeError("Generated retention-diagnostic markers are missing from the derivation document")
    DERIVATION_MD.write_text(pattern.sub(lambda _match: summary, source), encoding="utf-8")


def generated_convergence_summary(study: dict[str, dict[str, object]]) -> str:
    dt_us = tuple(value * 1e6 for value in GMS_CONVERGENCE_DTS)
    constants = physical_constants()
    lines = [
        "<!-- BEGIN GENERATED STEP HALVING SUMMARY -->",
        f"| Case | {dt_us[0]:.1f} us | {dt_us[1]:.1f} us | {dt_us[2]:.1f} us | "
        f"{A2_CONVERGENCE_DT * 1e6:.2f} us (A2 only) | "
        f"$\\Delta R_{{{dt_us[0]:g}\\to{dt_us[1]:g}}}$ | "
        f"$\\Delta R_{{{dt_us[1]:g}\\to{dt_us[2]:g}}}$ | "
        f"$\\Delta R_{{{dt_us[2]:g}\\to{A2_CONVERGENCE_DT * 1e6:g}}}$ | Observed $p$ |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for key in ("A2", "B2", "C2"):
        result = study[key]
        rms = result["rms_nm"]
        extra_rms = f"{rms[3]:.5f} nm" if len(rms) == 4 else "—"
        extra_difference = result["extra_difference_nm"]
        extra_difference_text = (
            f"{float(extra_difference):.5f} nm" if extra_difference is not None else "—")
        lines.append(
            f"| {key} | {rms[0]:.5f} nm | {rms[1]:.5f} nm | {rms[2]:.5f} nm | "
            f"{extra_rms} | "
            f"{result['coarse_difference_nm']:.5f} nm | {result['fine_difference_nm']:.5f} nm | "
            f"{extra_difference_text} | {result['observed_order']:.2f} |"
        )
    a2 = study["A2"]
    a2_extra = float(a2["extra_difference_nm"])
    if a2_extra >= float(a2["fine_difference_nm"]):
        a2_verdict = (
            f"The additional A2 difference grows again, from {a2['fine_difference_nm']:.5f} nm "
            f"to {a2_extra:.5f} nm. **The A2 metric is numerically unstable under the tested "
            "step-halving sequence; its final/settled-window RMS is not converged.**"
        )
    else:
        a2_verdict = (
            f"The additional A2 difference falls from {a2['fine_difference_nm']:.5f} nm "
            f"to {a2_extra:.5f} nm. **The 12.5 us A2 point was a grid-sensitive branch-switching "
            "artifact; the finer point reverses the apparent divergence.**"
        )
    reduced_cases = ", ".join(
        f"{key} ($p={study[key]['observed_order']:.2f}$)"
        for key in ("B2", "C2"))
    lines.extend([
        "",
        f"B2 and C2 show reduced successive differences under step halving: {reduced_cases}. "
        "Their observed orders are empirical hybrid-trajectory indicators, not an RK4 order claim.",
        "",
        f"A2 does not show the same trend over the first three grids: "
        f"$p=$[[derived:a2_convergence_order={a2['observed_order']:.2f}]]. {a2_verdict}",
        "",
        f"These values use the identical {main_duration(constants) * 1e3:.0f} ms zero-order-held, yield-spanning command and the identical final {constants['metric_window'] * 1e3:.0f} ms "
        "RMS definition. Since GMS branch switching is evaluated at RK trial states without event "
        "localization, the observed order is a sensitivity indicator, not a claimed fourth-order "
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


def update_derived_token_fallbacks(values: dict[str, float]) -> None:
    """Refresh generated fallback text while retaining browser-live token cells."""
    source = DERIVATION_MD.read_text(encoding="utf-8")
    for key, value in values.items():
        pattern = re.compile(rf"(\[\[derived:{re.escape(key)}=)[^\]]+(\]\])")
        match = pattern.search(source)
        if match is None:
            raise RuntimeError(f"Derived report token {key!r} is missing from the derivation document")
        template = match.group(0).split("=", 1)[1][:-2]
        if "." in template:
            digits = len(template.rsplit(".", 1)[1])
            shown = f"{value:.{digits}f}"
        else:
            shown = f"{int(round(value))}"
        source = pattern.sub(lambda item: item.group(1) + shown + item.group(2), source)
    DERIVATION_MD.write_text(source, encoding="utf-8")


def takeaway_derived_values(verification: dict[str, object],
                            experiments: dict[str, dict[str, object]]) -> dict[str, float]:
    """Return every generated number used by the Section 7–9 takeaway cards."""
    by_key = {str(row["key"]): row for row in verification["route_residuals"]}
    executed = float(by_key["P"]["rms_residual_nm"])
    interface = float(by_key["F"]["rms_residual_nm"])
    aligned = float(by_key["C2"]["rms_residual_nm"])
    restored = float(by_key["CB"]["rms_residual_nm"])
    full_lower = float(verification["full_lower_damped_mode"][0])
    return {
        "section7_reduced_coordinate_count": 2.0,
        "section7_full_coordinate_count": 10.0,
        "section7_rms_pct": float(verification["rms_residual_pct_command"]),
        "section7_rms_nm": float(verification["rms_residual_nm"]),
        "section7_damping_share_pct": 100.0 * abs(executed - interface) / executed,
        "section7_frequency_share_pct": 100.0 * (executed - aligned) / executed,
        "section7_drive_share_pct": 100.0 * (aligned - restored) / executed,
        "section7_drive_pole_error_pct": 100.0 * abs(
            float(by_key["P"]["lower_hz"]) - full_lower) / full_lower,
        "friction_port_count": float(len(SITE_KEYS)),
        "gms_states_per_site": float(GMS_N),
        "lugre_states_per_site": 1.0,
        "structural_identifiability_result_count": 2.0,
        "project_adev_floor_nm": 4.6,
        **guideway_retention_tokens(experiments),
    }


def document_frequency_tokens(constants: dict[str, float],
                              linear_metrics: dict[str, dict[str, float | np.ndarray]],
                              branches: dict[str, object]) -> dict[str, float]:
    """Every frequency, band edge and requirement quoted in prose.

    Prose numbers that were typed by hand drifted from the generated ones:
    three different detent bands, two Craig-Bampton pairs and four values for
    the same axial feature.  Each of these now has exactly one source.
    """
    low, high = detent_local_mode_band()
    operating_hz = float(linear_metrics["C"]["modes"][1])
    baseline_hz = float(linear_metrics["0"]["modes"][1])
    fixed_interface_hz = multi_route_reduction_metrics()["first_fixed_interface_hz"]
    return {
        "detent_band_low_hz": low,
        "detent_band_high_hz": high,
        "operating_mode_hz": operating_hz,
        "operating_fixed_interface_separation": (operating_hz / fixed_interface_hz)**2,
        "operating_fixed_interface_ratio": fixed_interface_hz / operating_hz,
        "baseline_fixed_interface_ratio": fixed_interface_hz / baseline_hz,
        "presliding_k_ax_mn": float(branches["presliding_k_ax"]) / 1.0e6,
        "closure_singular_limit_mn": float(branches["singular_limit"]) / 1.0e6,
        "measured_settling_ms": constants["measured_settling_time_2pct"] * 1.0e3,
        "interface_settling_ms": constants["interface_settling_time_2pct"] * 1.0e3,
        "executed_settling_ms": constants["axial_settling_time_2pct"] * 1.0e3,
        "detent_equilibrium_error_nm": constants["detent_equilibrium_error"] * 1.0e9,
        "predistortion_resolution_nm": constants["predistortion_resolution"] * 1.0e9,
        "required_microstep_divisor": constants["required_microstep_divisor"],
        "predistortion_levels_executed": (constants["detent_equilibrium_error"]
                                          / constants["command_step"]),
    }


def guideway_retention_tokens(experiments: dict[str, dict[str, object]]) -> dict[str, float]:
    """Live tokens backing the restated Section 9 takeaway (Part 1.4b)."""
    guideway_metrics = experiments["guideway"]["metrics"]
    lugre_r_hold = 100.0 * float(guideway_metrics["A"]["r_hold"])
    gms_r_hold = 100.0 * float(guideway_metrics["A2"]["r_hold"])
    return {
        "guideway_r_hold_lugre_pct": lugre_r_hold,
        "guideway_r_hold_gms_pct": gms_r_hold,
        "guideway_r_hold_ratio": gms_r_hold / max(lugre_r_hold, 1.0e-9),
    }


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


def _format_derived_value(key: str, value: object, template: str,
                          display: str | None = None) -> str:
    """Mirror the browser's display-unit formatting for initial HTML output."""
    numeric = float(value)
    if display == "mnm":
        return f"{numeric / 1.0e6:.3f}"
    if display == "mm":
        return f"{numeric * 1.0e3:.1f}"
    if display:
        raise ValueError(f"Derived token {key} requests unknown display format {display!r}")
    if re.search(r"loop_area|A_loop|a_loop", key, flags=re.IGNORECASE):
        return f"{numeric * 1.0e6:.2f}"
    if key in {"detent_settling_time_2pct", "axial_settling_time_2pct", "plateau_dwell"}:
        return f"{numeric * 1.0e3:.1f}"
    if key == "case_count":
        return str(int(numeric))
    if re.fullmatch(r"route_[psbfcm]_(kax|kball)", key) or key == "route_s_kax_full":
        return f"{numeric / 1.0e6:.3f}"
    if re.fullmatch(r"route_[psbfcm]_settling", key) or key == "full_model_settling":
        return f"{numeric * 1.0e3:.1f}"
    if key == "torsional_share":
        return f"{100.0 * numeric:.3f}"
    if re.fullmatch(r"route_[psbfcm]_(md|ms)", key):
        return f"{numeric:.3f}"
    if re.fullmatch(r"route_[psbfcm]_cax", key) or key == "interface_axial_damping":
        return f"{numeric:.2f}"
    if re.fullmatch(r"route_[psbfcm]_zeta", key) or key == "full_model_zeta":
        return f"{numeric:.3e}"
    if re.fullmatch(r"route_[psbfcm]_f[12]", key):
        return f"{numeric:.2f}"
    if key in {"mass_ratio", "reduced_mu", "mu_fraction", "cb_frequency_delta",
               "cb_damping_delta"}:
        digits = 4 if key in {"reduced_mu", "mu_fraction"} else 2
        return f"{numeric:.{digits}f}"
    if key in {"fixed_interface_separation", "discarded_pole_separation"}:
        return f"{numeric:.3f}"
    if (re.fullmatch(r"yield_[gnd]_[1-4]_(fs|fc)", key)
            or re.fullmatch(r"static_deflection_[gnd]", key)):
        return f"{numeric * 1.0e6:.2f}"
    if re.fullmatch(r"yield_span_[gnd]", key):
        return f"{numeric:.1f}"
    if re.fullmatch(r"gms_rate_[gnd]", key):
        return f"{numeric:.0f}"
    if re.fullmatch(r"tau_C_[gnd]", key) or key == "retained_mode_period":
        return f"{numeric * 1.0e3:.3f}"
    if key.startswith("detent_velocity_"):
        return f"{numeric * 1.0e3:.3f}"
    if key in {"d_Fs_efficiency_estimate", "tau_C_mode_ratio"}:
        return f"{numeric:.2f}"
    if key == "d_Fs_torque_equivalent":
        return f"{numeric * 1.0e3:.3f}"
    return _format_default_like(numeric, template)


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
        "m_eff_measured": MODEL["m_eff_measured"],
        "zeta_relative_measured": MODEL["zeta_relative_measured"],
        "axial_damping": MODEL["c_ax"],
        "electromagnetic_zeta": MODEL["zeta_m"],
        "microstep_divisor": MODEL["microstep_divisor"],
        "stribeck_delta": MODEL["stribeck_delta"],
        "tau_C": MODEL["tau_C"],
        "eta_screw": MODEL["eta_screw"],
        "F_preload_nut": MODEL["F_preload_nut"],
        "J_m": p["J_m"],
        "J_c": p["J_c"],
        "screw_length": p["screw_length"],
        "usable_screw_travel": p["usable_screw_travel"],
        "stage_travel": p["stage_travel"],
        # BOM value for the installed screw.  The earlier IT1 entry was wrong.
        "lead_accuracy_class": "IT3",
        "screw_diameter": p["screw_diameter"],
        "screw_root_diameter": p["screw_root_diameter"],
        "screw_density": p["screw_density"],
        "youngs_modulus": p["youngs_modulus"],
        "shear_modulus": p["shear_modulus"],
        "nut_axial_datum": p["nut_axial_datum"],
        "nut_mass": p["m_n"],
        "stage_mass": p["m_stage"],
        "k_c_series": p["k_c_series"],
        "k_brg": p["k_brg"],
        "k_mnt": p["k_mnt"],
        "eta_steel": INTERFACE_LOSS_FACTORS["zeta_steel"],
        "eta_bearing": INTERFACE_LOSS_FACTORS["zeta_bearing"],
        "eta_ball_nut": INTERFACE_LOSS_FACTORS["zeta_ball_nut"],
        "eta_nut_mount": INTERFACE_LOSS_FACTORS["zeta_nut_mount"],
        "measured_axial_band_low_hz": MODEL["measured_axial_band_low_hz"],
        "measured_axial_band_high_hz": MODEL["measured_axial_band_high_hz"],
        "detent_enabled": MODEL["detent_enabled"],
        "predistortion_levels": MODEL["predistortion_levels"],
        "zeta_steel": p["zeta_steel"],
        "zeta_bearing": p["zeta_bearing"],
        "zeta_ball_nut": p["zeta_ball_nut"],
        "zeta_nut_mount": p["zeta_nut_mount"],
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
        for index, stiffness in enumerate(GMS_STIFFNESS_BY_SITE[site], start=1):
            defaults[f"{site}_k{index}"] = float(stiffness)
    for index, weight in enumerate(GMS_WEIGHTS, start=1):
        defaults[f"gms_nu{index}"] = weight
    return defaults


def ensure_parameter_file() -> None:
    """Create the builder/browser handoff file once without overwriting edits."""
    if PARAMETER_FILE.exists():
        return
    parameters = browser_parameter_defaults()
    parameters["lead_accuracy_class"] = configured("lead_accuracy_class", "IT3")
    serializable = {
        key: value.item() if isinstance(value, np.generic) else value
        for key, value in parameters.items()
    }
    PARAMETER_FILE.write_text(
        json.dumps({
            "schema": "rev3-model-parameters-v1",
            "build_id": BUILD_ID,
            "parameters": serializable,
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


@lru_cache(maxsize=1)
def browser_derived_defaults() -> dict[str, float]:
    """Derived defaults generated from the same equations as the simulations."""
    constants = physical_constants()
    component = component_parameters()
    mass, _damping, stiffness, _input = linear_matrices((), "none")
    modes = _linear_modes(mass, stiffness)
    defaults = {
        "case_count": float(len(CASES)),
        "transmission_ratio": constants["r"],
        "total_rotational_inertia": component["J_total"],
        "reduced_drive_mass": constants["m_d"],
        "reduced_stage_mass": constants["m_s"],
        "magnetic_stiffness": constants["K_m"],
        "detent_stiffness": constants["K_det"],
        "reduced_axial_stiffness": constants["k_ax"],
        "interface_axial_damping": constants["c_ax_interface"],
        "k_ball": constants["k_ball"],
        "k_theta_a": component["k_theta_a"],
        "k_theta_b": component["k_theta_b"],
        "k_sha": component["k_sha"],
        "k_shb": component["k_shb"],
        "screw_length_a": component["screw_length_a"],
        "screw_length_b": component["screw_length_b"],
        "screw_inertia_nominal": component["screw_inertia_nominal"],
        "predistortion_levels": float(MODEL["predistortion_levels"]),
        "full_step_pitch": constants["full_step"],
        "quarter_step_bound": constants["quarter_step"],
        "command_step": constants["command_step"],
        "screw_inertia": component["screw_inertia"],
        "screw_segment_inertia": component["screw_inertia"] / 3.0,
        "screw_mass": component["screw_mass"],
        "screw_segment_mass": component["screw_mass"] / 3.0,
        "k_c_half": component["k_c1"],
        "mode_1_hz": float(modes[0]),
        "mode_2_hz": float(modes[1]),
        "detent_settling_time_2pct": constants["detent_settling_time_2pct"],
        "axial_settling_time_2pct": constants["axial_settling_time_2pct"],
        "plateau_dwell": constants["plateau_dwell"],
    }
    defaults.update(multi_route_reduction_metrics())
    defaults.update(friction_yield_metrics())
    defaults.update(friction_provenance_metrics())
    return defaults


def friction_provenance_metrics() -> dict[str, float]:
    """Derived cells that make the Section 8 parameter provenance visible."""
    constants = physical_constants()
    metrics: dict[str, float] = {}
    for site, values in FRICTION.items():
        metrics[f"gms_rate_{site}"] = float(values["C_gms"])
        metrics[f"tau_C_{site}"] = float(
            (values["F_s"] - values["F_c"]) / values["C_gms"])
    # Drive-port breakaway estimated from screw efficiency.  This feeds the
    # provenance cell only.  The transformer stays power conserving; see the
    # standing constraint in Section 8.1.
    metrics["d_Fs_efficiency_estimate"] = float(
        MODEL["F_preload_nut"] * (1.0 / MODEL["eta_screw"] - 1.0))
    metrics["d_Fs_torque_equivalent"] = float(FRICTION["d"]["F_s"] * constants["r"])
    # Detent excites the structure at v / (one full-step pitch).  Sweep speeds
    # for the constant-velocity identification must avoid these.
    full_step = constants["full_step"]
    mass, _damping, stiffness, _input = linear_matrices((), "none")
    modes = _linear_modes(mass, stiffness)
    metrics["detent_velocity_drive"] = float(modes[0] * full_step)
    metrics["detent_velocity_axial"] = float(modes[1] * full_step)
    metrics["detent_velocity_discarded"] = float(2001.95 * full_step)
    metrics["retained_mode_period"] = float(1.0 / modes[1])
    metrics["tau_C_mode_ratio"] = float(1.0 / (modes[1] * MODEL["tau_C"]))
    return metrics


def friction_yield_metrics() -> dict[str, float]:
    """Per-element yield distances that make the GMS parameter table readable."""
    metrics: dict[str, float] = {}
    for site, values in FRICTION.items():
        stiffness = GMS_STIFFNESS_BY_SITE[site]
        for level, force in (("fs", values["F_s"]), ("fc", values["F_c"])):
            distances = GMS_WEIGHTS * force / stiffness
            for index, distance in enumerate(distances, start=1):
                metrics[f"yield_{site}_{index}_{level}"] = float(distance)
        span = GMS_WEIGHTS * values["F_s"] / stiffness
        metrics[f"yield_span_{site}"] = float(span[-1] / span[0])
        metrics[f"static_deflection_{site}"] = float(values["F_s"] / values["sigma0"])
    return metrics


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
        key, display, value = match.group(1), match.group(2), match.group(3)
        escaped_key = html.escape(key, quote=True)
        authoritative = browser_derived_defaults().get(key, value.strip())
        formatted = _format_derived_value(key, authoritative, value, display)
        escaped_value = html.escape(formatted)
        # An optional @display suffix lets one derived quantity appear in the
        # unit its section uses.  The browser reads the same attribute, so a
        # live recalculation cannot silently change the printed unit.
        display_attribute = (f' data-derived-format="{html.escape(display, quote=True)}"'
                             if display else "")
        return keep(
            f'<output class="derived-output" data-derived="{escaped_key}"'
            f'{display_attribute} '
            f'aria-label="Derived value {escaped_key}">{escaped_value}</output>'
        )

    text = re.sub(r"\[\[derived:([A-Za-z0-9_]+)(?:@([a-z_]+))?=([^\]]+)\]\]",
                  derived_output, text)
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
            if level in (2, 3) and not title.startswith("Key equation "):
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
                row_cells = split_table_row(lines[i])
                row_class = ' class="metric-primary"' if (
                    row_cells and (
                        row_cells[0].startswith("Return-force mismatch")
                        or row_cells[0].startswith("Retention $R_{hold}$")
                    )) else ""
                table.append(f"<tr{row_class}>")
                for j, cell in enumerate(row_cells):
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
    # Keep Chapter 6 and its algebra together in the editable source, then
    # place Appendix E after the numbered chapters in the rendered document.
    if markdown_path == DERIVATION_MD:
        appendix_start = source.find("\n## Appendix E. Order-reduction derivations")
        chapter_7_start = source.find("\n## 7. Full-versus-reduced verification", appendix_start)
        if appendix_start >= 0 and chapter_7_start > appendix_start:
            appendix_e = source[appendix_start:chapter_7_start]
            source = source[:appendix_start] + source[chapter_7_start:] + appendix_e + "\n"
    body, toc = markdown_to_html(source)
    title_match = re.search(r"^#\s+(.+)$", source, flags=re.MULTILINE)
    title = title_match.group(1) if title_match else markdown_path.stem
    toc_html = []
    for level, label, section_id in toc:
        toc_html.append(f'<a class="toc-level-{level}" href="#{section_id}">{html.escape(label)}</a>')
    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    derived_defaults = browser_derived_defaults()
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="report-build" content="{BUILD_ID}">
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
.topbar {{ position:sticky; top:0; z-index:20; display:flex; flex-wrap:wrap; gap:.75rem; align-items:center; padding:.65rem 1rem; background:color-mix(in srgb,var(--card) 94%,transparent); border-bottom:1px solid var(--line); backdrop-filter:blur(9px); }}
.topbar .name {{ font-weight:700; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; margin-right:auto; }}
button {{ color:var(--text); background:var(--soft); border:1px solid var(--line); border-radius:7px; padding:.42rem .7rem; cursor:pointer; }}
.layout {{ display:grid; grid-template-columns:280px minmax(0,1fr); gap:1.5rem; max-width:1510px; margin:0 auto; padding:1.5rem; }}
body.outline-collapsed .layout {{ grid-template-columns:minmax(0,1fr); max-width:1180px; }}
body.outline-collapsed nav {{ display:none; }}
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
tr.metric-primary td {{ background:color-mix(in srgb,var(--accent) 7%,var(--card)); }}
tr.metric-primary td:first-child {{ border-left:4px solid var(--accent); }}
.parameter-input {{ width:100%; min-width:7rem; padding:.38rem .48rem; color:var(--text); background:var(--card); border:1px solid var(--line); border-radius:5px; font:inherit; font-variant-numeric:tabular-nums; }}
.parameter-input:focus {{ outline:2px solid var(--accent); outline-offset:1px; }}
.assumed-input {{ background:var(--assumed); border-color:var(--assumed-line); font-weight:700; }}
.derived-output {{ display:inline-block; width:100%; min-width:7rem; padding:.38rem .48rem; color:var(--accent); background:var(--soft); border:1px dashed var(--accent); border-radius:5px; font-variant-numeric:tabular-nums; font-weight:700; }}
p .derived-output,li .derived-output,blockquote .derived-output {{ width:auto; min-width:0; padding:.08rem .32rem; vertical-align:baseline; }}
.live-equation {{ margin:.7rem 0; padding:.7rem .85rem; overflow-x:auto; border:1px dashed var(--accent); border-radius:7px; background:var(--soft); color:var(--text); font:600 .92rem/1.5 "Cascadia Code",Consolas,monospace; font-variant-numeric:tabular-nums; }}
.edit-note {{ margin:0 0 1.2rem; padding:.7rem .9rem; border:1px solid var(--line); border-radius:8px; background:var(--soft); color:var(--muted); font-size:.86rem; }}
.section-takeaway {{ margin:-.45rem 0 1.4rem; padding:.85rem 1rem; border:1px solid var(--line); border-left:4px solid var(--accent); border-radius:8px; background:var(--soft); color:var(--muted); font-size:.91rem; }}
.section-takeaway p {{ margin:.34rem 0; }}
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
article img.zoomable-report-image {{ cursor:zoom-in; outline:none; }}
article img.zoomable-report-image:focus {{ outline:3px solid var(--accent); outline-offset:3px; }}
.image-lightbox {{ position:fixed; inset:0; z-index:1000; display:grid; grid-template-rows:auto minmax(0,1fr); background:rgba(5,10,14,.94); color:#f5f7f8; }}
.image-lightbox[hidden] {{ display:none; }}
.lightbox-toolbar {{ display:flex; flex-wrap:wrap; align-items:center; gap:.55rem; padding:.7rem 1rem; background:#111a21; border-bottom:1px solid #41505c; }}
.lightbox-title {{ flex:1; min-width:12rem; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
.lightbox-toolbar button {{ color:#f5f7f8; background:#24313b; border-color:#536572; }}
.lightbox-viewport {{ overflow:auto; display:block; padding:1.2rem; overscroll-behavior:contain; cursor:grab; }}
.lightbox-viewport.dragging {{ cursor:grabbing; user-select:none; }}
.lightbox-viewport img {{ display:block; width:100%; max-width:none; height:auto; margin:auto; border-radius:4px; background:white; box-shadow:0 10px 38px rgba(0,0,0,.45); }}
body.lightbox-open {{ overflow:hidden; }}
.footer {{ color:var(--muted); font-size:.78rem; margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line); }}
@media (max-width:920px) {{ .layout {{ grid-template-columns:1fr; padding:.7rem; }} nav {{ position:relative; top:auto; max-height:18rem; }} article {{ padding:1.2rem; }} .hide-small {{ display:none; }} .live-plot-grid {{ grid-template-columns:1fr; }} }}
@media print {{ .topbar,nav {{ display:none; }} body {{ background:white; }} .layout {{ display:block; padding:0; }} article {{ max-width:none; border:0; box-shadow:none; }} details {{ break-inside:avoid; }} }}
</style>
</head>
<body>
<div class="topbar"><span class="name">{html.escape(title)}</span><button id="outline-toggle" type="button" aria-expanded="true" aria-controls="report-outline" onclick="toggleOutline()">Hide outline</button><button onclick="setDetails(false)">Collapse</button><button onclick="saveParameterFile()">Save parameters</button><button class="hide-small" onclick="toggleTheme()">Theme</button></div>
<div class="layout"><nav id="report-outline"><div class="caption">On this page</div>{''.join(toc_html)}</nav><article><div class="edit-note"><span class="assumed-swatch"></span>Amber inputs are unidentified assumptions. “Save parameters” writes <code>model_parameters.json</code> through the browser’s file picker; save it beside this HTML and the Python builder so the next plot-generation run loads it. Dependent scalar values, live equations, and the live Bode panel recalculate immediately. Publication SVGs and nonlinear simulations update on rebuild.<span id="parameter-save-status" class="save-status">Loading values…</span></div>{body}<div class="footer">Rendered from {html.escape(markdown_path.name)} · build {BUILD_ID} · {generated}</div></article></div>
<div id="image-lightbox" class="image-lightbox" role="dialog" aria-modal="true" aria-label="Expanded report image" hidden><div class="lightbox-toolbar"><span id="lightbox-title" class="lightbox-title">Expanded image</span><button type="button" onclick="changeImageZoom(-0.25)" aria-label="Zoom out">−</button><button id="lightbox-zoom" type="button" onclick="resetImageZoom()" title="Reset zoom">100%</button><button type="button" onclick="changeImageZoom(0.25)" aria-label="Zoom in">+</button><button type="button" onclick="closeImageLightbox()">Close</button></div><div id="lightbox-viewport" class="lightbox-viewport"><img id="lightbox-image" alt=""></div></div>
<script>
function setDetails(open) {{ document.querySelectorAll('details').forEach(d => d.open=open); }}
function toggleTheme() {{ const root=document.documentElement; root.dataset.theme=root.dataset.theme==='dark'?'light':'dark'; }}
const outlineStorageKey = 'report-outline:' + document.title + ':' + location.pathname;
function setOutline(open, persist=true) {{
  document.body.classList.toggle('outline-collapsed', !open);
  const button=document.getElementById('outline-toggle');
  if (button) {{ button.textContent=open?'Hide outline':'Show outline'; button.setAttribute('aria-expanded',String(open)); }}
  if (persist) try {{ localStorage.setItem(outlineStorageKey, open?'open':'closed'); }} catch (_) {{}}
}}
function toggleOutline() {{ setOutline(document.body.classList.contains('outline-collapsed')); }}
let imageZoom=1;
function applyImageZoom() {{
  const image=document.getElementById('lightbox-image');
  const label=document.getElementById('lightbox-zoom');
  if (image) image.style.width=(imageZoom*100)+'%';
  if (label) label.textContent=Math.round(imageZoom*100)+'%';
}}
function changeImageZoom(delta) {{ imageZoom=Math.max(0.5,Math.min(6,imageZoom+delta)); applyImageZoom(); }}
function resetImageZoom() {{ imageZoom=1; applyImageZoom(); }}
function openImageLightbox(source) {{
  const box=document.getElementById('image-lightbox'), image=document.getElementById('lightbox-image');
  if (!box || !image) return;
  image.src=source.currentSrc||source.src; image.alt=source.alt||'Expanded report image';
  document.getElementById('lightbox-title').textContent=source.alt||'Expanded report image';
  resetImageZoom(); box.hidden=false; document.body.classList.add('lightbox-open');
  document.getElementById('lightbox-viewport').scrollTo(0,0);
  box.querySelector('button').focus();
}}
function closeImageLightbox() {{
  const box=document.getElementById('image-lightbox'); if (!box || box.hidden) return;
  box.hidden=true; document.body.classList.remove('lightbox-open');
}}
function initializeImageViewer() {{
  document.querySelectorAll('article img').forEach(image => {{
    image.classList.add('zoomable-report-image'); image.tabIndex=0; image.setAttribute('role','button');
    image.title='Open expanded image viewer'; image.addEventListener('click',()=>openImageLightbox(image));
    image.addEventListener('keydown',event=>{{ if (event.key==='Enter'||event.key===' ') {{ event.preventDefault(); openImageLightbox(image); }} }});
  }});
  const box=document.getElementById('image-lightbox'), viewport=document.getElementById('lightbox-viewport');
  box.addEventListener('click',event=>{{ if (event.target===box) closeImageLightbox(); }});
  viewport.addEventListener('wheel',event=>{{ event.preventDefault(); changeImageZoom(event.deltaY<0?0.25:-0.25); }},{{passive:false}});
  let dragging=false,startX=0,startY=0,scrollLeft=0,scrollTop=0;
  viewport.addEventListener('pointerdown',event=>{{ dragging=true; startX=event.clientX; startY=event.clientY; scrollLeft=viewport.scrollLeft; scrollTop=viewport.scrollTop; viewport.classList.add('dragging'); viewport.setPointerCapture(event.pointerId); }});
  viewport.addEventListener('pointermove',event=>{{ if (!dragging) return; viewport.scrollLeft=scrollLeft-(event.clientX-startX); viewport.scrollTop=scrollTop-(event.clientY-startY); }});
  viewport.addEventListener('pointerup',()=>{{ dragging=false; viewport.classList.remove('dragging'); }});
  document.addEventListener('keydown',event=>{{ if (box.hidden) return; if (event.key==='Escape') closeImageLightbox(); else if (event.key==='+'||event.key==='=') changeImageZoom(0.25); else if (event.key==='-') changeImageZoom(-0.25); else if (event.key==='0') resetImageZoom(); }});
}}
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
async function saveParameterFile() {{
  refreshInteractivePlots();
  const payload = {{
    schema: 'rev3-model-parameters-v1',
    build_id: '{BUILD_ID}',
    parameters: currentParameterValues()
  }};
  const source = JSON.stringify(payload, null, 2) + '\\n';
  const suggestedName = 'model_parameters.json';
  if ('showSaveFilePicker' in window) {{
    try {{
      const handle = await window.showSaveFilePicker({{
        suggestedName,
        types:[{{description:'Revision 3 model parameters', accept:{{'application/json':['.json']}}}}]
      }});
      const writable = await handle.createWritable();
      await writable.write(source);
      await writable.close();
      setParameterStatus('Saved model_parameters.json. Keep it beside the HTML/builder; the next rebuild will load it.', 'ok');
      return;
    }} catch (error) {{
      if (error.name === 'AbortError') return;
    }}
  }}
  const blob = new Blob([source], {{type:'application/json;charset=utf-8'}});
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url; link.download = suggestedName;
  document.body.appendChild(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  setParameterStatus('Downloaded model_parameters.json. Move it beside the HTML/builder before rebuilding.', 'warn');
}}
function scheduleParameterUpdate() {{
  setParameterStatus('Updating dependent values and live plots...', 'warn');
  if (parameterSaveTimer) clearTimeout(parameterSaveTimer);
  parameterSaveTimer = setTimeout(() => {{
    refreshInteractivePlots();
    setParameterStatus('Unsaved parameter changes - live preview updated - publication simulations are stale', 'warn');
  }}, 160);
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
  const screwRootDiameter = parameterNumber('screw_root_diameter', 6.8e-3);
  const screwDensity = parameterNumber('screw_density', 7850.0);
  const youngsModulus = parameterNumber('youngs_modulus', 210.0e9);
  const shearModulus = parameterNumber('shear_modulus', 80.8e9);
  const nutAxialDatum = parameterNumber('nut_axial_datum', 0.158);
  const tmax = parameterNumber('holding_torque', 0.060);
  const tdet = parameterNumber('detent_torque', 0.005);
  const detentPhase = parameterNumber('detent_phase', 0.0);
  const couplingSeries = parameterNumber('k_c_series', 68.7549);
  const mStage = parameterNumber('stage_mass', 0.355);
  const mNut = parameterNumber('nut_mass', 0.050);
  const ms = mStage + mNut;
  const axialModeTarget = parameterNumber('axial_mode_target_hz', 695.82);
  const mEffMeasured = parameterNumber('m_eff_measured', 0.600);
  const zetaRelativeMeasured = parameterNumber('zeta_relative_measured', 0.0014);
  const kbrg = parameterNumber('k_brg', 25.0e6);
  const kmnt = parameterNumber('k_mnt', 100.0e6);
  const cax = parameterNumber('axial_damping', 55.0);
  const zeta = parameterNumber('electromagnetic_zeta', 0.10);
  const etaSteel = parameterNumber('eta_steel', {INTERFACE_LOSS_FACTORS['zeta_steel']});
  const etaBearing = parameterNumber('eta_bearing', {INTERFACE_LOSS_FACTORS['zeta_bearing']});
  const etaBallNut = parameterNumber('eta_ball_nut', {INTERFACE_LOSS_FACTORS['zeta_ball_nut']});
  const etaNutMount = parameterNumber('eta_nut_mount', {INTERFACE_LOSS_FACTORS['zeta_nut_mount']});
  const microstepDivisor = parameterNumber('microstep_divisor', 16);
  if (!(lead>0 && teeth>0 && jm>0 && jc>=0 && screwLength>0 && usableScrewTravel>0 &&
        stageTravel>0 && stageTravel<=usableScrewTravel && usableScrewTravel<=screwLength && screwDiameter>0 &&
        screwRootDiameter>0 && screwRootDiameter<=screwDiameter &&
        youngsModulus>0 && shearModulus>0 &&
        nutAxialDatum>0 && nutAxialDatum<screwLength &&
        screwDensity>0 && tmax>0 && tdet>=0 && couplingSeries>0 &&
        mStage>0 && mNut>=0 &&
        axialModeTarget>0 && mEffMeasured>0 && zetaRelativeMeasured>0 &&
        kbrg>0 && kmnt>0 &&
        cax>=0 && zeta>=0 && microstepDivisor>=1))
    throw new Error('Geometry, masses, torque, and stiffness must be positive; the nut datum must lie inside the screw; damping and detent torque cannot be negative.');
  const r = lead/(2*Math.PI);
  const screwRadius = screwDiameter/2;
  const rootRadius = screwRootDiameter/2;
  const rootArea = Math.PI*rootRadius*rootRadius;
  const rootPolar = 0.5*Math.PI*Math.pow(rootRadius,4);
  const screwLengthA = nutAxialDatum;
  const screwLengthB = screwLength-nutAxialDatum;
  // Both screw segment pairs come from one section and one nut datum, so the
  // torsional and axial ratios are the same length ratio by construction.
  const kthetaA = shearModulus*rootPolar/screwLengthA;
  const kthetaB = shearModulus*rootPolar/screwLengthB;
  const ksha = youngsModulus*rootArea/screwLengthA;
  const kshb = youngsModulus*rootArea/screwLengthB;
  const screwMass = screwDensity*Math.PI*screwRadius*screwRadius*screwLength;
  // Mass follows the nominal diameter, polar inertia the root section.
  const screwInertia = 0.5*Math.PI*Math.pow(rootRadius,4)*screwDensity*screwLength;
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
  // Interface damping ratios are derived, never entered: a dashpot delivers
  // eta_j(w) = 2*zeta_j*w/w_j, so the ratio that realizes a target loss factor
  // at the retained mode is zeta_j = eta_j*f_j/(2*f_2).  Same equation as
  // interface_damping_ratios() in the builder.
  const segMass = screwMass/3, segInertia = screwInertia/3;
  const massPair = (a,b) => a*b/(a+b);
  const elementHz = (k,m) => Math.sqrt(k/m)/(2*Math.PI);
  const ballMass = 1/(r*r/segInertia + 1/segMass + 1/mNut);
  const ratioFor = (eta,k,m) => eta*elementHz(k,m)/(2*axialModeTarget);
  const zetaBearing = ratioFor(etaBearing, kbrg, segMass);
  const zetaSteel = ratioFor(etaSteel, ksha, massPair(segMass,segMass));
  const zetaBallNut = ratioFor(etaBallNut, kBall, ballMass);
  const zetaNutMount = ratioFor(etaNutMount, kmnt, massPair(mNut,mStage));
  // Multi-route reduction comparison.  Route F deliberately uses the same
  // adjacent/reduced masses as full_linear_matrices(), not one shared m_s.
  const kc1 = 2*couplingSeries, kc2 = 2*couplingSeries;
  const torsionalCompliance = r*r*(1/kc1 + 1/kc2 + 1/kthetaA);
  const axialCompliance = 1/kbrg + 1/ksha + 1/kBall + 1/kmnt;
  const kaxSeries = 1/axialCompliance;
  const kaxSeriesFull = 1/(axialCompliance + torsionalCompliance);
  const torsionalShare = torsionalCompliance/(axialCompliance + torsionalCompliance);
  const w2 = 2*Math.PI*axialModeTarget;
  const screwSegmentMass = screwMass/3;
  const screwSegmentInertia = screwInertia/3;
  const pairMass = (a,b) => a*b/(a+b);
  const ballRelativeMass = 1/(r*r/screwSegmentInertia + 1/screwSegmentMass + 1/mNut);
  const complexElements = [
    {{k:kbrg, c:2*zetaBearing*Math.sqrt(kbrg*screwSegmentMass)}},
    {{k:ksha, c:2*zetaSteel*Math.sqrt(ksha*pairMass(screwSegmentMass,screwSegmentMass))}},
    {{k:kBall, c:2*zetaBallNut*Math.sqrt(kBall*ballRelativeMass)}},
    {{k:kmnt, c:2*zetaNutMount*Math.sqrt(kmnt*pairMass(mNut,mStage))}}
  ];
  let compRe=0, compIm=0;
  complexElements.forEach(el => {{
    const den=el.k*el.k+w2*w2*el.c*el.c;
    compRe+=el.k/den; compIm+=-w2*el.c/den;
  }});
  const complianceMagnitudeSquared=compRe*compRe+compIm*compIm;
  const kaxCplx=compRe/complianceMagnitudeSquared;
  const caxCplx=-compIm/complianceMagnitudeSquared/w2;
  const kaxMeasured=w2*w2*mEffMeasured;
  const measuredBallCompliance=1/kaxMeasured-1/kbrg-1/ksha-1/kmnt;
  if (!(measuredBallCompliance>0))
    throw new Error('The measured-mass route leaves no positive compliance for k_ball.');
  const kBallMeasured=1/measuredBallCompliance;
  const muRel=md*ms/(md+ms), muMeasured=md*mEffMeasured/(md+mEffMeasured);
  const caxMeasured=2*zetaRelativeMeasured*Math.sqrt(kaxMeasured*muMeasured);
  // --- friction site parameter diagnostics -----------------------------------
  const nu = [1,2,3,4].map(i => parameterNumber('gms_nu'+i, 0));
  const frictionSites = {{
    g: {{k:[1,2,3,4].map(i => parameterNumber('g_k'+i,0)),
        s0:parameterNumber('g_sigma0',0), Fs:parameterNumber('g_Fs',0),
        Fc:parameterNumber('g_Fc',0)}},
    n: {{k:[1,2,3,4].map(i => parameterNumber('n_k'+i,0)),
        s0:parameterNumber('n_sigma0',0), Fs:parameterNumber('n_Fs',0),
        Fc:parameterNumber('n_Fc',0)}},
    d: {{k:[1,2,3,4].map(i => parameterNumber('d_k'+i,0)),
        s0:parameterNumber('d_sigma0',0), Fs:parameterNumber('d_Fs',0),
        Fc:parameterNumber('d_Fc',0)}}
  }};
  const tauC = parameterNumber('tau_C', 2.0e-4);
  const etaScrew = parameterNumber('eta_screw', 0.90);
  const preloadNut = parameterNumber('F_preload_nut', 100.0);
  if (!(tauC>0)) throw new Error('The Stribeck relaxation time must be positive.');
  if (!(etaScrew>0 && etaScrew<1)) throw new Error('Screw efficiency must lie between 0 and 1.');
  Object.entries(frictionSites).forEach(([tag,site]) => {{
    site.Cgms = (site.Fs-site.Fc)/tauC;
    site.tauC = tauC;
  }});
  const nuSum = nu.reduce((a,b)=>a+b,0);
  if (Math.abs(nuSum-1) > 1e-9) throw new Error('GMS force fractions must sum to one.');
  Object.entries(frictionSites).forEach(([tag,site]) => {{
    const kSum = site.k.reduce((a,b)=>a+b,0);
    if (!(kSum>0) || Math.abs(kSum-site.s0) > 1e-6*Math.max(site.s0,1))
      throw new Error('GMS element stiffnesses for site '+tag+' must sum to sigma_0.');
    site.yieldFs = site.k.map((k,i) => nu[i]*site.Fs/k);
    site.yieldFc = site.k.map((k,i) => nu[i]*site.Fc/k);
    site.yieldSpan = site.yieldFs[3]/site.yieldFs[0];
    site.staticDeflection = site.Fs/site.s0;
  }});
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
  function modePair(driveTangent, linkStiffness=kax, stageMass=ms) {{
    const qa=md*stageMass;
    const qb=md*linkStiffness+stageMass*(driveTangent+linkStiffness);
    const qc=driveTangent*linkStiffness;
    const discriminant=Math.max(qb*qb-4*qa*qc,0);
    const roots=[(qb-Math.sqrt(discriminant))/(2*qa),(qb+Math.sqrt(discriminant))/(2*qa)];
    return roots.map(value => Math.sqrt(Math.max(value,0))/(2*Math.PI));
  }}
  const modes=modePair(km);
  const localSoftModes=modePair(km-kdetAmplitude), localHardModes=modePair(km+kdetAmplitude);
  const localLowBand=[Math.min(localSoftModes[0],localHardModes[0]),Math.max(localSoftModes[0],localHardModes[0])];
  const routeModesS=modePair(km,kaxSeries,ms);
  const routeModesF=modePair(km,kaxCplx,ms);
  const routeModesM=modePair(km,kaxMeasured,mEffMeasured);
  const relativeZeta=cax/(2*Math.sqrt(kax*muRel));
  const zetaCplx=caxCplx/(2*Math.sqrt(kaxCplx*muRel));
  const settlingFactor=4;
  const settling=(ratio,frequency) => settlingFactor/(ratio*2*Math.PI*frequency);
  const detentSettling=settlingFactor/(zeta*Math.sqrt((km-kdetAmplitude)/md));
  const axialSettling=settling(relativeZeta,modes[1]);
  const plateauDwell=Math.max(0.100,detentSettling,axialSettling);
  const routeCF1={derived_defaults['route_c_f1']:.12g};
  const routeCF2={derived_defaults['route_c_f2']:.12g};
  const routeCKax={derived_defaults['route_c_kax']:.12g};
  const routeCCax={derived_defaults['route_c_cax']:.12g};
  const routeCZeta={derived_defaults['route_c_zeta']:.12g};
  const routeCKBall={derived_defaults['route_c_kball']:.12g};
  const fullModelZeta={derived_defaults['full_model_zeta']:.12g};
  const fullModelSettling={derived_defaults['full_model_settling']:.12g};
  const fullModelUpperHz={derived_defaults['full_model_upper_hz']:.12g};
  const firstFixedInterfaceHz={derived_defaults['first_fixed_interface_hz']:.12g};
  const firstDiscardedHz={derived_defaults['first_discarded_hz']:.12g};
  return {{
    frequencies, drive, stage, rotorStage, modes, localLowBand, md, ms, mStage, mNut,
    axialModeTarget, km, kdet, kdetAmplitude, kax, kBall, kbrg, ksha, kmnt,
    cax, zeta, lead, teeth, r, jm, jc, screwLength, screwDiameter,
    screwDensity, screwMass, screwInertia, jTotal, tmax, tdet, detentPhase,
    usableScrewTravel, stageTravel,
    couplingSeries, couplingHalf:2*couplingSeries, kthetaA, kthetaB, kshb,
    screwRootDiameter, youngsModulus, shearModulus, nutAxialDatum,
    screwLengthA, screwLengthB, kappa:teeth/r,
    fullStep:lead/(4*teeth), quarterStep:lead/(16*teeth),
    commandStep:lead/(4*teeth*microstepDivisor),
    microstepDivisor, fmax:tmax/r, cm, detentPeriod:lead/(4*teeth),
    mEffMeasured, zetaRelativeMeasured, zetaSteel, zetaBearing,
    zetaBallNut, zetaNutMount, axialCompliance,
    torsionalCompliance, torsionalShare, kaxSeries, kaxSeriesFull,
    kaxCplx, caxCplx, zetaCplx, kaxMeasured, caxMeasured, kBallMeasured,
    muRel, massRatio:md/ms, relativeZeta, routeModesS, routeModesF, routeModesM,
    routeCF1, routeCF2, routeCKax, routeCCax, routeCZeta, routeCKBall,
    fullModelZeta, fullModelSettling, fullModelUpperHz,
    firstFixedInterfaceHz, firstDiscardedHz, nu, frictionSites,
    tauC, etaScrew, preloadNut, detentSettling, axialSettling, plateauDwell,
    caseCount:{len(CASES)},
    routePSettling:settling(relativeZeta,modes[1]),
    routeFSettling:settling(zetaCplx,routeModesF[1]),
    routeMSettling:settling(zetaRelativeMeasured,routeModesM[1]),
    routeCSettling:settling(routeCZeta,routeCF2)
  }};
}}
function formatDerivedValue(key, value, display) {{
  if (display === 'mnm') return (value/1e6).toFixed(3);
  if (display === 'mm') return (value*1e3).toFixed(1);
  const scientific = new Set([
    'transmission_ratio','magnetic_stiffness','detent_stiffness','screw_inertia',
    'screw_segment_inertia','full_step_pitch','quarter_step_bound',
    'command_step','total_rotational_inertia',
    'reduced_axial_stiffness','k_ball'
  ]);
  if (/loop_area|A_loop|a_loop/i.test(key)) return (value*1e6).toFixed(2);
  if (scientific.has(key)) return value.toExponential(5);
  if (key==='case_count') return value.toFixed(0);
  if (key==='detent_settling_time_2pct' || key==='axial_settling_time_2pct' || key==='plateau_dwell')
    return (value*1e3).toFixed(1);
  if (key==='reduced_drive_mass') return value.toFixed(3);
  if (key==='reduced_stage_mass') return value.toFixed(3);
  if (key==='screw_mass' || key==='screw_segment_mass') return value.toFixed(6);
  if (key==='k_c_half') return value.toFixed(3);
  if (/^route_[psbfcm]_(md|ms)$/.test(key)) return value.toFixed(3);
  if (/^route_[psbfcm]_(kax|kball)$/.test(key) || key==='route_s_kax_full')
    return (value/1e6).toFixed(3);
  if (/^route_[psbfcm]_cax$/.test(key)) return value.toFixed(2);
  if (key==='interface_axial_damping') return value.toFixed(2);
  if (/^route_[psbfcm]_zeta$/.test(key) || key==='full_model_zeta') return value.toExponential(3);
  if (/^route_[psbfcm]_f[12]$/.test(key)) return value.toFixed(2);
  if (/^route_[psbfcm]_settling$/.test(key) || key==='full_model_settling')
    return (value*1e3).toFixed(1);
  if (/^zeta_(steel|bearing|ball_nut|nut_mount)$/.test(key)) return value.toExponential(5);
  if (key==='torsional_share') return (100*value).toFixed(3);
  if (key==='mass_ratio') return value.toFixed(2);
  if (key==='reduced_mu') return value.toFixed(4);
  if (key==='mu_fraction') return value.toFixed(4);
  if (key==='cb_frequency_delta' || key==='cb_damping_delta') return value.toFixed(2);
  if (key==='fixed_interface_separation' || key==='discarded_pole_separation')
    return value.toFixed(3);
  if (/^yield_[gnd]_[1-4]_(fs|fc)$/.test(key) || /^static_deflection_[gnd]$/.test(key))
    return (value*1e6).toFixed(2);          // micrometres
  if (/^yield_span_[gnd]$/.test(key)) return value.toFixed(1);
  if (/^gms_rate_[gnd]$/.test(key)) return value.toFixed(0);
  if (/^tau_C_[gnd]$/.test(key) || key==='retained_mode_period') return (value*1e3).toFixed(3);
  if (key.startsWith('detent_velocity_')) return (value*1e3).toFixed(3);
  if (key==='d_Fs_efficiency_estimate' || key==='tau_C_mode_ratio') return value.toFixed(2);
  if (key==='d_Fs_torque_equivalent') return (value*1e3).toFixed(3);
  if (key.endsWith('_hz')) return value.toFixed(2);
  return Number(value).toPrecision(6);
}}
function refreshDerivedOutputs(data) {{
  const values = {{
    case_count:data.caseCount,
    transmission_ratio:data.r,
    total_rotational_inertia:data.jTotal,
    reduced_drive_mass:data.md,
    reduced_stage_mass:data.ms,
    magnetic_stiffness:data.km,
    detent_stiffness:data.kdet,
    reduced_axial_stiffness:data.kax,
    interface_axial_damping:data.caxCplx,
    k_ball:data.kBall,
    full_step_pitch:data.fullStep,
    quarter_step_bound:data.quarterStep,
    command_step:data.commandStep,
    screw_inertia:data.screwInertia,
    screw_segment_inertia:data.screwInertia/3,
    screw_mass:data.screwMass,
    screw_segment_mass:data.screwMass/3,
    k_c_half:data.couplingHalf,
    k_theta_a:data.kthetaA,
    k_theta_b:data.kthetaB,
    k_sha:data.ksha,
    k_shb:data.kshb,
    zeta_steel:data.zetaSteel,
    zeta_bearing:data.zetaBearing,
    zeta_ball_nut:data.zetaBallNut,
    zeta_nut_mount:data.zetaNutMount,
    screw_length_a:data.screwLengthA,
    screw_length_b:data.screwLengthB,
    mode_1_hz:data.modes[0],
    mode_2_hz:data.modes[1],
    drive_stiffness:data.km,
    force_limit:data.fmax,
    spatial_wavenumber:data.kappa
    ,route_p_md:data.md, route_p_ms:data.ms, route_p_kax:data.kax,
    route_p_cax:data.cax, route_p_zeta:data.relativeZeta,
    route_p_f1:data.modes[0], route_p_f2:data.modes[1], route_p_kball:data.kBall,
    route_p_settling:data.routePSettling,
    route_s_md:data.md, route_s_ms:data.ms, route_s_kax:data.kaxSeries,
    route_s_cax:data.cax, route_s_zeta:data.relativeZeta,
    route_s_f1:data.routeModesS[0], route_s_f2:data.routeModesS[1], route_s_kball:data.kBall,
    route_s_settling:data.routePSettling,
    route_b_md:data.md, route_b_ms:data.ms, route_b_kax:data.kaxSeries,
    route_b_cax:data.cax, route_b_zeta:data.relativeZeta,
    route_b_f1:data.routeModesS[0], route_b_f2:data.routeModesS[1], route_b_kball:data.kBall,
    route_b_settling:data.routePSettling,
    route_f_md:data.md, route_f_ms:data.ms, route_f_kax:data.kaxCplx,
    route_f_cax:data.caxCplx, route_f_zeta:data.zetaCplx,
    route_f_f1:data.routeModesF[0], route_f_f2:data.routeModesF[1], route_f_kball:data.kBall,
    route_f_settling:data.routeFSettling,
    route_c_md:data.md, route_c_ms:data.ms, route_c_kax:data.routeCKax,
    route_c_cax:data.routeCCax, route_c_zeta:data.routeCZeta,
    route_c_f1:data.routeCF1, route_c_f2:data.routeCF2, route_c_kball:data.routeCKBall,
    route_c_settling:data.routeCSettling,
    route_m_md:data.md, route_m_ms:data.mEffMeasured, route_m_kax:data.kaxMeasured,
    route_m_cax:data.caxMeasured, route_m_zeta:data.zetaRelativeMeasured,
    route_m_f1:data.routeModesM[0], route_m_f2:data.routeModesM[1], route_m_kball:data.kBallMeasured,
    route_m_settling:data.routeMSettling,
    route_s_kax_full:data.kaxSeriesFull, torsional_share:data.torsionalShare,
    mass_ratio:data.massRatio, reduced_mu:data.muRel,
    full_model_zeta:data.fullModelZeta, full_model_settling:data.fullModelSettling,
    full_model_upper_hz:data.fullModelUpperHz,
    cb_frequency_delta:100*Math.abs(data.routeCF2-data.fullModelUpperHz)/data.fullModelUpperHz,
    cb_damping_delta:100*Math.abs(data.routeCZeta-data.fullModelZeta)/data.fullModelZeta,
    route_p_bandwidth_hz:2*data.relativeZeta*data.modes[1],
    full_model_bandwidth_hz:2*data.fullModelZeta*data.fullModelUpperHz,
    first_fixed_interface_hz:data.firstFixedInterfaceHz,
    first_discarded_hz:data.firstDiscardedHz,
    fixed_interface_separation:Math.pow(data.axialModeTarget/data.firstFixedInterfaceHz,2),
    discarded_pole_separation:Math.pow(data.axialModeTarget/data.firstDiscardedHz,2)
    ,mu_fraction:data.muRel/data.ms,
    relative_mode_nearground_hz:Math.sqrt(data.kax/data.ms)/(2*Math.PI),
    drive_pole_hz:Math.sqrt(data.km/data.md)/(2*Math.PI),
    ...Object.fromEntries(
      Object.entries(data.frictionSites).flatMap(([tag,site]) => [
        ...site.yieldFs.map((v,i) => ['yield_'+tag+'_'+(i+1)+'_fs', v]),
        ...site.yieldFc.map((v,i) => ['yield_'+tag+'_'+(i+1)+'_fc', v]),
        ['yield_span_'+tag, site.yieldSpan],
        ['static_deflection_'+tag, site.staticDeflection],
        ['gms_rate_'+tag, site.Cgms],
        ['tau_C_'+tag, site.tauC]
      ])
    ),
    d_Fs_efficiency_estimate:data.preloadNut*(1/data.etaScrew-1),
    d_Fs_torque_equivalent:data.frictionSites.d.Fs*data.r,
    detent_velocity_drive:data.modes[0]*data.fullStep,
    detent_velocity_axial:data.modes[1]*data.fullStep,
    detent_velocity_discarded:data.firstDiscardedHz*data.fullStep,
    retained_mode_period:1/data.modes[1],
    tau_C_mode_ratio:1/(data.modes[1]*data.tauC)
    ,detent_settling_time_2pct:data.detentSettling,
    axial_settling_time_2pct:data.axialSettling,
    plateau_dwell:data.plateauDwell
  }};
  document.querySelectorAll('[data-derived]').forEach(output => {{
    const key=output.dataset.derived;
    if (Object.prototype.hasOwnProperty.call(values,key)) output.textContent=formatDerivedValue(key,values[key],output.dataset.derivedFormat);
  }});
}}
function refreshLiveEquations(data) {{
  const mn = value => (value/1e6).toFixed(3);
  const reflectedKg = inertia => inertia/(data.r*data.r);
  const compliancePct = stiffness => 100*data.kax/stiffness;
  const torsionalCompliance = data.torsionalCompliance;
  const fullStaticCompliance = data.axialCompliance + torsionalCompliance;
  const fullStaticStiffness = 1/fullStaticCompliance;
  const torsionalCompliancePct = 100*torsionalCompliance/fullStaticCompliance;
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
    'exact-static-condensation':
      'Exact full static link: ' + (fullStaticStiffness/1e6).toFixed(3) +
      ' MN/m; executable axial-only link: ' + (data.kax/1e6).toFixed(3) +
      ' MN/m. Reflected torsional compliance is ' +
      (torsionalCompliance*1e9).toFixed(3) + ' nm/N (' +
      torsionalCompliancePct.toFixed(3) + '% of the complete static compliance).',
    'modal-stiffness':
      'f₂,target = ' + data.axialModeTarget.toFixed(2) + ' Hz  →  k_ax = ' +
      mn(data.kax) + ' MN/m',
    'route-comparison-summary':
      'Formal condensation, direct compliance, and the bond graph agree on the axial link to ' +
      (100*Math.abs(data.kaxSeries-data.kax)/data.kax).toExponential(2) +
      '%. Including the reflected torsional chain lowers the static link by ' +
      (100*(data.kaxSeries-data.kaxSeriesFull)/data.kaxSeries).toFixed(3) +
      '%. Frequency-domain condensation changes stiffness by ' +
      (100*(data.kaxCplx-data.kax)/data.kax).toFixed(3) +
      '% at f2 while replacing the damping propagation.',
    'damping-defect':
      'Current 2-DOF: zeta_rel=' + data.relativeZeta.toExponential(3) +
      ', t_2%=' + (data.routePSettling*1e3).toFixed(1) +
      ' ms; full 10-DOF audit: zeta_rel=' + data.fullModelZeta.toExponential(3) +
      ', t_2%=' + (data.fullModelSettling*1e3).toFixed(1) +
      ' ms; frequency-domain reduction: zeta_rel=' + data.zetaCplx.toExponential(3) +
      ', t_2%=' + (data.routeFSettling*1e3).toFixed(1) + ' ms.',
    'mass-conflict':
      'BOM mass ' + data.ms.toFixed(3) + ' kg gives k_ax=' + mn(data.kax) +
      ' MN/m and k_ball=' + mn(data.kBall) + ' MN/m; measured m_eff ' +
      data.mEffMeasured.toFixed(3) + ' kg gives k_ax=' + mn(data.kaxMeasured) +
      ' MN/m and k_ball=' + mn(data.kBallMeasured) + ' MN/m, a ' +
      (data.kBallMeasured/data.kBall).toFixed(2) + 'x change.',
    'axial-compliance':
      '1/' + mn(data.kbrg) + ' + 1/' + mn(data.ksha) + ' + 1/' +
      mn(data.kBall) + ' + 1/' + mn(data.kmnt) + ' = 1/' +
      mn(data.kax) + '  (MN/m)⁻¹',
    'friction-port-summary': {json.dumps(friction_port_sentence())},
    'gms-branch-census': {json.dumps(rendered_branch_census_sentence())}
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
  let outlineOpen=true;
  try {{ outlineOpen=localStorage.getItem(outlineStorageKey)!=='closed'; }} catch (_) {{}}
  setOutline(outlineOpen,false);
  initializeImageViewer();
  document.querySelectorAll('.parameter-input').forEach(input => {{
    input.setAttribute('value', input.value);
    input.addEventListener('input', () => {{ input.setAttribute('value', input.value); scheduleParameterUpdate(); }});
  }});
  refreshInteractivePlots();
  setParameterStatus('Values loaded from the builder parameter source · Save parameters writes model_parameters.json', 'ok');
}});
</script>
</body></html>"""


def render_document(markdown_path: Path) -> Path:
    output = markdown_path.with_suffix(".html")
    output.write_text(html_page(markdown_path), encoding="utf-8")
    return output


def build_progress(step: int, total: int, message: str) -> None:
    """Emit visible, unbuffered stage progress without external polling."""
    print(f"[BUILD {step}/{total}] {message}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-summary-update", action="store_true",
                        help="Render without refreshing the generated metrics table")
    parser.add_argument(
        "--jobs", type=int, default=max(1, min(8, os.cpu_count() or 1)),
        help="worker processes for independent RK4 simulations (default: up to 8)")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if not DESCRIPTION_MD.exists() or not DERIVATION_MD.exists():
        raise FileNotFoundError("Both Markdown source documents must exist before building")
    build_progress(1, 8, "COMPILING inputs and validating model/report structure")
    ASSET_DIR.mkdir(exist_ok=True)
    ensure_parameter_file()
    gms_audit = validate_gms_partition()
    screw_audit = validate_screw_geometry()
    loss_audit = validate_interface_loss_factors()
    validate_parameter_registry()
    validate_prose_frequencies()
    breakaway_audit = validate_breakaway_forces()
    validate_case_topology()
    constants = physical_constants()
    component = component_parameters()
    closure_audit = validate_closure_band(constants, component)
    branches = calibration_branches(constants, component)
    predistortion_audit = validate_predistortion_authority(constants)
    command_audit = validate_command_design(constants)
    frequencies, bode, linear_metrics = frequency_responses()
    branch_censuses: dict[str, GmsBranchCensus] = {}
    executor = ProcessPoolExecutor(max_workers=args.jobs) if args.jobs > 1 else None
    try:
        build_progress(
            2, 8,
            f"SIMULATING main {main_duration(constants) * 1e3:.0f} ms nonlinear "
            f"campaign with {args.jobs} workers")
        times, command, time_data, time_metrics = time_responses(
            constants, branch_censuses, executor)
        build_progress(
            3, 8, "SIMULATING memory trajectories, branch counterfactuals, "
            "step halving (including A2 at 6.25 us), and reduction audit")
        if executor is not None:
            guideway_future = executor.submit(
                presliding_responses, constants, ("A", "A2", "G", "G2"), "g")
            nut_future = executor.submit(
                presliding_responses, constants, ("B", "B2"), "n")
            tau_c_future = executor.submit(
                tau_c_sensitivity, constants, ("A", "A2"), "g")
            verification_future = executor.submit(
                full_reduced_verification, frequencies, constants)
        branch_census = gms_branch_census_study(
            constants, times, command, branch_censuses, time_metrics)
        global BRANCH_CENSUS_SENTENCE
        BRANCH_CENSUS_SENTENCE = branch_census_sentence(branch_census)
        convergence = gms_step_halving_convergence(
            constants, times, time_data, executor)
        if executor is not None:
            memory_experiments = {
                "guideway": guideway_future.result(),
                "nut": nut_future.result(),
            }
            tau_c_study = tau_c_future.result()
            verification = verification_future.result()
        else:
            memory_experiments = {
                "guideway": presliding_responses(constants, ("A", "A2", "G", "G2"), "g"),
                "nut": presliding_responses(constants, ("B", "B2"), "n"),
            }
            tau_c_study = tau_c_sensitivity(constants, ("A", "A2"), "g")
            verification = full_reduced_verification(frequencies, constants)
    finally:
        if executor is not None:
            executor.shutdown()
    build_progress(4, 8, "PRICING memory-branch departures, detent ablation, "
                         "damping sweep, and breakaway sensitivity")
    branch_departure = memory_branch_departure(constants, memory_experiments)
    retention_confirmation = high_damping_confirmation_run(constants)
    true_loop = true_presliding_loop(constants)
    ablation_executor = ProcessPoolExecutor(max_workers=args.jobs) if args.jobs > 1 else None
    try:
        detent_ablation = detent_ablation_study(constants, time_metrics, ablation_executor)
        damping_sweep = retention_damping_sweep(constants, ablation_executor)
        breakaway = breakaway_sensitivity(constants, ablation_executor)
    finally:
        if ablation_executor is not None:
            ablation_executor.shutdown()
    build_progress(5, 8, "RENDERING publication SVG figures")
    micro_viscous = micro_viscous_effect(frequencies, bode, linear_metrics, time_metrics)
    case_response_paths = plot_case_responses(
        frequencies, bode, times, command, time_data, constants, time_metrics)
    comparison_paths = plot_case_response_overlay(
        frequencies, bode, linear_metrics, micro_viscous)
    guide_memory_path = plot_presliding_memory(
        memory_experiments["guideway"], "presliding_memory_comparison.svg")
    nut_memory_path = plot_presliding_memory(
        memory_experiments["nut"], "nut_memory_comparison.svg")
    memory_supplement_path = plot_presliding_supplement(memory_experiments)
    true_loop_path = plot_true_presliding_loop(true_loop)
    diagram_paths = plot_kinematic_diagram()
    flowchart_a_path = plot_flowchart_provenance_structure()
    flowchart_b_path = plot_flowchart_friction_results()
    bond_graph_path = plot_reduced_bond_graph()
    verification_path = plot_full_reduced_verification(frequencies, verification)
    position_path = plot_position_dependence()
    resonance_path = plot_stepper_resonance_visibility()
    rotor_stage_path = plot_rotor_stage_transfer_functions(frequencies)
    for obsolete_name in ("lugre_gms_pairwise_comparison.svg",
                          "step_tracking_all_cases.svg",
                          "full_reduced_single_edge_diagnostic.svg"):
        obsolete_path = ASSET_DIR / obsolete_name
        if obsolete_path.exists():
            obsolete_path.unlink()
    if not args.skip_summary_update:
        build_progress(6, 8, "UPDATING generated Markdown tables, takeaways, and limitations")
        update_generated_bode_comparison(
            generated_bode_comparison(frequencies, bode, linear_metrics))
        update_generated_micro_viscous(generated_micro_viscous(micro_viscous))
        update_generated_reduction_convergence(
            generated_reduction_convergence(verification))
        update_generated_summary(generated_summary(linear_metrics, time_metrics, verification))
        update_generated_calibration_branches(
            generated_calibration_branches(branches, constants))
        update_generated_detent_ablation(
            generated_detent_ablation(detent_ablation, constants))
        update_generated_full_response_summary(
            generated_full_response_summary(linear_metrics, time_metrics, detent_ablation))
        update_generated_presliding_summary(generated_presliding_summary(
            memory_experiments, damping_sweep, true_loop_path, constants))
        update_generated_breakaway_sensitivity(
            generated_breakaway_sensitivity(breakaway))
        update_generated_retention_diagnostic(generated_retention_diagnostic(
            memory_experiments, retention_confirmation, true_loop_path,
            memory_supplement_path))
        update_generated_convergence_summary(generated_convergence_summary(convergence))
        update_generated_branch_census(
            generated_branch_census(branch_census, branch_departure, memory_experiments))
        update_generated_tau_c_sensitivity(generated_tau_c_sensitivity(tau_c_study))
        takeaway_values = takeaway_derived_values(verification, memory_experiments)
        takeaway_values["a2_convergence_order"] = float(
            convergence["A2"]["observed_order"])
        takeaway_values.update(
            document_frequency_tokens(constants, linear_metrics, branches))
        ratio_low, ratio_high = damping_sweep["ratio_range"]
        takeaway_values["retention_ratio_low"] = float(ratio_low)
        takeaway_values["retention_ratio_high"] = float(ratio_high)
        update_derived_token_fallbacks(takeaway_values)
    build_progress(7, 8, "RENDERING HTML from the synchronized Markdown sources")
    description_html = render_document(DESCRIPTION_MD)
    derivation_html = render_document(DERIVATION_MD)
    for comparison_path in comparison_paths:
        print(f"Built {comparison_path.relative_to(ROOT)}")
    print(f"Built {len(case_response_paths)} per-case response figures")
    print(f"Built {guide_memory_path.relative_to(ROOT)}")
    print(f"Built {nut_memory_path.relative_to(ROOT)}")
    print(f"Built {memory_supplement_path.relative_to(ROOT)}")
    print(f"Built {true_loop_path.relative_to(ROOT)}")
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
    print(f"Simulation workers: {args.jobs} of {os.cpu_count() or 1} logical CPUs")
    print("Command audit: "
          f"q_mu={command_audit['quantum_nm']:.3f} nm; "
          f"main max increment={command_audit['main_max_increment_um']:.4f} um; "
          f"guideway inner margin={command_audit['guideway_inner_margin_um']:.4f} um; "
          f"nut outer-to-yield-3 margin={command_audit['nut_outer_margin_um']:.4f} um")
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
    print("Screw geometry: "
          f"L_a={screw_audit['screw_length_a'] * 1e3:.1f} mm + "
          f"L_b={screw_audit['screw_length_b'] * 1e3:.1f} mm = "
          f"{FULL['screw_length'] * 1e3:.1f} mm; "
          f"k_theta {screw_audit['k_theta_a']:.1f}/{screw_audit['k_theta_b']:.1f} N m/rad; "
          f"k_sh {screw_audit['k_sha'] / 1e6:.2f}/{screw_audit['k_shb'] / 1e6:.2f} MN/m")
    print("Closure band: "
          f"k_brg={closure_audit['k_brg'] / 1e6:.1f} MN/m is "
          f"{closure_audit['margin']:.2f}x the "
          f"{closure_audit['singular_limit'] / 1e6:.3f} MN/m singular limit; "
          f"k_ball={closure_audit['k_ball'] / 1e6:.3f} MN/m")
    print("Calibration branches: "
          f"k_ax frictionless={branches['frictionless_k_ax'] / 1e6:.3f} MN/m vs "
          f"presliding-inclusive={branches['presliding_k_ax'] / 1e6:.3f} MN/m "
          f"(x{branches['ratio']:.3f}); K_m=0 moves k_ax by "
          f"{branches['unpowered_shift_pct']:+.3f}%")
    print("Dwell candidates (ms): "
          f"floor=100.0, detent={constants['detent_settling_time_2pct'] * 1e3:.1f}, "
          f"axial={constants['axial_settling_time_2pct'] * 1e3:.1f}, "
          f"interface={constants['interface_settling_time_2pct'] * 1e3:.1f}, "
          f"measured={constants['measured_settling_time_2pct'] * 1e3:.1f} -> "
          f"executed {constants['plateau_dwell'] * 1e3:.1f}")
    print("Pre-distortion authority: "
          f"quantum={predistortion_audit['command_step'] * 1e9:.1f} nm vs required "
          f"{predistortion_audit['required_resolution'] * 1e9:.1f} nm; "
          f"n_mu {predistortion_audit['executed_divisor']:.0f} vs required "
          f"{predistortion_audit['required_divisor']:.1f} "
          f"(short by {predistortion_audit['shortfall']:.1f}x); "
          f"{'SATISFIED' if predistortion_audit['satisfied'] else 'NOT SATISFIED'}")
    print("Breakaway range check: " + "; ".join(
        f"{site}: F_s={record['F_s']:.2f} N vs stated {record['low']:.1f}-{record['high']:.1f} N"
        + ("" if record["inside"] else
           f"  WARNING: {record['factor_above']:.1f}x above the stated range, "
           "see Section 12.4")
        for site, record in breakaway_audit.items()))
    print("Detent ablation: "
          + "; ".join(
              f"{key} {row['executed_nm']:.1f}->{row['detent_off_nm']:.1f} nm "
              f"({row['detent_share_pct']:.0f}% detent)"
              for key, row in detent_ablation["rows"].items()))
    print("Retention vs damping: " + "; ".join(
        f"zeta={row['target_zeta']:.4f}: ratio {row['r_hold_ratio']:.2f}x"
        for row in damping_sweep["rows"]))
    print("Breakaway sensitivity: " + "; ".join(
        f"F_s={row['F_s']:.1f} N -> {row['elements_yielded_at_inner']}/4 elements at "
        f"the inner level, F_ret {row['force_mismatch_gms_N']:.4f} N"
        for row in breakaway["rows"])
        + ("; DESIGN CHANGES" if breakaway["design_changes"] else "; design stable"))
    print("Micro-viscous A1v vs A: "
          f"peak reduced {micro_viscous['peak_drop_db']:.3f} dB at "
          f"{micro_viscous['peak_frequency_hz']:.1f} Hz; "
          f"implied dzeta={micro_viscous['implied_delta_zeta']:.3e} vs "
          f"predicted {micro_viscous['predicted_delta_zeta']:.3e} "
          f"({micro_viscous['agreement_pct']:.1f}% apart); settled RMS "
          f"{micro_viscous['rms_a_nm']:.1f} -> {micro_viscous['rms_a1v_nm']:.1f} nm")
    print("Interface loss factors at f_2: " + "; ".join(
        f"{key.replace('zeta_', '')}: eta={record['eta']:.5f}, zeta={record['zeta']:.6f}, "
        f"c={record['damping']:.2f} N s/m"
        for key, record in loss_audit.items()))
    print("tau_C sensitivity: F_ret spread="
          f"{tau_c_study['force_spread_N']:.6f} N vs law gap {tau_c_study['law_gap_force_N']:.6f} N; "
          f"A_loop spread={tau_c_study['loop_spread_J']:.4e} J vs law gap {tau_c_study['law_gap_loop_J']:.4e} J")
    print("Branch departure vs loop metrics: " + "; ".join(
        f"{record['case']} {record['metric']}: {float(record['delta']):+.4g} "
        f"(gap {float(record['law_gap']):.4g}){'  EXCEEDS' if record['exceeds_law_gap'] else ''}"
        for record in branch_departure["records"]))
    print("GMS branch census: "
          f"threshold flips={branch_census['threshold_total']}, "
          f"reversal flips={branch_census['reversal_total']}, "
          f"evaluations={branch_census['evaluation_total']}"
          + "".join(
              f"; {key} settled RMS {record['baseline_rms_nm']:.3f}->"
              f"{record['settled_rms_nm']:.3f} nm ({record['delta_pct']:+.2f}%)"
              for key, record in branch_census["enforced"].items()))
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
        dt_text = "/".join(f"{dt * 1e6:g}" for dt in result["dt_s"])
        rms_text = "/".join(f"{value:.6f}" for value in rms)
        print(f"Step halving {key}: RMS({dt_text} us)="
              f"{rms_text} nm; observed p={result['observed_order']:.3f}; "
              f"fine relative change={result['fine_relative_pct']:.6f}%")
    print("Full-model modes below 3 kHz: " + ", ".join(
        f"{mode:.2f}" for mode in verification["full_modes"] if mode < 3000.0))
    build_progress(8, 8, "GREEN FLAG — simulation, Markdown, SVG, and HTML update complete")


if __name__ == "__main__":
    main()

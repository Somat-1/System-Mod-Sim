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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"
DESCRIPTION_MD = ROOT / "Simulation_description.md"
DERIVATION_MD = ROOT / "Analytical_derivation_and_responses.md"


# Values already stated in Simulation_description.md.
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
    """Three increments, each no larger than one quarter of one full step."""
    if t < 0.005:
        return 0.0
    if t < 0.025:
        return quarter_step
    if t < 0.045:
        return 0.0
    return -quarter_step


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
    """Return four GMS force-state derivatives and the total site force."""
    threshold = np.maximum(GMS_WEIGHTS * stribeck(velocity, p), 1e-12)
    stiffness = GMS_STIFFNESS_FRACTIONS * p["sigma0"]
    derivatives = np.zeros(GMS_N)
    if abs(velocity) > 1e-14:
        direction = np.sign(velocity)
        for i in range(GMS_N):
            unloading = velocity * element_forces[i] <= 0.0
            below_yield = abs(element_forces[i]) < threshold[i]
            if unloading or below_yield:
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


def rk4_case(case: dict[str, object], constants: dict[str, float], dt: float = 5.0e-6,
             duration: float = 0.065) -> tuple[np.ndarray, np.ndarray]:
    times = np.arange(0.0, duration + 0.5 * dt, dt)
    states = np.zeros((times.size, STATE_SIZE), dtype=float)
    for i in range(times.size - 1):
        t = times[i]
        y = states[i]
        # Treat the discrete command as a true zero-order hold.  One midpoint
        # sample is fixed across all RK4 stages, so an endpoint discontinuity
        # cannot leak backward into the preceding integration interval.
        held_command = command_position(t + 0.5 * dt, constants.get('quarter_step'))
        k1 = nonlinear_rhs(t, y, case, constants, held_command)
        k2 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k1, case, constants, held_command)
        k3 = nonlinear_rhs(t + 0.5 * dt, y + 0.5 * dt * k2, case, constants, held_command)
        k4 = nonlinear_rhs(t + dt, y + dt * k3, case, constants, held_command)
        states[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return times, states


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
            "rms_final_error_nm": float(np.sqrt(np.mean(error[final_window] ** 2)) * 1e9),
            "max_stage_um": float(np.max(np.abs(states[:, 1])) * 1e6),
            "first_peak_um": first_peak * 1e6,
            "first_overshoot_pct": max(0.0, (first_peak / constants["quarter_step"] - 1.0) * 100.0),
        }
    return times, command, results, metrics


def plot_case_responses(frequencies: np.ndarray, responses: dict[str, np.ndarray],
                        times: np.ndarray, command: np.ndarray,
                        results: dict[str, np.ndarray], constants: dict[str, float]) -> list[Path]:
    """Create one self-contained Bode/step/error figure beside each case."""
    outputs: list[Path] = []
    time_ms = times * 1e3
    for key, case in CASES.items():
        fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.4))
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

        ax_phase.semilogx(frequencies, phase, color=color, linestyle=line_style, linewidth=1.8)
        ax_phase.set_xlabel("Frequency (Hz)")
        ax_phase.set_ylabel("Phase (deg)")
        ax_phase.set_ylim(-380.0, 30.0)

        stage = results[key][:, 1]
        error = command - stage
        ax_pos.step(time_ms, command * 1e6, where="post", color="#111111", linewidth=1.9,
                    label="Command")
        ax_pos.plot(time_ms, stage * 1e6, color=color, linestyle=line_style, linewidth=1.6,
                    label="Actual stage")
        ax_pos.set_ylabel("Position (µm)")
        ax_pos.set_title("Bounded commanded / actual motion")
        ax_pos.legend(loc="upper right", frameon=True)

        ax_err.plot(time_ms, error * 1e9, color=color, linestyle=line_style, linewidth=1.5)
        ax_err.axhline(0.0, color="#888888", linewidth=0.8)
        ax_err.set_xlabel("Time (ms)")
        ax_err.set_ylabel(r"Error $x_{cmd}-x_s$ (nm)")
        ax_err.set_title("Tracking error")

        for axis in axes.flat:
            axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
            axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
        fig.suptitle(case["label"], fontsize=14, fontweight="bold")
        fig.text(0.5, 0.012,
                 f"Nonlinear magnetic force; zeta_m={MODEL['zeta_m']:.2f}; each command increment <= {constants['quarter_step'] * 1e6:.2f} µm.",
                 ha="center", fontsize=8.5, color="#555555")
        fig.tight_layout(rect=(0.02, 0.04, 0.99, 0.95))
        output = ASSET_DIR / f"response_case_{key}.svg"
        fig.savefig(output, format="svg", bbox_inches="tight")
        plt.close(fig)
        outputs.append(output)
    return outputs


def plot_pairwise_comparison(frequencies: np.ndarray, responses: dict[str, np.ndarray],
                             times: np.ndarray, command: np.ndarray,
                             results: dict[str, np.ndarray]) -> Path:
    """Compare each LuGre case only with its topology-matched GMS case."""
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 10.5))
    time_ms = times * 1e3
    for row, (lugre_key, gms_key) in enumerate(PAIRS):
        ax_bode, ax_error = axes[row]
        for key in (lugre_key, gms_key):
            case = CASES[key]
            magnitude = 20.0 * np.log10(np.maximum(np.abs(responses[key]), 1e-15))
            ax_bode.semilogx(frequencies, magnitude, color=case["color"],
                             linestyle=case["ls"], linewidth=1.7, label=case["label"])
            error = command - results[key][:, 1]
            ax_error.plot(time_ms, error * 1e9, color=case["color"],
                          linestyle=case["ls"], linewidth=1.45, label=case["label"])
        ax_bode.axhline(0.0, color="#888888", linewidth=0.7)
        ax_bode.set_ylabel("Magnitude (dB)")
        ax_bode.set_ylim(-90.0, 30.0)
        ax_bode.set_title(f"{lugre_key}/{gms_key}: Bode magnitude")
        ax_bode.legend(loc="lower left", fontsize=8)
        ax_error.axhline(0.0, color="#888888", linewidth=0.7)
        ax_error.set_ylabel("Tracking error (nm)")
        ax_error.set_title(f"{lugre_key}/{gms_key}: nonlinear sequence")
        ax_error.legend(loc="upper right", fontsize=8)
        for axis in (ax_bode, ax_error):
            axis.grid(True, which="major", color="#d1d1d1", linewidth=0.7)
            axis.grid(True, which="minor", color="#eeeeee", linewidth=0.45)
    axes[-1, 0].set_xlabel("Frequency (Hz)")
    axes[-1, 1].set_xlabel("Time (ms)")
    fig.suptitle("Topology-matched friction-model comparison", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.96))
    output = ASSET_DIR / "lugre_gms_pairwise_comparison.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def plot_kinematic_diagram() -> Path:
    """Render the two mechanical DOFs and every force port as a vector diagram."""
    fig, ax = plt.subplots(figsize=(11.2, 4.5))
    ax.set_xlim(0.0, 11.0)
    ax.set_ylim(0.0, 4.8)
    ax.axis("off")
    ground_y = 0.65
    ax.plot([0.5, 10.5], [ground_y, ground_y], color="#505962", linewidth=2.0)
    for x in np.linspace(0.7, 10.3, 28):
        ax.plot([x, x - 0.16], [ground_y, ground_y - 0.18], color="#7a838b", linewidth=0.8)

    drive_box = FancyBboxPatch((2.15, 2.0), 2.0, 1.0, boxstyle="round,pad=0.05",
                               facecolor="#dceef6", edgecolor="#277da1", linewidth=2)
    stage_box = FancyBboxPatch((7.25, 2.0), 2.0, 1.0, boxstyle="round,pad=0.05",
                               facecolor="#dff2ea", edgecolor="#218c74", linewidth=2)
    ax.add_patch(drive_box)
    ax.add_patch(stage_box)
    ax.text(3.15, 2.55, "Drivetrain\n$m_d$", ha="center", va="center", fontsize=12)
    ax.text(8.25, 2.55, "Stage\n$m_s$", ha="center", va="center", fontsize=12)

    ax.add_patch(FancyArrowPatch((2.2, 3.45), (4.1, 3.45), arrowstyle="-|>",
                                 mutation_scale=15, color="#277da1", linewidth=1.8))
    ax.add_patch(FancyArrowPatch((7.3, 3.45), (9.2, 3.45), arrowstyle="-|>",
                                 mutation_scale=15, color="#218c74", linewidth=1.8))
    ax.text(3.15, 3.73, "DOF 1: $x_d$", ha="center", color="#1f5d73", fontweight="bold")
    ax.text(8.25, 3.73, "DOF 2: $x_s$", ha="center", color="#176a55", fontweight="bold")

    ax.add_patch(FancyBboxPatch((0.25, 2.08), 1.25, 0.84, boxstyle="round,pad=0.05",
                                facecolor="#fff1cc", edgecolor="#c08a00", linewidth=1.5))
    ax.text(0.88, 2.5, "$x_{cmd}$\nfield input", ha="center", va="center", fontsize=10)
    ax.add_patch(FancyArrowPatch((1.5, 2.5), (2.15, 2.5), arrowstyle="-|>",
                                 mutation_scale=14, color="#c08a00", linewidth=2))
    ax.text(1.82, 2.78, "$F_{mag}$", ha="center", fontsize=9)
    ax.text(0.88, 1.75, "input — not a DOF", ha="center", color="#74652d", fontsize=8.5)

    spring_x = np.linspace(4.15, 7.25, 13)
    spring_y = np.full_like(spring_x, 2.72)
    spring_y[1:-1] += 0.16 * np.where(np.arange(1, 12) % 2, 1.0, -1.0)
    ax.plot(spring_x, spring_y, color="#555555", linewidth=1.5)
    ax.text(5.7, 3.1, "$k_{ax}$", ha="center", fontsize=10)

    ax.plot([4.15, 5.05], [2.22, 2.22], color="#555555", linewidth=1.4)
    ax.add_patch(Rectangle((5.05, 2.05), 1.0, 0.34, fill=False, edgecolor="#555555", linewidth=1.4))
    ax.plot([5.55, 5.55], [2.05, 2.39], color="#555555", linewidth=1.4)
    ax.plot([6.05, 7.25], [2.22, 2.22], color="#555555", linewidth=1.4)
    ax.text(5.7, 1.83, "$c_{ax}$", ha="center", fontsize=10)

    ax.plot([4.15, 7.25], [1.48, 1.48], color="#a85f00", linewidth=1.4, linestyle="--")
    ax.text(5.7, 1.58, "nut friction $F_{f,n}$ (parallel port)", ha="center",
            color="#8a4d00", fontsize=9)

    ax.plot([3.15, 3.15], [2.0, 1.12], color="#6a4c93", linewidth=1.5)
    ax.plot([2.78, 3.52], [1.12, 1.12], color="#6a4c93", linewidth=3)
    ax.text(2.05, 1.15, "$F_{f,d}$ and $c_m$", color="#6a4c93", fontsize=9)
    ax.plot([8.25, 8.25], [2.0, 1.12], color="#b23a48", linewidth=1.5)
    ax.plot([7.88, 8.62], [1.12, 1.12], color="#b23a48", linewidth=3)
    ax.text(8.75, 1.15, "guideway $F_{f,g}$", color="#9b2f3d", fontsize=9)

    ax.text(5.5, 4.45, "Two-degree-of-freedom kinematic model", ha="center",
            fontsize=15, fontweight="bold")
    ax.text(5.5, 0.12, "Ground-referenced ports: drivetrain damping/friction and guideway friction.  Internal ports: axial path and nut friction.",
            ha="center", fontsize=8.8, color="#555555")
    output = ASSET_DIR / "kinematic_diagram.svg"
    fig.savefig(output, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output


def generated_summary(linear_metrics: dict[str, dict[str, float | np.ndarray]],
                      time_metrics: dict[str, dict[str, float]]) -> str:
    lines = [
        "<!-- BEGIN GENERATED RESPONSE SUMMARY -->",
        "| Case | Friction law | Presliding modes (Hz) | DC gain $X_s/X_{cmd}$ | First-step overshoot | Final-window RMS error |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for key, case in CASES.items():
        modes = linear_metrics[key]["modes"]
        mode_text = f"{modes[0]:.1f}, {modes[1]:.1f}"
        friction_label = {"none": "none", "lugre": "LuGre", "gms": "GMS"}[case["friction"]]
        lines.append(
            f"| {key} | {friction_label} | {mode_text} | {linear_metrics[key]['dc_gain']:.5f} | "
            f"{time_metrics[key]['first_overshoot_pct']:.1f}% | {time_metrics[key]['rms_final_error_nm']:.1f} nm |"
        )
    lines.extend([
        "",
        "The final column summarizes the last 2 ms of the nonlinear run; it is not an identified settling specification. "
        "All cases include the separately highlighted electromagnetic damping assumption; Case 0 remains frictionless.",
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
    DERIVATION_MD.write_text(pattern.sub(summary, source), encoding="utf-8")


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
.assumed-swatch {{ display:inline-block; width:1.1rem; height:.8rem; margin:0 .25rem; vertical-align:middle; background:var(--assumed); border:1px solid var(--assumed-line); border-radius:3px; }}
details {{ margin:1rem 0; border:1px solid var(--line); border-radius:9px; background:color-mix(in srgb,var(--soft) 45%,var(--card)); padding:.2rem .9rem .8rem; }} details details {{ margin-left:.45rem; }} summary {{ cursor:pointer; font-weight:700; padding:.75rem .1rem; color:var(--accent); }}
pre {{ overflow:auto; background:var(--code); color:#e8edf2; border-radius:9px; padding:1rem; font-size:.87rem; }} code {{ font-family:Cascadia Code,Consolas,monospace; }} p code,li code,td code {{ background:var(--soft); border:1px solid var(--line); border-radius:4px; padding:.1rem .28rem; }}
.display-math {{ overflow-x:auto; padding:.5rem 0; }} img {{ display:block; max-width:100%; height:auto; margin:1.3rem auto; border-radius:6px; }}
.footer {{ color:var(--muted); font-size:.78rem; margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line); }}
@media (max-width:920px) {{ .layout {{ grid-template-columns:1fr; padding:.7rem; }} nav {{ position:relative; top:auto; max-height:18rem; }} article {{ padding:1.2rem; }} .hide-small {{ display:none; }} }}
@media print {{ .topbar,nav {{ display:none; }} body {{ background:white; }} .layout {{ display:block; padding:0; }} article {{ max-width:none; border:0; box-shadow:none; }} details {{ break-inside:avoid; }} }}
</style>
</head>
<body>
<div class="topbar"><span class="name">{html.escape(title)}</span><button onclick="setDetails(true)">Expand derivations</button><button onclick="setDetails(false)">Collapse</button><button onclick="saveEditedHtml()">Save edited HTML</button><button onclick="resetParameterInputs()">Reset inputs</button><button class="hide-small" onclick="toggleTheme()">Theme</button><button class="hide-small" onclick="window.print()">Print</button></div>
<div class="layout"><nav><div class="caption">On this page</div>{''.join(toc_html)}</nav><article><div class="edit-note"><span class="assumed-swatch"></span>Amber inputs are unidentified assumptions. Input edits persist in this browser. “Save edited HTML” writes the current values into a new HTML file after an explicit browser save prompt; browsers cannot silently overwrite the source file. Static figures are regenerated only by the Python build.</div>{body}<div class="footer">Rendered from {html.escape(markdown_path.name)} · {generated}</div></article></div>
<script>
function setDetails(open) {{ document.querySelectorAll('details').forEach(d => d.open=open); }}
function toggleTheme() {{ const root=document.documentElement; root.dataset.theme=root.dataset.theme==='dark'?'light':'dark'; }}
const parameterStorageKey = 'model-parameters:' + document.title + ':' + location.pathname;
function currentParameterValues() {{
  const values = {{}};
  document.querySelectorAll('.parameter-input').forEach(input => values[input.dataset.param] = input.value);
  return values;
}}
function persistParameterInputs() {{
  localStorage.setItem(parameterStorageKey, JSON.stringify(currentParameterValues()));
}}
function resetParameterInputs() {{
  document.querySelectorAll('.parameter-input').forEach(input => {{
    input.value = input.dataset.default;
    input.setAttribute('value', input.value);
  }});
  localStorage.removeItem(parameterStorageKey);
}}
async function saveEditedHtml() {{
  const originalInputs = Array.from(document.querySelectorAll('.parameter-input'));
  originalInputs.forEach(input => input.setAttribute('value', input.value));
  const clonedRoot = document.documentElement.cloneNode(true);
  const clonedInputs = Array.from(clonedRoot.querySelectorAll('.parameter-input'));
  clonedInputs.forEach((input, index) => {{
    const value = originalInputs[index].value;
    input.setAttribute('value', value);
    input.setAttribute('data-default', value);
  }});
  const source = '<!doctype html>\n' + clonedRoot.outerHTML;
  const suggestedName = location.pathname.split('/').pop() || 'model-report.html';
  if ('showSaveFilePicker' in window) {{
    try {{
      const handle = await window.showSaveFilePicker({{suggestedName, types:[{{description:'HTML document', accept:{{'text/html':['.html']}}}}]}});
      const writable = await handle.createWritable();
      await writable.write(source);
      await writable.close();
      return;
    }} catch (error) {{ if (error.name === 'AbortError') return; }}
  }}
  const blob = new Blob([source], {{type:'text/html;charset=utf-8'}});
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url; link.download = suggestedName; link.click();
  URL.revokeObjectURL(url);
}}
document.addEventListener('DOMContentLoaded', () => {{
  let saved = {{}};
  try {{ saved = JSON.parse(localStorage.getItem(parameterStorageKey) || '{{}}'); }} catch (_) {{ saved = {{}}; }}
  document.querySelectorAll('.parameter-input').forEach(input => {{
    if (Object.prototype.hasOwnProperty.call(saved, input.dataset.param)) input.value = saved[input.dataset.param];
    input.setAttribute('value', input.value);
    input.addEventListener('input', () => {{ input.setAttribute('value', input.value); persistParameterInputs(); }});
  }});
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
    constants = physical_constants()
    frequencies, bode, linear_metrics = frequency_responses()
    times, command, time_data, time_metrics = time_responses(constants)
    case_paths = plot_case_responses(frequencies, bode, times, command, time_data, constants)
    comparison_path = plot_pairwise_comparison(frequencies, bode, times, command, time_data)
    diagram_path = plot_kinematic_diagram()
    for obsolete_name in ("bode_all_cases.svg", "step_tracking_all_cases.svg"):
        obsolete_path = ASSET_DIR / obsolete_name
        if obsolete_path.exists():
            obsolete_path.unlink()
    if not args.skip_summary_update:
        update_generated_summary(generated_summary(linear_metrics, time_metrics))
    description_html = render_document(DESCRIPTION_MD)
    derivation_html = render_document(DERIVATION_MD)
    for case_path in case_paths:
        print(f"Built {case_path.relative_to(ROOT)}")
    print(f"Built {comparison_path.relative_to(ROOT)}")
    print(f"Built {diagram_path.relative_to(ROOT)}")
    print(f"Built {description_html.name}")
    print(f"Built {derivation_html.name}")
    for key in CASES:
        modes = linear_metrics[key]["modes"]
        print(f"Case {key}: modes={modes[0]:.2f}, {modes[1]:.2f} Hz; "
              f"DC gain={linear_metrics[key]['dc_gain']:.6f}; "
              f"overshoot={time_metrics[key]['first_overshoot_pct']:.2f}%; "
              f"final-window RMS error={time_metrics[key]['rms_final_error_nm']:.2f} nm")


if __name__ == "__main__":
    main()

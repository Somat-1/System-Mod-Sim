#!/usr/bin/env python3
"""Build the Rev 4 frictionless Bode response by Lagrangian assembly.

The structural matrices are assembled independently as sums of rank-one
element contributions k*J.T@J and c*J.T@J. Parameters are loaded from the
Rev 4.2 parameter file because it contains the complete current structural
set, but no LuGre parameter or state is used in this model.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REV4 = ROOT.parent
PARAMETER_FILE = REV4 / "lugre_friction" / "Rev 4.2" / "model_parameters.json"
ASSET_DIR = ROOT / "rendered_assets"
STATE_LABELS = ("theta_m", "theta_c", "theta_s", "theta_sb", "x_s", "x_n")


def load_parameters() -> dict[str, float]:
    payload = json.loads(PARAMETER_FILE.read_text(encoding="utf-8"))
    return payload["parameters"]


def outer(row: np.ndarray) -> np.ndarray:
    return np.outer(row, row)


def build_lagrange_matrices(
    p: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return M, C, K and command vector from scalar energy functions."""
    lead_ratio = p["L"] / (2.0 * np.pi)
    k_em = p["N_r"] * p["T_hold"]
    k_detent = 4.0 * p["N_r"] * p["T_d"]

    mass = np.diag([
        p["I_m"], p["I_c"], p["I_s"], p["I_sb"],
        p["M_screw"], p["M_s"],
    ])

    j_c = np.array([1.0, -1.0, 0.0, 0.0, 0.0, 0.0])
    j_s1 = np.array([0.0, 1.0, -1.0, 0.0, 0.0, 0.0])
    j_s2 = np.array([0.0, 0.0, 1.0, -1.0, 0.0, 0.0])
    j_nut = np.array([0.0, 0.0, lead_ratio, 0.0, 1.0, -1.0])
    j_brg = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
    j_motor = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    stiffness = (
        p["k_c"] * outer(j_c)
        + p["k_s1"] * outer(j_s1)
        + p["k_s2"] * outer(j_s2)
        + p["k_nut"] * outer(j_nut)
        + p["k_brg"] * outer(j_brg)
        + (k_em + k_detent) * outer(j_motor)
    )

    damping = (
        p["c_c"] * outer(j_c)
        + p["c_s1"] * outer(j_s1)
        + p["c_s2"] * outer(j_s2)
        + p["c_nut"] * outer(j_nut)
        + p["c_brg"] * outer(j_brg)
        + p["c_EM"] * outer(j_motor)
    )

    command = k_em * j_motor
    return mass, damping, stiffness, command


def load_module(name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def response(
    frequencies_hz: np.ndarray,
    mass: np.ndarray,
    damping: np.ndarray,
    stiffness: np.ndarray,
    command: np.ndarray,
) -> np.ndarray:
    output = np.empty(frequencies_hz.size, dtype=complex)
    for index, frequency in enumerate(frequencies_hz):
        omega = 2.0 * np.pi * frequency
        dynamic_stiffness = stiffness - omega**2 * mass + 1j * omega * damping
        output[index] = np.linalg.solve(dynamic_stiffness, command)[5]
    return output


def modal_data(
    mass: np.ndarray, damping: np.ndarray, stiffness: np.ndarray
) -> tuple[list[float], list[float]]:
    inverse_mass = np.diag(1.0 / np.diag(mass))
    state_matrix = np.block([
        [np.zeros((6, 6)), np.eye(6)],
        [-inverse_mass @ stiffness, -inverse_mass @ damping],
    ])
    eigenvalues = np.linalg.eigvals(state_matrix)
    modes = sorted(
        (value for value in eigenvalues if value.imag > 1.0e-6),
        key=lambda value: abs(value),
    )
    frequencies = [float(abs(value) / (2.0 * np.pi)) for value in modes]
    damping_ratios = [float(-value.real / abs(value)) for value in modes]
    return frequencies, damping_ratios


def write_bode_svg(
    path: Path,
    frequencies_hz: np.ndarray,
    magnitude_db: np.ndarray,
    phase_deg: np.ndarray,
) -> None:
    """Write a dependency-free two-panel SVG Bode chart."""
    width, height = 1200, 820
    left, right = 105.0, 35.0
    plot_width = width - left - right
    panel_height = 285.0
    top_magnitude, top_phase = 95.0, 475.0
    decimation = max(1, frequencies_hz.size // 10000)
    frequency = frequencies_hz[::decimation]

    def points(values: np.ndarray, top: float, low: float, high: float) -> str:
        clipped = np.clip(values[::decimation], low, high)
        x = left + plot_width * frequency / frequencies_hz[-1]
        y = top + panel_height * (high - clipped) / (high - low)
        return " ".join(f"{x_i:.2f},{y_i:.2f}" for x_i, y_i in zip(x, y))

    magnitude_low = float(np.floor(np.min(magnitude_db) / 20.0) * 20.0)
    magnitude_high = float(np.ceil(np.max(magnitude_db) / 20.0) * 20.0)
    phase_low = float(np.floor(np.min(phase_deg) / 90.0) * 90.0)
    phase_high = float(np.ceil(np.max(phase_deg) / 90.0) * 90.0)
    if magnitude_high <= magnitude_low:
        magnitude_high = magnitude_low + 20.0
    if phase_high <= phase_low:
        phase_high = phase_low + 90.0

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.axis{stroke:#333;stroke-width:1}.grid{stroke:#d6d6d6;stroke-width:1}.curve{fill:none;stroke-width:1.5}</style>',
        '<text x="600" y="38" text-anchor="middle" font-size="22">Rev 4 frictionless Lagrangian model</text>',
        '<text x="600" y="66" text-anchor="middle" font-size="17">x_n(s) / theta_cmd(s)</text>',
    ]

    for top, low, high, label in (
        (top_magnitude, magnitude_low, magnitude_high, "Magnitude (dB re 1 m/rad)"),
        (top_phase, phase_low, phase_high, "Phase (deg)"),
    ):
        for index in range(6):
            x = left + plot_width * index / 5.0
            frequency_tick = frequencies_hz[-1] * index / 5.0
            svg.append(f'<line class="grid" x1="{x}" y1="{top}" x2="{x}" y2="{top + panel_height}"/>')
            if top == top_phase:
                svg.append(f'<text x="{x}" y="{top + panel_height + 25}" text-anchor="middle" font-size="13">{frequency_tick:.0f}</text>')
        for index in range(6):
            y = top + panel_height * index / 5.0
            value = high - (high - low) * index / 5.0
            svg.append(f'<line class="grid" x1="{left}" y1="{y}" x2="{left + plot_width}" y2="{y}"/>')
            svg.append(f'<text x="{left - 12}" y="{y + 5}" text-anchor="end" font-size="13">{value:.0f}</text>')
        svg.append(f'<rect x="{left}" y="{top}" width="{plot_width}" height="{panel_height}" fill="none" class="axis"/>')
        svg.append(f'<text transform="translate(25 {top + panel_height / 2}) rotate(-90)" text-anchor="middle" font-size="15">{label}</text>')

    svg.append(f'<polyline class="curve" stroke="#2b6cb0" points="{points(magnitude_db, top_magnitude, magnitude_low, magnitude_high)}"/>')
    svg.append(f'<polyline class="curve" stroke="#c05621" points="{points(phase_deg, top_phase, phase_low, phase_high)}"/>')
    svg.append(f'<text x="{left + plot_width / 2}" y="{height - 18}" text-anchor="middle" font-size="15">Frequency (Hz)</text>')
    svg.append('</svg>')
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def main() -> None:
    parameters = load_parameters()
    mass, damping, stiffness, command = build_lagrange_matrices(parameters)

    parent_payload = json.loads(
        (REV4 / "model_parameters.json").read_text(encoding="utf-8")
    )
    parent_parameters = parent_payload["parameters"]
    mass_n, damping_n, stiffness_n, command_n = build_lagrange_matrices(
        parent_parameters
    )
    root_k_detent = 4.0 * parent_parameters["N_r"] * parent_parameters["T_d"]
    command_n = command_n + root_k_detent * np.array(
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )

    root_matrix_residuals = {
        "M": float(np.max(np.abs(mass - mass_n))),
        "C": float(np.max(np.abs(damping - damping_n))),
        "K": float(np.max(np.abs(stiffness - stiffness_n))),
        "G": float(np.max(np.abs(command - command_n))),
    }

    rev42 = load_module(
        "lugre_model_rev42",
        REV4 / "lugre_friction" / "Rev 4.2" / "scripts" / "lugre_model_rev42.py",
    )
    mass_42, damping_42, stiffness_42, command_42 = rev42.build_structural_matrices(
        parameters
    )
    j_motor = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    k_detent = 4.0 * parameters["N_r"] * parameters["T_d"]
    tangent_stiffness_42 = stiffness_42 + k_detent * outer(j_motor)
    rev42_residuals = {
        "M": float(np.max(np.abs(mass - mass_42))),
        "C": float(np.max(np.abs(damping - damping_42))),
        "K_tangent": float(np.max(np.abs(stiffness - tangent_stiffness_42))),
        "G": float(np.max(np.abs(command - command_42))),
    }
    if any(value > 1.0e-12 for value in rev42_residuals.values()):
        raise AssertionError(f"Lagrange/Rev 4.2 mismatch: {rev42_residuals}")

    frequencies_hz = np.linspace(0.0, 8000.0, 80001)
    lagrange_response = response(
        frequencies_hz, mass, damping, stiffness, command
    )
    rev42_response = response(
        frequencies_hz, mass_42, damping_42, tangent_stiffness_42, command_42
    )
    response_residual = float(np.max(np.abs(lagrange_response - rev42_response)))

    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(lagrange_response), 1.0e-30))
    phase_deg = np.unwrap(np.angle(lagrange_response)) * 180.0 / np.pi
    modes_hz, damping_ratios = modal_data(mass, damping, stiffness)
    dc_gain = float(np.real(lagrange_response[0]))
    ideal_gain = parameters["L"] / (2.0 * np.pi)

    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(
        ASSET_DIR / "bode_lagrange_frictionless_data.npz",
        frequencies_hz=frequencies_hz,
        response=lagrange_response,
        magnitude_db=magnitude_db,
        phase_deg=phase_deg,
        mass=mass,
        damping=damping,
        stiffness=stiffness,
        command=command,
    )

    summary = {
        "model": "Rev 4 frictionless six-coordinate Lagrangian model",
        "parameter_source": str(PARAMETER_FILE.relative_to(REV4)),
        "nonlinear_friction_states": 0,
        "detent_treatment": "small-signal tangent k_d=4*N_r*T_d at theta_m=0",
        "frequency_range_hz": [0.0, 8000.0],
        "matrix_max_abs_residual_vs_rev42_frictionless_tangent": rev42_residuals,
        "response_max_abs_residual_vs_rev42_frictionless_tangent": response_residual,
        "root_rev4_matrix_residuals": root_matrix_residuals,
        "root_rev4_command_note": (
            "Root build_bode_rev4.py uses G[0]=k_EM+k_d; the Lagrange and "
            "Rev 4.2 grounded-detent convention uses G[0]=k_EM."
        ),
        "dc_gain_m_per_rad": dc_gain,
        "ideal_lead_ratio_m_per_rad": ideal_gain,
        "modes_hz": modes_hz,
        "damping_ratios": damping_ratios,
    }
    (ASSET_DIR / "bode_lagrange_frictionless_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    write_bode_svg(
        ASSET_DIR / "bode_lagrange_frictionless.svg",
        frequencies_hz,
        magnitude_db,
        phase_deg,
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

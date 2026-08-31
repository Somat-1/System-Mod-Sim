#!/usr/bin/env python3
"""Generate the rest-point small-signal Bode for Guyan + nonlinear friction."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from guyan_friction_model import GuyanFrictionModel, N_STATES


ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "rendered_assets"


def response(A, B, C, D, frequencies):
    identity = np.eye(A.shape[0])
    values = np.empty(frequencies.size, dtype=complex)
    for i, frequency in enumerate(frequencies):
        s = 2j * np.pi * frequency
        values[i] = (C @ np.linalg.solve(s * identity - A, B) + D)[0, 0]
    return values


def main():
    ASSETS.mkdir(parents=True, exist_ok=True)
    model = GuyanFrictionModel()
    equilibrium = np.zeros(N_STATES)
    residual = np.max(np.abs(model.rhs(0.0, equilibrium, 0.0)))
    if residual > 1e-12:
        raise AssertionError(f"Zero state is not an equilibrium: {residual:.3e}")
    A, B, C, D = model.analytical_linearization(equilibrium)
    frequencies = np.logspace(-2, np.log10(8000.0), 2400)
    values = response(A, B, C, D, frequencies)
    magnitude = 20.0 * np.log10(np.maximum(np.abs(values), 1e-300))
    phase = np.unwrap(np.angle(values)) * 180.0 / np.pi

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 8), sharex=True, constrained_layout=True)
    axes[0].semilogx(frequencies, magnitude, color="#b2182b", lw=1.25)
    axes[1].semilogx(frequencies, phase, color="#2166ac", lw=1.25)
    axes[0].set_ylabel(r"Magnitude $|x_n/\theta_{cmd}|$ [dB re m/rad]")
    axes[1].set_ylabel("Phase [deg]"); axes[1].set_xlabel("Frequency [Hz]")
    for axis in axes: axis.grid(True, which="both", alpha=.25)
    fig.suptitle("Two-master Guyan model with parallel LuGre friction\n"
                 "Small-signal tangent at rest; exact periodic detent linearized only for this Bode")
    figure = ASSETS / "guyan_friction_bode.png"
    fig.savefig(figure, dpi=180); plt.close(fig)

    np.savez_compressed(ASSETS / "guyan_friction_bode_data.npz",
                        frequency_hz=frequencies, response=values,
                        magnitude_db=magnitude, phase_deg=phase,
                        A=A, B=B, C=C, D=D,
                        M=model.mass, damping=model.damping,
                        stiffness_without_detent=model.stiffness,
                        J_way=model.jacobians["way"], J_nut=model.jacobians["nut"],
                        J_sb=model.jacobians["sb"])
    summary = {
        "model": "two-master Guyan + exact periodic detent + parallel Rev 4.2 LuGre",
        "coordinates": ["theta_m", "x_n"],
        "nonlinear_state_count": N_STATES,
        "friction_ports": list(model.jacobians),
        "nut_k_c_retained": True,
        "detent_in_nonlinear_rhs": "T_d*sin(4*N_r*theta_m)",
        "bode_definition": "analytical local linearization at zero rest equilibrium",
        "frequency_range_hz": [float(frequencies[0]), float(frequencies[-1])],
        "equilibrium_rhs_max_abs": float(residual),
        "figure": figure.name,
    }
    (ASSETS / "guyan_friction_bode_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

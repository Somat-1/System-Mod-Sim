#!/usr/bin/env python3
"""Bode overlay: full 15-state Rev 4.2 LuGre tangent vs Rev 4 baseline."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lugre_model_rev42 import LuGreModelRev42, N_STATES, PORTS


ROOT = Path(__file__).resolve().parent.parent
REV4_DIR = ROOT.parents[1]
ASSET_DIR = ROOT / "rendered_assets"
NPZ_DIR = ASSET_DIR / "npz"
STAGE_VELOCITY = 5.0e-3
FREQUENCY_HZ = np.linspace(0.0, 8000.0, 80001)

sys.path.insert(0, str(REV4_DIR / "scripts"))
from build_bode_rev4 import (  # noqa: E402
    INPUT_LABELS,
    build_matrices as build_baseline_matrices,
    build_state_space as build_baseline_state_space,
    load_parameters as load_baseline_parameters,
)


def response(A: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    identity = np.eye(A.shape[0])
    result = np.empty(FREQUENCY_HZ.size, dtype=complex)
    for index, frequency in enumerate(FREQUENCY_HZ):
        result[index] = (
            c @ np.linalg.solve(1j * 2.0 * np.pi * frequency * identity - A, b)
        )[0]
    return result


def bode_arrays(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    magnitude = 20.0 * np.log10(np.maximum(np.abs(values), 1.0e-300))
    phase = np.unwrap(np.angle(values)) * 180.0 / np.pi
    return magnitude, phase


def modal_data(A: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    eigenvalues = np.linalg.eigvals(A)
    oscillatory = eigenvalues[np.imag(eigenvalues) > 1.0e-5]
    oscillatory = oscillatory[np.argsort(np.abs(oscillatory))]
    frequency = np.abs(oscillatory) / (2.0 * np.pi)
    damping = -np.real(oscillatory) / np.maximum(np.abs(oscillatory), 1.0e-300)
    return eigenvalues, frequency, damping


def main() -> None:
    model = LuGreModelRev42(enforce_interface_power=False)
    operating_state = model.cruise_state(STAGE_VELOCITY)
    A_lugre, b_lugre, c_lugre = model.analytical_linearization(operating_state)

    complex_step = 1.0e-30
    complex_jacobian = np.column_stack([
        np.imag(model.rhs(
            0.0,
            operating_state.astype(complex)
            + 1j * complex_step * np.eye(N_STATES)[column],
            0.0,
        )) / complex_step
        for column in range(N_STATES)
    ])
    jacobian_absolute_error = float(np.max(np.abs(A_lugre - complex_jacobian)))
    jacobian_relative_error = float(np.max(
        np.abs(A_lugre - complex_jacobian) / np.maximum(np.abs(complex_jacobian), 1.0)
    ))
    if jacobian_relative_error > 1.0e-10:
        raise RuntimeError(
            f"Analytical Jacobian check failed: relative error={jacobian_relative_error:.6e}"
        )

    baseline_parameters = load_baseline_parameters()
    M0, C0, K0, Bu0 = build_baseline_matrices(baseline_parameters)
    A0, B0, c0 = build_baseline_state_space(M0, C0, K0, Bu0)
    b0 = B0[:, INPUT_LABELS.index("theta_cmd")]

    lugre_response = response(A_lugre, b_lugre, c_lugre)
    baseline_response = response(A0, b0, c0)
    lugre_magnitude, lugre_phase = bode_arrays(lugre_response)
    baseline_magnitude, baseline_phase = bode_arrays(baseline_response)

    eig_lugre, modes_lugre, damping_lugre = modal_data(A_lugre)
    eig_baseline, modes_baseline, damping_baseline = modal_data(A0)
    if np.max(np.real(eig_lugre)) >= 0.0:
        raise RuntimeError(
            f"Rev 4.2 tangent is not asymptotically stable: "
            f"max(real(lambda))={np.max(np.real(eig_lugre)):.6e}"
        )

    zero_friction_parameters = dict(model.p)
    for port in PORTS:
        for coefficient in ("sigma0", "sigma1", "sigma2"):
            zero_friction_parameters[f"{coefficient}_{port}"] = 0.0
    zero_friction_model = LuGreModelRev42(
        zero_friction_parameters, enforce_interface_power=False
    )
    zero_state = zero_friction_model.cruise_state(STAGE_VELOCITY)
    A_zero, _, _ = zero_friction_model.analytical_linearization(zero_state)
    _, zero_modes, _ = modal_data(A_zero)
    zero_friction_mode_error = float(np.max(np.abs(zero_modes - modes_baseline)))
    if zero_friction_mode_error > 1.0e-8:
        raise RuntimeError(
            f"Zero-friction eigenvalue regression failed: {zero_friction_mode_error:.6e} Hz"
        )

    test_velocity = np.array([0.31, -0.27, 0.19, -0.13, 0.07, -0.03])
    test_force = np.array([0.23, -0.17, 0.11])
    stacked_jacobian = np.vstack([model.jacobians[port] for port in PORTS])
    generalized_reaction = -stacked_jacobian.T @ test_force
    virtual_power_error = float(abs(
        test_velocity @ generalized_reaction
        + test_force @ (stacked_jacobian @ test_velocity)
    ))

    observations = model.port_observables(operating_state)
    interface_power = float(sum(
        observations[port]["force"] * observations[port]["velocity"]
        for port in PORTS
    ))
    if interface_power < -1.0e-12:
        raise RuntimeError(f"Negative interface power at cruise: {interface_power:.6e}")

    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    data_path = NPZ_DIR / "bode_rev42_vs_frictionless.npz"
    np.savez(
        data_path,
        frequency_hz=FREQUENCY_HZ,
        baseline_response=baseline_response,
        lugre_response=lugre_response,
        baseline_magnitude_db=baseline_magnitude,
        lugre_magnitude_db=lugre_magnitude,
        baseline_phase_deg=baseline_phase,
        lugre_phase_deg=lugre_phase,
        baseline_eigenvalues=eig_baseline,
        lugre_eigenvalues=eig_lugre,
        baseline_modes_hz=modes_baseline,
        lugre_modes_hz=modes_lugre,
        baseline_damping_ratios=damping_baseline,
        lugre_damping_ratios=damping_lugre,
        operating_state=operating_state,
    )

    fig, (magnitude_axis, phase_axis) = plt.subplots(
        2, 1, figsize=(12.0, 8.5), sharex=True
    )
    magnitude_axis.plot(
        FREQUENCY_HZ, baseline_magnitude, color="#2b6cb0", linewidth=1.25,
        label="Rev 4 frictionless baseline",
    )
    magnitude_axis.plot(
        FREQUENCY_HZ, lugre_magnitude, color="#c0392b", linewidth=1.25,
        label="Rev 4.2 full 15-state LuGre tangent",
    )
    magnitude_axis.set_ylabel("Magnitude (dB)")
    magnitude_axis.set_title(
        r"Command-to-stage Bode overlay: $x_n(s)/\theta_{cmd}(s)$ "
        f"at V_stage={STAGE_VELOCITY * 1e3:.1f} mm/s"
    )
    magnitude_axis.grid(True, linewidth=0.4, color="#cccccc")
    magnitude_axis.legend(loc="best")

    phase_axis.plot(
        FREQUENCY_HZ, baseline_phase, color="#2b6cb0", linewidth=1.25
    )
    phase_axis.plot(
        FREQUENCY_HZ, lugre_phase, color="#c0392b", linewidth=1.25
    )
    phase_axis.set_xlabel("Frequency (Hz)")
    phase_axis.set_ylabel("Phase (deg)")
    phase_axis.set_xlim(FREQUENCY_HZ[0], FREQUENCY_HZ[-1])
    phase_axis.grid(True, linewidth=0.4, color="#cccccc")

    figure_path = ASSET_DIR / "bode_rev42_vs_frictionless.png"
    fig.tight_layout()
    fig.savefig(figure_path, dpi=160)
    plt.close(fig)

    magnitude_delta = lugre_magnitude - baseline_magnitude
    phase_delta = lugre_phase - baseline_phase
    summary = {
        "formulation": "15-state LuGre; v=J*qdot; generalized friction=-J.T*F",
        "stage_operating_velocity_m_per_s": STAGE_VELOCITY,
        "frequency_range_hz": [float(FREQUENCY_HZ[0]), float(FREQUENCY_HZ[-1])],
        "k_nut_N_per_m": model.p["k_nut"],
        "c_nut_N_s_per_m": model.p["c_nut"],
        "sigma1_target_zeta": model.p["sigma1_target_zeta"],
        "baseline_modes_hz": modes_baseline.tolist(),
        "rev42_modes_hz": modes_lugre.tolist(),
        "baseline_damping_ratios": damping_baseline.tolist(),
        "rev42_damping_ratios": damping_lugre.tolist(),
        "maximum_real_eigenvalue_rev42_per_s": float(np.max(np.real(eig_lugre))),
        "interface_power_at_cruise_W": interface_power,
        "analytical_jacobian_max_absolute_error": jacobian_absolute_error,
        "analytical_jacobian_max_relative_error": jacobian_relative_error,
        "zero_friction_maximum_mode_error_hz": zero_friction_mode_error,
        "jacobian_virtual_power_identity_error": virtual_power_error,
        "dc_gain_baseline_m_per_rad": float(np.real(baseline_response[0])),
        "dc_gain_rev42_m_per_rad": float(np.real(lugre_response[0])),
        "maximum_absolute_magnitude_delta_db": float(np.max(np.abs(magnitude_delta))),
        "frequency_of_maximum_magnitude_delta_hz": float(
            FREQUENCY_HZ[np.argmax(np.abs(magnitude_delta))]
        ),
        "maximum_absolute_phase_delta_deg": float(np.max(np.abs(phase_delta))),
        "figure": figure_path.name,
        "data": str(data_path.relative_to(ROOT)),
    }
    summary_path = ASSET_DIR / "bode_rev42_vs_frictionless.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {figure_path}")
    print(f"Wrote {data_path}")


if __name__ == "__main__":
    main()

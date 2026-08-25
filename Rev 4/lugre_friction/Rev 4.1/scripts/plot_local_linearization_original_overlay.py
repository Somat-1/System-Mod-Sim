#!/usr/bin/env python3
"""Overlay the current and pre-change locally linearized LuGre responses.

The current parameter file is read without modification.  The original
support-bearing parameters are the values stored in the repository version
of model_parameters.json before the current parameter edits.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lugre_model import load_parameters
from run_local_linearization_bode import V_STAGE, build_linearized_matrices


ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = ROOT / "rendered_assets"
NPZ_DIR = ASSET_DIR / "npz"
FREQUENCY_HZ = np.linspace(0.0, 2000.0, 20001)

ORIGINAL_SUPPORT_BEARING = {
    "sigma0_sb": 500.0,
    "sigma2_sb": 0.45,
    "Tc_sb": 1.5e-3,
    "Ts_sb": 2.2e-3,
    "vs_sb": 2.3e-4,
}


def frequency_response(parameters: dict[str, float]) -> np.ndarray:
    mass, stiffness, damping, command = build_linearized_matrices(parameters)
    response = np.empty(FREQUENCY_HZ.size, dtype=complex)
    for index, frequency_hz in enumerate(FREQUENCY_HZ):
        omega = 2.0 * np.pi * frequency_hz
        dynamic_stiffness = stiffness - omega**2 * mass + 1j * omega * damping
        response[index] = np.linalg.solve(dynamic_stiffness, command)[5]
    return response


def bode(response: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(response), 1.0e-300))
    phase_deg = np.unwrap(np.angle(response)) * 180.0 / np.pi
    return magnitude_db, phase_deg


def main() -> None:
    current = load_parameters()
    original = dict(current)
    original.update(ORIGINAL_SUPPORT_BEARING)

    current_response = frequency_response(current)
    original_response = frequency_response(original)
    current_magnitude, current_phase = bode(current_response)
    original_magnitude, original_phase = bode(original_response)

    ASSET_DIR.mkdir(exist_ok=True)
    NPZ_DIR.mkdir(exist_ok=True)
    data_path = NPZ_DIR / "local_linearization_original_current_overlay.npz"
    np.savez(
        data_path,
        frequency_hz=FREQUENCY_HZ,
        original_response=original_response,
        current_response=current_response,
        original_magnitude_db=original_magnitude,
        current_magnitude_db=current_magnitude,
        original_phase_deg=original_phase,
        current_phase_deg=current_phase,
        original_support_bearing=json.dumps(ORIGINAL_SUPPORT_BEARING),
        current_support_bearing=json.dumps(
            {key: current[key] for key in ORIGINAL_SUPPORT_BEARING}
        ),
    )

    fig, (magnitude_axis, phase_axis) = plt.subplots(
        2, 1, figsize=(11.0, 8.0), sharex=True
    )
    magnitude_axis.plot(
        FREQUENCY_HZ,
        original_magnitude,
        color="#2b6cb0",
        linewidth=1.35,
        label=(r"Original: $\sigma_{0,sb}=500$ N m/rad, "
               r"$\sigma_{2,sb}=0.45$ N m s/rad"),
    )
    magnitude_axis.plot(
        FREQUENCY_HZ,
        current_magnitude,
        color="#c0392b",
        linewidth=1.35,
        label=(r"Current: $\sigma_{0,sb}=0.076$ N m/rad, "
               r"$\sigma_{2,sb}=10^{-5}$ N m s/rad"),
    )
    magnitude_axis.set_ylabel("Magnitude (dB)")
    magnitude_axis.set_title(
        r"Rev 4 locally linearized LuGre: original vs current, "
        r"$x_n(s)/\theta_{cmd}(s)$"
        + f" at V_stage={V_STAGE * 1e3:.1f} mm/s"
    )
    magnitude_axis.grid(True, linewidth=0.4, color="#cccccc")
    magnitude_axis.legend(loc="best")

    phase_axis.plot(
        FREQUENCY_HZ, original_phase, color="#2b6cb0", linewidth=1.35
    )
    phase_axis.plot(
        FREQUENCY_HZ, current_phase, color="#c0392b", linewidth=1.35
    )
    phase_axis.set_xlabel("Frequency (Hz)")
    phase_axis.set_ylabel("Phase (deg)")
    phase_axis.set_xlim(0.0, 2000.0)
    phase_axis.grid(True, linewidth=0.4, color="#cccccc")

    figure_path = ASSET_DIR / "local_linearization_original_current_overlay.png"
    fig.tight_layout()
    fig.savefig(figure_path, dpi=160)
    plt.close(fig)

    max_delta_index = int(np.argmax(np.abs(current_magnitude - original_magnitude)))
    print(f"Maximum magnitude difference: "
          f"{abs(current_magnitude[max_delta_index] - original_magnitude[max_delta_index]):.4f} dB "
          f"at {FREQUENCY_HZ[max_delta_index]:.1f} Hz")
    print(f"Wrote {figure_path}")
    print(f"Wrote {data_path}")


if __name__ == "__main__":
    main()

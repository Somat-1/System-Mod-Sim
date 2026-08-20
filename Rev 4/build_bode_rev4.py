#!/usr/bin/env python3
"""Command-to-stage Bode response for the Rev 4 6-DOF state-space model.

Pipeline: assemble M/C/K/B_u from model_parameters.json -> first-order state
space (A, B, C_y) per state_space_6dof.md Sec. 6 -> Laplace-domain transfer
function G(s) = C_y (sI - A)^-1 B_col for the theta_cmd -> x_n channel,
with the friction-port inputs held at zero (Sec. 7: the plant is linear only
with friction outputs frozen) -> frequency sweep 0-2000 Hz -> magnitude (dB)
and phase (deg) -> Bode plot.

As of the 2026-08-17 recoupling, the screw-nut interface (k_nut, c_nut) is
embedded directly in K/C as a linear reflected coupling between theta_s and
the axial coordinates (x_s, x_n), scaled by the lead ratio L/2*pi, rather
than injected through a frozen T_fric,nut friction port. u is now
[theta_cmd, T_fric_sb, F_fric_way] (3 inputs, was 4).

Depends only on NumPy and Matplotlib, matching the Rev 3 build script.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
PARAMETER_FILE = ROOT / "model_parameters.json"
ASSET_DIR = ROOT / "rendered_assets"

STATE_LABELS = ["theta_m", "theta_c", "theta_s", "theta_sb", "x_s", "x_n"]
INPUT_LABELS = ["theta_cmd", "T_fric_sb", "F_fric_way"]


def load_parameters() -> dict[str, float]:
    payload = json.loads(PARAMETER_FILE.read_text(encoding="utf-8"))
    return payload["parameters"]


def build_matrices(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assemble M, C, K, B_u per Sec. 5 of state_space_6dof.md (recoupled model,
    2026-08-18 sign correction).

    The screw-nut interface (k_nut, c_nut) is reflected onto theta_s through
    the lead ratio L/2*pi and enters K and C directly, coupling the
    rotational block (theta_m, theta_c, theta_s, theta_sb) to the axial
    block (x_s, x_n). T_fric,nut is no longer an input. The cross-term sign
    is fixed so positive theta_s drives positive x_n (Sec. 1 convention);
    F_nut = k_nut*(x_n - x_s - (L/2pi)*theta_s) + damping (Sec. 3.4).

    C[0,0] includes c_EM (Sec. 5.4): plain viscous damping of theta_m against
    a fixed frame, applied by default because mode 1 (176.7 Hz) is otherwise
    almost undamped (c_c, c_s1, c_s2 barely engage it). No theta_dot_cmd
    feedforward column is added to B_u -- see Sec. 5.4 for why.
    """
    k_EM = p["N_r"] * p["T_hold"]
    k_d = 4.0 * p["N_r"] * p["T_d"]
    lead_ratio = p["L"] / (2.0 * np.pi)

    M = np.diag([p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"]])

    k_c, k_s1, k_s2 = p["k_c"], p["k_s1"], p["k_s2"]
    k_nut, k_brg = p["k_nut"], p["k_brg"]
    K = np.array([
        [k_c + k_EM + k_d, -k_c, 0.0, 0.0, 0.0, 0.0],
        [-k_c, k_c + k_s1, -k_s1, 0.0, 0.0, 0.0],
        [0.0, -k_s1, k_s1 + k_s2 + lead_ratio**2 * k_nut, -k_s2, lead_ratio * k_nut, -lead_ratio * k_nut],
        [0.0, 0.0, -k_s2, k_s2, 0.0, 0.0],
        [0.0, 0.0, lead_ratio * k_nut, 0.0, k_brg + k_nut, -k_nut],
        [0.0, 0.0, -lead_ratio * k_nut, 0.0, -k_nut, k_nut],
    ])

    c_c, c_s1, c_s2 = p["c_c"], p["c_s1"], p["c_s2"]
    c_nut, c_brg, c_EM = p["c_nut"], p["c_brg"], p["c_EM"]
    C = np.array([
        [c_c + c_EM, -c_c, 0.0, 0.0, 0.0, 0.0],
        [-c_c, c_c + c_s1, -c_s1, 0.0, 0.0, 0.0],
        [0.0, -c_s1, c_s1 + c_s2 + lead_ratio**2 * c_nut, -c_s2, lead_ratio * c_nut, -lead_ratio * c_nut],
        [0.0, 0.0, -c_s2, c_s2, 0.0, 0.0],
        [0.0, 0.0, lead_ratio * c_nut, 0.0, c_brg + c_nut, -c_nut],
        [0.0, 0.0, -lead_ratio * c_nut, 0.0, -c_nut, c_nut],
    ])

    B_u = np.array([
        [k_EM + k_d, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
    ])

    return M, C, K, B_u


def build_state_space(M: np.ndarray, C: np.ndarray, K: np.ndarray, B_u: np.ndarray):
    """First-order form per Sec. 6: A, B, and the x_n output row C_y."""
    n = M.shape[0]
    M_inv = np.diag(1.0 / np.diag(M))

    A = np.zeros((2 * n, 2 * n))
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -M_inv @ K
    A[n:, n:] = -M_inv @ C

    B = np.zeros((2 * n, B_u.shape[1]))
    B[n:, :] = M_inv @ B_u

    C_y = np.zeros((1, 2 * n))
    C_y[0, STATE_LABELS.index("x_n")] = 1.0

    return A, B, C_y


def laplace_transfer_function(A: np.ndarray, b_col: np.ndarray, C_y: np.ndarray, s: complex) -> complex:
    """G(s) = C_y (sI - A)^-1 b_col, from sZ(s) = A Z(s) + b_col U(s)."""
    n = A.shape[0]
    z = np.linalg.solve(s * np.eye(n) - A, b_col)
    return complex((C_y @ z)[0])


def main() -> None:
    params = load_parameters()
    M, C, K, B_u = build_matrices(params)
    A, B, C_y = build_state_space(M, C, K, B_u)
    b_col = B[:, INPUT_LABELS.index("theta_cmd")]

    frequencies_hz = np.linspace(0.0, 2000.0, 20001)
    response = np.array([
        laplace_transfer_function(A, b_col, C_y, 1j * 2.0 * np.pi * f)
        for f in frequencies_hz
    ])

    magnitude = np.abs(response)
    magnitude_db = 20.0 * np.log10(np.maximum(magnitude, 1e-15))
    phase_deg = np.unwrap(np.angle(response)) * 180.0 / np.pi

    ASSET_DIR.mkdir(exist_ok=True)
    npz_dir = ASSET_DIR / "npz"
    npz_dir.mkdir(exist_ok=True)
    np.savez(
        npz_dir / "bode_rev4_data.npz",
        frequencies_hz=frequencies_hz,
        magnitude=magnitude,
        magnitude_db=magnitude_db,
        phase_deg=phase_deg,
    )

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True)

    ax_mag.plot(frequencies_hz, magnitude_db, color="#2b6cb0", linewidth=1.4)
    ax_mag.axhline(0.0, color="#999999", linewidth=0.8, zorder=0)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title(r"Rev 4: command-to-stage Bode, $x_n(s) / \theta_{cmd}(s)$")
    ax_mag.grid(True, which="both", linewidth=0.4, color="#cccccc")

    ax_phase.plot(frequencies_hz, phase_deg, color="#c05621", linewidth=1.4)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_xlim(0.0, 2000.0)
    ax_phase.grid(True, which="both", linewidth=0.4, color="#cccccc")

    fig.tight_layout()

    out_path = ASSET_DIR / "bode_rev4.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    print(f"Wrote {out_path}")
    print(f"Wrote {npz_dir / 'bode_rev4_data.npz'}")


if __name__ == "__main__":
    main()

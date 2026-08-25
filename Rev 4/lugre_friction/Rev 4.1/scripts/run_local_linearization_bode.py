#!/usr/bin/env python3
"""Dynamic Local Linearization of the LuGre friction ports, then a proper
Laplace-domain Bode plot of the resulting linear "extended" state space.

Procedure:
  1. Pick a single stage cruising speed V_STAGE (m/s) and derive the
     consistent operating velocity v0 at each of the three interfaces:
       v0_way = V_STAGE                          (guideway sees stage speed)
       v0_sb  = V_STAGE * 2*pi/L                  (screw/bearing rotation
                                                    rate for ideal no-slip
                                                    tracking at V_STAGE)
       v0_nut = 0                                  (v_nut is a SLIP/tracking-
                                                    error term, not a bulk
                                                    speed -- at an idealized
                                                    no-slip steady cruise it
                                                    is exactly zero, not
                                                    V_STAGE; see README.md)
  2. Freeze each bristle at its steady-state deflection z0 = z_ss(v0), then
     take the Jacobian of F(z,v) = sigma0*z*(1-sigma1*|v|/g(v)) + v*(sigma1+
     sigma2) [F with z_dot already substituted in] at (z0, v0):
       K_eq = dF/dz |(z0,v0)   C_eq = dF/dv |(z0,v0)
     Computed numerically (finite difference) rather than by hand-derived
     formula, so this is correct regardless of the sigma1 value. Because
     sigma1=0 at all three interfaces in this sub-branch's parameter set,
     F(z,v) is already exactly linear, so K_eq=sigma0 and C_eq=sigma2
     independent of v0 -- disclosed explicitly below rather than silently
     making the cruising-speed choice look consequential when it isn't
     for this specific parameterization.
  3. Drop K_eq/C_eq into the 6x6 K/C matrices, replacing the nonlinear
     LuGre ports entirely:
       - nut:  same (L/2*pi)-lever-arm cross-coupling structure as the
               baseline's (corrected-sign) k_nut/c_nut recoupling
       - support bearing: NEW grounding term on theta_sb (K_eq_sb, C_eq_sb)
       - guideway: NEW grounding term on x_n (K_eq_way, C_eq_way)
     T_fric,sb and F_fric,way are no longer separate inputs -- they are now
     linear feedback through K/C, exactly like the nut. B_u shrinks to a
     single column (theta_cmd only): this is a genuine single-input LTI
     system, solved the same way as build_bode_rev4.py.
  4. Sweep a dense linear 0-2000 Hz grid and plot the resulting Bode
     (magnitude dB, phase deg), matching the Hz-axis convention used
     throughout this sub-branch.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lugre_model import load_parameters

ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = ROOT / "rendered_assets"

STATE_LABELS = ["theta_m", "theta_c", "theta_s", "theta_sb", "x_s", "x_n"]

V_STAGE = 5.0e-3   # m/s -- chosen representative cruising speed (placeholder)


def lugre_output(z: float, v: float, sigma0: float, sigma1: float, sigma2: float,
                  Fc: float, Fs: float, vs: float) -> float:
    g = Fc + (Fs - Fc) * np.exp(-(v / vs) ** 2)
    zdot = v - sigma0 * abs(v) / g * z
    return sigma0 * z + sigma1 * zdot + sigma2 * v


def steady_z(v0: float, sigma0: float, Fc: float, Fs: float, vs: float) -> float:
    if v0 == 0.0:
        return 0.0
    g0 = Fc + (Fs - Fc) * np.exp(-(v0 / vs) ** 2)
    return np.sign(v0) * g0 / sigma0


def equivalent_stiffness_damping(v0: float, sigma0: float, sigma1: float, sigma2: float,
                                  Fc: float, Fs: float, vs: float):
    """Numerical Jacobian of F(z,v) at the frozen operating point (z0, v0)."""
    z0 = steady_z(v0, sigma0, Fc, Fs, vs)
    dz = max(abs(z0), 1e-9) * 1e-6
    dv = max(abs(v0), vs) * 1e-6
    args = (sigma0, sigma1, sigma2, Fc, Fs, vs)
    K_eq = (lugre_output(z0 + dz, v0, *args) - lugre_output(z0 - dz, v0, *args)) / (2 * dz)
    C_eq = (lugre_output(z0, v0 + dv, *args) - lugre_output(z0, v0 - dv, *args)) / (2 * dv)
    F0 = lugre_output(z0, v0, *args)
    return K_eq, C_eq, z0, F0


def build_linearized_matrices(p: dict[str, float]):
    k_EM = p["N_r"] * p["T_hold"]
    k_d = 4.0 * p["N_r"] * p["T_d"]
    lead_ratio = p["L"] / (2.0 * np.pi)

    omega0 = V_STAGE * 2.0 * np.pi / p["L"]
    v0_way = V_STAGE
    v0_sb = omega0
    v0_nut = 0.0

    K_eq_nut, C_eq_nut, z0_nut, F0_nut = equivalent_stiffness_damping(
        v0_nut, p["sigma0_nut"], p["sigma1_nut"], p["sigma2_nut"], p["Fc_nut"], p["Fs_nut"], p["vs_nut"])
    K_eq_sb, C_eq_sb, z0_sb, F0_sb = equivalent_stiffness_damping(
        v0_sb, p["sigma0_sb"], p["sigma1_sb"], p["sigma2_sb"], p["Tc_sb"], p["Ts_sb"], p["vs_sb"])
    K_eq_way, C_eq_way, z0_way, F0_way = equivalent_stiffness_damping(
        v0_way, p["sigma0_way"], p["sigma1_way"], p["sigma2_way"], p["Fc_way"], p["Fs_way"], p["vs_way"])

    print(f"V_STAGE = {V_STAGE*1e3:.2f} mm/s  ->  omega0 (screw/bearing) = {omega0:.4f} rad/s")
    print(f"  nut : v0={v0_nut:.6e} m/s   z0={z0_nut:.4e}  F0={F0_nut:.4e} N   "
          f"K_eq={K_eq_nut:.4e} N/m   C_eq={C_eq_nut:.4e} N*s/m")
    print(f"  sb  : v0={v0_sb:.6e} rad/s  z0={z0_sb:.4e}  F0={F0_sb:.4e} N*m  "
          f"K_eq={K_eq_sb:.4e} N*m/rad C_eq={C_eq_sb:.4e} N*m*s/rad")
    print(f"  way : v0={v0_way:.6e} m/s   z0={z0_way:.4e}  F0={F0_way:.4e} N   "
          f"K_eq={K_eq_way:.4e} N/m   C_eq={C_eq_way:.4e} N*s/m")

    M = np.diag([p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"]])

    k_c, k_s1, k_s2, k_brg = p["k_c"], p["k_s1"], p["k_s2"], p["k_brg"]
    K = np.array([
        [k_c + k_EM + k_d, -k_c, 0.0, 0.0, 0.0, 0.0],
        [-k_c, k_c + k_s1, -k_s1, 0.0, 0.0, 0.0],
        [0.0, -k_s1, k_s1 + k_s2 + lead_ratio**2 * K_eq_nut, -k_s2, lead_ratio * K_eq_nut, -lead_ratio * K_eq_nut],
        [0.0, 0.0, -k_s2, k_s2 + K_eq_sb, 0.0, 0.0],
        [0.0, 0.0, lead_ratio * K_eq_nut, 0.0, k_brg + K_eq_nut, -K_eq_nut],
        [0.0, 0.0, -lead_ratio * K_eq_nut, 0.0, -K_eq_nut, K_eq_nut + K_eq_way],
    ])

    c_c, c_s1, c_s2, c_brg, c_EM = p["c_c"], p["c_s1"], p["c_s2"], p["c_brg"], p["c_EM"]
    C = np.array([
        [c_c + c_EM, -c_c, 0.0, 0.0, 0.0, 0.0],
        [-c_c, c_c + c_s1, -c_s1, 0.0, 0.0, 0.0],
        [0.0, -c_s1, c_s1 + c_s2 + lead_ratio**2 * C_eq_nut, -c_s2, lead_ratio * C_eq_nut, -lead_ratio * C_eq_nut],
        [0.0, 0.0, -c_s2, c_s2 + C_eq_sb, 0.0, 0.0],
        [0.0, 0.0, lead_ratio * C_eq_nut, 0.0, c_brg + C_eq_nut, -C_eq_nut],
        [0.0, 0.0, -lead_ratio * C_eq_nut, 0.0, -C_eq_nut, C_eq_nut + C_eq_way],
    ])

    B_em = np.array([k_EM, 0.0, 0.0, 0.0, 0.0, 0.0])  # single input: theta_cmd only
    return M, K, C, B_em


def build_state_space(M, K, C, B_em):
    n = M.shape[0]
    M_inv = np.diag(1.0 / np.diag(M))
    A = np.zeros((2 * n, 2 * n))
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -M_inv @ K
    A[n:, n:] = -M_inv @ C
    B = np.zeros(2 * n)
    B[n:] = M_inv @ B_em
    C_y = np.zeros((1, 2 * n))
    C_y[0, STATE_LABELS.index("x_n")] = 1.0
    return A, B, C_y


def main():
    p = load_parameters()
    M, K, C, B_em = build_linearized_matrices(p)
    A, B, Cy = build_state_space(M, K, C, B_em)

    freq_hz = np.linspace(0.0, 2000.0, 20001)
    freq_hz[0] = 0.0  # DC is fine here: matrix inversion, no division by period
    n = A.shape[0]
    response = np.empty(len(freq_hz), dtype=complex)
    for i, f in enumerate(freq_hz):
        s = 1j * 2.0 * np.pi * f
        z = np.linalg.solve(s * np.eye(n) - A, B)
        response[i] = (Cy @ z)[0]

    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(response), 1e-300))
    phase_deg = np.unwrap(np.angle(response)) * 180.0 / np.pi

    ASSET_DIR.mkdir(exist_ok=True)
    npz_dir = ASSET_DIR / "npz"
    npz_dir.mkdir(exist_ok=True)
    np.savez(npz_dir / "local_linearization_bode_data.npz",
             freq_hz=freq_hz, magnitude_db=magnitude_db, phase_deg=phase_deg)

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True)
    ax_mag.plot(freq_hz, magnitude_db, color="#2b6cb0", linewidth=1.2)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title(r"Rev 4 LuGre sub-branch: dynamic local linearization, "
                      rf"$x_n(s)/\theta_{{cmd}}(s)$ @ V_stage={V_STAGE*1e3:.1f} mm/s")
    ax_mag.grid(True, linewidth=0.4, color="#cccccc")

    ax_phase.plot(freq_hz, phase_deg, color="#c05621", linewidth=1.2)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_xlim(0.0, 2000.0)
    ax_phase.grid(True, linewidth=0.4, color="#cccccc")

    fig.tight_layout()
    out_path = ASSET_DIR / "local_linearization_bode.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()

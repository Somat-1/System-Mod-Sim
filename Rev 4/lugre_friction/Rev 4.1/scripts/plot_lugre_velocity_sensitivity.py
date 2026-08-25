#!/usr/bin/env python3
"""Does the frozen-linearization reference velocity (V_STAGE) matter?

Requested check: run_local_linearization_bode.py freezes the LuGre bristle
at a single operating point (V_STAGE = 5 mm/s) to get K_eq/C_eq for the
Co-MAC Set B matrices. This script tests whether a different V_STAGE choice
would have given different K_eq/C_eq (and therefore a different Bode plot
and different Co-MAC numbers), by sweeping v0 directly through
equivalent_stiffness_damping() for all three ports, and by re-running the
full Bode computation at several different V_STAGE values.

For this parameter set sigma1_nut = sigma1_sb = sigma1_way = 0.0 exactly
(confirmed below, not assumed), which makes the LuGre force law
F(z, v) = sigma0*z + sigma2*v -- already perfectly linear in (z, v), with no
term that couples the tangent stiffness/damping to the operating point v0 at
all. So K_eq = dF/dz = sigma0 and C_eq = dF/dv = sigma2 are mathematically
invariant to v0 for this parameterization; what actually varies with v0 is
the frozen bristle deflection z0 itself (the Stribeck steady-state curve),
not the tangent stiffness evaluated there. This is disclosed in
run_local_linearization_bode.py's own module docstring; this script verifies
it numerically (sweep + a from-scratch Bode re-run at multiple V_STAGE
values) rather than just re-quoting that claim.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lugre_model import load_parameters
from run_local_linearization_bode import (
    STATE_LABELS, equivalent_stiffness_damping, V_STAGE as V_STAGE_REF,
)

LUGRE_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp"

PORTS = [
    ("nut", "sigma0_nut", "sigma1_nut", "sigma2_nut", "Fc_nut", "Fs_nut", "vs_nut"),
    ("support bearing", "sigma0_sb", "sigma1_sb", "sigma2_sb", "Tc_sb", "Ts_sb", "vs_sb"),
    ("guideway", "sigma0_way", "sigma1_way", "sigma2_way", "Fc_way", "Fs_way", "vs_way"),
]
PORT_COLORS = {"nut": "#2b6cb0", "support bearing": "#c05621", "guideway": "#2f855a"}


def build_matrices_at_v_stage(p: dict, v_stage: float):
    """Same assembly as run_local_linearization_bode.build_linearized_matrices,
    parameterized by v_stage explicitly instead of reading the module-level
    V_STAGE constant, so it can be swept here without touching that module."""
    k_EM = p["N_r"] * p["T_hold"]
    k_d = 4.0 * p["N_r"] * p["T_d"]
    lead_ratio = p["L"] / (2.0 * np.pi)

    omega0 = v_stage * 2.0 * np.pi / p["L"]
    v0_way = v_stage
    v0_sb = omega0
    v0_nut = 0.0

    K_eq_nut, C_eq_nut, _, _ = equivalent_stiffness_damping(
        v0_nut, p["sigma0_nut"], p["sigma1_nut"], p["sigma2_nut"], p["Fc_nut"], p["Fs_nut"], p["vs_nut"])
    K_eq_sb, C_eq_sb, _, _ = equivalent_stiffness_damping(
        v0_sb, p["sigma0_sb"], p["sigma1_sb"], p["sigma2_sb"], p["Tc_sb"], p["Ts_sb"], p["vs_sb"])
    K_eq_way, C_eq_way, _, _ = equivalent_stiffness_damping(
        v0_way, p["sigma0_way"], p["sigma1_way"], p["sigma2_way"], p["Fc_way"], p["Fs_way"], p["vs_way"])

    M = np.diag([p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"]])
    k_c, k_s1, k_s2, k_brg = p["k_c"], p["k_s1"], p["k_s2"], p["k_brg"]
    K = np.array([
        [k_c + k_EM + 4.0 * p["N_r"] * p["T_d"], -k_c, 0.0, 0.0, 0.0, 0.0],
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
    B_em = np.array([k_EM, 0.0, 0.0, 0.0, 0.0, 0.0])
    return M, K, C, B_em


def bode(M, K, C, B_em, freq_hz):
    n = M.shape[0]
    M_inv = np.diag(1.0 / np.diag(M))
    A = np.zeros((2 * n, 2 * n))
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -M_inv @ K
    A[n:, n:] = -M_inv @ C
    B = np.zeros(2 * n)
    B[n:] = M_inv @ B_em
    Cy = np.zeros((1, 2 * n))
    Cy[0, STATE_LABELS.index("x_n")] = 1.0

    response = np.empty(len(freq_hz), dtype=complex)
    for i, f in enumerate(freq_hz):
        s = 1j * 2.0 * np.pi * f
        z = np.linalg.solve(s * np.eye(2 * n) - A, B)
        response[i] = (Cy @ z)[0]
    return response


def main() -> None:
    p = load_parameters()

    for name, s0k, s1k, *_ in PORTS:
        print(f"{name}: sigma0={p[s0k]:.3e}  sigma1={p[s1k]:.3e}")

    # ---- Panel A/B: K_eq (normalized) and z0 vs v0, swept independently of V_STAGE ----
    v0_sweep = np.logspace(-6, 0, 200)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
    ax_keq, ax_z0, ax_bode, ax_resid = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    v_stage_markers = {
        "nut": 0.0,
        "support bearing": V_STAGE_REF * 2.0 * np.pi / p["L"],
        "guideway": V_STAGE_REF,
    }

    for name, s0k, s1k, s2k, ak, bk, vsk in PORTS:
        color = PORT_COLORS[name]
        K_eq_arr = np.array([equivalent_stiffness_damping(
            v0, p[s0k], p[s1k], p[s2k], p[ak], p[bk], p[vsk])[0] for v0 in v0_sweep])
        z0_arr = np.array([equivalent_stiffness_damping(
            v0, p[s0k], p[s1k], p[s2k], p[ak], p[bk], p[vsk])[2] for v0 in v0_sweep])

        ax_keq.plot(v0_sweep, K_eq_arr / p[s0k], color=color, linewidth=1.6, label=name)
        ax_z0.plot(v0_sweep, z0_arr * 1e6, color=color, linewidth=1.6, label=name)

        v0_marker = v_stage_markers[name]
        if v0_marker > 0:
            K_eq_m, _, z0_m, _ = equivalent_stiffness_damping(
                v0_marker, p[s0k], p[s1k], p[s2k], p[ak], p[bk], p[vsk])
            ax_keq.plot(v0_marker, K_eq_m / p[s0k], "o", color=color, markersize=7, zorder=5)
            ax_z0.plot(v0_marker, z0_m * 1e6, "o", color=color, markersize=7, zorder=5)

    ax_keq.set_xscale("log")
    ax_keq.set_ylim(0.0, 1.5)
    ax_keq.axhline(1.0, color="#999999", linewidth=0.8, linestyle=":", zorder=0)
    ax_keq.set_xlabel("v0 (m/s or rad/s)")
    ax_keq.set_ylabel("K_eq / sigma0 (dimensionless)")
    ax_keq.set_title("Equivalent tangent stiffness vs. frozen operating velocity")
    ax_keq.legend(fontsize=8)
    ax_keq.grid(True, which="both", linewidth=0.4, color="#cccccc")

    ax_z0.set_xscale("log")
    ax_z0.set_xlabel("v0 (m/s or rad/s)")
    ax_z0.set_ylabel("z0, frozen bristle deflection (um or urad)")
    ax_z0.set_title("What actually moves with v0: the Stribeck deflection point")
    ax_z0.legend(fontsize=8)
    ax_z0.grid(True, which="both", linewidth=0.4, color="#cccccc")

    # ---- Panel C/D: Bode at several V_STAGE choices, and the residual vs. the
    # 5 mm/s reference used throughout the Co-MAC derivations ----
    v_stage_choices_mm_s = [0.05, 0.5, 5.0, 50.0, 500.0]
    freq_hz = np.linspace(0.0, 2000.0, 4001)
    freq_hz[0] = 1e-3

    responses = {}
    for v_mm_s in v_stage_choices_mm_s:
        M, K, C, B_em = build_matrices_at_v_stage(p, v_mm_s * 1e-3)
        responses[v_mm_s] = bode(M, K, C, B_em, freq_hz)

    cmap = plt.get_cmap("viridis")
    for i, v_mm_s in enumerate(v_stage_choices_mm_s):
        mag_db = 20.0 * np.log10(np.maximum(np.abs(responses[v_mm_s]), 1e-300))
        ax_bode.plot(freq_hz, mag_db, color=cmap(i / (len(v_stage_choices_mm_s) - 1)),
                     linewidth=1.4, label=f"V_stage = {v_mm_s:g} mm/s")
    ax_bode.set_xlabel("Frequency (Hz)")
    ax_bode.set_ylabel("Magnitude (dB)")
    ax_bode.set_title(r"$x_n(s)/\theta_{cmd}(s)$ at different frozen reference velocities")
    ax_bode.legend(fontsize=8)
    ax_bode.grid(True, linewidth=0.4, color="#cccccc")

    ref_mag_db = 20.0 * np.log10(np.maximum(np.abs(responses[5.0]), 1e-300))
    for i, v_mm_s in enumerate(v_stage_choices_mm_s):
        mag_db = 20.0 * np.log10(np.maximum(np.abs(responses[v_mm_s]), 1e-300))
        resid = np.abs(mag_db - ref_mag_db)
        ax_resid.plot(freq_hz, np.maximum(resid, 1e-16),
                      color=cmap(i / (len(v_stage_choices_mm_s) - 1)), linewidth=1.2,
                      label=f"{v_mm_s:g} mm/s vs. 5 mm/s")
    ax_resid.set_yscale("log")
    ax_resid.set_xlabel("Frequency (Hz)")
    ax_resid.set_ylabel("|magnitude_dB - magnitude_dB @ 5 mm/s|")
    ax_resid.set_title("Residual vs. the Co-MAC reference velocity (floating-point noise floor)")
    ax_resid.legend(fontsize=8)
    ax_resid.grid(True, which="both", linewidth=0.4, color="#cccccc")

    fig.suptitle(
        "LuGre local-linearization: does the frozen reference velocity (V_STAGE) matter?\n"
        "sigma1 = 0 at all three ports for this parameter set -> K_eq/C_eq are exactly "
        "velocity-independent, so no -- verified both directly (left) and via a from-scratch "
        "Bode re-run (right)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.90])

    out_path = OUT_DIR / "lugre_velocity_sensitivity.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    print(f"\nMax residual across all V_STAGE choices and all frequencies: "
          f"{max(np.max(np.abs(20*np.log10(np.maximum(np.abs(responses[v]),1e-300)) - ref_mag_db)) for v in v_stage_choices_mm_s):.3e} dB")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()

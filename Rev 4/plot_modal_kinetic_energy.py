#!/usr/bin/env python3
"""Modal kinetic energy decomposition of the frictionless baseline, requested
2026-08-20 -- the exact mirror of plot_modal_strain_energy.py: mode shape
says what moves, this says which inertia carries it, strain energy says
which spring resists it.

KE_frac[i,j] = M_eq[i] * phi_scaled[i,j]^2, using the SAME phi_scaled
(rotational DOFs converted to equivalent axial displacement by L/(2*pi))
already computed in plot_mode_shapes.py, paired with an equivalent mass
M_eq so the identity carries through exactly:

    M_eq[i] = I_ii / (L/2*pi)^2   for the four rotational DOFs (same
                                    inertia-reflection used for I_eff in
                                    the broadband-ID torque-limit check)
    M_eq[i] = M_ii                 for x_s, x_n (already a mass)

so that M_eq[i]*phi_scaled[i,j]^2 = M_ii*phi[i,j]^2 exactly (the unscaled,
raw mass-normalized identity), and therefore
    sum_i KE_frac[i,j] = phi_j.T @ M @ phi_j = 1
automatically, by mass-normalization alone -- same "MANDATORY CHECK true
by construction" property the strain-energy fractions have, just verified
numerically below rather than asserted structurally.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from build_bode_rev4 import STATE_LABELS, build_matrices, load_parameters
from plot_colors import DOF_COLORS

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "rendered_assets" / "temp"

ROTATIONAL_DOFS = {"theta_m", "theta_c", "theta_s", "theta_sb"}


def main() -> None:
    params = load_parameters()
    M, C, K, B_u = build_matrices(params)
    n = M.shape[0]
    lead_ratio = params["L"] / (2.0 * np.pi)

    lam, phi = eigh(K, M)
    omega = np.sqrt(lam)
    freq_hz = omega / (2.0 * np.pi)

    scale = np.array([lead_ratio if lbl in ROTATIONAL_DOFS else 1.0 for lbl in STATE_LABELS])
    phi_scaled = phi * scale[:, None]
    for j in range(n):
        idx_max = np.argmax(np.abs(phi_scaled[:, j]))
        if phi_scaled[idx_max, j] < 0:
            phi_scaled[:, j] *= -1.0

    M_diag = np.diag(M)
    M_eq = np.array([
        M_diag[i] / scale[i] ** 2 if lbl in ROTATIONAL_DOFS else M_diag[i]
        for i, lbl in enumerate(STATE_LABELS)
    ])
    print("Equivalent mass M_eq (kg for x_s/x_n; I_ii/(L/2pi)^2, reflected to kg, for rotational DOFs):")
    for lbl, m in zip(STATE_LABELS, M_eq):
        print(f"  {lbl:10s} {m:.4e}")

    KE = M_eq[:, None] * phi_scaled ** 2   # (n_dof, n_modes)
    col_sums = KE.sum(axis=0)
    print(f"\n{'mode':>4s} {'f (Hz)':>10s}  " + "  ".join(f"{lbl:>10s}" for lbl in STATE_LABELS) + f"  {'sum':>8s}")
    for j in range(n):
        vals = "  ".join(f"{KE[i, j]:10.4f}" for i in range(n))
        print(f"{j+1:4d} {freq_hz[j]:10.2f}  {vals}  {col_sums[j]:8.4f}")
    assert np.allclose(col_sums, 1.0, atol=1e-6), "Kinetic-energy fractions do not sum to 1 per mode."

    fig, ax = plt.subplots(figsize=(10.0, 6.5))
    x_pos = np.arange(n)
    bottom = np.zeros(n)
    for i, lbl in enumerate(STATE_LABELS):
        ax.bar(x_pos, KE[i], bottom=bottom, color=DOF_COLORS[i], label=lbl, width=0.6)
        bottom += KE[i]

    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"mode {j+1}\n{freq_hz[j]:.0f} Hz" for j in range(n)], fontsize=8)
    ax.set_ylabel("Fraction of modal kinetic energy")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Modal kinetic energy decomposition, frictionless baseline\n"
                 r"KE$[i,j] = M_{eq,i} \cdot \phi_{scaled,i,j}^2$, stacked per mode (sums to 1)")
    ax.grid(True, axis="y", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=8, ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.12))

    fig.tight_layout()
    out_path = OUT_DIR / "modal_kinetic_energy.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "modal_kinetic_energy_data.npz",
              freq_hz=freq_hz, KE=KE, M_eq=M_eq, state_labels=np.array(STATE_LABELS))

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'modal_kinetic_energy_data.npz'}")


if __name__ == "__main__":
    main()

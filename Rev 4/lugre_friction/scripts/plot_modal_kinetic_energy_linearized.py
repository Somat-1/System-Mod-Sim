#!/usr/bin/env python3
"""Modal kinetic energy decomposition of the frozen-linearized (LuGre local-
linearization) system, requested 2026-08-21 -- mirror of
../../scripts/plot_modal_kinetic_energy.py (itself the mirror of
plot_modal_strain_energy.py): mode shape says what moves, strain energy says
which spring resists it, this says which inertia carries it.

KE_frac[i,j] = M_eq[i] * phi_scaled[i,j]^2, same M_eq convention as the
frictionless version (M_eq[i] = I_ii/(L/2*pi)^2 for the four rotational
DOFs, M_eq[i] = M_ii for x_s/x_n) -- friction ports add stiffness/damping
only, never mass, so M itself and this M_eq are IDENTICAL to the
frictionless baseline's.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from lugre_model import load_parameters
from run_local_linearization_bode import STATE_LABELS, build_linearized_matrices

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parent
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from plot_colors import DOF_COLORS  # noqa: E402

ROTATIONAL_DOFS = {"theta_m", "theta_c", "theta_s", "theta_sb"}


def main() -> None:
    params = load_parameters()
    M, K, C, B_em = build_linearized_matrices(params)
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

    KE = M_eq[:, None] * phi_scaled ** 2
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
    ax.set_title("Modal kinetic energy decomposition, frozen-linearized LuGre system\n"
                 r"KE$[i,j] = M_{eq,i} \cdot \phi_{scaled,i,j}^2$, stacked per mode (sums to 1)")
    ax.grid(True, axis="y", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=8, ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.12))

    fig.tight_layout()
    out_path = OUT_DIR / "modal_kinetic_energy_linearized.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "modal_kinetic_energy_linearized_data.npz",
              freq_hz=freq_hz, KE=KE, M_eq=M_eq, state_labels=np.array(STATE_LABELS))

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'modal_kinetic_energy_linearized_data.npz'}")


if __name__ == "__main__":
    main()

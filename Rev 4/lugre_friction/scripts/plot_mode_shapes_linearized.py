#!/usr/bin/env python3
"""Mode shapes of the frozen-linearized (LuGre local-linearization) system,
requested 2026-08-21 -- exact mirror of ../../scripts/plot_mode_shapes.py,
applied to run_local_linearization_bode.py's K/C/M (frozen bristle
equivalent stiffness/damping dropped into the 6x6 matrices at
V_STAGE=5 mm/s) instead of the frictionless baseline's.

scipy.linalg.eigh(K, M) -- same generalized symmetric solver, same
mass-normalization (phi.T @ M @ phi = I) and units caveat (m/sqrt(kg), not
m) as the frictionless version. Sign convention and per-mode normalization
+ annotation policy are identical, for direct visual comparability against
the frictionless mode_shapes.png.
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
    phi_scaled = phi * scale[:, None]   # equivalent axial displacement, every row, m/sqrt(kg)

    for j in range(n):
        idx_max = np.argmax(np.abs(phi_scaled[:, j]))
        if phi_scaled[idx_max, j] < 0:
            phi_scaled[:, j] *= -1.0
            phi[:, j] *= -1.0

    print(f"{'mode':>4s} {'f (Hz)':>10s}  " + "  ".join(f"{lbl:>12s}" for lbl in STATE_LABELS))
    for j in range(n):
        vals = "  ".join(f"{phi_scaled[i, j]:12.4e}" for i in range(n))
        print(f"{j+1:4d} {freq_hz[j]:10.2f}  {vals}")

    phi_norm = phi_scaled / np.max(np.abs(phi_scaled), axis=0, keepdims=True)

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0))
    x_pos = np.arange(n)
    for j, ax in enumerate(axes.flat):
        bars = ax.bar(x_pos, phi_norm[:, j], color=DOF_COLORS)
        ax.axhline(0.0, color="#333333", linewidth=0.8)
        for i, b in enumerate(bars):
            val = phi_scaled[i, j]
            y = b.get_height()
            va = "bottom" if y >= 0 else "top"
            offset = 0.03 if y >= 0 else -0.03
            ax.text(b.get_x() + b.get_width() / 2.0, y + offset, f"{val:.2e}",
                    ha="center", va=va, fontsize=6, rotation=90, clip_on=False)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(STATE_LABELS, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"Mode {j+1}: {freq_hz[j]:.1f} Hz", fontsize=10)
        ax.set_ylabel("Normalized shape (max |.| = 1)", fontsize=7.5)
        ax.set_ylim(-1.35, 1.35)
        ax.grid(True, axis="y", linewidth=0.4, color="#cccccc")

    fig.suptitle(
        "Mode shapes, frozen-linearized LuGre system (mass-normalized, sign-fixed)\n"
        "normalized to max |.| = 1 per mode; rotational DOFs scaled by L/(2π) to equivalent axial "
        "displacement, units m/√kg\n"
        "numbers above bars are the actual (unnormalized) values",
        fontsize=9.5,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.87])

    out_path = OUT_DIR / "mode_shapes_linearized.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "mode_shapes_linearized_data.npz",
              freq_hz=freq_hz, omega=omega, phi=phi, phi_scaled=phi_scaled, phi_norm=phi_norm,
              state_labels=np.array(STATE_LABELS), lead_ratio=lead_ratio)

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'mode_shapes_linearized_data.npz'}")


if __name__ == "__main__":
    main()

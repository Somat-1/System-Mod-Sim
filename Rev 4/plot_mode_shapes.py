#!/usr/bin/env python3
"""Mode shapes of the frictionless baseline, requested 2026-08-20; revised
2026-08-20 per review: correct units, per-mode normalization, and numeric
annotations, after the linear scale was found to hide mode 2's motor
participation (see below).

scipy.linalg.eigh(K, M) -- the generalized SYMMETRIC solver (M diagonal
positive-definite, K symmetric), not eig(A) on the 12-state first-order
system: it returns real, mass-normalized eigenvectors (phi.T @ M @ phi = I)
directly and avoids the state_space_6dof.md Sec. 9 item 3 conditioning
problem (rotational vs. translational entries differ by ~2*pi/L ~ 6.28e3).

Unit scaling. The four rotational DOFs (theta_m, theta_c, theta_s,
theta_sb) are in rad; the two translational DOFs (x_s, x_n) are in m.
Plotted raw, the translational entries are invisible next to the
rotational ones. Every mode shape below is converted to "equivalent axial
displacement" by multiplying the rotational components by L/(2*pi). Units
are m/sqrt(kg), not m: mass-normalization (phi.T @ M @ phi = I) makes phi
itself a mass-normalized shape, not a physical displacement -- physical
displacement is recovered as phi times a modal coordinate carrying the
matching sqrt(kg) factor. Stated on the figure, since getting this wrong
(labelling it "m") is the difference between a readable plot and a
misleading one.

Sign convention. Mass-normalization fixes magnitude, not sign -- phi_j and
-phi_j are equally valid. Each mode's sign is fixed here by forcing its
largest-magnitude (scaled) component positive, so mode shapes stay
comparable across repeated runs and against the LuGre version later.

Per-mode normalization + annotations (2026-08-20 fix). A linear axis
auto-scaled to each mode's own largest bar hides small-but-dynamically-
critical components: mode 2's scaled theta_m component is 0.0080, which on
mode 2's raw axis (dominated by x_s/x_n ~1.3-1.5) was 0.53% of full scale
-- invisible -- even though mode 2's residue |R|/|R1|=1.130 is LARGER than
mode 1's, and it produces the clear 746 Hz peak in the Bode. Read off the
old plot, mode 2 looked like it had no motor participation and therefore
couldn't be excited by theta_cmd -- backwards. Every mode's bars are now
normalized to max(|.|)=1 (a consistent 0-1 axis everywhere) AND the actual
signed value is printed above each bar, so the small numbers that matter
are legible regardless of bar height.
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
    phi_scaled = phi * scale[:, None]   # equivalent axial displacement, every row, m/sqrt(kg)

    # Sign fix: force each mode's largest-magnitude scaled component positive.
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
        "Mode shapes, frictionless baseline (mass-normalized, sign-fixed)\n"
        "normalized to max |.| = 1 per mode; rotational DOFs scaled by L/(2π) to equivalent axial "
        "displacement, units m/√kg\n"
        "numbers above bars are the actual (unnormalized) values",
        fontsize=9.5,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.87])

    out_path = OUT_DIR / "mode_shapes.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "mode_shapes_data.npz",
              freq_hz=freq_hz, omega=omega, phi=phi, phi_scaled=phi_scaled, phi_norm=phi_norm,
              state_labels=np.array(STATE_LABELS), lead_ratio=lead_ratio)

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'mode_shapes_data.npz'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Modal strain energy decomposition of the frictionless baseline, requested
2026-08-20.

K = sum_e K_e, one term per physical spring, each an outer product of its
own connectivity vector b_e and its stiffness: K_e = k_e * outer(b_e, b_e).
STATE_LABELS order: [theta_m, theta_c, theta_s, theta_sb, x_s, x_n].

  k_EM, k_d   grounded on theta_m           b = e_theta_m
  k_c         theta_m <-> theta_c           b = e_theta_m - e_theta_c
  k_s1        theta_c <-> theta_s           b = e_theta_c - e_theta_s
  k_s2        theta_s <-> theta_sb          b = e_theta_s - e_theta_sb
  k_brg       grounded on x_s               b = e_x_s
  k_nut       x_n - x_s - (L/2*pi)*theta_s  b = [0,0,-L/2*pi,0,-1,+1]
              (backlog issue #1's corrected sign convention -- getting this
              wrong silently corrupts every fraction below)

MANDATORY CHECK: sum(K_e) must equal the assembled K from build_matrices()
to numerical tolerance. This assert is the entire validity condition for
everything downstream of it -- if a term is missing, every fraction is
wrong, so it is not a warning, it is a hard failure.

Fractions: with the mass-normalized modes from plot_mode_shapes.py's
eigenanalysis, SE[e,j] = phi_j.T @ K_e @ phi_j / lambda_j. Each column
(each mode) sums to 1, since sum_e K_e = K and phi_j.T @ K @ phi_j =
lambda_j for mass-normalized phi_j.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from build_bode_rev4 import STATE_LABELS, build_matrices, load_parameters
from plot_colors import ELEMENT_COLORS as _ELEMENT_COLOR_LIST

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "rendered_assets" / "temp"

IDX = {lbl: i for i, lbl in enumerate(STATE_LABELS)}
N = len(STATE_LABELS)

_ELEMENT_ORDER = ["k_EM", "k_d", "k_c", "k_s1", "k_s2", "k_brg", "k_nut"]
ELEMENT_COLORS = dict(zip(_ELEMENT_ORDER, _ELEMENT_COLOR_LIST))


def unit(label: str) -> np.ndarray:
    v = np.zeros(N)
    v[IDX[label]] = 1.0
    return v


def build_element_stiffnesses(p: dict[str, float]) -> dict[str, np.ndarray]:
    k_EM = p["N_r"] * p["T_hold"]
    k_d = 4.0 * p["N_r"] * p["T_d"]
    lead_ratio = p["L"] / (2.0 * np.pi)

    b_EM = unit("theta_m")
    b_d = unit("theta_m")
    b_c = unit("theta_m") - unit("theta_c")
    b_s1 = unit("theta_c") - unit("theta_s")
    b_s2 = unit("theta_s") - unit("theta_sb")
    b_brg = unit("x_s")
    b_nut = np.zeros(N)
    b_nut[IDX["theta_s"]] = -lead_ratio
    b_nut[IDX["x_s"]] = -1.0
    b_nut[IDX["x_n"]] = 1.0

    elements = {
        "k_EM": (k_EM, b_EM), "k_d": (k_d, b_d), "k_c": (p["k_c"], b_c),
        "k_s1": (p["k_s1"], b_s1), "k_s2": (p["k_s2"], b_s2),
        "k_brg": (p["k_brg"], b_brg), "k_nut": (p["k_nut"], b_nut),
    }
    return {name: k * np.outer(b, b) for name, (k, b) in elements.items()}


def main() -> None:
    params = load_parameters()
    M, C, K, B_u = build_matrices(params)

    K_elements = build_element_stiffnesses(params)
    K_sum = sum(K_elements.values())
    max_err = np.max(np.abs(K_sum - K))
    print(f"MANDATORY CHECK: max|sum(K_e) - K| = {max_err:.3e} "
          f"(K scale ~{np.max(np.abs(K)):.3e})")
    assert np.allclose(K_sum, K, atol=1e-6, rtol=1e-9), (
        "Element stiffness decomposition does not reconstruct K -- "
        "a term is missing or a sign is wrong; every fraction below would be invalid."
    )
    print("  -> PASS: element stiffnesses reconstruct K exactly.\n")

    lam, phi = eigh(K, M)
    omega = np.sqrt(lam)
    freq_hz = omega / (2.0 * np.pi)
    n_modes = len(lam)

    element_names = list(K_elements.keys())
    SE = np.zeros((len(element_names), n_modes))
    for e_idx, name in enumerate(element_names):
        K_e = K_elements[name]
        for j in range(n_modes):
            SE[e_idx, j] = (phi[:, j] @ K_e @ phi[:, j]) / lam[j]

    col_sums = SE.sum(axis=0)
    print(f"{'mode':>4s} {'f (Hz)':>10s}  " + "  ".join(f"{nm:>8s}" for nm in element_names) + f"  {'sum':>8s}")
    for j in range(n_modes):
        vals = "  ".join(f"{SE[e_idx, j]:8.4f}" for e_idx in range(len(element_names)))
        print(f"{j+1:4d} {freq_hz[j]:10.2f}  {vals}  {col_sums[j]:8.4f}")
    assert np.allclose(col_sums, 1.0, atol=1e-6), "Strain-energy fractions do not sum to 1 per mode."

    fig, ax = plt.subplots(figsize=(10.0, 6.5))
    x_pos = np.arange(n_modes)
    bottom = np.zeros(n_modes)
    for e_idx, name in enumerate(element_names):
        ax.bar(x_pos, SE[e_idx], bottom=bottom, color=ELEMENT_COLORS[name], label=name, width=0.6)
        bottom += SE[e_idx]

    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"mode {j+1}\n{freq_hz[j]:.0f} Hz" for j in range(n_modes)], fontsize=8)
    ax.set_ylabel("Fraction of modal strain energy")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Modal strain energy decomposition, frictionless baseline\n"
                 r"SE$[e,j] = \phi_j^T K_e \phi_j / \lambda_j$, stacked per mode (sums to 1)")
    ax.grid(True, axis="y", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=8, ncol=len(element_names), loc="upper center", bbox_to_anchor=(0.5, -0.12))

    fig.tight_layout()
    out_path = OUT_DIR / "modal_strain_energy.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "modal_strain_energy_data.npz",
              freq_hz=freq_hz, SE=SE, element_names=np.array(element_names))

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'modal_strain_energy_data.npz'}")


if __name__ == "__main__":
    main()

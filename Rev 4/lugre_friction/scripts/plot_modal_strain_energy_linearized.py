#!/usr/bin/env python3
"""Modal strain energy decomposition of the frozen-linearized (LuGre local-
linearization) system, requested 2026-08-21 -- mirror of
../../scripts/plot_modal_strain_energy.py, with the frictionless baseline's
k_nut element replaced by three LuGre port elements (sigma0_nut, sigma0_sb,
sigma0_way), the frozen equivalent stiffnesses run_local_linearization_
bode.py drops into K at V_STAGE=5 mm/s.

K = sum_e K_e: six structural springs (k_EM, k_d, k_c, k_s1, k_s2, k_brg --
same connectivity vectors as the frictionless version, MINUS k_nut, which
this sub-branch drops entirely -- see lugre_model.py/README.md) plus three
port-stiffness elements K_eq_port * outer(b_port, b_port), where b_nut is
IDENTICAL to the frictionless version's b_nut (same physical interface, same
sign convention), and b_sb/b_way are pure grounding vectors on theta_sb/x_n
(verified against build_linearized_matrices' own K assembly: its nut block
is exactly K_eq_nut*outer(b_nut,b_nut), the K_eq_sb term appears only on the
theta_sb diagonal, and K_eq_way only on the x_n diagonal).

MANDATORY CHECK: sum(K_e) must equal the assembled K from
build_linearized_matrices() to numerical tolerance -- same hard-failure
policy as the frictionless version.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from lugre_model import load_parameters
from run_local_linearization_bode import (
    STATE_LABELS,
    V_STAGE,
    build_linearized_matrices,
    equivalent_stiffness_damping,
)

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parent
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from plot_colors import ELEMENT_COLORS as _FRICTIONLESS_ELEMENT_COLORS  # noqa: E402

IDX = {lbl: i for i, lbl in enumerate(STATE_LABELS)}
N = len(STATE_LABELS)

_STRUCTURAL_ORDER = ["k_EM", "k_d", "k_c", "k_s1", "k_s2", "k_brg"]
STRUCTURAL_COLORS = dict(zip(_STRUCTURAL_ORDER, _FRICTIONLESS_ELEMENT_COLORS[:6]))
PORT_ELEMENT_NAMES = {"nut": "sigma0_nut", "sb": "sigma0_sb", "way": "sigma0_way"}
PORT_COLORS = {"nut": "#6b3fa0", "sb": "#2b6cb0", "way": "#c05621"}


def unit(label: str) -> np.ndarray:
    v = np.zeros(N)
    v[IDX[label]] = 1.0
    return v


def build_structural_elements(p: dict[str, float]) -> dict[str, np.ndarray]:
    k_EM = p["N_r"] * p["T_hold"]
    k_d = 4.0 * p["N_r"] * p["T_d"]

    b_EM = unit("theta_m")
    b_d = unit("theta_m")
    b_c = unit("theta_m") - unit("theta_c")
    b_s1 = unit("theta_c") - unit("theta_s")
    b_s2 = unit("theta_s") - unit("theta_sb")
    b_brg = unit("x_s")

    elements = {
        "k_EM": (k_EM, b_EM), "k_d": (k_d, b_d), "k_c": (p["k_c"], b_c),
        "k_s1": (p["k_s1"], b_s1), "k_s2": (p["k_s2"], b_s2), "k_brg": (p["k_brg"], b_brg),
    }
    return {name: k * np.outer(b, b) for name, (k, b) in elements.items()}


def build_port_elements(p: dict[str, float], K_eq: dict[str, float]) -> dict[str, np.ndarray]:
    lead_ratio = p["L"] / (2.0 * np.pi)
    b_nut = np.zeros(N)
    b_nut[IDX["theta_s"]] = -lead_ratio
    b_nut[IDX["x_s"]] = -1.0
    b_nut[IDX["x_n"]] = 1.0
    b_sb = unit("theta_sb")
    b_way = unit("x_n")
    b = {"nut": b_nut, "sb": b_sb, "way": b_way}
    return {PORT_ELEMENT_NAMES[port]: K_eq[port] * np.outer(b[port], b[port]) for port in ("nut", "sb", "way")}


def main() -> None:
    p = load_parameters()
    M, K, C, B_em = build_linearized_matrices(p)

    omega0 = V_STAGE * 2.0 * np.pi / p["L"]
    v0 = {"way": V_STAGE, "sb": omega0, "nut": 0.0}
    K_eq = {}
    K_eq["nut"], _, _, _ = equivalent_stiffness_damping(
        v0["nut"], p["sigma0_nut"], p["sigma1_nut"], p["sigma2_nut"], p["Fc_nut"], p["Fs_nut"], p["vs_nut"])
    K_eq["sb"], _, _, _ = equivalent_stiffness_damping(
        v0["sb"], p["sigma0_sb"], p["sigma1_sb"], p["sigma2_sb"], p["Tc_sb"], p["Ts_sb"], p["vs_sb"])
    K_eq["way"], _, _, _ = equivalent_stiffness_damping(
        v0["way"], p["sigma0_way"], p["sigma1_way"], p["sigma2_way"], p["Fc_way"], p["Fs_way"], p["vs_way"])

    K_elements = build_structural_elements(p)
    K_elements.update(build_port_elements(p, K_eq))

    K_sum = sum(K_elements.values())
    max_err = np.max(np.abs(K_sum - K))
    print(f"MANDATORY CHECK: max|sum(K_e) - K| = {max_err:.3e} (K scale ~{np.max(np.abs(K)):.3e})")
    assert np.allclose(K_sum, K, atol=1e-6, rtol=1e-9), (
        "Element stiffness decomposition (structural + 3 frozen friction ports) "
        "does not reconstruct K -- a term is missing or a sign is wrong; "
        "every fraction below would be invalid."
    )
    print("  -> PASS: element stiffnesses (structural + frozen ports) reconstruct K exactly.\n")

    lam, phi = eigh(K, M)
    omega = np.sqrt(lam)
    freq_hz = omega / (2.0 * np.pi)
    n_modes = len(lam)

    element_names = _STRUCTURAL_ORDER + [PORT_ELEMENT_NAMES[port] for port in ("nut", "sb", "way")]
    element_colors = dict(STRUCTURAL_COLORS)
    element_colors.update({PORT_ELEMENT_NAMES[port]: PORT_COLORS[port] for port in ("nut", "sb", "way")})

    SE = np.zeros((len(element_names), n_modes))
    for e_idx, name in enumerate(element_names):
        K_e = K_elements[name]
        for j in range(n_modes):
            SE[e_idx, j] = (phi[:, j] @ K_e @ phi[:, j]) / lam[j]

    col_sums = SE.sum(axis=0)
    print(f"{'mode':>4s} {'f (Hz)':>10s}  " + "  ".join(f"{nm:>10s}" for nm in element_names) + f"  {'sum':>8s}")
    for j in range(n_modes):
        vals = "  ".join(f"{SE[e_idx, j]:10.4f}" for e_idx in range(len(element_names)))
        print(f"{j+1:4d} {freq_hz[j]:10.2f}  {vals}  {col_sums[j]:8.4f}")
    assert np.allclose(col_sums, 1.0, atol=1e-6), "Strain-energy fractions do not sum to 1 per mode."

    fig, ax = plt.subplots(figsize=(11.0, 6.5))
    x_pos = np.arange(n_modes)
    bottom = np.zeros(n_modes)
    for e_idx, name in enumerate(element_names):
        hatch = "//" if name in PORT_ELEMENT_NAMES.values() else None
        ax.bar(x_pos, SE[e_idx], bottom=bottom, color=element_colors[name], label=name, width=0.6,
               hatch=hatch, edgecolor="#333333" if hatch else "none", linewidth=0.4 if hatch else 0)
        bottom += SE[e_idx]

    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"mode {j+1}\n{freq_hz[j]:.0f} Hz" for j in range(n_modes)], fontsize=8)
    ax.set_ylabel("Fraction of modal strain energy")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Modal strain energy decomposition, frozen-linearized LuGre system\n"
                 r"SE$[e,j] = \phi_j^T K_e \phi_j / \lambda_j$, stacked per mode (sums to 1); "
                 "hatched = LuGre port stiffness ($\\sigma_0$), replacing the frictionless baseline's k_nut",
                 fontsize=9.5)
    ax.grid(True, axis="y", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=8, ncol=len(element_names), loc="upper center", bbox_to_anchor=(0.5, -0.12))

    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)

    out_path = OUT_DIR / "modal_strain_energy_linearized.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "modal_strain_energy_linearized_data.npz",
             freq_hz=freq_hz, SE=SE, element_names=np.array(element_names))

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'modal_strain_energy_linearized_data.npz'}")


if __name__ == "__main__":
    main()

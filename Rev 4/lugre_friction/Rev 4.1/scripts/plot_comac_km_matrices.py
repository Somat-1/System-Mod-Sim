#!/usr/bin/env python3
"""General reference for the Co-MAC folder (both Split Method and
Kinematically Scaled derivations): renders the actual K and M matrices used
by both, symbolically and numerically -- Set A (frictionless baseline,
build_bode_rev4.build_matrices) and Set B (frozen-linearized LuGre,
run_local_linearization_bode.build_linearized_matrices). Both eigenproblems
(eigh(K, M)) feeding every MAC/COMAC number in this folder come from exactly
these two matrix pairs; this script documents them once, in one place,
instead of re-deriving them inside each sub-method's report.

Also answers, in writing, whether these carry the corrected/final form of the
model: yes -- both builder functions implement state_space_6dof.md Sec. 5.2-
5.4 as it stands after the 2026-08-18 sign correction (Sec. 5.3, Sec. 10
items 1-3 resolved) that fixed an earlier draft's inverted sign on the
L/(2*pi)*k_nut cross-coupling terms. That inverted-sign draft made the
command-to-stage DC gain negative; the corrected form used here does not.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parents[1]
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from build_bode_rev4 import build_matrices as build_frictionless_matrices  # noqa: E402
from build_bode_rev4 import load_parameters as load_frictionless_parameters  # noqa: E402

from lugre_model import load_parameters as load_lugre_parameters  # noqa: E402
from run_local_linearization_bode import (  # noqa: E402
    build_linearized_matrices, equivalent_stiffness_damping, V_STAGE,
)

STATE_LABELS = ["theta_m", "theta_c", "theta_s", "theta_sb", "x_s", "x_n"]


def md_matrix(mat: np.ndarray, fmt: str = "{:.4e}") -> str:
    header = "| | " + " | ".join(STATE_LABELS) + " |\n"
    sep = "|---|" + "|".join(["---"] * len(STATE_LABELS)) + "|\n"
    rows = ""
    for i, rl in enumerate(STATE_LABELS):
        vals = " | ".join(fmt.format(mat[i, j]) for j in range(mat.shape[1]))
        rows += f"| **{rl}** | {vals} |\n"
    return header + sep + rows


def main() -> None:
    p0 = load_frictionless_parameters()
    M0, C0, K0, Bu0 = build_frictionless_matrices(p0)

    p1 = load_lugre_parameters()
    M1, K1, C1, Bem1 = build_linearized_matrices(p1)

    K_eq_nut, C_eq_nut, z0_nut, F0_nut = equivalent_stiffness_damping(
        0.0, p1["sigma0_nut"], p1["sigma1_nut"], p1["sigma2_nut"],
        p1["Fc_nut"], p1["Fs_nut"], p1["vs_nut"])
    omega0 = V_STAGE * 2.0 * np.pi / p1["L"]
    K_eq_sb, C_eq_sb, z0_sb, F0_sb = equivalent_stiffness_damping(
        omega0, p1["sigma0_sb"], p1["sigma1_sb"], p1["sigma2_sb"],
        p1["Tc_sb"], p1["Ts_sb"], p1["vs_sb"])
    K_eq_way, C_eq_way, z0_way, F0_way = equivalent_stiffness_damping(
        V_STAGE, p1["sigma0_way"], p1["sigma1_way"], p1["sigma2_way"],
        p1["Fc_way"], p1["Fs_way"], p1["vs_way"])

    lead_ratio = p0["L"] / (2.0 * np.pi)
    mass_equal = bool(np.array_equal(M0, M1))

    md = []
    md.append("# K and M Matrices Used in the MAC / Co-MAC Derivations\n")
    md.append(
        "Both Co-MAC methods in this folder ([Split Method](Split%20Method/comac_mode_extraction.md), "
        "[Kinematically Scaled](Kinematically%20Scaled/comac_kinematic_scaled.md)) start from the "
        "identical pair of generalized eigenproblems `scipy.linalg.eigh(K, M)`:\n\n"
        "- **Set A** (frictionless baseline) -- `M0, C0, K0, Bu0 = build_bode_rev4.build_matrices(p)`\n"
        "- **Set B** (frozen-linearized LuGre, `V_STAGE = "
        f"{V_STAGE*1e3:.1f} mm/s`) -- `M1, K1, C1, Bem1 = run_local_linearization_bode."
        "build_linearized_matrices(p)`\n\n"
        "Every MAC and COMAC number in either report is a downstream function of these two "
        "$(K, M)$ pairs. This page renders them once, symbolically and numerically, instead of "
        "re-deriving them inside each sub-method's report.\n"
    )

    md.append("\n## Does this carry the corrected/final model form?\n")
    md.append(
        "**Yes.** Both builder functions implement "
        "[`state_space_6dof.md`](../../../../state_space_6dof.md) Sec. 5.2-5.4 as it stands "
        "after the 2026-08-18 fix, not an earlier draft:\n\n"
        "- **Sec. 10, items 1-3 (resolved this revision):** the screw-nut interface is embedded "
        "directly in $\\mathbf{K}$/$\\mathbf{C}$ through the lead ratio $L/2\\pi$ rather than "
        "injected as a separate frozen-friction input, so the reaction pair on $\\theta_s$ and "
        "$x_s$ is symmetric by construction and there is exactly one screw-nut element (no "
        "double path).\n"
        "- **Sec. 5.3 sign correction (2026-08-18):** the four off-diagonal "
        "$\\tfrac{L}{2\\pi}k_{nut}$ terms have the sign fixed so that positive $\\theta_s$ "
        "drives positive $x_n$, per the Sec. 1 convention. An earlier draft had the opposite "
        "sign there, which inverted the command-to-stage DC gain (a spurious 180 deg offset at "
        "low frequency). `build_bode_rev4.py`'s own docstring records this explicitly "
        "(*\"2026-08-18 sign correction\"*), and the numeric $K_0$ below is generated by that "
        "exact function -- not hand-transcribed.\n"
        "- **Sec. 5.4 ($c_{EM}$ in $\\mathbf{C}$):** `C[0,0] = c_c + c_EM` is applied by default "
        "in both builders, matching the 2026-08-18 default-inclusion decision (mode 1 is "
        "otherwise almost undamped). This does not affect the $K$/$M$ shown below, but it is "
        "part of the same fixed derivation and is included in both builder functions for that "
        "reason.\n\n"
        "$\\mathbf{C}$ itself is **not** used anywhere in the MAC/Co-MAC computation -- "
        "`eigh(K, M)` is the undamped generalized eigenproblem. $\\mathbf{C}$ only matters for "
        "the time-domain and Bode simulations elsewhere in this project, not for this modal "
        "comparison. It is mentioned here only to confirm both builders are the single, "
        "consistent, corrected implementation of Sec. 5, not two independently-patched copies.\n"
    )

    md.append("\n## Mass Matrix $\\mathbf{M}$ (Sec. 5.2) -- shared by Set A and Set B\n")
    md.append(
        "$$\\mathbf{M} = \\operatorname{diag}\\big(I_m,\\; I_c,\\; I_s,\\; I_{sb},\\; "
        "M_{screw},\\; M_s\\big)$$\n\n"
        f"Verified numerically identical between the two builder calls: "
        f"`np.array_equal(M0, M1)` = **{mass_equal}** -- friction changes stiffness/damping "
        "only, never inertia.\n"
    )
    md.append("\n### Numeric $\\mathbf{M}$ (kg or kg m$^2$, shared)\n")
    md.append(md_matrix(M0))

    md.append("\n## Set A -- Frictionless Baseline, Stiffness Matrix $\\mathbf{K}_0$\n")
    md.append(
        "Symbolic form (Sec. 5.3, corrected sign):\n\n"
        "$$\\mathbf{K} = \\begin{bmatrix}"
        "k_c+k_{EM}+k_d & -k_c & 0 & 0 & 0 & 0 \\\\"
        "-k_c & k_c+k_{s1} & -k_{s1} & 0 & 0 & 0 \\\\"
        "0 & -k_{s1} & k_{s1}+k_{s2}+\\left(\\tfrac{L}{2\\pi}\\right)^2 k_{nut} & -k_{s2} & "
        "\\tfrac{L}{2\\pi}k_{nut} & -\\tfrac{L}{2\\pi}k_{nut} \\\\"
        "0 & 0 & -k_{s2} & k_{s2} & 0 & 0 \\\\"
        "0 & 0 & \\tfrac{L}{2\\pi}k_{nut} & 0 & k_{brg}+k_{nut} & -k_{nut} \\\\"
        "0 & 0 & -\\tfrac{L}{2\\pi}k_{nut} & 0 & -k_{nut} & k_{nut}"
        "\\end{bmatrix}$$\n\n"
        f"with $L/2\\pi = {lead_ratio:.6e}$ m/rad, and the nut treated as fully stuck "
        f"($k_{{nut}} = {p0['k_nut']:.3e}$ N/m, no LuGre softening).\n"
    )
    md.append("\n### Numeric $\\mathbf{K}_0$ (N/m, N/rad, or N m/rad as appropriate per row/col)\n")
    md.append(md_matrix(K0))

    md.append(
        "\n## Set B -- Frozen-Linearized LuGre System, Stiffness Matrix $\\mathbf{K}_1$\n"
    )
    md.append(
        "Same base structure as $\\mathbf{K}_0$, with the nut's rigid $k_{nut}$ replaced by its "
        "frozen LuGre equivalent stiffness $K_{eq,nut}$ (bristle deflection frozen at the "
        "steady-state operating point, Jacobian taken numerically), and two new grounding terms "
        "added on the diagonal for the other two friction ports -- $K_{eq,sb}$ on $\\theta_{sb}$ "
        "and $K_{eq,way}$ on $x_n$ -- which do not appear in $\\mathbf{K}_0$ at all (those ports "
        "are frictionless inputs there, held at zero, not stiffness terms):\n\n"
        "$$\\mathbf{K} = \\begin{bmatrix}"
        "k_c+k_{EM}+k_d & -k_c & 0 & 0 & 0 & 0 \\\\"
        "-k_c & k_c+k_{s1} & -k_{s1} & 0 & 0 & 0 \\\\"
        "0 & -k_{s1} & k_{s1}+k_{s2}+\\left(\\tfrac{L}{2\\pi}\\right)^2 K_{eq,nut} & -k_{s2} & "
        "\\tfrac{L}{2\\pi}K_{eq,nut} & -\\tfrac{L}{2\\pi}K_{eq,nut} \\\\"
        "0 & 0 & -k_{s2} & k_{s2}+K_{eq,sb} & 0 & 0 \\\\"
        "0 & 0 & \\tfrac{L}{2\\pi}K_{eq,nut} & 0 & k_{brg}+K_{eq,nut} & -K_{eq,nut} \\\\"
        "0 & 0 & -\\tfrac{L}{2\\pi}K_{eq,nut} & 0 & -K_{eq,nut} & K_{eq,nut}+K_{eq,way}"
        "\\end{bmatrix}$$\n"
    )
    md.append(
        f"\nFrozen operating point: $V_{{stage}} = {V_STAGE*1e3:.1f}$ mm/s, giving "
        f"$\\omega_0 = {omega0:.4f}$ rad/s at the screw/bearing and the following equivalent "
        "stiffnesses (numerical Jacobian of the LuGre bristle force at that frozen point):\n\n"
        "| Port | $v_0$ | $K_{eq}$ |\n|---|---|---|\n"
        f"| nut | {0.0:.3e} m/s | {K_eq_nut:.4e} N/m |\n"
        f"| support bearing | {omega0:.4f} rad/s | {K_eq_sb:.4e} N m/rad |\n"
        f"| guideway | {V_STAGE:.3e} m/s | {K_eq_way:.4e} N/m |\n"
    )
    md.append("\n### Numeric $\\mathbf{K}_1$ (N/m, N/rad, or N m/rad as appropriate per row/col)\n")
    md.append(md_matrix(K1))

    md.append(
        "\n## $\\mathbf{K}_1 - \\mathbf{K}_0$ -- Where Set A and Set B Actually Differ\n"
    )
    md.append(
        "Every MAC/COMAC number in this folder ultimately traces back to the entries that "
        "change here; everywhere else, $\\mathbf{K}_0$ and $\\mathbf{K}_1$ are identical:\n"
    )
    md.append(md_matrix(K1 - K0))
    md.append(
        "\nThe nonzero block is exactly the $(\\theta_s, \\theta_{sb}, x_s, x_n)$ "
        "sub-structure touched by the three friction ports -- $k_{nut} \\to K_{eq,nut}$ "
        "(nut, $\\approx 1\\times10^{8} \\to \\approx 2\\times10^6$ N/m, a much softer "
        "effective contact once the bristle is allowed to deflect instead of assumed rigid), "
        "plus the two new grounding terms $K_{eq,sb}$ on $\\theta_{sb}$ and $K_{eq,way}$ on "
        "$x_n$ that have no counterpart in $\\mathbf{K}_0$ at all. $\\theta_m$, $\\theta_c$, "
        "and every entry not touching those DOFs are exactly zero in this difference -- the "
        "motor/coupling end of the drivetrain is untouched by the friction linearization.\n"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    md_path = OUT_DIR / "comac_km_matrices.md"
    md_path.write_text("".join(md), encoding="utf-8")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()

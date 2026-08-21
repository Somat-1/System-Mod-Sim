#!/usr/bin/env python3
"""Co-MAC (Kinematically Scaled) -- alternative to the Split Method
(../rendered_assets/temp/Co-MAC/Split Method/): instead of computing two
separate 4-DOF/2-DOF MAC matrices, rescale the 4 rotational rows of both raw
eigenvector matrices to equivalent linear displacement using the lead-screw
ratio R = L/(2*pi), then run a single 6-DOF Global MAC, mode-pairing, sign
alignment, and COMAC over all 6 coordinates together.

Same two eigensolves as the Split Method (scipy.linalg.eigh(K, M) for the
frictionless baseline and the frozen-linearized LuGre system) -- not repeated
in full here; see Split Method/comac_mode_extraction.md for the raw
eigenvalue/eigenvector tables. This report picks up from those raw phi0/phi1
and documents only the scaling -> MAC -> pairing -> sign-fix -> COMAC path.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from lugre_model import load_parameters as load_lugre_parameters
from run_local_linearization_bode import STATE_LABELS, build_linearized_matrices

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parent
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC" / "Kinematically Scaled"
NPZ_DIR = OUT_DIR / "npz"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from build_bode_rev4 import build_matrices as build_frictionless_matrices  # noqa: E402
from build_bode_rev4 import load_parameters as load_frictionless_parameters  # noqa: E402

N = len(STATE_LABELS)
ROT_LABELS = STATE_LABELS[:4]
TRANS_LABELS = STATE_LABELS[4:6]


def md_matrix(mat: np.ndarray, row_labels: list[str], col_labels: list[str], fmt: str = "{:.4e}") -> str:
    header = "| | " + " | ".join(col_labels) + " |\n"
    sep = "|---|" + "|".join(["---"] * len(col_labels)) + "|\n"
    rows = ""
    for i, rl in enumerate(row_labels):
        vals = " | ".join(fmt.format(mat[i, j]) for j in range(mat.shape[1]))
        rows += f"| **{rl}** | {vals} |\n"
    return header + sep + rows


def compute_mac(phi_a: np.ndarray, phi_b: np.ndarray) -> np.ndarray:
    numer = np.abs(phi_a.T @ phi_b) ** 2
    norm_a = np.sum(phi_a * phi_a, axis=0)
    norm_b = np.sum(phi_b * phi_b, axis=0)
    return numer / np.outer(norm_a, norm_b)


def main() -> None:
    p0 = load_frictionless_parameters()
    M0, C0, K0, B_u0 = build_frictionless_matrices(p0)
    lam0, phi0 = eigh(K0, M0)
    omega0 = np.sqrt(lam0)
    freq0_hz = omega0 / (2.0 * np.pi)

    p1 = load_lugre_parameters()
    M1, K1, C1, B_em1 = build_linearized_matrices(p1)
    lam1, phi1 = eigh(K1, M1)
    omega1 = np.sqrt(lam1)
    freq1_hz = omega1 / (2.0 * np.pi)

    assert p0["L"] == p1["L"], "Set A and Set B must share the same lead-screw L for R to be meaningful."
    L = p0["L"]
    R = L / (2.0 * np.pi)

    mode_cols = [f"mode{j+1}" for j in range(N)]
    a_mode_labels = [f"A-mode{i+1}" for i in range(N)]
    b_mode_labels = [f"B-mode{j+1}" for j in range(N)]

    md = []
    md.append("# Co-MAC (Kinematically Scaled)\n")
    md.append(
        "Alternative to the [Split Method](../Split%20Method/comac_mode_extraction.md), which "
        "computed two separate MAC matrices (4-DOF rotational, 2-DOF translational) to sidestep "
        "the unit mismatch between $\\theta$ (rad) and $x$ (m). Here the mismatch is removed at "
        "the source instead: the 4 rotational rows are rescaled to equivalent linear "
        "displacement using the lead-screw kinematic relationship, then all 6 DOFs are "
        "carried through MAC, mode-pairing, sign-alignment, and COMAC together, in one pass.\n\n"
        "Set A (frictionless baseline) and Set B (frozen-linearized LuGre, `V_STAGE = 5 mm/s`) "
        "are the identical two `scipy.linalg.eigh(K, M)` solves used in the Split Method -- not "
        "reproduced in full here (see that report for the raw $\\lambda$, $\\omega$, $f$, and "
        "unscaled $\\phi$ tables). This report starts from those same raw, unfixed-sign, "
        "mass-normalized eigenvectors $\\phi_0$ (Set A) and $\\phi_1$ (Set B), state order "
        "`theta_m, theta_c, theta_s, theta_sb, x_s, x_n`.\n"
    )

    md.append("\n## Step 1 -- Kinematic Scaling of the Rotational Rows\n")
    md.append(
        "The lead screw ideal kinematic relationship (ball-screw/nut, no slip) ties axial "
        "travel to screw rotation by $x = \\dfrac{L}{2\\pi}\\theta$, where $L$ is the screw "
        "lead. Both models share the same physical lead screw, so both use the same $L$:\n\n"
        f"$$L = {L:.3e}\\ \\text{{m}}, \\qquad R = \\frac{{L}}{{2\\pi}} = {R:.6e}\\ "
        "\\text{m/rad}$$\n\n"
        "Multiply rows 1-4 ($\\theta_m, \\theta_c, \\theta_s, \\theta_{sb}$) of both raw "
        "eigenvector matrices by $R$; rows 5-6 ($x_s, x_n$) are already in meters and are left "
        "untouched:\n\n"
        "$$\\vec{\\phi}_{A,\\text{scaled},l} = \\begin{bmatrix}R\\,\\theta_m \\\\ R\\,\\theta_c "
        "\\\\ R\\,\\theta_s \\\\ R\\,\\theta_{sb} \\\\ x_s \\\\ x_n\\end{bmatrix}_{A,l}, "
        "\\qquad \\vec{\\phi}_{B,\\text{scaled},k} = \\begin{bmatrix}R\\,\\theta_m \\\\ "
        "R\\,\\theta_c \\\\ R\\,\\theta_s \\\\ R\\,\\theta_{sb} \\\\ x_s \\\\ "
        "x_n\\end{bmatrix}_{B,k}$$\n\n"
        "Why this works: every one of the 6 rows now reports the same physical quantity "
        "(equivalent linear displacement, m), in the same numerical range, so a 6-element dot "
        "product weighs all 6 DOFs on their actual structural contribution instead of letting "
        "whichever DOF group has larger raw magnitude dominate.\n"
    )

    scale = np.array([R, R, R, R, 1.0, 1.0])
    phi0_s = phi0 * scale[:, None]
    phi1_s = phi1 * scale[:, None]

    md.append("\n### Scaled eigenvector matrix phi0_scaled -- Set A\n")
    md.append(md_matrix(phi0_s, STATE_LABELS, mode_cols))
    md.append("\n### Scaled eigenvector matrix phi1_scaled -- Set B\n")
    md.append(md_matrix(phi1_s, STATE_LABELS, mode_cols))

    mac_scaled = compute_mac(phi0_s, phi1_s)
    md.append("\n## Step 2 -- Single 6-DOF Global MAC\n")
    md.append(
        "With every row in the same units, one standard MAC matrix now covers all 6 DOFs "
        "(same formula as the Split Method, applied to the scaled vectors instead of a "
        "sub-vector):\n\n"
        "$$\\mathrm{MAC}_{\\text{scaled}}(l,k) = \\frac{\\left|\\vec{\\phi}_{A,\\text{scaled},l}^"
        "{\\mathsf{T}}\\vec{\\phi}_{B,\\text{scaled},k}\\right|^2}"
        "{\\left(\\vec{\\phi}_{A,\\text{scaled},l}^{\\mathsf{T}}\\vec{\\phi}_{A,\\text{scaled},l}"
        "\\right)\\left(\\vec{\\phi}_{B,\\text{scaled},k}^{\\mathsf{T}}"
        "\\vec{\\phi}_{B,\\text{scaled},k}\\right)}$$\n"
    )
    md.append("\n### MAC_scaled (rows: Set A, cols: Set B)\n")
    md.append(md_matrix(mac_scaled, a_mode_labels, b_mode_labels, fmt="{:.4f}"))
    md.append("\nSee `comac_kinematic_mac_matrix.png` for this table as a heatmap (green border = row-wise argmax, same convention as the Split Method figure).\n")

    master_pairs = []
    for i in range(N):
        j = int(np.argmax(mac_scaled[i, :]))
        master_pairs.append((i, j))

    b_used = [b for _, b in master_pairs]
    duplicated_b = sorted(set(b for b in b_used if b_used.count(b) > 1))

    md.append("\n## Step 3 -- Master Pairing List\n")
    md.append(
        "Row-wise argmax of MAC_scaled -- one Set-B mode picked per Set-A mode:\n"
    )
    rows = "| Set A mode | Set B mode | MAC |\n|---|---|---|\n"
    for i, j in master_pairs:
        rows += f"| A-mode{i+1} | B-mode{j+1} | {mac_scaled[i, j]:.4f} |\n"
    md.append(rows)
    if duplicated_b:
        dup_str = ", ".join(f"B-mode{b+1}" for b in duplicated_b)
        md.append(
            f"\n**Not a clean bijection**: {dup_str} is claimed as the best match by more than "
            "one Set-A mode above (A-mode5 and A-mode6 both peak at B-mode4 here). Preserving "
            "all 6 DOFs together removes the severe aliasing seen in the Split Method's "
            "translational block (that whole block was >= 0.93 with no separation at all), but "
            "it does not by itself guarantee a strict one-to-one assignment -- resolving that "
            "fully would need an assignment algorithm (e.g. Hungarian/linear-sum-assignment) on "
            "top of the row-wise argmax, which is not run here. Both rows are still carried "
            "through Steps 4-5 below as given, matching how the Split Method handled its own "
            "duplicate (B-mode6 claimed by two rows in the rotational sub-MAC).\n"
        )

    signs = []
    for i, j in master_pairs:
        s = np.sign(np.dot(phi0_s[:, i], phi1_s[:, j]))
        signs.append(1.0 if s == 0 else float(s))

    md.append("\n## Step 4 -- Sign Alignment\n")
    md.append(
        "For each pair (l, k) above, compute S = the dot product of the full 6-element scaled "
        "vectors phi_A_scaled_l and phi_B_scaled_k. If S is negative, eigh returned mode k of "
        "Set B with the opposite overall sign from mode l of Set A; multiply column k of "
        "phi1_scaled by -1 before using it in Step 5 (applied per pair, not per column, since "
        "B-mode4 above is reused by two different pairs and is not guaranteed to need the same "
        "sign both times):\n"
    )
    rows = "| Set A mode | Set B mode | S | sign applied |\n|---|---|---|---|\n"
    for (i, j), s in zip(master_pairs, signs):
        S = float(np.dot(phi0_s[:, i], phi1_s[:, j]))
        rows += f"| A-mode{i+1} | B-mode{j+1} | {S:.4e} | {'+' if s > 0 else '-'} |\n"
    md.append(rows)

    n_dof = N
    comac = np.zeros((n_dof, 4))
    for i in range(n_dof):
        num_sum = sum(phi0_s[i, a] * s * phi1_s[i, b] for s, (a, b) in zip(signs, master_pairs))
        den_a = sum(phi0_s[i, a] ** 2 for a, b in master_pairs)
        den_b = sum(phi1_s[i, b] ** 2 for a, b in master_pairs)
        comac[i] = [num_sum, den_a, den_b, num_sum ** 2 / (den_a * den_b)]

    md.append("\n## Step 5 -- Unified 6-DOF COMAC\n")
    md.append(
        "One loop over all 6 coordinates (no rotational/translational split), summing over "
        "all L=6 master-paired, sign-aligned modes at once:\n\n"
        "$$\\mathrm{COMAC}_i = \\frac{\\left(\\sum_{m=1}^{L}(\\phi_{A,\\text{scaled}})_{i,l_m}\\,"
        "(\\phi_{B,\\text{scaled}})_{i,k_m}\\right)^2}"
        "{\\sum_{m=1}^{L}(\\phi_{A,\\text{scaled}})_{i,l_m}^2 \\cdot "
        "\\sum_{m=1}^{L}(\\phi_{B,\\text{scaled}})_{i,k_m}^2}\\;\\in[0,1]$$\n"
    )
    md.append("\n### Unified COMAC\n")
    header = ("| DOF | sum(phiA*phiB) | sum(phiA^2) | sum(phiB^2) | "
              "COMAC |\n|---|---|---|---|---|\n")
    body = ""
    for lbl, (num_sum, den_a, den_b, val) in zip(STATE_LABELS, comac):
        body += f"| **{lbl}** | {num_sum:.4e} | {den_a:.4e} | {den_b:.4e} | {val:.4f} |\n"
    md.append(header + body)

    md.append(
        "\nReading these against the Split Method result: theta_sb (support bearing ring) is "
        f"clearly the lowest of all 6 DOFs at {comac[3, 3]:.3f} -- and it is the one DOF the "
        "LuGre bearing friction torque acts on directly, so this is the cleanest physical signal "
        f"either method has produced. theta_s is next-lowest ({comac[2, 3]:.3f}), with x_s close "
        f"behind ({comac[4, 3]:.3f}). Both translational DOFs move well away from the Split "
        f"Method's near-zero reading (x_s = {comac[4, 3]:.3f}, x_n = {comac[5, 3]:.3f}) -- "
        "consistent with that earlier near-zero result being an artifact of leaving the "
        "eigenvector sign unresolved before summing (Step 4 above is exactly the fix), not a "
        f"genuine finding about x_s/x_n. x_n and theta_m score highest ({comac[5, 3]:.3f}, "
        f"{comac[0, 3]:.3f}), i.e. least reshaped by the added friction; theta_c sits in between "
        f"({comac[1, 3]:.3f}). See `comac_kinematic_comac.png` for these six values plotted "
        "against the same COMAC >= 0.9 threshold used in the Split Method figure.\n"
    )

    md.append(
        "\n## Notes\n"
        "- R is a single scalar shared by both models (same physical lead screw, same L), "
        "applied identically to Set A and Set B -- it rescales units, it does not introduce any "
        "new information relating the two systems.\n"
        "- This MAC_scaled is not comparable cell-by-cell with the Split Method's "
        "MAC_rot/MAC_trans: those were computed on sub-vectors with theta still in raw radians "
        "(mass-normalized units, not physically scaled), so the resulting best-match pairs "
        "differ -- each method's pairing is only self-consistent within itself.\n"
        "- As in the Split Method, sign is otherwise arbitrary per raw mode (eigh guarantees "
        "phi.T @ M @ phi = I, not a consistent sign); Step 4 above resolves it only for the 6 "
        "pairs actually used in Step 5, not for the full raw phi0/phi1 matrices.\n"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)

    md_path = OUT_DIR / "comac_kinematic_scaled.md"
    md_path.write_text("".join(md), encoding="utf-8")

    np.savez(NPZ_DIR / "comac_kinematic_scaled_data.npz",
             state_labels=np.array(STATE_LABELS), R=R, L=L,
             phi0=phi0, phi1=phi1, phi0_scaled=phi0_s, phi1_scaled=phi1_s,
             lam0=lam0, freq0_hz=freq0_hz, lam1=lam1, freq1_hz=freq1_hz,
             mac_scaled=mac_scaled,
             master_pairs=np.array(master_pairs), signs=np.array(signs),
             comac=comac)

    print(f"Wrote {md_path}")
    print(f"Wrote {NPZ_DIR / 'comac_kinematic_scaled_data.npz'}")


if __name__ == "__main__":
    main()

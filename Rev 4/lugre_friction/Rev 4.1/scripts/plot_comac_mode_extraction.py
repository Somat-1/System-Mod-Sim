#!/usr/bin/env python3
"""Co-MAC Stage 1 -- raw normal-mode extraction, requested 2026-08-21;
revised 2026-08-21 to drop the figure and write a Markdown report instead
(the npz alone isn't readable without running code).

Extracts the two sets of normal modes (scipy.linalg.eigh(K, M), the
generalized SYMMETRIC eigenproblem) for:
  Set A: the frictionless baseline (../../scripts/build_bode_rev4.py)
  Set B: the frozen-linearized LuGre system (run_local_linearization_bode.py,
         frozen bristle equivalent stiffness/damping at V_STAGE=5 mm/s)

and presents both AS EXTRACTED -- no sign-fixing, no lead_ratio rescaling to
"equivalent axial displacement", no cross-matching between the two sets.
That post-processing is what plot_mode_shapes(_linearized).py already do
for readability, and mode correspondence (which baseline mode became which
LuGre mode) is a separate question for an actual MAC comparison -- this
script only extracts and verifies, hence "Co-MAC": the raw material a MAC
computation would consume, not the MAC itself.

Mass normalization: eigh(K, M) with M passed as the second argument returns
eigenvectors satisfying phi.T @ M @ phi = I by construction (LAPACK's
generalized symmetric-definite solver, driver dsygv under the hood). This is
verified here numerically for both sets (full phi.T @ M @ phi matrix
reported, not just asserted).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from lugre_model import load_parameters as load_lugre_parameters
from run_local_linearization_bode import STATE_LABELS, build_linearized_matrices

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parents[1]
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC" / "Split Method"
NPZ_DIR = OUT_DIR / "npz"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from build_bode_rev4 import build_matrices as build_frictionless_matrices  # noqa: E402
from build_bode_rev4 import load_parameters as load_frictionless_parameters  # noqa: E402

N = len(STATE_LABELS)


def md_matrix(mat: np.ndarray, row_labels: list[str], col_labels: list[str], fmt: str = "{:.4e}") -> str:
    header = "| | " + " | ".join(col_labels) + " |\n"
    sep = "|---|" + "|".join(["---"] * len(col_labels)) + "|\n"
    rows = ""
    for i, rl in enumerate(row_labels):
        vals = " | ".join(fmt.format(mat[i, j]) for j in range(mat.shape[1]))
        rows += f"| **{rl}** | {vals} |\n"
    return header + sep + rows


def md_eigen_table(lam: np.ndarray, omega: np.ndarray, freq_hz: np.ndarray) -> str:
    header = "| mode | lambda (rad^2/s^2) | omega (rad/s) | f (Hz) |\n|---|---|---|---|\n"
    rows = ""
    for j in range(len(lam)):
        rows += f"| {j+1} | {lam[j]:.6e} | {omega[j]:.4f} | {freq_hz[j]:.4f} |\n"
    return header + rows


def compute_mac(phi_a: np.ndarray, phi_b: np.ndarray) -> np.ndarray:
    """Standard MAC matrix between two sets of (possibly sub-partitioned) mode-shape
    column vectors. phi_a is (n_dof x p), phi_b is (n_dof x q); returns (p x q) with
    entries in [0, 1]. Sign- and scale-invariant (both drop out of the ratio), so it
    is safe to apply directly to the raw, unfixed-sign eigenvectors above.
    """
    numer = np.abs(phi_a.T @ phi_b) ** 2
    norm_a = np.sum(phi_a * phi_a, axis=0)
    norm_b = np.sum(phi_b * phi_b, axis=0)
    denom = np.outer(norm_a, norm_b)
    return numer / denom


def md_best_matches(mac: np.ndarray) -> str:
    rows = ""
    for i in range(mac.shape[0]):
        j = int(np.argmax(mac[i, :]))
        rows += f"| A-mode{i+1} | B-mode{j+1} | {mac[i, j]:.4f} |\n"
    return "| Set A mode | best-matching Set B mode | MAC |\n|---|---|---|\n" + rows


# Valid (unaliased) mode pairs for the Co-MAC step, 0-indexed (Set-A mode, Set-B mode).
# Excludes exactly the slots flagged in comac_mac_matrices.png:
#   - (A-mode2, B-mode1) in BOTH sub-MACs -- this is each row's naive argmax, but it is a
#     spurious low-frequency subspace-leakage artifact of splitting the rows (loses the
#     mass-matrix orthogonality weighting), not a genuine shape match.
#   - (A-mode3..6, B-mode5/B-mode6) in the TRANSLATIONAL sub-MAC only -- x_s/x_n carry too
#     little modal energy at these frequencies to discriminate between B-mode5/6, so every
#     value there is a near-tied >0.93 (aliased); A-mode4/5/6's only strong translational
#     candidates all fall inside this block, so those three rows have no valid translational
#     pair and are dropped rather than forced onto an unreliable match.
PAIRS_ROT = [(0, 0), (2, 2), (3, 4), (4, 5), (5, 5)]
PAIRS_TRANS = [(0, 0), (1, 1), (2, 2)]


def paired_signs(phi_a_full: np.ndarray, phi_b_full: np.ndarray,
                  pairs: list[tuple[int, int]]) -> list[float]:
    """Per-pair +/-1 sign so that phi_a_full[:,a] and sign*phi_b_full[:,b] point the
    same way (full 6-DOF eigenvector dot product, sign chosen once per pair -- not
    per Set-B column, since the same Set-B mode can be reused across two pairs, e.g.
    B-mode6 above, with no guarantee its required sign is the same both times).
    Required because eigh returns each mode with an arbitrary overall sign; COMAC's
    row-sum happens BEFORE squaring, so leaving that sign arbitrary would let terms
    cancel and understate a real correlation.
    """
    signs = []
    for a, b in pairs:
        s = np.sign(np.dot(phi_a_full[:, a], phi_b_full[:, b]))
        signs.append(1.0 if s == 0 else float(s))
    return signs


def compute_comac(phi_a_sub: np.ndarray, phi_b_sub: np.ndarray,
                   pairs: list[tuple[int, int]], signs: list[float]) -> np.ndarray:
    """COMAC per row (DOF) of phi_a_sub/phi_b_sub, summed over the given paired
    modes only. Returns an (n_dof, 4) array of [numerator_sum, denom_a, denom_b, comac]."""
    n_dof = phi_a_sub.shape[0]
    out = np.zeros((n_dof, 4))
    for i in range(n_dof):
        num_sum = sum(phi_a_sub[i, a] * s * phi_b_sub[i, b] for s, (a, b) in zip(signs, pairs))
        den_a = sum(phi_a_sub[i, a] ** 2 for a, b in pairs)
        den_b = sum(phi_b_sub[i, b] ** 2 for a, b in pairs)
        out[i] = [num_sum, den_a, den_b, num_sum ** 2 / (den_a * den_b)]
    return out


def md_comac_table(dof_labels: list[str], comac_data: np.ndarray) -> str:
    header = ("| DOF | sum($\\phi_A\\phi_B$) | sum($\\phi_A^2$) | sum($\\phi_B^2$) | "
              "COMAC |\n|---|---|---|---|---|\n")
    rows = ""
    for lbl, (num_sum, den_a, den_b, comac) in zip(dof_labels, comac_data):
        rows += f"| **{lbl}** | {num_sum:.4e} | {den_a:.4e} | {den_b:.4e} | {comac:.4f} |\n"
    return header + rows


def md_pairs_table(pairs: list[tuple[int, int]], signs: list[float], mac: np.ndarray) -> str:
    header = "| Set A mode | Set B mode | sign applied | MAC value |\n|---|---|---|---|\n"
    rows = ""
    for (a, b), s in zip(pairs, signs):
        rows += f"| A-mode{a+1} | B-mode{b+1} | {'+' if s > 0 else '-'} | {mac[a, b]:.4f} |\n"
    return header + rows


def main() -> None:
    # ---- Set A: frictionless baseline ----
    p0 = load_frictionless_parameters()
    M0, C0, K0, B_u0 = build_frictionless_matrices(p0)
    lam0, phi0 = eigh(K0, M0)
    omega0 = np.sqrt(lam0)
    freq0_hz = omega0 / (2.0 * np.pi)

    # ---- Set B: frozen-linearized LuGre system ----
    p1 = load_lugre_parameters()
    M1, K1, C1, B_em1 = build_linearized_matrices(p1)
    lam1, phi1 = eigh(K1, M1)
    omega1 = np.sqrt(lam1)
    freq1_hz = omega1 / (2.0 * np.pi)

    # ---- Mass-normalization check: phi.T @ M @ phi should be I ----
    norm_check0 = phi0.T @ M0 @ phi0
    norm_check1 = phi1.T @ M1 @ phi1
    max_off0 = np.max(np.abs(norm_check0 - np.eye(N)))
    max_off1 = np.max(np.abs(norm_check1 - np.eye(N)))
    pass0 = np.allclose(norm_check0, np.eye(N), atol=1e-8)
    pass1 = np.allclose(norm_check1, np.eye(N), atol=1e-8)

    mode_cols = [f"mode{j+1}" for j in range(N)]

    md = []
    md.append("# Co-MAC Stage 1 -- Raw Normal-Mode Extraction\n")
    md.append(
        "Two independent runs of `scipy.linalg.eigh(K, M)` (generalized symmetric "
        "eigenproblem), presented **as extracted** -- no sign-fixing, no lead_ratio "
        "rescaling to \"equivalent axial displacement\", no cross-matching between the "
        "two sets. That post-processing is what `plot_mode_shapes(_linearized).py` "
        "already do for readability; mode correspondence (which baseline mode became "
        "which LuGre mode) is a separate question for an actual MAC comparison -- this "
        "is the raw material a MAC computation would consume, not the MAC itself.\n"
    )
    md.append(
        "**State order** (rows below): `theta_m, theta_c, theta_s, theta_sb, x_s, x_n` "
        "-- rotational DOFs in rad, x_s/x_n in m.\n"
    )

    md.append("\n## Set A -- Frictionless Baseline\n")
    md.append(f"`{REV4_DIR.name}/scripts/build_bode_rev4.py` -- `M, C, K, B_u = build_matrices(p)`, "
              "`eigh(K, M)`.\n")
    md.append("\n### Eigenvalues\n")
    md.append(md_eigen_table(lam0, omega0, freq0_hz))
    md.append("\n### Raw eigenvector matrix phi (sign as returned by `eigh`, NOT fixed)\n")
    md.append(md_matrix(phi0, STATE_LABELS, mode_cols))

    md.append("\n## Set B -- Frozen-Linearized LuGre System\n")
    md.append("`run_local_linearization_bode.py` -- `M, K, C, B_em = build_linearized_matrices(p)` "
              "(frozen bristle equivalent stiffness/damping at V_STAGE=5 mm/s), `eigh(K, M)`.\n")
    md.append("\n### Eigenvalues\n")
    md.append(md_eigen_table(lam1, omega1, freq1_hz))
    md.append("\n### Raw eigenvector matrix phi (sign as returned by `eigh`, NOT fixed)\n")
    md.append(md_matrix(phi1, STATE_LABELS, mode_cols))

    # ---- Partition each mode shape into rotational (DOFs 1-4) / translational (DOFs 5-6) ----
    ROT_LABELS = STATE_LABELS[:4]     # theta_m, theta_c, theta_s, theta_sb
    TRANS_LABELS = STATE_LABELS[4:6]  # x_s, x_n
    phi0_rot, phi0_trans = phi0[:4, :], phi0[4:6, :]
    phi1_rot, phi1_trans = phi1[:4, :], phi1[4:6, :]

    md.append("\n## Partitioned Mode Shapes\n")
    md.append(
        "For every mode $l$, the 6-element raw eigenvector $\\vec{\\phi}_l$ (columns above) "
        "split into two sub-vectors:\n\n"
        "- Rotational: $\\vec{\\phi}_{\\text{rot},l} = "
        "\\begin{bmatrix}\\theta_m & \\theta_c & \\theta_s & \\theta_{sb}\\end{bmatrix}^{\\mathsf{T}}$ "
        "(DOFs 1-4, rad)\n"
        "- Translational: $\\vec{\\phi}_{\\text{trans},l} = \\begin{bmatrix}x_s & x_n\\end{bmatrix}^{\\mathsf{T}}$ "
        "(DOFs 5-6, m)\n\n"
        "Same raw (unsigned, unscaled) values as the full eigenvector tables above, just split out.\n"
    )

    md.append("\n### Set A -- Frictionless Baseline\n")
    md.append("\n**Rotational sub-vectors** $\\vec{\\phi}_{\\text{rot},l}$\n")
    md.append(md_matrix(phi0_rot, ROT_LABELS, mode_cols))
    md.append("\n**Translational sub-vectors** $\\vec{\\phi}_{\\text{trans},l}$\n")
    md.append(md_matrix(phi0_trans, TRANS_LABELS, mode_cols))

    md.append("\n### Set B -- Frozen-Linearized LuGre System\n")
    md.append("\n**Rotational sub-vectors** $\\vec{\\phi}_{\\text{rot},l}$\n")
    md.append(md_matrix(phi1_rot, ROT_LABELS, mode_cols))
    md.append("\n**Translational sub-vectors** $\\vec{\\phi}_{\\text{trans},l}$\n")
    md.append(md_matrix(phi1_trans, TRANS_LABELS, mode_cols))

    # ---- MAC between Set A and Set B, evaluated separately on each sub-vector ----
    mac_rot = compute_mac(phi0_rot, phi1_rot)
    mac_trans = compute_mac(phi0_trans, phi1_trans)

    a_mode_labels = [f"A-mode{i+1}" for i in range(N)]
    b_mode_labels = [f"B-mode{j+1}" for j in range(N)]

    md.append("\n## Modal Assurance Criterion (MAC) on Rotational and Translational Sub-Vectors\n")
    md.append(
        "The standard Modal Assurance Criterion between two real mode-shape vectors "
        "$\\vec{u}$, $\\vec{v}$ is\n\n"
        "$$\\mathrm{MAC}(\\vec{u}, \\vec{v}) = \\frac{\\left|\\vec{u}^{\\mathsf{T}}\\vec{v}\\right|^2}"
        "{\\left(\\vec{u}^{\\mathsf{T}}\\vec{u}\\right)\\left(\\vec{v}^{\\mathsf{T}}\\vec{v}\\right)}"
        "\\;\\in [0, 1],$$\n\n"
        "which is invariant to the sign and scale of either vector (both cancel between "
        "numerator and denominator), so it can be evaluated directly on the raw, "
        "unfixed-sign eigenvectors extracted above -- no sign-fixing or `lead_ratio` "
        "rescaling is needed first. A value of 1 means the two shapes are perfectly "
        "collinear; 0 means they are orthogonal (uncorrelated).\n\n"
        "Rather than evaluate this once on the full 6-element $\\vec{\\phi}_l$ (which mixes "
        "rotational DOFs, ~1e2-1e3 rad, with translational DOFs, ~1e-3-1e0 m, so the dot "
        "products would be dominated by whichever DOF group happens to have larger raw "
        "magnitude), the sub-vectors defined above are substituted for $\\vec{u}$, "
        "$\\vec{v}$ directly, giving two separate, dimensionally-consistent MAC matrices:\n\n"
        "- **Torsional global MAC** -- uses only the 4-element rotational sub-vectors:\n\n"
        "$$\\mathrm{MAC}_{\\text{rot}}(i,j) = \\frac{\\left|\\vec{\\phi}_{\\text{rot},A,i}^{\\mathsf{T}}"
        "\\vec{\\phi}_{\\text{rot},B,j}\\right|^2}"
        "{\\left(\\vec{\\phi}_{\\text{rot},A,i}^{\\mathsf{T}}\\vec{\\phi}_{\\text{rot},A,i}\\right)"
        "\\left(\\vec{\\phi}_{\\text{rot},B,j}^{\\mathsf{T}}\\vec{\\phi}_{\\text{rot},B,j}\\right)}$$\n\n"
        "with $\\vec{\\phi}_{\\text{rot}} = [\\theta_m\\ \\theta_c\\ \\theta_s\\ \\theta_{sb}]^{\\mathsf{T}}$ "
        "(4 elements).\n\n"
        "- **Linear/translational global MAC** -- uses only the 2-element translational "
        "sub-vectors:\n\n"
        "$$\\mathrm{MAC}_{\\text{trans}}(i,j) = \\frac{\\left|\\vec{\\phi}_{\\text{trans},A,i}^{\\mathsf{T}}"
        "\\vec{\\phi}_{\\text{trans},B,j}\\right|^2}"
        "{\\left(\\vec{\\phi}_{\\text{trans},A,i}^{\\mathsf{T}}\\vec{\\phi}_{\\text{trans},A,i}\\right)"
        "\\left(\\vec{\\phi}_{\\text{trans},B,j}^{\\mathsf{T}}\\vec{\\phi}_{\\text{trans},B,j}\\right)}$$\n\n"
        "with $\\vec{\\phi}_{\\text{trans}} = [x_s\\ x_n]^{\\mathsf{T}}$ (2 elements).\n\n"
        "In both cases $i$ indexes Set A (frictionless baseline) modes and $j$ indexes Set B "
        "(frozen-linearized LuGre) modes, so each matrix is a **cross**-MAC between the two "
        "systems, not a within-set self-MAC (which would just be the identity given "
        "mass-orthogonality). Per the Notes above, Set A mode $i$ and Set B mode $j$ are not "
        "presumed to correspond by index -- that correspondence is exactly what the MAC value "
        "at $(i,j)$ is testing for.\n"
    )

    md.append("\n### Torsional Global MAC -- MAC$_{\\text{rot}}$ (rows: Set A, cols: Set B)\n")
    md.append(md_matrix(mac_rot, a_mode_labels, b_mode_labels, fmt="{:.4f}"))
    md.append("\n**Best-matching Set B mode per Set A mode (by MAC$_{\\text{rot}}$)**\n")
    md.append(md_best_matches(mac_rot))

    md.append("\n### Linear/Translational Global MAC -- MAC$_{\\text{trans}}$ (rows: Set A, cols: Set B)\n")
    md.append(md_matrix(mac_trans, a_mode_labels, b_mode_labels, fmt="{:.4f}"))
    md.append("\n**Best-matching Set B mode per Set A mode (by MAC$_{\\text{trans}}$)**\n")
    md.append(md_best_matches(mac_trans))

    # ---- Co-MAC (Coordinate MAC): per-DOF row correlation over a filtered, paired,
    # sign-aligned subset of modes -- separately for the rotational and translational
    # sub-systems, per PAIRS_ROT / PAIRS_TRANS above ----
    signs_rot = paired_signs(phi0, phi1, PAIRS_ROT)
    signs_trans = paired_signs(phi0, phi1, PAIRS_TRANS)
    comac_rot = compute_comac(phi0_rot, phi1_rot, PAIRS_ROT, signs_rot)
    comac_trans = compute_comac(phi0_trans, phi1_trans, PAIRS_TRANS, signs_trans)

    md.append("\n### Reading the Global MAC Matrices -- Flagged Slots\n")
    md.append(
        "`comac_mac_matrices.png` plots both tables above as heatmaps (same palette used "
        "throughout this report). Two markings on it are not explained on the figure itself "
        "(no legend was added, to keep it uncluttered) and are the direct input to Stage 3's "
        "Step 2, so they are recorded here instead:\n\n"
        "- **Green border** -- the row-wise argmax (naive best match per Set-A mode), i.e. "
        "the *Best-matching Set B mode* rows tabulated above -- not a one-to-one assignment, "
        "just the single highest value in that row.\n"
        "- **Black cell, boxed value flanked by X** -- a slot flagged as unreliable and "
        "excluded from the Stage 3 pairing (below), even where it happens to also be the "
        "row's argmax. Two distinct reasons:\n\n"
        "  1. *(A-mode2, B-mode1), both sub-MACs.* Rotationally this is A-mode2's argmax "
        "(0.9127), but the runner-up in the same row, B-mode2 (0.8840), is barely behind it -- "
        "and B-mode2 is also A-mode2's actual best match translationally (0.9571 vs. only "
        "0.6931 for B-mode1, and B-mode2's 909 Hz sits far closer to A-mode2's 746 Hz than "
        "B-mode1's 391 Hz does). Splitting the mass-orthonormal 6-DOF eigenvector into rotational "
        "and translational halves discards the $M$-weighted orthogonality that keeps unrelated "
        "modes from cross-correlating in the full vector, so the rotational sub-MAC's slight "
        "edge for B-mode1 reads as leakage from that lost weighting, not a genuine shape match. "
        "B-mode2 is the physically consistent pairing; B-mode1 is not used for A-mode2 in "
        "either sub-system.\n"
        "  2. *(A-mode3..6, B-mode5/B-mode6), translational sub-MAC only.* Every entry in that "
        "2x2-plus block is $\\geq 0.93$ (see the Linear/Translational Global MAC table above) -- "
        "with only 2 elements ($x_s$, $x_n$) to distinguish 4 higher-frequency modes, the "
        "translational sub-vector cannot separate B-mode5 from B-mode6, or tell A-mode4/5/6 "
        "apart from each other against either. That is spatial aliasing (too few coordinates "
        "for the number of modes being discriminated), not a resolvable correspondence, so "
        "A-mode4/5/6 are dropped from the translational pairing entirely rather than assigned "
        "an arbitrary winner.\n"
    )

    md.append("\n## Co-MAC Stage 3 -- Coordinate Modal Assurance Criterion (COMAC)\n")
    md.append(
        "Global MAC (previous section) pairs whole *modes* across the two systems; COMAC "
        "goes the other way -- for a fixed set of already-paired modes, it correlates one "
        "*coordinate* (DOF) at a time, across all those modes, to show which specific "
        "rows of the drivetrain are most disrupted by adding friction. The standard "
        "definition, for DOF $i$ over $L$ paired modes $(l \\to k_l)$:\n\n"
        "$$\\mathrm{COMAC}_i = \\frac{\\left|\\sum_{l=1}^{L}(\\phi_A)_{i,l}\\,(\\phi_B)_{i,k_l}"
        "\\right|^2}{\\sum_{l=1}^{L}(\\phi_A)_{i,l}^2\\;\\cdot\\;\\sum_{l=1}^{L}(\\phi_B)_{i,k_l}^2}"
        "\\;\\in[0,1]$$\n\n"
        "**Step 1 -- Partition rows**: the rotational (4-row) and translational (2-row) "
        "blocks already split out above.\n\n"
        "**Step 2 -- Filter and pair columns**: take each sub-MAC's row-wise argmax "
        "(Torsional/Linear-Translational Global MAC tables above) and drop any pair that "
        "lands in a flagged/aliased slot (see `comac_mac_matrices.png`) rather than force a "
        "match. This also requires fixing the relative sign of each pair before summing -- "
        "`eigh` returns each mode with an arbitrary overall sign, and COMAC sums "
        "$(\\phi_A)_{i,l}(\\phi_B)_{i,k_l}$ **before** squaring, so an unresolved sign flip "
        "would let a genuinely good match cancel itself out. The sign used per pair (from "
        "the full 6-DOF dot product, applied consistently to both sub-systems) is listed "
        "alongside the pairs below.\n\n"
        "*Rotational pairs* (5 of 6 Set-A modes; A-mode2 dropped -- its only candidate, "
        "B-mode1, is the flagged subspace-leakage slot):\n"
    )
    md.append(md_pairs_table(PAIRS_ROT, signs_rot, mac_rot))
    md.append(
        "\n*Translational pairs* (3 of 6 Set-A modes; A-mode4/5/6 dropped -- their only "
        "candidates, B-mode5/B-mode6, fall inside the flagged aliased block):\n"
    )
    md.append(md_pairs_table(PAIRS_TRANS, signs_trans, mac_trans))
    md.append(
        "\n**Steps 3-5 -- Numerator, denominator, and normalize**: for each DOF row, sum the "
        "sign-corrected cross-products over the paired modes and square it (numerator), "
        "separately sum the squared elements of Set A and of Set B over those same modes and "
        "multiply them (denominator), then divide:\n"
    )
    md.append("\n### Rotational COMAC\n")
    md.append(md_comac_table(ROT_LABELS, comac_rot))
    md.append("\n### Translational COMAC\n")
    md.append(md_comac_table(TRANS_LABELS, comac_trans))
    md.append(
        "\nReading these: $\\theta_{sb}$ (support bearing ring, where the LuGre bearing "
        "friction torque acts) and $\\theta_c$ come out lowest among the rotational DOFs "
        "(0.63 and 0.44), meaning their relative modal contribution is the most reshaped by "
        "the added friction. Both translational DOFs collapse to near zero (<0.01) -- across "
        "the only 3 modes that could be validly paired, $x_s$ and $x_n$ barely correlate "
        "row-wise at all between the two systems, even though the *combined* 2-vector "
        "direction was highly correlated in the Global MAC table above; that combined score "
        "is apparently carried by whichever of $x_s$/$x_n$ has the larger raw magnitude in "
        "each mode, not by both coordinates moving consistently together across modes. With "
        "only $L=3$ pairs this translational result rests on a small sample and should be "
        "read as indicative, not conclusive. See `comac_stage3.png` for these six values "
        "plotted against the conventional COMAC >= 0.9 \"unaffected\" threshold.\n"
    )

    md.append("\n## Mass-Normalization Check\n")
    md.append(
        "`eigh(K, M)` with `M` passed as the second (generalized) argument returns "
        "eigenvectors satisfying `phi.T @ M @ phi = I` by construction (LAPACK's "
        "generalized symmetric-definite driver). Verified numerically below, not just "
        "cited.\n"
    )
    md.append("\n### Set A: `phi0.T @ M0 @ phi0`\n")
    md.append(md_matrix(norm_check0, mode_cols, mode_cols, fmt="{:.2e}"))
    md.append(f"\nmax \\|off-identity\\| = `{max_off0:.3e}`, diagonal range "
              f"`[{np.min(np.diag(norm_check0)):.9f}, {np.max(np.diag(norm_check0)):.9f}]` "
              f"-> **{'PASS' if pass0 else 'FAIL'}** (tolerance 1e-8)\n")

    md.append("\n### Set B: `phi1.T @ M1 @ phi1`\n")
    md.append(md_matrix(norm_check1, mode_cols, mode_cols, fmt="{:.2e}"))
    md.append(f"\nmax \\|off-identity\\| = `{max_off1:.3e}`, diagonal range "
              f"`[{np.min(np.diag(norm_check1)):.9f}, {np.max(np.diag(norm_check1)):.9f}]` "
              f"-> **{'PASS' if pass1 else 'FAIL'}** (tolerance 1e-8)\n")

    md.append(
        "\n## Notes\n"
        "- `M` differs in scale per DOF between rotational (`I_ii`, kg·m², ~1e-7-1e-6) "
        "and translational (`M_ii`, kg, ~0.1-0.4) rows -- mass-normalization is w.r.t. "
        "**this** `M`, not the identity, so raw `|phi|` entries are not directly "
        "comparable across DOF types without the `lead_ratio` rescaling "
        "`plot_mode_shapes(_linearized).py` apply for readability. Not applied here "
        "deliberately -- this report shows the extraction as `scipy` actually returns it.\n"
        "- Sign is arbitrary (`phi_j` and `-phi_j` are equally valid solutions); not "
        "fixed here, unlike the mode-shapes plots.\n"
        "- Set A and Set B mode indices are **not** claimed to correspond to each other "
        "-- e.g. Set A mode 1 (176.7 Hz) is not asserted to be \"the same mode\" as Set B "
        "mode 1 (391.5 Hz). Establishing that correspondence is a MAC computation, not "
        "done here.\n"
        "- The two MAC matrices above are independent of each other by construction -- "
        "a mode can correlate strongly in $\\mathrm{MAC}_{\\text{rot}}$ and weakly in "
        "$\\mathrm{MAC}_{\\text{trans}}$ (or vice versa) since they draw on disjoint DOF "
        "subsets. Neither is a substitute for a full-vector MAC; combining them (e.g. a "
        "weighted product) would be a further step beyond this Stage 1/2 extraction, not "
        "attempted here.\n"
    )

    assert pass0, "Set A eigenvectors are not mass-normalized."
    assert pass1, "Set B eigenvectors are not mass-normalized."

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)

    md_path = OUT_DIR / "comac_mode_extraction.md"
    md_path.write_text("".join(md), encoding="utf-8")

    np.savez(NPZ_DIR / "comac_mode_extraction_data.npz",
             state_labels=np.array(STATE_LABELS),
             rot_labels=np.array(ROT_LABELS), trans_labels=np.array(TRANS_LABELS),
             lam0=lam0, omega0=omega0, freq0_hz=freq0_hz, phi0=phi0, M0=M0, K0=K0,
             phi0_rot=phi0_rot, phi0_trans=phi0_trans, norm_check0=norm_check0,
             lam1=lam1, omega1=omega1, freq1_hz=freq1_hz, phi1=phi1, M1=M1, K1=K1,
             phi1_rot=phi1_rot, phi1_trans=phi1_trans, norm_check1=norm_check1,
             mac_rot=mac_rot, mac_trans=mac_trans,
             comac_rot=comac_rot, comac_trans=comac_trans,
             pairs_rot=np.array(PAIRS_ROT), pairs_trans=np.array(PAIRS_TRANS),
             signs_rot=np.array(signs_rot), signs_trans=np.array(signs_trans))

    print(f"Set A: {'PASS' if pass0 else 'FAIL'} (max off-identity {max_off0:.3e})")
    print(f"Set B: {'PASS' if pass1 else 'FAIL'} (max off-identity {max_off1:.3e})")
    print(f"\nWrote {md_path}")
    print(f"Wrote {NPZ_DIR / 'comac_mode_extraction_data.npz'}")


if __name__ == "__main__":
    main()

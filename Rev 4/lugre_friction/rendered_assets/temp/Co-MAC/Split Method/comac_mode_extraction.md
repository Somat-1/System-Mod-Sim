# Co-MAC Stage 1 -- Raw Normal-Mode Extraction
Two independent runs of `scipy.linalg.eigh(K, M)` (generalized symmetric eigenproblem), presented **as extracted** -- no sign-fixing, no lead_ratio rescaling to "equivalent axial displacement", no cross-matching between the two sets. That post-processing is what `plot_mode_shapes(_linearized).py` already do for readability; mode correspondence (which baseline mode became which LuGre mode) is a separate question for an actual MAC comparison -- this is the raw material a MAC computation would consume, not the MAC itself.
**State order** (rows below): `theta_m, theta_c, theta_s, theta_sb, x_s, x_n` -- rotational DOFs in rad, x_s/x_n in m.

## Set A -- Frictionless Baseline
`Rev 4/scripts/build_bode_rev4.py` -- `M, C, K, B_u = build_matrices(p)`, `eigh(K, M)`.

### Eigenvalues
| mode | lambda (rad^2/s^2) | omega (rad/s) | f (Hz) |
|---|---|---|---|
| 1 | 1.232441e+06 | 1110.1535 | 176.6864 |
| 2 | 2.196815e+07 | 4687.0194 | 745.9623 |
| 3 | 1.075091e+08 | 10368.6612 | 1650.2237 |
| 4 | 4.642834e+08 | 21547.2369 | 3429.3493 |
| 5 | 1.686554e+09 | 41067.6744 | 6536.1234 |
| 6 | 1.859705e+09 | 43124.2962 | 6863.4449 |

### Raw eigenvector matrix phi (sign as returned by `eigh`, NOT fixed)
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **theta_m** | 5.7747e+02 | -4.5000e+01 | 8.7195e+02 | 1.2373e+02 | -7.3050e-01 | -3.1859e+00 |
| **theta_c** | 5.9839e+02 | -3.4415e+01 | -3.0957e+02 | -6.2179e+02 | 1.5360e+01 | 7.4209e+01 |
| **theta_s** | 6.0109e+02 | -2.6738e+01 | -5.0842e+02 | 7.4976e+02 | -1.2427e+02 | -6.7237e+02 |
| **theta_sb** | 6.0161e+02 | -2.7162e+01 | -5.5049e+02 | 1.1191e+03 | 6.2455e+02 | 2.0877e+03 |
| **x_s** | 4.4415e-03 | 1.3330e+00 | 8.4478e-02 | -1.3429e-01 | 3.2621e+00 | -8.6556e-01 |
| **x_n** | 1.0061e-01 | 1.4585e+00 | 6.3072e-03 | 1.6992e-02 | -5.5609e-01 | 1.4890e-01 |

## Set B -- Frozen-Linearized LuGre System
`run_local_linearization_bode.py` -- `M, K, C, B_em = build_linearized_matrices(p)` (frozen bristle equivalent stiffness/damping at V_STAGE=5 mm/s), `eigh(K, M)`.

### Eigenvalues
| mode | lambda (rad^2/s^2) | omega (rad/s) | f (Hz) |
|---|---|---|---|
| 1 | 6.050317e+06 | 2459.7393 | 391.4797 |
| 2 | 3.264042e+07 | 5713.1792 | 909.2807 |
| 3 | 1.542755e+08 | 12420.7700 | 1976.8269 |
| 4 | 1.775571e+08 | 13325.0564 | 2120.7486 |
| 5 | 7.089187e+08 | 26625.5273 | 4237.5843 |
| 6 | 4.858072e+09 | 69699.8737 | 11093.0794 |

### Raw eigenvector matrix phi (sign as returned by `eigh`, NOT fixed)
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **theta_m** | 3.3494e+00 | 8.2910e+02 | 6.4865e+02 | 1.4552e+01 | -5.2250e+01 | -1.3224e-01 |
| **theta_c** | 3.2595e+00 | 5.1825e+02 | -6.2740e+02 | -1.8510e+01 | 4.2991e+02 | 8.2709e+00 |
| **theta_s** | 3.1199e+00 | 3.2236e+02 | -5.0187e+02 | -1.0903e+01 | -1.1174e+03 | -2.1370e+02 |
| **theta_sb** | 9.2705e-01 | 9.6329e+01 | -1.5395e+02 | -3.3614e+00 | -3.8992e+02 | 2.5459e+03 |
| **x_s** | 2.4222e-01 | -1.0809e-02 | 9.3125e-02 | -3.6229e+00 | -8.8830e-03 | -1.9189e-04 |
| **x_n** | 1.5678e+00 | -7.7436e-03 | -4.4371e-04 | 1.0483e-01 | 1.3133e-03 | 3.4816e-05 |

## Partitioned Mode Shapes
For every mode $l$, the 6-element raw eigenvector $\vec{\phi}_l$ (columns above) split into two sub-vectors:

- Rotational: $\vec{\phi}_{\text{rot},l} = \begin{bmatrix}\theta_m & \theta_c & \theta_s & \theta_{sb}\end{bmatrix}^{\mathsf{T}}$ (DOFs 1-4, rad)
- Translational: $\vec{\phi}_{\text{trans},l} = \begin{bmatrix}x_s & x_n\end{bmatrix}^{\mathsf{T}}$ (DOFs 5-6, m)

Same raw (unsigned, unscaled) values as the full eigenvector tables above, just split out.

### Set A -- Frictionless Baseline

**Rotational sub-vectors** $\vec{\phi}_{\text{rot},l}$
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **theta_m** | 5.7747e+02 | -4.5000e+01 | 8.7195e+02 | 1.2373e+02 | -7.3050e-01 | -3.1859e+00 |
| **theta_c** | 5.9839e+02 | -3.4415e+01 | -3.0957e+02 | -6.2179e+02 | 1.5360e+01 | 7.4209e+01 |
| **theta_s** | 6.0109e+02 | -2.6738e+01 | -5.0842e+02 | 7.4976e+02 | -1.2427e+02 | -6.7237e+02 |
| **theta_sb** | 6.0161e+02 | -2.7162e+01 | -5.5049e+02 | 1.1191e+03 | 6.2455e+02 | 2.0877e+03 |

**Translational sub-vectors** $\vec{\phi}_{\text{trans},l}$
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **x_s** | 4.4415e-03 | 1.3330e+00 | 8.4478e-02 | -1.3429e-01 | 3.2621e+00 | -8.6556e-01 |
| **x_n** | 1.0061e-01 | 1.4585e+00 | 6.3072e-03 | 1.6992e-02 | -5.5609e-01 | 1.4890e-01 |

### Set B -- Frozen-Linearized LuGre System

**Rotational sub-vectors** $\vec{\phi}_{\text{rot},l}$
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **theta_m** | 3.3494e+00 | 8.2910e+02 | 6.4865e+02 | 1.4552e+01 | -5.2250e+01 | -1.3224e-01 |
| **theta_c** | 3.2595e+00 | 5.1825e+02 | -6.2740e+02 | -1.8510e+01 | 4.2991e+02 | 8.2709e+00 |
| **theta_s** | 3.1199e+00 | 3.2236e+02 | -5.0187e+02 | -1.0903e+01 | -1.1174e+03 | -2.1370e+02 |
| **theta_sb** | 9.2705e-01 | 9.6329e+01 | -1.5395e+02 | -3.3614e+00 | -3.8992e+02 | 2.5459e+03 |

**Translational sub-vectors** $\vec{\phi}_{\text{trans},l}$
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **x_s** | 2.4222e-01 | -1.0809e-02 | 9.3125e-02 | -3.6229e+00 | -8.8830e-03 | -1.9189e-04 |
| **x_n** | 1.5678e+00 | -7.7436e-03 | -4.4371e-04 | 1.0483e-01 | 1.3133e-03 | 3.4816e-05 |

## Modal Assurance Criterion (MAC) on Rotational and Translational Sub-Vectors
The standard Modal Assurance Criterion between two real mode-shape vectors $\vec{u}$, $\vec{v}$ is

$$\mathrm{MAC}(\vec{u}, \vec{v}) = \frac{\left|\vec{u}^{\mathsf{T}}\vec{v}\right|^2}{\left(\vec{u}^{\mathsf{T}}\vec{u}\right)\left(\vec{v}^{\mathsf{T}}\vec{v}\right)}\;\in [0, 1],$$

which is invariant to the sign and scale of either vector (both cancel between numerator and denominator), so it can be evaluated directly on the raw, unfixed-sign eigenvectors extracted above -- no sign-fixing or `lead_ratio` rescaling is needed first. A value of 1 means the two shapes are perfectly collinear; 0 means they are orthogonal (uncorrelated).

Rather than evaluate this once on the full 6-element $\vec{\phi}_l$ (which mixes rotational DOFs, ~1e2-1e3 rad, with translational DOFs, ~1e-3-1e0 m, so the dot products would be dominated by whichever DOF group happens to have larger raw magnitude), the sub-vectors defined above are substituted for $\vec{u}$, $\vec{v}$ directly, giving two separate, dimensionally-consistent MAC matrices:

- **Torsional global MAC** -- uses only the 4-element rotational sub-vectors:

$$\mathrm{MAC}_{\text{rot}}(i,j) = \frac{\left|\vec{\phi}_{\text{rot},A,i}^{\mathsf{T}}\vec{\phi}_{\text{rot},B,j}\right|^2}{\left(\vec{\phi}_{\text{rot},A,i}^{\mathsf{T}}\vec{\phi}_{\text{rot},A,i}\right)\left(\vec{\phi}_{\text{rot},B,j}^{\mathsf{T}}\vec{\phi}_{\text{rot},B,j}\right)}$$

with $\vec{\phi}_{\text{rot}} = [\theta_m\ \theta_c\ \theta_s\ \theta_{sb}]^{\mathsf{T}}$ (4 elements).

- **Linear/translational global MAC** -- uses only the 2-element translational sub-vectors:

$$\mathrm{MAC}_{\text{trans}}(i,j) = \frac{\left|\vec{\phi}_{\text{trans},A,i}^{\mathsf{T}}\vec{\phi}_{\text{trans},B,j}\right|^2}{\left(\vec{\phi}_{\text{trans},A,i}^{\mathsf{T}}\vec{\phi}_{\text{trans},A,i}\right)\left(\vec{\phi}_{\text{trans},B,j}^{\mathsf{T}}\vec{\phi}_{\text{trans},B,j}\right)}$$

with $\vec{\phi}_{\text{trans}} = [x_s\ x_n]^{\mathsf{T}}$ (2 elements).

In both cases $i$ indexes Set A (frictionless baseline) modes and $j$ indexes Set B (frozen-linearized LuGre) modes, so each matrix is a **cross**-MAC between the two systems, not a within-set self-MAC (which would just be the identity given mass-orthogonality). Per the Notes above, Set A mode $i$ and Set B mode $j$ are not presumed to correspond by index -- that correspondence is exactly what the MAC value at $(i,j)$ is testing for.

### Torsional Global MAC -- MAC$_{\text{rot}}$ (rows: Set A, cols: Set B)
| | B-mode1 | B-mode2 | B-mode3 | B-mode4 | B-mode5 | B-mode6 |
|---|---|---|---|---|---|---|
| **A-mode1** | 0.8698 | 0.7159 | 0.1013 | 0.1307 | 0.2053 | 0.2147 |
| **A-mode2** | 0.9127 | 0.8840 | 0.0197 | 0.0417 | 0.1061 | 0.1334 |
| **A-mode3** | 0.0007 | 0.0788 | 0.7827 | 0.6865 | 0.1621 | 0.1814 |
| **A-mode4** | 0.0433 | 0.0071 | 0.0025 | 0.0012 | 0.6807 | 0.4979 |
| **A-mode5** | 0.0043 | 0.0017 | 0.0044 | 0.0039 | 0.0149 | 0.9869 |
| **A-mode6** | 0.0000 | 0.0001 | 0.0002 | 0.0004 | 0.0001 | 0.9481 |

**Best-matching Set B mode per Set A mode (by MAC$_{\text{rot}}$)**
| Set A mode | best-matching Set B mode | MAC |
|---|---|---|
| A-mode1 | B-mode1 | 0.8698 |
| A-mode2 | B-mode1 | 0.9127 |
| A-mode3 | B-mode3 | 0.7827 |
| A-mode4 | B-mode5 | 0.6807 |
| A-mode5 | B-mode6 | 0.9869 |
| A-mode6 | B-mode6 | 0.9481 |

### Linear/Translational Global MAC -- MAC$_{\text{trans}}$ (rows: Set A, cols: Set B)
| | B-mode1 | B-mode2 | B-mode3 | B-mode4 | B-mode5 | B-mode6 |
|---|---|---|---|---|---|---|
| **A-mode1** | 0.9881 | 0.3815 | 0.0015 | 0.0002 | 0.0105 | 0.0182 |
| **A-mode2** | 0.6931 | 0.9571 | 0.4504 | 0.4264 | 0.3129 | 0.2830 |
| **A-mode3** | 0.0510 | 0.7294 | 0.9937 | 0.9893 | 0.9518 | 0.9369 |
| **A-mode4** | 0.0008 | 0.5379 | 0.9854 | 0.9906 | 0.9996 | 0.9971 |
| **A-mode5** | 0.0002 | 0.4949 | 0.9733 | 0.9805 | 0.9995 | 0.9999 |
| **A-mode6** | 0.0003 | 0.4934 | 0.9728 | 0.9801 | 0.9994 | 0.9999 |

**Best-matching Set B mode per Set A mode (by MAC$_{\text{trans}}$)**
| Set A mode | best-matching Set B mode | MAC |
|---|---|---|
| A-mode1 | B-mode1 | 0.9881 |
| A-mode2 | B-mode2 | 0.9571 |
| A-mode3 | B-mode3 | 0.9937 |
| A-mode4 | B-mode5 | 0.9996 |
| A-mode5 | B-mode6 | 0.9999 |
| A-mode6 | B-mode6 | 0.9999 |

### Reading the Global MAC Matrices -- Flagged Slots
`comac_mac_matrices.png` plots both tables above as heatmaps (same palette used throughout this report). Two markings on it are not explained on the figure itself (no legend was added, to keep it uncluttered) and are the direct input to Stage 3's Step 2, so they are recorded here instead:

- **Green border** -- the row-wise argmax (naive best match per Set-A mode), i.e. the *Best-matching Set B mode* rows tabulated above -- not a one-to-one assignment, just the single highest value in that row.
- **Black cell, boxed value flanked by X** -- a slot flagged as unreliable and excluded from the Stage 3 pairing (below), even where it happens to also be the row's argmax. Two distinct reasons:

  1. *(A-mode2, B-mode1), both sub-MACs.* Rotationally this is A-mode2's argmax (0.9127), but the runner-up in the same row, B-mode2 (0.8840), is barely behind it -- and B-mode2 is also A-mode2's actual best match translationally (0.9571 vs. only 0.6931 for B-mode1, and B-mode2's 909 Hz sits far closer to A-mode2's 746 Hz than B-mode1's 391 Hz does). Splitting the mass-orthonormal 6-DOF eigenvector into rotational and translational halves discards the $M$-weighted orthogonality that keeps unrelated modes from cross-correlating in the full vector, so the rotational sub-MAC's slight edge for B-mode1 reads as leakage from that lost weighting, not a genuine shape match. B-mode2 is the physically consistent pairing; B-mode1 is not used for A-mode2 in either sub-system.
  2. *(A-mode3..6, B-mode5/B-mode6), translational sub-MAC only.* Every entry in that 2x2-plus block is $\geq 0.93$ (see the Linear/Translational Global MAC table above) -- with only 2 elements ($x_s$, $x_n$) to distinguish 4 higher-frequency modes, the translational sub-vector cannot separate B-mode5 from B-mode6, or tell A-mode4/5/6 apart from each other against either. That is spatial aliasing (too few coordinates for the number of modes being discriminated), not a resolvable correspondence, so A-mode4/5/6 are dropped from the translational pairing entirely rather than assigned an arbitrary winner.

## Co-MAC Stage 3 -- Coordinate Modal Assurance Criterion (COMAC)
Global MAC (previous section) pairs whole *modes* across the two systems; COMAC goes the other way -- for a fixed set of already-paired modes, it correlates one *coordinate* (DOF) at a time, across all those modes, to show which specific rows of the drivetrain are most disrupted by adding friction. The standard definition, for DOF $i$ over $L$ paired modes $(l \to k_l)$:

$$\mathrm{COMAC}_i = \frac{\left|\sum_{l=1}^{L}(\phi_A)_{i,l}\,(\phi_B)_{i,k_l}\right|^2}{\sum_{l=1}^{L}(\phi_A)_{i,l}^2\;\cdot\;\sum_{l=1}^{L}(\phi_B)_{i,k_l}^2}\;\in[0,1]$$

**Step 1 -- Partition rows**: the rotational (4-row) and translational (2-row) blocks already split out above.

**Step 2 -- Filter and pair columns**: take each sub-MAC's row-wise argmax (Torsional/Linear-Translational Global MAC tables above) and drop any pair that lands in a flagged/aliased slot (see `comac_mac_matrices.png`) rather than force a match. This also requires fixing the relative sign of each pair before summing -- `eigh` returns each mode with an arbitrary overall sign, and COMAC sums $(\phi_A)_{i,l}(\phi_B)_{i,k_l}$ **before** squaring, so an unresolved sign flip would let a genuinely good match cancel itself out. The sign used per pair (from the full 6-DOF dot product, applied consistently to both sub-systems) is listed alongside the pairs below.

*Rotational pairs* (5 of 6 Set-A modes; A-mode2 dropped -- its only candidate, B-mode1, is the flagged subspace-leakage slot):
| Set A mode | Set B mode | sign applied | MAC value |
|---|---|---|---|
| A-mode1 | B-mode1 | + | 0.8698 |
| A-mode3 | B-mode3 | + | 0.7827 |
| A-mode4 | B-mode5 | - | 0.6807 |
| A-mode5 | B-mode6 | + | 0.9869 |
| A-mode6 | B-mode6 | + | 0.9481 |

*Translational pairs* (3 of 6 Set-A modes; A-mode4/5/6 dropped -- their only candidates, B-mode5/B-mode6, fall inside the flagged aliased block):
| Set A mode | Set B mode | sign applied | MAC value |
|---|---|---|---|
| A-mode1 | B-mode1 | + | 0.9881 |
| A-mode2 | B-mode2 | - | 0.9571 |
| A-mode3 | B-mode3 | + | 0.9937 |

**Steps 3-5 -- Numerator, denominator, and normalize**: for each DOF row, sum the sign-corrected cross-products over the paired modes and square it (numerator), separately sum the squared elements of Set A and of Set B over those same modes and multiply them (denominator), then divide:

### Rotational COMAC
| DOF | sum($\phi_A\phi_B$) | sum($\phi_A^2$) | sum($\phi_B^2$) | COMAC |
|---|---|---|---|---|
| **theta_m** | 5.7399e+05 | 1.1091e+06 | 4.2348e+05 | 0.7015 |
| **theta_c** | 4.6423e+05 | 8.4627e+05 | 5.7860e+05 | 0.4401 |
| **theta_s** | 1.2650e+06 | 1.6495e+06 | 1.5918e+06 | 0.6095 |
| **theta_sb** | 7.4268e+06 | 6.6659e+06 | 1.3139e+07 | 0.6298 |

### Translational COMAC
| DOF | sum($\phi_A\phi_B$) | sum($\phi_A^2$) | sum($\phi_B^2$) | COMAC |
|---|---|---|---|---|
| **x_s** | 2.3351e-02 | 1.7841e+00 | 6.7460e-02 | 0.0045 |
| **x_n** | 1.6903e-01 | 2.1374e+00 | 2.4581e+00 | 0.0054 |

Reading these: $\theta_{sb}$ (support bearing ring, where the LuGre bearing friction torque acts) and $\theta_c$ come out lowest among the rotational DOFs (0.63 and 0.44), meaning their relative modal contribution is the most reshaped by the added friction. Both translational DOFs collapse to near zero (<0.01) -- across the only 3 modes that could be validly paired, $x_s$ and $x_n$ barely correlate row-wise at all between the two systems, even though the *combined* 2-vector direction was highly correlated in the Global MAC table above; that combined score is apparently carried by whichever of $x_s$/$x_n$ has the larger raw magnitude in each mode, not by both coordinates moving consistently together across modes. With only $L=3$ pairs this translational result rests on a small sample and should be read as indicative, not conclusive. See `comac_stage3.png` for these six values plotted against the conventional COMAC >= 0.9 "unaffected" threshold.

## Mass-Normalization Check
`eigh(K, M)` with `M` passed as the second (generalized) argument returns eigenvectors satisfying `phi.T @ M @ phi = I` by construction (LAPACK's generalized symmetric-definite driver). Verified numerically below, not just cited.

### Set A: `phi0.T @ M0 @ phi0`
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **mode1** | 1.00e+00 | -1.35e-17 | -3.02e-16 | -1.78e-17 | -3.67e-18 | -5.15e-17 |
| **mode2** | -9.96e-18 | 1.00e+00 | 1.48e-18 | -4.81e-18 | -1.04e-16 | 2.13e-17 |
| **mode3** | -3.37e-16 | 3.70e-18 | 1.00e+00 | -4.56e-17 | -7.13e-17 | 9.50e-17 |
| **mode4** | -1.03e-17 | -5.56e-18 | -5.26e-17 | 1.00e+00 | -7.59e-17 | 1.74e-16 |
| **mode5** | -1.38e-17 | -1.19e-16 | -6.85e-17 | -5.23e-17 | 1.00e+00 | 2.85e-17 |
| **mode6** | -6.05e-17 | 2.17e-17 | 7.83e-17 | 2.34e-16 | 4.37e-17 | 1.00e+00 |

max \|off-identity\| = `4.441e-16`, diagonal range `[1.000000000, 1.000000000]` -> **PASS** (tolerance 1e-8)

### Set B: `phi1.T @ M1 @ phi1`
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **mode1** | 1.00e+00 | -8.88e-18 | 6.68e-18 | -7.81e-17 | -1.20e-18 | 2.43e-19 |
| **mode2** | -8.48e-18 | 1.00e+00 | -4.77e-17 | 1.68e-18 | -1.60e-16 | -9.58e-18 |
| **mode3** | 7.46e-18 | -3.10e-17 | 1.00e+00 | 5.71e-18 | 3.67e-16 | -2.54e-17 |
| **mode4** | -7.90e-17 | 1.11e-18 | 6.54e-18 | 1.00e+00 | 7.41e-18 | -7.68e-19 |
| **mode5** | -1.17e-18 | -1.51e-16 | 3.64e-16 | 7.40e-18 | 1.00e+00 | 8.24e-18 |
| **mode6** | 2.80e-19 | -3.13e-18 | -2.02e-17 | -6.20e-19 | 2.61e-17 | 1.00e+00 |

max \|off-identity\| = `6.661e-16`, diagonal range `[1.000000000, 1.000000000]` -> **PASS** (tolerance 1e-8)

## Notes
- `M` differs in scale per DOF between rotational (`I_ii`, kg·m², ~1e-7-1e-6) and translational (`M_ii`, kg, ~0.1-0.4) rows -- mass-normalization is w.r.t. **this** `M`, not the identity, so raw `|phi|` entries are not directly comparable across DOF types without the `lead_ratio` rescaling `plot_mode_shapes(_linearized).py` apply for readability. Not applied here deliberately -- this report shows the extraction as `scipy` actually returns it.
- Sign is arbitrary (`phi_j` and `-phi_j` are equally valid solutions); not fixed here, unlike the mode-shapes plots.
- Set A and Set B mode indices are **not** claimed to correspond to each other -- e.g. Set A mode 1 (176.7 Hz) is not asserted to be "the same mode" as Set B mode 1 (391.5 Hz). Establishing that correspondence is a MAC computation, not done here.
- The two MAC matrices above are independent of each other by construction -- a mode can correlate strongly in $\mathrm{MAC}_{\text{rot}}$ and weakly in $\mathrm{MAC}_{\text{trans}}$ (or vice versa) since they draw on disjoint DOF subsets. Neither is a substitute for a full-vector MAC; combining them (e.g. a weighted product) would be a further step beyond this Stage 1/2 extraction, not attempted here.

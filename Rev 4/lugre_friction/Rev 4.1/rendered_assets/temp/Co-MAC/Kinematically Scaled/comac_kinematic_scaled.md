# Co-MAC (Kinematically Scaled)
Alternative to the [Split Method](../Split%20Method/comac_mode_extraction.md), which computed two separate MAC matrices (4-DOF rotational, 2-DOF translational) to sidestep the unit mismatch between $\theta$ (rad) and $x$ (m). Here the mismatch is removed at the source instead: the 4 rotational rows are rescaled to equivalent linear displacement using the lead-screw kinematic relationship, then all 6 DOFs are carried through MAC, mode-pairing, sign-alignment, and COMAC together, in one pass.

Set A (frictionless baseline) and Set B (frozen-linearized LuGre, `V_STAGE = 5 mm/s`) are the identical two `scipy.linalg.eigh(K, M)` solves used in the Split Method -- not reproduced in full here (see that report for the raw $\lambda$, $\omega$, $f$, and unscaled $\phi$ tables). This report starts from those same raw, unfixed-sign, mass-normalized eigenvectors $\phi_0$ (Set A) and $\phi_1$ (Set B), state order `theta_m, theta_c, theta_s, theta_sb, x_s, x_n`.

## Step 1 -- Kinematic Scaling of the Rotational Rows
The lead screw ideal kinematic relationship (ball-screw/nut, no slip) ties axial travel to screw rotation by $x = \dfrac{L}{2\pi}\theta$, where $L$ is the screw lead. Both models share the same physical lead screw, so both use the same $L$:

$$L = 1.000e-03\ \text{m}, \qquad R = \frac{L}{2\pi} = 1.591549e-04\ \text{m/rad}$$

Multiply rows 1-4 ($\theta_m, \theta_c, \theta_s, \theta_{sb}$) of both raw eigenvector matrices by $R$; rows 5-6 ($x_s, x_n$) are already in meters and are left untouched:

$$\vec{\phi}_{A,\text{scaled},l} = \begin{bmatrix}R\,\theta_m \\ R\,\theta_c \\ R\,\theta_s \\ R\,\theta_{sb} \\ x_s \\ x_n\end{bmatrix}_{A,l}, \qquad \vec{\phi}_{B,\text{scaled},k} = \begin{bmatrix}R\,\theta_m \\ R\,\theta_c \\ R\,\theta_s \\ R\,\theta_{sb} \\ x_s \\ x_n\end{bmatrix}_{B,k}$$

Why this works: every one of the 6 rows now reports the same physical quantity (equivalent linear displacement, m), in the same numerical range, so a 6-element dot product weighs all 6 DOFs on their actual structural contribution instead of letting whichever DOF group has larger raw magnitude dominate.

### Scaled eigenvector matrix phi0_scaled -- Set A
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **theta_m** | 9.1908e-02 | -7.1619e-03 | 1.3877e-01 | 1.9692e-02 | -1.1626e-04 | -5.0706e-04 |
| **theta_c** | 9.5237e-02 | -5.4773e-03 | -4.9269e-02 | -9.8961e-02 | 2.4445e-03 | 1.1811e-02 |
| **theta_s** | 9.5666e-02 | -4.2555e-03 | -8.0917e-02 | 1.1933e-01 | -1.9778e-02 | -1.0701e-01 |
| **theta_sb** | 9.5750e-02 | -4.3230e-03 | -8.7613e-02 | 1.7812e-01 | 9.9400e-02 | 3.3227e-01 |
| **x_s** | 4.4415e-03 | 1.3330e+00 | 8.4478e-02 | -1.3429e-01 | 3.2621e+00 | -8.6556e-01 |
| **x_n** | 1.0061e-01 | 1.4585e+00 | 6.3072e-03 | 1.6992e-02 | -5.5609e-01 | 1.4890e-01 |

### Scaled eigenvector matrix phi1_scaled -- Set B
| | mode1 | mode2 | mode3 | mode4 | mode5 | mode6 |
|---|---|---|---|---|---|---|
| **theta_m** | 5.3307e-04 | 1.3196e-01 | 1.0324e-01 | 2.3160e-03 | -8.3158e-03 | -2.1047e-05 |
| **theta_c** | 5.1876e-04 | 8.2482e-02 | -9.9854e-02 | -2.9460e-03 | 6.8422e-02 | 1.3163e-03 |
| **theta_s** | 4.9654e-04 | 5.1306e-02 | -7.9875e-02 | -1.7352e-03 | -1.7784e-01 | -3.4011e-02 |
| **theta_sb** | 1.4755e-04 | 1.5331e-02 | -2.4502e-02 | -5.3498e-04 | -6.2058e-02 | 4.0519e-01 |
| **x_s** | 2.4222e-01 | -1.0809e-02 | 9.3125e-02 | -3.6229e+00 | -8.8830e-03 | -1.9189e-04 |
| **x_n** | 1.5678e+00 | -7.7436e-03 | -4.4371e-04 | 1.0483e-01 | 1.3133e-03 | 3.4816e-05 |

## Step 2 -- Single 6-DOF Global MAC
With every row in the same units, one standard MAC matrix now covers all 6 DOFs (same formula as the Split Method, applied to the scaled vectors instead of a sub-vector):

$$\mathrm{MAC}_{\text{scaled}}(l,k) = \frac{\left|\vec{\phi}_{A,\text{scaled},l}^{\mathsf{T}}\vec{\phi}_{B,\text{scaled},k}\right|^2}{\left(\vec{\phi}_{A,\text{scaled},l}^{\mathsf{T}}\vec{\phi}_{A,\text{scaled},l}\right)\left(\vec{\phi}_{B,\text{scaled},k}^{\mathsf{T}}\vec{\phi}_{B,\text{scaled},k}\right)}$$

### MAC_scaled (rows: Set A, cols: Set B)
| | B-mode1 | B-mode2 | B-mode3 | B-mode4 | B-mode5 | B-mode6 |
|---|---|---|---|---|---|---|
| **A-mode1** | 0.2184 | 0.5201 | 0.0557 | 0.0001 | 0.1580 | 0.1674 |
| **A-mode2** | 0.6930 | 0.0070 | 0.1081 | 0.4264 | 0.0005 | 0.0000 |
| **A-mode3** | 0.0085 | 0.0517 | 0.8162 | 0.1641 | 0.1221 | 0.1513 |
| **A-mode4** | 0.0002 | 0.0105 | 0.0778 | 0.2438 | 0.4809 | 0.3757 |
| **A-mode5** | 0.0002 | 0.0031 | 0.2307 | 0.9796 | 0.0023 | 0.0009 |
| **A-mode6** | 0.0003 | 0.0031 | 0.2049 | 0.8463 | 0.0014 | 0.1298 |

See `comac_kinematic_mac_matrix.png` for this table as a heatmap (green border = row-wise argmax, same convention as the Split Method figure).

## Step 3 -- Master Pairing List
Row-wise argmax of MAC_scaled -- one Set-B mode picked per Set-A mode:
| Set A mode | Set B mode | MAC |
|---|---|---|
| A-mode1 | B-mode2 | 0.5201 |
| A-mode2 | B-mode1 | 0.6930 |
| A-mode3 | B-mode3 | 0.8162 |
| A-mode4 | B-mode5 | 0.4809 |
| A-mode5 | B-mode4 | 0.9796 |
| A-mode6 | B-mode4 | 0.8463 |

**Not a clean bijection**: B-mode4 is claimed as the best match by more than one Set-A mode above (A-mode5 and A-mode6 both peak at B-mode4 here). Preserving all 6 DOFs together removes the severe aliasing seen in the Split Method's translational block (that whole block was >= 0.93 with no separation at all), but it does not by itself guarantee a strict one-to-one assignment -- resolving that fully would need an assignment algorithm (e.g. Hungarian/linear-sum-assignment) on top of the row-wise argmax, which is not run here. Both rows are still carried through Steps 4-5 below as given, matching how the Split Method handled its own duplicate (B-mode6 claimed by two rows in the rotational sub-MAC).

## Step 4 -- Sign Alignment
For each pair (l, k) above, compute S = the dot product of the full 6-element scaled vectors phi_A_scaled_l and phi_B_scaled_k. If S is negative, eigh returned mode k of Set B with the opposite overall sign from mode l of Set A; multiply column k of phi1_scaled by -1 before using it in Step 5 (applied per pair, not per column, since B-mode4 above is reused by two different pairs and is not guaranteed to need the same sign both times):
| Set A mode | Set B mode | S | sign applied |
|---|---|---|---|
| A-mode1 | B-mode2 | 2.5532e-02 | + |
| A-mode2 | B-mode1 | 2.6096e+00 | + |
| A-mode3 | B-mode3 | 3.5720e-02 | + |
| A-mode4 | B-mode5 | -3.7994e-02 | - |
| A-mode5 | B-mode4 | -1.1876e+01 | - |
| A-mode6 | B-mode4 | 3.1514e+00 | + |

## Step 5 -- Unified 6-DOF COMAC
One loop over all 6 coordinates (no rotational/translational split), summing over all L=6 master-paired, sign-aligned modes at once:

$$\mathrm{COMAC}_i = \frac{\left(\sum_{m=1}^{L}(\phi_{A,\text{scaled}})_{i,l_m}\,(\phi_{B,\text{scaled}})_{i,k_m}\right)^2}{\sum_{m=1}^{L}(\phi_{A,\text{scaled}})_{i,l_m}^2 \cdot \sum_{m=1}^{L}(\phi_{B,\text{scaled}})_{i,k_m}^2}\;\in[0,1]$$

### Unified COMAC
| DOF | sum(phiA*phiB) | sum(phiA^2) | sum(phiB^2) | COMAC |
|---|---|---|---|---|
| **theta_m** | 2.6613e-02 | 2.8145e-02 | 2.8150e-02 | 0.8940 |
| **theta_c** | 1.9516e-02 | 2.1466e-02 | 2.1473e-02 | 0.8263 |
| **theta_s** | 3.2742e-02 | 4.1799e-02 | 4.0645e-02 | 0.6310 |
| **theta_sb** | 1.4543e-02 | 1.6887e-01 | 4.6872e-03 | 0.2672 |
| **x_s** | 1.5283e+01 | 1.3193e+01 | 2.6318e+01 | 0.6728 |
| **x_n** | 2.3598e+00 | 2.4691e+00 | 2.4801e+00 | 0.9094 |

Reading these against the Split Method result: theta_sb (support bearing ring) is clearly the lowest of all 6 DOFs at 0.267 -- and it is the one DOF the LuGre bearing friction torque acts on directly, so this is the cleanest physical signal either method has produced. theta_s is next-lowest (0.631), with x_s close behind (0.673). Both translational DOFs move well away from the Split Method's near-zero reading (x_s = 0.673, x_n = 0.909) -- consistent with that earlier near-zero result being an artifact of leaving the eigenvector sign unresolved before summing (Step 4 above is exactly the fix), not a genuine finding about x_s/x_n. x_n and theta_m score highest (0.909, 0.894), i.e. least reshaped by the added friction; theta_c sits in between (0.826). See `comac_kinematic_comac.png` for these six values plotted against the same COMAC >= 0.9 threshold used in the Split Method figure.

## Notes
- R is a single scalar shared by both models (same physical lead screw, same L), applied identically to Set A and Set B -- it rescales units, it does not introduce any new information relating the two systems.
- This MAC_scaled is not comparable cell-by-cell with the Split Method's MAC_rot/MAC_trans: those were computed on sub-vectors with theta still in raw radians (mass-normalized units, not physically scaled), so the resulting best-match pairs differ -- each method's pairing is only self-consistent within itself.
- As in the Split Method, sign is otherwise arbitrary per raw mode (eigh guarantees phi.T @ M @ phi = I, not a consistent sign); Step 4 above resolves it only for the 6 pairs actually used in Step 5, not for the full raw phi0/phi1 matrices.

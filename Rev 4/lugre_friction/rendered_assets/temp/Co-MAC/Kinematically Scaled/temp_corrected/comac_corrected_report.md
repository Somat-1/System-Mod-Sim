# Corrected Kinematically Scaled Co-MAC Audit

This directory is a separate corrected result set. The original Kinematically
Scaled files are unchanged.

## Corrections implemented

1. **Coordinate-invariant MAC metric.** The two systems have exactly the same
   physical mass matrix. Global mode correlation therefore uses the common-mass
   inner product instead of Euclidean MAC after an arbitrary coordinate scaling:

   $$MAC_M(a,b)=\frac{|\phi_{A,a}^T M\phi_{B,b}|^2}
   {(\phi_{A,a}^T M\phi_{A,a})(\phi_{B,b}^T M\phi_{B,b})}.$$

   This is invariant to a consistent coordinate transformation and does not
   claim that equal units imply equal DOF weighting.

2. **Bijective pairing.** `linear_sum_assignment(..., maximize=True)` selects
   one unique Set-B mode for every Set-A mode. No mode column is reused.

3. **Consistent polarity.** Each unique pair receives one sign from
   $\operatorname{sign}(\phi_A^T M\phi_B)$. The standard signed-sum COMAC is
   then evaluated. The per-term-absolute result is reported only as a sensitivity
   comparison, not as the primary definition.

4. **COMAC is evaluated on raw coordinates.** For fixed pairs and signs, scaling
   a row in both datasets cancels identically from COMAC. The lead-screw scale is
   relevant to the superseded Euclidean pairing method, not to the row-wise
   COMAC normalization itself.

## Mass-weighted Global MAC

| Set A / Set B | B-mode1 | B-mode2 | B-mode3 | B-mode4 | B-mode5 | B-mode6 |
|---|---|---|---|---|---|---|
| A-mode1 | 0.0048 | 0.8513 | 0.0916 | 0.0000 | 0.0275 | 0.0249 |
| A-mode2 | 0.9030 | 0.0043 | 0.0003 | 0.0923 | 0.0000 | 0.0001 |
| A-mode3 | 0.0000 | 0.1253 | 0.8211 | 0.0000 | 0.0318 | 0.0218 |
| A-mode4 | 0.0001 | 0.0157 | 0.0771 | 0.0022 | 0.7998 | 0.1051 |
| A-mode5 | 0.0860 | 0.0001 | 0.0012 | 0.8450 | 0.0028 | 0.0649 |
| A-mode6 | 0.0061 | 0.0034 | 0.0087 | 0.0604 | 0.1381 | 0.7833 |

See `comac_corrected_mac_matrix.png`. Green borders show the optimal bijection.

## Optimal bijective mode pairs

| Set A | f_A [Hz] | Set B | f_B [Hz] | MAC_M | sign |
|---|---|---|---|---|---|
| A-mode1 | 176.686 | B-mode2 | 909.281 | 0.8513 | + |
| A-mode2 | 745.962 | B-mode1 | 391.480 | 0.9030 | + |
| A-mode3 | 1650.224 | B-mode3 | 1976.827 | 0.8211 | + |
| A-mode4 | 3429.349 | B-mode5 | 4237.584 | 0.7998 | - |
| A-mode5 | 6536.123 | B-mode4 | 2120.749 | 0.8450 | - |
| A-mode6 | 6863.445 | B-mode6 | 11093.079 | 0.7833 | + |

Total assigned mass-weighted MAC = `5.003532`.

## Corrected COMAC

| DOF | sum(phiA phiB) | sum(phiA^2) | sum(phiB^2) | signed COMAC | absolute-product check | original duplicate |
|---|---|---|---|---|---|---|
| theta_m | 1.050697e+06 | 1.111111e+06 | 1.111111e+06 | 0.894211 | 0.894724 | 0.893961 |
| theta_c | 7.724396e+05 | 8.474576e+05 | 8.474576e+05 | 0.830794 | 0.831276 | 0.826257 |
| theta_s | 1.428948e+06 | 1.650165e+06 | 1.650165e+06 | 0.749856 | 0.752878 | 0.630995 |
| theta_sb | 5.896194e+06 | 6.666667e+06 | 6.666667e+06 | 0.782215 | 0.782228 | 0.267206 |
| x_s | 1.214780e+01 | 1.319261e+01 | 1.319261e+01 | 0.847878 | 0.848225 | 0.672765 |
| x_n | 2.344202e+00 | 2.469136e+00 | 2.469136e+00 | 0.901364 | 0.902601 | 0.909356 |

The corrected result does not support interpreting the original low
`theta_sb = 0.267` value as a direct friction signature. With unique paired
modes it becomes approximately `0.782`. The result should still be interpreted
as coordinate correlation conditional on this mode map, rather than as a direct
damage or friction-localization measurement.

## Scaling audit

The superseded kinematic method used
`R = 1.591549430919e-04 m/rad` and `R^2 = 2.533029591058e-08`. Its fraction of scaled Euclidean modal norm in the two
translational rows was:

| Mode | Set A translational fraction | Set B translational fraction |
|---|---|---|
| mode1 | 0.220579 | 1.000000 |
| mode2 | 0.999970 | 0.006486 |
| mode3 | 0.166558 | 0.239033 |
| mode4 | 0.246035 | 0.999999 |
| mode5 | 0.999062 | 0.002000 |
| mode6 | 0.863447 | 0.000000 |

These fractions show that equalizing units did not equalize modal weighting.
The corrected Global MAC therefore uses the common physical mass metric.

## Split Method sign clarification

The Split Method already applied signs `(+,-,+)` to its three translational
pairs. Its near-zero `x_s` and `x_n` results were caused by the selected
three-pair amplitude correspondence and the small sample size, not by unresolved
eigenvector signs. No claim to the contrary is carried into this result set.

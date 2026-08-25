# Rev 4.2 — Jacobian-injected LuGre friction

Rev 4.2 implements three LuGre ports without hand-written reaction signs.
For every port, `v_p = J_p qdot` and the mechanical equation receives
`-J_p.T F_p`. The guideway and support-bearing ports are grounded. The nut
port spans the corrected screw/nut relative coordinate.

Unlike Rev 4.1, the baseline `k_nut` and `c_nut` contact path remains in the
structural matrices. The nut LuGre element is therefore parasitic pre-rolling
drag in parallel with the load-bearing rolling contact; it cannot cap the
transmitted thrust at the LuGre breakaway force.

The state is `[q(6), qdot(6), z_way, z_nut, z_sb]`. The LuGre law uses
`sqrt(v^2 + epsilon^2)` instead of `abs(v)`, and the model supplies an
analytical 15-state Jacobian. `sigma1` is sized with
`2*zeta*sqrt(sigma0*m_eff)` for `zeta=0.7`. The motor detent is evaluated as
`T_d*sin(4*N_r*theta_m)` rather than embedded as `k_d` in the structural K
matrix.

## Bode comparison

`scripts/build_bode_rev42.py` linearizes the complete 15-state model at a
frozen 5 mm/s cruise point and compares `x_n/theta_cmd` with the Rev 4
frictionless baseline over 0–8 kHz. It saves the overlay, response arrays and
a JSON summary under `rendered_assets`.

## Parameter caveats

The friction levels are reused from Rev 4.1 rather than newly identified.
The baseline `c_nut=101 N*s/m` is also retained, even though LuGre
`sigma1/sigma2` now supplies parallel dissipation. That split must be fitted
before the combined damping is treated as physical rather than provisional.

## Nonlinear stepping

`scripts/generate_stepping_rev42.py` mirrors the Rev 4 full-step/16x and
fast/settled command cases. It uses piecewise Radau integration with the
analytical Jacobian, generates the montage, single-step diagnostic and axial
spectrum, and adds a four-case tracking-error overlay against the frictionless
baseline. The accepted trajectory's interface power and cumulative work are
saved with convergence and solver statistics in `stepping_rev42_summary.json`.

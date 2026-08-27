# Rev 4 Frictionless-Plant Backlog

This backlog tracks only the structural, frictionless Rev 4 reference plant
derived in `state_space_6dof.md` and implemented by
`scripts/build_bode_rev4.py`. Nonlinear friction findings are maintained with
the model that produced them:

- `lugre_friction/Rev 4.1/README.md` documents the replacement approach, in
  which the nut LuGre element replaces `k_nut` and `c_nut`.
- `lugre_friction/Rev 4.2/README.md` documents the parallel approach, in
  which `k_nut` and `c_nut` remain structural and LuGre adds a separate
  pre-rolling friction branch.

The root Rev 4 model is the common comparison baseline for both approaches;
it does not select between them.

---

## Resolved — 2026-08-18

### 1. Screw–nut structural coupling had the wrong sign

**Issue.** After assembling the linear screw–nut contact directly in
`K`/`C`, the command-to-stage transfer function
`x_n(s)/theta_cmd(s)` had a negative DC gain. A unit `theta_cmd` produced
`x_n = -1.326e-4`, while the ideal relation
`x_n = (L/2*pi)*theta_s` requires the same magnitude with positive sign.
The Bode phase consequently started at 180 degrees rather than zero.

**Cause.** The four off-diagonal `(L/2*pi)*k_nut` and `c_nut` cross-terms
between `theta_s` and `x_s`/`x_n` represented

```text
x_n - x_s + (L/2*pi)*theta_s
```

instead of the convention-consistent structural deformation

```text
x_n - x_s - (L/2*pi)*theta_s.
```

Positive screw rotation must drive positive stage translation.

**Fix.** Negated the four cross-terms in `K` and `C` and made the same sign
correction in the `F_nut` constitutive relation. A direct solution of
`K*q = B_u*theta_cmd` then changed `x_n` from `-1.326e-4` to
`+1.326e-4`, matching the ideal ratio. The correction is a coordinate-sign
similarity transformation, so it does not alter the modal eigenvalues.

**Where:** `state_space_6dof.md` Sections 3.4, 4.2, 5.3 and 5.4;
`scripts/build_bode_rev4.py`, `build_matrices()`.

---

### 2. Mode 1 was almost undamped

**Issue.** The 176.7 Hz mode had `zeta_1 = 4.89e-5`, producing an extremely
narrow resonance and an approximately 74 s two-percent settling time.
Frequency-grid changes therefore produced unstable-looking peak estimates,
and the settling time was impractical for the stepping simulations.

**Cause.** In this mode, all four rotational coordinates move almost in
unison against the electromagnetic spring `k_EM`. The local dampers `c_c`,
`c_s1` and `c_s2` act only on relative velocity, which is nearly zero for
that mode shape. Increasing those dampers therefore has little leverage on
the first mode.

**Fix.** Added `c_EM = 1.3283e-4 N m s/rad` to `C[0,0]` as viscous damping
of `theta_m` against a fixed frame. It was selected by bisection to give the
first mode a damping ratio of 0.02. The `theta_cmd_dot` half of a literal
relative electromagnetic damper is intentionally omitted because the
stepping command is a staircase and its derivative would be an impulse
train on the finite simulation grid.

| Mode | Frequency | Damping before | Damping after |
|---|---:|---:|---:|
| 1 | 176.69 Hz | 4.89e-5 | **0.020000** |
| 2 | 745.97 Hz | 7.24e-3 | 7.27e-3 |
| 3 | 1650.23 Hz | 1.69e-2 | 2.18e-2 |
| 4 | 3429.32 Hz | 1.99e-2 | 1.99e-2 |
| 5 | 6536.22 Hz | 2.40e-2 | 2.40e-2 |
| 6 | 6863.23 Hz | 2.19e-2 | 2.19e-2 |

The first-mode settling time fell from approximately 74 s to 176 ms, which
sets the scale for the 250 ms settled dwell in
`scripts/generate_stepping_trajectory.py`.

**Where:** `model_parameters.json`, `state_space_6dof.md` Section 5.4 and
`scripts/build_bode_rev4.py`, `build_matrices()`.

---

## Open / next — frictionless plant

- Identify the placeholder structural parameters, especially `k_nut`,
  `c_nut`, `k_s1`, `k_s2` and their damping counterparts.
- Confirm whether `c_EM` represents a physical motor loss or should remain
  an explicitly numerical/modal closure.
- Confirm the screw lead and geometric datum used for the numerical
  evaluation; the BOM indicates `L = 1 mm`.
- Retain the root model as a zero-friction comparison case when either
  LuGre sub-revision is evaluated.

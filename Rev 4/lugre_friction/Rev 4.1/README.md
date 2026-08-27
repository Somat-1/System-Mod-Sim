# Rev 4.1 — LuGre replacement approach

Sub-branch of `Rev 4/`. Does not modify anything in the parent directory
(`state_space_6dof.md`, `model_parameters.json`, `build_bode_rev4.py`,
`generate_stepping_trajectory.py` are all untouched). The parent is retained
as the common frictionless comparison plant. This directory is self-contained:
its own `model_parameters.json` is a standalone copy, not an import.

## Model scope

Rev 4.1 is the **replacement experiment**. At the screw–nut interface, the
LuGre element replaces the structural `k_nut` and `c_nut` branch. The root
Rev 4 documents do not prescribe this topology; they define only the common
frictionless comparison plant. Rev 4.2 tests the distinct alternative in
which the structural contact remains and LuGre acts in parallel.

## What changed vs. the frictionless baseline

The baseline (`../../model_parameters.json`, `../../scripts/build_bode_rev4.py`) embeds
the screw-nut interface as a linear structural spring/damper (`k_nut`,
`c_nut`) directly in `K`/`C`. The frictionless reference makes no claim
about slip; selecting a nonlinear interface interpretation is the purpose of
the separate Rev 4.1 and Rev 4.2 experiments.

This sub-branch removes that assumption. `k_nut` and `c_nut` are dropped
entirely because Rev 4.1 deliberately interprets the LuGre bristle as the
complete screw–nut contact compliance as well as its friction law. In their
place are three independent LuGre friction ports: screw–nut, support bearing
and guideway. This replacement is a Rev 4.1 modelling hypothesis, not a
general requirement and not the topology later adopted by Rev 4.2.

State vector grows from 12 to 15:

```
x = [theta_m, theta_c, theta_s, theta_sb, x_s, x_n,
     theta_m_dot, ..., x_n_dot,
     z_sb, z_nut, z_way]
```

Everything else from the baseline carries over unchanged: the rotational
chain (`k_c/c_c`, `k_s1/c_s1`, `k_s2/c_s2`, `k_EM/k_d/c_EM` on `theta_m`),
and the axial bearing preload spring `k_brg/c_brg` on `x_s` (not a friction
interface -- it's the thrust-bearing pair's own elastic/damping behavior,
distinct from the support bearing's *rotational* friction `T_fric,sb`).

## Relative velocities

```
v_sb  = theta_sb_dot
v_way = x_n_dot
v_nut = (L/2*pi)*theta_s_dot - (x_n_dot - x_s_dot)
```

## LuGre port (generic form, one per interface)

```
g(v)    = Fc + (Fs - Fc)*exp(-(v/vs)^2)
z_dot   = v - sigma0*|v|/g(v)*z
F       = sigma0*z + sigma1*z_dot + sigma2*v
```

`F` is `T_fric,sb`, `F_fric,way`, or `F_fric,nut` depending on which port's
`(v, z, sigma0, sigma1, sigma2, Fc, Fs, vs)` are used (`lugre_model.py`
`lugre_force()`).

## Sign convention: v_nut and the F_fric,nut injection (read this before
touching the equations of motion)

The baseline's K/C sign fix (backlog.md issue #1, 2026-08-18) established
`x_n - x_s - (L/2*pi)*theta_s` as the physically correct screw-nut coupling
variable -- the one where positive `theta_s` drives positive `x_n`, per
Sec. 1's stated convention. Call its time derivative the "closing velocity"
`u_dot = x_n_dot - x_s_dot - (L/2*pi)*theta_s_dot`.

`v_nut` as specified for this LuGre port is `(L/2*pi)*theta_s_dot -
(x_n_dot - x_s_dot)`, which is exactly `-u_dot`. Because a LuGre force is
odd in its velocity argument (reversing the input velocity reverses the
force), `F_fric,nut` computed from `v_nut` is the negative of what the
baseline's `F_nut` meant. Injecting it into the equations of motion with
the baseline's signs (`+F_nut` on `x_s`, `-F_nut` on `x_n`, from Sec. 4 of
`state_space_6dof.md`) would silently reintroduce a sign inversion.

Instead, the injection signs here are re-derived directly from `v_nut` by
virtual work, so they are self-consistent with the port as specified rather
than with the old labels. For a relative velocity `v = sum(c_i * qdot_i)`
and an interface force `F` that opposes it, the generalized force on
coordinate `q_i` is `Q_i = -F * c_i`. Here `c_theta_s = +L/2*pi`,
`c_x_s = +1`, `c_x_n = -1`, giving:

```
theta_s row: -(L/2*pi) * F_fric,nut   (same sign as Sec. 4 -- see note below)
x_s row:     -F_fric,nut               (flipped vs. Sec. 4's "+F_nut")
x_n row:     +F_fric,nut               (flipped vs. Sec. 4's "-F_nut")
```

The `theta_s` coefficient did not flip: both `v_nut`'s sign and `F_nut`'s
meaning flipped relative to the baseline, and those two flips cancel for
that one row. This was checked numerically (`lugre_model.py`, smoke test):
a positive commanded step produces positive `theta_s` and positive `x_n`
(correct macro direction), and slightly negative `x_s` (the screw
compressing into the thrust bearing, matching the physical description of
`F_fric,nut`'s reaction).

`T_fric,sb` and `F_fric,way` are unambiguous by comparison: each is driven
directly by its own DOF's velocity (`v_sb = theta_sb_dot`, `v_way =
x_n_dot`), so each is simply subtracted from that DOF's own equation, same
as `T_fric,sb` already was in the frictionless baseline before it was
frozen at zero.

## Parameters

See `model_parameters.json` `parameter_notes` for full provenance of every
value. Summary: mechanical parameters (inertias, the rotational-chain
stiffnesses/dampings, `k_brg/c_brg`, `c_EM`) are unchanged copies from the
baseline. The eighteen new LuGre parameters (`sigma0/1/2`, `Fc`or`Tc`,
`Fs`or`Ts`, `vs`, times three interfaces) are placeholders: `sigma0`,
`sigma1`, `sigma2`, and `vs` are reused from Rev 3's GMS friction parameter
sets (`n_*` for nut, `g_*` for guideway, `d_*` for support bearing by
elimination); `Fc`/`Fs`/`Tc`/`Ts` are derived or guessed where no workspace
reference existed.

One correction worth flagging explicitly: `sigma0_sb` was originally reused
directly from Rev 3's `d_sigma0 = 3.0e6 N m/rad`. With
`I_sb = 1.5e-7 kg m^2`, that gives a bristle natural frequency of roughly
712 kHz, far beyond every structural mode and severe enough to make the
solver trials fail. An intermediate trial used `500 N m/rad`; the executed
Rev 4.1 parameter was subsequently reduced to **`0.076 N m/rad`**. The
chronology matters: 500 is a superseded diagnostic value, not the current
Rev 4.1 model parameter. All of these values remain unmeasured placeholders.

## Solver

`scipy.integrate.solve_ivp` on the full nonlinear 15-state RHS
(`lugre_model.LuGreModel.rhs`), not `scipy.signal.lsim` (linear-only).
Tried Radau, BDF, and LSODA on a representative case (A=5e-3 rad,
omega=1 rad/s -- the slowest, most stiction-cycling point in the sweep).
Radau and BDF both diverge to NaN in under a second, unaffected by
tightening `atol`/`rtol` or capping `max_step`. LSODA completes reliably.
This is consistent with LuGre's non-smooth kink at `v=0` (the
`|v|/g(v)*z` term is not differentiable there) defeating Radau/BDF's
Newton-iteration and finite-difference-Jacobian machinery, while LSODA's
on-the-fly explicit/implicit switching tolerates it. `run_sinusoidal_sweep.py`
defaults to `method="LSODA"`; revisiting Radau/BDF would need an analytical
Jacobian (or an event-based re-start at `v=0` crossings) rather than just
tolerance tuning.

## Dynamic Local Linearization (`run_local_linearization_bode.py`)

An alternative to the nonlinear time-domain sweep above: freeze each
bristle at its steady-state deflection for a chosen cruising speed, and
take a Jacobian to collapse the LuGre port into a static equivalent
stiffness/damping pair, then run a proper Laplace-domain Bode (same
approach as `../../scripts/build_bode_rev4.py`) on the resulting linear model.

Substituting `z_dot` into `F` gives `F` as an explicit function of `(z,v)`
alone: `F(z,v) = sigma0*z*(1 - sigma1*|v|/g(v)) + v*(sigma1+sigma2)`. The
Jacobian of this at the frozen operating point `(z0, v0)` -- `K_eq =
dF/dz`, `C_eq = dF/dv`, computed numerically (finite difference) rather
than by hand-derived formula -- is the "frozen bristle" limit: since the
bristle doesn't have time to relax for fast/structural-frequency
perturbations, it behaves like a fixed elastic element, giving a genuine
static stiffness+damping pair rather than the frequency-dependent transfer
function a full (non-frozen) linearization of the coupled `(z,v)` dynamics
would produce.

**All three interfaces have `sigma1=0`** in this sub-branch's parameters,
which makes `F(z,v)` already exactly linear -- so `K_eq = sigma0` and
`C_eq = sigma2` come out **independent of the chosen operating velocity
v0**. Disclosed rather than hidden: the cruising-speed choice is still
conceptually necessary (it's what justifies treating the port as "gross
sliding" rather than presliding in the first place, and it sets the steady
bias force `F0`), but for this specific parameter set it doesn't change
the numbers. It would matter if `sigma1` were ever made nonzero.

**Operating velocities**, derived from one chosen stage speed
`V_STAGE = 5 mm/s` (placeholder): `v0_way = V_STAGE`; `v0_sb = V_STAGE *
2*pi/L` (screw/bearing rotation rate under ideal no-slip tracking); `v0_nut
= 0`. The last one is deliberate, not an oversight: `v_nut` is defined as a
slip/tracking-error term (`(L/2*pi)*theta_s_dot - (x_n_dot-x_s_dot)`), not
a bulk speed, and at an idealized no-slip steady cruise it is exactly zero.
This isn't a singularity -- it just means the nut's equivalent
stiffness/damping reduce to its raw presliding values.

**A parameter red flag surfaced by this run**: the steady bias torque at
the support bearing, `F0_sb = sigma0_sb*z0 + sigma2_sb*v0_sb`, comes out to
~14.1 N*m at `v0_sb` = 31.4 rad/s -- about 235x the motor's own holding
torque (0.06 N*m), dominated by `sigma2_sb*v0_sb = 0.45*31.4`. Not
physically plausible for a real bearing. It doesn't affect the Bode plot
itself (`F0` is a constant bias, not part of the linearized dynamics
matrix), but `sigma2_sb = 0.45 N*m*s/rad` (the other value reused from Rev
3's `d_sigma2`) is likely too large once actually evaluated at a rotational
operating speed, the same kind of cross-domain-reuse issue that made
`sigma0_sb` numerically fatal earlier (see the Solver section above). Not
fixed here since it doesn't corrupt this particular result; worth
revisiting before trusting `F0_sb` for anything.

**Result**: the DC magnitude of this locally-linearized Bode
(~-114 dB) closely matches the small-amplitude nonlinear sweep's
stiction-trapped value (-113.83 dB from `sinusoidal_bode.svg`) -- a
reassuring cross-check between the two independently-built models, since
both are dominated by the same frozen/presliding stiffness at low
frequency.

## Files

- `model_parameters.json` -- self-contained parameter set (baseline copy +
  LuGre additions, `k_nut`/`c_nut` dropped).
- `lugre_model.py` -- `LuGreModel`: builds the linear baseline (everything
  except the three friction ports) once, and evaluates the full nonlinear
  RHS per `solve_ivp` call.
- `run_sinusoidal_sweep.py` -- the nonlinear time-domain sweep described
  above; writes `rendered_assets/sinusoidal_sweep_montage.svg` (+`.png`)
  and the raw sweep data (`sinusoidal_sweep_data.npz`).
- `plot_sinusoidal_bode.py` -- standalone magnitude/phase Bode chart built
  from the saved sweep data (`rendered_assets/sinusoidal_bode.svg`).
- `run_local_linearization_bode.py` -- the Dynamic Local Linearization
  method described above; writes `rendered_assets/local_linearization_bode.svg`
  (+`.png`) and `local_linearization_bode_data.npz`.

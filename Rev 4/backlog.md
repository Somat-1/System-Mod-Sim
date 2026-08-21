# Rev 4 Backlog

Tracks issues found during model build-out and their resolution, so the
reasoning survives past the chat session that found it. Add new entries
above old ones within a section; don't delete resolved entries.

---

## Resolved — 2026-08-18 (LuGre sub-branch)

Full detail lives in `lugre_friction/README.md`; this is a pointer entry so
the backlog stays the single place that lists what's been found.

- **Sign convention on `v_nut`.** The LuGre `v_nut = (L/2*pi)*theta_s_dot -
  (x_n_dot - x_s_dot)` given for this port is the *negative* of the
  baseline's "closing velocity" convention (issue #1 above). Re-derived the
  `F_fric,nut` injection signs by virtual work rather than reusing the
  baseline's `+F_nut`/`-F_nut` labels: `x_s` row flips to `-F_fric,nut`,
  `x_n` row flips to `+F_fric,nut`, `theta_s` row unchanged. Verified by a
  step-input smoke test: positive command gives positive `theta_s` and
  `x_n`, slightly negative `x_s` (screw compressed into the thrust
  bearing), matching the physical description.
- **Phase-measurement window bug.** The zero-crossing-based phase
  calculation started its "steady state" window exactly at
  `N_TRANSIENT_PERIODS*T`, which could clip the very first target crossing
  (no valid "before" sample), causing that one term to match against a
  crossing a full period away and silently biasing the averaged phase by
  about `-360/N_PERIODS` degrees -- surfaced as a suspicious, frequency-
  independent `-90 deg` at every low/mid frequency, contradicted by the raw
  trace (which showed near-zero lag). Fixed by padding the window start
  back by half a period and rejecting any matched pair with `|dt| > T/2`.
- **`sigma0_sb` three orders of magnitude too stiff.** Reused directly from
  Rev 3's `d_sigma0` = 3.0e6 N*m/rad; paired with `I_sb` = 1.5e-7 kg*m^2
  that's a bristle natural frequency of ~712 kHz, versus a ~6863 Hz max
  across every real structural mode in this system. Confirmed as the cause
  of Radau and BDF diverging to NaN (tolerance/max_step tuning didn't
  help). Revised to 500 N*m/rad (~9.2 kHz, comparable to `k_s1`/`k_s2`);
  LSODA's runtime on the same case dropped from 230s to 71s afterward.
- **Radau/BDF don't work on this model; LSODA does.** Even after the
  `sigma0_sb` fix, Radau and BDF still diverge to NaN on the same
  (A=5e-3 rad, omega=1 rad/s) case in under a second, unaffected by
  tolerance tuning -- consistent with LuGre's non-smooth kink at `v=0`
  defeating their Newton-iteration/Jacobian machinery. LSODA integrates it
  reliably. `run_sinusoidal_sweep.py` uses `method="LSODA"` by default;
  revisiting Radau/BDF would need an analytical Jacobian, not just tuning.

## Resolved — 2026-08-18

### 1. Screw-nut coupling had the wrong sign (K, C — Sec. 5.3/5.4)

**Issue.** After recoupling the screw-nut interface directly into `K`/`C`
(replacing the frozen `T_fric,nut` friction port with a linear reflected
coupling between `theta_s` and the axial coordinates, scaled by the lead
ratio `L/2*pi`), the command-to-stage transfer function `x_n(s)/theta_cmd(s)`
had a **negative DC gain**: a unit `theta_cmd` produced `x_n = -1.326e-4`,
while the ideal no-slip kinematic relation from Sec. 1 (`x_n = (L/2*pi)*theta_s`)
predicts `+1.326e-4` — same magnitude, opposite sign. This showed up as the
Bode phase plot starting at 180 deg instead of 0 deg.

**Cause.** The four off-diagonal `(L/2*pi)*k_nut` (and `c_nut`) cross-terms
between `theta_s` and `x_s`/`x_n` had the sign that corresponds to the
internal coupling variable `x_n - x_s + (L/2*pi)*theta_s`. The physically
correct variable, matching Sec. 1's "positive rotation drives positive
translation" convention, is `x_n - x_s - (L/2*pi)*theta_s` (the nut's actual
position minus its ideal no-slip position `x_s + (L/2*pi)*theta_s`) —
opposite sign on the `theta_s` term.

**Fix.** Negated the four cross-terms in `K` and `C`
(`K[theta_s,x_s]=K[x_s,theta_s]`, `K[theta_s,x_n]=K[x_n,theta_s]`, and the
same positions in `C`), and the matching sign inside `F_nut` in Sec. 3.4.
Verified by direct solve of `K*q = B_u*e_1`: `x_n` flips from `-1.326e-4` to
`+1.326e-4`, matching the ideal ratio exactly. No other matrix entries
changed, and this sign flip does not alter the eigenvalues of `A` (it is
equivalent to relabeling `x_s -> -x_s`, `x_n -> -x_n`, a diagonal similarity
transform), so the modal frequencies and damping ratios found in issue #2
below are unaffected by it.

**Where:** `state_space_6dof.md` Sec. 3.4, 4.2, 5.3, 5.4; `build_bode_rev4.py`
`build_matrices()`.

---

### 2. Mode 1 (176.7 Hz) was almost undamped, producing a runaway resonance peak and a ~74 s settling time

**Issue.** Two symptoms traced back to the same root cause:
- The Bode magnitude near 176.7 Hz behaved like a near-singular spike rather
  than a finite resonance: an initial coarse sweep (0.1 Hz steps, 0-2000 Hz)
  happened to sample `|G| = 0.777` (-2.19 dB) there; a later fine sweep
  (0.05 Hz steps around the peak) found the *true* peak at `|G| = 1.446`
  (+3.21 dB) — both numbers are grid-dependent artifacts of sampling a very
  narrow, very tall peak, not a stable finite-Q resonance.
- The mode's 2%-settling time was **~74 seconds**, which makes a
  "let the system settle between steps" time-domain run impractical (the
  stepping-trajectory task needed a settled dwell on the order of a few
  hundred ms, not tens of seconds).

**Cause.** With only the original placeholder dampers (`c_c`, `c_s1`,
`c_s2`, `c_nut`, `c_brg` — each sized to ~2% of critical damping for its own
*local* two-body reduced-mass pair, in isolation), mode 1 came out at
`zeta_1 = 4.89e-5`. Root cause: mode 1's eigenvector has all four rotational
DOFs (`theta_m`, `theta_c`, `theta_s`, `theta_sb`) moving in near-unison —
nearly equal amplitude, same phase — i.e. the whole rotational chain rocking
almost rigidly against the electromagnetic spring `k_EM`. The dampers
`c_c`, `c_s1`, `c_s2` only act on *relative* velocities between adjacent
rotational DOFs, and those relative velocities are all close to zero in this
mode shape, so those dampers barely dissipate any energy from it — no
placeholder value assigned to them meaningfully changes `zeta_1`.
(Confirmed by scanning `c_c` up to 450x its placeholder value: `zeta_1` only
reached ~0.005.)

**Fix.** Added `c_EM` (Sec. 3.1/5.4) to `C[0,0]` — plain viscous damping of
`theta_m` against a fixed frame, which *does* act directly on mode 1's
dominant motion. This is a deliberate simplification of the literal Sec. 3.1
relation (`T_EM` damping term depends on `theta_dot_cmd - theta_dot_m`, not
just `-theta_dot_m`): the paired `theta_dot_cmd` feedforward column in `B_u`
is intentionally **not** implemented, because `theta_cmd` is a step
staircase in the time-domain sims this model feeds, and its derivative is a
train of Dirac impulses — not representable on a finite time grid.

The value, `c_EM = 1.3283e-4`, was found by bisection to put `zeta_1` at
exactly 0.02 (2%, per direct request), and barely disturbs the other five
modes (they were already close to their own ~2% placeholders):

| mode | frequency | zeta before | zeta after |
|---|---|---|---|
| 1 | 176.69 Hz | 4.89e-5 | **0.020000** |
| 2 | 745.97 Hz | 7.24e-3 | 7.27e-3 |
| 3 | 1650.23 Hz | 1.69e-2 | 2.18e-2 |
| 4 | 3429.32 Hz | 1.99e-2 | 1.99e-2 |
| 5 | 6536.22 Hz | 2.40e-2 | 2.40e-2 |
| 6 | 6863.23 Hz | 2.19e-2 | 2.19e-2 |

**Verification.** With `c_EM` applied, the true peak at 176.7 Hz dropped from
`|G| = 1.446` (+3.21 dB) to `|G| = 0.00354` (**-49.03 dB**) — a 52 dB drop,
matching the ~1/(2*zeta) scaling expected from a ~408x increase in `zeta_1`.
Mode 1's 2%-settling time dropped from ~74 s to ~176 ms (now the longest of
the six modes, but practical — this is what set `DWELL_SETTLE = 250 ms` in
`generate_stepping_trajectory.py`).

**Where:** `model_parameters.json` (`c_EM`, plus `parameter_notes.c_EM`);
`state_space_6dof.md` Sec. 5.4; `build_bode_rev4.py` `build_matrices()`.

---

## Resolved — 2026-08-21

### 3. Does the frozen LuGre linearization reference velocity (`V_STAGE`) affect the Bode/Co-MAC results?

**Question.** `run_local_linearization_bode.py` freezes the LuGre bristle at a
single operating point, `V_STAGE = 5 mm/s`, to get `K_eq`/`C_eq` for the
Co-MAC Set B matrices (`lugre_friction/rendered_assets/temp/Co-MAC/`). Raised
while reviewing the Co-MAC derivation: would a different reference velocity
choice change `K_eq`/`C_eq`, and therefore the Bode plot and the Co-MAC
numbers?

**Finding.** No. `sigma1_nut = sigma1_sb = sigma1_way = 0.0` for every port in
the current parameter set (confirmed directly from
`lugre_model.load_parameters()`), which makes the LuGre force law
`F(z, v) = sigma0*z + sigma2*v` — exactly linear in `(z, v)`, with no term
coupling the tangent stiffness/damping to the operating point. So
`K_eq = dF/dz = sigma0` and `C_eq = dF/dv = sigma2` are mathematically
invariant to `v0`, for every port, at every velocity including zero.

**Verification.** `plot_lugre_velocity_sensitivity.py` (writes
`rendered_assets/temp/lugre_velocity_sensitivity.png`):
- Swept `v0` directly through `equivalent_stiffness_damping()` from 1e-6 to
  1 m/s for all three ports: `K_eq/sigma0` stays at exactly 1.0 throughout
  (what actually moves with `v0` is `z0`, the frozen bristle deflection,
  following the Stribeck curve — not the tangent stiffness evaluated there).
- Rebuilt the full `K`, `C`, and Bode response from scratch at five
  different `V_STAGE` choices (0.05, 0.5, 5, 50, 500 mm/s): residual against
  the 5 mm/s reference tops out at 2.7e-5 dB across 0-2000 Hz — floating-
  point noise, not a real dependence.

**Implication for Co-MAC.** The `V_STAGE = 5 mm/s` choice used throughout
`Co-MAC/` is provably inconsequential — any velocity, including 0, gives
bit-identical `K1`/`C1`. The substantive assumption isn't "which velocity to
freeze at," it's `sigma1 = 0`: this removes any sliding-regime softening of
the *tangent* stiffness from the model entirely, so Set B's `K_eq` is the
same "stuck/presliding" tangent stiffness (`sigma0`) at every velocity,
including full sliding. Not an error in the Co-MAC derivation as done, but a
scope caveat worth flagging — if the intended physical picture requires
sliding friction to present a different dynamic stiffness than presliding
contact, that behavior isn't present in this parameterization at all.

**Where:** `lugre_friction/scripts/plot_lugre_velocity_sensitivity.py` (new);
`lugre_friction/scripts/run_local_linearization_bode.py` (already disclosed
this in its own module docstring, not previously verified numerically);
`lugre_friction/rendered_assets/temp/Co-MAC/comac_km_matrices.md`.

---

## Open / next

- **LuGre friction is implemented as a sub-branch, not yet merged.**
  `lugre_friction/` (own `model_parameters.json`, `lugre_model.py`,
  `run_sinusoidal_sweep.py`, `README.md`) implements all three ports
  (screw-nut, support bearing, guideway) as a 15-state nonlinear ODE and
  runs a sinusoidal describing-function sweep. It does not touch anything
  in `Rev 4/` itself — the frictionless recoupled baseline (this file's
  first two entries) stays the reference/comparison point. Folding the
  LuGre ports back into the main Rev 4 line (replacing the frozen-at-zero
  ports and the rigid `k_nut`/`c_nut` embedding) is still open.
- **Radau/BDF don't work on the LuGre model.** See the LuGre entry above
  and `lugre_friction/README.md` — only LSODA is verified reliable. Would
  need an analytical Jacobian (or event-based re-starts at `v=0`) to make
  the implicit fixed-Jacobian methods viable, which matters if the sweep
  needs to go faster than LSODA's ~5 minutes for 40 points.
- **Re-check `c_EM` once LuGre is the live friction model.** It was sized
  (issue #2 above) against zero friction damping at any port. If the
  support-bearing or nut LuGre ports end up dissipating real energy out of
  mode 1 once they're live in the main line (not just the sub-branch),
  `c_EM`'s value should be revisited rather than assumed to still be
  correct.
- **Revisit `sigma1 = 0` if sliding-regime stiffness softening matters.**
  Confirmed 2026-08-21 (Resolved section above) that the current
  parameterization gives velocity-invariant `K_eq`/`C_eq` at all three
  friction ports specifically because `sigma1 = 0` everywhere. If the thesis
  needs the Co-MAC/Bode comparison to reflect a genuinely different dynamic
  stiffness during sliding vs. presliding contact, `sigma1` would need a
  nonzero value at the relevant port(s) — the "frozen at `V_STAGE`"
  linearization step would then actually matter, unlike today.

# Parameter Optimization

Staging area for identifying the Rev 4.2 model's friction parameters
against the v3 hardware campaign. Nothing here runs an optimizer yet --
`model_parameters_optimization.json` is a review document: it partitions
every Rev 4.2 parameter into three groups so it's clear what's being held
fixed, what's derived per run from a physical relation, and what's
actually free for a fitting routine to move.

## Structure

- **`fixed_parameters`**: geometry, inertia, and structural stiffness/
  damping -- treated as known and not touched by the optimizer.
- **`current_dependent_parameters`**: currently just `T_hold`, the peak
  holding-torque amplitude. Rev 4.2's baseline model applies one constant
  value to every simulation; v3's six runs actually alternate between two
  drive/hold currents (200 mA and 400 mA peak), so `T_hold` is recomputed
  per run from `T_hold(I) = sqrt(2) * K_t * I` instead. See the main
  `README.md`, "Current-dependent holding torque", for the derivation and
  its effect on the model's state-space equations.
- **`variable_parameters`**: the three LuGre friction ports' parameters
  (`way` = linear guideway, `nut` = leadscrew-nut parasitic pre-rolling
  drag, `sb` = support bearing), each with its current Rev 4.2 value and a
  proposed `bounds: [min, max]` for an optimizer to search within.

## Before running anything

The bounds are a first-pass proposal (roughly 0.1x-10x the current value,
or a physically-motivated velocity/force range) -- not derived from a
sensitivity study. Go through `variable_parameters` and adjust:

- Any bound that's too tight (would clip the optimizer against a wall) or
  too loose (would let it wander somewhere non-physical).
- Whether `sigma1_*` should stay independently free, or be re-coupled to
  `sigma0_*` via the damping-ratio relation the parent model used to set
  it (`sigma1 = 2*zeta*sqrt(sigma0*m_eff)`, noted per-port in the JSON).
- Whether `Fc <= Fs` (`Tc_sb <= Ts_sb`) should be enforced as an explicit
  optimizer constraint rather than left to the bounds alone.
- Whether `T_d` (detent torque amplitude, currently fixed) should move to
  `variable_parameters` if the detent-torque montage panels don't fit
  well once real optimization starts.

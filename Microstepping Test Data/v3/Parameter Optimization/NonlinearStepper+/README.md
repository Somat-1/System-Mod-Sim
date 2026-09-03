# Nonlinear stepper drive law (Option A)

## Reference and why it matters here

MathWorks Simscape `Stepper Motor` block electrical torque
([https://www.mathworks.com/help/sps/ref/steppermotor.html](https://www.mathworks.com/help/sps/ref/steppermotor.html)):

```
Te = -Km*(iA - eA/Rm)*sin(Nr*theta) + Km*(iB - eB/Rm)*cos(Nr*theta) - Td*sin(4*Nr*theta)
```

Two things fall out of comparing this to `Rev 4/lugre_friction/Rev 4.2/scripts/lugre_model_rev42.py`:

1. **The detent term is independently validated.** `-Td*sin(4*Nr*theta)` is
   character-for-character what this repo already implements
   (`lugre_model_rev42.py:216`, `detent[0] = self.p["T_d"] * np.sin(4.0 * self.p["N_r"] * q[0])`).
   Citable as a reference-implementation match.
2. **The electromagnetic part is a current projection, not a spring.**
   Under an ideal-chopper assumption (commanded current tracks perfectly,
   back-EMF neglected — no new states), the two-phase sum collapses to
   `T_hold * sin(N_r * theta_err)`, `theta_err = theta_cmd - theta_m`. This
   *saturates* at `+-T_hold` when `N_r*theta_err = +-pi/2` (one full step
   of lag) and then **rolls back down** past that point — the actual
   pull-out/stall mechanism. There is no ceiling *parameter* anywhere;
   the ceiling emerges from the sine itself.

That is exactly the structural gap identified earlier in this campaign:
the model currently in use (`build_structural_matrices`, `command =
[k_em, 0, ...]`, `k_em = N_r*T_hold`) implements only the **small-angle
linearization** of this law, `T_hold*sin(N_r*theta_err) -> k_em*theta_err`
as `theta_err -> 0`. That linear spring has no ceiling: however large the
tracking error grows, it keeps supplying more restoring torque, so the
model structurally cannot represent a stall or lost step — a candidate
explanation for the large, parameter-insensitive D_70/D_200 model-vs-
measurement mismatch found in the sensitivity/ablation work
(`../calibration_bracketing/`).

## Why not port to Simscape directly

- You'd be rebuilding a 15-state model (three custom LuGre friction
  ports, a 6-DOF drivetrain, an analytical Jacobian, a Radau harness) to
  gain one equation that's an afternoon of work to write directly.
- The Simscape block **cannot model a hybrid stepper** (torque from both
  permanent-magnet and variable-reluctance effects) — and at 1.8°/step
  with `N_r=50`, this is a hybrid motor. Same limitation either way.
- Its Averaged mode's slip detection is approximate (torque-vector
  threshold sustained over one step period) and MathWorks' own docs warn
  results may be wrong once slipping starts. Stepping mode is accurate,
  but only because of the same sin/cos law implemented directly below.

## Two implementation options

**Option A — current-mode, no new states** (implemented here).
Replace `k_em*(theta_cmd - theta_m)` with `T_hold*sin(N_r*(theta_cmd -
theta_m))`. Reduces to the existing model for small lag (`sin(x) -> x`
gives `N_r*T_hold = k_em`). Ceiling at `T_hold`, reached at one full step
(`pi/(2*N_r)` rad = **1.8°** for `N_r=50`). Assumes the chopper delivers
commanded current perfectly (no winding dynamics).

**Option B — add the two winding currents as states** (not implemented
here). `diA/dt = (vA - R*iA - eA)/L`, `eA = -Km*omega*sin(N_r*theta)` —
17 states total. Captures torque roll-off with speed (electrical time
constant limiting current rise at high step rates), which Option A
cannot.

**Which is needed depends on where step loss lives, and the data can
discriminate**: if closure errors appear only at the fast blocks (D_70,
D_200), it's speed-dependent/electrical -> Option B. If they appear
across rates including the slow blocks, it's the static pull-out limit
-> Option A is enough. Start with A regardless: one-line change,
parameterized entirely from what's already in `model_parameters.json`,
turns the model from something that structurally cannot lose steps into
something that can. Add electrical states only if Option A still
under-predicts loss specifically at the fast blocks.

## What's implemented in this folder

- **`lugre_model_rev42_optionA.py`** — `LuGreModelRev42OptionA`, a variant
  of `LuGreModelRev42` with *only* the drive term changed:
  `build_structural_matrices_optionA` drops `k_em` from the `(0,0)`
  stiffness entry and drops the command vector entirely; `rhs()` computes
  `T_hold*sin(N_r*theta_err)` directly and applies it to the `theta_m` row;
  `analytical_linearization(state, theta_cmd)` uses the *local* tangent
  `N_r*T_hold*cos(N_r*theta_err)` in place of the constant `k_em` — this
  is why, unlike the linear model, it needs `theta_cmd` as an explicit
  argument (the local slope is operating-point-dependent, exactly like
  the detent's own linearization already is). Everything else (mass,
  damping, `k_c/k_s1/k_s2/k_nut/k_brg` chain, three LuGre ports, exact
  `sin()` detent) is untouched.
- **`step_compare_linear_vs_optionA.py`** — runs both models on run 2,
  blocks `D_0.125`, `D_3.5`, `D_70`, `D_200` (the same four used in the
  earlier detent-vs-friction ablation), and produces:
  - `torque_law_linear_vs_optionA.png` — static analytical comparison of
    the two torque laws vs. tracking error, no simulation needed, with
    the `+-T_hold` ceiling and `+-1.8°` pull-out lines marked.
  - `compare_D_<rate>.png` (for the three controller-paced blocks) — a
    10-detent-cycle window centered mid-cruise (same windowing as
    `../torque_diagnostics/plot_cruise_zoom.py`), tracking error and
    motor torque, both models overlaid.
  - `compare_linear_vs_optionA.json` — per-block, per-model
    `cost_um` (RMS vs. measurement, same metric as `step6`/`step7`) and
    `max_abs_err_deg` (peak tracking error reached over the whole block,
    directly comparable to the 1.8° pull-out threshold — the
    discriminating diagnostic for whether Option A actually predicts
    anything close to step loss at D_70/D_200).

## Open decision: Km / T_hold current-scaling convention

The reference law's constant is `Km`, the true torque-per-amp constant,
with `T_hold = Km*I`. This repo's existing convention is instead
`T_hold(I) = sqrt(2)*K_t*I` with `K_t = 0.10606601717798213` N·m/A
(`model_parameters.json`), which implies an effective `Km = sqrt(2)*K_t
~= 0.150` N·m/A — matching **neither** real datasheet variant cited
(0674A: 0.090 N·m/A, 0956A: 0.045 N·m/A; ~1.7x and ~3.3x off,
respectively). Resolving which variant this motor actually is would let
`T_hold = Km*I` be set directly from a datasheet constant instead of the
current scaling assumption, and would also resolve the still-open
motor-variant/anchor question from Step 1 of the calibration-bracketing
work.

**Deliberately deferred here**: this comparison keeps the existing
`T_hold(I)` convention unchanged on *both* sides (linear and Option A),
so that the only variable under test is the drive law's functional form
(linear spring vs. sine), not the current-scaling constant. Swapping to
a datasheet-accurate `Km` is a follow-up, and a one-parameter change once
the motor variant is confirmed.

## Results

Ran on run 2 (MRES 1/4, 100% I), `D_0.125`/`D_3.5`/`D_70`/`D_200`:

| Block | linear cost (µm) | optA cost (µm) | linear max\|err\| | optA max\|err\| | pull-out |
|---|---:|---:|---:|---:|---:|
| D_0.125 | 0.7744 | 0.7746 | 0.5407° | 0.5408° | 1.80° |
| D_3.5   | 21.2709 | 21.2711 | 0.1258° | 0.1260° | 1.80° |
| D_70    | 148.9909 | 148.9912 | 0.1835° | 0.1845° | 1.80° |
| D_200   | 891.2447 | 891.2451 | 0.3138° | 0.3188° | 1.80° |

**Regression check passed**: linear-model costs match the previously
computed step6/step7 values exactly (e.g. D_3.5 = 21.2709 µm), confirming
this driver introduces no discrepancy relative to the existing pipeline.

**Small-lag reduction confirmed**: Option A tracks the linear model to
within 0.0002–0.005 µm cost and <0.005° tracking error at every block —
exactly the `sin(x) -> x` behavior the derivation predicts.

**Discriminating-test verdict: Option A does NOT explain the D_70/D_200
anomaly.** Peak tracking error never exceeds 0.32° at any tested block —
roughly 1/6 of the 1.80° pull-out threshold, even at D_200 where the
model-vs-measurement mismatch is 891 µm. Cost differs between the two
drive laws by <0.001% everywhere; the two laws never meaningfully diverge
because the error never leaves the near-linear region of the sine. Under
Option A's ideal-chopper assumption (commanded current delivered
perfectly, no winding dynamics), this motor is nowhere close to stalling
at any rate tested here — the 891 µm/149 µm anomaly has some other cause.

Tracking error is also **not monotonic in step rate**: largest at the
slowest block (D_0.125, 0.54°) and smallest at D_3.5/D_70 (~0.13–0.18°),
only ticking back up slightly at D_200 (0.32°). D_0.125's error looks like
a quasi-static deflection against load/friction/detent through the stiff
drivetrain rather than a dynamic tracking lag, so "faster commanded rate"
does not simply mean "closer to pull-out" in this model.

**This is the "Option A still under-predicts loss at D_200 specifically"
signal** the decision rule above was watching for, and it lines up with
the electrical-time-constant argument raised earlier: at D_200 the
0674A-datasheet time constant (`tau = L/R ~= 0.61 ms`) is about half the
1.25 ms microstep period, so real winding current plausibly cannot fully
settle between microsteps — exactly the effect Option A's ideal-chopper
assumption excludes by construction, and exactly what Option B (winding
currents as states) would capture. **Next step, if pursued: implement
Option B** rather than concluding the stepper-drive nonlinearity is a
dead end — Option A ruled out the *static* pull-out limit at these rates,
it did not rule out a *speed-dependent electrical* torque shortfall.

Plots: `torque_law_linear_vs_optionA.png` (static law comparison),
`compare_D_3.5.png` / `compare_D_70.png` / `compare_D_200.png` (per-block
cruise-window tracking-error and motor-torque overlays — visually
confirm the two traces are near-indistinguishable at this scale, matching
the table above). Raw numbers in `compare_linear_vs_optionA.json`.

### Why the comparison came out flat: the command itself is smoothed

`step_plot_commanded_sequence.py` -> `commanded_sequence_smooth_vs_discrete.png`
plots the smooth `trapezoid_fraction()` reconstruction that every
simulation in this repo actually integrates against the **true discrete
microstep staircase** a real stepper receives (one step every
`1/(rate*MRES)` s, `MRES=4` for run 2). The gap between the two is a
constant **one microstep (0.45°)** at every rate — invisible against
D_3.5's 70° window, clearly visible at D_70, and the dominant feature of
the D_200 window (a ~0.45° sawtooth riding on the ramp, comparable in
size to the window's own ~12.5° span).

That constant 0.45° gap is the **same order of magnitude** as the
tracking errors found above (0.13°–0.54°). Since neither the linear nor
Option A drive law ever sees a step edge — both only ever integrate the
smooth blue line — this comparison was structurally incapable of
surfacing a discretization-driven tracking-error effect even if one
exists on the real hardware. It isn't just that Option A found nothing;
the input it was given couldn't have shown it something.

This sharpens the Option B case further: feeding a true discrete
staircase into the linear/Option A drive (both assume commanded current
is delivered instantly) would just produce artificial torque spikes at
every microstep edge — not physically real either. The combination that
would actually be faithful to the hardware is **discrete step commands +
Option B's winding-current dynamics**, since the coil inductance is
exactly what smooths a raw step edge into a continuous current rise in
reality — the same smoothing this repo's reconstruction currently does
by assumption rather than by physics. Implementing Option B without also
switching to a discrete command would only add electrical dynamics on
top of a command that's already artificially pre-smoothed; the two
changes belong together.

### Net drive torque (motor minus detent)

`step_plot_net_drive_minus_detent.py` -> `net_drive_minus_detent.png`
plots `net = T_motor - T_detent` (Option A; linear is indistinguishable),
the grouped term that's actually left over to drive the rest of the
drivetrain once the motor and detent have "fought it out" at the theta_m
node (`I_m*theta_m_ddot = (T_motor - T_detent) - k_c*(theta_m-theta_c) -
...`). **Sampling note**: the windowed segment's `t_eval` density had to
be bumped from a flat 200 samples/s to `max(200, rate_hz*40)` in
`simulate_variant` -- the flat rate was only ~2.9 samples/detent-cycle at
D_70 and exactly at the Nyquist limit at D_200, producing a visibly
jagged, artifact-looking waveform in the first render. The fix only
affects this segment's *output interpolation* density (cheap -- the
adaptive Radau integration itself, and therefore every `cost_um`/
`max_abs_err_deg` value reported above, is unaffected).

With correct sampling: `T_motor` and `T_detent` are in phase (same sign,
comparable shape) at all three rates tested — consistent with, not
contradicting, the motor/detent rows in
`../torque_diagnostics/figureB_cruise_zoom_D9.5_D200.png`. Because they
are in phase, subtracting one from the other is what produces
cancellation, not opposition: at D_3.5 their amplitudes nearly match and
`net` stays almost flat; at D_70 the match is partial and `net` shows a
real, visible ripple well below either raw component; at D_200 the
amplitude match breaks down and `net`'s swing becomes comparable to
`T_motor` itself — the near-total cancellation seen at low rate has
essentially disappeared. This tracks the detent forcing frequency (200 Hz
at D_200) approaching the ~183-211 Hz drivetrain resonance identified in
`../calibration_bracketing/step4c_ringdown_fit.png`, consistent with a
control-loop-style cancellation degrading near a structural mode.

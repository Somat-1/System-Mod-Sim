# Calibration bracketing — results

Cross-run investigation (run 2 unless noted; runs 5/6 for the MRES 1/1
tread analysis) run against the Rev 4.2 parallel-LuGre model, following up
on the torque-diagnostics work in `../run_02_mres_4_i_100pct/torque_diagnostics/`.
Scripts and plots referenced below all live in this same folder.

## Step 1 — T_hold anchor and rotor tangent stiffness

**Script/plot**: `step1_rotor_tangent_stiffness.py` / `step1_rotor_tangent_stiffness.png`

Confirmed `k_em = N_r * T_hold`:

| Current | T_hold | k_em |
|---|---:|---:|
| 200 mA (I_50pct) | 30.0 mN·m | 1.500 N·m/rad |
| 400 mA (I_100pct) | 60.0 mN·m | 3.000 N·m/rad |

Plotted the combined rotor tangent stiffness across one full mechanical step:

```
k_eff(theta_m) = k_em + 4*N_r*T_d*cos(4*N_r*theta_m),   k_d = 4*N_r*T_d = 1.0 N*m/rad
```

| Current | k_eff range | Stable? |
|---|---|---|
| 200 mA | [0.500, 2.500] N·m/rad | Yes (k_eff > 0 everywhere) |
| 400 mA | [2.000, 4.000] N·m/rad | Yes (k_eff > 0 everywhere) |

**Reading**: both currents stay stable (no local instability/cogging-jump
risk anywhere in a full step), but the 50%-current margin (min = 0.5) is
markedly thinner than the 400 mA case, and thinner than it was before the
T_d correction (was ~0.9 at T_d = 3.0 mN·m). Worth re-checking if T_d ever
moves further.

## Step 2 — epsilon sensitivity on the C-block dwell

**Script/plot**: `step2_epsilon_dwell_check.py` / `step2_epsilon_dwell_comparison.png`

Ran the full C block (run 2) once at `smooth_velocity_epsilon = 1e-9` and
once at `1e-12`, compared simulated position throughout.

- Max |difference| over the whole block: **0.0000 um (~0.002-0.004 nm at
  the two move transients; indistinguishable from zero during both 60 s
  dwells)**.
- Difference at block end: 0.00000 um.

**Reading**: epsilon is inert at either tested value, for this solver
setup. Most likely explanation: the ODE solver's own `ATOL` on the
velocity states (1e-7 to 1e-9 depending on state) already dominates over
epsilon's regularization floor at both 1e-9 and 1e-12, so the two choices
are numerically indistinguishable here. Conclusion: stop worrying about
epsilon for this model/solver combination -- use 1e-9 everywhere, no
special-casing needed for C-block fits.

## Step 3 — nut damping split (declared)

**Files changed**: `Rev 4/lugre_friction/Rev 4.2/model_parameters.json`,
`Parameter Optimization/model_parameters_optimization.json`

`c_nut`, `sigma1_nut`, and `sigma2_nut` all act on the same relative
velocity `v_nut` and are linearly indistinguishable at a single operating
point (identifiability problem, not a physics error). **Decision: `c_nut`
fixed to 0** (was 101 N·s/m); LuGre's `sigma1_nut`/`sigma2_nut` now carry
all nut-port dissipation. No other value needed to change --
`sigma1_nut = 500.29 N·s/m` was already set independently via the
target-zeta=0.7 rule (`sigma1 = 2*zeta*sqrt(sigma0*m_eff)`), with no
reference to `c_nut`. `k_nut` (the load-bearing spring) is untouched.

## Step 4a — sigma0_way / Fs_way bracket from D_0.125 (MRES 1/1)

**Script/plots**: `step4_bracket_sigma0_way.py` / `step4a_d0125_treads.png`

At MRES 1/1, each `D_0.125` MOVE_ACK is exactly one full 10 um step, and
pulses are ~8 s apart -- each step's transient is fully isolated.

| Run | t (s) | Step height (um) | Pre-jump creep (um) | Creep fraction |
|---|---:|---:|---:|---:|
| 5 | 8.22 | +10.232 | +2.251 | 0.220 |
| 5 | 16.23 | +9.295 | +0.857 | 0.092 |
| 5 | 28.23 | -8.914 | -2.760 | 0.310 |
| 5 | 36.22 | -10.444 | -1.548 | 0.148 |
| 6 | 8.22 | +10.110 | +0.567 | 0.056 |
| 6 | 16.23 | +9.570 | +0.777 | 0.081 |
| 6 | 28.23 | -9.392 | -1.058 | 0.113 |
| 6 | 36.22 | -10.179 | -1.646 | 0.162 |

- Median creep fraction: **0.130** (13%).
- Rev-3-placeholder prediction (`Fs_way/sigma0_way = 4.5/7.6e5 = 5.9 um`
  out of a 10 um step) would give creep_frac ~ 0.59.
- **Implied sigma0_way (holding Fs_way = 4.5 N fixed): ~3.45e6 N/m**,
  roughly 4.5x the Rev-3 carry-over value.

**Reading**: the measured treads are genuinely close to square (visible
directly in `step4a_d0125_treads.png` -- a sharp jump with a small
overshoot/ring-down, not a long compliant ramp). This confirms the
hypothesis that `sigma0_way = 7.6e5 N/m` is too soft; **7.6e5 N/m should
not be treated as a credible fixed point going into the optimizer** --
bracket around ~3.45e6 N/m (order-of-magnitude estimate, not a precise
fit) pending a joint sigma0/Fs fit.

*(Bug fixed along the way: the first version of this script indexed the
IDS position array using block-relative time directly instead of adding
the block's absolute start time back in, producing nonsense ~50 um "steps"
at unrelated points in the recording. Fixed by converting every window to
an absolute-time index before slicing `position_nm`.)*

## Step 4b — C-block approach presliding (inconclusive)

**Plot**: `step4b_c_approach_presliding.png`

Intended as a second, independent bracket on presliding distance, but the
approach move is a fast 150-full-steps/s, 4-step burst (~27 ms total), and
the plot shows most of the -40 um travel already complete *before* the
logged `MOVE_ACK` timestamp -- the same command-log-timestamp-lags-
completion effect identified earlier this session for fast bursts (the
ack is logged after a full command round-trip, not at physical motion
onset). This means the visible "ramp" before t=0 in the plot is largely
the ordinary multi-pulse execution profile, not isolated presliding
compliance. **No usable sigma0_way bracket from this measurement as
currently analyzed** -- would need the actual per-pulse commanded
staircase reconstructed and subtracted before any "excess compliance"
could be read off cleanly. Left as an open item, not force-fit.

## Step 5 — quasi-static detent transmission at 1.25 Hz

**Script/plot**: `step5_quasi_static_1p25hz.py` / `step5_quasi_static_1p25hz.png`

Compared measured vs. a **freshly simulated** (current T_d=5 mN·m,
current-dependent T_hold, c_nut=0 -- not the stale pipeline npz) D_1.25
cruise-phase ripple, isolating the 1.25 Hz component via a least-squares
sin/cos fit over the middle 80% of the positive-direction plateau
(27.9 s window, 34.8 cycles, both linearly detrended first):

| | Amplitude at 1.25 Hz |
|---|---:|
| Measured | 0.0629 um |
| Simulated | 0.0345 um |
| Ratio (sim/measured) | **0.548** |

**Reading**: the model under-predicts the quasi-static detent-to-stage
transmission by roughly **1.8x** at a frequency low enough (well below
the ~140 Hz drivetrain resonance found earlier) that damping and friction
parameters shouldn't be the main driver of this gap -- this is a
transmission-gain check, largely independent of the LuGre parameters
being bracketed elsewhere in this document. Two candidate explanations,
not distinguished by this check alone: T_d may still be too small, or the
static stiffness chain (`k_c`/`k_s1`/`k_s2`/`k_nut`) is transmitting less
rotor-side ripple to the stage than the real hardware does.

*(Bug fixed along the way: D_1.25 is software-paced -- every segment is a
`hold`, there is no `trapezoid` segment -- so the cruise-window logic
written for the fast/trapezoid rates didn't apply. Rewrote it to find the
positive-direction plateau span directly from the pulse train's jump
timestamps instead.)*

## Summary table

| Item | Status | Headline number |
|---|---|---|
| T_hold anchor | Confirmed | k_em = 1.5 / 3.0 N·m/rad (200/400 mA) |
| Rotor tangent stiffness | Checked | Stable at both currents; thin margin at 200 mA |
| Epsilon sensitivity | Resolved | Inert (differences ~picometers) -- use 1e-9 everywhere |
| Nut damping split | Decided | c_nut = 0; LuGre carries all nut dissipation |
| sigma0_way bracket (D_0.125) | Bracketed | ~3.45e6 N/m (vs. 7.6e5 placeholder) |
| sigma0_way bracket (C approach) | Inconclusive | Confounded by command-timestamp lag |
| Quasi-static detent transmission | Checked | Model at 0.55x of measured ripple at 1.25 Hz |

## Not yet done

**Step 6** (sensitivity table + staged fit within the Step 4 brackets) has
not been started -- scope to be agreed before running anything open-ended.

**Pipeline currency**: the T_d (3.0->5.0 mN·m) and c_nut (101->0) changes
have not yet been propagated into a full 66-block `simulate_block_responses.py`
rerun. `block_montage_simulated.png`, `torque_montage.png`, and
`tracking_error_montage.png` in every run folder are stale relative to the
current best parameters until that rerun happens.

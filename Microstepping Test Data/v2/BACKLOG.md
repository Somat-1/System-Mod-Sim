# Motion-sequence backlog

## Objective

Create a deterministic stepper-motor motion sequence for execution on the
physical stage. Record displacement with the interferometric sensor and
compare the measured response against the frictionless tangent, exact-detent
frictionless, and Rev 4.2 LuGre models.

No file in v2 currently issues hardware commands. Motion generation must
remain a preview/export step until the safety and pilot gates below pass.

## Revision layout

- data/: commanded event tables and imported measurement runs.
- scripts/: deterministic sequence builders, validators, and plotting code.
- rendered_assets/trajectory_visualization_plots/: command and measured-response
  trajectory figures, grouped by test purpose.
- rendered_assets/: run summaries, numerical exports, and comparison plots.
- docs/: derivations, controller protocol, and test notes.

## Run matrix

Keep MRES fixed during one execution and change it only on a full-step
boundary.

| MRES | u per driver pulse | Burst pulse rate |
|---|---:|---:|
| 1/16 | 1 | 4000 pulses/s |
| 1/4 | 4 | 1000 pulses/s |
| 1/2 | 8 | 500 pulses/s |
| 1/1 | 16 | 250 pulses/s |

All burst and doublet moves therefore use 250 full-steps/s. The outer loop is
four MRES settings by three holding-current settings, for 12 executions.
Every execution starts from the same fixed physical stage position; no
position sweep is included.
The numerical values of I_lo, I_mid, and I_hi remain to be supplied.

## Authoritative block definitions

### Block 0: reference fingerprint

Two seconds lead-in, then the following 12 moves with 1.0 s dwell after each,
then two seconds tail:

1. +16
2. -16
3. +4
4. -4
5. +1
6. -1
7. -16
8. +16
9. -4
10. +4
11. -1
12. +1

The block is net zero and lasts approximately 16 s including finite burst
durations. Store the residual after moves 2, 4, 6, 8, 10, and 12. Compare the
six start-of-run residuals against the end-of-run fingerprint for thermal
contamination.

### Conditioning

Move +4 full steps and then -4 full steps at 150 full-steps/s, followed by a
2 s settle. Conditioning is specified in physical full steps, not driver
pulses, so it is identical for all MRES settings.

### Block A1: unidirectional step-and-settle

For N in {1, 2, 4, 8, 16, 32} driver pulses: repeat +N and dwell, then repeat
-N and dwell. Execution uses 25 repetitions and the command preview shows 3,
annotated as execution x25. Ladder dwell is parameterized as dwell_ladder;
the preview and initial execution configuration use 0.4 s. Change it to
0.8 s everywhere if the pilot identifies zeta near 0.005.

### Block A2: alternating steps

For the same N set: repeat +N, dwell, -N, dwell. Execution uses 25
repetitions and the preview shows 3.

### Block B: nested reversal loops

Use 0.30 s dwell after each move and 10 repetitions of each pattern. The
preview shows 2 repetitions and is annotated as execution x10:

- Descending: +32, -16, +8, -4, +2, -1, +1, -2, +4, -8, +16, -32
- Asymmetric: +8, -3, +2, -5, +6, -8, +4, -4
- Minor loop: +64, -16, +2, -2, -16, +2, -2, -16, +2, -2, -16

All values in Block B are driver pulses at the loaded MRES.

### Block C: dwell and creep

- Arrive from positive: -4 full steps at 150 full-steps/s, dwell 1.0 s,
  +4 full steps at 150 full-steps/s, then immediately record a 60 s hold.
- Arrive from negative: +4 full steps at 150 full-steps/s, dwell 1.0 s,
  -4 full steps at 150 full-steps/s, then immediately record a 60 s hold.

Both approaches are net zero and specified in physical full steps.

### Block D: velocity plateaus

Use the following MRES-invariant full-step rates:

0.125, 0.375, 1.25, 3.5, 9.5, 27.5, 70, and 200 full-steps/s.

Run both directions. Plateau duration is
clip(200 full steps / f_fs, 5 s, 20 s), excluding ramps. This is one motor
revolution when the cap does not bind. Use a 0.5 s trapezoidal ramp in and
out only when f_fs is at least 10 full-steps/s; below 10 full-steps/s use no
ramp. Discard ramp intervals during identification. Pulse rate is the
full-step rate multiplied by the loaded microsteps per full step.

The slow plateaus below 10 full-steps/s are local friction measurements.
Execute them from the same fixed stage starting position as every other run.
They are therefore intentionally local measurements: position dependence is
not averaged and must be retained as a limitation when interpreting them.

### Block E: doublets

For N in {1, 2, 4, 8, 16} driver pulses: repeat +N at the MRES-specific burst
rate, immediately apply -N, then dwell 1.0 s. Execution uses 20 repetitions
and the preview shows 3, annotated as execution x20.

Run Block 0 again after Block E.

## Figure plan

Generate one figure per block:

- Columns: MRES = 1/16, 1/4, 1/2, 1/1.
- Rows: block sub-condition, such as N, loop pattern, or full-step rate.
- First figures are command-only previews in physical travel units. Current
  does not change the command trajectory, so current traces are not shown.
- Measured/model-response figures later overlay I_lo, I_mid, and I_hi in each
  cell while retaining the same row/column skeleton.
- Block 0 measurement result: one residual-versus-doublet-index acceptance
  chart containing all 12 executions.

Store trajectory figures under `rendered_assets/trajectory_visualization_plots`
in four groups: reference and conditioning, step ladders, reversal and creep,
and velocity and doublets.

## Pilot and safety gates

- Determine the final A1/A2 dwell from fitted damping ratio:
  0.4 s when zeta is at least 0.01; 0.8 s near zeta = 0.005.
- Determine the A1/A2 repeat count from the 200-repeat N=1 breakaway pilot.
  Use 25 only if coefficient of variation is at most 25 percent.
- Step-loss pilot: 100 repetitions of +32, -32 at MRES 1/16 and again at
  MRES 1/1. Reject monotonic interferometer endpoint drift.
- Verify conditioning produces no sustained low-speed resonance.
- Supply I_lo, I_mid, and I_hi.
- Confirm stage travel, direction convention, limit-switch behavior, pulse
  timing limits, enable sequencing, and emergency-stop procedure.
- Log actual pulse timestamps, MRES, current, enable state, and interferometer
  timestamps for every execution.

## Remaining inputs before hardware execution

1. Supply I_lo, I_mid, and I_hi.
2. Replace dwell_ladder if required by the damping-ratio pilot.
3. Revisit the execution repeat count if the 200-repeat breakaway pilot gives
   a coefficient of variation above 25 percent.

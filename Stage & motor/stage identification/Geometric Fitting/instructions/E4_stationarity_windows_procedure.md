# E4 — stationarity windows

## Purpose

Model C assumes that the screw-periodic error has the same harmonic amplitude
and phase over the full screw. E4 tests that assumption directly: it estimates
the periodic error locally at seven positions and tests whether the fitted
harmonics change with position.

Runout and eccentricity may vary along the screw, for example near mid-span or
the nut preload region. A stationary model would then over-compensate some
positions and under-compensate others. E1 can indicate that variation exists,
while the full-travel E2/E3 fits average it. E4 is designed to localise it.

## Experimental design

- Use seven windows centred at 10, 27, 44, 61, 78, 95, and 112 mm.
- Each window is 6 mm long: 7–13, 24–30, 41–47, 58–64, 75–81,
  92–98, and 109–115 mm.
- Every window start is an integer multiple of 64 microsteps.
- Sample at 0.12 mm = 192 microsteps, giving 50 intervals and 51 targets per
  window.
- Because one screw revolution is 3200 microsteps, the 192-microstep spacing
  advances three of 50 electrical screw-phase bins. The first 50 targets visit
  every phase bin exactly once; target 51 repeats the starting phase and closes
  the three-revolution interval.
- Traverse every window upward and downward using the E1 stop–settle–record
  routine: move, settle for 2.0 s, and record for 0.5 s.
- Use a 3 mm lower pre-approach and upper overrun. All E4 turn points lie inside
  the proven 2.000–120.695 mm soft envelope.
- Use SpreadCycle at 400 mA run and hold current, 2 mm/s, and 40 mm/s². Generate
  STEP with the buffered 1 MHz RMT path used by E1.

The core data set contains `7 × 51 × 2 = 714` analysis stops. With the fixed
reference records described below, the executable plan contains 726 recorded
stops and has a nominal duration of 36.0 minutes.

## Executable sequence

The stage is measurement-only during E4. Before arming it must be homed,
warmed by the established end-to-end routine, physically parked at
`x_ref = 61.35 mm`, still energised, and already being recorded by the IDS
interferometer.

The deterministic sequence is:

1. Verify the coordinate arm, driver, end switches, configuration, and logging.
2. Record three fixed-side datum visits at x_ref, always approached from
   58.35 mm.
3. Execute W1 upward and downward.
4. Record one x_ref visit, then execute W2.
5. Repeat the reference/window pattern through W7. There is one reference after
   W1–W6, not after W7.
6. After the W4 reference, project the remaining duration. If completion would
   exceed 58 minutes, skip W5–W7, mark stationarity coverage incomplete, record
   the end references, and park.
7. Record three end references, request the environmental snapshot, and park at
   x_ref using the fixed-side approach.

Abort, switch, driver, RMT, incomplete-motion, or soft-limit faults invalidate
completion. Parking is attempted only while the coordinate remains trustworthy;
otherwise the motor is released and manual recovery is required.

## Logging and acquisition products

The ESP32 logs the configuration hash, all integer target coordinates, screw
phase bins, movement and recording boundaries, post-window MSCNT, reference
visits, watchdog projection, and final completeness counts. The host terminal
stores the serial stream without interpretation.

The IDS/analysis side must create the normal `stops.csv`, `events.csv`,
`env.csv`, and `summary.json` products. For each stop, retain at least the
window, direction, target, absolute screw phase, mean displacement, standard
deviation, slope, sample count, validity flags, and reference correction.
Do not repeat rejected stops inline, because this would destroy the matched
timing and phase structure.

## Analysis

For each window and direction, fit the requested harmonics on absolute screw
angle:

\[
e_w(\theta)=a_{0,w}+\sum_{h=1}^{H}
\left[c_{h,w}\cos(h\theta)+s_{h,w}\sin(h\theta)\right].
\]

Compute amplitude and phase as

\[
A_{h,w}=\sqrt{c_{h,w}^2+s_{h,w}^2},\qquad
\phi_{h,w}=\operatorname{atan2}(-s_{h,w},c_{h,w}),
\]

with uncertainty propagated from the coefficient covariance or a window-level
bootstrap. Compare direction-averaged amplitudes and unwrapped absolute phases
across W1–W7. A stationary Model C requires the between-window variation to be
consistent with uncertainty and reference stability. A systematic amplitude
profile supports extending C with position-dependent harmonic amplitudes.
Abrupt phase discontinuities require investigation for lost steps or coupling
slip before any physical stationarity conclusion.

Fit the structured GP Model D on `E4 ∪ E2 ∪ E3` as an independent check of the
same position dependence.

## Acceptance and outputs

Require all seven windows for the primary stationarity verdict, valid paired
directions at at least 95% of targets, at least 95% valid analysis stops, no
unsafe or time-limit abort, consistent reference visits, total recording time
at or below 60 minutes, and air-temperature change below 0.5 K.

Report:

- harmonic amplitude and phase with uncertainty for every window;
- amplitude-versus-position and unwrapped phase-versus-position profiles;
- the formal stationary/non-stationary verdict for each retained harmonic;
- direction dependence, reference drift, and phase-continuity diagnostics;
- the Model-D cross-check and any proposed position-dependent extension to C.


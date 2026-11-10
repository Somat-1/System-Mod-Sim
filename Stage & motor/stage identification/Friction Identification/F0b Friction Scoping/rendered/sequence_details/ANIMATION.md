# L: illustrative friction loops

`f0b_L_friction_loops.gif` shows three representative ladder strokes
(0.195, 1.992 and 20 µm), first increasing and then decreasing. Each stroke
executes two complete cycles. The shortened selection carries model state
between strokes; it is not a simulation of every intervening F0b stroke.

The command, normalized lag loop, and cycle-summary plots are synchronized.
Gold dashed curves are cycle 1; teal solid curves are cycle 2. During the
decreasing pass, grey curves show the corresponding increasing-pass response.
Filled summary markers represent the increasing pass; open markers represent
the decreasing pass. Each summary point averages the two half-cycle losses
in that cycle: `loss = s - abs(delta_u + delta_e)`.

This is an illustrative, rate-independent Maxwell-slip model, not measured
carriage motion or a force estimate. Six displacement thresholds
`[0.05, 0.15, 0.4, 1, 2.5, 6] µm` have equal stiffness weights, scaled so the
total lag saturates at ±1.757/2 µm. Thresholds are placeholders. The initial
state is positive saturation following a wipe from below. Small loops can
retain history-dependent vertical offsets; closed second cycles need not
coincide between passes. Small command motion is not necessarily entirely
absorbed by this model. Use F0b measurements to identify actual behavior.

The half-cycle axis and playback speed are schematic. Two-cycle repetition
shows the inherited initial branch and the subsequent closed loop; the
decreasing revisit exposes memory of larger excursions. There is a brief
pause after each cycle. The GIF repeats automatically.

Regenerate from the experiment directory with:

```sh
python3 render_f0b_l_animation.py
```

The script exports only the GIF, using a shared 64-colour palette and reducing
resolution as needed to keep the file below 5,000,000 bytes.

# S: illustrative slow ramps

`f0b_S_slow_ramps.gif` compares a 20 µm ramp after a reversal with a 20 µm
ramp in continuation over the same commanded positions. The reversal ramp is S1
(wipe from below, then ramp down). The continuation ramp is S4 (wipe from
above, then ramp down). Both use 512 microsteps at 1/256.

The animation has three synchronised panels:

- **(a) Mechanism.** A cartoon of the drive chain. The handle is the command
  u, the spring is the drive-chain compliance, and the block is the carriage
  x. The bristles are the friction elements: they are grey while sticking and
  ochre while slipping. Lag and travel share one scale, and the view pans to
  keep the carriage in frame. Bristle tilt and spacing are not to scale.
- **(b) Lag vs commanded travel.** The lag is e in the ramp direction,
  measured from the ramp start, plotted on a linear axis with a log-travel
  inset. S1 is solid, S4 dashed, and Δ = S1 − S4 bold.
- **(c) First 30 microsteps.** The commanded staircase and the carriage
  position on a time axis at 6.4 µsteps/s.

The model is the quasi-static Maxwell-slip model from the S-ramp instructions.
It is written in displacement units as `u − x = Σ γᵢ zᵢ`, with `|zᵢ| ≤ aᵢ`,
and solved exactly from one slip event to the next.

- **Slip limits:** the placeholders `aᵢ = [0.02, 0.08, 0.25, 0.8, 2.5] µm`.
- **Stiffnesses:** γᵢ, scaled so that Σγa = β₀/2.
- **Behaviour after a reversal:** the initial slope is −Σγ/(1+Σγ) = −0.92. The
  lag reaches exactly −β₀ (fully flipped) within the 20 µm.
- **Regimes in S1:**
  - stuck (slope < −0.8) up to 0.55 µm of commanded travel;
  - pre-sliding up to 3.3 µm;
  - sliding after that.

  These boundaries come from the placeholder slip limits, not from
  measurement.

The model is rate-independent and has no stick-slip, so the S4 lag is exactly
flat. The thresholds are placeholders: re-render with F0b-identified
parameters once they exist.

Regenerate from the experiment directory with:

```sh
python3 render_f0b_s_animation.py
```

The script uses the same 64-colour palette approach and 5 MB limit as the L
animation. Its runtime log prints the regime boundaries and the final lag.

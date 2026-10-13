# Motor identification

**Not populated yet.** This folder exists to keep motor-side work out of
[`../stage identification/`](../stage%20identification/), because the two answer
different questions and mixing them is what makes a settling number impossible
to act on.

## The question this side answers

Does the rotor go where the commanded microstep says it should — before the
screw, the nut and the carriage get involved?

A stage-side measurement cannot separate that from the mechanism. If a 0.625 µm
commanded microstep produces 0.4 µm at the carriage, stage identification tells
you *that* it happened; it cannot tell you whether the rotor moved the full
microstep and the mechanism lost the rest, or the rotor never moved that far.
The fixes are completely different, so the measurement has to be.

## What belongs here

| | Measures | Needs |
|---|---|---|
| Static phase/torque curve | holding torque vs electrical angle, detent torque | rotor-side angle, or a torque arm |
| Current fidelity | commanded vs delivered phase current | current probe on a coil |
| Chopper comparison | SpreadCycle vs StealthChop microstep fidelity | rotor angle |
| Interpolation effect | MicroPlyer on vs off, per commanded microstep | rotor angle at high rate |
| Thermal behaviour | current derating and step accuracy vs winding temperature | temperature, long runs |

## What it needs that the rig does not yet have

A **rotor-side angle reference**. The encoder currently fitted is on the screw,
so everything it sees has already been through the coupling, the screw and the
nut — which is precisely the part this side is trying to exclude. Until a
motor-shaft encoder or an equivalent reference exists, the honest options are:

- do motor work on a **decoupled motor** on the bench, with its own encoder, or
- infer motor behaviour indirectly from the stage data, and **say that is what
  was done** — the 2 µm out-of-roundness and the 159.5 µm eccentricity already
  measured on this rig show how much mechanism sits in the way of such an
  inference.

The TMC2209 configuration background that applies to all of this is in
[`../../TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md`](../../TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md);
the shared wiring is in [`../CONNECTIONS.md`](../CONNECTIONS.md).

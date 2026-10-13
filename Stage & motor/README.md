# Stage & motor

Characterising the two things that sit between a commanded microstep and where
the stage actually ends up, separately rather than together.

```
Stage & motor/
├── README.md                ← this document
├── CONNECTIONS.md           ← wiring for the shared ESP32-S3; read this first
├── stage identification/    ← what the mechanism does with a correct command
└── motor identification/    ← what the motor and driver do with a correct current
```

## Why two parts

A settling error measured at the stage is the sum of at least two different
things, and they have different fixes:

- the **motor and driver** may not put the rotor where the commanded microstep
  says — detent torque, phase current accuracy, chopper mode, interpolation;
- the **mechanism** may not put the carriage where the rotor says — screw
  pitch error, eccentricity, backlash, stiction, compliance.

Measured together, the two are confounded, and tuning one while the other
dominates wastes the effort. The split exists so each can be measured against
its own reference.

| | reference | measured against |
|---|---|---|
| **stage identification** | commanded position | the screw-mounted encoder |
| **motor identification** | commanded phase | rotor angle / current |

Both share one ESP32-S3, one TMC2209 and one encoder — see
[`CONNECTIONS.md`](CONNECTIONS.md), which also records a **pin collision between
the two existing sketches** and how it is resolved.

## Where things are

| | |
|---|---|
| [`stage identification/`](stage%20identification/) | the 43-minute identification sequence, its firmware, and the jog/calibration tool |
| [`stage identification/jog/`](stage%20identification/jog/) | run-in oscillation between two set-points, and the microstep calibration oscillation |
| [`motor identification/`](motor%20identification/) | scope and plan; not yet populated |

## Scale, once

| Quantity | Value |
|---|---|
| Full steps per mm | 100 → 1 full step = 10 µm |
| Electrical period | 4 full steps = 40 µm |
| Screw revolution | 200 full steps = 2 mm |
| Microstep at MRES 16 | 0.625 µm |
| Stage travel | 44 mm measured |

Everything in both folders uses **full steps from a zero you set**, never an
absolute machine coordinate: there is no homing and no end stop on this rig.

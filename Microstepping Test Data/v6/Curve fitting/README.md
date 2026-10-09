# v6 curve fitting

Two independent analyses of `../EXPrun.csv`, kept as separate versions rather
than one merged deliverable since they answer different questions.

| | Question | |
|---|---|---|
| [`v1/`](v1/README.md) | *During each commanded move, what does the axis actually do* — including overshoot and individual microsteps a simulation model needs to reproduce? | [process spec](v1/CURVE_FITTING_PROCESS.md) |
| [`v2/`](v2/README.md) | *When nothing is commanded at all*, how much does the measured position vary — i.e. what's the noise floor? | — |

`v1/` fits two curves (idealized "reference" and transient-preserving
"dynamic") to every commanded oscillation and trajectory block, grouped by
MRES. `v2/` looks at the opposite case: the 153.65 s at the end of the
capture after the last commanded move finishes, where the true answer is
"nothing" and any measured variation is either noise or drift, not signal —
see `v2/README.md` for why that distinction matters and how the noise number
was isolated from a slow-wander component and a single transient at the
segment's own boundary.

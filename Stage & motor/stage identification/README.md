# Stage identification

What the mechanism does when the command is taken as correct: gain, static
error against electrical phase, hysteresis, settling, and bidirectional
repeatability — measured against the screw-mounted encoder.

```
stage identification/
├── README.md
├── scripts/
│   ├── sequence_design.py            the campaign, defined as data
│   └── render_sequence.py            draws it
├── sequence/
│   ├── stage_identification_sequence.png
│   ├── sequence_plan.csv             every move, in order
│   └── sequence_summary.json
├── firmware/stage_identification/    executes a plan uploaded by the host
├── chirp/                            mechanical frequency sweeps, 1 Hz - 1 kHz
├── microstepping/                    MRES selection test, 1/16 - 1/256
└── jog/                              run-in oscillation + microstep calibration
```

Two kinds of test live here. The **positioning** sequence below asks where the
stage comes to rest; [`chirp/`](chirp/) asks how it responds to excitation
against frequency. They share the rig and the wiring and nothing else.

> **[`chirp/`](chirp/) carries a safety constraint**: 150-250 Hz produced
> uncontrolled motor motion on this rig on 2026-09-15, and attenuating the band
> was not enough to stop it. The sweep omits the band entirely. Read
> [`chirp/README.md`](chirp/README.md) before running anything there.

---

## 1. The sequence

![the sequence](sequence/stage_identification_sequence.png)

**43 minutes, 230 commanded moves**, within ±10 mm of a zero you set.

| | Phase | What it gets you | How | Time |
|---|---|---|---|---|
| 1 | Gain calibration | the true step → distance gain | long moves of whole electrical periods *and* whole screw revolutions, both directions; fit the slope | 3.6 min |
| 2 | Microstep staircase | static error vs phase, and hysteresis | quasi-static staircase through 2 electrical periods forward then back, 0.75 s per microstep, at 3 locations | 10.1 min |
| 3 | Settling matrix | **the main metric** | trapezoidal moves, randomised, distance × speed × direction, 5 s dwell after each stop | 19.2 min |
| 4 | Bidirectional repeatability | ISO 230-2 style accuracy | 5 targets, approached from both sides with a run-up, 5 repeats | 8.4 min |
| D | Drift checks | separate drift from the effects | return to a fixed reference every ~5 min | interleaved, 7 checks |

Regenerate after changing anything:

```bash
python scripts/sequence_design.py          # -> sequence_plan.csv, summary.json
python scripts/render_sequence.py          # -> the figure
```

### Why the phases are built the way they are

**Gain: whole electrical periods *and* whole screw revolutions.** Both, because
they are different errors. A move that is not a whole number of electrical
periods starts and ends at different points in the motor's electrical cycle, so
the static phase error does not cancel and leaks into the fitted slope. A move
that is not a whole number of screw revolutions starts and ends at different
screw angles, so eccentricity does not cancel either — and on this rig that term
was *measured at 159.5 µm amplitude*, far larger than the gain error being looked
for. 200 full steps is 50 electrical periods and exactly one screw revolution, so
every multiple of 200 satisfies both at once.

**Staircase: forward then back over the same ground.** The difference between the
two passes at the same commanded position *is* the hysteresis; a one-way sweep
cannot recover it. Repeating at three locations separates error that follows the
electrical phase (repeats every 4 full steps wherever you are) from error that
follows the screw (does not).

**Settling: randomised.** The confounders here are all monotonic in time — the
motor warms, lubricant redistributes, the room drifts. Run in a tidy nested
order, any of those would masquerade as an effect of whichever factor varied
slowest.

**Repeatability: a run-up, not just a repeat.** The stage is taken past the
target in the direction of travel and brought back, so backlash is always taken
up the same way within a direction and the two directions can be compared.

**Drift checks never interrupt a block.** A detour to the reference in the middle
of a staircase, or between a run-up and its target, inserts a backlash reversal
and measures something else. Moves that must stay together are tagged with a
group, and a check is only inserted at a group boundary — which is why the first
two checks sit at 6.9 and 6.6 min rather than exactly 5.

---

## 2. Running it

The plan is **not compiled into the firmware**. The host uploads it move by move
and then starts it, so the campaign has one definition — `sequence_design.py` —
and changing it never means reflashing.

```
PLAN CLEAR
PLAN ADD <target> <speed> <dwell_ms> <mres> <micro_n> <micro_dwell_ms> <tag>
PLAN RUN
```

While running, the board streams `# S,<t_ms>,<commanded_full_steps>,<encoder>,<index>`
at 50 Hz, plus `# MOVE` and `# ARRIVED` at each boundary. `ABORT` stops between
pulses.

Before a first run, with the stage roughly centred in its 44 mm:

1. `ZERO` — everything is relative to this, and it is lost on power-cycle.
2. Check the direction sense and the scale with a short move (see `jog/`).
3. Confirm `# READY,driver=1,encoder=1`.

---

## 3. The jog tool

[`jog/`](jog/) is the run-in oscillation, migrated here from the repository root,
with the microstep calibration added:

```bash
python jog/jog_stepper.py --port COM11            # find board, prompt
jog> -200 200                                     # two set-points, full steps
jog> speed 70                                     # 0.7 mm/s, "moderate"
jog> go                                           # shuttle until a key

jog> mres 16                                      # 0.625 µm per microstep
jog> calib 30 2000                                # +1 / -1, 2 s dwells, 30 s
```

`calib` is the smallest commanded motion the rig can make, held still long
enough to be measured. A +1/−1 pair returns to the same commanded position, so
whatever the encoder reports at the end of a pair is the stage's error and not
the command's — lost motion, stiction and the dead-band around a reversal show
up here and nowhere else. `mres` sets the step size being probed.

---

## 4. State

**The sequence design and the figure are verified** — the renders are generated
from the same module the firmware is fed from, and building the figure is what
caught two real design faults: drift checks landing inside a staircase, and
deferred checks firing in a burst.

**Neither sketch has run on hardware.** No board has been attached at any point.
Specifically unexercised: the step timing under load, the direction sense, the
soft limits, the TMC2209 UART readback, and the PCNT encoder configuration on
the new pins. Treat the first run as a bring-up, not a measurement — and see
[`../CONNECTIONS.md`](../CONNECTIONS.md) §1, because the encoder has to be
rewired off GPIO5/6 before either sketch can work at all.

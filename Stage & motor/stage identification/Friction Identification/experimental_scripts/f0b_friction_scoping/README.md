# F0b: friction scoping

The governing specification is [`F0b_scoping_procedure.md`](F0b_scoping_procedure.md).
This folder holds the session plan built from it, the checks on that plan, and
the figure that explains the session. **Status: not run.** Nothing here is
measured; `data/` and `analysis/` arrive with the session.

F0b is the first block in the friction roadmap
([`../../rendered/fig_friction_roadmap.png`](../../rendered/fig_friction_roadmap.png)).
It does not identify a friction model. It sizes every effect the model will
need — how large, over what distance, over what time — so that F1, F2 and F3
can be dense where something happens and short where nothing does.

## The figure

![F0b sequence overview](rendered/f0b_sequence_overview.png)

Four panels, in the form of the E1–E5 sequence renders:

- **(a)** the eight blocks in the procedure's order, with their minutes;
- **(b)** the complete commanded trajectory. The axis is symmetric-log
  because this one session commands displacements from 39 nm to 45 mm; on a
  linear axis everything except the reference pairs is a flat line at x₀;
- **(c)** the same session by *amplitude* instead of by time. This is the axis
  F0b is really built on, and it is the panel that shows the design is sound:
  the ladder runs from one microstep to 76 800 and brackets β₀ by about a
  decade and a half either side, so the saturation level it measures can be
  compared with the frozen reversal constant rather than merely reported;
- **(d)** where the time goes, against the 58 min watchdog gate.

## What the plan says about the session

`f0b_plan.py` builds the command table from the procedure and times it on the
same trapezoidal model the firmware uses. The result:

| | |
|---|---|
| Segments | 343 |
| Predicted duration | **45.3 min** (the procedure estimates ≈ 46) |
| Margin to the 58 min gate | 12.7 min |
| Inside a 60 min scope record | yes, by 14.7 min — no repeat of the E5 stitch |
| Plan checks | 10 of 10 pass |

Per block, in minutes:

| R1 | N | L | R2 | S | D | V | R3 |
|---|---|---|---|---|---|---|---|
| 3.5 | 8.8 | 6.9 | 3.5 | 11.7 | 3.9 | 3.5 | 3.5 |

**The total matches the procedure; the split does not.** §3 budgets 1 min for
each reference pair and ≈ 3 min for transits — 6 min in all. The plan spends
10.4 min, 23 % of the session, on the three reference pairs, and takes it back
from L (6.9 against 10), D (3.9 against 5) and V (3.5 against 5).

The reason is in §2: holding **MRES at 1/256 for the whole session** caps the
step rate at 25.6 kHz, which is 1 mm/s. A reference pair starts at x₀ = 60 mm,
visits 15 mm from both directions, then 105 mm from both directions, and
returns — 198 mm of travel, and at 1 mm/s that is 3.3 min of it before a single
stop is recorded. E5 could transit at 2 mm/s; F0b cannot, because changing MRES
mid-session would move the step-to-position mapping and the electrical phase
reference, which is exactly what §2 forbids.

Nothing is broken by this: 45.3 min clears the gate with 12.7 min in hand. But
if the session needs to be shortened, the reference pairs are where the time
is, not the measurement blocks.

## Files

| Path | Purpose |
|---|---|
| `F0b_scoping_procedure.md` | the governing procedure, as written |
| `f0b_config.json` | every number from §1.1, §2 and §3, in one place |
| `f0b_plan.py` | builds the command table, times it, and checks it |
| `render_f0b_sequence.py` | draws the sequence figure from `sequence/` |
| `sequence/f0b_command_table.csv` | the table the firmware replays (§2.1 format) |
| `sequence/f0b_schedule.csv` | the same rows with their start and end times |
| `sequence/f0b_plan_summary.json` | durations, counts and the checks |
| `rendered/f0b_sequence_overview.png` | the figure above |

```sh
python f0b_plan.py              # rebuild sequence/ and print the checks
python render_f0b_sequence.py   # redraw the figure from sequence/
```

`f0b_plan.py` exits non-zero if any check fails, so it can gate the session.

## What the plan checks

All ten pass as written:

- x₀ is on the 1024-µstep electrical grid, and u₀ = 1 536 000 matches 60.000 mm
- every commanded position is inside the 0–120 mm soft limits
- no segment exceeds the 25.6 kHz driver limit
- the ladder's 3 mm stroke stays inside the 0.75 mm/s cap §3 sets for it
- the stroke list is monotonic, so each stroke really does wipe the last
- one microstep (39.1 nm) is well above the interferometer floor (≈ 2 nm)
- the ladder brackets β₀ = 1.757 µm from both sides
- the session fits the 58 min gate and the 60 min scope record

## Carry-alongs

Three things from the geometric close-out ride along with this session at
almost no cost, and are set out in [`../../README.md`](../../README.md): the
daily two-point block with a BME280 and two NTCs, logging the serial port end
to end, and the B3 re-analysis. The reference pairs this plan already spends
10 min on are the same measurement the two-point block wants, so the first of
the three is most of the way done by the procedure as written.

## A note on the folder name

The folder is named for the experiment it holds. In the roadmap, **F1** is the
reversal sweep — the main experiment — and **F0b** is this scoping session, so
calling this folder `f1_…` would collide with F1's own folder when it arrives.

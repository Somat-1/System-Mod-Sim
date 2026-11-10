# F2 — dwell (frictional ageing and creep)

F2 asks whether the time the stage rests before a move changes how that move
starts. Each trial is the same preparation — a 3 mm wipe from below — then a
rest of T_d, then a 6 µm ramp either against the wipe (mode R, a reversal) or
with it (mode C, a continuation). The dwell is the variable; everything else
is held fixed. The procedure is
[`F2_dwell_procedure.md`](F2_dwell_procedure.md).

**Status: sequence generated, overview figure rendered, command table and
firmware written. There is no data yet; the companion plots of procedure §9
are not rendered.**

Two things are still open, and both are recorded in the artefacts rather than
assumed away:

- **F2 is conditional.** The full grid is justified only if F0b's dwell probe
  (1 s vs 60 s) exceeded 3·d_min. F0b has not been run, so
  `f2_config.json` carries `dwell_probe_delta_um: null` and the full variant
  is generated on the procedure's defaults. `check_f2_sequence.py` turns that
  into a real check as soon as the number exists: fill it in and the checker
  insists the variant matches it.
- **The post-dwell ramp is provisional.** It is 2·s₉₀ long, and s₉₀ comes
  from F0b; the placeholder is 3 µm, so the ramp is 154 µsteps = 6.0156 µm.
  The generated header sets `CFG_PROVISIONAL = 1` and the firmware says so at
  boot. Replace `s90_um` and `d_min_um`, rerun both builders, reflash.

## Files

| Path | Contents |
|---|---|
| `F2_dwell_procedure.md` | the supplied procedure, including the render instructions (§8–§9) |
| `f2_config.json` | positions, dwell grid, ramp rule, pauses, seed, gate and colours; the F0b-dependent values sit under `f0b_dependent` |
| `f2_sequence.py` | builds the timed plan: every transit, wipe, settle, dwell, ramp, hold and reference stop as one segment |
| `f2_make_table.py` | builds the table the board replays, times it, checks it, writes the firmware header |
| `check_f2_sequence.py` | re-derives procedure §5 from the written files alone |
| `render_f2_sequence.py` | the overview figure of §8, drawn from the timed plan |
| `rendered/f2_sequence.png` | the figure, full variant |
| `rendered/f2_sequence_reduced.png` | the same figure for `--variant reduced` |
| `sequence/f2_segments.csv` | the plan, one row per segment, tagged by block, trial, repetition and T_d |
| `sequence/f2_plan_summary.json` | counts, durations, coverage, correlations and the trial order |
| `sequence/f2_command_table.csv` | the table |
| `sequence/f2_trial_order.csv` | the shuffled order, one row per trial, with the dwell and ramp times |
| `sequence/f2_table_summary.json` | durations, counts, hashes and the checks |
| `firmware/f2_dwell/f2_config.h` | generated; do not edit |
| `firmware/f2_dwell/f2_dwell.ino` | the sketch |

```sh
python f2_sequence.py                  # the timed plan, for the figures
python f2_make_table.py                # the table and the firmware header
python check_f2_sequence.py            # the pre-run list, from the files
python render_f2_sequence.py           # rendered/f2_sequence.png

python f2_sequence.py       --variant reduced    # procedure §4.3
python f2_make_table.py     --variant reduced
python render_f2_sequence.py --variant reduced
```

## Plan with the default values

**278 segments, 41 trials, 42.93 min**, of which **21.6 min is spent at
rest**. The gate is 52 min, so there are 9.1 minutes of margin.

| Block | Position | Contents | Time |
|---|---|---|---|
| R1 | — | reference pair, 15 / 105 mm, both directions | 2.72 min |
| B | 60 mm | 34 trials: R and C × T_d = 0.1, 0.5, 2, 10, 30, 120 s, 3 repetitions (120 s: 2) | 24.67 min |
| A | 25 mm | 6 trials: R only × T_d = 0.5, 10, 120 s, 2 repetitions | 6.99 min |
| L | 60 mm | one 300 s dwell, mode R | 5.83 min |
| R3 | — | reference pair | 2.72 min |

Time by what it is spent on: 21.61 min dwelling, 8.22 min on the post-dwell
ramps, 5.44 min on the reference pairs, 3.98 min wiping, 1.97 min in transit,
1.37 min holding after the ramps, 0.34 min settling.

### Against the estimates in procedure §4.2

- **B lands at 24.67 min against ≈ 25.** The dwelling comes out at exactly the
  12.26 min the procedure predicts; the overhead is 23.2 s per trial rather
  than 22, because the transit into a block is charged to the first trial's
  wipe run-up.
- **A is 6.99 min against 6.6, L 5.83 against 5.5.** Same reason: §4.2 lists
  the 0.6 min transits separately, the builders fold each into the run-up of
  the trial that needs it.
- **The reference pairs are 2.72 min, not 3.5.** Four stops with 3 mm
  approaches and transits at 1 mm/s, as measured in F1 for the same block.
- **Total 42.93 min against ≈ 45–46.** The difference is almost entirely the
  two reference pairs.
- **The reduced variant** (§4.3) is 14.00 min with both reference pairs and
  7.81 min without them, which brackets the procedure's ≈ 12. Appending it to
  another session means the host session supplies the references; set
  `reduced_includes_references` to `false` for that.

## The overview figure

![F2 command sequence](rendered/f2_sequence.png)

`render_f2_sequence.py` draws procedure §8 from the timed plan, not from a
second description of it, so the figure cannot disagree with the table the
board replays: both come out of `f2_config.json`, and the panel (a)
trajectory is literally the segment list. It re-runs the plan checks and
exits non-zero if any of them fail, so a figure that renders is a figure of a
plan that passes.

The three panels are §8's: **(a)** the ribbon and the commanded trajectory
relative to the current test position, dwells as flat bars on the zero line
shaded by duration; **(b)** the dwell times covered at each position, filled
for reversal and open for continuation, against the dwells already present in
F0b, F1 and E5; **(c)** the anatomy of one trial per mode, aligned at the end
of the dwell.

Where it departs from §8, and why:

- **The ribbon's lower tier shows the repetitions of every block, not only
  B's.** A has two and L has one; drawing them costs nothing and makes the
  ribbon's two tiers mean the same thing all the way across.
- **The dwell colour scale starts at a visible green**, not at near-white.
  With the palest shade a filled 0.1 s marker in (b) was not tellable from an
  open one. A mid stop is also given explicitly, because the two end colours
  interpolated alone put everything below 10 s in the palest quarter.
- **(c) is not to scale and says so on its axis.** The two drawn lengths are
  2 s and 120 s; both end at t = 0, the long one is broken with //, and the
  short one is drawn as the lighter bar nested inside it, so the reader sees
  two dwells ending at the same instant rather than one bar of ambiguous
  length. The 300 s dwell is in (a) and (b), where it is to scale.
- **No position map.** F1's panel (c) placed its three test positions along
  the travel against the first-harmonic wobble; F2 reuses two of those same
  positions and adds nothing to that picture, and §8 does not ask for it.
- **(a) reads the dwell times, not the dwell shades.** At 43 min across the
  page a 2 s dwell is about one pixel wide, so the shade is only legible on
  the 120 s and 300 s bars. That is the point the panel is making — the short
  dwells are invisible because they *are* short — and (b) carries the grid.

`--variant reduced` renders the §4.3 session to
`rendered/f2_sequence_reduced.png` (8 trials, 80 segments, 14.0 min).

## What is checked, and where

The post-dwell ramp is the measurement, so its direction has to be right for
every one of the 41 trials, and the dwell has to be the only thing that
differs between them. Three layers check this, and each can fail on its own:

1. **While the plan is built.** `wipe()` asserts that the approach ends moving
   up, and `ramp()` asserts `ramp_dir == −last_dir` for mode R and
   `== +last_dir` for mode C. The dwell between them is motionless, so
   `last_dir` still describes the wipe.
2. **After the table is written.** `recheck_directions()` and
   `check_ramp_types()` recompute every direction from the positions in the
   finished table, ignoring how it was built, and also assert that each dwell
   is motionless and preceded by an upward wipe.
3. **From the files alone.** `check_f2_sequence.py` imports neither builder.
   It reads the command table, the trial order, both summaries and the
   firmware header, and re-derives procedure §5: seed logged, order printed,
   every (T_d, mode) present the planned number of times and exactly once per
   repetition, the order actually shuffled, the dwell grid whole in
   milliseconds, the directions correct, the header matching the table row for
   row, and the session inside the gate. It reports **40/40** and returns a
   non-zero exit status otherwise.

`f2_make_table.py` additionally cross-checks itself against
`sequence/f2_plan_summary.json`: segment count, trial count, duration to
0.01 min, per-block times to 0.01 min, the ramp length and **the whole trial
order**. The two builders are written separately, including two separate
copies of the seeded shuffle, so agreement is evidence that the grid in
`f2_config.json` was read the same way twice. As in F1, that is the point of
having two.

### Randomisation

The (T_d, mode) list is shuffled within each repetition from a logged seed
(`shuffle_seed: 20261110`), with one sub-seed per (block, repetition):
`seed + 1000003·block + 101·rep`. The generator is SplitMix64 feeding a
downward Fisher–Yates, written out in both builders rather than taken from
`random`, so the order is reproducible from the seed alone on any
interpreter and can be re-derived by the analysis.

What shuffling is for is checked rather than assumed: the correlation between
log₁₀(T_d) and when the trial starts is **−0.165 at B and +0.132 at A**,
against a configured limit of 0.45. A trend in time therefore cannot be read
as an ageing effect, and a residual trend in the data is a thermal or drift
question for the R1/R3 pair, not a dwell question.

## Assumptions not fixed by the procedure

- **Block order** is §4.2's: R1 → B → A → L → R3. L is at 60 mm, so the
  session transits back from A; that transit is inside L's first run-up.
- **The 120 s dwell is dropped from repetition 3**, which is how "3 (120 s:
  2)" is read. B therefore has 12 + 12 + 10 = 34 trials.
- **Start and references.** The session starts homed and parked at
  x_B = 60.00 mm, and each reference pair visits the nearer reference first,
  as in F1.
- **Where t = 0 of the dwell is.** The procedure makes the 0.5 s settle part
  of every trial and puts t = 0 at its end, so the table holds the settle and
  the dwell as two separate segments and the board times them separately.
  Physically the rest begins when the wipe stops: the total rest is
  0.5 s + T_d, identical in every trial. For the creep fit of §6.2 the stop
  instant is the **start of the settle segment**, which the firmware log
  timestamps; the 0.2 s datum window of §6.1 sits inside that settle.
- **Ramp rounding.** 2·s₉₀ = 6 µm becomes 154 µsteps = 6.0156 µm, symmetric
  in the sense that mode R and mode C use the same length in opposite
  directions. At 0.5 µm/s that is 12.80 Hz and 12.03 s.
- **The slow path carries the data.** Every measurement ramp in F2 is at
  12.8 Hz, below the ~33 Hz floor where one period still fits an RMT symbol,
  so all 41 ramps run on the timed single-pulse path. F0b and F1 both use it,
  F1 down to 2 Hz.
- **The motor is never released between trials.** IRUN = IHOLD,
  TPOWERDOWN = 0 and IHOLDDELAY = 0, so the current during a 300 s dwell is
  the current during motion; `configureDriver()` reads this back and refuses
  to run if the driver disagrees, and a clean run leaves the motor energised.
  Releasing it would wipe the hold state the dwell is supposed to age.

## The dwell timer

Procedure §5 asks for the dwell timer to be verified to better than 10 ms
before the session. The sketch does three things about it:

- every DWELL and MARK is held against a single `micros()` deadline rather
  than a 2 ms `millis()` poll, so the error does not grow with T_d;
- the achieved duration is recorded per segment and printed in the `SEG`
  lines as `held_us` and `dwell_err_us`, and the run ends with a
  `DWELL_CHECK` line giving the worst error over the 41 trial dwells against
  `CFG_DWELL_TOL_MS = 10`;
- **`d` runs the check on its own**, without moving the stage: it times a
  100 ms, a 2 s, a 30 s and a 300 s hold, reports each error, and asserts the
  commanded position did not change. Run it once before the session and the
  pre-run box is ticked with a logged number.

The watchdog is tested between segments, so a dwell can overrun the 52 min
gate by at most its own length. Both builders check that
`gate + longest dwell ≤ scope record`: 52 + 5 = 57 min against the 60 min
record.

## Running the session

> **Provisional and conditional.** Re-run `f2_make_table.py` after F0b is
> processed and reflash before treating an F2 run as final. If F0b's dwell
> probe came out below 3·d_min, run `--variant reduced` instead: it bounds
> the effect in 14 min and can be appended to another session.

Preconditions the firmware prints at boot and expects to be true: E0 warm-up
complete, IDS recording at 1 kHz with a scope record of at least 60 min, F0b
processed, the dwell probe checked, and temperature logging running — long
dwells are where thermal drift and ageing are easiest to confuse, which is
why §6.1 asks for the thermistor log as well as the R1/R3 pair.

The segment types, the two step paths and the bench procedure are the same as
F0b's and F1's — see
[`../F0b Friction Scoping/README.md`](../F0b%20Friction%20Scoping/README.md)
and [`../F1 Reversal Sweep/README.md`](../F1%20Reversal%20Sweep/README.md).
F2 emits no `step` segments; that path is retained only so the motion engine
stays identical to the proven one.

## Not done here

- The **companion plots D1–D4** of §9 (creep during the dwell, reversal after
  dwell, breakaway when continuing, effect against dwell time). They are
  illustrative curves of what the data should look like, as F1's V1/M1/Q1
  are, and need the F1 pre-sliding solver extended with an ageing term before
  they would show anything honest. The plan and table carry what they need
  from the sequence side: `f2_segments.csv` has the trajectory with each
  segment's kind, trial, repetition and T_d, and `f2_trial_order.csv` has the
  shuffled order with the dwell and ramp timestamps.
- Any analysis script. The outputs of §7 (`f2_creep.csv`,
  `f2_reversal.csv`, `f2_breakaway.csv`, `f2_summary.json`) come after the
  session.

# F0b_pt2 — reversal vs microstep resolution, and speed

F0b measured a saturation lost motion of 2.30 µm at 1/256 with 0.5 s pauses.
E2–E5 measured β₀ = 1.76 µm at 1/16 with 2 s settles at 2 mm/s. F0b_pt2 varies
exactly those three factors — resolution, settle time, approach speed — one at
a time in one session, so the 30 % gap is attributed rather than argued about.
The procedure is [`F0b_pt2_procedure.md`](F0b_pt2_procedure.md).

**Status: run once on 12 Nov (complete, but the windows did not settle).
Sequence revised for a second run. Compiled; not yet run.**

## Revision for the second run

In the first run every arrival had a 2.0 s settle and a 0.5 s mark. The mark
ended at the moment the next move began, and the analysis read W1
(0.3–0.5 s) and W2 (2.0–2.5 s). Under the 20 nm/s settling screen, W2
passed at about 75 % in M and almost nowhere in V. There were three
problems:

1. **W1 was screened.** At 0.3–0.5 s after the stop the stage is still
   creeping at about 100–150 nm/s, and no pause length changes that, because
   W1 is fixed relative to the stop. W1 measures the creep. It is a
   measurement, not a settled reading.
2. **W2 was marginal.** The log creep fitted to F1 (x = a·ln t) still gives
   about 20 nm/s at 2.0–2.5 s after a saturated reversal.
3. **No guard before the next move.** The scope and ESP clocks are aligned by
   a single offset (`analysis/run_data.py`, `best_lag`), so a clock-rate
   error grows through the session. A window that ends exactly where the next
   move starts picks up that move once the error is a few tens of ms. This
   hits V hardest, because V is the late block.

The revised stop, all offsets measured from the last step edge:

| Window | Offset from last edge | Use |
|---|---|---|
| W1 | 0.3–0.5 s | creep measurement, **not screened** |
| W2 | 2.0–2.5 s | same window as the first run, for comparison, **not screened** |
| **W3** | **6.0–6.5 s** | **the settled reading; the only window screened** (the MARK) |
| guard | 6.5–6.8 s | nothing moves; absorbs up to ~110 ppm of clock-rate error over the session |

Under log creep, W3 has 2.8 times less slope than W2. A turn that read
20–55 nm/s at W2 should pass at W3. One stop contains all three windows, so
the run measures the creep curve of every arrival as well as settling it.

- **Session:** 1149 segments, 248 arrivals, **42.2 min**, inside the 52 min
  gate, 58 min watchdog and 60 min record. The checker passes 31/31.
- **Self-test `c` fixed:** the test now compares the measured time with the
  expected time of the same symbols, and reports the real cruise rate (50.0 kHz
  for a 51.2 kHz target, because of the 1 µs period rounding) and `pass=1/0`.
  The first run's bench reading gives a timing ratio of 0.991, a pass.
- **Still recommended in the analysis:** fit a clock rate as well as an
  offset when aligning (lag at the start and at the end of the run, both
  well defined by R1 and R3). The guard makes the sequence safe without this
  fit, but the windows land more accurately with it.

## Files

| Path | Contents |
|---|---|
| `F0b_pt2_procedure.md` | the supplied procedure |
| `f0b_pt2_config.json` | positions, MRES and speed orders, overruns, windows, gate |
| `f0b_pt2_make_table.py` | builds the table the board replays, times it, checks it, writes the firmware header |
| `check_f0b_pt2_sequence.py` | re-derives procedure §5 from the written files alone |
| `sequence/f0b_pt2_command_table.csv` | the table, one row per segment |
| `sequence/f0b_pt2_arrivals.csv` | the 248 settled arrivals the analysis reads |
| `sequence/f0b_pt2_summary.json` | counts, durations, hashes, checks |
| `firmware/f0b_pt2/f0b_pt2_config.h` | generated; do not edit |
| `firmware/f0b_pt2/f0b_pt2.ino` | the sketch |

```sh
python3 f0b_pt2_make_table.py        # table, arrivals, summary, firmware header
python3 check_f0b_pt2_sequence.py    # the pre-run list, from the files
```

## The session as first generated (12 Nov run)

**901 segments, 248 arrivals, 24.4 min** against a 52 min gate — 27.6 minutes
of margin, so the §4 fallback of dropping the 1/32 and 1/128 visits is not
needed. The estimate in §4 is ≈ 45 min; the difference is that the moves are
quicker than that estimate assumed, not that anything is missing.

| Block | Contents | Time |
|---|---|---|
| R1 | reference pair, 15 / 105 mm both directions at 1/16, E5 speed | 2.62 min |
| M | 10 MRES visits, ABBA 16·32·64·128·256·256·128·64·32·16 | 12.92 min |
| V | 8 speed visits at 1/256, ABBA 2·1·0.1·0.01·0.01·0.1·1·2 mm/s | 6.22 min |
| R3 | reference pair | 2.62 min |

Each M visit is one MRES switch, 7 saturated arrivals at x₀ alternating from
below and above with 3 mm overruns (6 reversal values), then anchored 5 µm and
20 µm loops, 3 cycles each. Each V visit is a 3 mm wipe then 5 arrivals with
100 µm overruns (4 reversal values). Every arrival is a 2.0 s settle and a
0.5 s record window.

Step rates run from **256 Hz** (0.01 mm/s at 1/256) to **51 200 Hz** (2 mm/s at
1/256). Everything is above the ~33 Hz RMT floor, so every measurement segment
uses the buffered path.

## The four firmware changes of procedure §0

**1. One internal unit.** `positionU` is in 1/256 counts at every MRES, and so
is every target in the table. A segment's travel is converted to driver steps
by dividing by `uStepPerDriverStep = 1 << mres_code`. A travel that does not
divide is a fault that stops the run, not a rounding — the builder asserts it
for every segment against the resolution actually in force at that point, and
the firmware asserts it again at run time. Rates are stored as milli-counts per
second, so **a velocity is the same number at every MRES** and only the step
rate changes, which is what makes MRES the only variable in block M.
Acceleration is converted the same way.

**2. MRES switch at rest, with the MSCNT check.** `applyMres()` reads MSCNT and
CHOPCONF, writes the MRES field, reads both back, and requires that the field
took and that **MSCNT did not move**. MSCNT is the microstep table position
0–1023 at every resolution, so if it shifts across a switch the coil-current
vector moved with it and the comparison between levels is void: the run aborts
with `mres_fault=1`. Both readings go into the log for every switch. There are
13 switches, one at the head of every M visit including the repeat at the same
level, because the check is the point of it.

**3. Windows timed from the last step edge.** `lastEdgeUs` is stamped when the
final pulse has actually been transmitted. `end_offset_ms` in the table is not
a duration from when the segment was issued but the offset from that edge at
which the segment ends: the settle ends 2.000 s after the stage stopped and the
mark 2.500 s after it. One dwell therefore contains both analysis windows,
0.3–0.5 s and 2.0–2.5 s, measured from the same instant. `t_last_edge` is also
logged per segment.

**4. Telemetry in RAM**, printed after `SEQ_DONE`, as before.

This is the fix for what F0b could not answer. It does not make mid-move
readings usable — that needs option A or B of §0, which are hardware — but it
removes the ambiguity from every settled reading, which is what this session is
built on.

## What is checked, and where

`check_f0b_pt2_sequence.py` imports no builder. It reads the table, the
arrivals, the summary and the generated header, and reports **30/30**:

- the config and table hashes agree across the summary and the header, and the
  header matches the table row for row;
- every move divides into whole steps at the resolution in force, every turning
  point is a multiple of 0.625 µm (one 1/16 step), and x₀ and the ±3 mm
  saturated points are on the 40 µm electrical grid;
- every declared approach side matches the direction the move actually goes;
- every settle ends at 2000 ms and every mark at 2500 ms from the last edge,
  and both analysis windows fit inside one dwell;
- the MRES order is the ABBA palindrome, every level is visited twice, every M
  visit carries its own switch, and every switch follows a dwell or a mark so
  it happens at rest;
- 7 saturated arrivals per MRES visit and 5 per speed visit, sides alternating
  below/above, 120 partial-stroke arrivals;
- the fastest rate is within the bench limit, the slowest is above the RMT
  floor, and the session fits the gate inside the watchdog inside the record.

Three things it cannot check, because they are bench facts:

- **that the RMT path really delivers 51.2 kHz.** The sketch's `c` command
  measures it without moving the stage: the driver output is disabled, a
  quarter-second train is sent and timed, the same number of pulses is sent
  back, and it reports the achieved rate and whether MSCNT returned.
- **that MSCNT survives a live MRES switch.** The run aborts if it does not.
- **that a 1 mm move measures 1 mm at every resolution.**

## If the 51.2 kHz check fails

Set `m_speed_mm_s` to 1.0 and `v_extra_coarse` to true in
`f0b_pt2_config.json`, then rerun both scripts and reflash. Block M then runs
at 1 mm/s throughout and block V gains a 1/16 at 2 mm/s visit, so the E5
condition is still covered somewhere in the session. The checker's
`fastest_rate_within_bench_limit` will then pass at 25 600 Hz.

## Deviations from the procedure, and what is not here

- **No second independent builder.** F0b, F1 and F2 each have a plan builder
  and a table builder written separately and cross-checked. Here there is one
  builder and a checker that reads only the written files. The properties that
  matter for this session are different from F1's — divisibility at each MRES,
  sides, switch-at-rest, window offsets — and the checker tests all of them
  from the artefacts rather than from the code that made them. Worth adding the
  second builder before the run if you want the same belt and braces as F1.
- **Not compiled.** No Arduino toolchain here. The sketch was checked
  structurally instead: brackets balanced, every `CFG_` symbol it uses defined
  in the generated header, every segment enum defined, and each of the four §0
  changes present and wired. The motion engine is F0b's, changed only where
  §0 requires.
- **No figures yet.** The §8 overview and results figures are not drawn.

# Settling offset against microstep size

**Written for: whoever picks up this folder — the campaign's rationale, what it
measures, and how the data turns into an answer.**

## The question

The stage is commanded to move a distance. It ends up somewhere. **Does the
microstep size change how far it actually ends up moving?**

That is the only question here. Not where the stage is relative to a datum, not
how repeatable it is, not how it gets there. The deliverable is one number per
(microstep size, distance) cell:

```
displacement = x_settled_after - x_settled_before
offset       = displacement - commanded
```

Both terms are settled positions a few seconds apart, so **every number in this
campaign is a difference**. Nothing needs an absolute reference, a datum, or a
drift correction. That is deliberate, and it is why this campaign is separate
from `../iso230-2/` — ISO 230-2 measures positioning accuracy against a target
in machine coordinates, which is a different quantity with different error
sources.

## Why the earlier microstepping run could not answer it

The executed MRES sequence (`terminal_2026102*.log`, now in git history) was
built around trajectory execution: staircases and ramps, sampled to show *how*
the stage moves. Block B4 — the settling matrix, the only block that produced a
settling offset per distance — was dropped from the run. So the recorded data
shows trajectory shape, not landing error. Hence this rebuild from scratch.

## Microstep switching was verified — the identical step sizes were by design

The earlier plots show **the same physical step size whatever MRES was set**,
which looks like MRES switching silently failing. It was not. Two separate
findings, both checkable:

**1. The driver really did change microstep size.** `MSCNT` is the driver's own
index into the sine table — 1024 entries per electrical period, 256 per full
step — so it can only land on the grid the driver is actually stepping on. A
`CHOPCONF` readback proves only that the register holds the code; `MSCNT`
proves the stepping followed it. Run on the recorded sessions
(`scripts/verify_mres_switching.py`), every block that executed a plan row
comes back clean: 99.0–99.9 % of samples on a single coset of the correct
spacing — 8 counts at 1/32, 4 at 1/64, 2 at 1/128.

The residual 0.1–1 % are the few samples straddling an MRES change, sitting on
a coset half a microstep away. That matters methodologically: testing the
**GCD** of `MSCNT` differences makes those few samples halve the apparent
spacing and report a confident false failure for every size at once. The script
counts residues instead, which separates "the grid is wrong" from "a few
samples straddle a phase change".

**2. The identical step sizes come from the old plan, not the driver.** The
previous design (`mres_design.py`) set

```python
B3_COARSE_INC = U // 16           # 1/16 FS, the same spacing for every MRES
```

so the coarse staircase tier commanded **0.625 µm at every MRES**, built from a
different number of microsteps each time — 2 at 1/32, 4 at 1/64, 8 at 1/128, 16
at 1/256. Equal physical increments were the intent. Only the fine tier scaled
with MRES, at one microstep: 0.312 / 0.156 / 0.078 / 0.039 µm. Any plot
dominated by the coarse tier shows four identical staircases, correctly.

This also explains an earlier misreading of a block as 1/16 from its ~625 nm
increments: 625 nm is the coarse-tier increment at *every* size, so it
identifies nothing.

**Carried into this campaign.** Register readback is not enough on its own, so
the firmware now proves the step size physically before it will run: menu item
`7`, and automatically ahead of the matrix. It pulses once and five times at
each size and checks `MSCNT` advanced by exactly `256/MRES` per pulse and
closed back to its starting phase. A failure refuses the run — comparing
microstep sizes is meaningless if the sizes are not real. It also catches
`intpol` being on, which would interpolate every size to 1/256.

## Design

`scripts/settling_design.py` is the campaign. It is the single source of truth;
the CSV, the summary and the firmware header are all generated from it.

| | |
|---|---|
| Microstep sizes | 1/16, 1/32, 1/64, 1/256 |
| Distances | 1/16 FS, 1/4 FS, 1 FS, 4 FS, 16 FS, 64 FS, 5 mm, 20 mm |
| Repeats | 4 per direction (2 for 20 mm, travel-capped) |
| Measured displacements | 240 — 60 per microstep size |
| Rows | 344 |
| Duration | 32 min |
| Speed | 78 FS/s = 0.78 mm/s, **identical for every size** |
| Dwell at every stop | 3.0 s; settled value = mean of the final 1.0 s |
| Travel used | ±20 mm about the centre of the stage |
| Order | shuffled within each block, seed 4201 |
| Current | 400 mA RMS, SpreadCycle, `intpol` off, no standstill reduction |

### Distances in microsteps

Every distance is a whole number of microsteps at every size — the design
refuses to emit a plan where it is not. The 1/16 FS distance is **one
microstep** at 1/16 and sixteen at 1/256, which is the extreme of the
comparison and the cell most likely to separate the sizes.

| distance | µm | 1/16 | 1/32 | 1/64 | 1/256 |
|---|---:|---:|---:|---:|---:|
| 1/16 FS | 0.625 | 1 | 2 | 4 | 16 |
| 1/4 FS | 2.5 | 4 | 8 | 16 | 64 |
| 1 FS | 10 | 16 | 32 | 64 | 256 |
| 4 FS | 40 | 64 | 128 | 256 | 1024 |
| 16 FS | 160 | 256 | 512 | 1024 | 4096 |
| 64 FS | 640 | 1024 | 2048 | 4096 | 16384 |
| 5 mm | 5000 | 8000 | 16000 | 32000 | 128000 |
| 20 mm | 20000 | 32000 | 64000 | 128000 | 512000 |

### Three decisions worth knowing about

**One speed for all four sizes.** Settling depends on approach velocity, so
letting each size run at its own top speed would compare speeds, not microstep
sizes. The finest size sets the speed for everything: 1/256 at the firmware's
20 kHz pulse ceiling gives 78 FS/s. This is **below the 5 mm/s operational
speed** — a real limitation, stated in the design output every time it runs.
Dropping 1/256 from the set would allow 312 FS/s. The trade was taken in favour
of covering 1/256, because the question is about microstep size.

**Walks, not returns.** Within each distance the stage walks out n steps and
back, straddling the centre, rather than returning to a start point between
moves. That gives a run of same-direction displacements with the reversal
isolated to the first move of each run (tagged `reversal` in the CSV), instead
of a reversal contaminating every measurement. Straddling rather than walking
out one way halves the travel the plan needs — 40 mm instead of 84 mm — and
changes nothing, because every number is a difference.

**Conditioning before each block.** Three ±2 mm strokes then an approach to the
start from below, so every block begins in the same friction state. Without it
the first few cells of each block carry whatever the stage was left in.

### Block structure

Four blocks, one per microstep size, 86 rows / 60 measured / 8.0 min each:

- **M1** — write `CHOPCONF.MRES` at a full step, read it back, log `MSCNT`,
  then a distinctive marker excursion (32/36/40/48 FS, one per size) so the
  block boundary is unmistakable in the recording.
- **M2** — conditioning strokes, then the start position from below.
- **M3** — the matrix: distances in shuffled order, each walked out and back.

## Files

```
scripts/settling_design.py      the campaign as code; prints the tables, writes the CSV
scripts/gen_settling_header.py  CSV -> firmware/settling_test/settling_plan.h
scripts/settling_terminal.py    host terminal; logs everything to data/
scripts/verify_mres_switching.py  MSCNT check that the step size followed the register
sequence/settling_plan.csv      344 rows, one per commanded action
sequence/settling_summary.json  the numbers above, machine-readable
firmware/settling_test/         the sketch and its generated plan header
data/                           terminal logs and the separately recorded encoder
analysis/                       results (empty until the run happens)
```

Change the design, rerun `settling_design.py`, then rerun
`gen_settling_header.py`. The header says it is generated; do not edit it.

## How the data becomes an answer

The board cannot see the stage — the encoder is recorded separately — so the
firmware's job is to command the moves and mark the recording. Each stop emits,
at the moment its dwell ends:

```
# STOP,<row>,<tag>,d=<u256>,dir=<+-1>,measured=<0|1>,pos=<u256>,mscnt=<n>
```

and position samples stream at 50 Hz as `# S,<ms>,<pos_u256>,<mscnt>,<row>`.

Analysis, once the data exists:

1. **Align clocks.** The export header's FILETIME and the terminal PC's clock
   disagree — by about 16 s in the previous sessions, two machines with two
   clocks. Cross-correlate the commanded position from `# S,` against the
   recorded encoder; the earlier sessions reached r = 0.994 this way. Do not
   trust the timestamps.
2. **Cut the blocks by contiguity**, not by microstep size alone. Grouping by
   size alone spliced separate 1/256 runs together last time. Split on a gap
   > 120 s or a size change.
3. **Settled value per stop** = mean of the encoder over the final 1.0 s of
   each 3.0 s dwell, located from the `# STOP,` line.
4. **Displacement and offset** per measured row, as above.
5. **Per cell** (size × distance × direction): mean offset, SD, n.
6. **The answer**: mean |offset| against distance, one line per microstep size.
   If the lines lie on top of each other, microstep size does not affect
   settling at this speed and current, and the coarsest size that resolves the
   needed increment is the right choice.

Keep the reversal-tagged rows separate in step 5 — they carry backlash and
should not be averaged with the same-direction moves.

## Status

The 344-row matrix was completed on 2026-10-22. The independent encoder export
and matching controller log are in `data/`; the dedicated analysis, PNG
figures, machine-readable results, and findings are in
`analysis/mres_debug/`. See
`analysis/mres_debug/MRES_DEBUG_RESULTS.md` for the result.

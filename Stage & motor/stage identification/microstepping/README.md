# Microstep-resolution test

Which MRES to run the identification sequence (Steps 0–6) at. The test measures
the stage's small-move response with the IDS 3010 at 1/16, 1/32, 1/64, 1/128
and 1/256 in SpreadCycle. The full procedure is in
[`PROCEDURE.md`](PROCEDURE.md). This folder turns the part of it that is kept
into a plan.

```
microstepping/
├── README.md
├── PROCEDURE.md                  the execution procedure, as originally written
├── scripts/
│   ├── run_mres.py               terminal menu: pick a mode and run it
│   ├── mres_design.py            the session, defined as data
│   └── render_mres.py            draws it
├── sequence/
│   ├── mres_test_sequence.png    the figure below (600 dpi, 160 mm)
│   ├── mres_test_sequence.pdf
│   ├── mres_plan.csv             every row, in order
│   └── mres_summary.json         §1 parameters, timing, checks
├── firmware/                     what the runner sketch has to support
├── data/                         one folder per run, one log per block: B<b>_MRES<m>.log
└── analysis/                     §6 metrics and plots, once there is data
```

```bash
python scripts/run_mres.py         # the menu
```

```
================================================================
 MRES test - microstep-resolution selection
 part 1 warm-up: not defined yet   |   part 2: 6 blocks, 30.4 min
================================================================
  1  Warm-up            part 1: homing to the end switches + warm-up
  2  Stepping sequence  part 2: x_ref baseline, six MRES blocks, end reference
  3  Full session       part 1 then part 2
  4  Single block       B1-B3 at one MRES, for bring-up or a repeat
  5  Preview            the plan and the exact commands, no hardware
  6  Simulate           the stepping sequence against a protocol model
  7  Regenerate         rebuild the plan CSV, summary and figure
  q  Quit
```

Any mode can also be started directly, e.g. `--mode stepping`, or
`--mode block --mres 64`. Add `--simulate` to use the protocol model instead
of a board.

- **1 and 3** report that part 1 isn't defined yet and return to the menu.
- **2 and 4** ask for the serial port, check the firmware and set and read back
  the §2 driver settings. They then wait for a typed `RUN` before anything
  moves.
- **Each segment** (reference, each block, end) is uploaded and run as its own
  plan.
- **Ctrl+C** sends `ABORT`.

Each run writes `data/<stamp>_<mode>/` with `session.json` (settings, §1
parameters, board readback), `plan.csv`, and one log per segment.

To regenerate the plan and figure by hand:

```bash
python scripts/mres_design.py      # -> mres_plan.csv, mres_summary.json
python scripts/render_mres.py      # -> the figure
```

---

## Scope: what is kept from the procedure

| Phase | Kept? | Why |
|---|---|---|
| B1 configure + re-zero | yes | MRES can't be switched safely without full-step alignment and the re-zero; ~5 s |
| B2 condition | yes | sets a defined friction state, so B3 doesn't depend on what happened before |
| B3 two-tier staircase | **in full** | the core of the test |
| B4 settling matrix | **dropped** | final positions sit on a common grid for every MRES, so B4 could only detect an MRES effect in the approach dynamics; the ISO baseline at the chosen MRES covers the distance range later |
| B5 creep / hold | **dropped** | depends on friction and hold current, at most weakly on MRES; measured once at the chosen MRES during the friction identification |

The x_ref visits are kept, as a baseline at the start of the stepping part and
as the end-of-session reference. §6 drift subtraction interpolates between them,
and with B4 gone there are no other reference visits.

PROCEDURE.md is left as written. This table is the record of what departs from it.

## One session, two parts

The procedure splits the blocks over two sessions because, with B4 and B5,
each block took ~16 min. Without them, all six blocks fit in **one 30.4 min
session**. The two anchor runs are kept, now as the first and last block:
1/256 brackets the other four, and the four between them run in shuffled
order (seed 7). §7.1's validity check then compares the start of the session
with the end. If the two 1/256 runs disagree, drift or warm-up over the session
is large enough to bias the comparison, and that has to be fixed first.

**Part 1, homing and warm-up, is not defined yet.** `warmup()` in
`mres_design.py` is an empty stub, waiting for the two end switches. It will
home to both ends, set x_start at the centre of the measured travel, and run
the ±5 mm warm-up. Its one obligation is to leave the stage **at x_start, on a
full step**, because part 2 starts from there.

**Part 2, the stepping sequence**, is defined:

| | Segment | What happens | Rows |
|---|---|---|---|
| REF | reference | MRES set for the first block, session marker, x_ref ×3 from below | 10 |
| B1 | configure | at a full step: write MRES, read it back, 2 s settle, IDS marker | 3 |
| B2 | condition | 3 strokes of ±2 mm, then x_start from below, 2 s | 7 |
| B3 | staircase | fine tier (single µsteps) + coarse tier (1/16 FS) forward, the same in reverse; 2 repeats, with a fixed-side re-approach between them | 9 |
| E | end | end marker, x_ref ×3 from below, back to x_start | 9 |

B1–B3 repeat for each of the six blocks. Each segment starts
and ends at x_start and is one plan upload (19 rows per block, well under the
firmware's 512).

---

## The sequence

![MRES test sequence](sequence/mres_test_sequence.png)

**(a)** The session. Part 1 is a placeholder, and the anchor 1/256 (bold) runs
first and last.
**(b)** One 1/256 block end to end. The position axis is symlog, so the ±2 mm
strokes and the 40 µm staircase fit on one axis, and the time axis is square-root
scaled so B1 and B2 stay visible next to B3. **(c)** The B3 staircase at every
MRES over identical displacements. **(d)** Microstep size against the
interferometer floor, labelled with the number of fine-tier microsteps.

| Block | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| MRES | **1/256** | 1/32 | 1/16 | 1/128 | 1/64 | **1/256** |
| duration [min] | 7.4 | 3.2 | 2.9 | 5.0 | 3.8 | 7.4 |

Part 2 is 30.4 min, including 0.3 min each for the start and end references.
B3 takes most of every block. Part 1 has roughly 30 min left in a one-hour
session.

### §1 parameters, as filled in

| Symbol | Value | Source |
|---|---|---|
| L | 2 mm | rig constant, 200 FS per revolution |
| s_FS | 10 µm | |
| δ | 625 / 312.5 / 156.25 / 78.125 / 39.06 nm | 1/16 … 1/256 |
| R | 1.71 µm | ISO 230-2 reversal value, `analysis/stage_results.json` |
| S_fine | max(2.5, 5.13) → **5.625 µm (9/16 FS)** | rounded up to the 1/16 FS grid |
| σ_floor | **3 nm — placeholder** | needs the Allan minimum from the IDS |
| v | 0.7 mm/s (70 FS/s) | the rig's MODERATE, same for every MRES |
| f_step at 1/256 | 17.9 kHz | under the firmware's 20 kHz ceiling |
| x_start / x_ref | 0 / −1 mm | x_start to come from homing (part 1) |
| fixed-side run-up | 0.1 mm, from below | ≫ R |

S_fine is rounded up to a whole number of 1/16 FS steps. That makes the fine
tier an integer number of microsteps at every MRES (9, 18, 36, 72, 144). The
coarse tier then continues in 1/16 FS steps without a remainder, 55 increments
to 4 FS. Because 3R is larger than ¼ FS, the fine tier is longer than the
procedure's worked example of 64 microsteps at 1/256.

Every position in the plan is an integer number of **u256 = 1/256 FS**. Every
MRES divides 256, so every target lies on every block's grid and the spans are
identical by construction. The design rejects an MRES change at any position
that is not a whole full step.

The furthest excursion is now 2.1 mm (a B2 stroke plus its run-up), so the
firmware's default ±11 mm soft limit is enough.

### Interferometer synchronisation

The IDS records separately, so the plan carries **separator signatures** in
the convention used for the Microstepping Test Data v2 runs. Each one is a rapid
negative leap, a 1 s dwell, a return and a 0.5 s settle, and its amplitude
says what follows:

| marker | amplitude |
|---|---|
| session start | 60 FS |
| block at 1/16, 1/32, 1/64, 1/128, 1/256 | 32, 36, 40, 44, 48 FS |
| end of session | 56 FS |

The marker in B1 comes before B2, so B2's conditioning wipes out any friction
history the marker leaves.

**B1.3 (re-zero the interferometer)** is done in analysis rather than on the
instrument. The block's zero is the mean x_IFM over the 2 s settle after the
MRES write, which is equivalent and means the session doesn't stop for an
operator action.

Set the IDS acquisition rate to at least 100 Hz. Earlier `IDSdata.txt` logs ran
at ~88 ms per sample, which gives only three samples in a 300 ms B3 window.

---

## State

**The plan and the figure are generated and checked.** Every row is continuous
with the previous one and lies on its block's microstep grid. MRES changes only
happen at full steps. Each segment starts and ends at x_start.

**Part 1 is empty** until the end switches are fitted.

**The runner is checked against its simulator only.** Every mode runs end to
end against the protocol model: upload acknowledgements, segment logs, the
x_start check after each segment, and refusal of the undefined part 1. It has
not talked to a real board.

**Nothing here can run yet.** The stage-identification sketch
(`../firmware/stage_identification`) takes targets in whole full steps, does
single-microstep staircases only, and doesn't report MSCNT. The required
protocol is in [`firmware/README.md`](firmware/README.md). `run_mres.py`
speaks that protocol, but no sketch implements it yet, so the hardware modes
stop at the firmware check.

**Before the first session:** measure σ_floor (an Allan deviation of a static
IDS record) and update `SIGMA_FLOOR_NM`. Confirm the operating I_HOLD and
TPOWERDOWN; the plan currently holds at run current with TPOWERDOWN 0.
`check_parameters()` flags any TPOWERDOWN that would end inside a B3 averaging
window.

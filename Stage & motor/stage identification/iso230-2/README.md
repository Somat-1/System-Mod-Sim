# ISO 230-2 repeatability & accuracy

The full-range baseline of positioning accuracy and repeatability (B-sequence
Step 0), measured with the IDS 3010 and evaluated in ISO 230-2 terms. The
specification is in [`PROCEDURE.md`](PROCEDURE.md). This folder implements it
as an ESP32-S3 sketch with a terminal menu.

```
iso230-2/
├── README.md
├── PROCEDURE.md                      the specification, as written
├── config/
│   └── iso_config.json               every setting, in one place
├── scripts/
│   ├── iso_design.py                 checks the config, previews the targets,
│   │                                 writes the firmware header
│   └── iso_terminal.py               control terminal (keys act at once)
├── sequence/
│   ├── iso_targets_seed<S>_rom<R>.csv    target preview for a given travel
│   └── iso_plan_summary.json         derived values, duration, checks
├── firmware/iso230_2/
│   ├── iso230_2.ino                  the runner
│   └── iso_config.h                  GENERATED from the config
└── data/                             terminal logs, one per session
```

> **Compiled, not flashed, never run.** It builds cleanly for the ESP32-S3,
> with no warnings. Bring it up the way the MRES firmware was brought up:
> switch check first, then a home with a hand near the power.

---

## Workflow

```bash
# 1. edit config/iso_config.json, then
python scripts/iso_design.py                 # checks, preview, writes iso_config.h

# 2. flash (close any terminal holding the port first)
arduino-cli compile --upload -p /dev/cu.usbmodem3101 \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" firmware/iso230_2

# 3. run
/usr/bin/python3 scripts/iso_terminal.py
```

The config is compiled in, so **any config change needs step 1 and a reflash**.
The board prints the config's hash at boot (`# BOOT,...,config=<hash>`) and in
every `# PLAN` line, so a log always says which config produced it.

## The two options

```
  1  Home + warm-up    home, place the targets, warm up (key = stop after cycle)
  2  ISO sequence      datum, n runs, end reference, park
  3  Settings          current, velocity, runs, warm-up max
  4  Show plan         targets, limits, x_ref
  5  Go to x_ref       fixed-side approach
  6  Home only         home and place the targets, then wait
  7  Driver readback   registers and status
  8  Switch state      read both end switches now
  0  Release motor     EN high, outputs off
```

**1 – Home + warm-up** runs INIT → HOME → WARMUP:

- **INIT** happens at boot. It checks that every velocity is possible at the
  configured MRES, writes the driver registers, reads them back, and refuses to
  enable the motor on any mismatch.
- **HOME** finds the home switch: a fast seek, a 2 mm back-off, then a slow
  re-approach. The slow trigger is **commanded zero**. With `measure_far_end`
  set, it then finds the far switch the same way to measure the travel, places
  the targets (next section) and parks at x_ref.
- **WARMUP** repeats the ISO run exactly, unrecorded. **Any key latches the
  stop.** The cycle in progress finishes first, so DATUM always starts from the
  same state. Press `i` as well to go straight into the ISO sequence;
  otherwise it asks when the cycle ends. At `warmup_t_max_min` it goes straight
  on by itself.

**2 – ISO sequence** runs DATUM (3 × x_ref) → `n_runs` ISO runs → END_REF
(3 × x_ref) → PARK. **Ctrl+C aborts.** Any other key prints progress and does
not stop the measurement. After an abort, a switch, soft-limit or driver fault
disables the motor. A plain Ctrl+C asks whether to park.

## How the targets are placed

They're placed **on the board, after homing, from the travel it measured**.
Nothing is hard-coded:

```
soft limits   margin (2 mm) inside each switch trigger
span          soft limits shrunk by the overrun (3 mm)  ->  5 .. ROM-5 mm
              the end is floored to 0.5 mm, so homing scatter does not move it
count m       span / p_nom (5 mm), rounded, + 1
interval p    span / (m - 1)      evenly distributed, close to p_nom
offsets       uniform +/- p/4 from a seeded generator; first and last pinned
rounding      0.001 mm, then the microstep grid of the configured MRES
x_ref         centre of the soft limits (or x_ref_mm), on a full step
```

Every planned position is then checked against the soft limits before
anything moves: targets ± overrun, and x_ref with its pre-approach.

| Travel | Targets | Span | Interval | x_ref |
|---|---|---|---|---|
| 140 mm (spec) | 27 | 5 – 135 mm | 5.000 mm | 70.00 mm |
| 122.69 mm (measured 2026-10-20) | 24 | 5 – 117.5 mm | 4.891 mm | 61.35 mm |

**Same list every session.** Same seed and the same travel (to 0.5 mm) give
the same targets, to the microstep. The "after compensation" run at the end of
the thesis then reuses them, as §4 asks. `iso_design.py` implements the same
algorithm with the same generator and the same rounding. To reproduce the
board's list offline, run it with the travel the board printed:

```bash
python scripts/iso_design.py --rom 122.6919
```

With a 2 mm lead the offsets spread the targets over the lead phase. For 140 mm
the largest gap in phase coverage is 0.12 of a lead period.

## Timing (defaults)

| | |
|---|---|
| dwell | 1.0 s: settle 0.5 s (t_hold 0 + 0.5 margin) + record 0.5 s |
| one run | 2.4 min (27 targets × 2 directions) |
| ISO sequence | 12.2 min for 5 runs, 270 measured stops |

The dwell is shorter than the procedure's "≈ 2 s" because the config holds at
run current with TPOWERDOWN 0, as the identification campaign does, so t_hold
is 0. Set `ihold_multiplier`, `tpowerdown` and `iholddelay` to the operating
values, and `iso_design.py` derives the settle time from them.

## What the board logs, and what it doesn't

The interferometer streams on its own and the board never reads it. The board
emits everything needed to cut the IDS record into stops offline: every event
and every stop's record window, with the board's millisecond clock and the
commanded position:

```
# STOP,<state>,<run>,<dir>,<i>,<x_u256>,<x_mm>,<t_arrive>,<t_rec0>,<t_rec1>,<enc_mean>
# EVENT,<t_ms>,<state>,<event>,<run>,<dir>,<i>,<x_u256>
# TARGET,<i>,<x_u256>,<x_mm>          # PLAN,...   # HOMED,...   # SWITCH,...
```

`iso_terminal.py` writes all of it, with host timestamps, to
`data/terminal_<stamp>.log`.

These parts of the specification depend on the interferometer or environment
sensors, so the board can't do them. They belong to the host or to analysis:

- the live drift display and the "READY" hint during warm-up (§3.3);
- the lost-steps and interferometer-signal watchdogs (§7);
- the environment log and the ΔT warning;
- the per-stop mean, std and slope of x_IFM (§5.2, §6.3), cut from the IDS
  record using the `# STOP` windows.

The watchdogs on the board are both limit switches, the soft limits, and the
TMC2209 fault flags (over-temperature and shorts abort; pre-warning and open
load are reported). Homing also fails if the **opposite** switch closes during
a seek, which catches swapped switch pins or a reversed DIR sense.

## Departures from the specification

- **The far switch is homed too** (`measure_far_end`), because the targets are
  distributed over the measured travel. Zero is still the home switch's slow
  trigger, homed on the same side with the same sequence every time.
- **Homing velocities default to 2.0 / 0.15 mm/s,** not 5 / 0.5. These are the
  speeds that homed reliably on this rig.
- **Each warm-up cycle starts with the run's own x_ref visit**, recorded as
  `WARMUP_REF`, so the warm-up move pattern is exactly the ISO run, rather than
  adding a second visit after it.
- **Parking after a plain abort asks first** rather than moving by itself.

## Before the first run

1. **The distance scale.** On 2026-10-20 homing measured **122.69 mm** between
   the switches at 100 full steps/mm (2 mm lead, 200 steps/rev), but the
   procedure gives the ROM as ≈ 140 mm. One of `lead_mm`, `steps_per_rev` or the
   ROM figure is wrong for this stage. Every mm in this test depends on it.
2. **MRES.** It is 16 until the microstep test chooses one. At the 5 mm/s
   approach velocity, 1/16 needs 8 kHz of step pulses. The step generator tops
   out at 20 kHz, so 1/64 allows at most 1.95 mm/s and 1/256 at most 0.49 mm/s.
   `iso_design.py` flags any velocity the configured MRES cannot reach.
3. **Hold current.** Set the operating `ihold_multiplier`, `tpowerdown` and
   `iholddelay`.
4. **Switches.** Run `8` and press each switch by hand: GPIO10 is the DIR-low
   end, and `home_side: min` homes to it.

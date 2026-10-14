# Bench changes — 2026-10-13

What was changed under `Stage & motor/` to run the jog tool and the
identification firmware on the **Stepper Isolated bench**, and what was added
along the way. Everything here has been flashed and run on that bench except
where marked.

---

## 1. Wiring: all sketches now use the bench's pin map

The sketches were written for the pin map in [`CONNECTIONS.md`](CONNECTIONS.md),
which is not how the bench is wired. On the bench, the jog firmware reported
`driver_not_ready` and would have pulsed STEP/DIR into the encoder's B/Z lines.
All four sketches now use the bench wiring, the same as `Stepper_Encoder_Slow_Test`
and `TMC2209_Communication_Check`:

| Signal | was (`CONNECTIONS.md`) | now (bench) |
|---|---|---|
| STEP | GPIO5 | **GPIO7** |
| DIR | GPIO6 | **GPIO8** |
| EN (TMC ENN, active low) | tied to GND, not driven | **GPIO9, driven** |
| TMC UART | TX 17 (via 1 kΩ) / RX 18 | unchanged |
| Encoder A / B / Z | 4 / 7 / 15 | **4 / 5 / 6** |
| R_sense | 0.11 Ω | **0.10 Ω** (module marked R100) |

EN is held HIGH (outputs off) from boot until the driver configuration has been
verified by register readback. Then it goes LOW and stays low, so the motor holds
position between moves.

Sketches changed:

- `stage identification/firmware/stage_identification/stage_identification.ino`
- `stage identification/chirp/firmware/stage_chirp_sweep/stage_chirp_sweep.ino`
- `motor identification/firmware/motor_identification/motor_identification.ino`
- `stage identification/jog/firmware/esp32s3_tmc2209_jog/esp32s3_tmc2209_jog.ino`

> **Not yet updated:** `CONNECTIONS.md`, the jog `README.md` and the stage
> identification `README.md` still describe the old pin map. This file is the
> current reference until they are.

---

## 2. Coil current: set by firmware alone, and held through every dwell

Two settings meant the actual current differed from the `CURRENT` value:

1. **`I_scale_analog` was left at its reset value of 1.** The driver scaled its
   current by the voltage on the module's VREF potentiometer. GCONF read back
   `0x1C5`.
2. **The current dropped to half at standstill.** TMCStepper sets
   `IHOLD = IRUN × 0.5` by default, and `TPOWERDOWN` was at its reset value
   (about 0.44 s). Every dwell, which is when the measurements are taken,
   therefore ran at half current after the first half-second.

Now, in the stage identification, chirp and jog firmware (motor identification
already did this):

| Register | Value |
|---|---|
| `I_scale_analog` | 0 (current from IRUN/IHOLD only) |
| `IHOLD` | = `IRUN` (`hold_multiplier(1.0)`, so it stays equal after a `CURRENT` change) |
| `IHOLDDELAY`, `TPOWERDOWN` | 0 |

The startup check now also refuses to run if `I_scale_analog` reads back as 1.
Verified on the bench: GCONF `0x1C4` (bit 0 clear), and `DRV_STATUS` reports
CS_ACTUAL = 9 at standstill, the same as IRUN.

At `CURRENT 360` with 0.10 Ω, IRUN = 9/31 with `vsense` = 1, which is about
330 mA RMS.

Drive mode, unchanged:

- SpreadCycle forced (`TPWMTHRS = 0`)
- MicroPlyer interpolation off
- `TOFF 5`, `TBL 24`, `HSTRT 5`, `HEND 0`
- all verified by readback before motion is allowed

---

## 3. Jog tool (`stage identification/jog/`)

### Firmware: `esp32s3_tmc2209_jog.ino`

| Command | What it does |
|---|---|
| `STEP <n>` | relative move of *n* microsteps at the present MRES |
| `JOGU <a> <b>` | oscillation ends in 1/256-full-step units, so ends between full steps are exact; refused if not on the present microstep grid |
| `SEQ [s] [dwell_ms] [jump] [m1,m2,…]` | MRES sweep (see below) |
| `HOLD <0..100>` | standstill current as % of run current (default 100) |
| `REINIT` | retry the driver configuration without a reset |

Other changes:

- **Position reports:** every report (`POS`, `STATUS`, `ARRIVED`, `REACHED`,
  `STOP`) also carries `pos_u256`, the exact position.
- **Driver ID at boot:** prints the driver's IOIN/version register
  (`# DRIVER_ID,…,version=0x21`), so a UART fault is visible immediately.
- **No default soft limit:** `LIMIT` is 0 (off) until you set one.

**MRES sweep (`SEQ`).** For each MRES in turn (default 1, 2, 4, 8, 16, 32, 64):

1. Set MRES.
2. **Marker jump:** out by *jump* full steps (default 4 = 40 µm), dwell, back,
   dwell. This makes each block boundary obvious in a trace.
3. **Oscillation:** +1 / −1 microstep with *dwell* ms each side (default
   2000 ms) for *s* seconds (default 30). A block always ends on a completed
   pair, so the stage is back at the start point.

A final marker jump closes the last block, and the original MRES is restored.
It must start on a full step. Defaults take about 4.3 min. Each event is printed
with the board's millisecond clock (`SEQ_BLOCK`, `SEQ_JUMP`, `SEQ_OSC`,
`SEQ_OSC_END`, `SEQ_DONE`).

### Host: `jog_stepper.py`

| At the `jog>` prompt | |
|---|---|
| `<n>` | move *n* microsteps (negative goes back) |
| `set a` / `set b` | take the present position as an oscillation end |
| `go` / `go <legs>` | oscillate between A and B (any key stops) |
| `seq [s] [dwell] [jump] [m1,m2,…]` | run the MRES sweep, logged to `runs/mres_sweep_<date>_<time>.csv` |
| `hold <pct>` | standstill current |
| `reinit` | retry driver setup |
| `limit <steps>` | optional soft limit (0 = off, the default) |

Also new:

- **`--seq [s dwell jump]`** runs the sweep from the command line.
- **Firmware check:** the script checks *which* firmware answers. If the board is
  running an identification sketch, it says so and gives the `--flash` command;
  before, it attached to it and every command came back `unknown_command`.
- **Bug fixed:** `Link.wait_for` dropped every line that arrived in the same read
  as the awaited reply. Short moves, whose `ARRIVED` lands with the
  acknowledgement, would hang.
- **Simulator** (`--simulate`) models all of the above.

Run on the bench: microstep moves, `set a/b`, `go`, `hold`, and the full default
MRES sweep (7 × 8 pairs, 254 s, log in
`jog/runs/mres_sweep_20261013_114701.csv`). That sweep was blind (no encoder
attached) and ran **before** the current fix in section 2. Its dwells were at
about half current.

---

## 4. Stage identification campaign

- **Firmware:** wiring and current as in sections 1–2. The sequence logic is
  unchanged.
- **Plan:** regenerated from `scripts/sequence_design.py`. It is identical to the
  committed `sequence/sequence_plan.csv`: 230 moves, 42.7 min, targets within
  ±10.5 mm. The acceleration-aware timing (4000 full steps/s², matching the
  firmware) was already in the working tree.
- **`scripts/run_campaign.py --blind`:** new. Skips the encoder preflight (which
  needs 200–5000 counts over 0.1 mm and so always fails with no encoder),
  records `"blind": true` in the log header, and says BLIND in the safety prompt.
  The safety prompt now names the bench pins (STEP/DIR/EN 7/8/9, encoder 4/5/6).

Run:

```bash
cd "Stage & motor/stage identification"
python3 scripts/run_campaign.py --port /dev/cu.usbmodem3101 --flash --blind
# then centre the stage and type the RUN confirmation it asks for
```

Logs go to `stage identification/recordings/`.

---

## 5. Prebuilt images

Rebuilt from the current sketches, and the manifests updated (sketch hash,
image hash, core version, pin map, current mode):

- `stage identification/firmware/prebuilt/`
- `stage identification/chirp/firmware/prebuilt/`
- `stage identification/jog/firmware/prebuilt/`

The images were built with `esp32:esp32@3.3.10`. That version pads merged
images to the full 4 MB of flash; the trailing erased bytes (0xFF) were trimmed,
which does not change what gets written to the chip.

---

## 6. Only one sketch can be on the board

The jog tool and each identification sketch replace one another. Flashing one
removes the other, which is why the jog script now checks which firmware it is
talking to. Close every serial session (`quit` at the jog prompt) before
flashing or starting another tool. Two programs with the port open each lose
some of the other's replies.

## 7. Encoder counter wrapped at ±30,000 counts (fixed)

The jog, stage identification and motor identification firmware configured the
encoder counter (PCNT) with limits of ±30,000 and `accum_count = 1`, but never
registered the limits as watch points. In this ESP-IDF version accumulation only
works with the watch points in place, so the reading silently restarted at 0
every 30,000 counts. On the 20,000-counts-per-revolution shaft encoder, that is
every 1.5 revolutions.

All three now call `pcnt_unit_add_watch_point(pcnt, ±30000)`, as
`Stepper_Encoder_Slow_Test` and the Stepper Isolated stream firmware already
did. Verified on the bench:

- 2 revolutions read 39,997 counts and returned to 1;
- the campaign's 6-revolution warm-up read 119,996 counts and returned to 0.

The stage identification run earlier today was blind, so it is unaffected.

## 8. Motor identification campaign

- **Encoder, measured at bring-up:**
  - 20,000 counts per revolution (5000 lines, quadrature)
  - positive command gives positive counts
  - one index pulse per revolution
  - 1-count return error over one revolution out and back
- **Firmware:**
  - `REINIT` retries driver setup without a reset
  - prints the driver's chip ID (`# DRIVER_ID`) so a dead UART shows as
    version `0x00` rather than as all-zero registers
  - wiring, EN and 0.10 Ω as in section 1
  - the counter fix in section 7
- **`scripts/run_campaign.py`:** zeroes the command and encoder just before
  `PLAN RUN` (as the README specifies), and writes the recording line-buffered
  so it can be followed live.
- **`scripts/monitor_run.py`** (new), a live progress window:
  - progress through the 83 actions, the current block, elapsed and remaining
    time
  - commanded vs encoder angle over the last 8 s
  - tracking error in degrees

  It only reads the recording, so it can be closed and reopened at any time
  without disturbing the run. It uses Homebrew's Python (`/opt/homebrew/bin/python3`)
  because the system Python's Tk is too old for this macOS, and the Homebrew
  matplotlib currently fails to import (a pyexpat/libexpat mismatch).

Run:

```bash
cd "Stage & motor/motor identification"
"../../Stepper Isolated/.venv/bin/python" scripts/run_campaign.py --port /dev/cu.usbmodem3101
/opt/homebrew/bin/python3 scripts/monitor_run.py      # in a second window
```

## 9. Status

| | compiled | flashed + run on bench |
|---|---|---|
| jog firmware + `jog_stepper.py` | yes | yes |
| stage identification | yes | flashed; campaign started blind |
| chirp sweep | yes | yes: 8 HIGH runs, no encoder; see `stage identification/chirp/data/2026-10-13_chirp/`. The amplitude was changed to a constant ±1.5–2.5 full steps at MRES 16/8, 400–1000 mA (see the chirp README), with an on-board generator check, a ±100-full-step pre-check and a 1 kHz turn-around marker |
| motor identification | yes | yes: encoder bring-up passed; full campaign started 15:50 |

None of these changes are committed.

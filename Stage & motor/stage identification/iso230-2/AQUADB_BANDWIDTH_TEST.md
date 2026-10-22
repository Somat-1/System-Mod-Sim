# AquadB bandwidth test — can the EL5101 keep up at 8.33 MHz?

## Why

The AquadB clock has been reduced from **400 ns to 120 ns**. At the unchanged
resolution of 1000 pm per increment that moves the metrology ceiling:

| | AquadB clock | max increment rate | ceiling |
|---|---|---|---|
| before | 400 ns | 2.500 MHz | **2.500 mm/s** |
| now | 120 ns | 8.333 MHz | **8.333 mm/s** |

The old ceiling is what truncated the v2 ISO recording: above 2.500 mm/s the
link queued increments, the recorded trace fell behind the stage, and a 10 s
hold left under 2 s of usable stationary data. Raising the clock removes that
limit *at the transmitter*. It does not prove the **receiver** can accept the
higher rate.

That is what this test establishes, and it is the only question it asks.

## The two failure modes, and how to tell them apart

They look similar in a trace and mean completely different things:

| | what happens | signature |
|---|---|---|
| **transmitter clips / queues** | increments wait their turn, none are lost | distance preserved, timing stretched |
| **receiver loses counts** | the EL5101 cannot decode the edge rate; counts vanish | **distance short, error accumulates** |

So **closure is the measurement**, not saturation. The old 400 ns behaviour was
the first kind — every transfer still closed to within 1.1 µm however hard it
saturated. If the EL5101 is out of its depth at 8.33 MHz it will be the second
kind, and the recorded distance will come up short.

## Sequence

Same structure as `OPEN_LOOP_TRAJECTORY_DIAGNOSTIC.md` — the switch baseline,
then each profile running the `2 mm → x_ref − 3 mm → x_ref → 2 mm` transfer
twice with 10 s holds, re-homing between profiles. Only the profiles changed:

| Profile | Velocity | Acceleration | Step rate at 1/16 | Increment rate | Per quadrature channel |
|---|---:|---:|---:|---:|---:|
| `v2` | 2 mm/s | 40 mm/s² | 3 200 Hz | 2.00 MHz | 0.50 MHz |
| `v5` | 5 mm/s | 40 mm/s² | 8 000 Hz | 5.00 MHz | 1.25 MHz |
| `v8` | 8 mm/s | 40 mm/s² | 12 800 Hz | 8.00 MHz | 2.00 MHz |

Acceleration is held at 40 mm/s² throughout — the earlier diagnostic already
showed acceleration has no bearing on the metrology path.

`v2` is the control: it ran clean under the old ceiling and must stay clean.
`v5` was the failing case before and should now be faithful. `v8` is the
question.

About **5.3 minutes** of motion and holds, plus the switch baseline and
re-homing.

### Two things to know before reading the result

**The 8 mm/s profile has only 4 % headroom.** 8.000 against 8.333 MHz. The
step-pulse train stalls every 20 ms for a telemetry write, the rotor rings,
and the instantaneous rate swings well beyond the mean — at 2 mm/s commanded
it reached the old 2.500 mm/s ceiling in both directions. So **expect `v8` to
touch the new ceiling even if the terminal is perfectly healthy.** Saturation
is not the failure criterion; lost counts are.

**Check the terminal's rating.** 8.33 M increments/s is 2.00 MHz per channel
after 4× quadrature decoding. Compare that against the EL5101-0011 datasheet's
maximum input frequency before trusting a pass — if the terminal is rated
below it, count loss at `v8` is expected behaviour rather than a fault, and
the usable ceiling is whatever rate it is rated for, not 8.333 mm/s.

## Running it

1. Regenerate and flash after any config edit:

   ```text
   python scripts/iso_design.py --rom 122.695
   ```

   then compile and flash `firmware/iso230_2/iso230_2.ino`.

2. Start `scripts/iso_terminal.py` **and keep the log**. Neither IsoTest.csv
   nor ISO230-2_v3.csv has one, which is why the profile boundaries in those
   had to be recovered from the trace. The `# DIAG_` records make this
   analysis direct instead of inferred.

3. Start the interferometer recording before motion, 1 kHz, and keep it
   running until `# DIAG_DONE`.

4. Let boot auto-home and park at the centre, then press `9`.

5. Save together: the raw CSV, the matching `terminal_*.log`, and the config
   hash printed in `# DIAG_BEGIN`.

## Reading the result

For each profile, from the recorded trace:

1. **Closure — the criterion.** Every transfer is the same physical distance.
   Compare them across profiles. Under the old ceiling all twelve agreed to
   **1.1 µm**. If `v8` comes up short of `v2` by more than a few µm, the
   EL5101 is dropping counts and 8 mm/s is above its usable rate. Also check
   the start-to-end closure over the whole run — lost counts accumulate,
   a queue does not.
2. **Achieved speed.** With the ceiling raised, the recording is faithful at
   `v2` and `v5` for the first time, so the move durations finally show what
   the firmware's step generator actually delivers. Compare against
   `L/v + v/a`. The earlier runs suggested moves run ~25 % long because of the
   20 ms telemetry stalls; this will confirm or refute it.
3. **Saturation.** Count samples sitting exactly on the new ceiling — 8333 nm
   per 1 ms sample. Expect some at `v8` (see above). What matters is whether
   it coincides with short closure.
4. **Deceleration.** A move that ends in far less time than its acceleration
   allows is a queue running dry, not a stage braking. 40 mm/s² needs 50 ms
   from 2 mm/s, 125 ms from 5, 200 ms from 8.

Run `analysis/aquadb_bandwidth_check.py`, which does all of this from the log
and the CSV:

```text
python analysis/aquadb_bandwidth_check.py     --log data/terminal_YYYYMMDD_HHMMSS.log     --csv ../<export>.csv --ceiling 8.333
```

It takes the profiles and per-move boundaries from the DIAG records rather
than inferring them, ties the clocks using the export's FILETIME header
against the terminal's line timestamps, and measures closure settled-to-
settled between consecutive holds.

Watch the `settled nm` column. It is the worst scatter in the 0.5 s before a
hold ended; if it is large the queue had not drained by then, so there is no
settled moment and that row's closure means nothing. Run against the old
400 ns IsoTest.csv it flags every row for exactly that reason - the backlog
there was ~12 s against a 10 s hold - and the residual closure error scales
with speed (1.06 mm at 1 mm/s rising to 3.48 mm at 5 mm/s), which is the
backlog signature rather than lost counts. On a healthy 120 ns run those
numbers should collapse to the noise floor.

`analysis/render_diagnostic_report.py` is a separate thing: it documents the
original four-profile IsoTest result and is hardcoded to those profiles.

## If it passes

Raise `v_mm_s` in `config/iso_config.json` from its present 2.0 to whatever
the test supports, regenerate, reflash, and re-run the ISO sequence. A faster
ISO run shortens the 25.8 min sequence and reduces the thermal drift that is
currently the largest systematic in the v3 dataset (+17 nm/min, common to all
targets).

Leave `settle_margin_s` at 2.0 s regardless. That margin is what kept v3 clean
— the queue always drained at least 1.14 s before a record window opened — and
it is cheap insurance at any speed.

## Status

Configured, **not yet run**. The firmware change that picks the gentlest
profile for the end-of-test park (rather than the last one listed, which is
now the fastest) has not been compiled.

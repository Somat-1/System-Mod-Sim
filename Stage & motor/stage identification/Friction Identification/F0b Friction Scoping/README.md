# F0b — friction scoping

**Status: run on 10 November 2026 and complete — 346/346 segments in
48.05 min, not aborted or time-limited.** The capture and the ESP terminal
log are in [`data/`](data/); the raw CSV is local-only, as every hardware
export in this repository is. A first look at it is in
[`../analysis/README.md`](../analysis/README.md), which puts s₉₀ near 7.5 µm,
β₀ near 2.28 µm and the measurement floor at 0.6 nm — all of them well away
from the placeholders F1 and F2 were planned on.

The [L friction-loop animation](rendered/sequence_details/f0b_L_friction_loops.gif)
uses three representative strokes and two cycles per stroke, in increasing
and decreasing order. It is an illustrative Maxwell-slip response, not
experimental data. See [model notes](rendered/sequence_details/ANIMATION.md).
Regenerate it with `python3 render_f0b_l_animation.py`. The export is limited
to below 5 MB and does not generate PNG snapshots.

The [S slow-ramp animation](rendered/sequence_details/f0b_S_slow_ramps.gif)
compares a 20 µm ramp after a reversal (S1) with one in continuation (S4). Each
frame shows a spring-block mechanism, the lag curve and the first 30
microsteps. It uses the same kind of illustrative model, with saturation
calibrated to β₀/2. Regenerate it with `python3 render_f0b_s_animation.py`.

The [D dwell probe](rendered/sequence_details/f0b_D_dwell_probe.png) and
[V speed probe](rendered/sequence_details/f0b_V_speed_probe.png) figures each
pair the commanded segments with an illustrative sketch of the expected signal.
The segments are drawn on an event axis in the S render's visual grammar.
Regenerate both with `python3 render_f0b_dv.py`.

- **D:**
  - Panel (a) shows the four trials in ABBA order: the wipe, the dwell
    (T_d), and the 5 µm reversal.
  - Panel (b) is a time ruler aligned at the reversal.
  - Panel (c) shows the expected lag after the reversal for 1 s and 60 s. The
    60 s case uses a placeholder static rise of +0.1 µm, which marks where a
    dwell effect would appear. An inset shows the creep during the dwell.
- **V:**
  - Panel (a) shows the eight passes in fast–slow–fast order.
  - Panels (b) and (c) show up and down lag curves that share one geometric
    error term. Their half-difference δ is the friction lag, and comparing δ
    between the two speeds tests for a speed effect. The δ values used are
    placeholders.

F0b is the first dedicated friction-identification session. It sizes the noise
floor, reversal ladder, slow-ramp, dwell and speed effects before the denser
F1–F3 experiments are designed.

- [`F0b_scoping_procedure.md`](F0b_scoping_procedure.md) is the bench and
  analysis procedure supplied for this experiment.
- [`f0b_config.json`](f0b_config.json) holds the values used by the sequence
  overview, including β₀ = 1.757 µm, the 58-minute gate and the 60-minute scope
  record limit.
- [`render_f0b_sequence.py`](render_f0b_sequence.py) regenerates
  [`rendered/f0b_sequence.png`](rendered/f0b_sequence.png).
- [`render_f0b_sequence_details.py`](render_f0b_sequence_details.py) is retained
  for regenerating static L/S/D/V figures if needed. Its PNG and PDF outputs
  have been removed; `rendered/sequence_details/` contains the animation and
  its model notes.

Regenerate the figure from this directory with:

```bash
python3 render_f0b_sequence.py
```


## Running the session

The board replays a host-generated table; nothing is computed on the ESP32
but step timing. `f0b_make_table.py` turns
[`f0b_sequence_config.json`](f0b_sequence_config.json) into the table, the
checks and the firmware header:

```bash
python3 f0b_make_table.py            # writes sequence/ and the firmware header
```

It reports **346 segments, 48.1 min**, 9.9 min inside the 58 min watchdog
gate and 11.9 min inside the 60 min scope record. All ten plan checks pass.
The segment count matches the 346 in `f0b_config.json` and the duration its
48.2 min, so the figure and the bench table agree.

| File | What it is |
|---|---|
| `f0b_sequence_config.json` | bench and firmware parameters (procedure sections 2 and 3) |
| `f0b_make_table.py` | builds the table, times it, checks it, writes the header |
| `sequence/f0b_command_table.csv` | the table, in the procedure's section 2.1 format |
| `sequence/f0b_plan_summary.json` | durations, counts, hashes and the checks |
| `firmware/f0b_friction_scoping/f0b_config.h` | generated; do not edit |
| `firmware/f0b_friction_scoping/f0b_friction_scoping.ino` | the sketch |

### Five segment types

| Type | Meaning |
|---|---|
| `move` | trapezoidal transit at 1 mm/s |
| `ramp` | trapezoidal profile peaking at the segment's own rate |
| `step` | n single microsteps, dwelling after each (the L ladder's 1, 2 and 5 µstep strokes) |
| `dwell` | hold in place, motor energised at hold current |
| `mark` | the recorded window of a reference stop |

`ramp` covers 3.25 Hz to 25.6 kHz in one mechanism. Below about 1 kHz the
acceleration phase is shorter than a single microstep, so a measurement ramp
is constant-rate in practice; at the top it is a real trapezoid, which is
what stops the stage being asked to start at 1 mm/s from standstill.

### Two step paths

Rates at or above **33.3 Hz** go through the buffered RMT path proven in
E1–E5. Below that one period no longer fits in an RMT symbol (30 ms maximum),
so pulses are emitted one at a time against a `micros()` deadline. This is
not a corner case: **block S runs at 6.4 Hz and the L ladder's 13 µstep
stroke at 3.25 Hz**, so the timed path carries real measurement segments.

### At the bench

1. E0 warm-up, then start the IDS recording and the scope (**≥ 60 min record**).
2. Flash `firmware/f0b_friction_scoping/`.
3. `t` to confirm the table hash matches `sequence/f0b_plan_summary.json`.
4. `g` to home, centre and run, or `p` then `r` if already homed and parked
   at 60.000 mm. Ctrl+C aborts.
5. Capture the serial output end to end — the `# SEG,...` dump after
   `SEQ_DONE` is the firmware log the analysis needs, and E5 lost four
   validity checks for want of it.

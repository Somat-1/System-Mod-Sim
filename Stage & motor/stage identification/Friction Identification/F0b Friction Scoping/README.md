# F0b — friction scoping

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


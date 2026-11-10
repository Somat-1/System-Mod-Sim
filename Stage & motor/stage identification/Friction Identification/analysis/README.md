# Preliminary check of the F0b and F1 runs

Both sessions were run on 10 November 2026 and both completed. This folder
answers two questions and stops there: **is the capture sound**, and **does
it contain what the experiments were built to measure**. It is not the
analysis of either procedure, and nothing here fits a friction model.

```sh
python3 run_data.py           # parse, align, report the alignment
python3 preliminary_check.py  # the two figures and preliminary_check.json
```

| File | Contents |
|---|---|
| `run_data.py` | loads a capture, parses the segment log, aligns the two |
| `preliminary_check.py` | the diagnostics, the two figures and the JSON |
| `preliminary_check.json` | every number quoted below |
| `rendered/preliminary_quality.png` | question 1 |
| `rendered/preliminary_physics.png` | question 2 |

The raw captures live beside each experiment, in
`F0b Friction Scoping/data/` and `F1 Reversal Sweep/data/`, with the ESP
terminal log from the same session. They are ~60 MB each and are local-only,
like every other hardware export in this repository; the npz parse caches in
`cache/` are rebuilt from them in about 10 s.

## The runs

| | F0b | F1 |
|---|---|---|
| Segments | 346 / 346 | 1119 / 1119 |
| Elapsed / planned | 48.05 / 48.06 min | 51.31 / 51.79 min |
| `complete` / `aborted` | 1 / 0 | 1 / 0 |
| Scope capture | 50.05 min, 3 003 197 samples | 53.61 min, 3 216 750 samples |
| Alignment R² | 0.999998 | 0.999999 |
| Encoder scale | 1.000122 count/nm | 1.000030 count/nm |

Neither run was aborted, time-limited or flagged unsafe, and no driver fault
appears in either log.

### The second capture was the F1 run, not an F2 run

It arrived named `F2_Experimental_Run.csv`. Correlated against the F1
command sequence it gives **R² = 0.9992**; against the F2 plan, **0.0015**.
The timing agrees as well: the F1 terminal log opens at 15:46:39, the run
takes 51.31 min, and the capture spans 15:54:03 to 16:47:39. It has been
renamed `F1_Experimental_Run.csv` and moved into the F1 folder. No F2 run
happened — F2 is still conditional on this very data.

## 1. Is the capture sound?

**Yes, and better than the procedures assumed.**

- **Quantisation is 1 nm.** The EL5101 reading the IDS gives one count per
  nanometre, 39 counts to a microstep.
- **At rest the stage is still to well under a nanometre.** Over the three
  N holds the Allan deviation bottoms out at **0.4–0.6 nm around 0.2 s**
  of averaging, against the 5 nm that `d_min` was set to as a placeholder.
  Raw 1 kHz scatter is 6–20 nm and drift is +4 nm/min at 60 and 105 mm,
  +32 nm/min at 15 mm (the motor end, as block N was designed to test).
- **The scale is right to 1.2 × 10⁻⁴.** Fitted against 96 mm of commanded
  travel, independently in both sessions.

**One thing had to be got right to make any of this true: the time
alignment.** The scope and the ESP share no clock. Correlating *positions*
is dominated by the 45 mm transits, which last 45 s, so a second of error
costs almost nothing there — that correlation lands up to 1.4 s out, and at
the ladder's 0.75 mm/s that is 1 mm of apparent lost motion. Correlating
*velocities* weights the starts and stops instead and fixes the offset to a
few ms. The check that it worked: the ladder's 3 mm rungs then come out
**0.92 µm** short of command rather than the ~220 µm a 300 ms error gives.

Two things to carry into the real analysis:

- **A ~24 ms pipeline delay** sits between command and capture. It is
  visible as a mid-move lag proportional to speed (24 µm at 1 mm/s, 0.2 µm
  at 10 µm/s) and it **cancels between two settled endpoints**, which is why
  the ladder is unaffected. Anything read mid-motion must remove it first.
- **Post-move creep inside the 0.5 s pause.** After the fastest ladder rungs
  the stage is still 2 µm short when the pause begins and catches up during
  it; after the slow rungs it is 0.3 µm. The lost motion you measure
  therefore depends on when you read it, and the 0.5 s pause is **not long
  enough above about 25 µm/s**. The 500 µm and 3 mm rungs are excluded from
  the figures for this reason — at 3 mm the catch-up cancels the lost motion
  almost exactly, which is why that rung reads ≈ 0.

## 2. Is the physics in there?

**Yes — and the two headline numbers are both about twice their
placeholders.**

| | Placeholder | Measured here |
|---|---|---|
| `d_min` | 5 nm | **0.6 nm** (0.5 s Allan deviation, worst of three positions) |
| `s90` | 3 µm | **≈ 7.5 µm** |
| `beta_0` | 1.757 µm | **≈ 2.28 µm** |

- **The ladder is resolved end to end.** Every rung from one microstep
  (39 nm) to 100 µm gives a lost motion well clear of the noise, with ℓ/s
  falling from 0.92 at one microstep to 0.09 at 50 µm. F0b's anchored
  ladder at 60 mm and F1's symmetric loops at 25, 60 and 100 mm trace the
  same curve, which is a real cross-check: two sessions, four hours apart,
  different loop geometry.
- **There is a mild position dependence.** The 25 mm loops sit about 10–15 %
  above the 60 mm ones at the same stroke, consistently across the range.
- **The four slow-ramp preparations separate exactly as the model predicts.**
  Reversals lose 2.3–2.8 µm over a 20 µm ramp, continuations 0.4–0.85 µm,
  with no overlap. Breakaway is detectable after ~60 nm of commanded travel,
  about 1.6 microsteps.
- **The dwell probe does not show ageing.** 60 s of rest changes the next
  5 µm reversal by **−63 nm** against a trial-to-trial scatter of **317 nm**.
  On two repetitions per condition the probe cannot resolve anything smaller
  than roughly 500 nm, and what it sees is consistent with zero. The limit
  is the trial scatter, not the metrology: the measurement floor is 0.6 nm.

## What this means for F1, F2 and F3

1. **F1 was run on a stroke grid built for `s90` = 3 µm, and the real value
   is about 7.5 µm.** The grid still spans the interesting region — the
   loops reach 2.8 µm of lost motion at a 15 µm stroke and the curve is
   still rising there — but its dense part sits a factor of two low, and the
   saturation plateau is only sampled by the large strokes whose 0.5 s pause
   is too short. Worth deciding whether to re-run F1 on a corrected grid
   before treating its loop set as final.
2. **F2's go/no-go reads "no" on this evidence, but weakly.** The procedure's
   test is whether the 1 s versus 60 s difference exceeds 3·`d_min`. With
   `d_min` now 0.6 nm that threshold is ~2 nm, which the probe has no hope of
   resolving against 317 nm of trial scatter. The honest reading is that the
   probe was underpowered, not that ageing is absent. Either accept the
   reduced F2 variant (14 min, bounds the effect) or repeat the dwell probe
   with more repetitions before committing to the full 43 min session.
3. **`f0b_dependent` can come off its placeholders** in `f1_config.json` and
   `f2_config.json` once the real F0b analysis confirms these numbers. Both
   `CFG_PROVISIONAL` flags clear at that point, and F2's post-dwell ramp
   (2·`s90`) grows from 6 µm to about 15 µm.
4. **The speed probe needs the trapezoidal command model** before it can say
   anything: at 1 mm/s the apparent lag is dominated by the 24 ms pipeline
   delay. At 10 µm/s, where that contributes 0.2 µm, the residual lag is of
   order β₀/2, as sliding friction should be.

All of this is preliminary. The numbers above come from settled endpoints
and simple differences, not from the fits each procedure specifies.

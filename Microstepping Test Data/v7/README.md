# v7 — settling-error campaign: the sequence design

**Question:** when the stage is commanded to a position, where does it actually
come to rest?

This folder is, so far, only the **proposed sequence and its render** — no
firmware and no measurement yet. The sequence is defined as data in
`scripts/sequence_design.py`, exported to CSV/JSON, and drawn by
`scripts/render_sequence.py`.

```
v7/
├── README.md
├── scripts/
│   ├── sequence_design.py         the plan, as data — no plotting, no hardware
│   └── render_sequence.py         -> sequence/experimental_sequence.png
└── sequence/
    ├── experimental_sequence.png  the render
    ├── experimental_sequence.svg  same, vector
    ├── sequence_plan.csv          every segment: t, x, phase, rate, measured
    ├── sequence_targets.csv       every measured position + its averaging window
    └── sequence_summary.json      every number below, machine-readable
```

Run from `scripts/`:

```powershell
python sequence_design.py
python render_sequence.py
```

## The design, and where each number comes from

The generic advice for this kind of measurement is sound, but several of its
numbers are placeholders that this rig can replace with measured ones. Every
constant below cites the measurement it came from rather than a rule of thumb.

| Choice | Value | Why this value |
|---|---|---|
| Full step | 10 µm | `FULL_STEPS_PER_MM = 100`, v6 firmware |
| Electrical cycle | 4 full steps = 40 µm | the cyclic error's period |
| Sweep span | ±18 mm (36 mm) | 44 mm measured travel, 2 mm guard each end |
| Move rate | 700 µm/s, **every measured move** | v6's MODERATE |
| Settling discard | **0.20 s** | v2.3 measured the free decay at 343.8–350.0 Hz, τ = 15–21 ms. 200 ms is ~10 time constants |
| Dwell | **1.5 s** → 1300 samples averaged | long enough to average down encoder noise, short enough that v6.1's wander does not accumulate inside one measurement |
| Back-off before re-approach | **500 µm** | v6 measured reversal errors up to **239 µm** at MRES 1; the back-off must clear the worst backlash with margin or the approach is not actually unidirectional |
| Comb spacing | ¼ full step = 2.5 µm | needs MRES ≥ 4 — see below |
| Warm-up | 20 min | v6.1 still measures +0.49 nm/s drift after 19 min of running |

### Four adaptations to this rig specifically

**1. The measurement is streamed, not triggered.** The EL5101 records
continuously at 1 kHz, so there is no indicator to trigger after settling.
Every dwell is captured in full and the position is extracted in
post-processing: discard the first 0.20 s, average the remaining 1.30 s.
`sequence_targets.csv` carries each measurement's window bounds so the
extraction is defined by the plan rather than re-derived later.

**2. The dwell has an upper bound, not just a lower one.** The usual advice is
"wait 100–500 ms for ringing to stop", which sets a floor. v6.1 adds a ceiling:
the standstill signal is dominated by slow wander, so past roughly a second a
longer dwell starts integrating drift rather than averaging noise away. 1.5 s
sits between the two.

**3. Drift witnesses are interleaved.** The sequence returns to the datum every
5 targets and between phases (23 witness visits, 3.3 min). Thermal drift over a
63-minute run would otherwise be read as position-dependent error — the
witnesses make it separable, because they are the same commanded position
measured repeatedly across the run.

**4. The comb respects what each MRES can command.** A ¼-full-step comb is
2.5 µm, which needs MRES ≥ 4; at MRES 1 the minimum commandable step *is* a
full step. Rather than emit a plan the hardware cannot execute, the comb widens
to the finest step each MRES allows:

| MRES | Comb step | Points per electrical cycle |
|---:|---:|---:|
| 1 | 10 µm | 4 |
| 4 / 16 / 32 | 2.5 µm | 16 |

## The phases

| Phase | What it isolates | Time |
|---|---|---:|
| Warm-up | thermal transient — not measured at all | 19.3 min |
| Datum | one direction-consistent origin | 0.4 min |
| Unidirectional sweep | absolute accuracy, backlash excluded by construction | 7.0 min |
| Bidirectional pass | reversal error — the one phase that deliberately lets backlash in | 11.9 min |
| Cyclic comb | the 4-full-step error, at ¼-step spacing | 13.8 min |
| Settling repeats | repeatability: 20× the identical move to one target | 7.1 min |
| Drift witness | the thermal drift that would otherwise look like error | 3.3 min |

**Budget: 62.8 min total** — 20 min warm-up plus four MRES blocks (9.0 min at
MRES 1, 11.5 min each at 4/16/32), **820 measured positions**.

## Target placement

19 global targets across ±18 mm, every one snapped to a whole full step (the
motor is accurate at full-step positions; a target between them confounds
placement error with microstep non-linearity).

They are **deliberately not evenly spaced**. A uniform grid over this span lands
on exactly 2000 µm spacing — the lead screw's pitch — so every target would
sample the same point in the screw's rotation and a once-per-revolution error
would read as a constant offset instead of a periodic one. A bounded
pseudo-random jitter (fixed seed, so the plan is reproducible) spreads the
sampled phase across the turn: the realised spread is **1960 µm of the 2000 µm
pitch, largest gap 620 µm**. Spacing ends up 790–3140 µm. The bottom-right
panel of the render shows this directly, against what a uniform grid would do.

## The render

`sequence/experimental_sequence.png` — six panels:

1. **The whole campaign** — 63 min end to end, and how much of it is warm-up.
2. **One MRES block** — the phase structure, with the three zoom windows marked.
3. **Zoom A — unidirectional**: four consecutive targets; move, arrive, hold.
4. **Zoom B — the reversal**: where the pass turns around and backlash acts.
5. **Zoom C — the ¼-full-step comb**: the 2.5 µm staircase, invisible at block
   scale.
6. **One measurement** — where the settling error is read off, with the
   ring-down and the discarded window; and **target design** — the lead-screw
   phase spreading.

The ring-down in panel 6 is a **model** drawn at v2.3's measured frequency and
damping, to show where the averaging window sits relative to the transient. Its
amplitude and the drawn error are illustrative. Everything else in the figure is
the actual commanded plan, read from the same code that exports the CSVs.

### Plot conventions

Phase colours are the dataviz reference categorical order, slots 1/2/3/5/7;
yellow (slot 4) is skipped because it is the slot that fails the all-pairs gate
beside orange, and warm-up takes neutral grey because it is the one phase that
is *not* measured. Every phase band is additionally direct-labelled or named in
the legend, so identity never rests on colour. The palette validator could not
be run here (no `node` on this machine), so the published slot order is used
unmodified rather than a re-stepped one.

## Open choices

Three things are set to defaults that are worth a decision before anything is
flashed:

- **MRES coverage.** All four (1, 4, 16, 32) are included, at 63 min. Dropping
  to MRES 16 alone is a ~32 min pilot.
- **Warm-up length.** 20 min is taken from the generic guidance; v6.1 suggests
  the stage is still drifting slightly at 19 min, so this could be worth
  measuring directly rather than assuming.
- **MRES ordering.** The campaign runs 1→4→16→32, the same order as v6 — which
  is exactly why v6 could not separate "MRES-1-specific" from "first few minutes
  of the run". Scrambling the order here would settle that, at no extra cost.

## A lesson from v6 worth carrying in

v6's analysis notes that **no serial event log accompanied its export**, so
every block boundary had to be inferred from the measurement — and one
inference bug produced a wrong result that survived a sanity check (see
[`../v6/analysis/README.md`](../v6/analysis/README.md)). This campaign has 820
measured positions whose windows must line up with the encoder stream. Logging
a timestamped event per segment, and keeping a distinctive marker move at phase
boundaries, would make that alignment something to verify rather than
reconstruct.

# E3 — global map, new seed, new positions

The E2 design, a different day, **different targets**. This is the session the
E3 procedure was written for, and the one the 30 Oct run turned out not to be:
that run re-derived E2's own seed and measured E2's own 400 positions, so it
tested the map's repeatability but never its ability to predict anywhere new.
Its analysis is in [`../e3_repeat_global_map/`](../e3_repeat_global_map/);
read its `FINDINGS.md` before this session runs, because one of its findings
should change how this one is run (see **Before you run it** below).

Status: **generated and ready to flash. No data yet.**

## What is different from E2

| | E2 | E3 |
|---|---|---|
| `seed_start` | 5200 | **9100** |
| seed accepted | 5315 | **9101** |
| config hash | `96579c689d4d` | **`b4e292911469`** |
| targets | 400 | 400 |
| grid | 64 µsteps | 64 µsteps |
| exclusions | E0 | **E0 and E2** |
| span of targets | 2.28–119.92 mm | 2.12–119.72 mm |
| estimated run time | 38.7 min | 38.6 min |

Everything else — mechanics, currents, chopper, velocity, settle and record
times, soft limits, reference protocol, the 120 ns AquadB clock — is byte-for-byte
E2's. The firmware is E2's state machine with the names and the target table
changed, so the only new code on the board is a different array of integers.

### Targets stay on the 64-µstep grid

64 µsteps is one complete motor electrical cycle: four full steps at mres = 16.
Putting every target on that grid commands every one of them at the same point
in the microstep table, so the motor's own microstep error is a **constant
offset shared by all targets** rather than a varying term aliased into the
screw harmonics. That is what keeps the geometry identifiable, and it is why
E0's residual still contains a motor-phase term that E2 and E3 cannot see. It
also keeps this session on the same electrical phase as E1/E2/E4, so all four
are directly comparable.

### The E2 exclusion is not redundant

A new seed alone is not enough. The design stratifies the travel into 400 bins
of 295 µm, and each bin holds only about seven grid points, so an independent
draw reproduces roughly one E2 target in seven — some fifty of them. Those
targets would contribute nothing to a test at new positions. `e3_config.json`
therefore excludes any grid point within ±1 grid point of an E2 target, and
the generator fails loudly if that empties a stratum rather than quietly
reopening to E2's position.

Result: **0 of 400 targets shared with E2.**

### How new the new positions actually are

| | min | median | max |
|---|---|---|---|
| distance to the nearest E2 target | 80 µm | 120 µm | 280 µm |
| screw-phase shift from that target | 14.4° | 21.6° | 50.4° |

**Every target lands in a different screw-phase bin from its nearest E2
neighbour** — 0 of 400 share one. At A₁ ≈ 1.5 µm that phase shift moves the
k = 1 prediction by a median 0.56 µm, against a residual of about 0.4 µm, so
the test has real discriminating power over the part of the model that can
overfit.

Be clear about what this does and does not test. The 295 µm stratum width is
the binding constraint — widening the exclusion to ±2 grid points empties
strata — so each target necessarily sits within ~150 µm of an E2 one in x.
That means:

- **the varying harmonics are genuinely tested**, at 400 screw phases none of
  which were in the training set; but
- **s(x) is tested at effectively the same x.** This is benign rather than a
  gap: s(x) varies on a ~20 mm scale and is fitted by two spline knots, so it
  has no fine structure to overfit, and a point 120 µm away is as honest a
  test of it as one 1 mm away. The quantity that *can* overfit is the
  coefficient surface, and that is the one being probed.

Getting large separations in x as well would mean abandoning stratification or
cutting the target count, both of which cost more in harmonic resolution than
they buy.

## Before you run it

Two things from the 30 Oct session are worth acting on first. Neither is
implemented here, because both change a working firmware you are about to
flash, and a broken state machine costs another session.

1. **Log air and stage temperature.** That session's reference moved 5.4 µm at
   −18 µm/h over its second half and nothing recorded temperature, so it could
   not even be tested against the 11.5 ppm/K of steel. The run emits
   `# ENV_SNAPSHOT_REQUEST` at the end and nothing answers it.
2. **Consider more reference visits.** Three anchors — start, turnaround, end —
   cannot distinguish a real change in s(x) from drift curvature within a
   half, and that ambiguity is what limits the 30 Oct result. The procedure
   anticipates a D3 variant for exactly this. The cheap version is extra
   visits at the turnaround, where the stage is already near the reference;
   mid-pass visits are not cheap, because the reference sits at mid-travel and
   reaching it means a reversal that would disturb the backlash state of every
   target after it. Say the word and I will add a `CFG_TURN_REF_VISITS` and
   regenerate — it is a small, contained firmware change, but it is a firmware
   change.

If the session drifts the way 30 Oct did, ȳ and the harmonics will still be
sound — the mirrored pass design protects them — but β₀ and the gain will
again be unusable.

## Running it

```sh
python e3_sequence.py --check-only     # regenerate nothing, just validate
python e3_sequence.py                  # writes sequence/ and the firmware header
```

Then flash `firmware/e3_global_map/` and drive it with:

```sh
python e3_terminal.py                  # logs to data/e3_terminal_<stamp>.log
```

The board requires the stage homed, warm, at `x_ref` and the interferometer
already recording before it will arm. Keys: `p` accept precondition and arm,
`g` home, centre and run without operator gates, `r` run, `t` target summary,
`s` status, `0` release.

Save the interferometer export as `data/E3_Experimental_Run.csv` and leave the
terminal log beside it. Then:

```sh
python e3_run.py                       # aligns the capture to the plan
```

which should report `matched 807/807 planned stops`.

## Layout

| Path | Contents |
|---|---|
| `e3_config.json` | the experiment definition; hashed byte-for-byte into the firmware header, so it is LF-only by `.gitattributes` |
| `e3_sequence.py` | generates and validates the targets, the plan and the header |
| `e3_terminal.py` | serial driver and logger |
| `e3_run.py` | capture → holds → planned stops, all of it imported from `e2_run` |
| `sequence/` | the target list, the full action plan and the validated summary |
| `firmware/e3_global_map/` | the sketch and the generated `e3_config.h` |
| `data/` | put the capture and the terminal log here; the alignment cache is gitignored |

## Validation

All thirteen checks in `sequence/e3_plan_summary.json` pass:

```
target_count                              unique_and_increasing
targets_on_64_ustep_grid                  phase_bins_at_least_five
no_e0_target_within_one_grid_point        no_e2_target_within_one_grid_point
no_target_shared_with_e2                  seed_differs_from_e2
all_motion_inside_proven_soft_limits      two_directional_stops_per_target
seven_reference_stops                     estimated_recording_below_60_min
aquadb_velocity_limit_respected
```

Phase-bin occupancy runs 5–15 against a minimum of 5. Target span 2.12–119.72 mm,
spacing 0.04–0.56 mm. 800 directional stops plus 7 reference stops = 807, the
same structure E2 produced, so the existing alignment and analysis code applies
unchanged.

# Microstepping Test Data v3

Reliable-identification campaign on the EVO dedicated controller: six
MRES/current configurations, each running the reference, conditioning,
Condition C (creep/settling), and Condition D (velocity-plateau sweep)
blocks documented in `../v2/docs/HARDWARE_RUNNERS.md`. This revision holds
the actual hardware run and its processing pipeline, not just the command
specification.

## Layout

- `data/raw_local/SteppingSequenceID.csv`: raw Dewesoft IDS linear-encoder
  export (nm/count, tab-separated). Local-only, .gitignored (>50 MB) --
  it is a direct hardware capture, not regenerable, so keep a backup of it
  outside this repo.
- `data/raw_local/identification_controller_log.csv`: the EVO controller's
  own event log (one row per command/event), the authoritative source of
  block boundaries and configuration metadata.
- `data/splice_index.csv`: master index of every spliced block (see
  "Segmentation" below).
- `data/reliable_identification_summary.json`: campaign-level metrics
  (duration, measurement polarity, per-run peak/residual figures).
- `data/processed_local/`: one npz per run (`run_0N_mres_M_i_XXXpct.npz`)
  plus `campaign_timeseries.npz` for the full aligned record.
- `data/splices_local/`: one npz per spliced block (156 total), named
  `run_<NN>_<sequence>_<block>.npz`.
- `scripts/process_reliable_identification.py`: parses the two raw files,
  aligns them, splices every block, and renders the montages below.
  Regenerate with `python scripts/process_reliable_identification.py`.
- `rendered_assets/`: `campaign_overview.png`, `configuration_montage.png`,
  `creep_c_montage.png`, `reference_repeatability_montage.png`, and
  `velocity_montages/run_0N_*.png`.

## What was run

Six configurations, run back to back in one ~39.6-minute session:

| Run | MRES | Current |
|---|---:|---:|
| 1 | 1/4 | 200 mA (50%) |
| 2 | 1/4 | 400 mA (100%) |
| 3 | 1/2 | 200 mA (50%) |
| 4 | 1/2 | 400 mA (100%) |
| 5 | 1/1 | 200 mA (50%) |
| 6 | 1/1 | 400 mA (100%) |

Each configuration runs the same block sequence:

```
BLOCK_0_START -> conditioning (before C) -> C -> conditioning (before D)
  -> D_0.125 -> D_0.375 -> D_1.25 -> D_3.5 -> D_9.5 -> D_27.5 -> D_70
  -> D_200 -> BLOCK_0_END
```

- **BLOCK_0_START / BLOCK_0_END**: an identical small reference move
  sequence run at the start and end of each configuration, to bound
  repeatability/drift over the run.
- **C**: a directional creep/settling probe -- small approach moves from
  both directions followed by a 60 s dwell, to expose any slow
  post-motion relaxation (presliding creep) and direction-dependent
  asymmetry.
- **D_<rate>**: eight constant-velocity plateaus from 0.125 to 200 full
  steps/s (both directions), a Stribeck-curve sweep from quasi-static to
  Coulomb/viscous-dominated motion.

## Segmentation: separating blocks out of the raw encoder trace

The raw IDS file (`SteppingSequenceID.csv`) carries no event labels of its
own -- just a timestamp column and a running encoder count. Two
independent mechanisms let you recover block boundaries from it:

1. **Controller log timestamps (what the processing script actually
   uses).** `identification_controller_log.csv` logs a UTC timestamp on
   every `BLOCK_START`/`BLOCK_END`/`RUN_CONFIG`/`RUN_COMPLETE` event.
   `process_reliable_identification.py` converts the IDS file's own
   `Starttime of export` (Windows FILETIME) and `SampleTime[ms]` header
   fields into the same epoch-time base, so every controller-log
   timestamp maps onto an exact sample index in the IDS record. Each
   `BLOCK_START`/`BLOCK_END` pair becomes one row in `splice_index.csv`
   and one file in `data/splices_local/`.

2. **Data-visible marker signatures (a hardware-level fallback that needs
   no log at all).** Immediately before every block, the controller fires
   a distinct, amplitude-coded blip: a rapid move of N full steps, a 1.0 s
   dwell, an equal move back, and a 0.5 s settle. The amplitude `N`
   identifies which block follows, so the blocks are identifiable directly
   from the position trace by eye or by amplitude-matching, independent of
   the log:

   | Marker amplitude (full steps) | Identifies |
   |---:|---|
   | 68, 72, 76, 80, 84, 88 | start of configuration (run 1-6) |
   | 12 | conditioning before C |
   | 16 | conditioning before D |
   | 20 | C |
   | 24 | D 0.125 |
   | 28 | D 0.375 |
   | 32 | D 1.25 |
   | 36 | D 3.5 |
   | 40 | D 9.5 |
   | 44 | D 27.5 |
   | 48 | D 70 |
   | 52 | D 200 |
   | 56 | final Block 0 reference |

   These markers are spliced out into their own `MARKER_*`/`COND_BEFORE_*`
   npz files alongside the real blocks (e.g.
   `run_01_005_marker_c.npz` immediately precedes `run_01_006_c.npz`) --
   do not mistake a marker for the block it announces. `find_block()` in
   the processing script anchors its regex to the full block name
   (`^C$`, `^D_9\.5$`, ...) specifically to avoid matching a `MARKER_*` or
   `COND_BEFORE_*` block instead of the real one; an earlier unanchored
   version of that pattern did make exactly that mistake for the creep and
   reference-repeatability plots.

`splice_index.csv` columns: `run_index, current, mres, block, start_s,
end_s, sequence, start_sample, end_sample_exclusive, duration_s, samples,
position_range_um, end_residual_um, local_npz`. `sequence` is the
1-based order of the block across the whole campaign (1-156); `block` is
the label from the controller log (the real block name or a
`MARKER_*`/`COND_BEFORE_*` separator).

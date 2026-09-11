# v4 analysis — EL5101 IDS captures

Rendering for the two encoder captures in `../data/hardware_runs/`. One
subfolder per run, plus shared loading and styling in `scripts/ids_common.py`.

```
analysis/
├── scripts/
│   ├── ids_common.py              loader, style, decimation, marker detection
│   ├── plot_tmc_stepping.py       -> tmc_stepping/
│   ├── plot_tmc2209_exp_run.py    -> tmc2209_exp_run/
│   └── analyse_ramps.py           -> tmc2209_exp_run/ (ramp decomposition)
├── tmc_stepping/                  MRES oscillation sweep
└── tmc2209_exp_run/               full 41 min campaign
```

Run either script from `analysis/scripts/`:

```powershell
python plot_tmc_stepping.py
python plot_tmc2209_exp_run.py
```

## Data format

Both exports are Beckhoff EL5101-0011 encoder traces from a TWinCAT "YT Scope
Project": a tab-separated header block, then `index<TAB>counter` at the rate in
`SampleTime[ms]` (1 kHz for both). The counter is a raw UINT32; the
repository's established scale is **1 nm per count**, so one full step of the
stage is 10 µm. Measured encoder noise floor is **3.2 nm** (1σ over the
quietest 1 s window).

---

## `tmc_stepping/` — MRES oscillation sweep

`TMCStepping.csv`, 187.25 s, recorded 18:09:30–18:12:37. The
`tmp_v4_oscillation_mres_sweep` firmware: **SpreadCycle, MicroPlyer disabled**,
MRES 1-2-4-8-16-32, 20 s per size, 2 s dwell at each end of the cycle.

| Plot | What it shows |
|---|---|
| `overview_full_capture.png` | the complete raw capture with marker onsets |
| `commanded_overlay.png` | measured against the exact commanded sequence |
| `oscillation_by_mres.png` | six small multiples, one per resolution |
| `step_size_vs_mres.png` | commanded vs realised step, and realised fraction |

The commanded trace is reconstructed from the firmware constants, so this is a
genuine comparison rather than an annotation.

### Alignment

Each block is anchored on its **marker return edge** — the sharp crossing back
up from −200 µm — plus the firmware's 0.5 s settle. An earlier version anchored
on the marker *onset* plus a nominal 1.7667 s covering two rate-dependent moves
and two dwells; that carried a systematic **−67 ms** bias, consistent across all
six blocks, and every downstream window inherited it. The return edge is only
0.5 s of nominal time away from the oscillation, so far less can go wrong.

The cycle period comes from the measured block spacing — 22.078 s, giving
4.062 s per cycle rather than the nominal 4.000 s, because the firmware spends
~0.3 s per block on UART register reads and logging.

**The first cycle is excluded** from the statistics. It sits on the stage's
settling transient from the marker's 200 µm move, which at fine MRES is larger
than the commanded step itself. It is still drawn, shaded, so the exclusion is
visible rather than silent.

### Result

Four cycles per resolution (first excluded), detrended against the return
plateaus:

| MRES | Commanded (µm) | Measured (µm) | Realised | t | Resolved? |
|---:|---:|---:|---:|---:|---|
| 1 | 10.0000 | 10.722 ± 0.049 | 1.07 | 221.0 | yes |
| 2 | 5.0000 | 3.471 ± 0.048 | 0.69 | 72.3 | yes |
| 4 | 2.5000 | 0.547 ± 0.017 | 0.22 | 32.3 | yes |
| 8 | 1.2500 | 0.117 ± 0.023 | 0.09 | 5.0 | **no** |
| 16 | 0.6250 | −0.000 ± 0.048 | −0.00 | −0.0 | **no** |
| 32 | 0.3125 | 0.023 ± 0.006 | 0.07 | 4.1 | **no** |

**Finer microstepping buys almost nothing on this stage.** The realised
fraction falls monotonically from 1.07 at full-step to 0.22 by MRES 4, and at
MRES 8 and beyond the response cannot be distinguished from zero.

Significance is a two-sided t-test at 99% over the retained cycles, not a fixed
multiple of the standard error. That distinction matters here: excluding the
first cycle leaves n=4, so the critical value is 5.841, and a "3 standard
errors" rule would have reported MRES 8 and MRES 32 as detections. They are
not — note that MRES 32 reads *larger* than MRES 16 while its commanded step is
half the size, which is the signature of scatter rather than signal.

The limit is repeatability, not encoder resolution: the noise floor is 3.2 nm,
two orders of magnitude below the 312.5 nm MRES-32 step.

Two further methodological points behind those numbers:

- Each block is **re-zeroed on its own pre-oscillation baseline**. The axis
  does not return to the same absolute position between blocks, so a global
  zero pushes the fine-MRES traces off-scale and invites a false reading.
- The per-cycle step is **detrended against the return plateaus**. Every block
  shows monotonic creep, and by MRES 8 that creep exceeds the commanded step;
  an un-detrended difference would absorb it and overstate the result.

---

## `tmc2209_exp_run/` — full v4 campaign

`TMC2209ExpRun.csv`, 2574.7 s (42.9 min), recorded 17:20:06–18:03:00. The
complete campaign: 4 MRES values × [oscillation + 6 × 25 mm out-and-back].

> ⚠️ **Provenance.** This run predates the SpreadCycle fix. MicroPlyer was
> disabled but the chopper was left at its power-on default, so it ran under
> **StealthChop**. See `TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md`.

> No serial event log exists for this execution — it was started from the
> board's own `RUN` command, not a host capture — so every block boundary is
> **inferred from the measurement** and labelled as such.

| Plot | What it shows |
|---|---|
| `overview_full_capture.png` | the complete raw capture |
| `trajectories_overview.png` | the 24 detected 25 mm legs highlighted |
| `trajectory_montage.png` | each leg individually |
| `oscillation_blocks.png` | the inferred one-microstep oscillation windows |
| `travel_and_drift.png` | per-leg travel and cumulative origin drift |
| `ramp_asymmetry_and_settling.png` | outbound-vs-return asymmetry and endpoint creep |

All 24 trajectories are detected. Detection thresholds against a low
percentile (the resting level) rather than the median, which on a triangular
profile sits mid-ramp.

### Result

**Travel is held; the origin is not.**

| | |
|---|---|
| Travel per leg | 24.94 mm mean (commanded 25.0 mm) |
| Cumulative origin drift | **−2.47 mm over 24 legs** |

The axis moves its full commanded distance essentially every time, but the
resting position marches downward by 2.47 mm across the campaign. The drift is
not gradual:

- **Legs 1–6 (MRES 1)** — resting level flat at 0.000 mm, no loss at all.
- **Leg 7** — the first MRES 4 leg. Travel is 23.15 mm, **1.8 mm short**, and
  the origin drops to −1.79 mm.
- **Legs 8–17** — progressively smaller losses, settling at −2.47 mm.
- **Legs 18–24** — flat again; no further loss.

Almost the entire error appears at the **MRES 1 → MRES 4 transition** and
during the MRES 4 block, then stops. The firmware's own `checkOrigin()` could
not see any of this: it tracks commanded position, which returns to zero by
construction.

### Why it drifts: the outbound strokes fall short

Splitting each leg into outbound and return travel settles the mechanism. The
returns deliver their full ~24.94 mm every time; the **outbound** strokes are
the ones that come up short, and the shortfall is never recovered.

| Leg | MRES | Mode | Rate | Out (mm) | Back (mm) | Out − Back |
|---:|---:|---|---|---:|---:|---:|
| 1–6 | 1 | both | all | 24.947–24.952 | 24.944–24.953 | −2 to +6 µm |
| **7** | 4 | DIRECT | slow | **23.148** | 24.943 | **−1795 µm** |
| 8 | 4 | DIRECT | moderate | 24.784 | 24.946 | −161 µm |
| 9 | 4 | DIRECT | fast | 24.614 | 24.850 | −236 µm |
| 10 | 4 | INDIVIDUAL | slow | 24.861 | 24.940 | −79 µm |
| 12 | 4 | INDIVIDUAL | fast | 24.822 | 24.940 | −118 µm |
| 16 | 16 | INDIVIDUAL | slow | 24.865 | 24.945 | −80 µm |
| others | 16, 32 | both | all | ≈24.94 | ≈24.94 | ≈0 |

Six legs account for the entire −2469 µm, and leg 7 alone for 1795 µm of it.
Leg 7 is the first leg after the MRES 1→4 transition. MRES 1, 16 and 32 are
clean; the loss is confined to the MRES 4 block plus one straggler.

### Settling error at the endpoints

Creep during the 1 s dwell, measured from the moment motion stops:

| Endpoint | Mean | Mean magnitude | Worst |
|---|---:|---:|---:|
| Top (25 mm) | **+0.53 µm** | 0.80 µm | 1.66 µm |
| Bottom (origin) | **−0.43 µm** | 0.71 µm | 1.12 µm |

Both sub-micron and **systematically signed**: the top creeps forward, the
bottom backward — relaxation continuing in the direction of the preceding
motion, the presliding-friction recovery signature relevant to the Rev 4 LuGre
work. Magnitude grows with MRES (~0.6 µm at MRES 1, ~0.96 µm at MRES 32), and
~0.8 µm of endpoint creep exceeds a commanded microstep at MRES 16 (0.625 µm)
and MRES 32 (0.3125 µm).

### DIRECT vs INDIVIDUAL step generation

| Metric | DIRECT | INDIVIDUAL |
|---|---:|---:|
| Mean travel error | 245.8 µm | **81.0 µm** |
| Mean \|net drift\| per leg | 183.6 µm | **24.3 µm** |
| Top settling | 0.832 µm | 0.777 µm |

**INDIVIDUAL is 3× more accurate in travel and 7.5× better on drift.** Settling
ties, as expected — that is mechanical relaxation after motion ends and does not
depend on how the pulses were generated.

But the advantage is not general. At MRES 1 and 32 the two modes are
indistinguishable (±5 µm); the entire gap comes from the MRES 4 block, where
DIRECT's continuous pulse train lost 1716 µm more than INDIVIDUAL on the slow
leg. At MRES 16 slow the sign even reverses, DIRECT winning by 78 µm. The
defensible statement is that DIRECT was more vulnerable to whatever went wrong
at MRES 4, not that INDIVIDUAL generates better pulses in general.

> **Before treating the MRES 4 collapse as a property of the mechanism, re-run
> this campaign under SpreadCycle.** Every number in this section comes from a
> StealthChop run. The stepping sweep — the only SpreadCycle dataset here —
> showed clean, repeatable behaviour with no comparable losses.

The resting level is measured from genuine dwells — samples whose local slope
is below 20 µm/s — because sampling a fixed interval before the threshold
crossing lands mid-ramp and reports the leg's speed instead of the axis
position.

---

## Plot conventions

Series colours are `#1769aa` (measured) and `#d94801` (commanded), checked in
OKLab rather than by eye: ΔE 32 for normal vision, 24 under the worst
dichromatic simulation, contrast ≥3:1 on both light and dark surfaces.
Identity never rests on colour alone — the commanded trace is dashed and both
series are named in the legend.

Both traces are positions in the same unit, so they share **one** axis; there
are no dual-scale plots here. Six resolutions are shown as small multiples
rather than six overlaid series. Long traces are decimated with a min/max
envelope so extremes survive rather than being dropped between samples.

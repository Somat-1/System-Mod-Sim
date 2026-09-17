# v6 analysis — EL5101 IDS capture

Rendering for `../EXPrun.csv`, the full v6 campaign capture. One shared
loader/style module, one rendering script, one output folder — same layout as
`../../v4/analysis/`.

```
analysis/
├── scripts/
│   ├── ids_common.py       loader, style, decimation, marker/ramp detection
│   └── plot_v6_exp_run.py  -> exp_run/
└── exp_run/                 full ~21.7 min campaign
```

Run from `analysis/scripts/`:

```powershell
python plot_v6_exp_run.py
```

## Data format

Beckhoff EL5101-0011 encoder trace from a TwinCAT "YT Scope Project": a
tab-separated header block, then `index<TAB>counter` at 1 kHz. The counter is
a raw UINT32; established scale is **1 nm per count**, so one full step of the
stage is 10 µm.

---

## `exp_run/` — full v6 campaign

`EXPrun.csv`, 1302.3 s (21.7 min), recorded 13:48:58–14:10:40 on 2026-09-15.
The complete campaign: 4 MRES values (1, 4, 16, 32) × [oscillation + 6 × 10 mm
out-and-back]. v6 is v4's campaign firmware forked with two corrections —
**R_SENSE = 0.11 Ω** (v4 used 0.03 Ω, a different board's value, which
under-delivered current ~3×) and a **10 mm leg** instead of v4's 25 mm (the
stage's real travel is ~44 mm total; a 25 mm one-way leg does not fit). Unlike
the v4 TMC2209ExpRun capture, `configureDriver()` here refuses to start unless
GCONF/CHOPCONF/DRV_STATUS all confirm SpreadCycle active and interpolation
off, so this run is not suspected of the StealthChop problem that affects v4's
pre-fix data. See `../README.md`.

No serial event log accompanies this export, so — as with v4 — every block
boundary is **inferred from the measurement** and labelled as such.

| Plot | What it shows |
|---|---|
| `overview_full_capture.png` | the complete raw capture |
| `trajectories_overview.png` | the 24 detected 10 mm legs highlighted |
| `trajectory_montage.png` | each leg individually, spliced 6-per-row by MRES group |
| `oscillation_blocks.png` | the inferred one-microstep oscillation windows |
| `travel_and_drift.png` | per-leg travel (outbound vs return) and cumulative origin drift |
| `rest_window_diagnostic.png` | why an earlier version of the above was wrong (see below) |
| `marker_return_residual.png` | does each leg's marker return to exactly where it started? |

All 24 trajectories are detected (4 MRES × 6 = D/I × slow/moderate/fast).

### Result

**Travel is accurate and repeatable across the whole campaign, not just the
finer MRES values.**

| | |
|---|---|
| Travel per leg | 9.9935 mm mean, **0.0014 mm sd** (commanded 10.0 mm) |
| Cumulative origin drift | **+26 µm over 24 legs** (+1.1 µm/leg, roughly linear) |

This is a much cleaner result than v4's pre-fix StealthChop run (which lost
2.47 mm of travel concentrated at the MRES 1→4 transition): with SpreadCycle
confirmed active and the correct sense resistor, travel tracks the 10 mm
command to within a few microns on every leg, and drift is small and gradual
rather than concentrated in one collapse.

### A correction: the "MRES 1 only" resting-level scatter was a bug, not a finding

An earlier version of this analysis reported the MRES 1 group's resting level
swinging by up to ~240 µm leg-to-leg (vs ~2 µm for MRES ≥ 4) and treated it as
real settling behaviour, having checked that widening/narrowing leg 1's lookback
window didn't change its answer. That check was insufficient — it only tested
whether the window's *size* mattered, not whether the underlying *method* was
sound. The method was: median every "still" (low-slope) sample across the
**whole gap** between one leg and the next. That gap is not a plain dwell — it
also contains the marker that precedes every leg (and, at the start of each
MRES group, the full 30 s oscillation block) — so the median could land on
whichever cluster of still samples happened to be larger: the genuine pre-ramp
rest, or an unrelated dwell/settling transient elsewhere in the gap. Which one
won varied leg to leg close to at random.

The fix, now in `plot_v6_exp_run.py`, computes each leg's resting level
separately from its own `ramp_extent` window (bounded to where motion has
genuinely stopped on *either side of that leg's own ramp*), never reaching into
an unrelated gap. Applied to all 24 legs, this drops travel's sd from
**0.067 mm to 0.0014 mm** — the scatter, including the apparent MRES 1
peculiarity, is gone. `rest_window_diagnostic.png` renders both the old and new
sampling windows on the same raw signal for the six MRES 1 legs: three (#3, #4,
#6) show ~0 µm difference because their gaps are short and clean, and three
(#1, #2, #5) show the old method's median visibly pulled onto a transient dip
elsewhere in the gap — 239 µm off, in leg #5's case.

**Takeaway:** don't trust a per-leg statistic computed from a window defined by
"time between two detected events" without checking what else can fall inside
that window. Bounding by *your own event's* local stillness, not the *previous
event's* endpoint, is the safer default here.

### But there is a real MRES-1-specific effect, just not in leg travel

Fixing the bug above also made it possible to measure something the bug had
been masking: whether each leg's **marker** — the small out-and-back move that
precedes every leg, used to segment the capture — actually returns the axis to
where it started. This is `rest_after[leg i]` compared against
`rest_before[leg i+1]`, both genuine ramp-local dwell measurements, so the gap
between them is entirely attributable to what happened inside that one marker.

`marker_return_residual.png` shows it cleanly: **the 5 marker returns within
the MRES 1 group land 40–239 µm away from where they started** (mean magnitude
~120 µm), while **all 17 marker returns from MRES 4 onward land within 0.3–2.2
µm** (sd 0.47 µm). Traced sample-by-sample for one case (leg 4 → leg 5, the
+239 µm outlier): the marker moves from −279.07 µm to −358.76 µm, holds ~1 s,
returns — and settles at −39.90 µm, not back at −279.07 µm. The commanded
move net cancels (the firmware's own `checkOrigin()` would have aborted the
run otherwise), so this is a genuine **measured** reversal error, not a
logging artefact.

This is confined to MRES 1 and to *reversal* specifically — each leg's own
outbound+return travel is accurate (see above) even within the MRES 1 group;
it's only the small reversal inside the marker that lands somewhere else. That
pattern — large reversal error at the coarsest step size, vanishing at every
finer one — is a textbook backlash / detent-torque signature, worth chasing
for the LuGre/GMS friction work. One caveat: MRES 1 is also chronologically
first in this campaign (the run always executes 1→4→16→32), so this dataset
alone cannot separate "MRES-1-specific" from "first few minutes of this
particular run" — a repeat with the MRES order scrambled would settle it.

### Oscillation step size (peak-to-peak, inferred window)

| MRES | Peak-to-peak measured [µm] |
|---:|---:|
| 1 | 14.48 |
| 4 | 6.94 |
| 16 | 0.87 |
| 32 | 0.65 |

These are inferred-window peak-to-peak readings (a quick per-block sanity
check), not the paired plateau-to-plateau statistics the v4/v5 `tmc_stepping`
and `microstepping_only` analyses compute — this script's job is the
trajectory campaign, not a step-fidelity study. For that, see
`../../v5/analysis/microstepping_only/` and `../../v4/analysis/tmc_stepping/`.

---

## Plot conventions

Same as v4: series colours `#1769aa` (measured) / `#d94801` (commanded),
identity never resting on colour alone (commanded is also dashed, both series
named in the legend). Long traces are decimated with a min/max envelope.

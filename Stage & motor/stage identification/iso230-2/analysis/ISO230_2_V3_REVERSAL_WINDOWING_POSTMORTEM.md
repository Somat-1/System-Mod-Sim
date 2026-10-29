# ISO 230-2 v3 reversal post-mortem: a record window with no guard band

## Executive conclusion

The reversal value published for the ISO 230-2 v3 run, `B_mean = -0.365 µm`, is
wrong. The correct value from the same capture is `-1.601 µm`, and the error was
introduced entirely in evaluation, not in the measurement.

`iso230_2_evaluate.py` maps the board's millisecond clock onto the scope clock by
aligning logged arrival instants onto ends of motion detected in the trace. The
motion-end detector smooths with a **centred** 51-sample boxcar, which cannot
report the end of motion until about 25 samples after it happens. Measured
against true rest, the detector is **+34 ms late**. The fit therefore places the
whole board timeline about **38 ms late**, and because the firmware starts the
next move at the same instant the record window closes, every one of the 240
record windows ends 38 ms *after* the stage has started moving again.

The next move always continues the pass in progress. The contamination is
therefore signed with the approach direction, it does not average out, and it
lands squarely in `B = x̄↑ - x̄↓`.

Nothing is wrong with the stage, the firmware, the interferometer or the
EL5101 link. The same defect does not affect E1, E2 or the ISO v4 run, because
those anchor the averaging window inside the measured hold.

## Symptom

E1 (October 2026) reported a mean reversal of about `-1.68 µm` across its three
windows, flat with position and tight (sd 0.3-0.5 µm). Step A, derived from the
ISO 230-2 v3 run with the same stage, same configuration and the same commanded
settle and record times, had reported `-0.365 µm` with `max|B_i| = 0.982 µm`.

A factor of 4.6 between two runs of one machine is not a tolerance question. One
of the two numbers had to be wrong before either could be quoted.

## How it was found

The investigation ran in four stages. Only the last one is conclusive; the first
three matter because they close off the alternatives.

**1. Eliminate the physical explanations.** The two runs differ in stop spacing
(0.12 mm with 101 stops between reversals in E1, 4.891 mm with 24 in ISO v3) and
in firmware generation (ISO v3 predates the RMT step generator and cruised at
1.581 mm/s, 79.1 % of command; E1 uses RMT at 2.001 mm/s). Each was checked and
rejected as a cause:

- *Take-up with approach length* - E1's `B_i` is flat with position, slope
  -1.09 nm/mm, and the ISO grid is 40x coarser, so a take-up effect would have
  to run the wrong way.
- *Drift between passes* - the reference visits bound it well below the
  discrepancy.
- *Motor electrical phase* - `fig2_electrical_cycle.csv` gives a
  forward-minus-reverse difference that is flat at `-0.491 ± 0.027 µm` across
  all 64 microstep phases, so electrical phase cannot produce a per-target
  difference of this size.
- *Per-pass datum correction* - removing it changes nothing in `B`, which is a
  difference of two means at one target and so is blind to any constant.

**2. Re-derive E1 without the board clock.** E1's own loader fits the ESP clock
rate and offset and takes its record windows from the firmware's logged
`record_start`/`record_end`. That is the one step a sceptic should want removed.
`check_e1_no_board_clock.py` finds the dwells in the trace, matches them to the
planned stop order by position alone, and averages inside the measured hold. It
returns `-1.681 µm` against the board-clock pipeline's `-1.680 µm`, and six
different averaging windows span 0.031 µm. E1's number is not a processing
artefact.

**3. Put the ISO v3 capture through that same code.** `check_iso_control.py`
reads the v3 log only for the ordered list of (target, direction) - plan
information, no timestamps - and otherwise uses the identical path. It returns
`-1.601 µm`, not `-0.365 µm`. This is the decisive step: a pipeline that
manufactured a large reversal would have manufactured one only for E1. Instead
it finds the same reversal in the data behind the published number.

**4. Locate the defect.** With measured hold boundaries available, the published
window can be compared directly against the hold it belongs to:

| quantity | median | 5-95 % |
|---|---:|---:|
| window start minus measured rest | +2030 ms | +2012 .. +2037 |
| window end minus measured hold end | **+38 ms** | +38 .. +39 |
| windows overrunning the hold | **240 of 240** | - |

A one-millisecond spread over 240 windows is a systematic offset, not jitter.

## Root cause

```python
# iso230_2_evaluate.py, fit_clock()
sp = np.convolve(np.abs(d), np.ones(51) / 51, mode="same") * 1000
moving = sp > 0.02
ends = np.where(moving[:-1] & ~moving[1:])[0] + 1
```

A centred boxcar of width 51 still contains moving samples for 25 samples after
the stage stops, so `ends` is late. Measured against the true rest instants:

| motion-end smoother | detected end minus true rest |
|---|---:|
| 51-sample (as published) | **+34.0 ms** |
| 5-sample | +0.0 ms |

`fit_clock` fits the offset so that mapped arrivals sit on these late `ends`, so
the offset inherits the lag. The firmware's own timings are exact to the
millisecond (`t0 - t_arr` = 2000 ms, `t1 - t0` = 500 ms), so the error passes
straight through to the window.

That the whole 38 ms is one offset error, rather than two effects, can be
checked by moving the offset back by 38 ms:

| | published offset | offset 38 ms earlier |
|---|---:|---:|
| `t_arr` minus measured rest | +29.6 ms | **-8.4 ms** |
| `t1` minus measured hold end | +38.1 ms | **+0.1 ms** |

At the corrected offset the window closes exactly on the measured hold end, and
the board's arrival instant falls 8.4 ms *before* physical rest - which is what
it should do, since the firmware logs arrival when the last step pulse is issued
while the stage is still in its deceleration tail. The measured dwell is a
consistent 2492 ms rest-to-rest either way.

## Why the error went into B and almost nowhere else

In the 38 ms of contaminated window the stage travels **22.9 µm**, signed with
the pass direction: `+22.8 µm` on up stops, `-22.9 µm` on down stops. Averaged
over the 500-sample window that drags each direction mean by `∓0.617 µm`.

`B_i = x̄↑ - x̄↓` subtracts two quantities biased in opposite senses, so the
biases add: `-1.235 µm`. The direction average `(x̄↑ + x̄↓)/2` adds two
quantities biased in opposite senses, so they cancel - measured shift
`+0.000 µm`.

That single fact explains the whole damage pattern. Everything built on
direction averages survived; everything built on the direction difference did
not.

## Corrected values

"Tail clipped" keeps the published window and truncates it at measured rest -
the minimal change. "Measured holds" is the independent blind path. They agree
on every quantity that the defect touches.

| quantity (µm) | published | tail clipped | measured holds |
|---|---:|---:|---:|
| `B`   reversal value | 0.982 | **2.264** | **2.264** |
| `B_mean` mean signed `B_i` | -0.365 | **-1.600** | **-1.601** |
| `R↑`  unidirectional repeatability, up | 2.098 | 2.013 | 2.012 |
| `R↓`  unidirectional repeatability, down | 1.850 | 1.891 | 1.891 |
| `R`   bidirectional repeatability | 2.427 | **3.616** | **3.617** |
| `E↑`  systematic deviation, up | 6.116 | 6.068 | 6.192 |
| `E↓`  systematic deviation, down | 6.611 | 6.610 | 6.735 |
| `E`   systematic deviation | 6.611 | **7.729** | **7.855** |
| `M`   mean bidirectional deviation range | 6.364 | 6.339 | 6.463 |
| `A`   accuracy | 7.949 | **8.901** | **9.151** |
| median `s_i` | 0.280 | 0.283 | 0.283 |
| max `s_i` | 0.524 | 0.503 | 0.503 |

Unaffected to within 0.1 µm: `R↑`, `R↓`, `E↑`, `E↓`, `M`, the per-target
direction averages, the repeatability spreads, and the Step A gain
(-24.17 ppm becomes -24.10 ppm). Affected: `B`, `R`, `E`, `A`.

## Why it survived review

`B_mean` as a function of the fitted clock offset is a cliff, and the published
offset sits one search-grid step past its edge:

| offset vs published | `B_mean` (µm) |
|---:|---:|
| -100 ms | -1.601 |
| -75 ms | -1.601 |
| -50 ms | -1.600 |
| -25 ms | -1.537 |
| **0 ms (published)** | **-0.365** |
| +25 ms | +3.497 |
| +50 ms | +11.125 |

Flat and correct for every offset up to -50 ms, then it collapses. Local
sensitivity is `+47 nm of B per millisecond`. `fit_clock` searches the offset on
a `np.arange(..., 25.0)` grid - **1.17 µm of `B` per grid step**, for a quantity
reported to three decimal places in micrometres.

Nothing in the evaluation's own output looked wrong. It reported a median
arrival residual of 8.6 ms, which reads as a good fit until one notices that the
quantity being measured needs the window placed to a few milliseconds.

Changing only the smoother width inside `fit_clock`, with everything else
untouched, walks the answer back:

| smoother | `B_mean` (µm) |
|---|---:|
| 51-sample (as published) | -0.365 |
| 21-sample | -1.194 |
| 11-sample | -1.537 |
| 5-sample | -1.537 |
| measured holds (blind) | -1.601 |

## Independent confirmation

Four runs, three grids, two firmware generations, all evaluated against measured
holds:

| run | grid | mean reversal (µm) | targets |
|---|---|---:|---:|
| E1 (spectral windows) | 0.12 mm | -1.681 | 303 |
| E2 (global map) | 0.04-0.56 mm | -1.605 | 400 |
| ISO 230-2 v3, re-derived | 4.891 mm | -1.601 | 24 |
| ISO 230-2 v4 | 4.891 mm | -1.819 | 24 |

The v4 run was commissioned specifically to test whether the larger reversal
repeats in the ISO sequence. It does. Per-target `B_i` correlates `r = +0.871`
between v3 and v4, so the reversal is a repeatable position-dependent property
of the stage rather than scatter.

## Corrective action

Done:

- `iso230_2_evaluate_blind.py` - ISO 230-2 evaluation from the capture alone.
  Dwells from the trace, approach direction from the sign of the last 150 ms of
  travel, targets from the planned grid, averaging window inside the measured
  hold. Used for v4, which again arrived without a board log.
- `check_e1_no_board_clock.py`, `check_iso_control.py`,
  `render_window_independence.py` - the cross-checks above.
- Results recorded in `iso230_2_v3_results_blind.json` and
  `iso230_2_v4_results_blind.json`.

Outstanding:

- `iso230_2_evaluate.py` is **unchanged**, so `iso230_2_v3_results.json`,
  `per_target_statistics.csv` and the `B`, `R`, `E`, `A` rows of `E0/TABLES.md`
  still carry the defective values. They must be regenerated or withdrawn before
  any of those four quantities is quoted again.
- `render_command_recorded_overlay.py` and `analyze_open_loop_diagnostic.py`
  have not been audited for the same window-placement assumption.

Recommended fixes, in order of value:

1. **Assert on the window, do not trust the clock.** Reject any record window
   whose peak-to-peak exceeds about 0.3 µm, and trim trailing motion. This
   catches the entire class of fault regardless of what the clock fit does, and
   would have failed loudly on all 240 windows here.
2. **Give the window a guard band in firmware.** Wait about 100 ms after
   `record_end` before commanding the next move, so that a timing error of a few
   tens of milliseconds cannot reach the average.
3. **Sharpen the fit.** A 5-sample motion-end smoother and a 1 ms offset grid in
   `fit_clock`.

Fix 1 is the important one. Fixes 2 and 3 reduce the probability of the error;
fix 1 makes it detectable when it happens anyway.

## Reproducing

```bash
cd "stage identification/iso230-2/analysis"

# blind evaluation of either ISO run (needs the raw capture, kept local-only)
python3 iso230_2_evaluate_blind.py --csv ../../ISO230-2_v4.csv \
        --json iso230_2_v4_results_blind.json
python3 iso230_2_evaluate_blind.py --csv ../../ISO230-2_v3.csv \
        --json iso230_2_v3_results_blind.json

cd "../../Geometric Fitting/experimental_scripts/e1_spectral_windows"
python3 check_e1_no_board_clock.py      # E1 with the log out of the chain
python3 check_iso_control.py            # ISO v3 through the same code
python3 render_window_independence.py   # the figure, all three runs
```

## Operational lessons

1. **A derived quantity needs a validity assertion, not just a plausible
   pipeline.** Every intermediate here looked healthy - the clock fit residual,
   the per-target spreads, the repeatability. The only thing that would have
   caught it is a check on the thing actually being averaged.
2. **Never let an analysis window end where motion begins.** A window whose edge
   coincides with a commanded event has no margin for timing error in either
   direction.
3. **A centred smoother cannot mark the time of an edge.** Use a short window,
   or detect the edge and correct for the known group delay.
4. **Check the sensitivity of a result to its own free parameters.** A one-line
   sweep of the clock offset would have shown that `B` moves 1.17 µm per search
   grid step, which is itself disqualifying.
5. **A disagreement between two measurements of one machine is a finding, not a
   nuisance.** The 4.6x gap was the only reason any of this was looked at.

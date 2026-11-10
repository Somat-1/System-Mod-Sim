# MRES settling matrix — results

## Result

The completed matrix does **not show a consistent improvement in settled
landing accuracy as MRES is increased from 1/16 to 1/256**. The four curves
cross repeatedly with distance: the best setting is 1/64 at 0.625 µm, 1/32 at
40 µm, 1/256 at 640 µm, 1/64 at 5 mm, and 1/16 at 20 mm. That is not the
monotonic ordering expected if microstep resolution were the dominant landing
error source.

However, the run also exposes a critical limitation: **the executed command
velocity was not held constant across MRES**. The long-move logical command
rates were 0.625, 0.625, 0.608, and 0.540 mm/s at 1/16, 1/32, 1/64, and 1/256,
respectively, rather than the planned 0.78 mm/s. The apparent aggregate
advantage at 1/256 is therefore partly confounded with a 14% slower approach
than 1/16 and a 31% shortfall from the planned rate.

Across the 176 same-direction moves used for the main comparison, the median
absolute error is **0.268 µm**, the 95th percentile is **2.138 µm**, and the
maximum is **2.587 µm**. Mean absolute error by setting is 0.622, 0.642, 0.583,
and 0.534 µm for MRES 1/16, 1/32, 1/64, and 1/256 respectively. Although 1/256
has the lowest aggregate value, its 0.108 µm separation from the worst setting
is not reproduced consistently across the individual distances and should not
be treated as a general resolution advantage.

The observed data do not support selecting 1/256 on the expectation of lower
settled-position error. They are compatible with using the coarsest setting
that resolves the required increment, but the experiment cannot make a clean
causal statement about MRES until the velocities are matched in execution.

## Most important figures

### Same trajectories, overlaid by MRES

The curves are direction-normalized and referred to each move's own pre-move
settled position. Solid lines are the median of the same-direction
repetitions; translucent bands span the interquartile range. The dashed line
is the intended endpoint.

![Direction-normalized trajectory montage](plots/01_trajectory_montage.png)

### Offset from the instantaneous commanded trajectory

This montage subtracts the terminal log's instantaneous logical command from
the encoder trajectory. It is useful for comparing the executed shapes, but
the transient residual includes the finite synchronization accuracy between
the 20 Hz terminal telemetry and the 1 kHz TwinCAT recording. It is not a
closed-loop following-error measurement.

![Instantaneous command tracking error](plots/02_command_tracking_error_overlays.png)

### Critical execution limitation: achieved velocity

The long trajectories make the velocity mismatch visible. The table below is
computed from the controller's logical-position telemetry between 10% and 90%
of the 5 mm and 20 mm moves.

| MRES | Planned [mm/s] | Achieved [mm/s] | Shortfall |
|---:|---:|---:|---:|
| 1/16 | 0.780 | 0.625 | 19.9% |
| 1/32 | 0.780 | 0.625 | 19.9% |
| 1/64 | 0.780 | 0.608 | 22.1% |
| 1/256 | 0.780 | 0.540 | 30.8% |

The firmware explains the trend. In
[`moveToU256()`](../../firmware/settling_test/settling_test.ino), each pulse
is emitted and then the code starts a full nominal-period wait from the
current `micros()` value. Per-pulse guard, GPIO, telemetry, and serial-service
work therefore extends the rising-edge interval rather than being included in
it. The fractional penalty grows as the requested pulse period becomes
shorter, so fine MRES runs slower.

This is a test-method limitation, not evidence that the driver ignored MRES:
the independent MSCNT pulse checks passed at all four settings.

### Settled offset from the intended endpoint

This is the direct view of the campaign's main quantity. Faint curves are
individual same-direction landings; bold curves are the MRES medians. Each
trace covers the final 1.5 s of the dwell and is expressed relative to the
intended displacement.

![Settled endpoint offset overlays](plots/02_endpoint_offset_overlays.png)

### Numerical summary

![Mean absolute error against distance](plots/03_offset_summary.png)

## Mean absolute landing error

Only same-direction moves are included. Values are in micrometres. Each cell
contains six moves (three per direction), except 20 mm, which contains two
(one per direction) because travel limited the design.

| Command | 1/16 | 1/32 | 1/64 | 1/256 |
|---|---:|---:|---:|---:|
| 0.625 µm | 0.168 | 0.197 | 0.147 | 0.152 |
| 2.5 µm | 0.236 | 0.271 | 0.191 | 0.214 |
| 10 µm | 0.197 | 0.278 | 0.268 | 0.202 |
| 40 µm | 0.191 | 0.107 | 0.183 | 0.142 |
| 160 µm | 0.293 | 0.370 | 0.135 | 0.224 |
| 640 µm | 1.228 | 1.283 | 1.288 | 0.858 |
| 5 mm | 2.134 | 1.999 | 1.864 | 1.924 |
| 20 mm | 0.357 | 0.617 | 0.595 | 0.600 |

The fact that 20 mm performs better than 5 mm for every MRES setting rules out
a simple distance-proportional scale error as the main explanation. The
landing error depends on trajectory history and direction, consistent with
friction/backlash effects being more important than the selected microstep
resolution.

## Aggregate statistics by MRES

| MRES | Mean \(|e|\) [µm] | Median \(|e|\) [µm] | P95 \(|e|\) [µm] | Gain error [ppm] | Gain-fit residual RMS [µm] |
|---:|---:|---:|---:|---:|---:|
| 1/16 | 0.622 | 0.260 | 2.114 | +12.8 | 0.959 |
| 1/32 | 0.642 | 0.356 | 2.176 | +27.4 | 0.944 |
| 1/64 | 0.583 | 0.260 | 2.116 | +24.9 | 0.909 |
| 1/256 | 0.534 | 0.256 | 1.910 | +14.6 | 0.849 |

The fitted gains are all within 28 ppm of unity and again have no monotonic
ordering with MRES.

## Direction and reversal

![Forward and reverse landing offsets](plots/04_directional_offset_montage.png)

The first move after each reversal is intentionally tagged separately. Its
mean absolute error is 1.202 µm, versus 0.595 µm for the same-direction data;
including it in the MRES mean would therefore mix backlash/approach history
into the resolution comparison. The main tables and conclusions exclude all
64 reversal-tagged moves.

## Acquisition and analysis settings

| Item | Value |
|---|---|
| Scope source | `data/MRESdebug.csv` |
| Matching terminal log | `data/terminal_20261022_171207.log` |
| Scope interval | 2026-10-22 17:20:03.551–18:00:19.900 |
| Scope sampling | 1.000 ms; 2,416,339 samples |
| Counter scale | 1 count = 1 nm |
| Counter rollovers | 0 |
| Matrix | 344 rows; 240 measured moves |
| Main comparison | 176 same-direction moves |
| MRES order | 1/16 → 1/32 → 1/64 → 1/256 |
| Physical MRES check | Passed before the run at all four settings |
| Nominal motion speed | 0.78 mm/s (78 full steps/s) |
| Achieved logical speed | 0.625 / 0.625 / 0.608 / 0.540 mm/s for 1/16 / 1/32 / 1/64 / 1/256 |
| Dwell | 3.0 s on measured rows |
| Settled value | Mean of final 1.0 s, ending 0.10 s before the STOP boundary |
| Reversal pre-window | 0.8 s where the setup row has only a 1.0 s dwell |
| Driver current | 400 mA RMS |
| Chopper | SpreadCycle |
| StealthChop | Off |
| Interpolation / microPlyer | Off |

The ESP32 board clock was aligned to the independent TwinCAT trace from the
full commanded trajectory. The fitted mapping is
`scope_time = board_time − 477.946000 s`. The full-trajectory
encoder-versus-command RMS is 15.98 µm; this number includes motion transients
and telemetry timing, and is not the settled endpoint uncertainty.

![Recording and clock alignment](plots/05_recording_and_alignment.png)

For the same-direction rows used in the comparison, the worse of the two
endpoint-window standard deviations has a median of **9.0 nm**, a 95th
percentile of **20.1 nm**, and a maximum of **37.4 nm**. The reported
sub-micrometre landing differences are therefore resolved by the encoder
recording. No counter rollover, clipping, or saturation signature was found.

## Run integrity

The matching log explicitly records:

- firmware `settling_test`, 344 rows and 240 measured moves;
- physical one-pulse/five-pulse MSCNT checks passing at 1/16, 1/32, 1/64,
  and 1/256;
- MRES switches at all four planned block boundaries;
- `RUN_COMPLETE,pos=0` after row 343;
- no `ERR`, abort, or failed MRES-check record during the completed run.

This removes the earlier ambiguity that existed before the matching log was
available: block identity and MRES switching are now confirmed by the
controller log, not inferred from the scope shape alone.

## Limitations

- The planned common speed was not achieved. MRES is confounded with the
  executed approach velocity, most strongly at 1/256. A causal MRES comparison
  requires a rerun with a deadline-based pulse scheduler or calibrated
  per-MRES speed commands.
- The result otherwise applies to 400 mA RMS, SpreadCycle, interpolation off,
  and StealthChop off.
- MRES blocks were executed once in a fixed order rather than repeatedly
  interleaved. A slow thermal or friction drift is therefore partially
  confounded with MRES.
- The 20 mm same-direction cells have only two observations, so their apparent
  ranking is weak evidence.
- The instantaneous command-offset montage is synchronization-sensitive.
  Settled differences are the primary result because they are insensitive to
  the absolute datum and much less sensitive to sub-sample clock alignment.

## Reproducibility

Run:

```powershell
$env:MPLBACKEND = 'Agg'
python .\analysis\mres_debug\analyze_mres_debug.py
```

The script writes only PNG figures, plus the two result CSVs and JSON summary.
No PDF output is generated.

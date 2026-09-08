# v2 chirp -- sweep configuration log

Generated 2026-09-08 12:19 UTC by `generate_sweep_config_log.py` from `chirp_v2_schedule.py`. **Do not hand-edit** -- change the schedule module and regenerate instead, or this stops being trustworthy as a processing reference.

## Mechanics / drive

| Parameter | Value |
|---|---:|
| Lead | 2 mm/rev |
| Full steps/rev | 200 |
| Full step size | 10 um (0.01 mm) |
| MRES | 16 (1/16 step) |
| 1 microstep | 0.0625 full steps (0.625 um) |
| Run current | 400 mA RMS |

## Plant placeholders -- REPLACE with same-day ringdown measurements

| Parameter | Placeholder value | Source |
|---|---:|---|
| f_n | 176.7 Hz | Prior analytical model; not measured |
| Q | 20 | Prior analytical model; not measured |
| zeta | 0.025 | = 1/(2*Q) |

## Amplitude profile

`A_cmd(f) = min(cruise * notch(f), stroke_clamp, step_rate_clamp / (2*pi*f*MRES))`, `notch(f)` a raised-cosine window centered on f_n.

| Parameter | Value |
|---|---:|
| Cruise levels | 3 full steps |
| Notch half-width | 56.7 Hz |
| Notch band (at f_n placeholder) | 120.0 - 233.4 Hz |
| Notch depth ratio | 0.1 (amplitude at f_n = 10% of cruise) |
| Stroke clamp | 0.5 mm (50 full steps) |
| Step-rate clamp | 200,000 Hz |
| Velocity/presliding clamp | not applied by default (see README Section 3) |

Minimum amplitude actually reached per cruise level (from `amplitude_profile`, swept 1-1000 Hz):

| Cruise (full steps) | Minimum amplitude | At frequency | Above 1 microstep? |
|---:|---:|---:|:---:|
| 3 | 4.8 microsteps | 176.8 Hz | yes |

Stribeck-velocity reference (not applied as a default clamp -- see README Section 3):

- Nut: 0.2 mm/s
- Guideway: 0.25 mm/s

## Sweep frequency law

| Segment | Range | Duration |
|---|---|---:|
| Log | 1 -> 60 Hz | 120 s |
| Linear | 60 -> 1000 Hz | 120 s (7.833 Hz/s) |

Down-sweep is the exact mirror (linear 1000->60 Hz, then log 60->1 Hz).

## Sync marker

| Parameter | Value |
|---|---:|
| Frequency | 750 Hz |
| Amplitude | 2 full steps |
| Burst count | 3 |
| Burst duration | 0.1 s |
| Gap duration | 0.1 s |
| Total marker duration | 0.5 s |

## Level timing breakdown (one amplitude level, one direction pass)

| Segment | Start (s) | End (s) | Duration (s) |
|---|---:|---:|---:|
| idle pre-roll | 0.00 | 30.00 | 30.00 |
| sync marker 1 | 30.00 | 30.50 | 0.50 |
| up-sweep: log segment | 30.50 | 150.50 | 120.00 |
| up-sweep: linear segment | 150.50 | 270.50 | 120.00 |
| mid-dwell | 270.50 | 275.50 | 5.00 |
| sync marker 2 | 275.50 | 276.00 | 0.50 |
| down-sweep: linear segment | 276.00 | 396.00 | 120.00 |
| down-sweep: log segment | 396.00 | 516.00 | 120.00 |
| tail (idle) | 516.00 | 546.00 | 30.00 |
| **Total** | | | **546.00 s (9.10 min)** |

Resonance crossing (f_n placeholder = 176.7 Hz): up-sweep at t=165.40 s, down-sweep at t=381.10 s (both relative to level start).

## Full campaign (all cruise levels)

1 cruise level x 546.0 s = 9.1 min of sweeping, plus ringdown characterization time (Section 1) not counted here.


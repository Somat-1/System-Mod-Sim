# ISO 230-2 initial-run post-mortem

## Executive conclusion

The large discrepancy in the first ISO 230-2 render was produced by two independent effects:

1. The original 5 mm/s motion exceeded the approximately 2.5 mm/s counting or output-rate ceiling visible in the interferometer/EL5101 recording path. The recorded counter therefore lagged the physical trajectory and continued to advance during the short stop window.
2. The ESP32 and TwinCAT scope did not share a trigger or clock. The first overlay treated their time axes as though they were already synchronized, which displaced the trajectories in time. A later point-by-point time warp made the short reference ramps look unrealistically fast. Only a single affine clock mapping preserves the commanded ramp durations.

The diagnostic data do **not** support the hypothesis that the controller simply stopped issuing steps before each target. Endpoint agreement across all tested profiles and repeatable home-switch closures instead show that the stage completed the commanded transfers, while the original measurement path and overlay made them appear incomplete.

## What was observed in the first run

- The raw interferometer trace appeared to cover less distance than the controller command during the large ramps.
- The small fixed-side reference moves appeared displaced relative to the main ISO ramps.
- In the early synchronized render, the dashed reference ramps appeared much faster than the corresponding controller ramps.
- Some intended stop windows did not contain a clearly settled plateau.

These symptoms initially resembled missed steps, insufficient move time, or a reference-position failure. They cannot be diagnosed from an unsynchronized overlay alone.

## Root cause 1: recorded velocity ceiling

The follow-up `IsoTest.csv` diagnostic repeated the same 59.35 mm transfer with four profiles:

| Profile | Commanded speed | Acceleration | Recorded full-transfer speed | Stable hold visible |
|---|---:|---:|---:|---:|
| Original | 5 mm/s | 40 mm/s² | 2.472 mm/s | about 1.2-1.5 s |
| Low acceleration | 5 mm/s | 10 mm/s² | 2.462 mm/s | about 1.7-2.3 s |
| Low speed | 2 mm/s | 40 mm/s² | 1.582 mm/s | about 9.4-9.5 s |
| Conservative | 1 mm/s | 5 mm/s² | 0.810 mm/s | about 9.5-9.6 s |

For both 5 mm/s profiles, 97.8-98.4% of the samples during a representative ramp were pinned at an exact 2500-count increment. With the nominal 1 nm/count scale and 1 kHz export rate, that is an exact 2.500 mm/s derivative ceiling. Reducing acceleration while retaining 5 mm/s did not remove it. Reducing the maximum speed did.

The important consequence is temporal: the counter output continued catching up after the stage had completed its pulse train. The original stop consisted of only 0.5 s settling plus a 0.5 s record window. That was too short to expose a clean plateau after a rate-limited ramp, so the following stop could appear offset or incomplete.

This ceiling belongs to the observed metrology/output/acquisition chain. The available CSV alone cannot identify whether its exact hardware origin is the interferometer quadrature output, EL5101 input path, or their configuration. Calling it a measurement-path ceiling is therefore more defensible than assigning it to one device without a direct register or bandwidth test.

## Root cause 2: independent clocks and misleading synchronization

The ESP32 reports command time using its millisecond clock. TwinCAT reports the interferometer samples on the scope clock. There was no common trigger. Serial host timestamps were also unsuitable as a precise trigger because buffered terminal output can arrive after the corresponding motion.

Three different time representations must therefore remain distinct:

- **Raw scope time:** the interferometer samples exactly as exported.
- **Controller time:** the commanded trajectory in the ESP32 clock.
- **Corrected overlay time:** controller time mapped to scope time.

The reliable correction is:

1. detect stationary plateaus independently in the raw trace;
2. order-match them to controller stops using position and sequence order;
3. fit one clock rate and one offset;
4. apply that single affine map to the complete command trace.

Using a separate local warp between successive stops is not acceptable because it changes individual ramp durations. That was why the dashed reference ramps in an earlier panel appeared implausibly fast even though the controller panel showed longer ramps.

## Evidence against missed steps as the primary cause

- All four diagnostic profiles converged to essentially the same transfer span: 59.3372-59.3384 mm.
- The two repeats within each profile differed by at most about 1.03 µm.
- Ten home-switch-like closures covered only 4.99 µm peak-to-peak.
- Lowering acceleration from 40 to 10 mm/s² while retaining 5 mm/s did not change the clipped behavior.
- The 2 mm/s profile produced long, directly visible stationary plateaus.
- The complete v3 trace spans 122.706 mm, consistent with full travel.

The common approximately 12 µm offset of the diagnostic transfer from the nominal 59.35 mm is profile-independent. It is therefore a small common scale/datum/end-point issue, not evidence that the high-speed profile alone lost motion.

This does not prove that an open-loop stage can never miss a step. It shows that accumulated lost steps or coupling slip were not the dominant mechanism in this event. The firmware can prove that pulses were emitted, while only the interferometer and end switches can verify mechanical following.

## Why the original short dwell was insufficient

The firmware does not use dwell time to decide whether a move is complete. `moveX()` emits the complete pulse train first; only after it returns does `dwellRecord()` perform settling and recording. The next move starts after both windows finish.

Consequently:

- a longer dwell cannot recover a genuinely missed motor step;
- it can reveal mechanical settling;
- it can prevent a rate-limited metrology output from contaminating the next record window;
- it provides a stable plateau that can be paired robustly during analysis.

The original 1.0 s total stop window was inadequate for the clipped 5 mm/s recording. The revised profile uses 2.0 s settling followed by a 0.5 s record window.

## Corrective action

The validated ISO configuration was changed to:

- velocity: **2.0 mm/s**;
- acceleration: **40 mm/s²**;
- settling before recording: **2.0 s**;
- record window: **0.5 s**;
- automatic startup: countdown, home, centre hold, then the ISO sequence.

The speed change addresses the observed rate limit. The acceleration was not reduced because the diagnostic isolated maximum speed, rather than acceleration, as the discriminating parameter. The longer settle window gives the acquisition and mechanical system margin before each formal record window.

## Validation from ISO230-2 v3

The corrected v3 acquisition produced:

- 2,321,285 raw samples at 1 kHz over 38.69 minutes;
- all **251/251** planned record stops paired to stable plateaus;
- all **240/240** formal ISO target stops resolved;
- ISO arrival error: median **-0.02 µm**, RMS **1.63 µm**, maximum absolute error **4.77 µm** after one constant absolute-position datum;
- no position-gain fitting, rescaling, filtering, or interpolation.

There are 27 missing 1 ms samples clustered near 80 s in the scope export. Their original time gaps are retained. They occur in the pre-ISO context and do not prevent any planned ISO stop from being resolved.

The v3 serial terminal log was not saved. Its command trace is therefore reconstructed from the flashed configuration and exact generated target table, then placed on the scope clock using the plateau-based affine alignment. This is sufficient for the full-trajectory validation above, but a matching serial log should be captured in future runs so the exact firmware record windows and driver messages can be audited directly.

## Operational lessons

- Start the interferometer and serial logger before resetting the ESP32.
- Preserve the raw scope CSV, terminal log, target table, and configuration hash together for every run.
- Do not infer missed steps from a raw unsynchronized overlay.
- Do not use local time warping that changes individual ramp duration.
- Keep end switches as safety and pre/post-run integrity checks. Do not home between formal ISO runs, because re-homing changes the datum and approach history and can hide accumulated error.
- If a future trace again shows a 2500-count-per-sample ceiling, inspect the interferometer quadrature-output rate and EL5101 input/filter configuration before changing the mechanical trajectory.

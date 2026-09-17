# Unnotched 250–1000 Hz Chirp

`FB_Chirp250To1000.st` is based on the safe high-frequency leg of ESP chirp v3. It does not contain the v3.2 notch and never traverses the lower-frequency resonance band.

## Sequence

1. Hold the captured origin for 5 s.
2. Linear up-chirp from 250 Hz to 1000 Hz over 96 s.
3. Hold the origin for 5 s.
4. Linear down-chirp from 1000 Hz to 250 Hz over 96 s.
5. Hold the origin for 5 s.

Default amplitude is 3 full steps peak. The phase is evaluated analytically on every PLC scan, avoiding accumulated phase-integration drift. Normal completion commands the captured origin.

## No notch

There is no amplitude notch, taper, or commanded frequency below 250 Hz. The requested band is generated directly. At the default 96 s duration, the sweep rate is 7.8125 Hz/s, matching the high-band leg of v3.

## Microstep resolution during the chirp

The microstep size does **not** change during the chirp. `uiTerminalMicrostepsPerFullStep` remains fixed at 64, so one terminal position count always represents 1/64 of a full step. The default ±3-full-step amplitude therefore remains ±192 terminal counts throughout both sweeps.

Only the commanded frequency changes. Consequently, the difference between consecutive 250 µs position targets becomes larger at higher frequency, but the underlying position-count size remains 1/64 step. At 1 kHz the sampled command is approximately `0, +192, 0, -192` counts per period.
## Timing limitation

The supplied EL70x1 manual states a 250 us internal cycle. With a 250 us PLC/EtherCAT task:

- 250 Hz receives 16 setpoint samples per period.
- 400 Hz receives 10 samples per period.
- 1000 Hz receives only 4 samples per period.

The 1 kHz endpoint is therefore a coarse four-point waveform. The function block requires the documented 250 us task period and requires explicit acknowledgment whenever the endpoint has fewer than ten samples per period.

## Speed check

For a sinusoidal position command, peak speed is `2*pi*f*A`. At 1000 Hz and 3 full steps peak, this is about 18,850 full steps/s. The code assumes the terminal speed range is configured for 32,000 full steps/s and rejects a request that consumes more than 90% of the declared range.

This speed check does not establish motor/load stability. The acceleration demand remains severe and open-loop pole slip is still possible. Treat this as an experimental sequence, start with reduced amplitude/bandwidth, and record the independent position sensor.

## Error IDs

- `16#0202`: task period is not 250 us.
- `16#0204`: requested frequencies leave the allowed 250–1000 Hz band.
- `16#0205`: fewer than four command samples per endpoint period.
- `16#0206`: marginal sampling was not acknowledged.
- `16#0208`: peak sine velocity exceeds the configured 90% speed margin.
## Rendered previews

![Complete chirp sequence](Rendered/chirp_sequence_overview.png)

![Chirp cyclic-sampling detail](Rendered/chirp_sampling_detail.png)

The source is `render_chirp_sequence.py`; generated files and regeneration instructions are in `Rendered/README.md`.
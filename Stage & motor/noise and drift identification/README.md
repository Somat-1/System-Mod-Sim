# Noise and drift identification

This area holds long-duration environmental and standstill measurements that
are used to quantify noise, drift, and disturbances independently of the
motion-identification campaigns.

## Layout

- `data/` — immutable raw measurement exports.
- `scripts/` — reproducible inspection and analysis code.
- `analysis/cache/` — compact display summaries rebuilt from raw data.
- `analysis/plots/` — rendered measurement sequences and event close-ups.
- `analysis/exclusion_windows.csv` — explicit time masks. Raw samples are never
  deleted; subsequent analyses must read this file and report which masks were
  applied.

## Current recording

`data/2026-10-15_seis/TvalSeis.csv` is a three-channel SEIS voltage recording
sampled at 20 kHz. The operator-confirmed table disturbance from 26:00 to
30:00 (`1560–1800 s`) is recorded as an active exclusion. The raw file remains
unchanged; processing scripts apply the time mask explicitly and never bridge
the resulting gap with an analysis window.

The table and floor frequency-domain results are consolidated in
`analysis/plots/seis_frequency_montage.png`. It presents direct FFT amplitude
spectral density and Welch PSD across the complete measured band. The SEIS/table
and floor recordings occupy separate column groups, while channels 1–3 occupy
separate rows. Every panel contains one measurement and one axis only, with a
faint estimate and a heavy smoothed per-axis trace. A raw
24-32 minute table context crop with the active 26-30 minute interval omitted
is rendered as `analysis/plots/seis_04_raw_exclusion_crop.png`.

`data/2026-10-15_seis_floor/TvalSeisfloor.csv` is a second three-channel,
20 kHz recording made with the same device on the floor. It contains
19,444,981 samples over 972.249 s (16.20 min). No interval is excluded. Its
raw overview is `analysis/plots/seis_floor_01_recorded_sequence.png`. Renderers
create PNG files by default; PDF output is opt-in with `--pdf` where supported.

A quantitative comparison of the two placements is provided in
`analysis/SEIS_TABLE_VS_FLOOR_COMPARISON.md`.


# Preemptive run visualization

- IDS source: `data/PreemptiveRundata.csv` (1,115,750 samples at 1 kHz; 1115.749 s)
- ESP source: `data/esp32_runs/full_campaign_marked_complete_20260828.txt`
- ESP campaign duration captured: 1058.160 s
- Estimated ESP `CAMPAIGN_START` in IDS time: 56.195 s
- Alignment derivative correlation: 0.2070
- Fitted encoder counts per commanded 1/16-full-step unit: 0.0336454
- Extracted labeled segments: 61 (60 complete, 1 partial)
- Run ended by commanded abort during run 3, A1, N2 negative; later tests are absent, not inferred.

The command overlays use a separate full-step axis; no encoder calibration is assumed.
Segment boundaries come from the ESP timestamps shifted by the marker-derived clock offset. `segment_manifest.csv`
contains both ESP-relative and IDS-relative boundaries and should be treated as the
machine-readable splice index.

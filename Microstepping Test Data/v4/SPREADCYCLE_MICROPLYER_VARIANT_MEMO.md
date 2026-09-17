# Memo — SpreadCycle + MicroPlyer variant of the v4 campaign (MRES 1/4/16)

## Scope

Per [SPREADCYCLE_VARIANT_MEMO.md](SPREADCYCLE_VARIANT_MEMO.md), driver-mode
changes are cloned into a separately versioned sketch rather than retrofitted
into the baseline. This is a second such clone, combining two changes
against the baseline `esp32_v4_mres_trajectory_campaign/`:

1. **MRES 32 is dropped.** This variant sweeps only MRES `1, 4, 16`. The
   baseline's four-MRES sweep left only ~6.4 min of headroom under the
   55 min recording limit / 2 min safety margin
   (`RECORDING_LIMIT_S - PLANNED_MARGIN_S` in
   [run_mres_trajectory_campaign.py](scripts/run_mres_trajectory_campaign.py)),
   not enough to add MicroPlyer runs on top of all four MRES values.
   Dropping MRES 32 frees enough time that MicroPlyer can be enabled for
   this entire three-MRES sequence instead of needing a separate reduced
   scope for just the added runs.
2. **SpreadCycle with MicroPlyer interpolation ON for the whole
   sequence.** `configureDriver()` calls `driver.en_spreadCycle(true)` and
   `driver.intpol(true)` — SpreadCycle chopping and MicroPlyer step
   interpolation both apply throughout, unlike
   `esp32_v4_mres_trajectory_campaign_spreadcycle/`, which disables
   interpolation (`intpol(false)`).

Every MRES 1/4/16 oscillation, marker, trajectory rate, dwell, and
separation is otherwise identical to the baseline, so the reduced-MRES
dataset stays comparable to the corresponding MRES 1/4/16 subset of the
baseline.

## Planned duration

`scripts/plot_planned_sequence_spreadcycle_microplyer.py` builds the same
plan model as the baseline's `plot_planned_sequence.py`, with `MRES_VALUES =
(1, 4, 16)` instead of `(1, 4, 16, 32)`. The complete planned sequence is
**2,090.2 s (34.84 min)**, about 20.16 min under the 55 min recording limit
(18.16 min under the 53 min effective budget) — comfortably inside budget,
unlike the ~6.4 min of headroom left in the four-MRES baseline.

- `rendered_assets/planned_mres134_spreadcycle_microplyer_campaign.png` —
  complete preview: three oscillation panels (MRES 1, 4, 16), the full
  18-trajectory plot (3 MRES x 2 modes x 3 rates), and the overall schedule.

## What changed vs. the baseline sketch

- `MRES_VALUES[] = {1, 4, 16}` (was `{1, 4, 16, 32}`).
- `configureDriver()`: added `driver.en_spreadCycle(true);` and
  `driver.intpol(true);`.
- Boot banner changed to
  `# ESP32_V4_MRES134_SPREADCYCLE_MICROPLYER_CAMPAIGN` and a
  `# DRIVER_MODE,spreadcycle=1,intpol=1` line is printed once the driver
  comes up.
- The `mode` field logged for `CAMPAIGN_START`/`RUN_CONFIG`/`RUN_COMPLETE`/
  `CAMPAIGN_COMPLETE` events is `SPREADCYCLE_MICROPLYER` (per-block `mode`
  fields for markers, oscillation, and the direct/individual trajectories
  are unchanged).
- Everything else — rates, dwells, separations, marker amplitudes,
  throughput preflight — is byte-for-byte the same as the baseline sketch.

## Validation status

The plan model's own internal invariants (every marker, oscillation, and
trajectory block returns to its start position; total duration under the
recording limit) are checked by `build_plan()` and pass. No ESP32 hardware
dry-run or live-run CSV exists yet for this variant — `latest_dry_run()` in
the plot script looks for
`data/hardware_runs/mres134_spreadcycle_microplyer_dry_run_*.csv` and, if
none is found, renders the plan model unvalidated against a captured run and
says so. Once this sketch is compiled, flashed, and run (or dry-run through
the ESP32's own serial output), save its CSV under that name to enable the
existing `validate_dry_run()` checks used by the baseline.

## Comparison notes

Compare against the MRES 1/4/16 subset of the baseline, and against
`esp32_v4_mres_trajectory_campaign_spreadcycle/`, only when mechanics,
supply voltage, current, origin, preload, and acquisition settings are
unchanged. Because both SpreadCycle *and* interpolation change together
here relative to the baseline, a difference between this dataset and the
baseline should not be attributed to SpreadCycle alone — the
`..._spreadcycle` (interpolation off) variant isolates that.

## Build and upload

```sh
arduino-cli compile \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  scripts/esp32_v4_mres134_spreadcycle_microplyer_campaign

arduino-cli upload \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  --port /dev/cu.usbmodem101 \
  scripts/esp32_v4_mres134_spreadcycle_microplyer_campaign
```

Status: sketch and plan/plot added, not yet compiled, flashed, or run.

# Memo — SpreadCycle variant of the v4 campaign

## Scope

The baseline v4 MRES/trajectory campaign
(`scripts/esp32_v4_mres_trajectory_campaign/`) leaves the TMC2209 at its
power-on driver-mode defaults (StealthChop enabled, MicroPlyer interpolation
on) and does not touch chopper mode, per its own header comment. Per
[STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md](STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md),
driver-mode changes are not retrofitted into that baseline sketch — they are
cloned into a separate, separately versioned sketch so the baseline dataset
stays reproducible and comparable.

`scripts/esp32_v4_mres_trajectory_campaign_spreadcycle/` is that clone. It
repeats exactly the same measurement plan as the baseline — MRES `1, 4, 16,
32`, the 15-cycle/30 s oscillations, and both 25 mm command-generation
approaches (direct and individual full-step) at 27.5, 70, and 200 full
steps/s — with one difference: `configureDriver()` explicitly calls
`driver.en_spreadCycle(true)` and `driver.intpol(false)`, so SpreadCycle
chopping and disabled step interpolation apply for the entire sequence
instead of the power-on default.

MicroPlyer-enabled (`intpol(true)`) 25 mm runs are **not** included in this
sketch or the baseline — they did not fit the ~6.4 min of remaining headroom
under the 55 min recording limit / 2 min safety margin
(`RECORDING_LIMIT_S - PLANNED_MARGIN_S` in
[run_mres_trajectory_campaign.py](scripts/run_mres_trajectory_campaign.py))
alongside the existing baseline content, and are being run as a separate,
additional sequence.

## What changed vs. the baseline sketch

- `configureDriver()`: added `driver.en_spreadCycle(true);` and
  `driver.intpol(false);` after the existing current/UART setup, matching the
  `en_spreadCycle(true)` / `intpol(false)` pattern already used for
  SpreadCycle-only mode in `ESP32_TMC2209_StepSize_Sweep` and
  `ESP32S3_TMC2209_Chirp_Test` elsewhere in this repo.
- Boot banner changed to `# ESP32_V4_MRES_TRAJECTORY_CAMPAIGN_SPREADCYCLE`
  and a `# DRIVER_MODE,spreadcycle=1,intpol=0` line is printed once the
  driver comes up, so the raw serial capture is self-identifying.
- The `mode` field logged for `CAMPAIGN_START`/`RUN_CONFIG`/`RUN_COMPLETE`/
  `CAMPAIGN_COMPLETE` events changed from `BASELINE` to `SPREADCYCLE` as a
  dataset-level tag (per-block `mode` fields for markers, oscillation, and
  the direct/individual trajectories are unchanged, so the CSV schema and
  block/event counts are identical to the baseline and validate against the
  existing `plot_planned_sequence.py` dry-run checks unchanged).
- Everything else — timing constants, MRES order, rates, dwells,
  separations, marker amplitudes, throughput preflight — is byte-for-byte
  the same as the baseline sketch.

## Comparison notes

Compare the baseline and SpreadCycle datasets only when mechanics, supply
voltage, current, origin, preload, and acquisition settings are unchanged, as
with the StealthChop/StallGuard repeat. Because interpolation is explicitly
disabled here (`intpol(false)`) while the baseline leaves it at its power-on
default (on), a difference in microstep smoothness between the two datasets
is expected and should not be read as a SpreadCycle-vs-StealthChop effect on
its own.

## Build and upload

```sh
arduino-cli compile \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  scripts/esp32_v4_mres_trajectory_campaign_spreadcycle

arduino-cli upload \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  --port /dev/cu.usbmodem101 \
  scripts/esp32_v4_mres_trajectory_campaign_spreadcycle
```

Status: sketch added, not yet compiled, flashed, or run.

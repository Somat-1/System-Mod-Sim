# v6 — MRES trajectory campaign, corrected current and distance

`scripts/esp32_v6_mres_trajectory_campaign` is v4's full MRES trajectory
campaign (`esp32_v4_mres_trajectory_campaign`), forked with two corrections:

1. **`R_SENSE_OHM` is 0.11**, the resistor actually fitted on this board, not
   the 0.03 v4 used (a different board's value). 0.03 under-delivered current
   by roughly 3x (~122 mA actual against the 360 mA requested); 0.11 actually
   delivers close to 360 mA. See
   `../v4/scripts/tmp_v4_manual_range_sweep` and the v5 work for how this was
   found.
2. **`TRAJECTORY_FULL_STEPS` is 10 mm each way**, not v4's 25 mm. The stage's
   actual travel was since measured (`../v4/scripts/tmp_v4_manual_range_sweep`)
   at 44 mm total — a 25 mm one-way leg does not fit at all.

Everything else — pin map, MRES ladder (1, 4, 16, 32), marker scheme, the
three speeds, DIRECT vs INDIVIDUAL command generation, the throughput
preflight — is unchanged from v4. See the
[v4 README](../v4/README.md) for that detail.

## Safety: no range-finding, starts moving at boot

Like v4, this sketch has **no range-finding or centering step** — it starts
moving automatically once the throughput preflight passes. A 10 mm leg each
way only fits safely if the stage already has at least 10 mm of clearance on
both sides of wherever it happens to be sitting when flashed. Confirm that
first (e.g. with `../v4/scripts/tmp_v4_manual_range_sweep`, or by checking
against a previously-established center) rather than assuming it.

If that clearance isn't already known, prefer
`../v5/scripts/esp32_v5_centered_mres_campaign` instead — it establishes and
centers on the actual range interactively (jog + `SETMIN`/`SETMAX`, `RUN`
only once centered) before anything moves automatically, at a shorter 12 mm
leg split into both a positive and a negative excursion per speed/mode,
rather than v6/v4's single one-way leg.

## Build and upload

```sh
arduino-cli compile \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  scripts/esp32_v6_mres_trajectory_campaign

arduino-cli upload \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  --port /dev/cu.usbmodem101 \
  scripts/esp32_v6_mres_trajectory_campaign
```

Upload resets the ESP, runs preflight, and starts the campaign automatically.

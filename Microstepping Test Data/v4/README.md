# v4 — ESP32 MRES oscillation and 25 mm trajectory campaign

## Current run

The live v4 firmware is
`scripts/esp32_v4_mres_trajectory_campaign/esp32_v4_mres_trajectory_campaign.ino`.
It was compiled with USB CDC enabled, flashed to the ESP32-S3, and started
automatically on 2026-09-07.

This is a baseline motion/IDS run. It does **not** configure StealthChop,
SpreadCycle, StallGuard, CoolStep, or interpolation. Only the UART interface,
motor current, and required MRES field are configured. A separate future repeat
with StealthChop and StallGuard is specified in
[STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md](STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md).

## Verified hardware configuration

| Signal | ESP32-S3 connection |
|---|---:|
| STEP | GPIO 5 |
| DIR | GPIO 6 |
| TMC UART RX at ESP | GPIO 18 |
| TMC UART TX from ESP | GPIO 17 |
| Onboard RGB LED | GPIO 48 |
| EN/ENN | Externally grounded |

The TMC2209 returned chip version `0x21` immediately before the live campaign.
The requested run/hold current is 360 mA RMS. Because EN/ENN is externally
grounded, firmware can stop STEP pulses but cannot de-energize the rotor; remove
motor power for a physical emergency stop.

## Measurement sequence

MRES is tested in the order `1, 4, 16, 32`. For every MRES, the firmware runs:

1. A unique configuration marker.
2. A unique oscillation marker.
3. Fifteen complete one-microstep forward/return cycles in 30 seconds: one
   microstep forward, 1 s dwell, one microstep back to origin, 1 s dwell.
4. A 25 mm out-and-back trajectory at slow, moderate, and fast speed using one
   direct distance call per leg. Every speed has its own preceding marker.
5. The same three trajectories using 2,500 individually submitted full-step
   calls per leg. Every speed again has its own preceding marker.

Every 25 mm leg returns to the same origin after a 1 s endpoint dwell. The lead
is 2 mm/revolution and the motor has 200 full steps/revolution, so one full step
is 0.010 mm and 25 mm is 2,500 full steps or 12.5 revolutions.

The selected speeds come from v3:

| Label | Rate (full steps/s) | One-way 25 mm time |
|---|---:|---:|
| slow | 27.5 | 90.91 s |
| moderate | 70 | 35.71 s |
| fast | 200 | 12.50 s |

27.5 full steps/s is the slowest v3 rate that permits all four MRES values,
both generation approaches, markers, and oscillations to fit within the
55-minute recording limit.

### Direct versus individual generation

The TMC2209 itself accepts STEP/DIR pulses, not a millimetre-distance UART
command. The two ESP implementations therefore distinguish where the move is
expanded:

- **Direct:** the firmware receives one 25 mm request for a leg, calculates
  `2,500 × MRES` pulses once, and emits the whole continuous pulse train.
- **Individual:** the firmware submits one full-step movement, waits for its
  `MRES` pulses to finish, and then submits the next full step. This repeats
  2,500 times per leg.

The physical distance and nominal speed are the same. Sending one software
command per microstep at MRES 32 would require 80,000 commands per leg and would
not test the requested one-full-step-at-a-time approach.

## Measured individual-command limit

Before the automatic run, the firmware performed zero-net alternating
full-step commands at MRES 32. A rate passes when measured throughput reaches at
least 95% of the requested rate. The live result was:

- 400 commands/s: passed at 387.011 commands/s;
- 500 commands/s: failed at 468.219 commands/s;
- highest tested reliable command frequency: **400 commands/s**;
- required fast trajectory: 200 commands/s, measured at 196.610 commands/s.

The preflight is logged as `THROUGHPUT_PREFLIGHT_UNRECORDED` and returns to the
origin before `CAMPAIGN_START`.

## Separation and processing markers

Each configuration, oscillation block, command mode, and speed is preceded by
a globally unique negative marker. Marker amplitudes increase from 12 to 136
full steps in four-step increments and return to origin. The serial CSV stream
also records:

- `run_index` and `mres`;
- exact `block`, `mode`, speed, and iteration;
- marker amplitude and ideal position in 1/32-full-step units;
- every oscillation forward/return event;
- trajectory endpoints and achieved individual-command rate.

This provides both visible IDS segmentation and machine-readable event labels.
The onboard LED is white once at boot, yellow during preflight, purple during
markers, green during oscillation, blue during direct trajectories, orange
during individual trajectories, cyan at completion, and red only on failure.

## Duration and render

The complete planned measurement trajectory is **2,477.7 s (41.30 min)**,
including all markers and dwells. This leaves about 13.7 minutes below the
55-minute recording limit.

`scripts/plot_planned_sequence.py` validates a full campaign dry-run and renders
the complete trajectory. The top plot shows all 24 25 mm out-and-back blocks;
four zoom panels show all 15 oscillation cycles at each MRES.

- `rendered_assets/planned_mres_trajectory_campaign.png` — complete preview.
- `rendered_assets/planned_settling_sequence_preview.png` — updated to the same
  complete preview, replacing the old representative-only figure.

## Files

- `scripts/esp32_v4_mres_trajectory_campaign/` — compiled and flashed live ESP
  firmware; the first run starts automatically after preflight.
- `scripts/run_mres_trajectory_campaign.py` — dry-run/host-controller mirror
  used to validate block count, order, timing, and CSV schema without ESP
  motion.
- `scripts/plot_planned_sequence.py` — full-sequence renderer.
- `data/hardware_runs/mres_trajectory_dry_run_*.csv` — validation log.
- `scripts/run_settling_dedicated_controller.py` — retained legacy settling
  experiment; not the current v4 campaign.

## Build and upload

```sh
arduino-cli compile \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  scripts/esp32_v4_mres_trajectory_campaign

arduino-cli upload \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  --port /dev/cu.usbmodem101 \
  scripts/esp32_v4_mres_trajectory_campaign
```

Upload resets the ESP, runs preflight, and starts the campaign automatically.
USB commands remain available: `STATUS`, `ABORT`, and `RUN` to repeat after the
automatic pass. Opening some serial monitors toggles DTR/RTS and may reset the
ESP, which starts a fresh automatic run.

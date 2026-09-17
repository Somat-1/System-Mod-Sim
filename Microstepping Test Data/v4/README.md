# v4 — ESP32 MRES oscillation and 25 mm trajectory campaign

## Current run

The live v4 firmware is
`scripts/esp32_v4_mres_trajectory_campaign/esp32_v4_mres_trajectory_campaign.ino`.
It was compiled with USB CDC enabled, flashed to the ESP32-S3, and started
automatically on 2026-09-07.

The checked-in source has since been revised for the next run: every 25 mm
trajectory now dwells 2 s at +25 mm and 2 s after returning to the origin,
measured experiments have explicit 10 s stationary separations, and each
microstep oscillation block verifies that its net position is unchanged. This
revised source still needs to be compiled and flashed before the next live run.

This is a baseline motion/IDS run. It does **not** configure StallGuard or
CoolStep. The UART interface, motor current, and required MRES field are
configured, and the chopper is forced to **SpreadCycle with MicroPlyer
interpolation disabled**. A separate future repeat with StealthChop and
StallGuard is specified in
[STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md](STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md).
A SpreadCycle-for-the-whole-sequence variant (interpolation off, otherwise
identical to this baseline) is specified in
[SPREADCYCLE_VARIANT_MEMO.md](SPREADCYCLE_VARIANT_MEMO.md) and implemented in
`scripts/esp32_v4_mres_trajectory_campaign_spreadcycle/`. A SpreadCycle +
MicroPlyer (interpolation on) variant, MRES 1/4/16 only, is specified in
[SPREADCYCLE_MICROPLYER_VARIANT_MEMO.md](SPREADCYCLE_MICROPLYER_VARIANT_MEMO.md)
and implemented in
`scripts/esp32_v4_mres134_spreadcycle_microplyer_campaign/`.

### Chopper mode and interpolation are forced, not left alone

Both of the TMC2209's relevant power-on defaults are wrong for this
measurement, so neither may be left unconfigured:

| Register bit | Resets to | Meaning if left alone |
|---|---|---|
| `GCONF` bit 2 `en_spreadCycle` | 0 | **StealthChop** — voltage-mode PWM with automatic amplitude regulation, lower microstep positional fidelity than SpreadCycle |
| `CHOPCONF` bit 28 `intpol` | 1 | **MicroPlyer on** — each commanded microstep expanded into 256 sub-steps and smeared across the interval to the next step |

"Configure nothing to keep the run neutral" therefore selects exactly the two
settings this experiment must not use. Both defects were live in the
2026-09-10 16:59 run, which was verified on the hardware to be running
StealthChop (`GCONF=0x000001C0`, `DRV_STATUS` bit 30 `stealth=1`) with
interpolation active.

Measured with the driver's own `MSCNT` register at MRES=1, a single commanded
full step produced a 13-count move followed by the remaining 243 counts
creeping out over the next second. With `intpol` forced off the same command
produces one clean 256-count move.

`configureDriver()` now refuses to start the campaign unless `GCONF`,
`CHOPCONF` and `DRV_STATUS` all confirm SpreadCycle active and interpolation
off, and logs `CHOPPER_CONFIG` / `CHOPPER_ACTIVE` lines so the mode is recorded
in the data rather than inferred.

**Data recorded before 2026-09-10 ran under StealthChop with interpolation on
and is not comparable with data recorded after.** This affected every block,
not only the oscillations. The partial run in
`data/hardware_runs/mres_trajectory_live_20260910_165920.csv` still carries the
StealthChop configuration and should be treated as superseded.

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
   microstep forward, 1 s dwell, one microstep back to the exact block-start
   position, 1 s dwell. Both host and ESP implementations reject a block whose
   final net position differs from its start.
4. A 25 mm out-and-back trajectory at slow, moderate, and fast speed using one
   direct distance call per leg. The stage dwells 2 s at +25 mm and 2 s after
   returning to the origin. Every speed has its own preceding marker.
5. The same three trajectories using 2,500 individually submitted full-step
   calls per leg, with the same two 2 s endpoint dwells. Every speed again has
   its own preceding marker.

Every new velocity experiment is separated from the preceding measured block
by a 10 s stationary dwell; the unique negative marker follows that gap and
still immediately identifies the experiment that comes next. Adjacent MRES
runs also have a 10 s separation before the next configuration marker.

Every 25 mm out-and-back block returns to the same origin. The lead
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
The additional 10 s zero-position gaps make consecutive measured experiments
visually distinct without replacing or weakening the marker signature.
The onboard LED is white once at boot, yellow during preflight, purple during
markers, green during oscillation, blue during direct trajectories, orange
during individual trajectories, cyan at completion, and red only on failure.

## Duration and render

The complete planned measurement trajectory is **2,795.7 s (46.60 min)**,
including all markers, 2 s endpoint dwells, and 10 s separations. This leaves
about 8.4 minutes below the
55-minute recording limit.

`scripts/plot_planned_sequence.py` validates a full campaign dry-run and renders
two views. In the complete view, the four oscillation panels come first, each
showing coincident start/end markers; the full 24-trajectory plot and overall
timeline are below them. The second view isolates the MRES 4 segment with blue
direct-command and orange individual-command position and signed-velocity
traces.

- `rendered_assets/planned_mres_trajectory_campaign.png` — complete preview.
- `rendered_assets/planned_mres4_velocity_detail.png` — MRES 4 position,
  velocity, and schedule detail.

## Files

- `scripts/esp32_v4_mres_trajectory_campaign/` — compiled and flashed live ESP
  firmware; the first run starts automatically after preflight.
- `scripts/esp32_v4_mres_trajectory_campaign_spreadcycle/` — SpreadCycle
  variant of the same firmware; see
  [SPREADCYCLE_VARIANT_MEMO.md](SPREADCYCLE_VARIANT_MEMO.md).
- `scripts/esp32_v4_mres134_spreadcycle_microplyer_campaign/` — SpreadCycle +
  MicroPlyer variant, MRES 1/4/16 only; see
  [SPREADCYCLE_MICROPLYER_VARIANT_MEMO.md](SPREADCYCLE_MICROPLYER_VARIANT_MEMO.md).
- `scripts/plot_planned_sequence_spreadcycle_microplyer.py` — plan/plot
  renderer for that variant (MRES 1/4/16 instead of 1/4/16/32).
- `scripts/run_mres_trajectory_campaign.py` — dry-run/host-controller mirror
  used to validate block count, order, timing, and CSV schema without ESP
  motion.
- `scripts/dedicated_controller_support.py` — shared guarded transport,
  logging, marker, origin, and shutdown infrastructure used by the active host
  runner.
- `scripts/plot_planned_sequence.py` — full-sequence renderer.
- `data/hardware_runs/mres_trajectory_dry_run_*.csv` — validation log.

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

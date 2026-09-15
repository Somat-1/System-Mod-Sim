# v5 — centered, bidirectional MRES campaign

## What's different from v4

v4's trajectory blocks moved 25 mm one-way from wherever the firmware
happened to boot, with no known safe range (no end-of-travel sensors). That
firmware is superseded here by two things this project has since established:

1. **The actual sense resistor is 0.11 Ω**, not the 0.03 Ω every earlier
   sketch (including v4) used. That earlier value under-delivered current by
   roughly 3x (~122 mA actual against a 360 mA request) — plausibly enough to
   explain degraded low-torque behavior at finer MRES by itself. This sketch
   uses 0.11 Ω and actually delivers close to the requested 360 mA RMS.
2. **The safe travel range has now been measured by hand**, using
   `../v4/scripts/tmp_v4_manual_range_sweep`: **44 mm wide** (min/max 22 mm
   either side of center).

A 25 mm one-way excursion does not fit inside a 44 mm range at all, so v5
does not reuse v4's trajectory distance. Unlike v4 it does not auto-run the
campaign at boot — but it *does* auto-center at boot (see below), which is
itself a real, automatic move worth understanding before flashing.

## The range is reused, not re-measured — this is an assumption

There is still no homing sensor, and the position counter always resets to 0
at boot while the physical stage does not. Rather than require re-jogging
both limits from scratch every session, this firmware hard-codes the range
last measured with `../v4/scripts/tmp_v4_manual_range_sweep`
(`REUSED_RANGE_MIN_FULL_STEPS = 2050`, `REUSED_RANGE_MAX_FULL_STEPS = 6450`,
44 mm wide) and assumes **boot position 0 is physically the same spot as that
session's confirmed max** (the stage was sitting at max when that range was
last read back). Flashing doesn't move the motor, so this holds as long as
nothing jogged the stage, on any firmware, between that session and this
boot — but it is an assumption, not a measurement, and it fails silently if
that's not true.

Because of that, `setup()`:

1. Applies the reused range and prints it (`# REUSED_RANGE_ASSUMED,...` plus
   a `# WARNING` line) immediately after driver configuration.
2. Gives a **5-second, `ABORT`-able countdown** (`# AUTO_CENTER_STARTING,...`)
   before moving.
3. Drives to the computed center automatically, at MRES 32 for exact
   reachability, then redefines position 0 as that center.

Watch the terminal (and the stage) right after flashing — that countdown is
the only window to catch a wrong assumption before the stage moves on it. If
anything looks off, send `ABORT`; the firmware falls back to idle with the
full manual toolkit (`J+`/`J-`/`NUDGE+`/`NUDGE-`/`SETMIN`/`SETMAX`/
`CLEARRANGE`) available to re-establish the range by hand, same as
`tmp_v4_manual_range_sweep`.

## Sequence

Auto-centering happens once at boot, before any commands are needed. Once
`RUN` is sent:

1. **Throughput preflight** — same alternating full-step command-rate test as
   v4, at MRES 32, before anything else moves.
2. **Re-center** — a no-op unless the stage was jogged since boot; otherwise
   this is where centering actually happens, at MRES 32 for exact
   reachability regardless of where a limit was last set. Position 0 is
   redefined here if it moves; every block after this returns to this same
   origin, exactly like v4's `CAMPAIGN_START` origin.
3. For each MRES in **1, 4, 16, 32**: a config marker, an oscillation marker,
   the unchanged 15-cycle one-microstep oscillation, then trajectory blocks
   for **both directions** (`POS`, `NEG`) × both command modes (`DIRECT`,
   `INDIVIDUAL`) × the same three v3/v4 speeds (slow/moderate/fast) — twelve
   trajectory blocks per MRES instead of v4's six, each preceded by its own
   marker.

## Trajectory distance and margin

`TRAJECTORY_FULL_STEPS = 1200` (12 mm) each way, against a measured 22 mm
half-range: **10 mm of physical clearance to either hard stop (~45%
margin)**. `SETMIN`/`SETMAX` sets `# WARNING,TRAJECTORY_FULL_STEPS_DOES_NOT_FIT_IN_RANGE`
if the range turns out narrower than that on a given session — check for it
before sending `RUN`. Edit the constant and reflash for a different distance.

## Duration

Estimated **~42 minutes** total (4 MRES × 2 directions × 2 modes × 3 speeds
of 12 mm legs, plus markers and oscillations), under the 55-minute budget v4
used.

## Commands

Same jog/range vocabulary as `tmp_v4_manual_range_sweep`, plus `RUN`. The
range starts pre-set (reused) and centered by the time these are usable:

- `J+ [n]` / `J- [n]` — jog n full steps (default 5), clamped to the range.
- `NUDGE+` / `NUDGE-` — single microstep, for the final approach to a stop.
- `SETMIN` / `SETMAX` — overwrite the reused limit with the current position.
- `CLEARRANGE` — forget the range entirely (not just override one side).
- `RUN` — preflight, re-center (no-op unless jogged since boot), then run the
  full MRES campaign. Refused until a range is set.
- `PAUSE` / `RESUME` — pause any in-progress motion or run; resume it.
- `ABORT` — stop now, not resumable. Also cancels the boot-time auto-center
  countdown/move.
- `STATUS` / `HELP`.

Only `ABORT`/`PAUSE`/`RESUME`/`STATUS`/`HELP` are accepted while busy.

## Build and upload

```sh
arduino-cli compile \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  scripts/esp32_v5_centered_mres_campaign

arduino-cli upload \
  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  --port /dev/cu.usbmodem101 \
  scripts/esp32_v5_centered_mres_campaign
```

Upload resets the ESP. **It then auto-centers within ~10-15 seconds** (driver
config, then a 5-second `ABORT`-able countdown, then the move) based on the
reused range above — watch the terminal and the stage. After that it idles,
waiting for jog/range-override/`RUN` commands over USB CDC serial
(115200 baud).

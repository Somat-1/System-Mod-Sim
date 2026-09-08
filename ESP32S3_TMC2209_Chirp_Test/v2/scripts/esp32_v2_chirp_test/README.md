# esp32_v2_chirp_test firmware

**Unflashed, untested.** A fresh implementation (not a port of v1) that
runs the cruise-amplitude + smooth-notch chirp designed in
`../chirp_v2_schedule.py` / `../../README.md` / `../../SWEEP_CONFIG_LOG.md`.
Two things must happen before this touches real hardware:

1. **Confirm the pins.** `STEP_PIN`, `DIR_PIN`, `EN_PIN`, the UART pins, and
   the trigger pins are placeholders reused from v1's documented wiring —
   not re-verified for this file. Check against the actual rig.
2. **Bench-test the achievable STEP pulse rate** before trusting
   `STEP_RATE_CLAMP_HZ` (200 kHz). v1's original `digitalWrite()` busy-wait
   approach implied only ~111 kHz was reliably achievable; this firmware
   switches the hot path to direct `GPIO.out_w1ts`/`out_w1tc` register
   writes specifically to buy headroom above that, but that is reasoning,
   not a measurement. A simple standalone test — drive `STEP_PIN` alone,
   alternating `DIR_PIN` occasionally, and count edges over a timed window
   — would give a real number, the same way the v4 campaign benchmarks
   throughput before trusting a plan
   (`benchmark_individual_rate()` in
   `../../../Microstepping Test Data/v4/scripts/run_mres_trajectory_campaign.py`).

I was also unable to compile-check this file (no `arduino-cli` available in
this environment) or verify a few TMCStepper API calls
(`hstrt()`/`hend()`/`tbl()`) against the installed library version from
memory alone — both are flagged with `NOT COMPILE-VERIFIED` comments at
their call sites in the `.ino`. Compile it before flashing, obviously, but
treat those specific lines as the first place to look if it fails to build.

## What's different from v1

- **MRES = 16** (register code 4), not native 1/256 — see the top-level
  README's "Why the notch comes out, and why MRES changed."
- **Position-tracking control loop**, not a fixed ±1 microstep alternation.
  `commandedStateAt(elapsed_s)` analytically evaluates the exact ideal
  microstep target for any elapsed time (closed-form phase integration, so
  no accumulated drift even if a loop iteration is late), covering idle,
  both sync markers, both sweep directions, the mid-dwell, and the tail as
  one continuous function. The main loop just steps the real position
  toward whatever that function currently says, one microstep at a time,
  rate-limited to `MIN_STEP_PERIOD_US` (3 µs, ~333 kHz physical ceiling)
  regardless of how fast the loop itself iterates.
- **Direct GPIO register writes** (`GPIO.out_w1ts`/`out_w1tc`) for the STEP
  pulse specifically, instead of `digitalWrite()`, to remove the timing
  ceiling that approach implied (see point 2 above). DIR still uses
  `digitalWrite()` since it only changes at direction reversals, not every
  step, so it is not on the hot path.
- **Smooth notch**, not a hard exclusion band — STEP commands are never
  suppressed; the amplitude is tapered by the raised-cosine window from
  `chirp_v2_schedule.resonance_notch_taper()`.
- A defensive sanity ceiling (`MAX_ABS_TARGET_MICROSTEPS`) aborts the run if
  the computed target ever exceeds 1.2x cruise amplitude, as a guard
  against a firmware bug commanding something the design never intended —
  independent of the two clamps already baked into `amplitudeMicrosteps()`.

## Constants must be kept in sync by hand

Firmware cannot import `chirp_v2_schedule.py`. Every constant in the
"Mirrors chirp_v2_schedule.py" block at the top of the `.ino` is a
hand-copied snapshot — re-check it against `../../SWEEP_CONFIG_LOG.md`
(regenerate that first if the schedule module has changed) before every
flash, not just the first one.

## Serial protocol

115200 baud. `CHECK` configures the driver and prints the plan without
moving. `RUN` executes one complete level (idle → marker → up-sweep →
dwell → marker → down-sweep → tail, ~9.1 min) and returns to the origin
before disabling the motor, same as on `ABORT`. Log lines are
`timestamp_us,event,segment,frequency_hz,position_microsteps` CSV, written
only at segment transitions (not every microstep — that would collide with
the timing this file exists to protect).

## Build and upload

```powershell
arduino-cli compile --fqbn esp32:esp32:esp32s3 .\esp32_v2_chirp_test
arduino-cli upload --fqbn esp32:esp32:esp32s3 --port COM_PORT .\esp32_v2_chirp_test
```

Requirements: Arduino-ESP32 3.x, TMCStepper. Replace `COM_PORT`; open the
serial console at 115200 baud afterward.

Perform the first run with the mechanism unloaded and observed directly,
exactly as v1's own README insists — nothing about a smoother notch or a
coarser MRES changes that.

# esp32_v2_chirp_test firmware

**Validated on real hardware.** Runs the cruise-amplitude + smooth-notch
chirp designed in `../chirp_v2_schedule.py` / `../../README.md` /
`../../SWEEP_CONFIG_LOG.md`, using dual-core step generation. A full `RUN`
(pre-roll, both sync markers, both sweep directions, mid-dwell, tail)
completed segment-by-segment on schedule with no watchdog reset, and the
TMC2209 connected and configured correctly over the real UART link
(`test_connection()` returns 0, MRES readback matches).

## Why dual-core

The original single-core design (computing `sin()`/notch/clamp math and
issuing the STEP pulse in the same loop iteration) was bench-tested on this
exact board (scratch-pin tests, no TMC2209 involved, not committed to the
repo) and could not sustain anywhere near `STEP_RATE_CLAMP_HZ` (200 kHz):

1. **Polled single-core loop**: real edge rate flatlined around ~55 kHz
   regardless of demand, because the trig/notch/clamp math ran inline with
   GPIO issuance.
2. **Single-core hardware-timer ISR**: decoupled the math from edge
   issuance, but sharing one core meant a firing rate high enough to help
   throughput instead starved the main loop's target computation
   (main-loop rate collapsed 5x) — trading spatial staircasing (too few
   microsteps) for temporal staircasing (too few position updates per
   cycle).
3. **Dual-core (this file)**: core 0 runs a dedicated step-issuing task
   (`stepTaskFn`) that only ever compares two integers and writes a GPIO
   register — no floating point anywhere near it. Core 1 runs the normal
   Arduino `setup()`/`loop()`, doing all the `sin()`/notch/clamp math and
   publishing a target microstep position into a shared `volatile`. Neither
   starves the other. Verified against the mathematically correct
   reference (a perfectly-tracked sine wave's time-averaged edge rate is
   exactly 2/π of its peak rate) that this reproduces the intended waveform
   accurately across the full 1-1000 Hz sweep at the design's original
   MRES=16 — MRES did not need to be lowered to make this work.

## Pin corrections from the original design

The pins in this file are confirmed against the working v4 MRES trajectory
campaign sketch's wiring for this specific rig, not reused from v1 as
placeholders:

- **STEP=GPIO5, DIR=GPIO6, UART_RX=GPIO18, UART_TX=GPIO17.** The original
  draft had STEP=6/DIR=7/EN=5, which did not match this rig at all.
- **No `EN_PIN`.** EN/ENN is physically grounded on this rig (always
  enabled), matching v4 — there is no GPIO controlling it, and no
  software-controlled disable between runs.
- **No `TRIG_OUT`/`TRIG_ECHO`.** The confirmed-correct v4 wiring for this
  rig has no such pins, so they were almost certainly never wired here
  either. If this rig does get a DAQ sync line, add it back deliberately
  with a confirmed pin.

## A real firmware bug found and fixed on hardware

The first full-sequence `RUN` attempt **crashed and rebooted the chip ~15 s
into `PRE_ROLL`**: `stepTaskFn` never yields on core 0, and real chirp
segments include long silent stretches (30 s pre-roll/tail, 5 s mid-dwell)
where the task is "running" but has nothing to step — starving core 0's
IDLE task long enough to trip the ESP-IDF task watchdog. This was invisible
in all prior bench-scale tests (1-1.5 s each), which never held that state
long enough to trigger it.

`disableCore0WDT()` was tried first and did not fully work on this IDF
version — it removes IDLE0 from the watchdog's registry, but IDLE0's own
idle hook still unconditionally calls `esp_task_wdt_reset()` afterward,
producing a continuous stream of `task not found` errors instead of a
crash. Since core 0 is permanently, deliberately dedicated to step
generation by design, the fix is `esp_task_wdt_deinit()` in `setup()` —
deiniting the whole Task Watchdog Timer subsystem rather than fighting a
partial per-task removal.

A second, unrelated usability bug was also found and fixed: typing
commands into a raw serial terminal (e.g. `screen`) at normal human speed
got truncated into unrecognized single characters, because
`Serial.readStringUntil('\n')` only waited 20 ms per call. `readSerialLine()`
now accumulates characters across as many calls as it takes, with no
timeout, so it works regardless of typing speed.

## What's different from v1

- **MRES = 16** (register code 4), not native 1/256 — see the top-level
  README's "Why the notch comes out, and why MRES changed."
- **Position-tracking control loop**, not a fixed ±1 microstep alternation.
  `commandedStateAt(elapsed_s)` analytically evaluates the exact ideal
  microstep target for any elapsed time (closed-form phase integration, so
  no accumulated drift even if a loop iteration is late), covering idle,
  both sync markers, both sweep directions, the mid-dwell, and the tail as
  one continuous function. Core 1 evaluates this and publishes it; core 0
  continuously steps the real position toward whatever the shared target
  currently is, rate-limited to `MIN_STEP_PERIOD_US` (3 µs) — benchmarking
  showed this is not the binding constraint once split across both cores.
- **Direct GPIO register writes** (`GPIO.out_w1ts`/`out_w1tc`) for the STEP
  pulse specifically, instead of `digitalWrite()`. DIR still uses
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
before leaving the driver enabled (EN is hardwired, not software-disabled),
same as on `ABORT`. Log lines are
`timestamp_us,event,segment,frequency_hz,position_microsteps` CSV, written
only at segment transitions (not every microstep — that would collide with
the timing this file exists to protect).

## Build and upload

```
arduino-cli compile --fqbn esp32:esp32:esp32s3:CDCOnBoot=cdc esp32_v2_chirp_test
arduino-cli upload --fqbn esp32:esp32:esp32s3:CDCOnBoot=cdc --port PORT esp32_v2_chirp_test
```

Requirements: Arduino-ESP32 3.x, TMCStepper. Replace `PORT`. The
`CDCOnBoot=cdc` board option is required on ESP32-S3 for `Serial` to route
over the native USB port instead of UART0 — without it, nothing appears on
the USB serial console at all. Open the serial console at 115200 baud
afterward; typing is line-buffered with no timeout, so normal typing speed
works.

`f_n=176.7 Hz` and `Q=20` are still the prior analytical-model
placeholders, not measured values — the notch is centered on a model, not
on this mechanism's actual resonance. Replace both with same-day ringdown
measurements (top-level README Section 1) before treating a run as
scientifically meaningful.

Perform the first run with the mechanism unloaded and observed directly,
exactly as v1's own README insists — nothing about a smoother notch, a
different pin set, or dual-core step generation changes that.

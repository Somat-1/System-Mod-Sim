# What was commanded to the motor — 2026-10-13 rotary-encoder runs

Command-side reference for the two encoder recordings in [`data/`](data/). It
documents **what the TMC2209 was told to do**: the register state in force, the
STEP/DIR waveform of every action, and the exact order in which the actions
ran. It deliberately says almost nothing about how the rotor responded — that
is the subject of
[`analysis/motor_encoder_analysis_report.pdf`](analysis/motor_encoder_analysis_report.pdf).

Reconstructed from the executed plans (`data/*/plan.csv`), the firmware
[`firmware/motor_identification/motor_identification.ino`](firmware/motor_identification/motor_identification.ino),
the host runner [`scripts/run_campaign.py`](scripts/run_campaign.py), and the
`#` header and boundary lines inside the recordings themselves.

## The two runs

| | Run A — identification | Run B — microstep resolution |
|---|---|---|
| Data | [`data/2026-10-13_motor_identification/`](data/2026-10-13_motor_identification/) | [`data/2026-10-13_microstep_resolution/`](data/2026-10-13_microstep_resolution/) |
| Host start (UTC) | 2026-10-13 13:50:05 | 2026-10-13 14:47:22 |
| Actions | 83 | 13 |
| Planned duration | 37.10 min | 4.75 min |
| STEP edges issued | 339,840 | 89,660 |
| Telemetry rows | 1,903,061 | 650,110 |
| Terminated | `# RUN_COMPLETE` | `# RUN_COMPLETE` |

A third recording was made at 15:45 local the same afternoon with firmware
lacking the PCNT watch-point fix; its encoder count restarted every 30,000
counts and it was discarded. It is not in `data/`.

## Units used throughout

| Quantity | Value |
|---|---|
| Motor | 200 full steps/rev, 50 electrical cycles/rev (4 full steps = 1 electrical cycle) |
| Encoder | 20,000 counts/rev (5,000 lines × 4), index pulse once per rev |
| Scale | 1 full step = 100 counts; 1 count = 0.018° = 64.8 arcsec |
| 1/16 microstep | 16/256 full step = 6.25 counts |
| `command_u256` | commanded phase in 1/256 full step; 256 = 1 full step. Invariant across MRES changes |

`command_u256` is pure bookkeeping: the firmware increments it by `256/MRES`
per STEP edge. It is what was *asked for*, never what was measured.

## Signal path

The ESP32-S3 drove the TMC2209 two ways, and the distinction matters when
reading the data:

1. **STEP/DIR GPIO** — carried all motion. Every commanded increment is one
   rising edge on STEP with DIR held at the intended level.
2. **UART (single-wire, 115200 baud, address 0)** — carried only
   *configuration*: microstep resolution (`CHOPCONF.MRES`) and current
   (`IHOLD_IRUN`), plus register readback for verification. No motion was ever
   commanded over UART.

Bench wiring as recorded in the run header (`# PINS`):

| Signal | GPIO |
|---|---|
| STEP | 7 |
| DIR | 8 |
| EN (TMC `ENN`, active low) | 9 |
| UART TX / RX | 17 / 18 |
| Encoder A / B / Z | 4 / 5 / 6 |

`R_sense` = 0.10 Ω. `EN` was held HIGH (driver disabled) from reset until
register readback confirmed the required mode, then pulled LOW and **left LOW
for the entire run** — the motor was energised continuously, including during
every standstill.

### STEP edge timing

Produced by `pulseOne()` in the firmware, identical for every action:

- DIR is written only when the direction differs from the previous edge,
  followed by a **20 µs** setup delay before the STEP edge.
- STEP is driven HIGH for **3 µs**, then LOW.
- `CHOPCONF.dedge = 0`, so only the **rising** edge counts: one pulse = one
  microstep.
- DIR HIGH = positive command = increasing encoder count (polarity confirmed
  on this bench).
- Pulse rates were capped at 20 kHz; the fastest actually used was 6.41 kHz.

Timing came from busy-wait loops on `micros()` interleaved with telemetry
emission, not from a hardware timer or the RMT peripheral, so intervals carry
software jitter. Pulse periods are integer microseconds, truncated from
`1e9 / (speed_milli_fs × MRES)` — e.g. 200 full-step/s at 1/16 gives 312 µs
(3205.1 Hz) rather than an exact 3200 Hz.

## TMC2209 register state

The firmware refused to run unless SpreadCycle, interpolation-off, inactive
StealthChop, and `IHOLD == IRUN` were all confirmed by readback, and it
re-verified before **every** action (`verifyDriverMode`), aborting the run on
any change. The readback line recorded at the start of Run A:

```
# DRIVER,gconf=0x000001C4,chopconf=0x04028055,drv_status=0x80090000,
  spreadcycle=1,intpol=0,stealth=0,ihold=9,irun=9,coolstep=0
```

Decoded:

| Register | Field | Value | Meaning |
|---|---|---|---|
| `GCONF` | `I_scale_analog` | 0 | current set digitally; VREF pin ignored |
| | `internal_Rsense` | 0 | external 0.10 Ω sense resistors |
| | `en_spreadCycle` | 1 | **SpreadCycle chopper, always** |
| | `shaft` | 0 | direction not inverted inside the driver |
| | `pdn_disable` | 1 | PDN_UART used as UART; standstill power-down off |
| | `mstep_reg_select` | 1 | **MRES from UART, not the MS1/MS2 pins** |
| | `multistep_filt` | 1 | step-pulse input filter on (library default, not a deliberate choice) |
| `CHOPCONF` | `TOFF` | 5 | chopper slow-decay time |
| | `HSTRT` / `HEND` | 5 / 0 | hysteresis start / end |
| | `TBL` | 1 | blank time 24 clocks |
| | `VSENSE` | 1 | 0.180 V full-scale sense — sets the current scale (below) |
| | `MRES` | 4 | 1/16 microstep at the moment of readback |
| | `dedge` | 0 | rising STEP edge only |
| | `INTPOL` | 0 | **MicroPlyer interpolation disabled** |
| `DRV_STATUS` | `CS_ACTUAL` | 9 | effective current scale |
| | `stealth` | 0 | **StealthChop inactive** |
| | `stst` | 1 | standstill at readback (expected) |

Written at configuration and never changed during a run:

| Setting | Value | Effect |
|---|---|---|
| `TPWMTHRS` | 0 | never switches to StealthChop at any velocity |
| `TCOOLTHRS` / `SEMIN` | 0 / 0 | CoolStep and StallGuard disabled |
| `TPOWERDOWN` / `IHOLDDELAY` | 0 / 0 | no standstill current reduction |
| `IHOLD` | = `IRUN` | hold current equals run current |
| `pwm_autoscale` | false | StealthChop PWM autoscaling off (moot) |

So across both runs the chopper, interpolation, and current-regulation
behaviour were **constant**. The only driver parameters that ever changed
mid-run were `CHOPCONF.MRES` and `IHOLD_IRUN`.

### What the requested currents actually became

`IRUN`/`IHOLD` are 5-bit codes, so requested milliamps quantise. With
`VSENSE = 1` and `R_sense = 0.10 Ω`:

> I_RMS = ((CS+1)/32) × 0.180 V / (0.10 + 0.02) Ω / √2 = (CS+1) × 33.15 mA

One CS code is **33 mA**, which is coarse relative to the sweep steps.

| Requested | CS | I_RMS | Fraction of nominal |
|--:|--:|--:|--:|
| 360 mA | 9 (read back) | 331 mA | 1.00 |
| 270 mA | 7 (computed) | 265 mA | 0.80 |
| 180 mA | 4 (computed) | 166 mA | 0.50 |

Only `CS = 9` is confirmed — it appears in the `# DRIVER` readback. The other
two are computed from the same formula, assuming the truncating integer
conversion in `TMCStepper::rms_current()`. The firmware logs the *requested*
milliamps in the `current_ma` column and never records `CS_ACTUAL` per action.

**This matters for the current sweep.** It was designed as 50 / 75 / 100 % of
nominal, but the realised levels are closer to **0.50 / 0.80 / 1.00**. Any fit
separating detent torque (angle error scaling roughly as 1/I) from
phase-current waveform error (roughly I-independent) should use the realised
currents, not the labels. Logging `CS_ACTUAL` per action would remove the
assumption entirely.

## Command primitives

Every action is one of six kinds, dispatched by `executeAction()`. Before each
one the firmware applies that action's current and MRES, then re-verifies the
mode. The pattern column in the tables below follows directly from these
definitions.

| Kind | `kind` | STEP/DIR pattern |
|---|:-:|---|
| `hold` | 0 | No edges at all. DIR untouched. Energised standstill for `dwell_ms`, logged at `sample_hz`. |
| `move` | 1 | `abs(microsteps)` edges at a uniform period, one direction, then a `dwell_ms` settle. Positioning, not measurement. |
| `stair` | 2 | One edge, then `dwell_ms` of settling — repeated `abs(microsteps)` times. The settled-position primitive. |
| `step_response` | 3 | `repeats` × (excite `abs(microsteps)` edges → `dwell_ms` → return the same edges the other way → `dwell_ms`). Net zero. |
| `velocity` | 4 | Identical code path to `move`: constant-period train, then dwell. |
| `reversal` | 5 | `repeats` × (+A, −2A, +A edges) at constant period — a triangle centred on the start point, peak-to-peak 2A. Net zero. |

Three consequences worth keeping in mind when reading the data:

- **There is no acceleration ramp anywhere.** `move`, `velocity`, and
  `reversal` begin pulsing at their final rate from the first edge and stop
  dead at the last. The 400 full-step/s action starts at 6.41 kHz from
  standstill.
- **`move` and `velocity` are the same firmware branch.** The distinction is
  intent only — positioning vs. measurement — not behaviour.
- **Edge marking is selective.** `stair`, `step_response`, and `reversal` emit
  an extra telemetry row stamped `event = +1/−1` at each edge; `move` and
  `velocity` do not. Run A contains exactly **12,736** marked edges — 6,368
  positive and 6,368 negative — out of 339,840 issued. The other 327,104
  (warm-up, velocity plateaus, dynamics repositioning) are reconstructable only
  from `command_u256`.

## Run startup, in order

The host deliberately does not reset the board — native USB CDC may or may not
reset on open — so it queries live state instead:

1. `PING` → `# PONG`, requires `ready=1`.
2. `DRIVER` → `# DRIVER,...`, must contain `spreadcycle=1`, `intpol=0`,
   `stealth=0`, else abort.
3. `PLAN CLEAR` → `# PLAN_CLEARED`.
4. 83 (or 13) × `PLAN ADD <kind> <microsteps> <speed_milli_fs> <dwell_ms>
   <repeats> <current_ma> <mres> <sample_hz> <tag>`, each acknowledged with
   `# PLAN_ADD`. The whole plan is resident on the board before any motion.
5. Operator confirmation that the motor is decoupled or known-load, the encoder
   is on the shaft, and rotation is safe.
6. `ZERO` → `# ZERO_OK`: sets `command_u256 = 0` and clears the encoder
   counter. **Both runs therefore begin at encoder 0 by definition, not by
   measurement.**
7. `PLAN RUN`. The board streams `# ACTION_BEGIN`, `D,...` rows, and
   `# ACTION_END` per action, finishing with `# RUN_COMPLETE`.

`Ctrl+C` on the host sends `ABORT`; neither run needed it.

---

# Run A — identification campaign

83 actions, 37.10 min planned, 339,840 STEP edges, 1,903,061 telemetry rows.
The plan is built and validated by
[`scripts/sequence_design.py`](scripts/sequence_design.py), which enforces that
the campaign returns to `command_u256 = 0`, that the static map is exactly one
revolution out and back, and that no pulse rate exceeds 20 kHz.

## Blocks

| block | actions | STEP edges | duration | MRES used | I req (mA) |
|---|--:|--:|--:|---|---|
| `warmup` | 13 | 230,400 | 4.10 min | 1/16 | 360 |
| `alignment` | 2 | 128 | 0.26 min | 1/16 | 360 |
| `static_map` | 2 | 6,400 | 16.00 min | 1/16 | 360 |
| `current_sweep` | 9 | 1,536 | 3.27 min | 1/16 | 180, 270, 360 |
| `dynamics` | 20 | 5,184 | 2.43 min | 1/1, 1/16 | 360 |
| `velocity` | 22 | 91,904 | 4.90 min | 1/16 | 360 |
| `reversal` | 12 | 4,032 | 5.00 min | 1/16 | 360 |
| `close` | 3 | 256 | 1.14 min | 1/16 | 360 |
| **total** | **83** | **339,840** | **37.10 min** | | |

## What each block commanded, and why

**`warmup`** — one 60 s energised standstill (the opening reference), then six
out-and-back pairs of six shaft revolutions each at 80 full-step/s. 230,400
edges, 68 % of the campaign's total, purely to bring the motor and driver to
thermal steady state and to run in the bearings. No edges are marked; this
block is not intended for position analysis.

**`alignment`** — one electrical cycle (64 microsteps at 1/16 = 4 full steps)
out and back as a single-microstep staircase with 120 ms dwells. Establishes
the electrical-phase reference against which the static map is read.

**`static_map`** — the core measurement. Every one of the 3,200 native 1/16
table entries over one complete shaft revolution, one microstep at a time with
a **150 ms settling dwell** at each, then the identical 3,200 positions in
reverse. 16.00 min, 6,400 marked edges. The forward/reverse pair is what makes
the separation possible: their mean is the repeatable static angle error, half
their difference is directional hysteresis.

**`current_sweep`** — at each of three current levels, a 4 s settle followed by
four electrical cycles (256 microsteps) out and back with 120 ms dwells. The
only block where `IHOLD_IRUN` changes. Intended to scale detent torque against
waveform error — see the current-quantisation caveat above.

**`dynamics`** — single-increment step responses at four rotor angles (0, 50,
100, 150 full steps), reached by unmarked `move` actions at 40 full-step/s. At
each angle, four actions: one full step (`MRES = 1`) positive and negative,
then one 1/16 microstep positive and negative, each 12 times out and back with
350 ms either side. The excitation and reset edges are both marked, so each is
a usable natural-frequency and damping record. A final `move` of −2,400
microsteps returns to zero.

**`velocity`** — eleven log-spaced plateaus from 0.25 to 400 full-step/s, each
immediately reversed. Span grows with speed (8 to 600 full steps) so each
plateau lasts long enough to measure, capped at three revolutions. 1,500 ms
dwell after each. Unmarked edges; analyse only the constant-rate interior and
discard the first and last electrical cycle.

**`reversal`** — triangular reversals at amplitudes of 1, 2, 4, 8, 16, and 32
native 1/16 microsteps, eight cycles each, at 0.2 Hz and then 1 Hz. Net zero by
construction. Targets pre-sliding friction and backlash-like behaviour at the
smallest commandable amplitudes.

**`close`** — two electrical cycles (128 microsteps) out and back as a 150 ms
staircase, repeating the opening geometry, then a 30 s energised standstill.
Together with the `warmup` hold these bracket the run for thermal drift and
repeatability.

## Every action, in order

| # | block | label | MRES | I req (mA) | STEP edges | DIR | commanded STEP/DIR pattern | cmd span (full steps) |
|--:|---|---|:-:|--:|--:|:-:|---|---|
| 0 | `warmup` | energized hold | 1/16 | 360 | 0 | — | hold 60 s, no STEP edges | 0 |
| 1 | `warmup` | warmup + rep 1 | 1/16 | 360 | 19,200 | + | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 0 → 1200 |
| 2 | `warmup` | warmup - rep 1 | 1/16 | 360 | 19,200 | − | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 1200 → 0 |
| 3 | `warmup` | warmup + rep 2 | 1/16 | 360 | 19,200 | + | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 0 → 1200 |
| 4 | `warmup` | warmup - rep 2 | 1/16 | 360 | 19,200 | − | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 1200 → 0 |
| 5 | `warmup` | warmup + rep 3 | 1/16 | 360 | 19,200 | + | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 0 → 1200 |
| 6 | `warmup` | warmup - rep 3 | 1/16 | 360 | 19,200 | − | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 1200 → 0 |
| 7 | `warmup` | warmup + rep 4 | 1/16 | 360 | 19,200 | + | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 0 → 1200 |
| 8 | `warmup` | warmup - rep 4 | 1/16 | 360 | 19,200 | − | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 1200 → 0 |
| 9 | `warmup` | warmup + rep 5 | 1/16 | 360 | 19,200 | + | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 0 → 1200 |
| 10 | `warmup` | warmup - rep 5 | 1/16 | 360 | 19,200 | − | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 1200 → 0 |
| 11 | `warmup` | warmup + rep 6 | 1/16 | 360 | 19,200 | + | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 0 → 1200 |
| 12 | `warmup` | warmup - rep 6 | 1/16 | 360 | 19,200 | − | 19200 edges, uniform 781 µs (1280.4 Hz), no ramp; then 0.5 s dwell | 1200 → 0 |
| 13 | `alignment` | electrical cycle + | 1/16 | 360 | 64 | + | 64 × (1 edge, then 120 ms dwell) | 0 → 4 |
| 14 | `alignment` | electrical cycle - | 1/16 | 360 | 64 | − | 64 × (1 edge, then 120 ms dwell) | 4 → 0 |
| 15 | `static_map` | full revolution + | 1/16 | 360 | 3,200 | + | 3200 × (1 edge, then 150 ms dwell) | 0 → 200 |
| 16 | `static_map` | full revolution - | 1/16 | 360 | 3,200 | − | 3200 × (1 edge, then 150 ms dwell) | 200 → 0 |
| 17 | `current_sweep` | settle at 180 mA | 1/16 | 180 | 0 | — | hold 4 s, no STEP edges | 0 |
| 18 | `current_sweep` | 180 mA + | 1/16 | 180 | 256 | + | 256 × (1 edge, then 120 ms dwell) | 0 → 16 |
| 19 | `current_sweep` | 180 mA - | 1/16 | 180 | 256 | − | 256 × (1 edge, then 120 ms dwell) | 16 → 0 |
| 20 | `current_sweep` | settle at 270 mA | 1/16 | 270 | 0 | — | hold 4 s, no STEP edges | 0 |
| 21 | `current_sweep` | 270 mA + | 1/16 | 270 | 256 | + | 256 × (1 edge, then 120 ms dwell) | 0 → 16 |
| 22 | `current_sweep` | 270 mA - | 1/16 | 270 | 256 | − | 256 × (1 edge, then 120 ms dwell) | 16 → 0 |
| 23 | `current_sweep` | settle at 360 mA | 1/16 | 360 | 0 | — | hold 4 s, no STEP edges | 0 |
| 24 | `current_sweep` | 360 mA + | 1/16 | 360 | 256 | + | 256 × (1 edge, then 120 ms dwell) | 0 → 16 |
| 25 | `current_sweep` | 360 mA - | 1/16 | 360 | 256 | − | 256 × (1 edge, then 120 ms dwell) | 16 → 0 |
| 26 | `dynamics` | full step + at 0 fs | 1/1 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 0 |
| 27 | `dynamics` | full step - at 0 fs | 1/1 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 0 |
| 28 | `dynamics` | one microstep + at 0 fs | 1/16 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 0 |
| 29 | `dynamics` | one microstep - at 0 fs | 1/16 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 0 |
| 30 | `dynamics` | move to 50 fs | 1/16 | 360 | 800 | + | 800 edges, uniform 1562 µs (640.2 Hz), no ramp; then 1 s dwell | 0 → 50 |
| 31 | `dynamics` | full step + at 50 fs | 1/1 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 50 |
| 32 | `dynamics` | full step - at 50 fs | 1/1 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 50 |
| 33 | `dynamics` | one microstep + at 50 fs | 1/16 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 50 |
| 34 | `dynamics` | one microstep - at 50 fs | 1/16 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 50 |
| 35 | `dynamics` | move to 100 fs | 1/16 | 360 | 800 | + | 800 edges, uniform 1562 µs (640.2 Hz), no ramp; then 1 s dwell | 50 → 100 |
| 36 | `dynamics` | full step + at 100 fs | 1/1 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 100 |
| 37 | `dynamics` | full step - at 100 fs | 1/1 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 100 |
| 38 | `dynamics` | one microstep + at 100 fs | 1/16 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 100 |
| 39 | `dynamics` | one microstep - at 100 fs | 1/16 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 100 |
| 40 | `dynamics` | move to 150 fs | 1/16 | 360 | 800 | + | 800 edges, uniform 1562 µs (640.2 Hz), no ramp; then 1 s dwell | 100 → 150 |
| 41 | `dynamics` | full step + at 150 fs | 1/1 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 150 |
| 42 | `dynamics` | full step - at 150 fs | 1/1 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 150 |
| 43 | `dynamics` | one microstep + at 150 fs | 1/16 | 360 | 24 | ± | 12 × (+1 edge, 350 ms, −1 edge, 350 ms) | 150 |
| 44 | `dynamics` | one microstep - at 150 fs | 1/16 | 360 | 24 | ± | 12 × (−1 edge, 350 ms, +1 edge, 350 ms) | 150 |
| 45 | `dynamics` | return from dynamics | 1/16 | 360 | 2,400 | − | 2400 edges, uniform 1562 µs (640.2 Hz), no ramp; then 1 s dwell | 150 → 0 |
| 46 | `velocity` | 0.25 fs/s + | 1/16 | 360 | 128 | + | 128 edges, uniform 250000 µs (4.0 Hz), no ramp; then 1.5 s dwell | 0 → 8 |
| 47 | `velocity` | 0.25 fs/s - | 1/16 | 360 | 128 | − | 128 edges, uniform 250000 µs (4.0 Hz), no ramp; then 1.5 s dwell | 8 → 0 |
| 48 | `velocity` | 0.5 fs/s + | 1/16 | 360 | 128 | + | 128 edges, uniform 125000 µs (8.0 Hz), no ramp; then 1.5 s dwell | 0 → 8 |
| 49 | `velocity` | 0.5 fs/s - | 1/16 | 360 | 128 | − | 128 edges, uniform 125000 µs (8.0 Hz), no ramp; then 1.5 s dwell | 8 → 0 |
| 50 | `velocity` | 1 fs/s + | 1/16 | 360 | 192 | + | 192 edges, uniform 62500 µs (16.0 Hz), no ramp; then 1.5 s dwell | 0 → 12 |
| 51 | `velocity` | 1 fs/s - | 1/16 | 360 | 192 | − | 192 edges, uniform 62500 µs (16.0 Hz), no ramp; then 1.5 s dwell | 12 → 0 |
| 52 | `velocity` | 2 fs/s + | 1/16 | 360 | 384 | + | 384 edges, uniform 31250 µs (32.0 Hz), no ramp; then 1.5 s dwell | 0 → 24 |
| 53 | `velocity` | 2 fs/s - | 1/16 | 360 | 384 | − | 384 edges, uniform 31250 µs (32.0 Hz), no ramp; then 1.5 s dwell | 24 → 0 |
| 54 | `velocity` | 5 fs/s + | 1/16 | 360 | 960 | + | 960 edges, uniform 12500 µs (80.0 Hz), no ramp; then 1.5 s dwell | 0 → 60 |
| 55 | `velocity` | 5 fs/s - | 1/16 | 360 | 960 | − | 960 edges, uniform 12500 µs (80.0 Hz), no ramp; then 1.5 s dwell | 60 → 0 |
| 56 | `velocity` | 10 fs/s + | 1/16 | 360 | 1,920 | + | 1920 edges, uniform 6250 µs (160.0 Hz), no ramp; then 1.5 s dwell | 0 → 120 |
| 57 | `velocity` | 10 fs/s - | 1/16 | 360 | 1,920 | − | 1920 edges, uniform 6250 µs (160.0 Hz), no ramp; then 1.5 s dwell | 120 → 0 |
| 58 | `velocity` | 20 fs/s + | 1/16 | 360 | 3,840 | + | 3840 edges, uniform 3125 µs (320.0 Hz), no ramp; then 1.5 s dwell | 0 → 240 |
| 59 | `velocity` | 20 fs/s - | 1/16 | 360 | 3,840 | − | 3840 edges, uniform 3125 µs (320.0 Hz), no ramp; then 1.5 s dwell | 240 → 0 |
| 60 | `velocity` | 50 fs/s + | 1/16 | 360 | 9,600 | + | 9600 edges, uniform 1250 µs (800.0 Hz), no ramp; then 1.5 s dwell | 0 → 600 |
| 61 | `velocity` | 50 fs/s - | 1/16 | 360 | 9,600 | − | 9600 edges, uniform 1250 µs (800.0 Hz), no ramp; then 1.5 s dwell | 600 → 0 |
| 62 | `velocity` | 100 fs/s + | 1/16 | 360 | 9,600 | + | 9600 edges, uniform 625 µs (1600.0 Hz), no ramp; then 1.5 s dwell | 0 → 600 |
| 63 | `velocity` | 100 fs/s - | 1/16 | 360 | 9,600 | − | 9600 edges, uniform 625 µs (1600.0 Hz), no ramp; then 1.5 s dwell | 600 → 0 |
| 64 | `velocity` | 200 fs/s + | 1/16 | 360 | 9,600 | + | 9600 edges, uniform 312 µs (3205.1 Hz), no ramp; then 1.5 s dwell | 0 → 600 |
| 65 | `velocity` | 200 fs/s - | 1/16 | 360 | 9,600 | − | 9600 edges, uniform 312 µs (3205.1 Hz), no ramp; then 1.5 s dwell | 600 → 0 |
| 66 | `velocity` | 400 fs/s + | 1/16 | 360 | 9,600 | + | 9600 edges, uniform 156 µs (6410.3 Hz), no ramp; then 1.5 s dwell | 0 → 600 |
| 67 | `velocity` | 400 fs/s - | 1/16 | 360 | 9,600 | − | 9600 edges, uniform 156 µs (6410.3 Hz), no ramp; then 1.5 s dwell | 600 → 0 |
| 68 | `reversal` | A=1 us f=0.2 Hz | 1/16 | 360 | 32 | ± | 8 × (+1, −2, +1 edges) at 0.8 Hz; then 1 s dwell | 0 |
| 69 | `reversal` | A=2 us f=0.2 Hz | 1/16 | 360 | 64 | ± | 8 × (+2, −4, +2 edges) at 1.6 Hz; then 1 s dwell | 0 |
| 70 | `reversal` | A=4 us f=0.2 Hz | 1/16 | 360 | 128 | ± | 8 × (+4, −8, +4 edges) at 3.2 Hz; then 1 s dwell | 0 |
| 71 | `reversal` | A=8 us f=0.2 Hz | 1/16 | 360 | 256 | ± | 8 × (+8, −16, +8 edges) at 6.4 Hz; then 1 s dwell | 0 |
| 72 | `reversal` | A=16 us f=0.2 Hz | 1/16 | 360 | 512 | ± | 8 × (+16, −32, +16 edges) at 12.8 Hz; then 1 s dwell | 0 |
| 73 | `reversal` | A=32 us f=0.2 Hz | 1/16 | 360 | 1,024 | ± | 8 × (+32, −64, +32 edges) at 25.6 Hz; then 1 s dwell | 0 |
| 74 | `reversal` | A=1 us f=1 Hz | 1/16 | 360 | 32 | ± | 8 × (+1, −2, +1 edges) at 4.0 Hz; then 1 s dwell | 0 |
| 75 | `reversal` | A=2 us f=1 Hz | 1/16 | 360 | 64 | ± | 8 × (+2, −4, +2 edges) at 8.0 Hz; then 1 s dwell | 0 |
| 76 | `reversal` | A=4 us f=1 Hz | 1/16 | 360 | 128 | ± | 8 × (+4, −8, +4 edges) at 16.0 Hz; then 1 s dwell | 0 |
| 77 | `reversal` | A=8 us f=1 Hz | 1/16 | 360 | 256 | ± | 8 × (+8, −16, +8 edges) at 32.0 Hz; then 1 s dwell | 0 |
| 78 | `reversal` | A=16 us f=1 Hz | 1/16 | 360 | 512 | ± | 8 × (+16, −32, +16 edges) at 64.0 Hz; then 1 s dwell | 0 |
| 79 | `reversal` | A=32 us f=1 Hz | 1/16 | 360 | 1,024 | ± | 8 × (+32, −64, +32 edges) at 128.0 Hz; then 1 s dwell | 0 |
| 80 | `close` | repeat two cycles + | 1/16 | 360 | 128 | + | 128 × (1 edge, then 150 ms dwell) | 0 → 8 |
| 81 | `close` | repeat two cycles - | 1/16 | 360 | 128 | − | 128 × (1 edge, then 150 ms dwell) | 8 → 0 |
| 82 | `close` | final reference standstill | 1/16 | 360 | 0 | — | hold 30 s, no STEP edges | 0 |

---

# Run B — microstep resolution

13 actions, 4.75 min planned, 89,660 STEP edges, 650,110 telemetry rows. Ran
57 minutes after Run A on the same bench and wiring, generated by
[`scripts/microstep_resolution_design.py`](scripts/microstep_resolution_design.py).
It asks one question: how finely does the motor actually step at the fine end
of the MRES table?

## A leftover MRES, and why action 0 exists

Run B's startup readback was `chopconf=0x07028055` — MRES code 7, i.e. **1/2
microstep**, not the 1/16 that Run A started from. Because the runner does not
reset the board, `CHOPCONF` still held whatever the previously executed plan
left behind. This is harmless (the firmware sets MRES per action regardless),
but it means the `# DRIVER` line in a recording describes the state *inherited*
at `PLAN RUN`, not the state any action ran under. Read MRES from the `mres`
column of the `D` rows instead.

Action 0 is a 1 s hold at **`MRES = 1`** whose only purpose is to snap the
driver's internal microstep counter to a full-step boundary. The driver's
`MSCNT` survives ESP32 resets, and an earlier test that afternoon showed the
rotor jumping two full steps on an MRES change from an unaligned position. It
is an alignment action, not a measurement.

## Structure

For each of 1/32, 1/64, and 1/128, in that order:

1. **settle** — 1 s hold immediately after the MRES change.
2. **fwd** — one complete shaft revolution (200 × MRES microsteps) forward, one
   microstep at a time, with a short dwell: 7 ms at 1/32, 3 ms at 1/64, 2 ms at
   1/128. These dwells are *not* settling times — this is a slow quasi-static
   rotation, not a settled staircase, and must not be read as one.
3. **back** — the identical revolution in reverse.
4. **pm1** — ten ±1-microstep out-and-back responses, 200 ms either side.

Every edge in Run B is marked: the staircases and step responses are all
edge-marking kinds.

## Blocks

| block | actions | STEP edges | duration | MRES used | I req (mA) |
|---|--:|--:|--:|---|---|
| `align` | 1 | 0 | 0.02 min | 1/1 | 360 |
| `m032_settle` | 1 | 0 | 0.02 min | 1/32 | 360 |
| `m032_fwd` | 1 | 6,400 | 0.75 min | 1/32 | 360 |
| `m032_back` | 1 | 6,400 | 0.75 min | 1/32 | 360 |
| `m032_pm1` | 1 | 20 | 0.07 min | 1/32 | 360 |
| `m064_settle` | 1 | 0 | 0.02 min | 1/64 | 360 |
| `m064_fwd` | 1 | 12,800 | 0.64 min | 1/64 | 360 |
| `m064_back` | 1 | 12,800 | 0.64 min | 1/64 | 360 |
| `m064_pm1` | 1 | 20 | 0.07 min | 1/64 | 360 |
| `m128_settle` | 1 | 0 | 0.02 min | 1/128 | 360 |
| `m128_fwd` | 1 | 25,600 | 0.85 min | 1/128 | 360 |
| `m128_back` | 1 | 25,600 | 0.85 min | 1/128 | 360 |
| `m128_pm1` | 1 | 20 | 0.07 min | 1/128 | 360 |
| **total** | **13** | **89,660** | **4.75 min** | | |

## Every action, in order

| # | block | label | MRES | I req (mA) | STEP edges | DIR | commanded STEP/DIR pattern | cmd span (full steps) |
|--:|---|---|:-:|--:|--:|:-:|---|---|
| 0 | `align` | full-step alignment | 1/1 | 360 | 0 | — | hold 1 s, no STEP edges | 0 |
| 1 | `m032_settle` | 1/32 settle | 1/32 | 360 | 0 | — | hold 1 s, no STEP edges | 0 |
| 2 | `m032_fwd` | 1/32 rev +6400 | 1/32 | 360 | 6,400 | + | 6400 × (1 edge, then 7 ms dwell) | 0 → 200 |
| 3 | `m032_back` | 1/32 rev -6400 | 1/32 | 360 | 6,400 | − | 6400 × (1 edge, then 7 ms dwell) | 200 → 0 |
| 4 | `m032_pm1` | 1/32 +-1 x10 | 1/32 | 360 | 20 | ± | 10 × (+1 edge, 200 ms, −1 edge, 200 ms) | 0 |
| 5 | `m064_settle` | 1/64 settle | 1/64 | 360 | 0 | — | hold 1 s, no STEP edges | 0 |
| 6 | `m064_fwd` | 1/64 rev +12800 | 1/64 | 360 | 12,800 | + | 12800 × (1 edge, then 3 ms dwell) | 0 → 200 |
| 7 | `m064_back` | 1/64 rev -12800 | 1/64 | 360 | 12,800 | − | 12800 × (1 edge, then 3 ms dwell) | 200 → 0 |
| 8 | `m064_pm1` | 1/64 +-1 x10 | 1/64 | 360 | 20 | ± | 10 × (+1 edge, 200 ms, −1 edge, 200 ms) | 0 |
| 9 | `m128_settle` | 1/128 settle | 1/128 | 360 | 0 | — | hold 1 s, no STEP edges | 0 |
| 10 | `m128_fwd` | 1/128 rev +25600 | 1/128 | 360 | 25,600 | + | 25600 × (1 edge, then 2 ms dwell) | 0 → 200 |
| 11 | `m128_back` | 1/128 rev -25600 | 1/128 | 360 | 25,600 | − | 25600 × (1 edge, then 2 ms dwell) | 200 → 0 |
| 12 | `m128_pm1` | 1/128 +-1 x10 | 1/128 | 360 | 20 | ± | 10 × (+1 edge, 200 ms, −1 edge, 200 ms) | 0 |

---

# Reading the command record against the response

Three things about the *command* sequence shape how far the encoder data can be
trusted. All three are visible in the `# ACTION_END` lines without any
reanalysis; comparing `encoder − command_u256 × 100/256` at each boundary gives
the residual in counts.

### Both runs were commanded back to exactly zero

Run A's plan is validated to end at `command_u256 = 0` and did; Run B likewise.
So the closing residual is a direct measure of accumulated position error.

### Run A's residual is clean through the static map and wanders afterwards

| Through action | Residual at `command_u256 = 0` |
|---|---|
| 0–16 (`warmup`, `alignment`, `static_map`) | within ±6 counts (±0.06 full step) |
| 17–25 (`current_sweep`) | within ±3 counts at every zero crossing |
| 26–45 (`dynamics`) | reaches +111 counts |
| 46–67 (`velocity`) | swings between −60 and +170 counts |
| 68–79 (`reversal`) | ends at +267 counts |
| 80–82 (`close`) | +267 counts, stable |

The static map — the campaign's primary measurement — sits in the clean region.
That is the single most useful thing this table says.

### Where the +267 counts came from

The data README flags the closing +267 counts (2.67 full steps, 4.8°) as
unexplained and asks that it be located before the later blocks are trusted.
It resolves into two distinct contributions:

**1. MRES switching in the `dynamics` block.** The residual is +1 count at the
end of action 25 and +111 by action 45. The jumps line up with the
`MRES = 16 ↔ 1` transitions, and they are consistent in sign — entering
`MRES = 1` shifts the residual negative, returning to `MRES = 16` shifts it
positive, by tens of counts each time. This is the same MRES-change behaviour
that Run B's action 0 was added to work around. Consequence: the full-step and
1/16 step-response records within `dynamics` sit at *different* rotor offsets,
and absolute position bookkeeping across that block is not reliable. The step
responses themselves are local measurements and are unaffected.

**2. Monotonic creep during the final reversal action.** Action 79 (`A=32 us
f=1 Hz`, 1,024 edges at 128 Hz) alone accounts for +234 counts. Sampling the
encoder at each cycle zero-crossing through that action gives −8, +40, +13,
+42, +69, +83, +73, +83, +85, +136, +157, +189, +178, +227, +243, +268 — a
steady forward walk of roughly +30 counts per cycle. There is **no discrete
jump anywhere in the action**: the largest sample-to-sample change is under 8
counts. The two neighbouring actions barely move at all (action 78, the same
1 Hz at half the amplitude, contributes +5; action 73, the same amplitude at
0.2 Hz, +13), so it is specific to the largest amplitude at the higher
frequency.

Two candidate explanations, and the data here does not separate them:

- genuine directional creep — the rotor ratcheting forward under symmetric
  commanded excursions of ±2 full steps at 1 Hz;
- encoder over-counting on the reversals, the same mechanism already suspected
  for Run B's 1/32 revolution, where quadrature ringing at the largest and
  fastest steps is the leading candidate.

A gradual walk cannot be lost steps: a stepper locked to its field loses
position in jumps of whole electrical cycles (400 counts), not 30 counts per
cycle. Repeating action 79 on its own, ideally with a second independent angle
reference, would settle it.

### Implications for the data as it stands

- `static_map`, `alignment`, and `current_sweep` position data are supported by
  a residual that never leaves ±6 counts.
- `dynamics` step-response *waveforms* are usable; `dynamics` absolute
  positions are not.
- `velocity` lag figures should be taken from the constant-rate interior of
  each plateau, as designed, not from boundary positions — the largest
  single-boundary jump in the whole run is +415 counts at action 59.
- `reversal` results at A = 32 and 1 Hz are suspect; the smaller amplitudes
  move the residual by only a few counts each and look sound.

### Suggested firmware additions for the next run

Each of these closes a gap that this reconstruction had to work around:

1. Log `CS_ACTUAL` per action, not just the requested milliamps, so realised
   current is recorded rather than computed.
2. Log `MSCNT` at every `ACTION_BEGIN` and `ACTION_END`, which would turn the
   MRES-transition offsets from an inference into a measurement.
3. Re-print the `# DRIVER` readback after each MRES or current change, so the
   recording carries the state actually in force.
4. Mark edges in `move` and `velocity` too, or at least log the commanded edge
   count per action, so all 339,840 edges are accounted for in the file.

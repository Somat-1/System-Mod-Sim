# Stage identification: what is commanded to the motor

This document records the driver settings and the motion commands in each
sequence of the stage identification campaign. It is taken from
`scripts/sequence_design.py`, the generated `sequence/sequence_plan.csv`,
`firmware/stage_identification/stage_identification.ino` and
`scripts/run_campaign.py`.

## Units

| Quantity | Value |
|---|---|
| Position | full steps from the zero set at the start of the run. 100 full steps = 1 mm (1 full step = 10 µm). |
| Speed | full steps per second |
| Microstep resolution | MRES 16 for every move: 1 microstep = 1/16 full step = 0.625 µm, 16 STEP pulses per full step |
| Electrical period | 4 full steps |
| Screw revolution | 200 full steps |

## Driver settings (TMC2209)

Applied at boot over UART and confirmed by reading the registers back before
any motion is allowed:

| Setting | Value |
|---|---|
| Chopper | SpreadCycle (`en_spreadCycle = 1`, `TPWMTHRS = 0`), StealthChop off |
| Interpolation | off (`intpol = 0`): each pulse is one native 1/16 microstep |
| Microstep resolution | MRES 16 |
| Run current | 360 mA requested (IRUN 9/31, `vsense = 1`, 0.10 Ω sense resistor, about 330 mA RMS) |
| Hold current | equal to run current (IHOLD = IRUN), `IHOLDDELAY = 0`, `TPOWERDOWN = 0` |
| Current reference | internal only (`I_scale_analog = 0`), not the module's potentiometer |
| Chopper timing | `TOFF = 5`, blank time 24 clocks |
| Enable | EN (GPIO9) held high until the settings are confirmed, then low for the rest of the session |
| Pins | STEP GPIO7, DIR GPIO8, EN GPIO9, UART TX 17 / RX 18 |

## How a move is executed

Each plan row is one move: target, speed, dwell, MRES, microstep count and
microstep dwell.

**Ordinary move** (`micro_steps = 0`):

- one continuous pulse train from the current position to the target;
- symmetric acceleration-limited profile, acceleration 4000 full steps/s²;
- the rate for each pulse is `min(cruise speed, √(2 · 4000 · s))` full
  steps/s, where `s` is the distance to the nearer end of the move;
- pulse rate = rate × 16, at most 20 kHz;
- STEP high time 3 µs; DIR set 20 µs before the first pulse;
- then the motor is held still for `dwell` seconds.

**Staircase move** (`micro_steps > 0`): `micro_steps` single pulses in one
direction. Each pulse is followed by a 1 ms gap and then `micro_dwell` seconds
held still. No continuous motion.

## Run sequence

Before the plan runs:

1. Soft limit set to ±1100 full steps (±11 mm) about zero.
2. All 230 moves uploaded and validated by the board: target inside the limit,
   pulse rate at most 20 kHz, MRES a power of two.
3. `ZERO`: the present position becomes 0.
4. Preflight `CHECK`: +10 full steps and back at 25 full steps/s (160 pulses
   each way). Skipped with `--blind`.
5. `ZERO` again, then `PLAN RUN`. The plan must start at 0 and end at 0.

The board then executes the 230 moves in order, about 42.7 min in total.

| # | Sequence | Starts at | Moves | Duration |
|---|---|---|---|---|
| 1 | Gain calibration | 0 s | 21 | 3.6 min |
| 2 | Microstep staircase | 214 s | 10 | 10.0 min |
| 3 | Settling matrix | 833 s | 91 | 19.3 min |
| 4 | Bidirectional repeatability | 2029 s | 101 | 8.4 min |
| D | Drift checks | interleaved | 7 | 1.4 min |

### 1 · Gain calibration

- **Commanded:** targets ±L/2 for L = 200, 400, 600, 800, 1000 full steps,
  taken in the order +L/2 then −L/2 for each L. The whole set is done twice,
  followed by a return to 0.
- **Moves as commanded**, in full steps:

  ```
  +100  −200  +300  −400  +500  −600  +700  −800  +900  −1000
  +600  −200  +300  −400  +500  −600  +700  −800  +900  −1000   +500
  ```

- **Speed:** 70 full steps/s (1120 pulses/s) for every move.
- **Dwell:** 2 s after every move.
- **Range of targets:** −500 to +500 full steps.

### 2 · Microstep staircase

Done at three locations, in this order: −500, 0 and +500 full steps. At each:

| Move | Commanded | Speed | Hold |
|---|---|---|---|
| 1 | to the location (ordinary move) | 70 full steps/s | 2 s |
| 2 | 128 single microsteps forward (+8 full steps) | single pulses | 0.75 s after each (96 s in total) |
| 3 | 128 single microsteps back to the location | single pulses | 0.75 s after each (96 s in total) |

The phase ends with a return to 0 at 70 full steps/s, then a 2 s hold.

### 3 · Settling matrix

- **Commanded:** 90 relative moves covering every combination of:

  | Factor | Values |
  |---|---|
  | distance | 10, 100, 1000 full steps (0.1, 1, 10 mm) |
  | speed | 25, 70, 200 full steps/s (400, 1120, 3200 pulses/s) |
  | direction | + and − |
  | repeats | 5 |

- **Order:** shuffled with a fixed random seed (7).
- **Each move** goes from the previous target. If it would pass ±1000 full
  steps, it goes the same distance in the opposite direction instead.
- **Dwell:** 5 s after every move.
- **End:** a return to 0 at 70 full steps/s, then 2 s.
- **After a drift check** (below) the next move starts from 0, so the plan
  contains four settling moves of 400, 750, 770 and 980 full steps.

### 4 · Bidirectional repeatability

- **Targets:** −1000, −500, 0, +500 and +1000 full steps.
- **Approaches:** each target is reached from below and from above.
- **Repeats:** the full set of targets and both directions is done 5 times.
- **One approach is two moves:**

  | Move | Commanded | Speed | Hold |
  |---|---|---|---|
  | run-up | 50 full steps past the target, on the approach side | 70 full steps/s | 0.5 s |
  | target | the target itself | 70 full steps/s | 2.5 s |

- **Furthest point:** the run-ups reach ±1050 full steps.
- **End:** a return to 0 at 70 full steps/s, then 2 s.

### D · Drift checks

- **Commanded:** a move to 0 at 70 full steps/s, then a 3 s hold.
- **When:** inserted about every 300 s. Never inside a staircase, and never
  between a run-up and its target.
- **In this plan:** 7 checks, at roughly 415, 821, 1135, 1456, 1772, 2091 and
  2412 s.

# E5: compensated validation (ISO 230-2, paired) — procedure v2

The governing specification is [`E5_validation_procedure.md`](E5_validation_procedure.md).
This implementation applies the frozen geometric and reversal model in the
command path and compares uncompensated (U) and compensated (C) runs within one
session. Status: **model, sequence and firmware generated; no E5 data yet**.

## Implemented design

- Frozen model: E2 data only, D2 drift variant, H = 2, K_s = 2, K_C = 8,
  lambda = 0.32. After fixing H = 2, this is the minimum-CV E2 M_C candidate.
  E3 determines only the harmonic order; no E3 measurements enter the fit.
- Reversal: frozen E2 beta0 = -1.757 um.
- ISO block: ten 24-target bidirectional runs in `U C C U U C C U U C`
  order (480 stops).
- Reference blocks: R1-R4, each measuring 15 and 105 mm from both directions
  in the fixed order `15 up, 15 down, 105 up, 105 down` (16 stops).
- B3: 32 seeded targets from 5-110 mm, nominal offsets -32,-30,...,+30
  microsteps, two bidirectional repetitions (128 stops).
- Total: 624 recorded stops in one predicted 56.9 min session.

The command correction is

    x_c = x* - [G_hat(x*) - a0] - (d/2) beta0
    u_cmd = round(x_c * 1600)

The ESP32 evaluates no model; it executes the generated integer table.

## Files

| Path | Purpose |
|---|---|
| `E5_validation_procedure.md` | governing v2 experimental and analysis procedure |
| `e5_config.json` | hardware, timing and session definition |
| `e5_freeze_model.py` | freezes the E2-only H=2 model |
| `model/geom_model_E5.json` | frozen model and covariance, hashed before the run |
| `e5_make_table.py` | builds, validates and hashes the complete session |
| `sequence/e5_command_table.csv` | one row per recorded stop |
| `sequence/e5_measurement_plan.csv` | all state, move and stop actions in execution order |
| `sequence/e5_plan_summary.json` | hashes, checks, warnings and duration estimates |
| `sequence/e5_b3_targets.csv` | B3 nominal targets and electrical offsets |
| `firmware/e5_compensated_validation/e5_config.h` | generated flash-resident action table |
| `firmware/e5_compensated_validation/e5_compensated_validation.ino` | ESP32 RMT executor |
| `render_e5_sequence.py` | sequence, correction and expected-result renders |

## Rebuild

```sh
python3 e5_freeze_model.py
python3 e5_make_table.py
python3 render_e5_sequence.py
```

`e5_make_table.py --check-only` validates without rewriting generated files.
Always regenerate and record the model and table SHA-256 values before flashing.

## Firmware operation

The firmware uses the same ESP32-S3 pins, TMC2209 configuration, buffered RMT
trapezoidal motion and directional switch handling as E1-E4.

1. Perform the established E0 warm-up using `step_timing_test` (`h`, then `w`).
2. Flash `e5_compensated_validation.ino` with its generated `e5_config.h`.
3. Start the IDS recording and lossless serial logger.
4. Press `g` to re-home, centre and run. (`p`, then `r` is available only when
   the stage is already homed and physically at 61.35 mm.)
5. Do not close the logger after motion ends: the firmware prints `SEQ_DONE`
   first and then dumps all 624 RAM-buffered timing records.

Each dumped `# STOP` record includes the table row, block, run, mode, direction,
target, nominal and commanded microsteps, electrical offset, MSCNT, move-end and
record-window timestamps, and TMC status. The session header prints the full
model hash, table hash, B3 seed and home MSCNT.

## Pre-run result and explicit procedure warnings

All executable safety and consistency checks pass. Command correction spans
-3.75 to +4.38 um and |q| <= 0.311 um. The generated model hash begins
`09200766589f`.

Three requirements in the supplied procedure cannot simultaneously be met by
the existing E0 v4 assets and are therefore reported, not hidden:

1. **Timing margin:** the full plan is 56.9 min. It is below the firmware's
   58 min gate but does not leave an eight-minute margin. Reducing B3 to one
   repetition gives 52.4 min; change `b3_repetitions` to 1 and regenerate if
   that fallback is required before the run.
2. **U electrical grid:** the procedure asks both for unchanged E0 v4 commands
   and for all U commands to be multiples of 64 microsteps. The E0 v4 target
   list contains off-grid commands, so the implementation preserves the more
   fundamental paired-baseline requirement.
3. **B3 final phase:** nominal B3 targets cover -32,-30,...,+30 exactly.
   Direction-dependent compensation shifts `u_cmd`, so the final commanded
   offsets are not the same exact even sequence. The actual offset and MSCNT
   are logged at every stop for the registered B3 analysis.

Do not change these choices after viewing E5 data. If a warning is resolved,
regenerate the table, header, hashes and renders together before the run.

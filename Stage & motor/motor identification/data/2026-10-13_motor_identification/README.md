# Motor identification run — 2026-10-13

The full rotor-side identification campaign (`../../sequence/motor_identification_plan.csv`,
copied here as `plan.csv`): 83 actions, 37.1 min, started 15:50 local
(13:50 UTC), completed (`RUN_COMPLETE`).

| File | |
|---|---|
| `motor_identification_run.csv` | the recording as written by `scripts/run_campaign.py`: 1,903,061 `D` sample rows plus `#` event lines (86 MB, Git LFS) |
| `plan.csv` | the plan that was run |

## Conditions

| | |
|---|---|
| Bench | Stepper Isolated wiring: STEP 7, DIR 8, EN 9, UART 17/18, encoder A/B/Z 4/5/6 |
| Driver | TMC2209, R_sense 0.10 Ω, SpreadCycle, `intpol = 0`, `I_scale_analog = 0`, IHOLD = IRUN = 9 (360 mA requested, about 330 mA RMS), `TPOWERDOWN = 0` |
| Encoder | on the motor shaft, 20,000 counts/rev (5000 lines, x4), index once per rev; positive command = positive counts |
| Firmware | `firmware/motor_identification` with the PCNT watch-point fix (counts accumulate past ±30,000) |
| Load | as fitted on the bench on 2026-10-13 |

## Notes

- Columns: `record_type,seq,t_us,action,event,command_u256,encoder,index_count,current_ma,mres`.
  `command_u256` is in 1/256 full step; 1 full step = 100 encoder counts.
- **The run ended with the encoder at +267 counts while the command was 0**
  (2.67 full steps, 4.8°). The motor gained or lost position somewhere during
  the run. Find where before trusting the later blocks: compare
  `encoder - command` at each `# ACTION_END`.
- An earlier attempt the same afternoon (`recordings/motor_identification_20261013_154541.csv`,
  not copied here) used firmware without the watch-point fix. Its encoder
  reading restarted at 0 every 30,000 counts. Do not use it.

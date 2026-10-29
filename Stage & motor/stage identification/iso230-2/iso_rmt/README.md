# ISO 230-2 baseline — RMT sequence

The v3 ISO 230-2 sequence on the motion and logging path proven by the E1, E2
and E4 runs of 2026-10-27. It replaces the CPU-timed step generator and the
50 Hz in-motion telemetry of `../firmware/iso230_2/` with:

- buffered 1 MHz RMT step generation (the 2/5/8 mm/s timing-test path);
- no serial output and no TMC UART access during motion; MSCNT and driver
  status are read only after each record window;
- homing on the board (`g`/`h`), so no separate homing or manual arming step.

## Defaults

| setting | value |
|---|---|
| motion | 2.0 mm/s, 40 mm/s², trapezoidal |
| driver | 1/16, SpreadCycle, 400 mA, hold = run, TPOWERDOWN 0 |
| stop | 2.0 s settle + 0.5 s record |
| targets | 24, seed 2302, 5.000–117.500 mm (v3's positions exactly) |
| runs | 5 bidirectional, 3 mm overrun at both ends |
| references | x_ref 61.35 mm from 58.35 mm: 3 datum, 1 per run, 3 end |
| soft limits | 2.000–120.695 mm, checked against the homed travel |
| duration | about 25.8 min after homing (5.05 min per run) |

The targets come from `place_targets()` in `../scripts/iso_design.py` for a
122.77 mm travel. The generator checks them against the `# TARGET` lines of the
v3 log `../data/terminal_20261021_183115.log`, so the result is directly
comparable with `../analysis/iso230_2_v3_results.json`.

## Sequence

`HOME → INIT → DATUM ×3 → 5 × (ISO_REF, ISO up, ISO down) → END_REF ×3 → PARK`

Each run is v3's `isoRun()`: an x_ref visit, a move to first target − 3 mm,
24 ascending stops, a move to last target + 3 mm, 24 descending stops, and the
return to first target − 3 mm. After each run but the last, a time watchdog
projects the finish; above 58 min the remaining runs are skipped, the end
references are still taken, and the session is marked incomplete. 240 analysis
stops and 11 reference stops in total.

## Generate and render

```bash
python3 iso_rmt_sequence.py
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib python3 render_iso_rmt_sequence.py
```

The generator writes `sequence/iso_rmt_targets.csv`,
`sequence/iso_rmt_measurement_plan.csv`, `sequence/iso_rmt_plan_summary.json`
and `firmware/iso230_2_rmt/iso_rmt_config.h`. The config's SHA-256 is compiled
into the header and printed at boot, so edit the JSON, regenerate, and reflash.
The render is `rendered/iso_rmt_sequence.png`.

## Firmware

```bash
arduino-cli compile --upload -p /dev/cu.usbmodem3101 \
  --fqbn "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc" \
  firmware/iso230_2_rmt
```

Then run `python3 iso_rmt_terminal.py [port]`; logs go to `data/`. Commands:

- `g` — home both ends, check the travel holds the soft limits, centre at
  x_ref from 58.35 mm, then run the sequence;
- `h` — home and centre only (arms the coordinate);
- `p` — arm at x_ref manually, with the stage already homed and parked there;
- `r` — run the sequence from an armed coordinate;
- `s` — status; `t` — plan summary; `0` — release the motor and invalidate
  the coordinate; Ctrl+C — request abort.

The warm-up is not part of this sketch: run it beforehand with
`../firmware/step_timing_test` (`h`, then `w`), then flash this one and press
`g`.

Each `# STOP` line carries state, window, run, direction, target index,
target in µsteps and mm, the ESP move-start, arrival, record-start and
record-end times, and the post-window MSCNT. `sync_gpio` is −1, so the IDS
record is cut using the event clock map.

# Open-loop trajectory diagnostic

## Purpose

Run this test before another ISO 230-2 dataset. It reproduces the failed
reference transfer (`2 mm -> x_ref - 3 mm -> x_ref -> 2 mm`) under four motion
profiles, holds each endpoint for 10 seconds, and checks the home switch after
each profile. The existing ISO sequence is not changed.

The stage has no position feedback. The firmware can prove that every requested
STEP pulse was emitted, but it cannot prove that the rotor, screw, coupling, and
carriage followed those pulses. The limit switches can detect accumulated
position loss at the ends; the independently recorded interferometer trace is
still required to validate intermediate positions.

Longer dwell has a deliberately limited role: it reveals whether the carriage
continues to settle after the pulse train ends. It cannot make a stalled
stepper execute missing steps after the pulse train has stopped.

## Configured sequence

The parameters live in `config/iso_config.json` and are compiled into the same
traceable config hash as the ISO settings.

| Phase | Velocity | Acceleration | Repetitions | Endpoint hold | What it isolates |
|---|---:|---:|---:|---:|---|
| Switch baseline | 0.15 mm/s | current ramp | 5 closures | — | Home-switch repeatability/noise floor |
| `current` | 5.0 mm/s | 40 mm/s² | 2 | 10 s | Reproduces the current reference move |
| `low_accel` | 5.0 mm/s | 10 mm/s² | 2 | 10 s | Acceleration sensitivity at the same maximum speed |
| `low_speed` | 2.0 mm/s | 40 mm/s² | 2 | 10 s | Maximum-speed/resonance sensitivity |
| `conservative` | 1.0 mm/s | 5 mm/s² | 2 | 10 s | Strongly derated candidate fix |

For each profile the controller performs:

1. Move from the lower soft limit (normally 2 mm) to `x_ref - 3 mm`.
2. Execute the final 3 mm fixed-side approach to `x_ref`.
3. Hold for 10 seconds.
4. Return directly to the lower soft limit.
5. Hold for 10 seconds.
6. Repeat steps 1–5 once.
7. Approach the home switch at 0.15 mm/s, log the closure correction, re-zero
   there, and back off by 2 mm before the next profile.

At the end, the stage parks at `x_ref` using the conservative profile. Allow
approximately 18–20 minutes; the exact duration is logged rather than assumed.

## How to run it

1. Inspect the mechanics and wiring using the pre-run checklist below.
2. Regenerate the header after any config edit:

   ```text
   python scripts/iso_design.py --rom 122.695
   ```

3. Compile and flash `firmware/iso230_2/iso230_2.ino`.
4. Start `scripts/iso_terminal.py` so the complete board stream is logged.
5. Start the interferometer recording before motion. Use the normal 1 kHz
   acquisition and keep it running until `# DIAG_DONE`.
6. Let boot auto-home and park at the centre. If auto-home was cancelled, run
   menu option `6` first.
7. Press `9`. Do not touch the table or stage. `Ctrl+C` is the immediate abort.
8. Save together:
   - the raw interferometer CSV;
   - the matching `terminal_YYYYMMDD_HHMMSS.log`;
   - the config hash printed in `# DIAG_BEGIN`;
   - notes on motor/driver temperature, supply voltage, unusual sound, and the
     point at which any stall or vibration was heard.

Do not start menu option `2` after this test. The diagnostic deliberately
re-homes between profiles and is not an ISO measurement session.

## Log evidence

- `# DIAG_MOVE_START/END`: profile, leg, requested/emitted pulse counts, actual
  elapsed time, and completion flag.
- `# DIAG_HOLD_START/END`: exact 10-second stationary windows for settling
  analysis.
- `# DIAG_HOME_CHECK`: open-loop coordinate at which the physical home switch
  closed. A nonzero value is the correction accumulated since the preceding
  re-zero.
- `# DRIVER`, `# WARN`, `# ERR`: driver readback and electrical/thermal faults.
- `# DIAG_DONE`: complete test and final commanded position.

## Interpretation we will apply to the returned data

| Observation | Interpretation |
|---|---|
| Interferometer keeps moving during a 10 s hold and reaches the target late | Settling/compliance is relevant; select a dwell from the measured settling time |
| Interferometer becomes flat quickly but is short of target | More dwell cannot fix it; investigate lost steps, slip, binding, or the distance scale |
| `current` fails and `low_accel` passes | Acceleration/available torque is the main limit |
| Both 5 mm/s profiles fail and `low_speed` passes | Maximum speed, resonance, or pulse-rate/torque margin is the main limit |
| Only `conservative` passes | Use a derated profile and then confirm it with a longer repeat test |
| All profiles show the same displacement ratio error | Check lead, steps/rev, MRES readback, and mechanism ratio before tuning dwell |
| Home correction grows beyond the switch baseline | The open-loop position was lost through missed steps or mechanical slip |
| Home correction is small but an outward endpoint is wrong | Outward/return errors may cancel, or the scale/overlay is wrong; use the IFM trajectory |
| Every move has `emitted_pulses=requested_pulses` | The controller did not truncate the move; this does not prove mechanical following |
| A move has `complete=0`, an unexpected switch, or a driver fault | Abort the diagnosis; the dataset is not a valid profile comparison |

## Pre-run validation checklist

- [ ] Coupling, motor shaft clamp, screw clamp, and stage fasteners are tight;
      add witness marks so slip can be seen.
- [ ] With power safely removed, the usable travel has no hard spot, cable pull,
      collision, or abrupt preload change.
- [ ] TMC2209 and motor current ratings permit 400 mA continuously; cooling and
      supply voltage are stable.
- [ ] Menu `8` reports both switches open at mid travel and each switch changes
      state when pressed by hand.
- [ ] Driver readback confirms SpreadCycle, `intpol=0`, MRES 16, IRUN/IHOLD, and
      no overtemperature or short flags.
- [ ] Config hash is recorded and `encoder=0` is shown in `# READY`.
- [ ] Homing completes without an unexpected switch and the measured ROM is
      consistent with earlier sessions.
- [ ] Interferometer signal is valid, its counter is streaming, and acquisition
      begins before the diagnostic.
- [ ] Nobody touches the stage/table during the run; audible stalls are noted,
      not ignored.

## Ready-to-confirm checklist after analysis

These are gates for a follow-up confirmation test, not automatic claims made by
the open-loop controller.

- [ ] All diagnostic moves report `complete=1` and emitted pulses equal requested
      pulses.
- [ ] No unexpected limit event, TMC fault, reset, serial gap, or interferometer
      counter discontinuity occurred.
- [ ] The home-switch baseline is repeatable; post-profile closure corrections
      remain within `max(0.020 mm, 3 × baseline peak-to-peak)`.
- [ ] Both repetitions of the selected profile reach the reference displacement
      without a missed plateau, jump, or steady accumulated offset.
- [ ] Measured `low -> x_ref` displacement agrees with the commanded displacement
      to within `max(0.010 mm, 5 microsteps)` before formal ISO analysis.
- [ ] The last second of each 10-second hold is visibly stationary and has no
      meaningful residual slope relative to the expected repeatability.
- [ ] Outward and return trajectories have the same scale; no profile-dependent
      shortfall remains.
- [ ] The selected operating profile has margin. Do not use the fastest profile
      that only barely passes; derate velocity and acceleration from the observed
      failure boundary.
- [ ] A subsequent confirmation sequence completes at least five consecutive
      reference cycles with the selected profile and dwell, followed by a home
      closure check within the baseline envelope.

Only after every gate passes should the ISO firmware settings be updated and a
new ISO session started from a fresh home and warm-up. Do not home between ISO
runs: it changes the datum and approach history and can conceal lost position.
Use end switches as pre/post-session checks and as immediate safety limits, not
as a correction inserted into the formal target sequence.

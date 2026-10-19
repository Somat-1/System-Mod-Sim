# Mechanical chirp sweep

Two frequency sweeps of the stage, 1 Hz to 1 kHz, in both directions.

```
chirp/
├── README.md
├── scripts/
│   ├── chirp_design.py           both profiles, as data
│   └── render_chirp.py           draws them
├── sequence/
│   ├── chirp_sweeps.png
│   ├── chirp_plan.csv
│   └── chirp_summary.json
└── firmware/stage_chirp_sweep/   RUN FULL | RUN HIGH | CHECK | ABORT
```

---

## ⚠ 150–250 Hz is not driven, and that is deliberate

Two incidents on **2026-09-15** produced **uncontrolled motor motion**, reported
live at around 200 Hz. A smooth amplitude notch over the band did not prevent
them — not at 90 % suppression, and not widened to 65–285 Hz at 97 %.
**Attenuating that band was not sufficient, however aggressively it was done.**

So the FULL profile does not attenuate it. It does not visit it: the up sweep
runs 1 → 150 Hz, changes target to 250 Hz in a single step, and continues
250 → 1000 Hz. The down sweep mirrors that.

Do not narrow the gap to recover the missing data. If 150–250 Hz has to be
characterised, it needs a method that does not involve driving the motor through
it — a ring-down from an impulse, or a measurement with the motor de-energised.

Stay at the rig with the power within reach for any FULL run. The history is in
[`../../../ESP32S3_TMC2209_Chirp_Test/v3`](../../../ESP32S3_TMC2209_Chirp_Test/v3)
and the v3.2 sketch header.

---

## The two profiles

![chirp sweeps](sequence/chirp_sweeps.png)

| | Band | Duration | Why |
|---|---|---|---|
| **FULL** | 1–150 Hz + 250–1000 Hz | 4.47 min | cover everything it is safe to drive, both directions |
| **HIGH** | 250–1000 Hz, continuous | 3.93 min | resolve what is happening at 250 Hz |

Both keep v3's shared sweep rate of **7.8125 Hz/s**, **400 mA RMS**, and the
750 Hz sync marker bursts at each end for aligning an external recording.

### Amplitude: constant acceleration, capped (changed 2026-10-13)

```
A(f) = A0                      f <= 100 Hz        A0 = 1 full step
A(f) = A0 * (100 Hz / f)^2     above it           (256 microsteps at MRES 256)
```

| f [Hz] | 1–100 | 150 | 250 | 500 | 750 | 1000 |
|---|---|---|---|---|---|---|
| amplitude [microsteps @ 1/256] | 256 | 114 | 41 | 10.2 | 4.6 | 2.56 |

The previous ±3 full steps (48 microsteps at MRES 16) had two problems:

- **It wasn't constant.** The step-rate ceiling (200k microsteps/s) clamped it
  above 663 Hz, falling to 66% at 1 kHz.
- **It exceeded what a stepper tolerates.** A rotor that falls more than
  2 full steps behind the field slips a whole electrical cycle, and ±3 full
  steps lets that happen. It may be what produced the 2026-09-15 runaway.

The new law fixes both:

- **The command never exceeds ±1 full step**, so the rotor cannot slip even if
  it does not move at all.
- **Above 100 Hz the peak acceleration, and so the torque demanded, is
  constant.**
- **The peak step rate (161k/s, at the corner) stays 20% under the ceiling**,
  so the amplitude is never clamped.

MRES moved from 16 to 256 so the small high-frequency amplitudes are still
many microsteps.

The 150–250 Hz gap is kept. Whether the smaller amplitude makes that band safe
has not been tested, and should not be assumed.

```bash
python scripts/chirp_design.py     # -> chirp_plan.csv, chirp_summary.json
python scripts/render_chirp.py     # -> the figure
```

---

## About the thing at 250 Hz

Worth checking before treating it as a property of the stage: **250 Hz is
exactly where FULL jumps back in.** The sweep leaves at 150 Hz and resumes at
250 Hz, and that change of target is instantaneous. An instantaneous frequency
change is a step in the excitation, and a step is broadband — it deposits energy
well outside the frequency either side of it. A feature sitting precisely at the
re-entry frequency is therefore a suspect before it is a finding.

That is what HIGH is for, and why it is not simply "FULL with the low part
deleted":

- it never approaches the gap, so there is no jump anywhere in it;
- it **dwells at 250 Hz for 3 s before sweeping**, so the start frequency is
  reached by steady excitation rather than a discontinuity;
- it ends the same way, so the down sweep's arrival at 250 Hz is clean too.

Read the two together:

| HIGH shows | Then the 250 Hz feature is |
|---|---|
| the same feature | real — a property of the stage, and FULL's jump is incidental |
| nothing at 250 Hz | an artefact of the gap edge in FULL, not the mechanism |
| something different | worth a dwell sweep: hold a fixed frequency either side and compare |

If it turns out to be the gap edge, the fix is not to remove the gap — it is to
stop reading that part of the FULL record, and to use HIGH for anything in the
upper band.

---

## Running it

```
CHECK          configure and report, move nothing
RUN HIGH       the safer profile; start here
RUN FULL       7 s countdown, then the full sweep
ABORT          stop now, or cancel the countdown
```

With the host scripts (pyserial for the runner, Homebrew Python for the window):

```bash
python scripts/run_chirp.py --port /dev/cu.usbmodem3101 --profile high --out data/<date>_chirp
/opt/homebrew/bin/python3 scripts/monitor_chirp.py data/<date>_chirp     # live window
python scripts/analyse_chirp.py data/<date>_chirp
```

The board streams `# SEG,...` at each segment boundary and `# GAP_JUMP` when it
crosses the gap. Every 0.1 s window (at least one cycle) it also streams:

```
# W,t,f,cycles,f_issued,amp_target,amp_issued,amp_enc,phase_deg,lag_max,dc_err
```

This is demodulated on the board against the sweep's own phase. It shows that
the controller reached the frequency and amplitude it was asked for:

- `amp_issued` against `amp_target`, from the steps actually issued;
- `f_issued` against `f`, from zero crossings of the issued position;
- `lag_max`, the worst lag of issued steps behind the target.

`amp_enc` and `phase_deg` are the shaft encoder's response when one is fitted.
A sync check then aborts the run if the rotor is more than 1.75 full steps from
the command. With no encoder they read 0, the check is inactive, and the
response has to come from the external sensor.

`REINIT` retries the driver setup and prints the chip ID, for example after the
motor supply has been switched on. The LED
shows phase, because serial has proven unreliable to watch live on this rig:
yellow waiting, blue counting down, green running, cyan done, red failed.

**Run HIGH first.** It is entirely above the band that caused the incidents, it
answers the question you actually have, and it exercises the firmware before a
FULL run does.

---

## Two things to check before the first run

**`R_SENSE` differs between boards on this rig.** The sketch now uses
**0.10 Ω**, for the R100 module on the Stepper Isolated bench (2026-10-13); the
v6 campaign's board used 0.11 Ω. Getting it wrong mis-scales the coil current by ~10 %,
which at 400 mA is enough to change what the motor does near resonance. Read the
markings on the board that is actually fitted and set `R_SENSE_OHM` to match.

**No sweep has run on hardware yet.** The sketch was flashed on 2026-10-13 and
passed `CHECK` (driver ready, MRES 256, SpreadCycle, no interpolation), but has
not yet driven a sweep. Given what the 150–250 Hz band did to this rig once already,
treat the first run as a bring-up: `CHECK` first, then `RUN HIGH`, and watch it.

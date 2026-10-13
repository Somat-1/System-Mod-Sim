# Mechanical chirp sweep

Two frequency sweeps of the stage, 1 Hz to 1 kHz, in both directions.

```
chirp/
├── README.md
├── scripts/
│   ├── chirp_design.py           both profiles, as data
│   └── render_chirp.py           draws them
├── sequence/
│   ├── chirp_sweeps.png / .svg
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

Both use the parameters validated on hardware in v3: one shared sweep rate of
**7.8125 Hz/s**, a cruise amplitude of **3 full steps** (48 microsteps at
MRES 16), **400 mA RMS**, and 750 Hz sync marker bursts at each end so the
recording can be aligned to the sequence.

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

The board streams `# SEG,...` at each segment boundary, `# F,<t>,<f>,<pos>`
periodically during a sweep, and `# GAP_JUMP` when it crosses the gap. The LED
shows phase, because serial has proven unreliable to watch live on this rig:
yellow waiting, blue counting down, green running, cyan done, red failed.

**Run HIGH first.** It is entirely above the band that caused the incidents, it
answers the question you actually have, and it exercises the firmware before a
FULL run does.

---

## Two things to check before the first run

**`R_SENSE` differs between boards on this rig.** The sketch defaults to
**0.11 Ω**, matching the v6 campaign and everything else under `Stage & motor`.
The v3 chirp used **0.10 Ω** for a board marked R100, noting that its 0.11
predecessor had died. Getting it wrong mis-scales the coil current by ~10 %,
which at 400 mA is enough to change what the motor does near resonance. Read the
markings on the board that is actually fitted and set `R_SENSE_OHM` to match.

**Nothing here has run on hardware.** The sketch compiles, and its structure and
constants come from the hardware-validated v3, but this file itself has never
driven a motor. Given what the 150–250 Hz band did to this rig once already,
treat the first run as a bring-up: `CHECK` first, then `RUN HIGH`, and watch it.

# E3 findings — the first session at positions E2 never saw

Capture `E3_Experimental_Run.csv`, log `e3_terminal_20261103_110536.log`, 44.3 min, 2,661,000 samples.

## The headline

E2's model, frozen for E5 before this capture existed, was applied unchanged
to 399 positions drawn from a different seed.

| | µm rms |
|---|---|
| E3 measured error | 1.489 |
| after E2's model | **0.522** |
| E3's own fit (floor) | 0.369 |

The model removes a factor of 2.85 and explains 87.7% of the variance (r = 0.946) at positions it was never shown. The stationary
model M_B leaves 0.620 µm on the same targets, so the varying coefficients
are carrying real structure rather than fitting E2's noise.

What remains, 0.522 µm, sits above E3's own-fit floor of 0.369 µm.
In quadrature the irreducible day-to-day term is 0.369 µm — almost exactly
the floor itself, so the session-to-session change is the same size as the
within-session noise.

## Two harmonics, not four

| k | E2 amp [µm] | E3 amp [µm] | ratio | mean phase shift |
|---|---|---|---|---|
| 1 | 1.528 | 1.540 | 1.008 | -0.3° |
| 2 | 0.434 | 0.426 | 0.983 | -19.1° |
| 3 | 0.075 | 0.022 | 0.288 | -45.8° |
| 4 | 0.165 | 0.051 | 0.309 | -71.6° |

k = 1 repeats to 0.8% in amplitude and 0.3° in phase; k = 2 to 1.7% and 19°.
k = 3 and k = 4 do not repeat, and both sit at or below the own-fit residual.

Two independent lines of evidence agree:

1. Out of sample, truncating E2's model to H = 3 or H = 2 predicts E3 slightly
   *better* than the selected H = 4 (0.514 / 0.515 µm against 0.522 µm).
2. E3's own model selection, run on this capture alone, chooses H = 2, K_s = 4, λ = 3.2 (edf 32.31, CV 0.4062 µm) where E2 chose H = 4.

E2's H = 4 was session-specific structure. Within-session cross-validation
holds out positions, not sessions, so it could not see this.

## Reversal moved

| | β₀ [µm] | slope [µm/100 mm] |
|---|---|---|
| E2 | -1.757 ± 0.021 | +0.187 ± 0.061 |
| E3 | -1.623 ± 0.020 | -0.359 ± 0.058 |

β₀ changed by 0.135 µm, which is 4.7σ: small,
but not noise. The slope changed sign, so this is a change in how reversal
varies along the screw, not just in its level.

## Validity

All 13 of 13 checks pass.

- Drift datum to end: **+0.394 µm**, against 5.41 µm on 30 Oct. A clean session.
- Mirrored timing holds the mean visit time to 0.032 s across all 400 targets, so D0 and D1 agree
  to 0.005 nm. D2 moves ȳ by at most 0.196 µm.
- MSCNT is identical at all references (520) and at all
  400 map targets (264): no lost steps, and every target sits on
  one point of the microstep table, as the 64-µstep grid intends.

## Two things to carry forward

**The capture's axis was reversed.** The IDS exported this session with its sign
opposite to E2's: the stage climbed from 2.12 mm while the recorded displacement
fell. `e3_run.py` now decides the sign from the data — the wrong orientation
matches 158 of 807 stops and the right one all 807 — but the setup difference is
real and worth tracing before the next session.

**The motor's electrical phase differs between sessions.** All E2 map targets sat
at MSCNT 8, all E3 map targets at 264: a 256-count shift,
a quarter of an electrical cycle. Within a session that is a constant absorbed by
the intercept, which is why the harmonics are unaffected. Across sessions it moves
the constant — and E5 corrects absolute position, so it lands on the datum.

## What this means for E5

- The varying-coefficient model transfers. That was the open question and it is
  answered: **0.522 µm** out of sample against 1.489 µm uncorrected.
- Harmonics 3 and 4 should be dropped. They cost 7 effective parameters and buy
  nothing out of sample; E5's pre-correction currently applies them.
- E5's §6 pass lines (E ≤ 1.96, M ≤ 1.735 µm) remain doubtful. A floor of 0.522 µm
  on ȳ alone, before reversal and before quantisation, is consistent with the
  replay's predicted M of 2.04 µm and not with the 1.335 µm the lines are built on.
- β₀ moved 0.13 µm between sessions. E5 commands a fixed −β₀ from E2.

---

Figures: `rendered/findings/`. Tables: `analysis/`.
Regenerate with `python process_e3.py --select && python render_e3_findings.py`.

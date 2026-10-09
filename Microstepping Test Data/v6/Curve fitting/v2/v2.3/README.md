# v2.3 — is the ~360 Hz peak a real mode, and is it dominant?

[v2](../README.md) found a peak at 362.8 ± 3.0 Hz in the post-recording idle
tail; [v2.2](../v2.2/README.md) reproduced it at 357.9 Hz on a different
segment with a different detrend. Both measured it in the **same state**:
nothing commanded, axis at rest. A peak in one state's noise floor cannot
distinguish a mechanical mode of the stage from an artifact of the driver
holding position — or from an alias of something else entirely.

This version applies v2.2's method to the states the tail cannot reach.

```
v2.3/
├── README.md                         this file
├── spectral_summary.json             every number below, machine-readable
├── scripts/
│   └── spectral_motion_states.py     the whole analysis
└── rendered_assets/
    ├── psd_by_motion_state.png       idle vs holding vs moving, same Welch treatment
    ├── psd_ramp_grid.png             12 ramp spectra with commanded lines marked
    └── ringdown_fit.png              free decay after one microstep + damped fit
```

```powershell
cd "Microstepping Test Data/v6/Curve fitting/v2/v2.3/scripts"
python spectral_motion_states.py
```

## Short answer

| Claim | Verdict |
|---|---|
| It is a **real mechanical mode**, not an idle-state artifact | **Yes** — established here |
| It is the **dominant** frequency of the system | **No** — only in one state |
| We know its **absolute frequency** | **No** — still aliasing-limited |

## 1. It is real: four independent states agree

All four are cut from the same capture and given identical Welch treatment
(Hann, 512-sample segments, 50% overlap, cubic detrend), so their levels are
directly comparable.

| State | What it is | Peak in 330–370 Hz | vs 200–300 Hz floor |
|---|---|---:|---:|
| Idle tail | nothing commanded, after the run | 359.4 Hz | 59× |
| Oscillation dwell (MRES 1 / 4 / 16 / 32) | at rest, driver holding, mid-campaign | 349.6 / 353.5 / 359.4 / 351.6 Hz | 7× / 205× / 32× / 431× |
| Mid-ramp, fast leg | **moving** | — | 7× / 10× / 7× / 10× |
| **Ring-down after one microstep** | **free decay** | **343.8 / 350.0 / 350.0 / 343.8 Hz** | — |

**The ring-down is the decisive one.** A free decay rings at the system's own
damped natural frequency — that is a mechanical property by definition, not a
noise-floor curiosity. The same band appearing in the decay *and* while
holding *and* while moving *and* in the idle tail rules out an artifact
confined to the idle state (a chopper or holding-current effect would not be
excited by a mechanical step and then decay).

Damping, from a damped-sinusoid fit to the ensemble median decay:

| MRES | Ring-down peak | ζ | Q | τ | R² |
|---|---:|---:|---:|---:|---:|
| 1 | 343.8 Hz | — | — | — | 0.05 |
| 4 | 350.0 Hz | 0.0219 | 23 | 20.9 ms | 0.87 |
| 16 | 350.0 Hz | 0.0294 | 17 | 15.3 ms | 0.65 |
| 32 | 343.8 Hz | — | — | — | 0.55 |

ζ is quoted **only where the single-mode fit actually describes the decay**
(R² ≥ 0.60). At MRES 32 the ring amplitude is ~40 nm, close enough to the
noise that the fit reports near-zero damping without evidence. At MRES 1 the
single-mode model fails outright — see below. The trustworthy figures are
**ζ ≈ 0.022–0.029, Q ≈ 17–23, decaying in 15–21 ms**.

### MRES 1 is a different response, not a worse measurement

MRES 1's decay carries a strong **~181 Hz** component beside the ~344 Hz one,
at ~8× the amplitude of the finer-MRES rings (3908 nm vs 476 nm at MRES 4).
A single damped sinusoid cannot fit both, which is why its R² collapses.

This matters for how the fit is seeded, and the script says so: started from a
fixed 350 Hz guess, the optimiser converged onto the 181 Hz component and
reported MRES 1 as a *different mode* from the other three for the same
measurement. Seeding instead from the dominant measured peak inside
330–370 Hz makes all four measure the same feature, and the lower component is
reported separately rather than being picked silently by an optimiser.

The physical reading — offered as a hypothesis, not a conclusion — is that a
full 10 µm step is a large enough disturbance to excite the rotor's own
oscillation about its detent position in addition to the structural mode,
and that this is amplitude-dependent. Testing it needs the faster capture.

## 2. It is not dominant — that depends entirely on the state

`psd_by_motion_state.png` makes this immediate: the moving trace sits **three
to four orders of magnitude** above both at-rest traces across the whole
spectrum.

- **While moving** — commanded stepping dominates completely. The mode is
  present (7–10× the local floor on fast legs) but is a minor feature.
- **At rest** — the sub-1 Hz thermal wander carries 11.92 nm RMS against
  4.34 nm in 300–400 Hz ([v2.2](../v2.2/README.md)). The mode is the dominant
  **narrowband** feature, not the dominant contributor to position error.
- **Right after a step** — here it genuinely dominates: it *is* the ring-down,
  and it sets how long the axis takes to settle.

So "dominant/parasitic frequency" overstates it. The accurate statement is:
**a lightly damped structural mode near 350 Hz, excited by every commanded
step transient, decaying in 15–21 ms** — which is a real constraint on
settling time and the right thing to put in a simulation model, but not the
largest error source in any state except the moments just after a step.

### What actually dominates while moving

`psd_ramp_grid.png` marks two families of commanded line on all 12 ramp
spectra:

- **Harmonics of the FULL-STEP rate** (27.5 / 70 / 200 Hz) — *not* the
  microstep rate. The stage advances one full step at a time whatever MRES
  does; microstepping only subdivides how it gets there. At MRES 1 and 4 these
  run up to six orders of magnitude above the broadband floor.
- **The microstep pulse rate `rate_fsps × mres`, folded into 0–500 Hz.** Above
  Nyquist these do not vanish, they alias: 800 Hz → 200 Hz at MRES 4 fast,
  3200 Hz → 200 Hz at MRES 16 fast, 6400 Hz → 400 Hz at MRES 32 fast. The
  marked line lands on a real peak in essentially every panel, which is what
  confirms the identification.

**Why the fast legs are the clean test:** a 27.5 Hz fundamental puts a
harmonic within ~14 Hz of 350 Hz no matter what, and 70 × 5 = 350 exactly. So
the slow and moderate legs cannot separate the mode from commanded content.
Only the fast legs (harmonics at 200 and 400 Hz) leave 330–370 Hz free — and
there the band still sits 7–10× above the floor, which is the evidence that
the mode is genuinely excited during motion rather than merely buried.

## 3. The absolute frequency is still unknown

Nyquist is 500 Hz. A measured 350 Hz is indistinguishable from **650 Hz** or
**1350 Hz** at this sample rate, and every state above is sampled at 1 kHz, so
they would all alias identically — their agreement with each other is not
evidence about the absolute value. The PSD also still climbs into Nyquist,
the usual signature of folding.

What v2.3 establishes is that **the peak is a genuine free-decay frequency of
the mechanics**. What it cannot establish is **which** frequency that is.

### Seeing it four times is not four pieces of evidence about the frequency

Every state above is sampled by the same EL5101 chain at the same 1 kHz. If
the true mode were at 652 Hz, all four would fold to ~348 Hz *identically*.
The recurrence is strong evidence that **one real physical resonance is being
excited in every state** — it is not independent confirmation of the value
350 Hz, and no amount of re-measuring at 1 kHz can make it so. Decimating or
resampling the existing data cannot recover it either: the information was
lost at acquisition, so this is an information-theoretic limit rather than a
processing problem.

### Which fitted numbers survive the ambiguity, and which do not

This decides what is safe to put in a model today:

| Quantity | Robust to aliasing? | Why |
|---|---|---|
| **τ (decay time, 15–21 ms)** | **Yes** | read off the decay *envelope*, which is slow and properly sampled |
| f_d (~350 Hz) | No | the carrier is the aliased quantity |
| ζ, Q | **No** | both scale with ω_n: ζ = 1/(τ·ω_n), Q = τ·ω_n/2 |

Concretely, if the mode is really at ~650 Hz rather than ~350 Hz:

| MRES | τ | ζ at 348 Hz | Q at 348 Hz | ζ at 652 Hz | Q at 652 Hz |
|---|---:|---:|---:|---:|---:|
| 4 | 20.9 ms | 0.0219 | 23 | 0.0117 | 43 |
| 16 | 15.3 ms | 0.0294 | 17 | 0.0160 | 31 |

τ is identical either way; ζ and Q move by about 2×. **So quote the settling
time, not the damping ratio**, until the frequency is pinned down — τ is the
number that is true regardless.

This does sharpen the experiment needed, though. It no longer requires
re-running the campaign: capture a **single oscillation block at ≥5 kHz** —
a few seconds of data — and fit the ring-down. If the decay still rings at
350 Hz it is real; if it moves to 650 Hz it was an alias. That is now the
cheapest open question in the dataset, and it gates the filter design in v2,
the noise floor in v2.2, and the transient model in v1 simultaneously.

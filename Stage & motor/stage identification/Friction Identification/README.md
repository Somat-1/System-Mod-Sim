# Friction identification

The friction layer of the stage model: pre-sliding, sliding and the
reversal error that the geometric layer could only handle as one constant.

Nothing is measured here yet. This folder holds the roadmap for the
F-series and the two sessions that are planned so far.

## Layout

| Folder | What it is |
|---|---|
| [`F0b Friction Scoping/`](F0b%20Friction%20Scoping/) | the scoping session: procedure, config, renders and the sequence overview |
| [`F1 Reversal Sweep/`](F1%20Reversal%20Sweep/) | the reversal sweep: procedure, sequence generator and companion figures |
| [`Preliminary/`](Preliminary/) | the existing stage-identification material, kept as an exploratory baseline while the friction workflow is developed — not a final identification campaign |

The roadmap itself lives at this level:

| File | What it is |
|---|---|
| `ROADMAP_SPEC.md` | the render instructions the figure works to, plus what was dropped and why |
| `render_friction_roadmap.py` | draws the roadmap figure |
| `rendered/fig_friction_roadmap.png` | **the roadmap figure**, 180 × 75 mm at 600 dpi |

```sh
python render_friction_roadmap.py
```

The run prints two checks: any run of text wider than the box it sits in
(currently none), and what the fills are worth in grayscale.

## Why this layer is next

E5 closed the geometric layer. The pooled model predicting an unseen
session leaves **0.383 µm** rms against a same-day floor of **0.369 µm**, so
the whole remaining geometric model error is 0.102 µm in quadrature. A
*perfect* geometric model would move the interior compensated rms from
0.424 to 0.412 µm — three per cent.

What it would not move is what is left in the ISO table after compensation
and a two-point gain calibration: E = 2.520, M = 1.735, **B = 0.959**,
R↑ = 1.930, R↓ = 1.635 µm. B is identical before and after the gain fix,
because reversal is not a gain. The geometric layer carries reversal as one
frozen scalar, β₀ = −1.757 µm, which put this session's *mean* reversal
39 nm from zero and left 0.959 µm of position-dependent spread it has no
structure to represent. That spread, and the two repeatabilities that set
the floor, are friction's.

E5's B3 block also bounds the layer below: the motor's electrical cycle
contributes 0.30 ± 0.13 µm at 2.3σ, under the geometric residual. The next
micron is here, not in the motor LUT and not in a finer screw map.

See [`../Geometric Fitting/experimental_scripts/e5_compensated_validation/FINDINGS.md`](../Geometric%20Fitting/experimental_scripts/e5_compensated_validation/FINDINGS.md)
and [`../Geometric Fitting/experimental_scripts/geometric_model/README.md`](../Geometric%20Fitting/experimental_scripts/geometric_model/README.md).

## The roadmap figure

Three rows, read left to right: the lab sequence on top, the one test that
can run alongside it in the middle, the desk work underneath. The rows are
not labelled - position and the arrows say which is which.

```
lab        Setup → F0b → F1 → F3 → F2 ............... F5
alongside                  F4
desk       F0a ..................... Identification → Compensation
```

F0a reads what the existing captures already say, so F0b can aim at the
gaps rather than re-measuring them. F0b sizes the effect and sets the
amplitudes, speeds and dwells of everything after it. F1 is the central
experiment: the reversal sweep, which identifies the pre-sliding elements
and δ_C, with the old single-step test folded into it as slow 1/256 ramps
at the small-amplitude end. F3 gives the Stribeck curve. F2 runs only if
F0b finds a dwell effect above the noise. F4 adds the force scale and the
guide/nut split and is optional. Identification fits the GMS model with a
hold-out check; compensation turns it into a host-side command table; F5
validates it over small and large moves in four modes — geometry only,
constant reversal, the friction model, and a unidirectional approach.

### Caption

> **Fig. X.** Friction identification roadmap. A desk review of existing
> data (F0a) and a short scoping session (F0b), including the measurement
> noise floor, set the amplitudes, speeds and dwells of the main
> experiments. The reversal sweep (F1, central) identifies pre-sliding, the
> velocity sweep (F3) the sliding behaviour, and the dwell test (F2) is run
> only if scoping shows a dwell effect above the noise. The guide-friction
> test (F4) runs in parallel and adds only the force scale. Identification
> and compensation produce the command table for the validation (F5), which
> compares geometry only, constant reversal, the friction model and a
> unidirectional approach.

### What the figure deliberately does not say

The first draft followed `ROADMAP_SPEC.md` literally: five lines in every
box, three decision diamonds with their questions, a feed-forward arc over
the whole row, a dashed bypass branch, and the route-B footnote. All of it
was accurate and none of it was readable — a roadmap that needs the thesis
open beside it is not doing its job.

So the structure stayed and the apparatus went. Each box now carries its
name and one plain sentence. The gates are gone; the only branch a reader
has to hold is F2's, and it is written above the box in words. Terms a
first-time reader would have to look up — pre-sliding elements, δ_C, the
Stribeck curve, GMS, route B — are in the paragraph above and in the
caption, where they can be read at the reader's own pace rather than
decoded from a diagram. `ROADMAP_SPEC.md` §7 lists every item dropped.

## Carry these into the sessions

Three things from the geometric close-out cost almost nothing if they ride
along with the F-series, and will not happen otherwise.

1. **The daily two-point block.** R-pair at 15 and 105 mm, both directions,
   ×2, after warm-up — two minutes at the start of every session, with a
   BME280 and two NTCs (one on the nut, one on the motor flange, ≈ €15).
   Over one to two weeks of November weather this decides whether the
   +11 ppm session gain is air, screw temperature or cosine error, and that
   verdict changes what the deployed stage needs. See `CLOSEOUT.md` §1.4,
   which already assumes these sessions as its host.
2. **Log the serial port.** E5 lost four validity checks — MSCNT across U
   and R stops, MSCNT against `el_off`, the echoed `u_cmd`, and the
   firmware summary — for want of a file that costs nothing to keep. Start
   the IDS logger early and stop it late as well.
3. **B3 re-analysis**, about an hour at the desk, whenever. It tightens the
   bound on the motor layer; it will not change the plan.

## Status

Nothing is measured. What exists is the roadmap and two planned sessions.
F0b scopes the effects and predicts 48.2 min against the 58 min watchdog
gate; F1 is the reversal sweep it feeds. `data/` and `analysis/` arrive
when the sessions are run; F2-F5 get folders of their own beside them.

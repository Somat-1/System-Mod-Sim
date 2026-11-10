# F1 — reversal-amplitude sweep

F1 is the central friction experiment. F0b measured how large the friction
effects are; F1 measures their shape. It covers virgin curves (V), symmetric
loops (P), nested loops (M) and a rate check (Q), at three positions along the
travel, so that the GMS pre-sliding elements can be identified and checked for
independence from position. The procedure is
[`F1_reversal_sweep_procedure.md`](F1_reversal_sweep_procedure.md).

**Status: sequence planned and overview rendered. There is no command table,
firmware or data yet.** The values marked ⟨…⟩ in the procedure depend on F0b
(s₉₀, d_min). [`f1_config.json`](f1_config.json) carries the procedure's
defaults, s₉₀ = 3 µm and d_min = 5 nm, under `f0b_dependent`. Replace those two
numbers once F0b is processed, then rerun both scripts.

## Files

| Path | Contents |
|---|---|
| `F1_reversal_sweep_procedure.md` | the supplied procedure, including the render instructions (§8) |
| `f1_config.json` | positions, stroke rule, speeds, pauses and colours; the F0b-dependent values sit under `f0b_dependent` |
| `f1_sequence.py` | builds the timed plan: every move, settle, turning pause and reference stop as one segment, with positions rounded to 1/256 µsteps |
| `sequence/f1_segments.csv` | the plan, one row per segment, with tags for group, sub-block, stroke, cycle and half |
| `sequence/f1_plan_summary.json` | segment count, duration, the stroke list as rounded, and time per group, sub-block and segment kind |
| `render_f1_sequence.py` | the overview figure, `rendered/f1_sequence.png` |

```sh
python3 f1_sequence.py
python3 render_f1_sequence.py   # needs the E5 frozen geometric model for A1(x)
```

## Plan with the default values

The plan has **1107 segments and lasts 51.8 min**. R2 is included, because the
procedure drops it only above 52 min.

| Group | Contents | Time |
|---|---|---|
| R1 | reference pair | 2.7 min |
| A (25 mm) | V · P (increasing) · M | 12.1 min |
| B (60 mm) | V · P (increasing and decreasing) · M · Q | 21.4 min |
| R2 | reference pair | 2.7 min |
| C (100.48 mm) | V · P (increasing) · M | 10.8 min |
| R3 | reference pair | 2.0 min |

These figures differ from the estimates in procedure §4:

- **B runs 21.4 min against an estimate of about 17.** Its two P passes take
  14.6 min. Each stroke also needs a move to x_p − s/2 and a move back, and the
  ≈ 6 min estimate per pass leaves those out.
- **The reference pairs are shorter than the 3.5 min estimate.** They take
  2.0–2.7 min, which is four stops with 3 mm approaches and transits at
  1 mm/s. This offsets most of B's overrun.
- **Margin:** the total sits only 0.2 min under the 52 min gate. A larger s₉₀
  from F0b raises the V, M and Q times, so if the regenerated plan exceeds the
  gate, R2 is dropped automatically.
- **Transit time:** moves between positions (about 1.9 min) are counted in the
  first wipe of the next group. That is why A's V block takes longer than C's.

## Reversal and continuation checks

Every segment stores `last_dir`, the direction of the last motion before it
(+1 up, −1 down). Declared measurement ramps also carry a `ramp_type` of
`reversal` or `continuation`. A wipe from below ends moving up, and a wipe
from above ends moving down; `wipe()` asserts this.

The ramp type is checked twice:

- **While the plan is built:** `Plan.ramp()` asserts `ramp_dir == -last_dir`
  for a reversal and `ramp_dir == +last_dir` for a continuation.
- **After the plan is written:** `check_ramp_types()` recomputes every
  declared type from the positions in the plan.

The virgin ramps (V) are declared as reversals for their slow first part and
as continuations for the 2 µm/s remainder. An earlier version of this plan
drew the V ramps in the wipe direction, which made them continuation ramps.
Both checks now reject that error. Carry the same checks into the F0b and F1
command-table generators.

## Assumptions not fixed by the procedure

- **Start position:** the session starts homed and centred at x_B. Each
  reference pair visits its nearer reference first.
- **M branch speeds:** each branch follows the 4 s half-cycle rule of its own
  loop. Outer branches run at S_o/4 s and the inner loop at S_i/4 s, which
  makes the inner loop comparable to a P loop of the same stroke. Under this
  rule M takes about 0.75 min per position, against the ≈ 1.5 min estimate.
- **Stroke rounding:** strokes are rounded symmetrically, as a whole number of
  µsteps on each side of x_p. With the defaults, the 0.15 µm stroke becomes
  0.156 µm (4 µsteps) and is the only stroke run as `step` segments.

## Figure

`rendered/f1_sequence.png` follows procedure §8.

- **(a)** A two-tier ribbon shows the position groups and the sub-blocks, with
  widths proportional to time. Below it is the commanded position relative to
  the current test position, on a symmetric log scale. P and Q appear as light
  envelopes with thin outlines. Reference pairs are purple dashed, and
  transits are grey dashed with an "→ position" label.
- **(b)** Amplitude coverage for each sub-block, with guides for d_min, one
  µstep, the s₉₀ band, β₀ and the electrical cycle.
- **(c)** A travel map showing the test positions with their screw angles, the
  reference positions, the weakly constrained zones and the known features. A
  strip above it shows the first-harmonic amplitude A₁(x) from the frozen E2
  geometric model (`geom_model_E5.json`).

The caption is embedded in the figure, and its segment count and duration come
from the plan.

## Companion figures

`render_f1_companions.py` writes the static companion figures to
`rendered/companion/`.

| File | Shows |
|---|---|
| `f1_V1_virgin_curve.png` | V1: the virgin-curve anatomy. Reversal curves in both directions on a log travel axis, with the stuck, pre-sliding and sliding regimes, s₉₀, the slow/fast ramp strip and a sketch of the V command |
| `f1_P1_loops_and_width.png` | P1: loops at 0.3 s₉₀, s₉₀ and 100 µm on a normalised position axis, and loop width against stroke at the 14 F1 strokes |
| `f1_P2_masing_overlay.png` | P2: a loop with a 6 µm stroke, overlaid with 2·f(Δu/2) from each turning point, and a residual strip against ±d_min |
| `f1_M1_nested_loop.png` | M1: the outer and inner loop against the uninterrupted outer branch, with a zoom on the branch-off and rejoin point |
| `f1_Q1_rate_check.png` | Q1: the s₉₀ loop at half-cycle times of 1.3, 4 and 12 s plus the repeat, and loop width against half-cycle time with a ±d_min band |

All the curves are illustrative and are labelled as such. They come from the
quasi-static Maxwell-slip model used for the F0b S animation, imported from
`../F0b Friction Scoping/render_f0b_s_animation.py`, with placeholder slip
limits and a saturation of β₀/2 per side. The model's own s₉₀ is 2.92 µm,
close to the 3 µm placeholder.

This model is rate-independent, obeys the Masing rule exactly and has exact
return-point memory. The overlays in P2, M1 and Q1 therefore coincide by
construction, and the script asserts both the Masing match and the zero
return-point error. These figures show what a GMS-type stage would produce;
they are not evidence for one.

X1 (position comparison and hold-out) is data-only and is not drawn. The
animations V2, P3 and M2 are not built yet.

## Running the session

> **The table is provisional.** The P, M and Q strokes are built from s₉₀ and
> d_min, which F0b measures. Until F0b is processed, `f1_config.json` carries
> the procedure's placeholder defaults (s₉₀ = 3 µm, d_min = 5 nm), the
> generated header sets `CFG_PROVISIONAL = 1`, and the firmware says so at
> boot and in the session header. **Re-run `f1_make_table.py` after F0b and
> reflash before treating an F1 run as final.**

`f1_sequence.py` plans the session in time, for the overview figure.
`f1_make_table.py` builds the table the board replays, in the same five
segment types as F0b:

```bash
python3 f1_make_table.py             # writes sequence/ and the firmware header
```

It reports **1119 segments, 51.8 min** with R2 included — which reproduces
`f1_sequence.py`'s independent 51.78 min and its per-group times to within
0.01 min. Two separately written builders agreeing is the check that the
sub-block definitions were read the same way twice.

**The margin is 0.2 min.** The procedure's gate is 52 min and the plan lands
at 51.8. `f1_make_table.py` drops R2 automatically if the session would
exceed the gate, exactly as `f1_sequence.py` does, but any growth in the
stroke grid when F0b's real s₉₀ arrives will trip that. Expect to run without
R2, or to shorten P, once the grid is regenerated.

| File | What it is |
|---|---|
| `f1_make_table.py` | builds the table, times it, checks it, writes the header |
| `sequence/f1_command_table.csv` | the table |
| `sequence/f1_table_summary.json` | durations, counts, hashes and the checks |
| `firmware/f1_reversal_sweep/f1_config.h` | generated; do not edit |
| `firmware/f1_reversal_sweep/f1_reversal_sweep.ino` | the sketch |

The segment types, the two step paths and the bench procedure are the same
as F0b's — see
[`../F0b Friction Scoping/README.md`](../F0b%20Friction%20Scoping/README.md).
F1's slowest measurement ramp is 2 Hz, at the smallest loop stroke, so it
leans on the timed path more than F0b does. The session parks at
x_B = 60.00 mm, the position F0b also uses.

# F1 — reversal-amplitude sweep

F1 is the central friction experiment. F0b measured how large the friction
effects are; F1 measures their shape. It covers virgin curves (V), symmetric
loops (P), nested loops (M) and a rate check (Q), at three positions along the
travel, so that the GMS pre-sliding elements can be identified and checked for
independence from position. The procedure is
[`F1_reversal_sweep_procedure.md`](F1_reversal_sweep_procedure.md).

**Status: run on 10 November 2026 and complete — 1119/1119 segments in
51.31 min, not aborted or time-limited.** The capture and the ESP terminal
log are in [`data/`](data/); the raw CSV is local-only, as every hardware
export in this repository is. A first look at it is in
[`../analysis/README.md`](../analysis/README.md).

**The run used the provisional stroke grid.** The values marked ⟨…⟩ in the
procedure depend on F0b (s₉₀, d_min), and
[`f1_config.json`](f1_config.json) still carries the procedure's defaults,
s₉₀ = 3 µm and d_min = 5 nm, under `f0b_dependent`; the firmware said so twice
in the log. The preliminary check puts the real s₉₀ near 7.5 µm, so the dense
part of the grid sits about a factor of two low. Replace those two numbers
once F0b is properly processed, rerun both scripts, and decide then whether
this loop set stands or the sweep is repeated on the corrected grid.

## Files

| Path | Contents |
|---|---|
| `F1_reversal_sweep_procedure.md` | the supplied procedure, including the render instructions (§8) |
| `f1_config.json` | positions, stroke rule, speeds, pauses and colours; the F0b-dependent values sit under `f0b_dependent` |
| `f1_sequence.py` | builds the timed plan: every move, settle, turning pause and reference stop as one segment, with positions rounded to 1/256 µsteps |
| `sequence/f1_segments.csv` | the plan, one row per segment, with tags for group, sub-block, stroke, cycle and half |
| `sequence/f1_plan_summary.json` | segment count, duration, the stroke list as rounded, and time per group, sub-block and segment kind |
| `render_f1_sequence.py` | the overview figure, `rendered/f1_sequence.png` |
| `render_f1_companions.py` | the shared illustrative solver and the V1 figure, as PNG |
| `f1_gif_common.py` | colours, canvas, palette and GIF writer shared by the three animations |
| `render_f1_P_gif.py` | the P animation, `rendered/companion/f1_P_loops.gif` |
| `render_f1_M_gif.py` | the M animation, `rendered/companion/f1_M_memory.gif` |
| `render_f1_Q_gif.py` | the Q animation, `rendered/companion/f1_Q_speed.gif` |
| `check_f1_companions.py` | the physical invariants of the shared illustrative solver |

```sh
python3 f1_sequence.py
python3 render_f1_sequence.py   # needs the E5 frozen geometric model for A1(x)

python3 render_f1_companions.py # the V1 figure
python3 render_f1_P_gif.py      # the three animations
python3 render_f1_M_gif.py
python3 render_f1_Q_gif.py
python3 check_f1_companions.py  # the solver they all share
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

Run python render_f1_companions.py to write the V1 PNG to
rendered/companion/, and the three GIF renderers to write the animations
beside it. The supplied design briefs are preserved in
[COMPANION_DESIGN.md](COMPANION_DESIGN.md),
[F1_P_loops_gif_brief.md](F1_P_loops_gif_brief.md) and
[F1_M_Q_gifs_brief.md](F1_M_Q_gifs_brief.md), each with a record of where
the build departs from it.

| Figure | Reading sequence |
|---|---|
| V1 | Wipe and reversal sketch, then positive lost motion with nearly-stuck, pre-sliding and sliding regions; a grey direction-dependent alternative |
| **P · GIF** | **What the motor is told to do, what the stage does, and the gap for each size of move: 1, 3 and 100 µm, drawn live, then summarised** |
| **M · GIF** | **A big move with a small detour halfway back, and a zoom on the moment it rejoins the path it left; a stage without memory would not** |
| **Q · GIF** | **The same small move done fast, medium and slow, each loop tracing over the last, and the gap for each speed** |

The three animations are meant to be shown in that order — P, M, Q — which
is how the brief sequences the three plain-language points: the gap depends
on the size of the move, the stage remembers its past, and speed does not
change it.

All plots are explicitly illustrative, not measured. The hardware sequence and
f1_config.json are unchanged. No processed F0b result was available here;
s90 = 3 um and d_min = 0.005 um remain configuration placeholders.

### The animations

Four of the five technical companion PNGs are now animations: P1 (loops and
width) and P2 (Masing overlay) became `f1_P_loops.gif`, M1 (nested loop)
became `f1_M_memory.gif`, and Q1 (rate check) became `f1_Q_speed.gif`. Those
four PNGs are deleted, as are the PDF copies of every companion figure: the
renderer writes PNG only, and `render_f1_companions.py` now holds the shared
solver and V1. What the three animations have in common lives in
`f1_gif_common.py`.

| Animation | Frames | Length | Size |
|---|---|---|---|
| `f1_P_loops.gif` | 573 | 28.6 s | 0.67 MB |
| `f1_M_memory.gif` | 390 | 19.5 s | 0.42 MB |
| `f1_Q_speed.gif` | 490 | 24.5 s | 0.44 MB |

All three are 1000 × 560 px at 20 fps, loop with a pause on the finished
picture, and are written for a reader with no background: no symbols, no
equations, numbers only where they help. Pillow merges identical consecutive
frames as it writes, so each file holds fewer stored frames carrying the same
total duration.

What survives from the technical figures is one panel each: P1's gap against
the size of the move, M1's return-point detail, and Q1's nested line weights.
What is no longer drawn is P2's Masing construction, Q1's ±d_min band and its
drift-check ring. The invariants behind them are still asserted by
`check_f1_companions.py`, whose solver all three animations share, and the
drift check is still in the Q block itself and in the procedure.

Each renderer asserts the claim its animation makes before a frame is drawn:

- **M**: the detour returns to the lag it left, and what follows it is the
  branch the big move was on rather than something close to it — **0 µm** and
  **4 × 10⁻¹⁶ µm** respectively. The grey alternative has to miss by enough to
  be worth drawing (0.090 µm) and then close most of it.
- **Q**: the three loops agree to **3 × 10⁻¹⁵ µm** at coincident commanded
  positions, and the three gaps are the same floating-point number,
  0.904447 µm, which is also `width_at_centre(3 µm)`.
- **P**: the two drawn cycles of each move agree to 1e-9, and each move's gap
  lands on the summary curve to better than 2 nm.

Choices worth knowing about, all recorded in the briefs:

- **Each move is settled before it is drawn.** The first back-and-forth after
  a change of stroke inherits the state the previous one left: for the 1 µm
  move its gap is 0.475 µm against the settled 0.037 µm. Drawing it would
  contradict the message on screen and make the gap arrow collapse halfway
  through the move, so one cycle is run through the solver before the two
  that are shown — as the bench does, with three cycles per stroke and the
  later ones read. The footer says so, and the two drawn cycles are asserted
  to agree to 1e-9.
  The same settling is run for M's big move and for Q's three passes.
- **The solver is imported, not re-implemented.** All three renderers take
  `run`, the element limits and the weights from `render_f1_companions.py`,
  so the animations and V1 cannot drift apart. The P brief's own §5 snippet
  is explicitly a placeholder and is not used.
- **M's "stage that forgets" is the same solver with the outer excursion
  erased**, not the exponential its brief sketches. That formula needs the
  lag at the end of the detour to differ from the big path's, and with
  return-point memory those are the same number, so it would draw nothing.
- **Q's panel (A) is drawn against animation time, not real time.** Against
  real time the slow pass would be twelve times as wide and its dot would
  cross the panel three times faster than the fast one's, which is backwards.
  The labels carry the real 1 s / 4 s / 12 s and the footer states that the
  drawing stretches time 1 : 2 : 4.
- **Q's gap panel starts at zero**, so the 10 percent rate-dependent
  alternative looks like the 10 percent it is rather than like a large
  effect.

### Model and technical choices

The new self-contained solver advances 30 Maxwell-slip elements exactly from
one slip event to the next. Limits are log-spaced from 0.01 to 20 um. The model
satisfies equilibrium, saturation beta0 = 1.757 um, sum(gamma) = 20, initial
lost-motion slope 20/21, and the illustrative s90 = 3 um.

The design brief's fixed log-normal centre near s90/3 and width near 1 cannot
also satisfy all these constraints. Both log-normal parameters are therefore
solved numerically instead; the actual values are recorded in
rendered/companion/illustrative_model.json. This is a broad, truncated
illustrative distribution, not an inferred friction-element distribution.

Masing's f is the initial-loading response from unstrained elements. The
saturated reversal response in V1 already contains the factor-two
transformation and is not doubled again. Small-loop vertical bias is retained
from the wipe history, rather than artificially recentering the loops; this
is why the small loop in the P animation sits low on its axis instead of in
the middle.

The M animation enlarges the detour to 3 um against a 10 um move, as its
brief asks, which is larger than the planned test; the footer says so. One
consequence: halfway down a 10 um move the stage is 5 um past the reversal
and s90 is 3 um, so the branch it rejoins is already most of the way to
saturation and only gently curved, which is why the rejoin is argued in the
zoom panel rather than left to the eye in the loop. Its no-memory example is
one way to lack memory, not a claim about every model. The Q rate ghost is
deliberately illustrative, 10 percent from fast to slow.

Run python check_f1_companions.py to check reversal symmetry, saturation,
equilibrium, loop-width bounds, Masing and return-point memory. The renderer
also checks equilibrium on every incremental solve. Those checks cover all
three animations, which use the same solver.

X1 remains data-only.

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

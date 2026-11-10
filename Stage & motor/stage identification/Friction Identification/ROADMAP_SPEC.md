# Friction identification roadmap: render instructions (v2)

A process figure for the thesis (no dates), with an optional dated variant for supervisor meetings. Style matches the E3/E5 figures.

---

## 1. Structure

Left-to-right flow in **three horizontal lanes**:

| Lane (top → bottom) | Boxes |
|---|---|
| **Lab sessions** | Setup → F0b → F1 → F3 → F2 → … → F5 |
| **Parallel / optional** | F4 |
| **Desk work** | F0a → … Identification → Compensation |

Changes from v1:
- **F0 is split** into a desk review of existing data (**F0a**, desk lane, left) and a short scoping session (**F0b**, lab lane), which now includes the noise floor.
- The original single-step test is **folded into F1** (slow 1/256 ramps at the small-amplitude end).
- **F2 is conditional:** drawn with a dashed border, entered only if F0b shows a dwell effect.
- A **feed-forward arrow** from F0b to F1, F3 and F2, labelled "tailors amplitudes, speeds, dwells".

---

## 2. Layout

- **Figure size:** about 7.2 × 3.6 in (text width); 16:9 for slides.
- **Lane labels** on the far left, small caps or rotated 90°: *lab sessions · parallel · desk*.
- **Lab lane (top):** seven positions, equal box widths (~0.9 in), ≥ 0.25 in gaps:
  `Setup → F0b → F1 → F3 → F2 → [gap] → F5`
  Leave a wider gap between F2 and F5, where the desk work sits underneath.
- **Desk lane (bottom):**
  - **F0a** at the far left, under Setup/F0b. A solid arrow up into F0b ("magnitudes known → scoping targets the gaps").
  - **Identification → Compensation** under the F2–F5 gap.
- **Parallel lane (middle):** **F4** as one box spanning roughly under F1–F3.

### Arrows

| From → to | Style | Label |
|---|---|---|
| Lab sequence (Setup → F0b → F1 → F3 → F2) | solid | — |
| F0a → F0b | solid, upward | "known magnitudes" |
| F0b → F1, F3, F2 | thin dashed, curved above the lab lane | "tailors amplitudes, speeds, dwells" (one shared label) |
| F1, F3, F2 → Identification | dashed, downward | "data" (once) |
| F4 → Identification | dashed | "force scale (optional)" |
| Identification → Compensation | solid | — |
| Compensation → F5 | solid, upward | "command table" |

### Decision gates (small diamonds on arrows; text beside, not inside)

| Position | Text |
|---|---|
| After F0b | "detectable? pre-sliding length?" |
| After F1 | "loops overlap across positions?" |
| Before F2 | "dwell effect > noise?" (if no: skip F2, dashed bypass arrow F3 → Identification) |
| Inside / under Identification | "route-B: refit Ĝ unchanged?" (small text, no diamond) |

---

## 3. Box contents

Three lines per box: **title** (bold, 9 pt), question (8 pt), output (7 pt, grey).

| Box | Title | Question | Output |
|---|---|---|---|
| Setup | **Setup** | 1/256, logging, references | — |
| F0a | **F0a · desk review** | what do existing data already show? | noise σ(τ), β₀, dither floor, creep |
| F0b | **F0b · scoping** | how big, how far, how slow? | noise floor, pre-sliding length, deadband |
| F1 | **F1 · reversal sweep** | lost motion vs distance? | pre-sliding elements, δ_C, small-step ramps |
| F3 | **F3 · velocity sweep** | lag vs speed? | Stribeck curve |
| F2 | **F2 · dwell** *(if needed)* | does rest change the next move? | static rise |
| F4 | **F4 · guide friction** | force of the guides alone? | force scale, guide/nut split |
| Identification | **Identification** | fit GMS, hold-out check | parameters (µm units) |
| Compensation | **Compensation** | host-side friction prediction | command table |
| F5 | **F5 · validation** | small + large moves, 4 modes | friction-chapter result |

Mark **F1** as the central experiment with a thicker border or a small "central" tag above the box.

---

## 4. Visual style

| Element | Fill | Border | Text |
|---|---|---|---|
| Lab boxes | light blue #e8f0f5 | dark blue #1f4e79, 1 pt | dark |
| F1 (central) | light blue #e8f0f5 | dark blue #1f4e79, 2 pt | dark |
| F2 (conditional) | light blue #e8f0f5 | dark blue, **dashed** | dark |
| F5 (goal) | dark blue #1f4e79 | — | white |
| F4 (parallel) | light ochre #f5ece0 | ochre #a0692a | dark |
| Desk boxes (F0a, Identification, Compensation) | light grey #f1f2f4 | grey #7a8088 | dark |
| Gates | white | purple #6c5a8e | purple, beside the diamond |

- **Lines:** solid = sequence; dashed = data or optional flow. Small filled arrowheads.
- **Type:** serif, as in the other figures. No figure title.
- **Rounded corners**, radius ~0.05 in.
- **Black-and-white safe:** the lanes are separated by position; solid vs dashed lines carry flow type; F5 stands out by its filled box; F2 by its dashed border. Colour is never the only cue.

---

## 5. Caption

> Friction identification roadmap. A desk review of existing data (F0a) and a short scoping session (F0b), including the measurement noise floor, set the amplitudes, speeds and dwells of the main experiments. The reversal sweep (F1, central) identifies pre-sliding, the velocity sweep (F3) the sliding behaviour, and the dwell test (F2) is run only if scoping shows a dwell effect above the noise. The guide-friction test (F4) runs in parallel and adds only the force scale. Identification and compensation produce the command table for the validation (F5), which compares geometry only, constant reversal, the friction model and a unidirectional approach.

---

## 6. Optional dated variant (planning)

- Add a thin strip under the desk lane with week ticks **W46 … W49**.
- Shade a narrow band at **~20 Nov** labelled "flexure decision".
- Lab boxes aligned approximately to their planned week; desk boxes span their weeks as horizontal bars.
- Use this version for supervisor meetings; keep the undated version for the thesis.

---

## 7. Implementation record

Everything above is the instruction as received. This section records what
the renderer actually does, and where it departs.

The figure is drawn in millimetres at final size, 1:1, with the same
primitives as `Stage & motor/paper_figures/render_identification_scheme.py`
— orthogonal arrows, one small filled head, serif type, and an automatic
check that no run of text is wider than the box it sits in. Two sizes:
180 × 75 mm, rendered to PNG at 600 dpi.

### 7.1 Simplified for a general reader

A first version followed §2–§4 literally. It was accurate and unreadable:
five lines of text in every box, three decision diamonds each with a
two-line question, a feed-forward arc spanning the lane, a dashed bypass,
and the route-B footnote. The brief for the figure changed to "simple
enough for anyone", so the structure was kept and the apparatus removed.

| Spec item | What happened |
|---|---|
| §3, three lines per box (title, question, output) | **Two.** Title and one plain sentence. The outputs — pre-sliding elements, δ_C, the Stribeck curve, the GMS parameters, the command table — moved to `README.md` and the caption. |
| §2, decision gates (three diamonds + questions) | **Dropped.** The only branch a reader must hold is F2's, and it is written above the box as *run if F0b finds an effect* - a criterion for running it, not a caveat. |
| §2, feed-forward arc F0b → F1, F3, F2 | **Dropped.** It crossed the whole row to say something the left-to-right order already says, and it was the single densest object on the page. The README states it. |
| §2, dashed bypass F3 → Identification | **Dropped** with the gate it branched from. F3 already feeds Identification; the bypass only restated the skip. |
| §2, route-B note under Identification | **Dropped.** Meaningless without the geometric chapter; it is in the README. |
| §3, "central" tag on F1 | **Kept, in words.** Heavier border *and* the note *the main experiment* above it. A border weight is not something a first-time reader decodes. |
| §3, box wording | **Rewritten in plain language.** "lost motion vs distance?" became "what changes when the stage reverses?"; "fit GMS, hold-out check" became "turn the measurements into a friction model". |
| §2, lane labels | **Dropped.** They were first renamed *IN THE LAB · IN PARALLEL · AT THE DESK* and set in a left gutter, then removed: the rows are already distinguished by what is in them and by where the arrows go, and the gutter cost 2.5 mm the drawing now uses. The three rows are named in the caption instead. |

### 7.2 Smaller departures

- **Type sizes.** §3 asks for 9 / 8 / 7 pt in boxes about 0.9 in wide. At
  9 pt bold, "F5 · validation" is 23 mm in a 20 mm box, so the spec's own
  two numbers cannot both hold. Type runs at 7.6 / 7.0 pt, with 6.4 pt for
  edge labels. Box width 23 mm, gaps 5.6 mm (§2 asks ≥ 6.35; the width went
  to the text instead).
- **Curved feed-forward.** Moot now that the arc is gone, but every
  connector in this figure and in the companion scheme is orthogonal, and
  one curve among them would read as a different kind of relation.
- **F4's hand-off** enters the side of Identification, not its top: the top
  already carries the measurement bus, and two parallel verticals into one
  box edge put the ochre label across the dashed one.
- **Colours** are §4's table unchanged. Its three pale fills sit within
  five grey levels of each other, which elsewhere would be a fault — here
  §4's own black-and-white clause applies, since row position separates
  the families and solid-against-dashed separates the flow types. The fills
  are decoration, not an encoding. The one fill that must survive a grey
  print is F5's, at luminance 71 against 239. `grayscale_report()` prints
  all of them on every run.

### 7.3 The dated variant (§6) — dropped

A dated variant was built: a week tag in each box and a W46–W49 ruler with
the flexure band under it. It is no longer rendered, and the renderer no
longer carries the code for it.

Two reasons. Its week assignment was never data — §6 names W46–W49 and a
flexure decision near 20 November but gives no per-experiment plan, so the
schedule in it was an assumption. And §6 asks for boxes aligned to their
week *and* equal box widths (§2), which cannot both hold: making x mean
time breaks the row layout. A schedule is better kept where it can be
edited — it is not a figure.

Only `rendered/fig_friction_roadmap.png` is produced. PDF output was
dropped with it; the drawing is still built in millimetres at final size,
so the PNG is 180 mm wide at 600 dpi.

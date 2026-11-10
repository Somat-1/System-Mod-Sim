# F1 companion plots: redesign instructions

The current plots are correct, but they show only the abstract friction quantity (lag against position). The reader can't see **what the stage was doing**, **which part of the motion** a curve belongs to, or **what a different outcome would look like**. The redesign fixes that with five rules, then applies them plot by plot.

---

## 0. Design rules (apply to every plot)

1. **Pair every friction plot with a mini command icon.** A small grey u-vs-time sketch of the motion (wipe, ramp, loop) with the measured segment highlighted in the block colour. Numbered markers (①, ②, ③ …) on the icon match the same numbers on the friction curve, so the reader can follow the motion from command to measurement.
2. **Use the most intuitive quantity.** Where possible, plot **lost motion** ℓ = commanded travel − carriage travel since the last reversal (positive, starts at 0, grows, saturates at β₀) instead of the signed lag e. Use the signed lag only for loops, where both branches have to be shown.
3. **Show the alternative.** Every test has a "pass" and a "fail" picture. Draw the expected outcome in the block colour, and the alternative (what would be seen if the hypothesis were wrong) as a faint grey dotted ghost, labelled. That makes the test visible even before data exist.
4. **Smooth model curves.** The current kinks come from using only 5 friction elements. Use about 30 elements so the illustrative curves look like physical behaviour (parameters in §6).
5. **Minimal in-plot text.** Direct labels on curves, numbered points, short tags. Explanations go in the caption. Arrows on loop branches large enough to read (≥ 4 pt).

Colours (as in the F1 overview): V slate #2b3640, P teal #1a6478, M ochre #c48a33, Q red #b5473a, positioning grey #b8bdc2, ghosts grey #9aa0a6 dotted. Black-and-white safe: line style plus numbering carries the meaning, colour is secondary.

---

## 1. V1 — Virgin curve (what a reversal does)

**Current issues:** the two mirrored lag curves cross in the middle, which is visually busy; the sign convention (lag going from −0.88 to +0.88) isn't intuitive; the kinks look synthetic.

**Layout:** one main panel (b) with a small command icon (a) at the top left, inside or beside it.

**(a) Command icon** (~1.5 × 1 in):
- u against time: grey wipe (down 3 mm, back up), then the slate reversal ramp going down.
- ① end of wipe (last motion up), ② reversal starts, ③ end of slow part, ④ end of ramp.
- Tag the two speed parts: "slow 0.25 µm/s" and "2 µm/s".

**(b) Main panel:**
- **x:** commanded travel since reversal s (log, 0.01–100 µm).
- **y:** **lost motion** ℓ(s) = s − Δx (µm), 0 to 2.
- **Curves:** reversal down (solid slate) and reversal up (dashed slate). They should lie on top of each other, which is itself the symmetry check. Label: "down and up coincide → no geometric leak".
- **Guides:**
  - faint grey line ℓ = s, labelled "carriage stuck: all motion lost";
  - dotted horizontal at **β₀ = 1.76 µm**, labelled "sliding: fixed loss";
  - vertical dashed line at **s₉₀**, labelled "pre-sliding length".
- **Regime bands** under the curve (light shading, labels at the bottom): "stuck" (up to the deadband) · "pre-sliding" (deadband to s₉₀) · "sliding" (beyond s₉₀).
- Place ②, ③, ④ on the curve where they occur.

**Optional (c), small:** transmission τ = Δx/s (0 to 1), the same information as "fraction of the command that arrives". Use only if space allows.

**Caption idea:** "After a reversal the carriage first stays put (all commanded motion lost), then follows partly as the friction contacts deform and slip, and finally slides with a fixed loss β₀. Both directions give the same curve."

---

## 2. P1 — Loops and loop width (what a stroke does)

> **Superseded.** P1 and P2 were replaced by one plain-language animation,
> `rendered/companion/f1_P_loops.gif`, built by `render_f1_P_gif.py` to the
> brief in [`F1_P_loops_gif_brief.md`](F1_P_loops_gif_brief.md). The two
> PNGs and all the companion PDFs were deleted. The sections below are kept
> as the record of what the technical P figures showed: the stroke-to-width
> relationship survives as panel (C) of the GIF, while the Masing check of
> §3 is no longer drawn — the invariant itself is still asserted by
> `check_f1_companions.py`.

**Current issues:** in the normalised overlay the small loop collapses into a line; three loop types on one axis are hard to separate; it isn't clear where a loop starts or what "width" means.

**Layout:** three small multiples on the left (one per stroke), one summary panel on the right.

**Small multiples (a1–a3), shared y, own x:**
- Strokes: small (≈ 0.3·s₉₀), medium (≈ s₉₀), large (100 µm).
- **x:** commanded position relative to the centre, in **absolute µm** (each panel its own range), so the real size of each stroke is visible.
- **y:** signed lag e (µm), shared range ±1.
- **Loop:** teal, with large arrows on the branches; ① start (bottom-left turning point), ② top turning point, ③ back.
- **Width bracket:** a vertical double arrow at the centre position, labelled **w**, with the definition once in (a1): "gap between going up and going down at the same commanded position".
- **Command icon** above each panel: a triangle wave of that stroke, 3 cycles, the same numbering.
- Short tag in each panel: "inside pre-sliding", "at s₉₀", "sliding".

**Summary (b):**
- **x:** stroke s (log), **y:** loop width w (µm), as now (this panel works well).
- Mark the three strokes of a1–a3 with letters on the curve, so the reader connects them.
- Shade the s₉₀ band; dotted line at β₀.
- **Ghost alternative:** a dotted grey step at β₀ (the "constant reversal" assumption β₀ used in E5), labelled "constant correction assumes this". The gap between the ghost and the curve at small strokes is exactly where the constant correction over-corrects.

**Caption idea:** "Small strokes produce thin loops because most of the motion stays in pre-sliding; large strokes produce full loops whose width saturates at β₀. The dotted step is what a constant reversal correction assumes."

---

## 3. P2 — Masing check (is the loop predicted by the virgin curve?)

> **Superseded**, with §2 above.

**Current issues:** the prediction is drawn exactly on top of the loop, so the check is invisible; the residual panel is a flat line that says nothing without data.

**Layout:** three steps from left to right, then the residual.

- **(a) Virgin curve f** from V1, small: lag change against travel since reversal (linear x, 0 to s). Label "measured after a reversal".
- **(b) Stretched ×2:** the same curve with both axes doubled, drawn as a thin dark dashed line, labelled "2·f(Δu/2)".
- **(c) Placed on the loop:** the medium loop as a **thick, light teal** line (so it's visible underneath); the prediction as a **thin dark dashed** line starting at each turning point (one for the down branch, one mirrored for the up branch). Arrows from (b) to the two turning points.
- **Ghost alternative in (c):** the prediction with ×1.5 instead of ×2, dotted grey, labelled "if friction did not follow Masing". It visibly misses the loop.
- **(d) Residual:** branch − 2f against position, with a ±d_min band. Illustrative: the "pass" case (inside the band, solid) and the ×1.5 case (dotted grey, outside the band), so the reader knows what failure would look like.

**Caption idea:** "Each loop branch should equal the virgin curve stretched by two from its turning point. If it does, the loops follow from the virgin curve alone, as a GMS-type model predicts."

For the paper, P2 can also be a small inset in P1(a2).

---

## 4. M1 — Nested loop (does the stage remember?)

> **Superseded.** M1 was replaced by a plain-language animation,
> `rendered/companion/f1_M_memory.gif`, built by `render_f1_M_gif.py` to the
> brief in [`F1_M_Q_gifs_brief.md`](F1_M_Q_gifs_brief.md). The PNG was
> deleted. The section below is kept as the record of what the technical
> figure showed; the return-point-memory invariant it rested on is still
> asserted by `check_f1_companions.py` and now also by the GIF's renderer.

**Current issues:** the inner loop is a thin spike that's hard to see; the order of motion is unclear; a large text block carries the message.

**Layout:** main panel with a command icon on top.

**Command icon:** u against time: outer up ① → ②, down to the centre ③, inner up ④, inner down back to ⑤ (= ③), rest of the outer down to ⑥.

**Main panel:**
- **x:** commanded position relative to centre (µm), **y:** lag e (µm).
- **Outer loop** ochre, with arrows; **inner loop** dark slate, thicker, with its own arrows.
- **Illustration scale:** draw the inner loop larger than in the real sequence (e.g. 0.4 × outer instead of 0.1 ×) so it's clearly visible. State "illustration not to scale" in the caption.
- Numbered points ① to ⑥ on the path, matching the icon. ③ and ⑤ coincide: draw one open circle with the label "branch-off = rejoin".
- **The uninterrupted outer branch:** a light grey band under the outer down-branch.
- **Ghost alternative:** after ⑤, a dotted grey path that does *not* rejoin the outer branch but approaches it gradually from a different curve, labelled "without return-point memory (single-state model)".
- **Zoom inset** around ③/⑤: the return-point error as a small double arrow, "≈ 0 with memory".
- Remove the explanatory text block; it goes in the caption.

**Caption idea:** "A small inner loop is inserted halfway down the outer loop. With return-point memory (GMS) the stage resumes the outer branch exactly where it left it (③ = ⑤); a single-state model would follow a different path (dotted)."

---

## 5. Q1 — Rate check (does speed matter?)

> **Superseded.** Q1 was replaced by a plain-language animation,
> `rendered/companion/f1_Q_speed.gif`, built by `render_f1_Q_gif.py` to the
> brief in [`F1_M_Q_gifs_brief.md`](F1_M_Q_gifs_brief.md). The PNG was
> deleted. The section below is kept as the record of what the technical
> figure showed. Its nested line weights survive in the GIF; the ±d_min
> band and the drift-check ring do not, and the drift check is left out of
> the animation on the brief's own instruction.

**Current issues:** the three speeds overlap exactly, so only one loop is visible; the width panel is zoomed to a tiny range with identical points, which conveys nothing.

**Layout:** two panels.

**(a) Loops at three speeds:**
- Draw the three loops with **decreasing line weight** so all stay visible when they overlap: 12 s thick and light red, 4 s medium, 1.3 s thin and dark red. The 4 s repeat as a thin black line.
- **Ghost alternative:** a dotted grey loop, slightly wider, labelled "if friction depended on speed (slow loop wider)".
- Command icon: three short triangle-wave segments of increasing period, colour-coded by line weight.

**(b) Loop width against half-cycle time:**
- **x:** half-cycle time (log, 1–15 s); **y:** loop width with a range that shows meaningful change, e.g. ±15 % around the expected width, not ±2 %.
- Expected: three points on a flat line within the ±d_min band.
- **Ghost alternative:** a dotted line rising by, say, 5 % per decade of time, labelled "rate-dependent".
- The repeated 4 s point as an open marker, labelled "repeat (drift check)".

**Caption idea:** "If pre-sliding friction doesn't depend on speed, the loops at three speeds coincide and the widths lie on a flat line; a speed effect would appear as a slope (dotted)."

---

## 6. Smooth illustrative model

Same quasi-static model as the S-ramp GIF, with more elements:

- **N = 30** elements, slip limits aᵢ log-spaced from 0.01 to 20 µm.
- Weights γᵢ from a log-normal distribution in aᵢ, centred at ≈ s₉₀/3 (σ ≈ 1 in ln), then scaled so that
  - Σ γᵢ aᵢ = β₀/2 = 0.88 µm (saturation), and
  - Σ γᵢ ≈ 20 (initial slope of lost motion ≈ 0.95, i.e. nearly stuck).
- Lost motion for V1: ℓ(s) = s − Δx(s), with Δx from the incremental solve.
- Loops for P1/P2/Q1: symmetric triangle waves through the same solver; the loop width is read at the centre position.
- Nested loop for M1: the solver gives return-point memory automatically. For the "no memory" ghost, use a single-element first-order model, e.g. e relaxing exponentially towards ±β₀/2 with a length constant of s₉₀/3.

Replace placeholder s₉₀ = 3 µm and d_min = 5 nm with F0b values when available.

---

## 7. Optional: one "F1 at a glance" figure for the presentation

Four rows, two columns:

| Row | Left: what the stage does (command icon) | Right: what we measure (expected, with ghost) |
|---|---|---|
| V | wipe → reversal ramp | lost motion vs travel: stuck → pre-sliding → sliding |
| P | triangle waves, three sizes | loop width vs stroke, with the constant-correction step |
| M | outer loop with inner excursion | loop with return point, plus the "no memory" ghost |
| Q | same loop at three speeds | width vs speed: flat, plus the "rate-dependent" ghost |

Each row ends with a one-line tag on the right edge: "→ model elements", "→ Masing / width", "→ memory", "→ rate independence". Use it as the slide that introduces F1 before going into the individual plots.

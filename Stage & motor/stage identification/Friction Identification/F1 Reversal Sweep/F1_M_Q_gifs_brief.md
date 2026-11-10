# F1 · M and Q segments: simple GIFs

Same style as the P GIF: three panels, short everyday words, no technical terms, no equations, numbers only where they help. Same colours: going up teal (#1a6478), coming back orange (#c48a33), old curves light grey (#d0d4d8), dots dark (#2b3640), "what if" curves grey dotted (#9aa0a6). Serif font, ≥ 12 px, white background, thick lines, large arrows. Footer: "illustration, not measured data".

---

# GIF M — "Does the stage remember where it was?"

**Idea in one sentence:** the stage makes a big back-and-forth, and halfway back it takes a small detour; afterwards it picks up exactly where it left off, as if the detour never happened.

## Layout (≈ 1000 × 560 px)

```
┌───────────────────────────┬───────────────────────────┐
│ (A) What the motor does   │ (B) What the stage does   │
│  big move + small detour  │  big loop, small detour   │
├───────────────────────────┴───────────────────────────┤
│ (C) Zoom: where the detour leaves and rejoins          │
└───────────────────────────────────────────────────────┘
 one-line message (changes per phase)
```

### (A) What the motor does
- Motor position over time, with a moving dot and a direction arrow.
- The path: **up** (big), **back down halfway**, **small up, small down** (the detour), **rest of the way down**.
- Mark the detour with a light box and the label **"small detour"**.

### (B) What the stage does
- **x:** "where the motor is"; **y:** "how far the stage is behind" (labels "going up" / "going down" at the bottom / top, as in the P GIF).
- The big loop drawn live: up in teal, down in orange.
- The **detour** drawn live in dark grey, thicker.
- Before the detour starts, show the rest of the big down-path as a **faint grey line** labelled **"where it was heading"**.
- Mark the point where the detour leaves the big path with a small open circle, **"leaves here"**.
- After the detour, the stage returns to the same circle and continues on the faint line. Change the label to **"back on track"**.
- **Illustration scale:** draw the detour about one third the size of the big move (larger than in the real test), so it is clearly visible.

### (C) Zoom
- A magnified view around the circle, appearing when the detour starts.
- Shows the detour leaving and returning to the same point.
- At the end, a dotted grey **"forgetful stage"** path: after the detour it does *not* return to the circle but drifts back onto the big path only gradually. Label: **"a stage that forgets would take this path"**.

## Sequence

| Phase | (A) | (B) | (C) | Message |
|---|---|---|---|---|
| 0 | empty | empty | hidden | "The stage makes a big move back and forth, with a small detour on the way back." |
| 1 · big move up | dot rises | teal path | hidden | "Going up: the stage follows, a little behind." |
| 2 · halfway down | dot falls to halfway | orange path; faint "where it was heading" appears | hidden | "Coming back down…" |
| 3 · detour (slow) | small up, small down | dark small loop | zoom appears, detour drawn | "A small detour: up a little, then back." |
| 4 · continue down | rest of the way down | orange continues on the faint line; label "back on track" | circle highlighted | "After the detour, it continues exactly where it left off." |
| 5 · what-if | — | — | dotted "forgetful" path | "A stage that forgets would take a different path." |
| 6 · hold | all | all | all | "The stage remembers where it was." |

## Timing (20 fps, ≈ 22 s)

| Phase | Frames |
|---|---|
| 0 | 40 |
| 1 | 60 |
| 2 | 40 |
| 3 (slowed down) | 80 |
| 4 | 50 |
| 5 | 60 |
| 6 | 60 |

## Curves (renderer only)
- Use the same smooth friction model as the P GIF.
- Big move: 10 µm; detour: 3 µm (illustration scale), starting halfway down.
- Path: −5 → +5 → 0 → +3 → 0 → −5 µm (positions relative to the centre).
- "Forgetful" alternative (dotted, illustration only): after the detour, let the lag approach the big down-path gradually, e.g.
  `lag_forget(u) = lag_big(u) + (lag_at_detour_end − lag_big(u_end_detour)) · exp(−|u − u_end_detour| / 1.0 µm)`.

## One-sentence caption
> "Halfway back, the stage takes a small detour. Afterwards it picks up exactly where it left off, as if the detour never happened. A stage that forgot its past would follow a different path."

---

# GIF Q — "Does speed matter?"

**Idea in one sentence:** the same back-and-forth is done fast, medium and slow, and the stage draws the same loop every time.

## Layout (≈ 1000 × 560 px)

```
┌───────────────────────────┬───────────────────────────┐
│ (A) What the motor does   │ (B) What the stage does   │
│  same move, three speeds  │  loops drawn on top       │
├───────────────────────────┴───────────────────────────┤
│ (C) Gap for each speed (fills in)                      │
└───────────────────────────────────────────────────────┘
 one-line message (changes per phase)
```

### (A) What the motor does
- Motor position over time, three back-and-forths of the **same size**, one after the other: **fast**, **medium**, **slow**.
- Each one labelled above it: **"fast"**, **"medium"**, **"slow"**, with a small clock icon or "1 s", "4 s", "12 s" per pass.
- **The dot moves at different animation speeds**: quick for fast, slow for slow (see timing), so the viewer *feels* the difference.

### (B) What the stage does
- Same axes as the P GIF ("where the motor is", start → end; "how far the stage is behind").
- Each speed draws its loop with a different line style, on top of the previous ones:
  - fast: thick light line;
  - medium: medium line;
  - slow: thin dark line.
- Because the loops coincide, each new loop traces exactly over the previous one. Add a short tag when a loop finishes: **"same loop"**.

### (C) Gap for each speed
- **x:** "speed" with three labels: **fast · medium · slow**. **y:** "gap", same scale as in (B).
- After each loop, its point drops in. The three points lie on a **flat line**.
- At the end, a dotted grey **"if speed mattered"** line, sloping up towards "slow".

## Sequence

| Phase | (A) | (B) | (C) | Message |
|---|---|---|---|---|
| 0 | empty | empty | empty | "The same small move, done at three speeds." |
| 1 · fast | quick back-and-forth ×2 | thick light loop | first point | "Fast: this is the loop." |
| 2 · medium | medium back-and-forth ×2 | medium loop traces over it | second point, same height | "Medium speed: the same loop." |
| 3 · slow | slow back-and-forth ×2 | thin dark loop traces over it | third point, same height | "Slow: still the same loop." |
| 4 · what-if | — | — | dotted sloping line | "If speed mattered, the gap would change like this." |
| 5 · hold | all | all | all | "Speed doesn't change what the stage does." |

## Timing (20 fps, ≈ 24 s)

| Phase | Frames |
|---|---|
| 0 | 40 |
| 1 · fast (2 back-and-forths, 10 frames per half) | 40 |
| hold | 20 |
| 2 · medium (20 frames per half) | 80 |
| hold | 20 |
| 3 · slow (40 frames per half) | 160 |
| hold | 20 |
| 4 | 50 |
| 5 | 60 |

The animation speeds keep the real 1 : 3 : 9 ratio only loosely (1 : 2 : 4); exact ratios would make the slow pass too long. The labels "1 s / 4 s / 12 s" carry the real values.

## Curves (renderer only)
- Same smooth friction model; stroke 3 µm (around the "loop opens" size from the P GIF).
- The model doesn't depend on speed, so all three loops are identical. That's the expected outcome.
- "If speed mattered" (dotted, illustration only): gap rising by about 10 % from fast to slow.

## Optional: the drift check
The real test repeats the medium speed at the end, to check nothing changed during the test. If you want to show it, add a fourth short pass labelled **"medium again (check)"** whose loop also traces over the others. Otherwise leave it out to keep the GIF simple.

## One-sentence caption
> "The same small move is done fast, medium and slow. The stage draws the same loop each time, so its behaviour doesn't depend on speed."

---

## Using the three GIFs together
In a presentation, show them in this order: **P** (bigger moves, bigger gap, up to a limit) → **M** (the stage remembers where it was) → **Q** (speed doesn't matter). Together they make the three plain-language points behind the friction model: the gap depends on the move size, the stage remembers its past, and speed doesn't change it.

---

## Implementation record

Everything above is the brief as received. This section records what
`render_f1_M_gif.py` and `render_f1_Q_gif.py` actually do, and where they
depart. The pieces the three GIFs share — colours, canvas, palette, frame
capture, GIF writer — were moved into `f1_gif_common.py` when the second and
third were written; the P GIF now imports them instead of carrying its own
copy, and its output is unchanged apart from one extra palette ramp.

### M

| Brief | Built | Why |
|---|---|---|
| the exponential "forgetful" formula in *Curves* | the same solver with the outer excursion erased | The brief's formula needs `lag_at_detour_end` to differ from `lag_big`, and in a model with return-point memory those are the same number to machine precision, so the formula gives a zero offset and draws nothing. What is drawn instead is a stage whose whole history is the detour: it responds from the detour's turning point with no memory of the big move, comes back onto the big path **0.090 µm low** and is still 0.013 µm off five microns later. It is a worked alternative rather than a drawn curve, but it is still one way to lack memory, not a statement about every model. |
| one big back-and-forth | one **drawn**, one run before it | As in the P GIF. Without it the loop starts from the wiped state, the up-branch is a flat saturated line with no bend, and the loop shown is not the repeatable one. The footer says the big move is repeated until the loop settles. |
| (C) "a magnified view around the circle" | magnified ×5 vertically, nearly the full width | The forgetful alternative closes its 0.090 µm over about three microns, so the panel has to show that whole run-out; there is nothing to gain by cropping horizontally. The title states the magnification. |
| "leaves here" / "back on track" beside the circle | the same, on a white patch | The faint "where it was heading" line runs through that spot. |
| "where it was heading" label on the faint line | the label in the empty middle of the panel, on a leader | The faint line runs along the top of the panel, where the two branch labels and the circle already are. |
| ≈ 22 s | 19.5 s | The brief's own frame table sums to 390 frames, which is 19.5 s at 20 fps. The frames are as tabled. |

The detour is drawn at 3 µm against a 10 µm move, as the brief asks, which
is larger than the planned test; the footer says so. One consequence worth
knowing: halfway down a 10 µm move the stage is 5 µm past the reversal and
s₉₀ is 3 µm, so the branch it rejoins is already most of the way to
saturation and only gently curved. The claim is carried by panel (C), where
the two paths are told apart, and `main()` asserts it numerically: the
detour closes to **0 µm** and what follows the detour is the uninterrupted
branch to **4 × 10⁻¹⁶ µm**.

### Q

| Brief | Built | Why |
|---|---|---|
| (A) three passes drawn in sequence | drawn against **animation** time, so the passes are 1 : 2 : 4 wide | Against real time the slow pass would be twelve times the width of the fast one and its dot would cross the panel three times *faster*, which is the opposite of what the panel is for. Against animation time the dot advances steadily and the zigzag is tight for fast and stretched for slow, which is what a viewer reads as speed. The labels carry the real 1 s / 4 s / 12 s and the footer states the 1 : 2 : 4 drawing. |
| (C) "same scale as in (B)" | its own zero-based scale | (B) is a signed lag of ±1.05 µm; (C) is an unsigned gap. Starting (C) at zero is what keeps the 10% alternative looking like 10% — on a scale cropped to the data it would look like a large effect. |
| a tag "same loop" when a loop finishes | the same, from the second pass on | After the first pass there is nothing for it to be the same as. |
| the optional drift-check pass | left out | The brief says to leave it out to keep the GIF simple. It is in the real Q block and in the procedure. |
| — | a "gap" arrow in (B), drawn once after the first pass and kept | (C)'s y axis is called "gap" and nothing else in this GIF says what that means. The P GIF introduces it, but Q has to stand on its own. |

The three loops are not merely close. Every pass samples the stroke at a
divisor of the slow pass's spacing, so `main()` compares them at coincident
commanded positions with no interpolation: the loops agree to
**3 × 10⁻¹⁵ µm** and the three gaps are the same floating-point number,
0.904447 µm, which is also `width_at_centre(3 µm)`. The only difference left
is that the fast pass is drawn with ten points per half-cycle, so it cuts
the corner at a reversal by 6.8 nm, about **0.8 px**.

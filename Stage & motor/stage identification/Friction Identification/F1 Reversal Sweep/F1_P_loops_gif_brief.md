# F1 · P segment: simple GIF ("back and forth, and what the stage does")

**Goal:** show anyone, without technical background, what happens when the stage is moved back and forth by different amounts: the motor is told to go back and forth, and the stage follows a little behind. Drawing *how far behind* against *where the motor is* gives a loop. Small moves give a thin loop, large moves a wide one, and the width stops growing at some point.

**Rule for all text:** short, everyday words. No technical terms (no "hysteresis", "pre-sliding", "lag e", "Masing", "β₀", "s₉₀"). Numbers only where they help (stroke size in µm).

---

## 1. Layout (≈ 1000 × 560 px)

```
┌───────────────────────────┬───────────────────────────┐
│ (A) What the motor does   │ (B) What the stage does   │
│  position ↕ over time     │  loop drawn live          │
│  moving dot               │  moving dot               │
├───────────────────────────┴───────────────────────────┤
│ (C) Loop size for each move size (fills in at the end)│
└───────────────────────────────────────────────────────┘
 one-line message under the panels (changes per phase)
```

### (A) What the motor does
- **x:** time (no numbers). **y:** "motor position" (no numbers, or "small / large" ticks).
- A grey line showing the planned back-and-forth for the current move, and a **dot** moving along it.
- Small arrow next to the dot showing the direction of motion (↑ or ↓).
- Title: **"What the motor is told to do"**.

### (B) What the stage does
- **x:** "where the motor is", from **start** to **end** of the move (labels at the two ends only).
- **y:** "how far the stage is behind", fixed scale for all moves, with two labels: **"behind (going up)"** at the bottom and **"behind (going down)"** at the top. No numbers, or just one tick "≈ 1.8 µm" for the full gap at the end.
- The loop is **drawn live** as the dot moves:
  - going up: one colour (e.g. teal), **"going up"** label with an arrow;
  - coming back: second colour (e.g. orange), **"coming back"** label with an arrow.
- A **double arrow** in the middle appears after the first full back-and-forth, labelled **"gap"**.
- Title: **"What the stage does"**.

### (C) Loop size for each move size
- **x:** "size of the move" (small → large; three labelled points: "1 µm", "3 µm", "100 µm").
- **y:** "gap", same scale as the gap in (B).
- After each move finishes, its point drops in. At the end, a smooth line joins them, rising and then flattening.
- A thin dotted horizontal line at the top, labelled **"maximum gap"**.
- Title: **"Bigger moves, bigger gap, up to a limit"**.

---

## 2. Sequence

Three moves, each going back and forth **twice**, then the summary.

| Phase | (A) | (B) | Message line |
|---|---|---|---|
| 0 | empty | empty | "The motor moves the stage back and forth. Watch how far the stage falls behind." |
| 1 · small move (≈ 1 µm) | small back-and-forth | very thin loop, almost a line; gap tiny | "Small move: the stage hardly moves, so the loop stays thin." |
| 2 · medium move (≈ 3 µm) | medium back-and-forth | loop opens up | "Medium move: the stage starts to follow, so the loop opens." |
| 3 · large move (≈ 100 µm) | large back-and-forth | wide loop with flat top and bottom; full gap | "Large move: the stage follows, but always a bit behind. The gap is now as big as it gets." |
| 4 · summary | all three moves shown faintly | the three loops faintly overlaid | "The gap grows with the size of the move, then stops growing." |

- **Between phases:** the previous loop fades to light grey and stays, so the reader can compare.
- **Panel (B) keeps the same width for every move** (x goes from "start" to "end" of the move), so only the *height* of the loop changes. In (A), the y-axis rescales per move. Show a small tag in (A): "×3 bigger", "×30 bigger", so the reader knows the moves differ in size.

---

## 3. Timing (20 fps, ≈ 25 s)

| Phase | Frames |
|---|---|
| 0 intro | 40 (2 s) |
| each move: 2 back-and-forths, 30 frames per half (up or down) | 120 (6 s) per move |
| hold after each move (gap arrow + point drops in (C)) | 20 (1 s) per move |
| 4 summary | 60 (3 s) |
| **Total** | ≈ 520 frames, ≈ 26 s |

Make it loop with a 1 s pause at the end.

---

## 4. Look

- **Background white; text dark grey; serif font, as in the other figures; ≥ 12 px.**
- **Colours:** going up teal (#1a6478), coming back orange (#c48a33), old loops light grey (#d0d4d8), dots dark (#2b3640).
- **Lines thick** (2.5–3 px); arrows large; no grid or a very faint one.
- **No equations, no axis numbers** except the three move sizes in (C) and optionally the one "≈ 1.8 µm" tick.
- Footer, small and light grey: "illustration, not measured data".

---

## 5. How to generate the curves (for the renderer only, not shown in the GIF)

Use the smooth friction model from the F1 companion plots (≈ 30 elements; maximum gap = 1.76 µm; most of the change happens between 1 and 5 µm).

```python
import numpy as np

# smooth friction model (placeholder values)
a = np.logspace(-2, np.log10(20), 30)                 # element limits (µm)
w = np.exp(-0.5 * (np.log(a / 1.0))**2)               # spread centred near 1 µm
g = w * (0.88 / np.sum(w * a))                        # max lag 0.88 µm per side
# check: g.sum() should be large (≥ ~10) so small moves barely move the stage;
# if not, shift the centre (1.0) of the spread to smaller values.

def run(path_u, z):
    """Move the motor along path_u (µm); return stage lag at each sample."""
    u_prev, x, lags = path_u[0], path_u[0], []
    for u in path_u[1:]:
        du = u - u_prev
        stick = np.abs(z + np.sign(du) * 1e-12) < a
        dx = du / (1 + np.sum(g[stick]))
        z = np.clip(z + dx, -a, a)
        x += dx; u_prev = u
        lags.append(x - u)
    return np.array(lags), z

def back_and_forth(s, cycles=2, n_half=30):
    up = np.linspace(-s/2, s/2, n_half)
    return np.concatenate([np.r_[up, up[::-1]] for _ in range(cycles)])

z = a.copy()                       # start: stage pushed up earlier (as after the wipe)
for s in [1.0, 3.0, 100.0]:
    u = back_and_forth(s)
    lag, z = run(np.r_[u[0], u], z)
    # panel B: x = (u + s/2) / s  (0 = start, 1 = end), y = lag (flip sign so "behind" reads naturally)
    # panel C: gap = lag on the way down minus lag on the way up, at the middle (x = 0.5)
```

Panel (B) plots the lag against the position within the move (0 = start, 1 = end). Panel (C) takes the gap at the middle of the move, from the second back-and-forth.

---

## 6. One-sentence caption (if used on a slide)

> "The motor moves the stage back and forth by small, medium and large amounts. The stage always trails a little; plotting how far it trails against where the motor is gives a loop. The loop's gap grows with the size of the move, then stops growing."

---

## 7. Implementation record

Everything above is the brief as received. This section records what
`render_f1_P_gif.py` actually does, and where it departs.

| Brief | Built | Why |
|---|---|---|
| §5 placeholder solver | the calibrated 30-element solver imported from `render_f1_companions.py` | The brief calls its own snippet a placeholder and asks for "the smooth friction model from the F1 companion plots". That model is already in the repository, already satisfies equilibrium, saturation β₀ = 1.757 µm and the illustrative s₉₀ = 3 µm, and is already checked by `check_f1_companions.py`. Importing it means the GIF and the three surviving figures cannot drift apart. |
| "×30 bigger" tag | "×100 bigger" | The moves are 1, 3 and 100 µm, as §2 and §5 specify, so the third tag is ×100. The tags are computed from the strokes rather than written out. |
| two back-and-forths per move | two **drawn**, one run before them | The first back-and-forth after a change of stroke inherits the previous state: for the 1 µm move its gap is 0.475 µm against the settled 0.037 µm, thirteen times wider. Drawn as-is, the opening frames of each move would contradict the message line, and the "gap" arrow would collapse halfway through the move. The bench handles this the same way, with three cycles per stroke and the later ones read. The footer says "each move is repeated until the loop settles" and `main()` asserts the two drawn cycles agree to 1e-9. |
| "going up" / "coming back" labels with arrows, inside the panel | a two-colour key above the panel, plus an arrowhead riding along each branch | The loops fill the panel, and no corner stays clear across all three moves: a floating label collided with the large loop's flat branches. The arrowheads sit on the lines, so they cannot collide, and they move with the trace. |
| y labels "behind (going up)" / "behind (going down)" | the same two labels, as coloured y ticks at the two saturation levels | Rotated axis text in the margin was cramped and read as one block. As tick labels they are horizontal, colour-matched to the branches, and placed exactly where the stage is furthest behind in each direction. |
| ≈ 520 frames, ≈ 26 s | 573 frames, 28.6 s | Ten frames per move are spent travelling into position before the first back-and-forth, which the brief's budget does not count, and the 1 s end pause is included. Pillow merges identical consecutive frames when writing, so the file holds 398 frames with their durations summed; the playback is unchanged. |
| one "≈ 1.8 µm" tick in (B) | that tick in (C) only | The number belongs next to the "maximum gap" line it labels. (B) stays free of numbers, as §4 asks. |

The gap arrow's head size is scaled to the gap it spans: at 1 µm the gap is
about 6 px, and fixed heads alone would have made it look several times
wider than it is.

# F1 — Reversal-amplitude sweep (pre-sliding identification)

> **Question:** How does the lost motion grow with the distance travelled after a reversal, does the stage remember nested reversals, and is this behaviour the same everywhere along the travel?

F1 is the central friction experiment. F0b measured the *scale* of the friction effects at one position. F1 measures their *shape* in detail, at three positions, with enough structure to identify the pre-sliding part of the friction model (the GMS elements) and to test that friction is separate from geometry.

**Values in ⟨angle brackets⟩ depend on F0b** (d_min = minimum detectable step, s₉₀ = pre-sliding length, ℓ_sat = saturation lost motion). Defaults are given for s₉₀ ≈ 3 µm and d_min ≈ 5 nm; regenerate the table once F0b is processed.

---

## 1. What F1 measures

| Sub-block | Measures | Identifies |
|---|---|---|
| **V** virgin curves | lag against travel after a full reversal, in fine resolution, both directions | the virgin (initial loading) curve → GMS element stiffnesses and slip limits; deadband; saturation |
| **P** symmetric loops | closed hysteresis loops of increasing stroke, 3 cycles each | loop width against stroke; Masing check (branches = 2 × virgin curve) |
| **M** nested loops | a small inner loop inside a larger loop, then continuation | return-point memory: does the stage rejoin the outer branch where it left it? (GMS yes, LuGre no) |
| **Q** rate check (middle position only) | the same loop at three speeds | whether pre-sliding is rate-independent, as GMS assumes |

**Three positions** test separation from geometry: the same friction curves at different positions and screw angles, after removing the geometric model, mean friction doesn't depend on position.

---

## 2. Configuration

| Item | Setting |
|---|---|
| Driver | TMC2209, **MRES 1/256 for the whole session**, SpreadCycle, interpolation off, IRUN = IHOLD (as F0b) |
| Step rates | ramps from ~6 steps/s up to the transit rate; transits 1 mm/s (25.6 kHz), or the F0b fallback rate if lower |
| **Test positions** (all on the 1024-µstep grid; different screw angles θ) | **x_A = 25.00 mm** (θ ≈ 180°, motor side) · **x_B = 60.00 mm** (θ ≈ 0°, same as F0b x₀) · **x_C = 100.48 mm** (θ ≈ 86°, far side) |
| Why these | inside the well-constrained span (5–110 mm); away from the 72–80 mm screw feature and the 95 mm reversal bump; x_B links to F0b; x_A and x_C differ in wobble amplitude and screw angle |
| Memory wipe | 3 mm approach in the stated direction at 1 mm/s, 2 s settle |
| Reference pair | 15 / 105 mm, both directions, start and end (R1, R3); R2 between positions if time allows |
| Recording | IDS continuous at 1 kHz, scope record ≥ 60 min; BME280 and thermistors if fitted |
| Command table | host-generated, board replays step counts; row, segment type, target, rate, dwell, tag (position, sub-block, stroke, cycle, half) |

---

## 3. Sub-block definitions (per position x_p)

### V — Virgin curves (≈ 3 min per position)
1. Wipe from **below** to x_p.
2. Reverse: ramp **down** 100 µm in two speeds: 0.25 µm/s for the first ⟨2·s₉₀ ≈ 6 µm⟩, then 2 µm/s to 100 µm (≈ 70 s total).
3. Wipe from **above** to x_p, then ramp **up** in the same way.

### P — Symmetric loops (≈ 6 min per position)
Loops are **centred on x_p**: x_p − s/2 ↔ x_p + s/2.
1. Wipe from below to x_p.
2. For each stroke s, **increasing**: move to x_p − s/2, then 3 full cycles between x_p − s/2 and x_p + s/2, return to x_p.
3. **At x_B only:** repeat the list **decreasing** (history check, as in F0b L).

**Stroke grid** (dense around s₉₀, anchors above):

| | Values |
|---|---|
| Rule | s = ⟨s₉₀⟩ × {0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5}, then 100 µm, 1 mm, 3 mm; every stroke ≥ max(1 µstep, 5·d_min); rounded to whole µsteps |
| Default (s₉₀ = 3 µm) | 0.15, 0.3, 0.6, 0.9, 1.5, 2.1, 3, 4.5, 6, 9, 15, 100, 1000, 3000 µm (14 strokes) |

Speed: each half-cycle ≈ 4 s (rate = s / 4 s, capped at 0.75 mm/s), with 0.5 s pauses at turning points. Strokes ≤ 5 µsteps as `step` segments with 0.5 s dwell per step.

Centred loops cover the same commanded positions in both directions, so the loop width at a given position is free of geometric and microstep error. They also give the standard loop shape for the Masing check.

### M — Nested loops (≈ 1.5 min per position)
Outer stroke S_o = ⟨3·s₉₀ ≈ 10 µm⟩, inner stroke S_i = ⟨0.3·s₉₀ ≈ 1 µm⟩.
1. Wipe from below to x_p − S_o/2, ramp up to x_p + S_o/2 (outer up-branch).
2. Ramp down to x_p (halfway down the outer down-branch).
3. **Inner loop:** up S_i, down S_i (back to x_p).
4. Continue down to x_p − S_o/2 (rest of the outer down-branch).
5. Repeat once.

Return-point memory: after the inner loop, the trace should rejoin the outer down-branch at x_p and continue on it as if the inner loop had not happened.

### Q — Rate check (x_B only, ≈ 2.5 min)
Centred loops of stroke ⟨s₉₀⟩, 3 cycles each, at half-cycle times of **1.3 s, 4 s and 12 s** (×3, ×1, ×⅓ speed), in the order 4 – 1.3 – 12 – 4 s (repeat of the middle speed at the end as a drift check).

---

## 4. Session

| Block | Content | Time (est.) |
|---|---|---|
| R1 | reference pair | 3.5 min |
| **A** (x_A = 25 mm) | V · P (increasing) · M | ≈ 11 min |
| **B** (x_B = 60 mm) | V · P (increasing + decreasing) · M · Q | ≈ 17 min |
| R2 (optional) | reference pair | 3.5 min |
| **C** (x_C = 100.48 mm) | V · P (increasing) · M | ≈ 11 min |
| R3 | reference pair | 3.5 min |
| — | transits between positions | ≈ 3 min |
| **Total** | | **≈ 49–53 min** (drop R2 if > 52 min) |

Order of positions: A → B → C, so B (the link to F0b and the richest block) sits in the middle of the session, closest to the reference anchors on both sides in time.

---

## 5. Pre-run checks
- [ ] Stroke grid regenerated from F0b (s₉₀, d_min); all strokes whole µsteps and ≥ 5·d_min
- [ ] x_A, x_B, x_C on the 1024-µstep grid; MSCNT logged at every block start
- [ ] Centred loops: start/end points x_p ± s/2 rounded symmetrically (same number of µsteps each side)
- [ ] Total duration ≤ 52 min; scope record ≥ 60 min
- [ ] Dry run of the first ~40 rows (echo of `u_cmd`)

---

## 6. Analysis

### 6.1 Pre-processing
- Drift and session gain from R1–R3.
- **Geometry:** subtract Ĝ along the path for strokes ≥ 100 µm and for the 100 µm virgin ramps; for smaller strokes the change is ≤ 0.1 µm, but subtract anyway for consistency.
- **Lag:** e(t) = x_meas − u_cmd (after drift and Ĝ), relative to a datum window at the start of each sub-block.
- **Loop width at fixed commanded position:** w(u) = e↓(u) − e↑(u). Microstep and geometric error cancel exactly, because both branches pass the same u.

### 6.2 Virgin curves (V)
- e against travel since reversal, 0–100 µm, per position and direction.
- Read the deadband, initial slope, s₉₀ and saturation level ℓ_sat/2 per side.
- Compare up vs down (symmetry) and across positions (separation).

### 6.3 Loops (P)
- Loop width against stroke: w_max(s), the width at the turning point.
- **Masing check:** scale the virgin curve by 2 and overlay it on the loop branches. GMS-type friction satisfies this; a systematic mismatch means the model needs more than rate-independent elements.
- Cycle 1 vs 2 vs 3: does the loop close after the first reversal?
- x_B increasing vs decreasing (history effect, as in F0b).

### 6.4 Nested loops (M)
- After the inner loop, does the trace rejoin the outer down-branch at x_p within d_min? Report the return-point error.

### 6.5 Rate check (Q)
- Loop width and shape at the three speeds; the repeated 4 s loop gives the drift bound.

### 6.6 Identification
- Fit N = 3–5 GMS elements (relative stiffness γᵢ, slip limit aᵢ) to the virgin curves of **x_A and x_B**, jointly with their loop widths.
- **Hold-out:** predict x_C's virgin curve, loops and nested loop with those parameters. A prediction within a few d_min confirms position independence; a systematic mismatch means a position term is needed (or the route-B iteration).
- Cross-check: fitted saturation Σγᵢaᵢ ≈ β₀/2 from E2–E5.

### 6.7 Separation checks (summary)

| Check | Pass means |
|---|---|
| Virgin curves and loop widths at x_A, x_B, x_C overlap (after Ĝ) | friction independent of position and screw angle |
| Up/down virgin curves mirror each other | no geometric leak |
| Saturation ≈ β₀/2 | the friction layer explains the reversal found in E0–E5 |
| Hold-out prediction of x_C within a few d_min | the identified model transfers along the travel |

---

## 7. Outputs
- `f1_command_table.csv` (+ hash), raw recording, firmware log
- `f1_virgin.csv` (position, direction, travel, e), `f1_loops.csv` (position, stroke, cycle, w_max, branch data), `f1_nested.csv`, `f1_rate.csv`
- `f1_gms_params.json`: γᵢ, aᵢ, fit residuals, hold-out residuals
- Figures: virgin curves (3 positions overlaid), loop width against stroke, Masing overlay, nested loop, rate check

---

## 8. Render instructions: F1 overview figure

Same visual grammar as the final F0b figure: ribbon + commanded trajectory (a), amplitude coverage (b), plus a small travel map (c) because F1 adds positions.

### Layout (full text width, ≈ 7.2 × 6.5 in)

```
(a) ribbon:  R1 | A: V P M | B: V P M Q | (R2) | C: V P M | R3
    trajectory: commanded position relative to the CURRENT test position, symlog
(b) amplitude coverage: rows V, P, M, Q
(c) travel map: 0–120 mm bar with x_A, x_B, x_C, reference positions, wobble profile
```

### (a) Ribbon + trajectory
- **Ribbon, two tiers:** top tier the position groups (R1 · **A 25 mm** · **B 60 mm** · R2 · **C 100.48 mm** · R3), bottom tier the sub-blocks inside each group (V, P, M, Q). Widths proportional to time.
- **y-axis:** commanded position relative to the *current* test position x_p, on a symmetric log scale with ±µm/mm tick labels and a shaded linear band (|Δ| < 1 µm), as in F0b. Using x_p per group keeps all three positions on the same axis.
- **Transits between positions** (75 mm moves) and reference pairs: grey dashed or purple dashed, as in F0b; draw the inter-position transits as short grey breaks with a label ("→ 60 mm").
- **Colours by sub-block:** V virgin dark slate (#2b3640), P loops teal (#1a6478), M nested ochre (#c48a33), Q rate red (#b5473a); positioning and wipes grey; reference pairs purple dashed.
- Expect P to appear as a filled envelope (as in the F0b ladder); draw it with a thin outline and lighter fill so it reads as oscillation.

### (b) Amplitude coverage
- x-axis: commanded amplitude (log, ~0.03 µm … 5 mm).
- Rows:
  - **V virgin:** a bar from 0.04 µm to 100 µm, with a darker segment over the slow part (0 … 2·s₉₀).
  - **P loops:** markers at the 14 strokes; three marker shapes or small offsets for x_A, x_B, x_C (all identical in amplitude, so a single row with "×3 positions" is enough).
  - **M nested:** two markers, inner and outer, joined by a bracket ("inner inside outer").
  - **Q rate:** one marker at s₉₀ with three small tick labels "1.3 / 4 / 12 s".
- **Guides (labels at the top, offset from the lines):** d_min (from F0b), 1 µstep = 39 nm, **s₉₀ band from F0b** (shaded, the region the grid is dense around), β₀ = 1.76 µm, electrical cycle 40 µm.

### (c) Travel map (small strip)
- Horizontal bar 0–120 mm (the stage travel), with:
  - x_A, x_B, x_C as filled markers with labels and screw angle (θ ≈ 180°, 0°, 86°);
  - reference positions 15 / 105 mm as open purple markers;
  - weakly constrained zones (< 5, > 110 mm) shaded grey;
  - known features as small ticks: 72–80 mm screw feature, 95 mm reversal bump;
  - optionally, the first-harmonic amplitude A₁(x) from the geometric model as a faint line above the bar, so the reader sees that the three positions sample different wobble amplitudes.

### Style
- Serif, panel letters only, no title; findings in the caption.
- Black-and-white safe: sub-blocks identified by ribbon letters and row labels, not colour alone; reference pairs dashed.
- Text ≥ 7 pt at final size.

### Caption (suggested)
> F1 command sequence. (a) Commanded position relative to the current test position on a symmetric logarithmic scale; the ribbon marks the three test positions and, within each, the virgin curves (V), symmetric loops (P), nested loops (M) and, at the middle position, the rate check (Q); reference pairs and transits are dashed/grey. (b) Amplitudes covered by each sub-block; the loop strokes are dense around the pre-sliding length s₉₀ measured in F0b and extend to 3 mm, bracketing the expected saturation lost motion β₀. (c) Test positions along the travel at different screw angles and wobble amplitudes, inside the well-constrained span and away from known local features. ⟨n⟩ segments, ⟨t⟩ min.

# F2 — Dwell (frictional ageing and creep)

> **Question:** Does the time the stage rests before a move change how that move starts, and if so, by how much and how does it grow with the rest time?

Static friction can grow while contacts are at rest: the contact spots creep and interlock (ageing). After a long rest, the next move may lose more motion after a reversal, or stick and then jump when continuing in the same direction. In operation, a stage often rests for unpredictable times between moves, so if this effect is larger than the measurement floor, the friction model needs a dwell term (or operation needs a fixed dwell). F2 measures the size of the effect and how it depends on rest time.

**F2 is conditional.** Run the full version only if F0b's dwell probe (1 s vs 60 s) showed a difference larger than ⟨3·d_min⟩. Otherwise run the reduced version (§4.3), which is enough to report an upper bound.

Values in ⟨angle brackets⟩ depend on F0b (s₉₀, d_min, dwell-probe result). Defaults: s₉₀ = 3 µm, d_min = 5 nm.

---

## 1. What F2 measures

| Observable | When | Shows |
|---|---|---|
| **Creep during the dwell** | while the stage rests | how far the carriage keeps drifting after stopping (pre-sliding relaxation, E5 hinted at ≈ 30 nm in 0.3 s) |
| **Lost motion after a reversal following the dwell** (mode R) | the post-dwell move reverses the last motion | whether ageing enlarges the deadband, the pre-sliding length or causes an overshoot ("stiction peak") before settling at β₀ |
| **Breakaway when continuing after the dwell** (mode C) | the post-dwell move continues in the same direction | whether the carriage sticks at restart and then jumps (extra lag, then a step) |

Each is measured against the **dwell time T_d**, on a logarithmic grid, because ageing typically grows with log(T_d).

---

## 2. Configuration

| Item | Setting |
|---|---|
| Driver | TMC2209, MRES 1/256 for the whole session, SpreadCycle, interpolation off, **IRUN = IHOLD** (the hold current during the dwell is the same as during motion; record it) |
| Positions | **x_B = 60.00 mm** (primary, full grid); **x_A = 25.00 mm** (check, reduced grid). Both on the 1024-µstep grid, as in F1 |
| Memory wipe | 3 mm approach **from below** at 1 mm/s for every trial, so the state before the dwell is always the same (last motion upward) |
| Dwell | the board holds position with no steps for exactly T_d (logged start and end) |
| Post-dwell ramp | ⟨2·s₉₀ ≈ 6 µm⟩ at **0.5 µm/s** (12 s), down (mode R) or up (mode C); then 2 s hold |
| Reference pair | 15 / 105 mm, both directions, at start and end (R1, R3) |
| Recording | IDS continuous at 1 kHz, scope record ≥ 60 min; BME280 and thermistors if fitted (temperature during long dwells matters) |
| Randomisation | within each repetition, the (T_d, mode) list is shuffled with a logged seed, so dwell time is not correlated with session time |

---

## 3. One trial

```
wipe from below (3 mm, 1 mm/s) → settle 0.5 s → DWELL T_d → ramp 6 µm at 0.5 µm/s → hold 2 s
                                                             ↓ mode R: ramp DOWN (reversal)
                                                             ↓ mode C: ramp UP   (continuation)
```

- The wipe is identical in both modes, so the only differences are T_d and the direction of the post-dwell move.
- The 0.5 s settle after the wipe is part of every trial and defines t = 0 of the dwell.
- **T_d = 0.1 s** acts as the "no dwell" baseline.
- Reference points for comparison: F1's virgin ramps followed a 2 s settle (≈ T_d = 2 s); E5's ISO stops were 2.5 s.

---

## 4. Session

### 4.1 Grids

| Block | Position | Modes | Dwell times T_d (s) | Repetitions |
|---|---|---|---|---|
| **B** | 60 mm | R and C | 0.1, 0.5, 2, 10, 30, 120 | 3 (120 s: 2) |
| **A** | 25 mm | R only | 0.5, 10, 120 | 2 |
| **L** (long dwell) | 60 mm | R | 300 | 1 |

### 4.2 Timing

| Block | Content | Time (est.) |
|---|---|---|
| R1 | reference pair | 3.5 min |
| **B** | 34 trials, shuffled per repetition | ≈ 25 min (12.3 min dwelling + ≈ 22 s overhead per trial) |
| transit | 60 → 25 mm | 0.6 min |
| **A** | 6 trials | ≈ 6.6 min |
| transit | 25 → 60 mm | 0.6 min |
| **L** | one 300 s dwell, mode R | ≈ 5.5 min |
| R3 | reference pair | 3.5 min |
| **Total** | | **≈ 45–46 min** |

### 4.3 Reduced version (if F0b showed no dwell effect)
B only, mode R only, T_d = 0.1, 2, 30, 120 s, 2 repetitions (≈ 12 min). It can be appended to another session (e.g. F3). This gives an upper bound on the dwell effect for the report.

---

## 5. Pre-run checks
- [ ] Seed logged; shuffled order printed in the table; every (T_d, mode) occurs the planned number of times
- [ ] Dwell timer verified (start/end timestamps in the firmware log match T_d to < 10 ms)
- [ ] All trials start from a wipe from below; mode R ramps down, mode C ramps up (direction check, as for F0b/F1)
- [ ] Total ≤ 52 min; scope record ≥ 60 min
- [ ] Temperature logging running (long dwells are sensitive to thermal drift)

---

## 6. Analysis

### 6.1 Pre-processing
- Drift and gain from R1/R3; **additionally** remove slow drift within long dwells using the thermistor / BME280 log if correlated (a 300 s dwell is long enough for thermal drift to show).
- Datum for each trial: the mean position over the 0.2 s before the dwell starts.

### 6.2 Creep during the dwell
- Carriage position x(t) during the dwell, against log(t − t₀).
- Fit x(t) = c₀ + c₁·log(t − t₀ + τ) (typical relaxation form) and report the creep after 1, 10, 100 s.
- Compare with the reference-pair drift: creep is a local relaxation right after stopping; drift is slow and global.

### 6.3 Mode R: reversal after the dwell
- Lost motion ℓ(s) = s − Δx against travel since the reversal, for each T_d, overlaid (log–log, as in the redesigned V1).
- Read for each T_d:
  - deadband (where the curve leaves ℓ = s),
  - pre-sliding length s₉₀,
  - **peak lost motion** (an overshoot above β₀ that relaxes back would indicate a stiction peak),
  - plateau (should stay at β₀ if only static friction ages).
- Effect size: Δℓ_peak(T_d) = ℓ_peak(T_d) − ℓ_peak(0.1 s).

### 6.4 Mode C: breakaway after the dwell
- Lag e against travel since restart, for each T_d.
- Read the **breakaway overshoot** (extra lag before the carriage starts following) and any **jump** (sudden step in x, i.e. stick-slip at restart).
- Effect size: Δe_break(T_d) relative to T_d = 0.1 s.

### 6.5 Dependence on dwell time
- Plot the effect sizes against log₁₀(T_d) with ±d_min band.
- Fit Δ(T_d) = B·log₁₀(T_d / 0.1 s) (logarithmic ageing) and, as an alternative, a saturating form Δ_∞(1 − exp(−T_d/τ)).
- **Significance:** effect over the operational dwell range (e.g. 0.1 to 120 s) compared with 3·d_min.
- Position check: does A (25 mm) follow the same law as B?

### 6.6 Decision

| Result | Consequence |
|---|---|
| Effect < 3·d_min over 0.1–120 s | dwell can be ignored; the friction model stays rate- and time-independent; report the upper bound |
| Effect significant, logarithmic | add a dwell-dependent static friction term to the model (static level growing with log T_d); the compensation needs the rest time since the last move, which the host knows from the move plan |
| Effect significant but saturates within seconds | a fixed minimum dwell in operation (e.g. always rest ≥ τ) removes the variability more simply than modelling it |
| Creep larger than the effect itself | the stop position keeps changing during the rest; settling time, not ageing, is the issue (compare with the E5 settling curves) |

---

## 7. Outputs
- `f2_command_table.csv` (+ seed, hash), raw recording, firmware log
- `f2_creep.csv` (trial, T_d, t, x), `f2_reversal.csv` (T_d, s, ℓ), `f2_breakaway.csv` (T_d, s, e)
- `f2_summary.json`: creep at 1/10/100 s, Δℓ_peak(T_d), Δe_break(T_d), fitted ageing law and significance, position comparison
- Figures: overview (§8), and the companion plots D1–D4 (§9)

---

## 8. Render instructions: F2 overview figure

Same visual grammar as F0b and F1, adapted to the fact that **time at rest**, not amplitude, is the variable.

### Layout (full text width, ≈ 7.2 × 6.5 in)

```
(a) ribbon:  R1 | B rep1 | B rep2 | B rep3 | A | L | R3
    trajectory: commanded position relative to x_p (symlog), dwells as flat segments
(b) dwell coverage: rows B·R, B·C, A·R, L
(c) trial anatomy: one trial of each mode, aligned at the end of the dwell
```

### (a) Ribbon + trajectory
- **Ribbon, two tiers:** top tier the blocks (R1, B, A, L, R3); bottom tier the repetitions of B (rep 1, 2, 3).
- **y:** commanded position relative to the current test position, symlog with ±µm/mm labels and the shaded linear band, as in F0b/F1. Wipes appear as grey dips to −3 mm; post-dwell ramps as short slopes to ±6 µm.
- **Dwell segments** drawn on the zero line, **coloured by T_d** with a sequential green scale (light = 0.1 s, dark = 300 s). Include a small colour key: "dwell time".
- **Post-dwell ramps:** mode R dark slate going down, mode C dark slate going up, so direction alone separates the modes.
- Reference pairs purple dashed; transits grey.
- Because the dwells dominate time, the shuffled order is visible as a random pattern of short and long flat segments. Note it in the caption ("order shuffled per repetition").

### (b) Dwell coverage
- **x:** dwell time T_d, log, 0.05–500 s.
- **Rows:** B · reversal, B · continuation, A · reversal, L · long dwell.
- **Markers:** one marker per planned trial at its T_d, stacked vertically by repetition (so 3 dots = 3 repetitions). Filled for mode R, open for mode C.
- **Guides (labels at the top, offset from the lines):**
  - F0b dwell probe points (1 s, 60 s), dotted;
  - F1 virgin-ramp settle (2 s) and E5 ISO stop (2.5 s), dotted: "dwells already present in earlier experiments";
  - optional shaded band for the expected operational dwell range (to be agreed, e.g. 0.5–60 s).

### (c) Trial anatomy
- Two horizontal timelines, **aligned at the end of the dwell (t = 0)**, one per mode: grey wipe and settle → green dwell bar (draw two lengths, e.g. 2 s and 120 s, the long one broken with a // gap) → dark ramp.
- Above the ramp: "down: reversal (R)" or "up: continuation (C)".
- Small labels under each part: "creep measured here" under the dwell; "response measured here" under the ramp.

### Style
- Serif, panel letters only, findings in the caption; text ≥ 7 pt.
- Black-and-white safe: dwell length is visible as segment length, mode as ramp direction, markers filled/open.

### Caption (suggested)
> F2 command sequence. (a) Commanded position relative to the test position; each trial is a 3 mm wipe from below, a dwell (green, shade by duration) and a 6 µm post-dwell ramp, down (reversal, R) or up (continuation, C). The order of dwell times and modes is shuffled per repetition. (b) Dwell times covered at each position, one marker per trial; dotted lines mark the dwells already present in F0b, F1 and E5. (c) Anatomy of one trial: creep is measured during the dwell, the friction response during the ramp. ⟨n⟩ trials, ⟨t⟩ min.

---

## 9. Companion plots (design rules as for the redesigned F1 plots)

Each plot: a **mini command icon** with numbered points, **intuitive quantity**, a **ghost alternative** (what "no ageing" would look like), smooth illustrative curves until data exist, serif text, findings in the caption.

### D1 — Creep during the dwell
- **Icon:** wipe → stop ① → dwell ② … ③ end.
- **x:** time since stop (log, 0.01–300 s); **y:** carriage position relative to the stop (nm).
- Illustrative: a log-like rise of a few tens of nm, flattening.
- **Ghost:** a flat line at 0 ("rigid stop, no creep").
- Band: ±d_min.

### D2 — Reversal after dwell (mode R)
- **Icon:** wipe → dwell (two lengths shown) → reversal ramp.
- **x:** travel since reversal (log); **y:** lost motion (log–log, as in V1).
- Curves for T_d = 0.1, 2, 30, 120, 300 s with the green shade scale (light to dark).
- Illustrative effect: longer dwell → curve stays on the stuck diagonal slightly longer and shows a small peak above β₀ before settling.
- **Ghost:** all curves on one line (light grey), labelled "no ageing: dwell doesn't matter".
- Guides: ℓ = s diagonal, β₀ plateau.

### D3 — Breakaway when continuing (mode C)
- **Icon:** wipe → dwell → continuation ramp.
- **x:** travel since restart (linear, 0–6 µm); **y:** extra lag relative to the steady sliding lag (nm).
- Illustrative: a brief overshoot growing with T_d, then return to 0; optionally a small step for long dwells (stick-slip at restart).
- **Ghost:** flat at 0 ("no breakaway").

### D4 — Effect against dwell time (summary)
- **x:** T_d (log); **y:** effect size (nm): Δℓ_peak (filled, mode R) and Δe_break (open, mode C), with ±d_min band; position A as a different marker shape.
- Illustrative: points along a logarithmic trend.
- **Ghost:** flat line at 0, labelled "no ageing".
- Shade the operational dwell range, so the reader sees directly whether the effect matters in practice.

For the paper, D4 alone (plus one panel of D2) may be enough; D1–D3 are for the presentation and the appendix.

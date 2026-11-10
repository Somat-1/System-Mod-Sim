# F0b — Friction scoping session

> **Question:** How large are the friction effects in this stage, over what distances and times do they act, and can the interferometer resolve them? The answers set the amplitudes, speeds and dwells of F1–F3.

F0b is deliberately coarse: one position, few repetitions, wide grids. It doesn't identify the friction model. It measures the **scale** of every effect the model will need, so the main experiments are dense where it matters and don't waste time where nothing happens.

---

## 1. What F0b measures

| Block | Measures | Decides |
|---|---|---|
| **N** noise floor | precision and drift at rest at the motor end, mid-travel and far end | the smallest effect worth measuring and whether it depends on stage position |
| **L** coarse reversal ladder | lost motion vs stroke, from 1 µstep to 3 mm | the pre-sliding length and the saturation level → F1 amplitude grid |
| **S** slow ramps | what happens microstep by microstep after a reversal and in continuation | the deadband, stick-slip, the usable resolution at 1/256 |
| **D** dwell probe | whether 1 s vs 60 s at rest changes the next move | whether F2 is needed, and its range |
| **V** speed probe | whether the steady lag differs between 0.01 and 1 mm/s | the F3 range, or whether sliding can be treated as speed-independent |

### 1.1 Resolution reference

At 1/256: **1 µstep = 2 mm / (200 · 256) = 39.1 nm**. One electrical cycle = 1024 µsteps = 40 µm. 1 mm = 25 600 µsteps.

| Stroke (µsteps) | 1 | 2 | 5 | 13 | 26 | 51 | 128 | 256 | 512 | 1280 | 2560 | 12 800 | 76 800 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Stroke (µm) | 0.04 | 0.08 | 0.20 | 0.51 | 1.0 | 2.0 | 5.0 | 10 | 20 | 50 | 100 | 500 | 3000 |

The interferometer floor (≈ 1–2 nm for a 0.5 s window mean) is well below one µstep. If single steps are invisible, the stage isn't moving: that is a mechanical result, not a metrology limit. Block N confirms the floor at three positions spanning the stage travel.

---

## 2. Configuration

| Item | Setting |
|---|---|
| Driver | TMC2209, **MRES 1/256 for the whole session** (MSCNT 0–1023), SpreadCycle, interpolation off, IRUN = IHOLD (as E5, 400 mA class) |
| Do not change MRES mid-session | it changes the step-to-position mapping and the electrical phase reference |
| Step generation | RMT, precomputed segments; constant-rate trains from ~1 step/s to 25.6 kHz (1 mm/s); verify ≥ 25.6 kHz before the session |
| Transit speed | 1 mm/s, 40 mm/s² (25.6 kHz at 1/256); transits are not analysed |
| Test position | **x₀ = 60.000 mm**, on the 1024-µstep grid (u₀ = 1 536 000), so every block starts at the same electrical phase |
| Memory wipe | 3 mm (76 800 µsteps) approach in the stated direction at 1 mm/s, then 2 s settle |
| Reference pair | 15 / 105 mm, both directions, at start, middle, end (as E5) |
| Recording | IDS continuous, TwinCAT scope at 1 kHz, **record length ≥ 60 min**; IDS logger file as backup |
| Logging | row index, segment type, start/end timestamps, `u_cmd`, MSCNT at every segment start; BME280 and thermistors if fitted |
| Warm-up | as E0 v4 / E5, then homing and centring, before the recording starts |

### 2.1 Command table

Same principle as E5: the host generates the full sequence as rows, the board replays them.

| column | meaning |
|---|---|
| `row` | execution index |
| `block` | R / N / L / S / D / V |
| `seg` | `move` (trapezoidal to target), `ramp` (constant step rate to target), `step` (n single steps with a dwell after each), `dwell` |
| `u_target` | target step count |
| `rate_hz` | step rate for `ramp`; ignored otherwise |
| `dwell_s` | dwell after the segment (for `step`: after each step) |
| `tag` | stroke, cycle, half (up/down), dwell class, speed class |

---

## 3. Sequence

Total recording ≈ 48.2 min, within the 58 min gate.

| Block | Content | Time |
|---|---|---|
| R1 | reference pair | 1 min |
| **N** | noise floor at three stage positions | 11.7 min |
| **L** | coarse reversal ladder, up then down the stroke list | 10 min |
| R2 | reference pair | 1 min |
| **S** | slow ramps | 12 min |
| **D** | dwell probe | 5 min |
| **V** | speed probe | 5 min |
| R3 | reference pair | 1 min |
| — | transits | ≈ 3 min |

### N — Noise floor (11.7 min)

Approach each position from below with the standard 3 mm memory wipe, settle,
then hold with the motor energised at hold current and no steps.

| Hold | Position | Duration | Purpose |
|---|---|---:|---|
| **N-end1** | 15 mm (motor end) | 2 min | floor near the motor: shortest beam path or motor-side stiffness, plus any motor heat effect |
| **N-mid** | 60 mm (x₀) | 5 min | full characterisation: Allan deviation up to 5 s, drift rate and spectrum |
| **N-end2** | 105 mm (far end) | 2 min | floor at the opposite end: longest or shortest beam path, depending on where the head sits |

Run them in the table order. The two 45 mm transfers and three approaches account
for the remaining ≈2.7 min of the block.

### L — Coarse reversal ladder (10 min)

Anchored loops above x₀: the stage oscillates between x₀ and x₀ + s.

1. Memory wipe from below to x₀.
2. For each stroke s in the list, **increasing**: `x₀ → x₀+s → x₀ → x₀+s → x₀` (the first up-half is a continuation after the wipe; the other three are reversals), with 0.5 s pauses at each turning point.
3. Repeat the whole list **decreasing**, without a new wipe.

Speeds:

- s ≤ 5 µsteps: `step` segments, 0.5 s dwell after each step.
- s > 5 µsteps: `ramp` with each half-cycle taking ≈ 4 s, i.e. rate = s / 4 s, capped at 0.75 mm/s (19.2 kHz) for the 3 mm stroke.

Increasing order matters: each larger stroke wipes the memory of the smaller ones, so every stroke starts from a defined state. The decreasing pass checks for history effects in the ladder itself.

### S — Slow ramps (12 min)

All ramps: 20 µm (512 µsteps) at **0.25 µm/s (6.4 µsteps/s)**, 80 s each. At 1 kHz that's about 150 samples per microstep.

| Ramp | Preparation | Then |
|---|---|---|
| S1 | wipe from below to x₀ | **reverse**: ramp down 20 µm |
| S2 | wipe from above to x₀ | **reverse**: ramp up 20 µm |
| S3 | wipe from below to x₀ | **continue**: ramp up 20 µm |
| S4 | wipe from above to x₀ | **continue**: ramp down 20 µm |

Run S1–S4, then S1–S4 again (two repetitions each).

### D — Dwell probe (5 min)

Each trial: wipe from below to x₀ → dwell T_d (recorded) → **reverse**: ramp down 5 µm at 0.25 µm/s (20 s).

Order **T_d = 1 s, 60 s, 60 s, 1 s** (ABBA).

### V — Speed probe (5 min)

1. 1 mm segment x₀ … x₀ + 1 mm, traversed up and down at **1 mm/s**, 2×.
2. The same segment up and down at **0.01 mm/s** (100 s per pass), 1×.
3. Again at 1 mm/s, 2× (ABA: fast–slow–fast).

---

## 4. Pre-run checks

- [ ] RMT produces the 25.6 kHz transit rate and the 6.4 Hz ramp rate cleanly (scope or pulse count check)
- [ ] x₀ and all segment start points on the 1024-µstep grid where required; MSCNT logged at each start
- [ ] Command table: stroke list, orders and dwells as specified; total duration ≤ 50 min
- [ ] Scope record length ≥ 60 min (avoid the E5 stitch)
- [ ] Dry run of the first ~30 rows: firmware echoes `u_cmd` = table

---

## 5. Analysis

### 5.1 Pre-processing

- **Drift:** piecewise-linear from R1–R3 (↑/↓ mean at the reference pair).
- **Per-segment datum:** e(t) = x_meas(t) − x_meas(t_start) − [u_cmd(t) − u_cmd(t_start)]/U, i.e. tracking error relative to the start of each segment. U = 25 600 µsteps/mm.
- **Geometry:** for L strokes ≥ 100 µm and the V block, subtract Ĝ along the path, or (V) use the half-difference of up and down at the same positions. Below 100 µm, geometry changes by < 0.6 µm over the stroke and < 0.06 µm below 10 µm, so it can be ignored for the small strokes.
- **Validity:** MSCNT equal at every return to x₀ (no lost steps); x at x₀ consistent across the L ladder within the drift.

### 5.2 N — Noise floor

- σ of window means for window lengths 0.05, 0.1, 0.5, 1, 2, 5 s (Allan deviation).
- **Minimum detectable step:** d_min ≈ 3 · √2 · σ(0.5 s).
- Drift rate at rest (nm/s) after removing the reference-pair drift.
- Compare N-end1, N-mid and N-end2: position dependence of the floor across the full travel.
- Use N-mid for the primary Allan-deviation, drift-rate and spectral estimates; use the two end holds as shorter cross-checks.
- Spectrum (Welch) for dominant vibration lines.

### 5.3 L — Lost motion vs stroke

For each stroke s and each half-cycle:

- **Measured travel** Δx over the half-cycle (end window mean minus start window mean).
- **Transmission** τ(s) = Δx / s: 0 = the carriage doesn't move, 1 = it follows.
- **Lost motion** ℓ(s) = s − Δx.

Plot τ(s) and ℓ(s) against s on a log axis, separately for the **continuation half** (first up after the wipe) and the **reversal halves**.

Expected shape:

- reversal halves: ℓ grows with s and saturates at **ℓ_sat ≈ β₀ ≈ 1.757 µm**;
- continuation half: ℓ ≈ 0 apart from small-stroke stick effects.

Read off:

- **Saturation level** ℓ_sat. Compare with β₀ from E2–E5.
- **Pre-sliding length** s₉₀: the stroke where ℓ reaches 90 % of ℓ_sat.
- **Smallest stroke with τ > 0** (beyond d_min): the effective resolution after a reversal.
- Increasing vs decreasing pass: any difference = history effect beyond what the stroke order explains.

### 5.4 S — Slow ramps

Plot e against commanded travel since the segment start (0–20 µm):

- **After reversal (S1, S2):** e starts at 0 and falls as the carriage lags: this is the **virgin curve** at fine resolution. Read the **deadband** (commanded travel before the carriage moves by more than d_min), the initial slope (pre-sliding stiffness, in displacement units), and whether e saturates within 20 µm.
- **Continuation (S3, S4):** e should stay near constant. Look for **stick-slip** (sawtooth: e drifts, then jumps), and for periodic structure with a 40 µm period (electrical cycle; only partly visible within 20 µm).
- **Step response within the ramp:** fold e on the microstep index: does each microstep produce a visible step, or only every few?
- Up vs down: symmetry of the deadband and virgin curve.

### 5.5 D — Dwell

Lost motion and initial deadband after 1 s vs 60 s (mean of 2 each, ABBA); creep during the 60 s dwell.

### 5.6 V — Speed

Steady lag = half-difference of up and down at the same positions, at 0.01 vs 1 mm/s. A difference of more than a few d_min means friction is speed-dependent in this range.

---

## 6. Decisions for F1–F3

| F0b result | Sets / decides |
|---|---|
| d_min (block N) | smallest F1 stroke and S-ramp resolution: ≥ 3–5 × d_min |
| s₉₀, ℓ_sat (block L) | F1 stroke grid: dense from about s₉₀/10 to 3 × s₉₀, sparse elsewhere; ℓ_sat as a cross-check with β₀ |
| Deadband, stick-slip (block S) | whether single 1/256 steps are usable; what "small move" can promise; whether the friction model needs a stick-slip term |
| 1 s vs 60 s difference (block D) | F2 run (if > 3 × d_min) or skipped; if run, a dwell range around the observed effect |
| 0.01 vs 1 mm/s difference (block V) | F3 velocity grid (dense where it changes) or a reduced F3 |
| Increasing vs decreasing ladder | whether F1 needs explicit memory wipes between strokes |

---

## 7. Outputs

- `f0b_command_table.csv` (+ hash), raw IDS recording, firmware log
- `f0b_noise.json`: σ(τ), d_min, drift rate
- `f0b_ladder.csv`: s, half, Δx, τ, ℓ
- `f0b_summary.json`: ℓ_sat, s₉₀, deadband ↑/↓, stick-slip flag, dwell effect, speed effect
- Figures:
  1. Noise: Allan deviation with d_min marked
  2. Ladder: τ(s) and ℓ(s) on log-s, reversal vs continuation, increasing vs decreasing
  3. Ramps: e vs commanded travel, S1–S4 overlaid
  4. Dwell and speed probes: two small panels
- A short "scales card" for the friction chapter:

| Quantity | Value |
|---|---|
| Minimum detectable step d_min | … nm |
| Saturation lost motion ℓ_sat | … µm (β₀ = 1.757 µm) |
| Pre-sliding length s₉₀ | … µm |
| Deadband after reversal ↑ / ↓ | … / … µm |
| Smallest stroke that moves after reversal | … µsteps |
| Stick-slip in continuation | yes / no, amplitude … |
| Dwell effect (60 s vs 1 s) | … µm |
| Speed effect (0.01 vs 1 mm/s) | … µm |

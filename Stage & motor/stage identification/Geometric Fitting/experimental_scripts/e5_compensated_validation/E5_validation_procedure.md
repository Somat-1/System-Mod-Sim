# E5 — Compensated Validation (ISO 230-2, paired) · v2

> **Question:** Does the frozen geometric + reversal model, applied as a command pre-correction, reduce the stage's positioning error on an independent session — and is what remains explained by the day-to-day slow-slope change seen in E3?

Everything validated so far (cross-validation on E2, replay on E0/E1/E4, transfer to E3) is **offline**: the model was subtracted from measured data. E5 is the first test where the correction goes into the **command**. It therefore also tests what offline replay cannot:

- microstep quantisation of the correction (±0.31 µm at 1/16),
- evaluating Ĝ at x\* instead of x_c,
- the motor microstep error that returns once commanded targets leave the 64-µstep electrical grid,
- whether the step-count origin (and so θ) survives a fresh homing — E3 says yes (k = 1 phase to 0.3°); E5 confirms it in closed form.

---

## 0. Changes since v1 (driven by E3)

| # | Change | E3 evidence |
|---|---|---|
| 1 | **Model = H = 2, varying coefficients, fitted on E2 only** (was H = 4) | k = 3, 4 not reproduced; E3's own selection picks H = 2; H2 transfers as well as H4 (0.515 vs 0.522 µm). H = 2 was a candidate in E2's original CV grid, so no E3 data enters the model. |
| 2 | **Two reference positions, each approached from both directions** (was one, unidirectional) | (a) E3 residual is dominated by a slow slope change between days (≈ −0.65 µm at 10 mm, +0.33 µm at 80 mm, ≈ 14 ppm) — a single reference cannot see a gain change. (b) Reversal in E2/E3 was confounded with time; measuring B at fixed positions at three moments separates time from position. |
| 3 | **Primary metric is the rms of the compensated direction-averaged residual**, compared to E3's transfer rms; ISO ranges kept as the standard report | rms is directly comparable to E3 (0.522 µm) and less sensitive to one bad target than a range. |
| 4 | **Acceptance criteria re-derived** from E3 transfer error + quantisation (were based on the offline E0 replay) | The E0 replay (M = 1.335 µm) was same-day-like; E3 shows the realistic cross-day floor. |
| 5 | **Secondary analysis: frozen geometry + session gain/offset from the reference pair** | Tests the operational model suggested by E3: fixed wobble, per-session slope. |
| 6 | **Edge targets (x < 5 mm or > 110 mm) flagged** and reported with/without | Model weakly constrained there; k = 1 diverged above 110 mm in E3. |
| 7 | β₀ stays at the E2 value (−1.757 µm), unchanged | E3's −1.623 µm is partly a time artefact of the mirrored passes; the time-clean far-end value in both sessions is ≈ −1.7 to −1.75 µm. Keeping the frozen value keeps E5 independent. Expect B_mean,C ≈ 0 to +0.13 µm. |

---

## 1. Design principles

| Choice | Reason |
|---|---|
| Same 24 targets, same ISO cycle, same settings as **E0 v4** | Direct before/after comparison against the existing baseline. |
| **Uncompensated (U) and compensated (C) runs in the same session**, interleaved | Day-to-day differences hit both modes equally; the comparison is paired. |
| **ABBA ordering** U C C U U C C U U C | Linear drift over the session cancels between modes. |
| **Model frozen and hashed before the run** | E5 stays independent. Any later refit (e.g. pooled E2 + E3) gets reported separately. |
| Corrections **precomputed on the host**, sent as a command table | No model evaluation on the ESP32; the exact commanded step count is logged per stop. |
| **Off-grid block (B3)** | ISO targets only sample a narrow band of electrical phase after correction; B3 covers the full cycle. |
| **Bidirectional reference pair** at 3 moments | Session gain/offset for the secondary analysis; reversal vs time at fixed positions. |

---

## 2. Model and compensation law

**Frozen model:** `geom_model_E5.json`

- Fitted on **E2 only** (seed 5315, D2 drift variant).
- H = 2, position-varying C₁(x), C₂(x); K_s and K_C as selected in E2's CV (K_s = 2, K_C = 1 unless the E2-only refit at H = 2 selects otherwise — decide once, record it).
- No constant k ≥ 3 terms.
- Record: coefficients, knot positions, x̄, fit span, and the file's SHA-256.

**Compensation law** for target x\* approached in direction d (+1 up, −1 down):

$$
x_c = x^* - \big[\hat G(x^*) - a_0\big] - \tfrac{d}{2}\,\beta_0, \qquad u_c = \operatorname{round}(x_c \cdot U), \qquad \beta_0 = -1.757\ \mu\text{m}
$$

- a₀ dropped (arbitrary datum).
- Quantisation residual q = u_c/U − x_c, |q| ≤ 0.3125 µm, logged per stop.

**Sign check:** with β₀ < 0 the up-approach command lands **higher** than the down-approach command for the same target, by ≈ 1.76 µm (± one µstep). Up-approaches stop short, so they must be sent further.

---

## 3. Session structure (one 60 min window)

| Block | Content | Stops | Time (est.) |
|---|---|---|---|
| **B0** | Homing → warm-up as in E0 v4 → re-home → centre | — | per E0 practice |
| **R1** | Reference pair, bidirectional (§3.1) | 4 | ~1 min |
| **B2a** | ISO runs 1–5: **U C C U U** | 240 | ≈ 240 × t_cycle |
| **R2** | Reference pair, bidirectional | 4 | ~1 min |
| **B2b** | ISO runs 6–10: **C C U U C** | 240 | ≈ 240 × t_cycle |
| **R3** | Reference pair, bidirectional | 4 | ~1 min |
| **B3** | Off-grid block, compensated only (§3.2) | 128 | ≈ 128 × t_cycle |
| **R4** | Reference pair, bidirectional | 4 | ~1 min |

`t_cycle` = move + settle + record of E0 v4 (≈ 2.9 s/stop in E3). Estimate: B2 ≈ 23 min, B3 ≈ 7 min, R1–R4 ≈ 4 min → **≈ 34 min**. Fill in the real `t_cycle` and keep ≥ 8 min margin to the 58 min gate. If tight: B3 down to one repetition first, never B2.

### 3.1 Reference pair (R blocks)

- Two positions inside the well-constrained span, far apart: **x_R1 ≈ 15 mm, x_R2 ≈ 105 mm**, both on the 64-µstep grid, both **uncompensated** (nominal command).
- At each: pre-approach 3 mm below → stop (up), overrun 3 mm above → stop (down). Same settle/record windows as the ISO stops.
- Order within a block: R1↑ R1↓ R2↑ R2↓ (same every block).
- Uses:
  - **Drift**: mean of ↑/↓ at each position vs time → D2 anchors (now four instead of three).
  - **Session gain/offset** for the secondary analysis (§5.6): the frozen model predicts Ĝ at both points; the measured-minus-predicted difference between R2 and R1, divided by their separation, is the session gain correction.
  - **Reversal vs time**: B_R(t) = e↑ − e↓ at each reference, four moments → separates time from position (the E2/E3 confound).

### 3.2 B3 — off-grid electrical-phase block

- 32 targets spread over **5–110 mm** (new seed, none coinciding with E0/E2/E3 targets).
- Target k at offset **o_k = −32 + 2k µsteps** from the 64-µstep grid → full electrical cycle in 2-µstep increments.
- Compensated with the same law, approached up and down (pre-approach/overrun as in E2/E3), **2 repetitions**.
- Purpose: residual vs electrical offset. Expected microstep contribution ≈ 0.14 µm rms (motor layer).

### 3.3 Motion details (unchanged from E0 v4 unless stated)

- 1/16 microstepping, SpreadCycle, interpolation off, IRUN/IHOLD as E0 v4, RMT step generation, AquadB 120 ns.
- Pre-approach and overrun computed from the **nominal** target. Corrections (≤ ~6 µm) never change the approach direction.
- Same velocity/acceleration, settle and record windows as E0 v4.
- Same homing routine and step-count origin as E2/E3.

---

## 4. Implementation

### 4.1 Host: `e5_make_table.py`

1. Load `geom_model_E5.json`, print its SHA-256.
2. Load the E0 v4 target list; generate B3 targets (record seed); define R positions.
3. For every stop compute `x_c`, `u_c`, `q`.
4. Write one row per stop **in execution order**:

| column | meaning |
|---|---|
| `row` | execution index (echoed by firmware) |
| `block` | R / B2 / B3 |
| `run` | ISO run (1–10), B3 rep, or R block (1–4) |
| `mode` | U / C |
| `dir` | +1 / −1 |
| `target_id` | ISO / B3 / R index |
| `u_nom` | nominal step count |
| `u_cmd` | commanded step count (= u_nom for U and R) |
| `corr_um` | (u_cmd − u_nom)/U in µm |
| `q_um` | quantisation residual (0 for U/R) |
| `theta_deg` | screw angle of u_cmd |
| `el_off` | signed offset of u_cmd from the 64-µstep grid (−32…31) |
| `edge` | 1 if x\* < 5 mm or > 110 mm |

5. Export as a C array / streamed file, same mechanism as E3.

### 4.2 Firmware

- Executes the table row by row; no model logic.
- Prints at start: model hash, B3 seed, table hash, `MSCNT` at home.
- At each stop logs `row`, `u_cmd` (echo), **`MSCNT`**, move-end and record-window timestamps.
- Timing records in RAM, printed after `SEQ_DONE`.

### 4.3 Pre-run checks (today)

- [ ] `corr_um` vs x for C rows looks like −Ĝ(x) (wobble on a slow trend, ≈ ±4 µm), with no k ≥ 3 content.
- [ ] Up-minus-down command per target = **+1.757 µm ± 0.625 µm** (sign check).
- [ ] `max |corr_um|` < ~6 µm; `max |q_um|` ≤ 0.3125 µm.
- [ ] Hand-check 3 targets against the offline Ĝ (one near each end, one mid-travel).
- [ ] U and R rows: `u_cmd mod 64 == 0`. U rows identical to the E0 v4 command sequence.
- [ ] B3 `el_off` covers −32…30 evenly; B3 targets inside 5–110 mm.
- [ ] Edge flags set for the E0 targets outside 5–110 mm (note which and how many).
- [ ] Predicted duration < 50 min.
- [ ] Dry run: firmware echoes `u_cmd` equal to the table for the first ~20 rows.

---

## 5. Analysis

Pre-registered: run as written; do not change the model or criteria after seeing E5.

### 5.1 Pre-processing (as E3)
- Blind measured-hold windowing (no ESP-log window edges).
- Drift: D2 piecewise-linear from the ↑/↓ mean at the reference pair, four anchors (R1–R4). Report D0/D1/D2 differences as in E3.
- Reject record windows > 0.3 µm p-p or with excessive slope; report counts.
- Deviation always relative to the **nominal** target x\*, for both modes.
- Validity table as in E3: stops matched, MSCNT constant across U/R stops, MSCNT matches `el_off` for C/B3 stops, no saturated increments, drift datum→end.

### 5.2 Primary: compensated residual rms
On the direction-averaged C means at the 24 targets (5 runs):

- **rms_C** (after removing one constant offset), all targets and interior-only (edge targets excluded).
- Compare with: E3 transfer rms **0.522 µm**, E3 own-fit floor **0.369 µm**, and rms_U of the same session.
- Expected: rms_C ≈ √(0.52² + q²) with q the direction-averaged quantisation (~0.1–0.15 µm) → **≈ 0.54 µm**.

### 5.3 ISO 230-2, U vs C
E, M, A, B, B_mean, R, R↑, R↓ for U and C separately.

Expected pattern:
- E, M, A, B drop strongly; B_mean moves from ≈ −1.75 µm to ≈ 0 … +0.13 µm.
- R drops (it contains |B_i|).
- **R↑, R↓ stay roughly unchanged** — compensation cannot remove random scatter.
- M_C is a range of 24 means, ≈ 4× rms_C → expect **M_C ≈ 2 µm**, not the 1.335 µm of the offline E0 replay.

Also: U session vs E0 v4 (is it the same stage as the baseline).

### 5.4 Achieved vs commanded correction
Per target and direction: Δ = x̄_C − x̄_U against `corr_um`. Fit Δ = k·corr + c.
- k ≈ 1: amplitude right. k < 1: overcorrection / partial phase error. k ≈ −1: sign error.
- Scatter about the line ≈ model error; compare with 0.5 µm.

### 5.5 Residual structure (C runs)
- vs **x**: fit a slope. Compare with E3's ≈ 14 ppm day-to-day slope. A slope of similar size → the same mechanism; check the session temperature log.
- vs **θ**: k = 1, 2 amplitudes C/U. Leftover k = 1 with the U phase → amplitude mismatch; shifted phase → θ origin moved (check home `MSCNT`, step count at switch).
- **Quantisation-corrected residual** (subtract `q_um`): if clearly better, finer microstepping for corrections is worth testing.

### 5.6 Secondary: session gain/offset recalibration
Operational model suggested by E3 — fixed geometry, per-session slope.

1. From R1–R4 (direction-averaged, drift-corrected): r_j = measured − Ĝ(x_Rj), j = 1, 2.
2. Session gain correction g = (r₂ − r₁)/(x_R2 − x_R1); offset from their mean.
3. Apply g to the C residuals **post hoc** (the command was not changed) → **rms_C,recal** and M_C,recal.
4. Report the gain correction in ppm and compare with E3's slope change.

If rms_C,recal approaches ≈ 0.40 µm, the per-session slope is the dominant cross-day term, and a two-point calibration at power-up is the practical recommendation.

### 5.7 Reversal
- B_C per target: mean and position trend. Expected ≈ 0 to +0.13 µm, flat.
- **B_R(t)** from R1–R4 at both references: does reversal change with session time at a fixed position? Same trend at both references → time effect (warm-up / friction state). Different → position.
- Compare B_U (ISO, short up/down separation) with the far-end E2/E3 value (≈ −1.7 to −1.75 µm).

### 5.8 B3 — electrical phase
Residual vs `el_off`; fit 1, 2, 4 harmonics per electrical cycle; compare with ≈ 0.14 µm rms from the motor layer. Flat → motor LUT not needed at this operating point.

### 5.9 Settling (thesis metric, secondary)
|error| at 0.25, 0.5, 1.0, 2.0 s after the last step pulse; median and p95; U vs C.

### 5.10 Figures (to produce)
1. ISO deviation vs x\*, U and C side by side, both directions, mean ± 2s, shared y-scale; edge zones shaded.
2. Achieved vs commanded correction with y = x and fitted line.
3. Compensated residual along the travel (frozen, and after session recalibration), binned means ± 2σ.
4. Error budget as a table: rms_U, rms_C, rms_C,recal vs E3 transfer and own-fit floor; ISO metrics U/C.
5. B3 residual vs electrical offset.
6. Reversal: B_C vs x and B_R vs session time.
7. Settling |error| vs time, U vs C.

Black-and-white safe: U open markers / dashed, C filled / solid.

---

## 6. Acceptance criteria (fixed before the run)

| Check | Pass | If failed |
|---|---|---|
| **rms_C (primary)** | ≤ 0.65 µm (interior) | Model doesn't transfer to the commanded stage beyond the known cross-day term |
| rms_C / rms_U | ≤ 0.45 | Compensation removes less than expected (E3: 0.35) |
| M_C | ≤ 2.2 µm | As above |
| E_C | ≤ 2.6 µm | As above, or β₀ off |
| \|B_mean,C\| | ≤ 0.3 µm | Reversal correction off (sign or day-to-day β₀) |
| R↑_C, R↓_C | within ±25 % of U | Correction adds scatter (quantisation, microstep) |
| Slope k (§5.4) | 0.9 – 1.1 | Amplitude/phase mismatch |
| k = 1 amplitude in C residual | < 25 % of U | θ origin / phase problem |
| U session vs E0 v4 | M within ±0.7 µm | Stage changed since baseline |
| *Secondary:* rms_C,recal | ≤ 0.50 µm | Cross-day term is not mainly slope |

Margins are derived from E3 (transfer rms 0.522 µm) plus quantisation; adjust once if needed, but before the run.

---

## 7. Failure modes — what to look at first

| Symptom | Likely cause | Check |
|---|---|---|
| Error roughly **doubles** after compensation | Sign error | Up/down command difference sign; `corr_um` vs Ĝ |
| k = 1 remains, same amplitude, phase shifted | Step-count origin moved | Home `MSCNT`, step count at switch, θ of `u_cmd` |
| Up and down split by a constant | β₀ changed | B_U this session vs −1.757; B_R(t) |
| Slow tilt remains | Session gain differs (thermal/air) | §5.6 gain; temperature log; compare with E3's ≈ 14 ppm |
| Large errors only near ends | Spline edge | Interior-only metrics vs all |
| R↑/R↓ worse in C | Off-grid microstep error or quantisation | B3 vs `el_off`; quantisation-corrected residual |
| One target bad in all C runs | Local feature (72–80 mm, 95 mm) | Compare with E4 windows at that x |
| MSCNT not constant across U/R stops | Lost steps | Discard the affected segment; check current/speed |

---

## 8. Outputs

- `geom_model_E5.json` + hash; `e5_command_table.csv` + hash
- Raw IDS recording + firmware log
- `e5_iso_results.json` (U and C), `e5_per_target.csv`, `e5_reference.csv` (R blocks)
- Figures 1–7
- Verdict paragraph: which criteria pass, primary rms vs E3, whether session recalibration closes the gap → decision to freeze the geometric layer and move to the friction layer.


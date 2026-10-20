# Microstep-Resolution Test — Execution Procedure

**Goal:** select the operating MRES (microstep resolution) for the full identification sequence (Steps 0–6), based on measured settling error, repeatability and small-move response of the stage.
**Driver mode:** SpreadCycle only (StealthChop already excluded).
**Resolutions tested:** 1/16, 1/32, 1/64, 1/128, 1/256 (1/256 = anchor, run in both sessions).
**Sensor:** IDS 3010 interferometer.

**Design principle:** every move and staircase span is defined in **physical units** (full steps, FS, or mm), never as a fixed microstep count. That way all MRES are compared over identical displacements.

---

## 1. Parameters to fill in before the first session

| Symbol | Meaning | Value |
|---|---|---|
| L | screw lead [mm] | ___ |
| s_FS = L/200 | stage travel per full step [µm] | ___ |
| δ_MRES = s_FS/MRES | travel per microstep [nm] | 1/16: ___ · 1/32: ___ · 1/64: ___ · 1/128: ___ · 1/256: ___ |
| R | reversal lost motion from existing data [µm] | ___ |
| S_fine = max(¼ s_FS, 3R) | B3a span (fine tier) [µm] | ___ |
| σ_floor | interferometer floor at a 300 ms average (Allan minimum) [nm] | ___ (≈ low single-digit nm) |
| v | approach velocity [mm/s], **identical for all MRES** | ___ |
| f_step = v / δ_1/256 | step rate at 1/256 — check against driver/controller limit | ___ |
| x_start | start position (centre region of travel) | ___ |
| x_ref | drift reference position | ___ |

**Resolvability check, per MRES:**
- δ_MRES ≥ 3 σ_floor: single microsteps are resolved directly.
- δ_MRES < 3 σ_floor: single steps are not resolved. Evaluate B3a statistically (slope and scatter over the staircase), not step by step.

> If f_step at 1/256 exceeds what the step generator can deliver, lower v for **all** MRES rather than only for 1/256.

## 2. Settings held fixed across all blocks

- SpreadCycle (`en_spreadCycle = 1`), identical TOFF / HSTRT / HEND / TBL.
- I_RUN, I_HOLD and TPOWERDOWN identical to the intended operating values. Record them.
  - The hold-current reduction changes motor stiffness and therefore the friction-induced offset. Always take the settled-position window **after** the TPOWERDOWN transition, so the measurement represents the real holding state.
- `intpol` (MicroPlyer interpolation): pick one setting and keep it for all MRES. Record it.
- Same physical target region, same approach velocity in mm/s, same acceleration profile.
- Interferometer environmental compensation ON; log T and p.

## 3. Session plan (≤ 60 min each)

| | Block 1 | Block 2 | Block 3 |
|---|---|---|---|
| **Session 1** | 1/64 | 1/16 | **1/256 (anchor)** |
| **Session 2** | 1/128 | **1/256 (anchor)** | 1/32 |

- Randomize the MRES order per session (the table is one valid example). The anchor must appear in both sessions.
- Timing per session: warm-up (~6 min) → 3 blocks (~16.5 min each) → end reference (~2 min), about 57.5 min in total.
- B4 duration depends on v, because of the 5 mm and 20 mm moves. If a block overruns, reduce the repeats of the **operational-length moves only** (from 4 to 3); keep all short-move conditions intact.

### Warm-up and checks (start of every session, ~6 min)
1. Power the driver at I_RUN and run the stage back and forth over ±5 mm at v for ~4 min, until the interferometer reading at x_ref is stable to within the noise floor.
2. Check interferometer signal quality and alignment; log T and p.
3. Record x_ref three times (approach from the same side) as the session's drift baseline.

## 4. Block procedure (repeat for each MRES)

### B1 — Configure and align (~1 min)
1. Move to a **full-step position** (MSCNT at a full-step value) **before** changing MRES. Switching off-grid shifts the microstep reference and causes a position jump.
2. Write MRES via UART (CHOPCONF.MRES) and read it back to verify.
3. Wait 2 s, then re-zero the position reference on the interferometer.
4. Log MRES, MSCNT, T, p and a timestamp as the block header.

### B2 — Condition (~1.5 min)
1. Run 3 strokes of ±2 mm at v to erase the friction history.
2. End with an approach to x_start **from the negative side**, always the same side in every block.
3. Dwell 2 s.

### B3 — Two-tier staircase (~4.5 min)
Each sweep has a **fine tier** where pre-sliding happens (at the start, and directly after the reversal), followed by a **coarse tier** in gross sliding.

**Forward sweep (+):**
1. **B3a, fine tier:** single microsteps until S_fine is covered.
   - Count = S_fine / δ_MRES. At S_fine = ¼ FS this is 4 (1/16), 8 (1/32), 16 (1/64), 32 (1/128), 64 (1/256).
   - Few points at coarse MRES are expected. That is the finding, not a flaw: coarse settings cannot command smaller moves.
2. **B3b, coarse tier:** continue in increments of **1/16 FS** (= MRES/16 microsteps per increment: 1 at 1/16 … 16 at 1/256) until the total sweep reaches **4 FS** (one full electrical cycle).
   - Every MRES therefore samples the electrical cycle at the same 1/16-FS spacing.

**Reverse sweep (−):** reverse direction immediately and repeat the same pattern. B3a fine tier **first** (this captures the reversal dead zone), then B3b coarse tier back to x_start.

**Timing and recording:**
- Dwell **0.5 s** after every command. Record the mean position over the **last 300 ms** of each dwell. Longer averaging gains little, because correlated vibration dominates above τ ≈ 0.1 s.
- **2 repeats** of the forward + reverse pair.
- Order of magnitude at 1/256: (64 fine + 60 coarse) dwells × 0.5 s ≈ 62 s per sweep, × 2 directions × 2 repeats ≈ 4.1 min.

*Purpose:* the achieved increment per commanded microstep (slope), the reversal dead zone (lost motion), periodicity within the electrical cycle, and small-move scatter. **This data is reused for the GMS pre-sliding identification.**

### B4 — Settling matrix (~7.5 min)
1. Distances, in physical units, identical for every MRES:
   - **1/16, 1/4, 1, 4, 16, 64 FS**, plus
   - **5 mm and 20 mm**, which are operational move lengths. Adjust them to your real use case and check that they stay inside the travel from x_start.
2. For each move: start at x_start (approached from the fixed side), move by ±D, dwell, record, then return to x_start.
3. 8 distances × 2 directions × 4 repeats = **64 measured moves**, in **randomized order** (fix the seed per block and log it).
4. Dwell 2 s at each target, longer than TPOWERDOWN. Record the full trajectory, plus the settled value = mean over the final 500 ms.
5. Every 20 moves, visit x_ref (fixed-side approach) as a drift check.

### B5 — Creep / hold (~2 min)
1. 2 moves per direction with D = 4 FS.
2. Hold 30 s at the target and record the full trajectory.

*Purpose:* post-settle drift (creep) as part of the settling error, and whether it depends on MRES.

### End of session (~2 min)
- Visit x_ref 3× (fixed-side approach) and log T and p. The difference from the start baseline is the session drift.

## 5. Data to log (per sample)

| Field | Notes |
|---|---|
| t | timestamp |
| MRES | per block |
| commanded step count | absolute, in microsteps of the current MRES |
| MSCNT | TMC2209 microstep counter → electrical phase θ_e |
| x_IFM | interferometer position |
| T, p | environment |
| phase / move ID | B1–B5, tier (B3a/B3b), move index, D, direction, repeat, random seed |

Store one file per block, named e.g. `S1_B2_MRES16.csv`, with a header holding all fixed settings and the §1 parameters.

## 6. Analysis

**Pre-processing:** subtract the drift (linear interpolation between x_ref visits) and convert everything to physical units.

**Per MRES:**

| Metric | From | Definition |
|---|---|---|
| Small-move slope k_fine | B3a | achieved / commanded increment over the fine tier (ideal = 1) |
| Reversal lost motion | B3a (reverse) | commanded displacement before the stage moves after reversal |
| Sliding slope k_coarse | B3b | achieved / commanded over the coarse tier (should be ≈ 1) |
| Electrical-cycle periodicity | B3b | error vs θ_e (from MSCNT): amplitude of the 1st, 2nd and 4th harmonics |
| Small-move σ | B3a | scatter of the per-step increment (or of the residual about the slope if not resolvable) |
| Settling bias b(D, dir) | B4 | mean of (x_settled − x_target) |
| Repeatability σ(D, dir) | B4 | std of the settled positions |
| Creep | B5 | drift over the 30 s hold |
| Session consistency | anchor | 1/256 results, session 1 vs session 2 |

**Plots:**
1. b and σ vs D (log axis), one curve per MRES; mark the operational-length moves separately.
2. Staircase: achieved vs commanded position for each MRES, with the fine tier shaded.
3. k_fine, k_coarse and reversal lost motion vs MRES.

## 7. Decision rule

1. **Validity first:** the two 1/256 anchor runs must agree within their uncertainty. If they don't, session-to-session effects dominate and must be fixed before any MRES comparison is meaningful.
2. Choose the **coarsest MRES at which σ and |b| stop improving** (no significant gain at the next finer setting), weighted towards the **operational-length moves**, **provided** that:
   - its microstep size δ is small enough to resolve the planned corrections (δ ≲ target settling error), and
   - its step rate at v is within limits.
3. If a finer MRES performs *worse* (higher σ, or k_fine well below 1), its microsteps are below what the motor and friction can resolve. Do not choose it.
4. Record the chosen MRES, together with all fixed settings, as the configuration for Steps 0–6.

## 8. Checklist

- [ ] Parameters in §1 filled in: lead, R, S_fine, σ_floor; resolvability checked per MRES
- [ ] Step rate at 1/256 verified
- [ ] Operational move lengths fit inside the travel from x_start
- [ ] Fixed settings in §2 written to the driver and logged
- [ ] Warm-up completed; x_ref stable
- [ ] MRES changed only at full-step positions
- [ ] Fixed-side approach used before every B3/B4 start
- [ ] B3 reverse sweep starts with the fine tier
- [ ] Random seed logged per block
- [ ] Anchor (1/256) run in both sessions
- [ ] End-of-session reference recorded

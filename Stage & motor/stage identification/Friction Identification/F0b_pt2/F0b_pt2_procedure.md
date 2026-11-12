# F0b_pt2 — Reversal vs microstep resolution, and speed (rerun)

> **Questions:** (1) Does the set microstep resolution (MRES) change the reversal lost motion? (2) Does the approach speed change where the stage stops after a reversal? (3) Is F0b's 2.3 µm saturation, against β₀ = 1.76 µm from E2–E5 (run at 1/16), explained by MRES, settle time or speed?

Question 3 is why this session matters: F0b ran at 1/256 with 0.5 s pauses, while E2–E5 ran at 1/16 with 2 s settles at 2 mm/s. F0b_pt2 varies exactly those factors, one at a time, in one session.

---

## 0. Software fix for the timing problem

**Root cause.** The firmware logged when segments were *issued*, not when step pulses actually started and stopped, and the capture was aligned to the ESP log by velocity correlation (a few ms). At 1 mm/s, every millisecond of uncertainty is 1 µm of apparent lag, so mid-move readings were unusable.

**Fixes, best first:**

| Option | What | Mid-move lag usable? | Effort |
|---|---|---|---|
| **A. Record the commanded position on the TwinCAT clock** | Feed STEP and DIR into a Beckhoff counter terminal in pulse/direction mode (a second EL5101 channel, or an EL1512 / EL1502 up/down counter) and log it in the same scope as the interferometer. Commanded and measured position are then sampled in the same EtherCAT cycle. | **yes**: skew < 1 cycle jitter (µs), even at 1 mm/s | one terminal + wiring, no firmware change |
| **B. Timestamp the step train + sync pulse** | Log `t_first_edge` (RMT transmit start, corrected by the measured fixed latency) and `t_last_edge` (RMT tx-done ISR, `esp_timer_get_time()` in µs) per segment. Toggle a GPIO at each segment start into a Beckhoff digital input in the scope (EL1252 with timestamps, or any EL1xxx at 1 kHz). The RMT schedule is deterministic (1 MHz resolution), so u(t) follows exactly from `t_first_edge` + the symbol table. | yes, to ~1 ms (EL1xxx) or µs (EL1252) | firmware + one input |
| **C. Settled endpoints only** | Compare readings at rest before and after each pass. | no (only where it stops) | analysis only |

**Recommendation:** implement **A** if a counter channel is free (it removes the alignment problem for every future session, including F3), otherwise **B**. F0b_pt2 itself is designed on **settled endpoints (C)**, so it runs correctly even without either fix. If A or B is in place, the same passes also give mid-move lag for free.

**Firmware changes needed for this session (independent of A/B):**
1. **One internal position unit:** 1/256 counts (1 mm = 25 600). Each segment carries its MRES; the firmware converts to steps at that MRES and **asserts divisibility** (abort if a target isn't representable).
2. **MRES switch segment:** at rest only. Write CHOPCONF.MRES over UART, read back CHOPCONF and MSCNT before and after, log both, and **abort if MSCNT changes**.
3. **Settle and record windows timed from `t_last_edge`**, not from segment issue.
4. Telemetry buffered in RAM, printed after `SEQ_DONE` (as now).

---

## 1. Design principles

| Choice | Reason |
|---|---|
| **All turning points representable at every MRES** (multiples of 0.625 µm = one 1/16 step) | the commanded endpoints are identical at every MRES |
| **Test position and saturated turning points on the 40 µm electrical grid** | the coil-current vector at the endpoints is identical at every MRES (with interpolation off, coarser MRES lands on the same table entries); only the path between endpoints differs |
| **Same velocity and acceleration at every MRES** (only the step rate changes) | MRES is the only factor varied in block M |
| **2 s settle + continuous recording** | read each stop at 0.3–0.5 s (F0b ladder timing) **and** 2.0–2.5 s (E2–E5 timing) from the same dwell: settle time tested at no extra cost |
| **ABBA ordering of MRES levels and speeds** | linear drift cancels between levels; the repeat difference bounds drift |
| **One level = E5 replica** (1/16, 2 mm/s, 40 mm/s², 3 mm overrun, 2 s settle) | links directly to β₀ = 1.76 µm |

---

## 2. Configuration

| Item | Setting |
|---|---|
| Driver | TMC2209, SpreadCycle, interpolation off, IRUN = IHOLD (as E5/F0b); **MRES per block** |
| MRES levels | 1/16, 1/32, 1/64, 1/128, 1/256 |
| Test position | **x₀ = 60.000 mm** (1500 electrical cycles, same as F0b) |
| Motion (block M) | **2 mm/s, 40 mm/s²** (E5 settings). Requires 51.2 kHz at 1/256: bench-check first. **Fallback:** 1 mm/s for all levels, and add a "1/16 at 2 mm/s" level to block V |
| Overrun (saturated reversal) | 3 mm (= 75 electrical cycles) |
| Partial strokes | 5 µm and 20 µm (= 8 and 32 steps at 1/16) |
| Settle / record | 2.0 s settle, recording continuous; analysis windows 0.3–0.5 s and 2.0–2.5 s after `t_last_edge` |
| Reference pairs | 15 / 105 mm, both directions, E5 settings (1/16), at start and end |
| Recording | IDS at 1 kHz, scope ≥ 60 min; commanded-position counter (option A) or sync input (option B) if available; BME280/thermistors if fitted |

---

## 3. Blocks

### R1, R3 — Reference pairs (E5 settings)
As E5: 15 / 105 mm, both directions. Also give this session's ISO-style reversal at two positions and the within-session gain (drift at the two ends).

### M — Reversal vs MRES (main block)
**Visit order (ABBA):** 1/16 → 1/32 → 1/64 → 1/128 → 1/256 → 1/256 → 1/128 → 1/64 → 1/32 → 1/16 (10 visits).

**Per visit (≈ 3 min):**
1. **Switch MRES** at rest at x₀ (MSCNT check).
2. **Saturated reversals, ISO-style:** a chain of 7 approaches to x₀, alternating from below and above with 3 mm overruns: below → above → below → … Each arrival: 2 s settle + 0.5 s record. Six consecutive differences = six reversal values (3 up→down, 3 down→up).
3. **Partial strokes:** anchored loops from x₀ (approached from below), x₀ ↔ x₀ + 5 µm, then x₀ ↔ x₀ + 20 µm, 3 cycles each, 2.5 s pause at every turning point.

### V — Speed (settled endpoints), at 1/256
**Speeds:** 0.01, 0.1, 1, 2 mm/s (2 mm/s only if the 51.2 kHz check passes), plus **1/16 at 2 mm/s** if block M used the 1 mm/s fallback.
**Order (ABBA):** 2 → 1 → 0.1 → 0.01 → 0.01 → 0.1 → 1 → 2 mm/s.

**Per speed visit:** 3 mm wipe from below at 1 mm/s → then a chain of 5 approaches to x₀ alternating from below/above with **100 µm overruns** (≫ s₉₀ = 7.5 µm, so saturated), all at the test speed and 40 mm/s² → four reversal values. 2 s settle + 0.5 s record per arrival.

(100 µm instead of 3 mm keeps the 0.01 mm/s visit to ≈ 1.5 min.)

If option A or B is in place: the overrun legs are constant-velocity segments, so the mid-move lag (half the up/down difference at the same commanded position) is extracted from the same passes.

---

## 4. Session

| Block | Content | Time (est.) |
|---|---|---|
| R1 | reference pair (1/16, E5 settings) | 3.5 min |
| **M** | 10 MRES visits × ≈ 3 min | ≈ 30 min |
| **V** | 8 speed visits | ≈ 6 min |
| R3 | reference pair | 3.5 min |
| — | MRES switches, transits | ≈ 2 min |
| **Total** | | **≈ 45 min** |

If over 52 min: drop the 1/32 and 1/128 visits in the second half of M (keeps 16, 64, 256 in ABBA).

---

## 5. Pre-run checks
- [ ] Bench: RMT produces 51.2 kHz cleanly at 1/256 (else use the fallback)
- [ ] MRES switch at rest: MSCNT unchanged before/after, at every level, on the bench
- [ ] 1 mm move at each MRES measures 1.000 mm on the interferometer (± gain), i.e. the unit conversion is right
- [ ] Every turning point divisible at 1/16; x₀ and the saturated turning points on the 40 µm grid
- [ ] Each approach's direction matches its declared side (below/above) in the table (the sign check from F0b/F1)
- [ ] Windows timed from `t_last_edge`; option A/B channel visible in the scope if used
- [ ] Total ≤ 52 min; scope ≥ 60 min

---

## 6. Analysis

### 6.1 Pre-processing
- Drift: R1/R3 plus the x₀ arrivals, interpolated in time (one position, so no position term).
- Each arrival: mean over the 0.3–0.5 s window **and** the 2.0–2.5 s window.
- Reject arrivals whose window std > 50 nm or |slope| > 20 nm/s (report count).
- **Detection limit:** d_min for consecutive arrivals (~10 s apart) from F0b's σ_Δ(T); state it.

### 6.2 Saturated reversal (M, step 2)
- Reversal B = arrival from below − arrival from above (consecutive pairs), per visit; mean ± 2σ over the six values.
- Per MRES: average the two ABBA visits; their difference = drift bound.
- **Model:** B against log₂(MRES), linear fit; slope ± 2σ.
- **Settle effect:** B(0.5 s window) − B(2 s window), per level.
- **E5 replica:** B at 1/16 (2 s window) against β₀ = 1.76 µm.

### 6.3 Partial strokes (M, step 3)
ℓ(5 µm) and ℓ(20 µm), cycles 2–3, settled turning points, up/down averaged (F0b definition), against MRES.

### 6.4 Speed (V)
- B against log₁₀(speed), per speed: mean of the two ABBA visits; drift bound from their difference.
- If option A/B: mid-move lag against speed (Stribeck preview for F3).

### 6.5 Decision

| Result | Conclusion |
|---|---|
| B rises with MRES (beyond 2σ and the drift bound), and 1/16 ≈ 1.76 µm, 1/256 ≈ 2.3 µm | **MRES explains the F0b/E5 difference.** β₀ must be identified at the operating MRES; the friction model needs MRES as a condition (or one MRES is fixed for operation) |
| B independent of MRES, but B(0.5 s) ≠ B(2 s) | **settle time explains it**: the ladder read the stage before it finished relaxing. Use 2 s windows in F1/F0b re-analysis where possible |
| B depends on speed | **speed explains it**: approach speed becomes a model input, or one approach speed is fixed in operation |
| None of these, and B at 1/16 ≈ 2.3 µm | **the stage changed between days** (temperature, lubrication, state): β₀ needs per-session identification, like the gain |
| Partial-stroke ℓ depends on MRES but saturated B doesn't | MRES changes the pre-sliding path (finer steps load the contacts differently) but not the saturation |

---

## 7. Outputs
- `f0b_pt2_command_table.csv` (+ hash), raw recording, firmware log (with MSCNT/CHOPCONF at each switch)
- `f0b_pt2_reversal.csv` (block, MRES, speed, visit, pair, window, B)
- `f0b_pt2_partial.csv` (MRES, stroke, cycle, ℓ)
- `f0b_pt2_summary.json`: B vs MRES (slope ± 2σ), settle effect, B vs speed, E5-replica value, decision

---

## 8. Figures

### Overview (sequence)
Same grammar as F0b/F1: ribbon (R1 · M with sub-ribbon of MRES levels 16 32 64 128 256 256 128 64 32 16 · V with speeds · R3) + commanded position relative to x₀ on the symlog axis. Colour the M visits with a sequential scale by MRES (light = 1/16, dark = 1/256); V visits red with line weight by speed. Caption: segments, minutes.

### Results (two figures + table)

**R-1 — Reversal against microstep resolution** (7.2 × 3.0 in, two panels)
- (a) Saturated reversal B against MRES (x: 16 … 256, log₂ ticks labelled 1/16 … 1/256): two ABBA visits as filled/open markers, mean ± 2σ; filled squares for the 2 s window, open squares for the 0.5 s window, dodged. Horizontal dotted line at β₀ = 1.76 µm ("E2–E5, 1/16") and at 2.30 µm ("F0b ladder, 1/256").
- (b) Partial-stroke lost motion ℓ(5 µm) and ℓ(20 µm) against MRES, same x axis.
- Caption: slope per doubling of MRES ± 2σ against the drift bound; the E5 replica value; the settle effect.

**R-2 — Reversal against approach speed** (3.5 × 2.8 in)
- B against speed (log x, 0.01–2 mm/s), ABBA visits filled/open, ±2σ; the same two reference lines.
- If option A/B: an inset of mid-move lag against speed.

**T1 —** one table: B per MRES (0.5 s / 2 s windows), B per speed, E5 replica, drift bound, d_min, decision row.

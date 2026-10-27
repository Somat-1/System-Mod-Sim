# Experiment E2 — Global Random Map: Sequence & Script Specification

**Purpose:** the main data set for geometric model B. It measures the direction-averaged error at a few hundred positions spread over the whole travel, so that the gain, the screw-periodic harmonics (structure from E1) and the slow cumulative lead trend can be fitted together with one model over the full range.

**Role in the identification chain:**
- **E0 (ISO):** gain (Step A) + final validation
- **E1:** which harmonics exist (H) and stationarity
- **E2:** fit model B over the full travel ← this experiment
- **E3:** same design, new seed, for C's priors and model selection

**Scope of this script:** measurement only. Homing and warm-up run beforehand in the separate warm-up script.
**Precondition at script start:** the stage is homed, warmed up, and parked at x_ref = 61.35 mm, with the interferometer recording.
**Time limit:** the complete recording must fit within **60 min** (target ≤ 50 min).

---

## 1. Fixed system parameters

| Quantity | Value |
|---|---|
| Screw lead L | 2 mm = 200 full steps = **3200 µsteps** per revolution |
| Microstep resolution | 1/16 → 0.625 µm/µstep, 1600 µstep/mm |
| Electrical cycle | 4 full steps = 40 µm = **64 µsteps** (50 per revolution) |

## 2. Settings (identical to ISO230-2_v2 and E1)

| Setting | Value |
|---|---|
| Chopper / currents | SpreadCycle, I_RUN = 400 mA, I_HOLD / TPOWERDOWN / IHOLDDELAY as in ISO v2 |
| Motion | 2 mm/s, 40 mm/s², one trapezoidal RMT pulse train per move |
| t_settle / t_record | as validated in the E1 pilot (≥ t_hold + 0.5 s / 0.5 s) |
| Pre-approach / overrun | 3 mm |
| x_ref | 61.35 mm, fixed-side approach from below (pre-move to 58.35 mm) |
| Datum convention | global mean d, as in ISO v2 |
| Span | **2.000–120.000 mm**, the same span as the ISO v2 baseline, so G covers exactly the range E0 validates |

---

## 3. Target design

### 3.1 Requirements and how the design meets them

| Requirement | Why | Design choice |
|---|---|---|
| Even coverage of the whole travel | the spline (slow lead trend) needs points everywhere; no gaps | **stratified random:** the span is split into N equal bins with one target per bin |
| All screw phases sampled evenly | the harmonics are fitted by folding on screw angle | random position within each bin, plus a **phase-balance check** (§3.3) |
| Motor error must not disturb the fit | microstep error (40 µm period) | **all targets on the 64-µstep grid** → identical electrical phase → microstep error = one constant, absorbed by the global offset a₀ (same principle as E1) |
| Independent of the validation set | ISO 230-2 §5.2: acceptance targets shall differ from compensation sampling points | exclude any target coinciding with an E0 target (±1 grid point) |
| Reproducible | E3 and the thesis | fixed seed, full target list logged |

**Consequence of the 64-µstep grid:** the targets can only fall on 50 screw phases (3200 / 64). That's the same 50 phases E1 uses, so E1 and E2 are directly comparable. As in E1, E2 resolves only the periodic and trend content of the error, not the absolute microstep offset.

### 3.2 Number of targets N

N is chosen from the time per stop measured in the E1 pilot (see the render, panel d):

$$N = \Big\lfloor \frac{(50 - 4.5)\ \text{min} \cdot 60}{2\, t_{stop}} \Big\rfloor \quad\text{(capped at 500)}$$

| t_stop | N per direction | Recording time |
|---|---|---|
| 1.6 s | 400 (default) | ≈ 25 min |
| 2.4 s | 400 | ≈ 37 min |
| 3.0 s | 400 | ≈ 45 min |
| 3.5 s | 380 | ≈ 49 min |

**Statistical reach with N = 400** (σ_ȳ ≈ 0.35 µm after up/down averaging):

| Quantity | Approximate precision |
|---|---|
| Harmonic amplitude | ≈ σ_ȳ · √(2/N) ≈ **0.025 µm** |
| Gain a₁ | ≈ σ_ȳ / (√N · std(x)) ≈ **0.5 ppm**, an independent cross-check of Step A (−25.09 ppm from E0) |
| Model size | ~ 2 + 2H + K ≈ 15–25 parameters, so ≥ 16 points per parameter |

More points add little precision. Coverage and independence matter more, which is why E3 (a new seed on a different day) is a separate experiment rather than a longer E2.

### 3.3 Generator

```python
# ── E2 target generator: stratified random on the 64-µstep grid, phase-balanced ──
import numpy as np

U, G, REV = 1600, 64, 3200              # µsteps/mm, electrical cycle, µsteps per revolution
LO_U, HI_U = 3200, 192000               # 2.000 mm, 120.000 mm (both multiples of 64)
N = 400
ISO_TARGETS_U = [...]                   # E0 targets in µsteps (exclusion list)

def make(seed):
    r = np.random.default_rng(seed)
    edges = np.linspace(LO_U, HI_U, N + 1)
    t = []
    for k in range(N):
        g0 = int(np.ceil(edges[k] / G))
        g1 = int(np.ceil(edges[k + 1] / G)) - (0 if k == N - 1 else 1)
        t.append(G * r.integers(g0, g1 + 1))
    return np.array(t)

def valid(t):
    phase_counts = np.bincount((t % REV) // G, minlength=50)
    clash = any(abs(t - iso).min() <= G for iso in ISO_TARGETS_U)
    return len(np.unique(t)) == N and phase_counts.min() >= 5 and not clash

seed = 5200
while not valid(t := make(seed)):
    seed += 1                           # log the seed actually used
targets_u = np.sort(t)                  # integer µsteps, all ≡ 0 (mod 64)
```

**Generator checks:**
- N unique targets, all multiples of 64 µsteps, sorted;
- every screw-phase bin has **≥ 5 targets** (mean = N/50 = 8);
- no target within ±1 grid point of an E0 target;
- consecutive spacing 0.04–0.6 mm (one grid point minimum);
- all positions including the 3 mm pre-approach/overrun inside the soft limits;
- the seed and the target list are saved with the session.

---

## 4. Sequence

```
[precondition: homed, warmed up, parked at x_ref, interferometer recording]

INIT → DATUM (3 × x_ref)
     → pre-approach to 2 mm − 3 mm
     → UP PASS:   N targets ascending, stop–settle–record at each
     → overrun to 120 mm + 3 mm
     → TURNAROUND REF (1 × x_ref, fixed side) → back to 120 mm + 3 mm
     → DOWN PASS: the same N targets descending, stop–settle–record at each
     → END_REF (3 × x_ref) → PARK → DONE
any fault → ABORT → PARK
```

| State | Actions |
|---|---|
| **INIT** | load config + target list; write and read back registers; check the stage is at x_ref; check the interferometer signal; write the session header |
| **DATUM** | 3 × x_ref (pre-move to 58.35 mm → 61.35 mm), settle, record |
| **UP PASS** | pre-approach to −1.0 mm, then for each target in ascending order: one trapezoidal move (≈ 0.1–0.4 s depending on spacing), settle, record (mean, std, slope); telemetry only after the record window |
| **TURNAROUND REF** | from the overrun at 123 mm: transit to 58.35 mm → x_ref, settle, record; transit back to 123 mm |
| **DOWN PASS** | for each target in descending order: move, settle, record |
| **END_REF** | 3 × x_ref, environment snapshot |
| **PARK** | stay at x_ref; stop recording; write the summary |

**Why the targets are visited as two sorted passes rather than in random order:**
- Each pass only moves one way, so every target is reached from a defined direction without a pre-approach per stop. Random order would need a 3 mm pre-approach before every target, roughly doubling the time.
- **Built-in drift cancellation:** the down pass mirrors the up pass in time. A target visited at time t on the way up is visited at about T − t on the way down, so its average visit time is ≈ T/2 for *every* target. A linear thermal drift therefore adds the same constant to every ȳᵢ, which is absorbed by a₀ and doesn't bias the shape of G. (In the reversal Bᵢ, drift appears as a linear trend in x. That's why the drift correction from the reference visits is still applied.)
- The turnaround reference sits at the symmetry point, so it doesn't break the symmetry.

**Move sizes:** consecutive targets are 0.04–0.6 mm apart (mean 0.295 mm). Moves below 0.1 mm are triangular profiles below 2 mm/s; longer moves reach 2 mm/s. All stay below the 2.5 mm/s encoder limit.

---

## 5. Time budget

| Block | Time |
|---|---|
| INIT + DATUM | ~0.5 min |
| Pre-approach + transits (incl. turnaround ref) | ~2.5 min |
| Up + down pass | 2N × t_stop |
| END_REF | ~0.5 min |
| **Total (N = 400)** | **≈ 25 min at 1.6 s · ≈ 37 min at 2.4 s · ≈ 45 min at 3.0 s** |

**Time watchdog:** if the projection at the turnaround exceeds 58 min, abort the down pass and flag it. An up-only data set can't be direction-averaged and would have to be repeated.

---

## 6. Checks (automatic, no operator pause)

The settle and record settings were validated in the E1 pilot, so E2 has no operator pause. The script checks automatically:
- per stop: record-window slope × t_record < ½ s and no saturated encoder samples. Flag the stop; don't repeat it inline (that would break the time symmetry);
- at every reference visit: step-count check (|measured − commanded − datum| < tolerance), otherwise ABORT;
- x_ref drift between consecutive visits > 0.5 µm → warning.

## 7. Logging

Same format as ISO and E1 (`stops.csv`, `events.csv`, `env.csv`, `summary.json`), with fields `state` (DATUM / UP / TURN_REF / DOWN / END_REF), `target_index`, `target_usteps`, `screw_phase_bin = (target_usteps mod 3200) / 64`, and MSCNT after the record window if available.

## 8. Acceptance criteria

- [ ] ≥ 95 % of stops pass the record-window check, both directions for ≥ 95 % of targets
- [ ] no watchdog abort; reference visits consistent within the drift tolerance
- [ ] recording ≤ 60 min
- [ ] ΔT_air < 0.5 K over the recording

---

## 9. Analysis (model B), fixed before the run

1. **Pre-processing:** deviation per stop as in ISO v2; drift correction from the 7 reference visits; per target ȳᵢ = ½(e↑ + e↓) and Bᵢ = e↑ − e↓. Targets missing one direction are excluded from ȳ.
2. **Model B:**

   $$\bar y(x) = a_0 + a_1 x + \sum_{h=1}^{H}\big[c_h\cos(2\pi h\,\theta) + s_h\sin(2\pi h\,\theta)\big] + \sum_k b_k B_k(x), \qquad \theta = \frac{\text{commanded µsteps mod } 3200}{3200}$$

   - **θ** is the absolute screw angle (from the commanded step count), the same phase reference as E1;
   - **H** from E1 (candidates H_E1 − 1 … H_E1 + 1);
   - **B-spline** basis with its linear part removed (so a₁ stays the pure gain), knot spacing ≥ 10 mm, K ∈ {0 … 8} interior knots;
   - fit by weighted least squares.
3. **Model selection:** **blocked** 10-fold cross-validation (10 contiguous blocks of ~12 mm) for H and K. Blocked rather than random folds, so the CV tests interpolation over real gaps, not just between neighbouring points.
4. **Residual checks:** the residual periodogram should be flat (no leftover peaks at k/L); residual rms ≈ σ_ȳ.

## 10. Cross-checks against independent data

| Check | Data | Expectation |
|---|---|---|
| Gain | a₁ from E2 vs Step A (−25.09 ppm, E0) | agree within the combined uncertainty |
| Harmonics | E2's fitted c_h, s_h vs the E1 window values (same θ reference) | amplitudes and phases agree; E2 averages over the travel, E1 is local |
| Prediction | model B evaluated at the E0 targets vs E0's x̄ᵢ | residual rms ≈ E0's repeatability level; this is the first independent validation of G |

## 11. Deliverables

- [ ] ȳ(x) over the full travel with the fitted B, and its components (gain, periodic, slow trend) plotted separately
- [ ] coefficient table: a₁ (ppm), c_h / s_h → amplitude and phase per harmonic, spline coefficients; all with CIs
- [ ] CV table (H, K) and the chosen structure
- [ ] residual vs x, and the residual periodogram
- [ ] B(x) over the full travel (handed to the friction work)
- [ ] cross-check table (§10)
- [ ] validity table (§8)

# E3 — repeat global map (30 Oct 2026)

> **Not the new-seed E3.** This session measured E2's own 400 targets; the
> session the procedure actually specifies was generated afterwards and lives
> in [`../e3_global_map/`](../e3_global_map/). Section 1 below explains how
> that happened and what it costs.

The session recorded on 30 Oct 2026, 13:56–14:42, processed against E2 of 27 Oct
under `instructions/E3_processing_instructions.md`.

Two things about it have to be said before any number is quoted, because both
change what the session can be used for.

## 1. It is a repeat at the same targets, not a new-seed session

The procedure specifies E3 as the E2 design with a **different seed**, so that
the map is tested at **new positions**. What ran was the unmodified
`e2_global_map` firmware, `config=96579c689d4d`. It derives its target list from
`seed_start: 5200` in `e2_config.json`, which was not changed, so it re-derived
the same accepted seed `5315` and emitted the same 400 targets. The 400
`# TARGET` lines of the two terminal logs diff clean; see
`sequence/PROVENANCE.md`.

| Question | Status |
|---|---|
| Q1, does the map predict a **new day**? | answered |
| Q1, does the map predict **new positions**? | **not answerable from this session** |
| Q2, are C_k(x), s(x) and the gain the same on both days? | answered for C_k(x); the gain and s(x) are limited by item 2 below |
| Q3, is the residual structure geometric or temporal? | answered, and more sharply than planned |
| Q4, is the screw-angle reference stable? | answered, with an independent cross-check |
| Q5, is friction leaking into the geometry? | answered, with the caveat in item 2 |
| Q6, final pooled model | identified, but it gains replication rather than new positional support |

Where the design loses Q1 it gains elsewhere. Two sessions at identical targets
pair one to one, so the between-day residual comparison needs no running mean,
no common grid and no circular-shift null — the smoothing that Q3's prescribed
test needs only because two seeds land in different places. Both tests are
reported: the paired one because it is the real answer, the smoothed one because
a future new-seed session will have to use it.

What the pooled fit does **not** get is new positional support. Leave-one-session-out
CV on these two sessions is a replication test, not a generalisation test, and
is reported as such.

## 2. The reference moved 5.4 µm during the session

The interferometer datum at `x_ref` reads, over the three reference groups:

| Visit group | t [s] | error [µm] |
|---|---|---|
| DATUM (3 visits) | 324–335 | +0.68, +0.79, +0.75 |
| TURN_REF | 1482 | +1.11 |
| END_REF (3 visits) | 2625–2642 | −4.61, −4.62, −4.78 |

against a ±0.5 µm validity criterion, so **this check fails**. Two features of
it matter more than the size:

- **It is not noise.** The three visits within each group agree to 0.17 µm. The
  datum genuinely moved.
- **It is not a lost step.** MSCNT is identical at all seven references, and
  every reference uses the same absolute command, so the driver did not lose
  steps. Whatever moved was mechanical, thermal or optical.
- **It is not linear.** The first half drifts +0.37 µm; the second half drifts
  −5.78 µm, a rate of −18 µm/h.

### Why the geometry survives it and the reversal does not

The mirrored pass design makes the **mean visit time of a target constant along
the travel**: going up, target *i* is reached after the gaps below it; coming
down, after the gaps above it, and the two sum to the same total for every *i*.
Measured spread of that mean time across the 400 targets: **0.013 s**. So a
drift that is *linear in time* adds the same constant to every ȳ and vanishes
when the mean is removed — which is why the D0 and D1 variants return ȳ shapes
agreeing to ~10⁻⁵ µm. Only the drift's **curvature** reaches ȳ.

B = e↑ − e↓ is the opposite case. It differences the two directions, so it
carries the whole drift accumulated over the up-to-down gap, and that gap runs
from 68 s at the turnaround to 2225 s at the ends. β₀ comes out at −1.16 µm
under D0 and −3.85 µm under D2.

Consequences, all carried in `analysis/e3_findings_summary.json` under
`drift_variant_robustness`:

| Quantity | Spread over D0/D1/D2 | Usable? |
|---|---|---|
| A₁, A₂ along x | 2·10⁻⁴ µm | yes, unaffected |
| s(x) range | 4·10⁻⁵ µm | yes, unaffected |
| gain | 24.5 ppm | **D2 only**, and see below |
| β₀ | 2.69 µm | **no** — E3 contributes no reversal estimate |

There is a clean reason the harmonics and s(x) come through untouched while the
gain does not. With a piecewise-linear drift and two mirrored monotone passes,
the correction to ȳ is **linear in x exactly**: on the up leg the time runs
with x, on the down leg against it, so the two segment slopes combine into one
slope and the curvature cancels. Measured: the D2 correction is +24.5 ppm of
gain and 0.003 µm of everything else. The drift this session suffered reaches
the gain and reaches nothing else.

That makes the gain comparison uninterpretable and leaves the rest standing.
The E3−E2 gain difference is +3.05 ppm, which is **12 % of the correction
applied** — a one-in-eight error in the drift model accounts for it entirely,
so it is not evidence about the screw. By the same token the 1.76 µm s(x)
difference is *not* a drift-correction artefact; it is either a real day-to-day
change in the lead error or drift curvature *within* a segment, and three
reference visits cannot separate those.

The spread is **not an error bar.** D0 applies no correction and D1 provably
removes only a constant, so with a drift this large neither is a candidate
answer: the D0-to-D2 distance measures how much contamination D2 took out, not
how uncertain the result is. What is genuinely unbounded is whether a straight
line between three references describes the real drift. Quantities whose spread
is negligible are safe on either reading, which is the useful half of section
3's rule.

**Two practical consequences.** The reversal applied in the validation comes
from E2, E4 and E0, with E3 excluded and its variant spread recorded. And the
split-direction leakage test of section 8.3 is read on **E2 alone**: s↑(x) −
s↓(x) is the smooth part of B(x), so a pooled test including this session
measures its drift and calls it friction. The three versions differ by 5×
(max |s↑−s↓| = 0.66 µm on E2, 3.37 µm on E3, 1.71 µm pooled), which is the
whole reason to run them separately.

### One excluded record window

One UP window of 807, at target 5.24 mm, fails the std (0.65 µm) and slope
(2.35 µm across the window) gates — the hold was short enough that the
averaging window overlapped settling. It is the session's largest single
residual at −1.78 µm. It is excluded, logged in
`analysis/e3_excluded_windows.csv`, and the same mask is applied to E2, where it
removes nothing, so both sessions still go through one code path. Model residual
rms without it: 0.415 → 0.402 µm.

## Headline results

Full tables and the section 13 readings are in [FINDINGS.md](FINDINGS.md).

- **The map transfers between days, at a cost.** E2's model applied to E3 with
  no refit and one offset: 1.707 µm raw → 0.689 µm, against 0.402 µm for E3's
  own fit and a 0.24 µm noise floor. 78 % of the gap closed. The reverse
  direction gives 76 %, so the difference is symmetric between the sessions
  rather than one session being at fault.
- **The harmonic profile is a stage property.** max |ΔC₁| = 0.217 µm,
  max |ΔC₂| = 0.230 µm, rms ≈ 0.09 µm for both, and the drift variants move
  the amplitudes by 2·10⁻⁴ µm. This is the part of the model that can be
  frozen with confidence.
- **The residual is about half geometric.** Paired over 399 identical targets,
  ρ = +0.479, which splits the ≈0.39 µm residual into **0.265 µm repeatable**
  and **0.279 µm independent** — the latter sitting on the 0.24 µm E4 local
  floor. The repeatable half is a candidate for a later refinement; the
  independent half is the accuracy floor of the geometric layer.
- **Homing moved by about two degrees, and three independent estimates agree.**
  δ₁ = −1.99 ± 0.36°, δ₂ = −2.62 ± 0.58° from the harmonic phases, and the
  driver's own MSCNT shifted 240 counts = 15 µsteps = 1.69° = 9.4 µm of travel.
  Two of those estimates know nothing of the other. Inside the 3° requirement,
  but homing should be fixed to a known µstep phase.
- **β₀ is translation-invariant across three sessions**: E2 −1.76, E4 −1.66,
  E0 −1.60 µm, range 0.156 µm. E3 contributes nothing (see item 2 above).
- **Pooled model**: H = 4, K_s = 2, K_C = 8, λ = 3.2; 59 parameters, 38.1
  effective; residual rms 0.478 µm; gain −28.25 ± 0.52 ppm. Empirical-Bayes
  hyperparameters τ² = 30.1, γ = 2.65 are frozen for the friction iteration.
- **Predicted compensated ISO 230-2 performance**, after G + B: E 8.52 → 2.30,
  M 6.88 → 1.96, A 9.25 → 3.50 µm on the v4 session (v3: 7.85 → 2.12,
  6.46 → 1.55, 9.15 → 3.51).

## Layout

| Path | Contents |
|---|---|
| `e3_run.py` | capture → holds → planned stops. All substantive logic imported from `e2_run`, so both sessions are reduced by one implementation |
| `process_e3.py` | sections 2–9: validity, drift, cross-prediction, coefficient and homing comparison, residual structure, pooled and empirical-Bayes identification, reversal and leakage, V2/V3 validation |
| `render_e3_findings.py` | figures F1–F6 and the appendix |
| `write_findings.py` | generates `FINDINGS.md` from the summary JSON, so no number in the document is hand-copied |
| `FINDINGS.md` | **the results**: tables T1–T5, the section 13 readings, and what to do next |
| `data/` | the capture, the ESP/TMC log, and the alignment cache (gitignored) |
| `sequence/PROVENANCE.md` | why there are no sequence files here |
| `analysis/` | tables T1–T5, the per-target CSV, the summary JSON and the fitted models |
| `rendered/` | `experimental_findings/` and `presentation/` |

The log is also committed at `../e2_global_map/data/e2_terminal_20261030_135726.log`,
where it first landed; the copy here is the one the scripts read.

## Reproducing

```sh
python e3_run.py          # rebuild the alignment cache from the CSV   (~25 s)
python process_e3.py      # sections 2-9; writes analysis/             (~8 min)
python write_findings.py  # FINDINGS.md from the summary JSON
python render_e3_findings.py
```

Almost all of `process_e3.py` is the pooled model selection: 270 candidate
models, each scored over five shifted 10-fold partitions plus
leave-one-session-out. Nothing else in the script takes meaningful time.

`process_e3.py` reads E2 through `e2_run` and `render_e2_findings`, so E2's
published numbers are recomputed rather than quoted: its M_C residual rms comes
back as 0.3685 µm against the 0.3685 µm in `e2_findings_summary.json`.

## Conventions

As E1/E2/E4. Deviation positive means measured travel exceeds command;
θ = 2π(commanded µsteps mod 3200)/3200; e = A cos(kθ + φ); pattern shift
Δθ_k = −(φ_k − φ̄_k)/k. The basis is E2's realised span, 2.280–119.920 mm, with
its reference position x̄ = 61.0022 mm, used unchanged for E3 and for the pooled
fit.

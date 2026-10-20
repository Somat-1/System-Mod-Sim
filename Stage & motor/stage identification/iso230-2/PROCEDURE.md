# ISO 230-2 Baseline Measurement — Sequence & Script Specification

**Purpose:** a full-range baseline of positioning accuracy and repeatability (B-sequence Step 0), evaluated in ISO 230-2 terms.
**This document does two things:**
1. it documents what the automated sequence does, for the thesis and the lab log, and
2. it is the specification for writing the script.

**Axis:** stepper-driven ball-screw stage, range of motion (ROM) ≈ 140 mm.
**Sensor:** IDS 3010 interferometer.
**Driver:** TMC2209, SpreadCycle, MRES as selected in the microstep test.
**Homing:** physical end switches. One is the home switch; both are hard limits.

---

## 1. Sequence overview

One script runs the whole session as a state machine:

```
INIT ──► HOME ──► WARMUP ──(operator stop or t_max)──► DATUM ──► ISO_CYCLE ──► END_REF ──► PARK ──► DONE
  │        │         │                                   │           │             │
  └────────┴─────────┴──────────── any fault ────────────┴───────────┴─────────────┴──► ABORT ──► PARK
```

| State | What it does | Exit |
|---|---|---|
| INIT | check hardware, write and verify driver registers, load config and targets, open log files | all checks pass |
| HOME | find the home switch, define commanded zero, park | home found |
| WARMUP | repeat the ISO move/dwell pattern **without recording targets**; visit x_ref periodically and display drift | operator presses STOP **or** t_max reached, then the current cycle is finished |
| DATUM | approach x_ref from the fixed side and record the interferometer/command pair as the measurement datum | datum recorded |
| ISO_CYCLE | n bidirectional runs over the target list, dwelling and recording at each target | all runs done |
| END_REF | final x_ref visits and environment snapshot | done |
| PARK | move to the park position and set the end motor state | done |
| ABORT | stop motion immediately, log the reason, then go to PARK if safe | — |

---

## 2. Configuration parameters

All parameters live in one config file (YAML or JSON), which is copied into the output folder at the start of every session.

| Key | Meaning | Default / value |
|---|---|---|
| `lead_mm` | screw lead | ___ |
| `mres` | microstep resolution | from microstep test |
| `i_run`, `i_hold` | driver currents (register values) | operating values |
| `tpowerdown`, `iholddelay` | hold-current transition registers | operating values |
| `t_hold_s` | resulting time to reach I_HOLD | ___ (derived) |
| `v_mm_s`, `a_mm_s2` | approach velocity and acceleration | e.g. 5 / ___ |
| `home_side` | which switch is home | `min` |
| `home_seek_v`, `home_slow_v` | fast and slow homing velocities | e.g. 5 / 0.5 mm/s |
| `home_backoff_mm` | back-off before the slow approach, and after homing | 2.0 |
| `soft_min_mm`, `soft_max_mm` | software travel limits | home + 2 / ROM − 2 |
| `x_start_mm`, `x_end_mm` | first and last target | 5.0 / 135.0 |
| `p_nom_mm` | nominal target interval | 5.0 |
| `seed` | random seed for target offsets | logged |
| `overrun_mm` | overrun past the first/last target | 3.0 |
| `x_ref_mm` | drift reference position | 70.0 |
| `ref_preapproach_mm` | fixed-side pre-move before every x_ref visit | 3.0 (from below) |
| `t_settle_s` | wait after arrival before recording | `t_hold_s + 0.5` |
| `t_record_s` | recording window | 0.5 |
| `n_runs` | ISO runs | 5 (8 if time allows) |
| `warmup_ref_every_cycles` | x_ref visit cadence during warm-up | 1 (after every warm-up cycle) |
| `warmup_t_max_min` | warm-up time cap | 30 |
| `warmup_drift_hint_um_min` | drift-rate threshold shown as "ready" to the operator | 0.1 |
| `ifm_rate_hz` | interferometer streaming rate | ≥ 1000 |
| `watchdog_tol_um` | allowed commanded-vs-measured jump at x_ref | ≈ 5 × microstep size + expected error span |
| `ifm_min_signal` | minimum interferometer signal quality | from IDS settings |
| `end_motor_state` | after PARK | `hold` or `disabled` (keep consistent) |

### Key numbers (140 mm ROM)

| Item | Value | Why |
|---|---|---|
| End safety margin | **2 mm** from each mechanical limit | protects the switches, end stops and nut |
| Overrun | **3 mm** | ≫ reversal lost motion (µm scale); leaves room to decelerate, reverse, and reach constant velocity before the target |
| Measured span | **5 → 135 mm** | ROM − 2 × (margin + overrun) |
| Target interval | **5 mm nominal → 27 targets**, offsets uniform ±1.25 mm | randomizes the phase relative to the screw lead |
| Dwell | **≈ 2 s** = t_settle + t_record | measure after the hold-current transition |
| Duration | **≈ 18 min** for 5 runs | plus homing (~1 min), warm-up (≤ 30 min), end (~2 min) |

---

## 3. State specifications

### 3.1 INIT
1. Create the output folder `YYYYMMDD_HHMM_iso230-2/` and copy the config file into it.
2. Connect: driver (UART), step generator, interferometer, environment sensors, end switches.
3. Write the driver registers (SpreadCycle, MRES, currents, TPOWERDOWN, IHOLDDELAY, intpol), **read them back and verify**. Save the register dump.
4. Check the interferometer is streaming, its signal is ≥ `ifm_min_signal`, and environmental compensation is in the expected state.
5. Load the target list (§4) or generate it from `seed`. Convert each target to an integer microstep count; store both the count and the corresponding commanded position.
6. Check that every planned position (targets ± overrun, x_ref ± pre-approach) lies inside the soft limits.
7. Start the continuous loggers (§6).

**Fail →** ABORT, with the motor never enabled.

### 3.2 HOME
1. Enable the motor at I_RUN.
2. Move towards the home switch at `home_seek_v` until it triggers, then stop.
3. Back off `home_backoff_mm`.
4. Re-approach at `home_slow_v` until it triggers, then stop. This is **commanded zero**.
5. Back off `home_backoff_mm` and move to x_ref via the fixed-side approach.
6. Log the switch trigger positions (fast and slow) in both commanded steps and interferometer reading. Over several sessions this gives the homing repeatability.

Home **always on the same side with the same sequence**.

**Fail** (switch not found within ROM + margin, or the opposite limit triggers) **→** ABORT.

### 3.3 WARMUP
**Pattern:** execute the ISO cycle of §5 exactly (same targets, overruns, velocities and dwells) but **without marking any targets as measurements**. The thermal load during warm-up then matches the thermal load during measurement. Continuous running would heat the motor more than the dwell-heavy ISO cycle and cause drift during the measurement.

**Loop:**
```
cycle = 0
while True:
    run one warm-up cycle (= one ISO run, unrecorded)
    cycle += 1
    if cycle % warmup_ref_every_cycles == 0:
        visit x_ref (fixed-side approach), dwell, record   # logged as WARMUP_REF
        update live display: x_ref vs time, drift rate over the last 3 visits
        show "READY" hint when |drift rate| < warmup_drift_hint_um_min
    if STOP_requested or elapsed >= warmup_t_max_min:
        break                                              # only at the end of a cycle
```

**Operator stop:** the operator ends the warm-up when it looks sufficient (a key press, GUI button, or a flag file, implementation's choice).
- The request is **latched**. The script finishes the current warm-up cycle and its x_ref visit, then transitions. It never stops mid-cycle, so the friction history and thermal state at the start of DATUM are the same in every session.
- **t_max** forces the same transition if no stop arrives.
- Log the stop reason (`operator` / `t_max`), the warm-up duration, and the last drift rate.

The x_ref log from the warm-up is kept: it doubles as a **thermal-settling curve** of the stage.

### 3.4 DATUM
1. Pre-move to x_ref − `ref_preapproach_mm`, then approach x_ref.
2. Dwell `t_settle_s`, then record for `t_record_s`.
3. Store the datum pair (commanded x_ref, mean interferometer reading). All ISO deviations are computed relative to this pair.
4. Repeat 3× and store all three (scatter = datum quality); use the mean as the datum.

**No re-homing after the warm-up.** The interferometer defines the datum; re-homing would add a switch error and an extra move for no benefit.

### 3.5 ISO_CYCLE
Run the cycle of §5 for `n_runs` runs, with the dwell and record procedure of §5.2 at every target.
- One x_ref visit **at the start of every run** (drift tracking).
- Live display: the deviation at each target as it is recorded, plus the x_ref drift.
- Watchdog checks (§7) run continuously.

### 3.6 END_REF
1. 3 × x_ref visits (fixed-side approach), dwell and record.
2. Environment snapshot.

### 3.7 PARK
1. Move to x_ref (or a defined park position).
2. Set `end_motor_state`.
3. Close the loggers, write the session summary (§6.3), and flush all files.

### 3.8 ABORT
1. Stop motion immediately (controlled deceleration if possible; hard stop if a limit switch triggered).
2. Log the state, reason, timestamp and the last valid data index.
3. If the abort wasn't caused by a limit switch or lost position, go to PARK; otherwise, disable the motor and stay put.
4. Mark the session folder `ABORTED` in the summary file.

---

## 4. Target positions

Generate them once, save the list, and reuse **the same list** for the "after compensation" run at the end of the thesis.

```python
# ── Target generator: ISO 230-2 targets with randomized spacing ──────────
# 27 targets over 5..135 mm, nominal 5 mm apart, each nudged by a random
# offset of up to ±p/4 so the targets don't all hit the same phase of the
# screw lead (a periodic error would otherwise be sampled at one phase only).
import numpy as np

X_START, X_END, P_NOM = 5.0, 135.0, 5.0      # mm
SEED = 2302                                  # log this; reuse it for re-runs
LEAD = None                                  # set the screw lead [mm] to run the phase check

rng = np.random.default_rng(SEED)
nominal = np.arange(X_START, X_END + 1e-9, P_NOM)
targets = nominal + rng.uniform(-P_NOM / 4, P_NOM / 4, nominal.size)
targets[0], targets[-1] = X_START, X_END     # pin the end targets to the span edges
targets = np.round(np.sort(targets), 3)

if LEAD:                                     # phase coverage relative to the screw lead
    phase = np.mod(targets, LEAD) / LEAD
    print("phase coverage (should spread over 0..1):", np.round(np.sort(phase), 2))

np.savetxt(f"iso_targets_seed{SEED}.csv", targets, fmt="%.3f", header="target_mm")
```

**Lead check:** if the nominal interval were an integer or half-integer multiple of the lead, the targets *without* offsets would sample only one or two phases of the periodic error. With L = 2 mm, for example, 5 mm = 2.5 L gives only 2 phases. The random offsets fix this; confirm with the phase-coverage printout.

**Rounding to microsteps:** the commanded position used in the evaluation is the **rounded** microstep count × microstep size, not the unrounded target.

---

## 5. ISO 230-2 standard linear cycle

### 5.1 One run
```
visit x_ref (fixed-side approach), dwell, record           # drift reference
move to P1 − overrun                                        # 3 mm below P1
for i = 1 … m:   move to Pi (↑), dwell, record
move to Pm + overrun                                        # 3 mm above Pm
for i = m … 1:   move to Pi (↓), dwell, record
move to P1 − overrun
```
- Overruns use the same velocity. No recording during overruns, but the interferometer keeps streaming.
- **No pauses between runs** apart from the x_ref visit, so the thermal state stays continuous.

### 5.2 Dwell and record at each stop
```
arrive ── settle wait t_settle = t_hold + 0.5 s ── record window t_record = 0.5 s ── move on
```
- **t_hold** = TPOWERDOWN + IHOLDDELAY ramp to I_HOLD. Recording after this transition measures the real holding state, which matters because the hold current changes stiffness and so the friction offset.
- The **+0.5 s margin** lets structural ringing (lowest mode ≈ 177 Hz) and early creep die out.
- **Settled value** = mean over the record window. Also store the std and the slope (a linear fit over the window) for the settling check (§8).
- Keep the dwell **identical for every stop** in warm-up and measurement.

---

## 6. Outputs

### 6.1 Continuous streams
| File | Content |
|---|---|
| `ifm_stream.*` | t, x_IFM, signal quality (≥ 1 kHz; binary or HDF5 recommended) |
| `env.csv` | t, T_air, p, RH, T_stage (≥ 1 Hz) |

### 6.2 Event log (`events.csv`)
One row per event: `t, state, run_j, direction, target_i, cmd_steps, cmd_pos_mm, event`, where `event` ∈ {move_start, arrive, settle_end, record_start, record_end, overrun, xref, warmup_ref, datum, stop_request, state_change, abort}.

### 6.3 Per-stop results (`stops.csv`) and summary (`summary.json`)
- `stops.csv`: `state, run_j, direction, target_i, cmd_pos_mm, x_mean, x_std, x_slope, t_record_start`. Covers ISO stops, x_ref visits, warm-up refs and datum.
- `summary.json`: config hash, register dump, seed, warm-up duration and stop reason, datum, number of valid runs, abort info if any, environment min/max.

---

## 7. Watchdogs (active in all motion states)

| Check | Condition | Action |
|---|---|---|
| Limit switch | either switch triggers outside HOME | ABORT (hard stop) |
| Soft limits | a commanded position is outside the soft limits | refuse the move → ABORT |
| Interferometer signal | quality < `ifm_min_signal` or stream gap | ABORT (incremental count may be lost) |
| Lost steps | at x_ref: \|measured − commanded − datum offset\| > `watchdog_tol_um` | ABORT and flag |
| Driver status | TMC2209 error flags (overtemperature, short, open load) via UART poll | ABORT |
| Environment | ΔT_air since DATUM > 0.5 K | warn and log (no abort) |

---

## 8. Post-processing and evaluation (ISO 230-2:2014)

**Validity checks first:**
- Settling: per stop, slope × t_record < ½ of the expected repeatability. ≥ 95 % of stops must pass.
- Drift: x_ref over the session < ½ R, or apply a linear drift correction between x_ref visits **and report it**.
- No watchdog events during ISO_CYCLE.
- Target phase coverage relative to the lead confirmed.

**Per stop:** x_ij↑/↓ = (x_mean − datum_IFM) − (cmd_pos − datum_cmd).

| Quantity | Formula |
|---|---|
| Mean unidirectional deviation | x̄_i↑ = (1/n) Σ_j x_ij↑ (likewise ↓) |
| Mean bidirectional deviation | x̄_i = (x̄_i↑ + x̄_i↓) / 2 |
| Reversal value | B_i = x̄_i↑ − x̄_i↓; B = max \|B_i\|; B̄ = mean B_i |
| Std. estimator | s_i↑ = sqrt( Σ_j (x_ij↑ − x̄_i↑)² / (n − 1) ) |
| Unidirectional repeatability | R_i↑ = 4 s_i↑; R↑ = max R_i↑ |
| Bidirectional repeatability | R_i = max(2s_i↑ + 2s_i↓ + \|B_i\|, R_i↑, R_i↓); R = max R_i |
| Systematic deviation | E↑ = max x̄_i↑ − min x̄_i↑; E = max(x̄↑, x̄↓) − min(x̄↑, x̄↓) |
| Mean bidirectional range | M = max x̄_i − min x̄_i |
| Accuracy | A = max(x̄_i↑ + 2s_i↑, x̄_i↓ + 2s_i↓) − min(x̄_i↑ − 2s_i↑, x̄_i↓ − 2s_i↓) |

Check these definitions against your copy of the standard before reporting.

**Standard plots:**
1. x̄_i↑ and x̄_i↓ vs position, with ±2s bands.
2. B_i vs position.
3. s_i↑ and s_i↓ vs position.
4. x_ref vs time over the whole session, with the warm-up and measurement phases marked.

**Hand-off to the B sequence:** x̄_i → Step 2 (geometric GP); B_i → Step 3 (backlash / friction split); s_i → repeatability floor for Step 5; A, R, E, M → baseline row of the final results table.

---

## 9. Pre-session checklist (manual, before starting the script)

- [ ] Interferometer head on the stage base, reflector on the carriage at the functional point. Abbe offset measured and entered in the config.
- [ ] Beam aligned: signal strength checked at **both** ends of the ROM.
- [ ] Dead path minimized.
- [ ] Environmental compensation configured; environment sensors connected; stage temperature sensor attached.
- [ ] End switches wired; both trigger correctly by hand; the home switch is on the configured side.
- [ ] Config file filled in (lead, t_hold, v, soft limits, seed).
- [ ] Dry run with the motor disabled (or at reduced current, without the interferometer): the state machine, soft-limit checks and logging all work.

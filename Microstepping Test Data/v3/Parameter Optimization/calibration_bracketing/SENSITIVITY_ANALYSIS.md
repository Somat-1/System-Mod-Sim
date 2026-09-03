# Step 6 — sensitivity table: how it's set up

**Script**: `step6_sensitivity_table.py` | **Live status**: `sensitivity_status.json`
(written by the script while it runs) | **Raw output**: `step6_sensitivity_table.json`
(written on completion)

Scope: run 2 only (MRES 1/4, 100% current), way-port LuGre parameters only.
Run 2 was chosen as the run this identification effort is focusing on; the
way port because it's the dominant friction contributor visible across the
D-rate campaign.

## What's being measured

Six way-port LuGre parameters, each perturbed independently:

| Parameter | Meaning |
|---|---|
| `sigma0_way` | bristle stiffness |
| `sigma1_way` | bristle damping |
| `sigma2_way` | viscous term |
| `Fc_way` | Coulomb (kinetic) friction level |
| `Fs_way` | static friction level |
| `vs_way` | Stribeck velocity |

...against eight D-rate blocks (`D_0.125` ... `D_200`, the controller step
rates in full-steps/s), each a distinct move profile in the run-2 log.

**Cost metric**, computed per (parameter, sign, block) job:

```
cost = RMS( measured_position_um(t) - simulated_position_um(t) )
```

over the full block window, with the simulated trace linearly interpolated
onto the measured time grid (`np.interp`). Both traces are baselined against
the median of the first 20 measured samples in the window before subtracting.

**Sensitivity**, per (parameter, block) cell, from the ±10% perturbation:

```
sensitivity = 0.5 * ( |cost(+10%) - cost(base)| + |cost(-10%) - cost(base)| ) / cost(base)
```

a normalized, symmetric relative sensitivity. `NOISE_FLOOR = 0.01` (1%): a
parameter row with no column above this floor is flagged `FIX` (not
identifiable from this data) rather than `FIT`.

A separate check correlates the `Fs_way` and `sigma0_way` sensitivity
columns (Pearson r across the 8 blocks) — if they move together
(`corr > 0.9`), the two are collinear in this dataset and only their ratio
is identifiable, not both values independently.

## Job structure

```
8 baseline jobs   (one per block, no perturbation)
+ 6 parameters × 2 signs (+10%/-10%) × 8 blocks = 96 perturbed jobs
= 104 total jobs
```

Each job is one full nonlinear time-domain simulation of one block with the
Rev 4.2 model (`lugre_model_rev42.LuGreModelRev42`, 15 states), solved with
`scipy.integrate.solve_ivp(method='Radau')` using the model's analytical
Jacobian, `rtol=1e-6`, and a per-state `ATOL` vector tuned to each state's
physical scale. Per-run `T_hold` (hence `k_em`) is set from the actual SC
peak current logged for run 2, matching every other script in this
campaign.

## Parallelization and the redundant-I/O fix

Jobs run under a `ProcessPoolExecutor` (`workers = min(18, cpu_count-2)`).

The IDS position log (`SteppingSequenceID.csv`, ~2.78M samples) and the
controller event log live on the `\\mult-fp01.hitdom.lan` network share and
take **~45-55s to parse from cold** (`parse_ids` + `load_log_rows`). The
first version of this script parsed them fresh inside every job — with up
to 18 worker processes doing that concurrently, that redundant/contended
network I/O was the actual cause of the wildly inconsistent per-job times
observed during development (23s to 260s+ for nominally the same block).

Fixed by parsing **once**, in the main process, before submitting any jobs,
and handing the parsed arrays to every worker via the pool's `initializer`
(`_worker_init`, populating a module-level `_CACHED` dict once per worker
process at pool startup) instead of re-parsing per job. Verified with a
smoke test: `D_0.125` jobs dropped from 200s+ to 2-4s each; the genuinely
compute-heavy `D_3.5` block still took ~168s on its own, unchanged (its cost
is real ODE integration work, not I/O).

Blocks are **not** all equally expensive — `D_0.375`, `D_1.25`, `D_3.5`,
`D_27.5` in particular take from tens of seconds up to several minutes each
per job, driven by stiff bristle dynamics near presliding/direction
reversals, not by parameter choice. With 13 jobs per block (1 baseline + 12
perturbed) sharing 18 workers, total wall time is dominated by the slow
blocks and is expected to run well beyond an hour; there is no reliable
upfront estimate (see status-file caveat below), which is why progress is
tracked live instead.

## Live status file (`sensitivity_status.json`)

Single-writer (main process only, never workers) to avoid write races;
written atomically via a `.json.tmp` file + `Path.replace()` after every job
submission and every job completion, so a reader never sees a
partially-written file. Schema:

```jsonc
{
  "schema": "step6-sensitivity-status-v1",
  "note": "...",                      // caveat text, see below
  "started_at": "...", "last_updated": "...",
  "workers": 18,
  "total_jobs": 104, "completed_jobs": N, "failed_jobs": N, "pending_jobs": N,
  "elapsed_s": 123.4,
  "avg_completed_job_wall_s": 45.6,   // diagnostic only, see caveat
  "estimated_remaining_s": 789.0,     // throughput-based, see caveat
  "jobs": {
    "baseline|None|3.5": {"status": "done", "cost_um": 21.27, "job_wall_s": 168.3},
    "plus|sigma0_way|9.5": {"status": "running"},
    ...
  }
}
```

Two deliberate design choices worth calling out:

- **`estimated_remaining_s` uses overall throughput** (`elapsed_s /
  finished_jobs * remaining_jobs`), not the average of individual
  `job_wall_s` values. All 104 jobs are submitted to the pool up front, so
  a job's own `job_wall_s` (submit-to-completion) includes however long it
  sat queued behind other jobs on a busy worker — using that directly as
  an ETA basis would inflate over the course of the run as later jobs queue
  longer. Overall elapsed-time-per-completed-job already nets out the
  parallelism actually being achieved, so it's the better ETA basis.
  `job_wall_s` is kept per-job purely as a diagnostic (e.g. to spot which
  block/parameter combination is unusually slow), not fed into the ETA.
- **The ETA is explicitly labeled low-confidence** (`note` field) given the
  order-of-magnitude per-job timing variance seen during development on
  this network share. Treat it as a loose bound, not a real countdown.

## Outputs

- Console: baseline cost per block, the 6×8 sensitivity table with a
  `FIT`/`FIX` verdict per row, and the `Fs_way`/`sigma0_way` collinearity
  correlation.
- `step6_sensitivity_table.json`: raw `cost_um` for all 104 jobs, keyed
  `"{tag}|{param}|{block}"` (`tag` ∈ `baseline|plus|minus`), for any
  downstream re-analysis without rerunning the simulations.

## Status as of this writing

The first full launch (2026-09-03, ~08:43 UTC) was interrupted after 9/104
jobs completed — the background shell it ran under was torn down (session
teardown), not a script failure. `sensitivity_status.json` currently
reflects that partial, stale run. It needs to be relaunched from a session
that will stay alive for the full duration (or via a persistent/detached
launch) before the table can be completed.

"""Initial-testing parameter optimization of the Rev 4.2 nonlinear LuGre model
against a v6 measured segment.

Continues directly from sensitivity_dynamic_fit.py's segment choices: reuses
the same MRES 1 oscillation block (`osc_MRES1`), the strongest above-noise
signal of the four segments tested there, as the optimization target. Ramp
legs are explicitly out of scope for this run -- see README note at the
bottom of this file for why.

Model: `Rev 4/lugre_friction/Rev 4.2/scripts/lugre_model_rev42.py`
(LuGreModelRev42, 15 states, analytical Jacobian, same model
`v3/Parameter Optimization/calibration_bracketing/step6_sensitivity_table.py`
uses). theta_cmd is a piecewise-constant target rotor angle (a spring target,
per that model's own convention in generate_stepping_rev42.py's edge_list/
simulate_rev42) -- built here directly from the ALREADY-EXTRACTED oscillation
seed times in MRES_1/fitted_blocks.npz (measured, refined edge timing), not
re-derived from nominal campaign rates, so the forcing matches the actual
capture rather than an idealized schedule.

Optimized: three of the six way-port LuGre parameters (sigma0_way, Fc_way,
Fs_way) -- the "dominant friction contributor" per
calibration_bracketing/SENSITIVITY_ANALYSIS.md -- starting from the current
Rev 4.2 baseline (model_parameters.json), against the first 3 of 15
oscillation cycles (6 edges). sigma1_way, sigma2_way, vs_way are held fixed
at baseline. Both cuts (3 of 6 parameters, 3 of 15 cycles) are purely for
per-evaluation cost: one evaluation of this model over 6 edges takes ~10s
even at loosened solver tolerances (measured directly, see the timing probe
comment below), so a naive 6-parameter/12-edge attempt was projected to run
for hours, not minutes -- this is deliberately the largest version of this
fit that stays runnable as an initial test.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
CURVE_FIT_V1 = HERE.parent
REV42_SCRIPTS = Path(
    r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim"
    r"\Rev 4\lugre_friction\Rev 4.2\scripts"
)
sys.path.insert(0, str(REV42_SCRIPTS))
from lugre_model_rev42 import N_STATES, LuGreModelRev42, load_parameters  # noqa: E402

FULL_STEP_RAD = np.deg2rad(1.8)
MRES = 1
# Timing probe: 12 edges (6 cycles) at the tight tolerances below the
# generate_stepping_rev42.py convention uses cost 49s/evaluation -- a 6-param
# least_squares fit at that rate would run for hours. Cut to 3 cycles
# (6 edges) and looser tolerances (10s/evaluation, verified against the
# tight-tolerance trace to within ~0.6 um) to keep this an actually-runnable
# initial test rather than a job that never finishes.
N_CYCLES_FIT = 3
WAY_PARAMS = ["sigma0_way", "Fc_way", "Fs_way"]  # sigma1/sigma2/vs_way held at
                                                    # baseline -- see README
                                                    # note: cut for per-eval
                                                    # cost, not identifiability
RTOL = 1.0e-4
ATOL = np.full(N_STATES, 1.0e-8)
ATOL[12:15] = 1.0e-9  # z states (bristle deflection) stay tighter


def load_target():
    npz = np.load(CURVE_FIT_V1 / "MRES_1" / "fitted_blocks.npz")
    t = npz["oscillation_time_s"]
    measured_um = npz["oscillation_measured_um"]
    seed_local = npz["oscillation_seed_times_s"]
    n_edges = 2 * N_CYCLES_FIT
    edges = seed_local[:n_edges]
    window_end = edges[-1] + 0.6
    keep = t <= window_end
    return t[keep], measured_um[keep], edges


def edge_values(n_edges, step_rad):
    """theta_cmd after each edge: +step, 0, +step, 0, ... (forward, return)."""
    return np.array([step_rad if i % 2 == 0 else 0.0 for i in range(n_edges)])


def simulate(params: dict, edges: np.ndarray, window_end: float, out_dt: float = 5.0e-3):
    model = LuGreModelRev42(parameters=params, enforce_interface_power=False)
    values = edge_values(edges.size, FULL_STEP_RAD / MRES)
    bounds_t = np.concatenate([[0.0], edges, [window_end]])
    commands = np.concatenate([[0.0], values])  # held BEFORE each edge in bounds_t[1:]

    state = np.zeros(N_STATES)
    t_parts, x_parts = [], []
    jac = lambda _t, y, m=model: m.analytical_linearization(y)[0]
    for lo, hi, cmd in zip(bounds_t[:-1], bounds_t[1:], commands):
        n_pts = max(2, int(round((hi - lo) / out_dt)) + 1)
        t_eval = np.linspace(lo, hi, n_pts)
        sol = solve_ivp(
            lambda t, y, c=cmd: model.rhs(t, y, c),
            (lo, hi), state, method="Radau", jac=jac,
            t_eval=t_eval, rtol=RTOL, atol=ATOL,
        )
        if not sol.success:
            raise RuntimeError(f"integration failed on [{lo:.3f},{hi:.3f}]: {sol.message}")
        state = sol.y[:, -1]
        keep = slice(0, -1) if hi < bounds_t[-1] else slice(0, None)
        t_parts.append(t_eval[keep])
        x_parts.append(sol.y[5, keep])
    return np.concatenate(t_parts), np.concatenate(x_parts) * 1.0e6  # x_n, um


def main():
    t_meas, measured_um, edges = load_target()
    window_end = t_meas[-1]
    base_params = load_parameters()

    print(f"target: MRES {MRES} oscillation, first {N_CYCLES_FIT} cycles "
          f"({edges.size} edges), window [0, {window_end:.3f}] s, "
          f"{t_meas.size} measured samples")

    t0 = time.perf_counter()
    t_sim, x_base_um = simulate(base_params, edges, window_end)
    one_eval_s = time.perf_counter() - t0
    print(f"one model evaluation: {one_eval_s:.2f} s ({t_sim.size} output points)")

    measured_on_sim_grid = np.interp(t_sim, t_meas, measured_um)
    base_cost = float(np.sqrt(np.mean((measured_on_sim_grid - x_base_um) ** 2)))
    print(f"baseline RMS-to-measured: {base_cost:.4f} um")

    budget_s = 1500.0
    max_evals = max(6, int(budget_s / max(one_eval_s, 0.05)))
    print(f"per-eval time suggests a budget of ~{max_evals} total evaluations "
          f"for an {budget_s:.0f} s run")

    # Round 1 pinned Fc_way/Fs_way at a 5x-baseline bound while still
    # climbing -- xtol termination there was a wall, not a real optimum.
    # Round 2: much wider bounds (0.05x-20x), warm-started from round 1's
    # result instead of the original baseline, so this run spends its
    # budget extending the search past the old wall instead of re-deriving
    # the direction it already found.
    round1_result_path = HERE / "model_optimization_results.json"
    if round1_result_path.exists():
        prior = json.loads(round1_result_path.read_text(encoding="utf-8"))
        x0 = np.array([prior["parameters"][p]["optimized"] for p in WAY_PARAMS])
        print(f"warm-starting from round 1 result: {x0}")
    else:
        x0 = np.array([base_params[p] for p in WAY_PARAMS])
    baseline_x0 = np.array([base_params[p] for p in WAY_PARAMS])
    lower = baseline_x0 * 0.05
    upper = baseline_x0 * 20.0

    history = []
    cache: dict[tuple, np.ndarray] = {}

    def residuals(x):
        key = tuple(round(float(v), 6) for v in x)
        if key in cache:
            return cache[key]
        params = dict(base_params)
        for name, val in zip(WAY_PARAMS, x):
            params[name] = float(val)
        t_sim, x_sim_um = simulate(params, edges, window_end)
        target = np.interp(t_sim, t_meas, measured_um)
        res = target - x_sim_um
        cost = float(np.sqrt(np.mean(res ** 2)))
        history.append({"x": x.tolist(), "cost_um": cost,
                        "elapsed_s": time.perf_counter() - opt_start})
        print(f"  eval {len(history):3d}: cost={cost:.4f} um  "
              f"x={np.array2string(x, precision=4)}")
        cache[key] = res
        return res

    opt_start = time.perf_counter()
    result = least_squares(
        residuals, x0, bounds=(lower, upper), method="trf",
        diff_step=0.05, max_nfev=max_evals, verbose=0,
    )
    opt_elapsed = time.perf_counter() - opt_start
    print(f"cache: {len(history)} real evaluations, "
          f"{result.nfev - len(history)} deduplicated")

    optimized_params = dict(base_params)
    for name, val in zip(WAY_PARAMS, result.x):
        optimized_params[name] = float(val)
    t_sim_opt, x_opt_um = simulate(optimized_params, edges, window_end)
    final_cost = float(np.sqrt(np.mean(
        (np.interp(t_sim_opt, t_meas, measured_um) - x_opt_um) ** 2)))

    print(f"\noptimization finished: success={result.success}, "
          f"nfev={result.nfev}, elapsed={opt_elapsed:.1f}s")
    print(f"cost: {base_cost:.4f} -> {final_cost:.4f} um "
          f"({100 * (1 - final_cost / base_cost):.1f}% reduction)")
    for name, before, warm, after in zip(WAY_PARAMS, baseline_x0, x0, result.x):
        print(f"  {name:12s} original_baseline={before:.6g}  "
              f"round1_warm_start={warm:.6g}  round2_final={after:.6g} "
              f"({100 * (after / before - 1):+.1f}% vs original baseline)")

    fig, ax = plt.subplots(figsize=(11, 4.4))
    ax.plot(t_meas, measured_um, color="#777777", lw=0.6, alpha=0.6, label="measured")
    ax.plot(t_sim, x_base_um, color="#9a9a9a", lw=1.3, ls="--",
            label=f"baseline params (RMS {base_cost:.3f} µm)")
    ax.plot(t_sim_opt, x_opt_um, color="#c0392b", lw=1.4,
            label=f"optimized params (RMS {final_cost:.3f} µm)")
    for e in edges:
        ax.axvline(e, color="#bd6b18", alpha=0.3, lw=0.8, ls=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Position (µm)")
    ax.set_title(f"MRES {MRES} oscillation, first {N_CYCLES_FIT} cycles — "
                f"Rev 4.2 way-port parameter fit (initial testing)", fontsize=10.5)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "plots" / "model_fit_overlay.png", dpi=160)
    plt.close(fig)

    out = {
        "schema": "v6-rev42-way-port-optimization-v1",
        "target": {"mres": MRES, "n_cycles_fit": N_CYCLES_FIT,
                   "n_edges": int(edges.size), "window_end_s": float(window_end),
                   "source": "MRES_1/fitted_blocks.npz (oscillation)"},
        "model": "Rev 4.2 LuGreModelRev42 (way-port parameters only)",
        "one_eval_s": one_eval_s,
        "optimizer": {"method": "trf", "diff_step": 0.05, "max_nfev": max_evals,
                     "actual_nfev": result.nfev, "success": bool(result.success),
                     "status": int(result.status), "message": result.message,
                     "elapsed_s": opt_elapsed},
        "cost_um": {"baseline": base_cost, "final": final_cost},
        "bounds": {"lower": lower.tolist(), "upper": upper.tolist(),
                  "relative_to_original_baseline": [0.05, 20.0]},
        "parameters": {
            name: {"original_baseline": float(b), "warm_start": float(w),
                  "optimized": float(a)}
            for name, b, w, a in zip(WAY_PARAMS, baseline_x0, x0, result.x)
        },
        "history": history,
    }
    (HERE / "model_optimization_results.json").write_text(
        json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote model_optimization_results.json, "
          f"plots/model_fit_overlay.png -> {HERE}")


if __name__ == "__main__":
    main()

# NOTE on scope: ramp legs are NOT attempted in this initial run. A 25 mm leg
# at MRES 1 needs TRAJECTORY_FULL_STEPS/MRES edges (2500+ for the v4 rig's
# geometry; v6 is expected to be the same order) -- 100-200x this block's 12
# edges -- and each edge here already costs a full Radau solve. Extrapolating
# this run's one_eval_s to that edge count is the honest way to decide
# whether a ramp-leg fit is even feasible before attempting one, rather than
# launching a multi-hour job blind.

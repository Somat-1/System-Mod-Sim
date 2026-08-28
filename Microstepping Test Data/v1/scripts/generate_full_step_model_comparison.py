#!/usr/bin/env python3
"""Reproduce the REAL hardware back-and-forth stepping sequence from each
"Microstepping Test Data/v1/StepSize{1,2,16}.csv" run (command edge times and
directions detected directly from the encoder trace, not an idealized
sequence), then run that exact per-file edge timing through all three Rev 4
model families and overlay the result against the real IDS interferometer
displacement for the same run:

  - Newton/free-body linear frictionless baseline, after verifying that its
    state-space realization is identical to the Lagrange realization
  - otherwise frictionless dynamics with the exact fixed-frame nonlinear
    detent torque T_d*sin(4*N_r*theta_m) replacing the constant k_d tangent
  - Rev 4.2 nonlinear parallel LuGre friction model (Rev 4/lugre_friction/
    Rev 4.2/scripts/lugre_model_rev42.py -- imported via an explicit
    sys.path insert below, since this script itself lives under
    Microstepping Test Data/scripts/, not next to that module)

StepSize8 dropped (2026-08-27): its encoder CSV was flagged as faulty and
removed from the data folder (git history still has it if ever needed).

Edge detection: each StepSize{N}.csv encoder trace turns out to be a simple
alternating +1/-1 step command at a consistent ~1.17 s cadence, confirmed by
peak-detecting the median-filtered counter derivative with a per-file
magnitude threshold (chosen to separate genuine command edges from small
sub-100-count post-reversal settling blips visible at the finer step sizes,
not to fit a particular expected pattern) -- all three remaining files show
exactly 26 edges (13 full there-and-back cycles) at a 1.170 s mean gap,
spanning the ~30 s test the file names/README describe. That detected edge
list (real times, real directions) is what is actually simulated here -- this
is not the idealized MOVES staircase used elsewhere in Rev 4/Rev 4.2.

Physical step size per file: StepSize{N} = one full mechanical step (5.000
um, verified elsewhere against state_space_6dof.md) divided by N.

Alignment (2026-08-27): the encoder and IDS interferometer are two
independent logging systems with their own, offset t=0 references (confirmed
by comparing StepSize1.csv's "Starttime of export" header against
IDSdata.txt's "Date:" header for the same run -- a few seconds apart), and
the IDS's own polarity convention is opposite the model's. Both are corrected
per step size before plotting: the real trace is time-shifted so its own
first detected step lines up with the model's first commanded edge, and
negated if its first step direction is opposite the model's.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import StateSpace, lsim, medfilt, find_peaks

VERSION_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = VERSION_ROOT.parents[1]                    # .../Sytem Mod & Sim
REV4 = PROJECT_ROOT / "Rev 4"
LUGRE_REV42_SCRIPTS = REV4 / "lugre_friction" / "Rev 4.2" / "scripts"
DATA_DIR = VERSION_ROOT
OUT_DIR = DATA_DIR / "rendered_assets"
FIGURE = OUT_DIR / "full_step_model_comparison.png"
DATA = OUT_DIR / "full_step_model_comparison.npz"

sys.path.insert(0, str(LUGRE_REV42_SCRIPTS))
from lugre_model_rev42 import N_Q, N_STATES, LuGreModelRev42  # noqa: E402

SUMMARY = OUT_DIR / "full_step_model_comparison_summary.json"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_microstepping_test_data import parse_encoder_csv, parse_ids_displacement  # noqa: E402

STEP_SIZES = [1, 2, 16]
FULL_STEP_UM = 5.0  # verified against state_space_6dof.md / Rev 3 doc elsewhere
EDGE_THRESHOLD_COUNTS = {1: 800, 2: 400, 16: 100}
COLORS = {1: "#2b6cb0", 2: "#c05621", 16: "#805ad5"}

OUTPUT_DT_S = 1.0e-3
RTOL = 1.0e-6
ATOL = np.array([
    1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-12, 1.0e-12,
    1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-9, 1.0e-9,
    1.0e-11, 1.0e-11, 1.0e-9,
])


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def detect_command_edges(t: np.ndarray, delta: np.ndarray, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    """Real command edge times and signs, detected from the encoder trace --
    not assumed. Median-filter to suppress single-sample noise, then peak-find
    the derivative with a magnitude threshold tuned per file to separate
    genuine command steps from small settling blips (verified against the
    raw traces, not fit to a target answer)."""
    smooth = medfilt(delta.astype(float), kernel_size=5)
    d = np.diff(smooth)
    peaks, _ = find_peaks(np.abs(d), height=threshold, distance=50)
    edge_times = t[peaks + 1]
    signs = np.sign(d[peaks])
    return edge_times, signs


def detect_ids_first_edge(t: np.ndarray, disp_um: np.ndarray) -> tuple[float, float]:
    """First genuine step edge (time, sign) in the real IDS displacement
    trace. Threshold is set relative to this trace's own excursion amplitude
    (25% of the 95th-percentile deviation from the median) rather than a
    hardcoded number, since the IDS signal amplitude varies by more than 10x
    across step sizes."""
    d = np.diff(disp_um)
    amplitude = np.percentile(np.abs(disp_um - np.median(disp_um)), 95)
    threshold = 0.25 * amplitude
    peaks, _ = find_peaks(np.abs(d), height=threshold, distance=3)
    if len(peaks) == 0:
        raise RuntimeError("No IDS step edge detected -- threshold too strict for this trace")
    return float(t[peaks[0] + 1]), float(np.sign(d[peaks[0]]))


def align_real_trace(real_t: np.ndarray, real_disp_um: np.ndarray,
                      model_first_edge_t: float, model_first_sign: float) -> tuple[np.ndarray, np.ndarray]:
    """Time-shift the real trace so its own first detected step lines up with
    the model's first commanded edge, and flip its sign if its first step
    direction is opposite the model's -- the encoder and IDS loggers have
    independent t=0 references and the IDS has its own polarity convention,
    neither of which is a modeling discrepancy."""
    ids_first_t, ids_first_sign = detect_ids_first_edge(real_t, real_disp_um)
    shift = model_first_edge_t - ids_first_t
    aligned_t = real_t + shift
    aligned_disp = real_disp_um if ids_first_sign == model_first_sign else -real_disp_um
    return aligned_t, aligned_disp


def build_segments(edge_times: np.ndarray, signs: np.ndarray, step_m: float,
                    lead: float, total_duration_s: float) -> list[tuple[float, float]]:
    """[(dwell_s, theta_command), ...] -- a leading zero-command pre-roll from
    t=0 to the first detected edge, one segment per detected edge (command
    holds at the new position until the next edge), and a trailing segment
    holding the final commanded position out to the end of the real test."""
    segments = [(float(edge_times[0]), 0.0)]
    travel = 0.0
    for i in range(len(edge_times)):
        travel += signs[i] * step_m
        theta_command = travel / lead
        end_t = edge_times[i + 1] if i + 1 < len(edge_times) else total_duration_s
        dwell = max(end_t - edge_times[i], OUTPUT_DT_S)
        segments.append((float(dwell), float(theta_command)))
    return segments


def simulate_lugre_segments(model: LuGreModelRev42, segments: list[tuple[float, float]], lead: float) -> dict:
    state = np.zeros(N_STATES)
    time_parts, position_parts = [], []
    t_cursor = 0.0
    nfev = njev = nlu = 0
    started = time.perf_counter()

    for i, (dwell, theta_command) in enumerate(segments):
        n_out = max(int(round(dwell / OUTPUT_DT_S)), 1)
        local_time = np.linspace(0.0, dwell, n_out + 1)
        solution = solve_ivp(
            lambda t, y, command=theta_command: model.rhs(t, y, command),
            (0.0, dwell), state, method="Radau",
            jac=lambda _t, y: model.analytical_linearization(y)[0],
            t_eval=local_time, rtol=RTOL, atol=ATOL,
        )
        if not solution.success:
            raise RuntimeError(f"LuGre integration failed at segment {i}: {solution.message}")
        state = solution.y[:, -1]
        nfev += solution.nfev
        njev += solution.njev
        nlu += solution.nlu
        time_parts.append(t_cursor + local_time[:-1])
        position_parts.append(solution.y[5, :-1])
        t_cursor += dwell
        print(f"  LuGre segment {i + 1}/{len(segments)}; elapsed {time.perf_counter() - started:.1f} s",
              flush=True)

    return {
        "time_s": np.concatenate(time_parts),
        "x_n_m": np.concatenate(position_parts),
        "elapsed_s": time.perf_counter() - started,
        "nfev": nfev, "njev": njev, "nlu": nlu,
    }


def simulate_full_detent_segments(
    model: LuGreModelRev42,
    segments: list[tuple[float, float]],
) -> dict:
    """No LuGre forces; retain the exact fixed-frame periodic detent."""
    state = np.zeros(2 * N_Q)
    time_parts, position_parts = [], []
    t_cursor = 0.0
    nfev = njev = nlu = 0
    started = time.perf_counter()

    def rhs(y: np.ndarray, theta_command: float) -> np.ndarray:
        q = y[:N_Q]
        velocity = y[N_Q:]
        detent = np.zeros(N_Q)
        detent[0] = model.p["T_d"] * np.sin(
            4.0 * model.p["N_r"] * q[0]
        )
        acceleration = model.mass_inverse @ (
            model.command * theta_command
            - model.damping @ velocity
            - model.stiffness @ q
            - detent
        )
        return np.concatenate([velocity, acceleration])

    def jacobian(y: np.ndarray) -> np.ndarray:
        position_block = -model.stiffness.copy()
        position_block[0, 0] -= (
            4.0 * model.p["N_r"] * model.p["T_d"]
            * np.cos(4.0 * model.p["N_r"] * y[0])
        )
        return np.block([
            [np.zeros((N_Q, N_Q)), np.eye(N_Q)],
            [
                model.mass_inverse @ position_block,
                -model.mass_inverse @ model.damping,
            ],
        ])

    for i, (dwell, theta_command) in enumerate(segments):
        n_out = max(int(round(dwell / OUTPUT_DT_S)), 1)
        local_time = np.linspace(0.0, dwell, n_out + 1)
        solution = solve_ivp(
            lambda _t, y, command=theta_command: rhs(y, command),
            (0.0, dwell), state, method="Radau",
            jac=lambda _t, y: jacobian(y),
            t_eval=local_time, rtol=RTOL, atol=ATOL[:2 * N_Q],
        )
        if not solution.success:
            raise RuntimeError(
                f"Full-detent integration failed at segment {i}: "
                f"{solution.message}"
            )
        state = solution.y[:, -1]
        nfev += solution.nfev
        njev += solution.njev
        nlu += solution.nlu
        time_parts.append(t_cursor + local_time[:-1])
        position_parts.append(solution.y[5, :-1])
        t_cursor += dwell
        print(
            f"  Full detent segment {i + 1}/{len(segments)}; "
            f"elapsed {time.perf_counter() - started:.1f} s",
            flush=True,
        )

    return {
        "time_s": np.concatenate(time_parts),
        "x_n_m": np.concatenate(position_parts),
        "elapsed_s": time.perf_counter() - started,
        "nfev": nfev, "njev": njev, "nlu": nlu,
    }


def command_trace(
    segments: list[tuple[float, float]], lead: float,
) -> tuple[np.ndarray, np.ndarray]:
    time_parts, travel_parts = [], []
    t_cursor = 0.0
    for dwell, theta_command in segments:
        n_out = max(int(round(dwell / OUTPUT_DT_S)), 1)
        local_time = np.linspace(0.0, dwell, n_out + 1)
        time_parts.append(t_cursor + local_time[:-1])
        travel_parts.append(
            np.full(local_time.size - 1, lead * theta_command)
        )
        t_cursor += dwell
    return np.concatenate(time_parts), np.concatenate(travel_parts)


def second_order_state_space(mass, damping, stiffness, command) -> StateSpace:
    n = mass.shape[0]
    inverse_mass = np.linalg.inv(mass)
    system = np.block([
        [np.zeros((n, n)), np.eye(n)],
        [-inverse_mass @ stiffness, -inverse_mass @ damping],
    ])
    input_vector = np.concatenate([np.zeros(n), inverse_mass @ command]).reshape(-1, 1)
    output = np.zeros((1, 2 * n))
    output[0, 5] = 1.0
    return StateSpace(system, input_vector, output, np.zeros((1, 1)))


def simulate_linear_segments(system: StateSpace, segments: list[tuple[float, float]]) -> dict:
    state = np.zeros(system.A.shape[0])
    time_parts, position_parts = [], []
    t_cursor = 0.0
    for dwell, theta_command in segments:
        n_out = max(int(round(dwell / OUTPUT_DT_S)), 1)
        local_time = np.linspace(0.0, dwell, n_out + 1)
        values = np.full(local_time.size, theta_command)
        _, output, states = lsim(system, U=values, T=local_time, X0=state)
        state = states[-1]
        time_parts.append(t_cursor + local_time[:-1])
        position_parts.append(np.asarray(output[:-1]).reshape(-1))
        t_cursor += dwell
    return {"time_s": np.concatenate(time_parts), "x_n_m": np.concatenate(position_parts)}


def full_step_travel_m(newton_parameters: dict) -> float:
    lead_ratio = newton_parameters["L"] / (2.0 * np.pi)
    full_step_angle = 2.0 * np.pi / (4.0 * newton_parameters["N_r"])
    return lead_ratio * full_step_angle


def baseline_systems():
    newton = load_module("newton_rev4", REV4 / "scripts" / "build_bode_rev4.py")
    newton_parameters = newton.load_parameters()
    mass_n, damping_n, stiffness_n, inputs_n = newton.build_matrices(newton_parameters)
    command_column = newton.INPUT_LABELS.index("theta_cmd")
    command_n = inputs_n[:, command_column].copy()
    command_n[0] = newton_parameters["N_r"] * newton_parameters["T_hold"]
    newton_system = second_order_state_space(mass_n, damping_n, stiffness_n, command_n)

    lagrange = load_module(
        "lagrange_rev4",
        REV4 / "Lagrange Derivation" / "scripts" / "build_bode_lagrange_frictionless.py",
    )
    lagrange_parameters = lagrange.load_parameters()
    mass_l, damping_l, stiffness_l, command_l = lagrange.build_lagrange_matrices(lagrange_parameters)
    lagrange_system = second_order_state_space(mass_l, damping_l, stiffness_l, command_l)

    return newton_system, lagrange_system, newton_parameters


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    newton_system, lagrange_system, newton_parameters = baseline_systems()
    realization_difference = max(
        float(np.max(np.abs(
            getattr(newton_system, name) - getattr(lagrange_system, name)
        )))
        for name in ("A", "B", "C", "D")
    )
    if realization_difference > 1.0e-12:
        raise AssertionError(
            "Newton and Lagrange frictionless realizations no longer match: "
            f"maximum matrix difference {realization_difference:.3e}"
        )
    print(
        "Newton/Lagrange frictionless state-space identity verified; "
        f"maximum matrix difference {realization_difference:.3e}.",
        flush=True,
    )
    step_travel_full_m = full_step_travel_m(newton_parameters)
    if not np.isclose(step_travel_full_m, 5.0e-6, rtol=1.0e-6):
        raise AssertionError(
            f"Full-step travel {step_travel_full_m*1e6:.4f} um drifted from the documented 5.000 um"
        )

    model = LuGreModelRev42(enforce_interface_power=False)
    lead = model.p["L"] / (2.0 * np.pi)

    ids_blocks = parse_ids_displacement(DATA_DIR / "IDSdata.txt")

    fig, axes = plt.subplots(1, len(STEP_SIZES), figsize=(5.0 * len(STEP_SIZES), 5.5))
    axes_flat = iter(axes)
    summary_per_size = {}
    data_payload = {}

    for n in STEP_SIZES:
        step_m = step_travel_full_m / n
        t_enc, delta = parse_encoder_csv(DATA_DIR / f"StepSize{n}.csv")
        edge_times, signs = detect_command_edges(t_enc, delta, EDGE_THRESHOLD_COUNTS[n])
        total_duration_s = float(t_enc[-1])
        print(f"\nStepSize{n}: {len(edge_times)} detected edges, step={step_m*1e6:.4f} um, "
              f"duration={total_duration_s:.2f} s", flush=True)

        segments = build_segments(edge_times, signs, step_m, lead, total_duration_s)

        print("  Running Newton linear frictionless baseline...", flush=True)
        frictionless = simulate_linear_segments(newton_system, segments)

        print(
            f"  Running frictionless full periodic detent "
            f"({len(segments)} segments, nonlinear Radau)...",
            flush=True,
        )
        full_detent = simulate_full_detent_segments(model, segments)

        print(f"  Running Rev 4.2 LuGre ({len(segments)} segments, nonlinear Radau)...", flush=True)
        lugre = simulate_lugre_segments(model, segments, lead)
        command_t, command_m = command_trace(segments, lead)

        ax = next(axes_flat)
        real_t, real_disp_um = ids_blocks[n]
        real_t, real_disp_um = align_real_trace(real_t, real_disp_um, edge_times[0], signs[0])
        ax.plot(
            command_t, command_m * 1.0e6, color="#888888",
            linewidth=0.8, linestyle=":", label="reconstructed command",
        )
        ax.plot(real_t, real_disp_um, color="#333333", linewidth=0.8, alpha=0.75,
                label="IDS interferometer (measured, aligned)")
        ax.plot(
            frictionless["time_s"], frictionless["x_n_m"] * 1.0e6,
            color="#2b6cb0", linewidth=1.1,
            label="frictionless linear tangent (Newton)",
        )
        ax.plot(
            full_detent["time_s"], full_detent["x_n_m"] * 1.0e6,
            color="#d97706", linewidth=1.05, linestyle="-.",
            label="frictionless full periodic detent",
        )
        ax.plot(lugre["time_s"], lugre["x_n_m"] * 1.0e6, color="#c0392b",
                linewidth=1.0, label="Rev 4.2 LuGre + detent")
        ax.set_title(f"StepSize{n} ({step_m*1e6:.4f} $\\mu$m/edge)", fontsize=10)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"Stage displacement ($\mu$m)")
        ax.grid(True, linewidth=0.4, color="#cccccc")
        if n == STEP_SIZES[0]:
            ax.legend(fontsize=7.5, loc="lower right")

        ids_first_t_raw, ids_first_sign_raw = detect_ids_first_edge(*ids_blocks[n])
        summary_per_size[f"StepSize{n}"] = {
            "step_um": step_m * 1.0e6,
            "detected_edges": int(len(edge_times)),
            "mean_edge_gap_s": float(np.mean(np.diff(edge_times))),
            "total_duration_s": total_duration_s,
            "model_first_edge_s": float(edge_times[0]),
            "model_first_sign": float(signs[0]),
            "ids_first_edge_s_raw": ids_first_t_raw,
            "ids_first_sign_raw": ids_first_sign_raw,
            "ids_time_shift_applied_s": float(edge_times[0] - ids_first_t_raw),
            "ids_sign_flipped": bool(ids_first_sign_raw != signs[0]),
            "linear_baseline": (
                "Newton; Lagrange state-space verified identical"
            ),
            "newton_lagrange_state_space_max_difference": (
                realization_difference
            ),
            "full_detent_elapsed_s": full_detent["elapsed_s"],
            "full_detent_nfev": full_detent["nfev"],
            "lugre_elapsed_s": lugre["elapsed_s"],
            "lugre_nfev": lugre["nfev"],
            "max_linear_vs_full_detent_um": float(
                np.max(np.abs(
                    frictionless["x_n_m"] - full_detent["x_n_m"]
                )) * 1.0e6
            ),
            "max_full_detent_vs_lugre_um": float(
                np.max(np.abs(
                    full_detent["x_n_m"] - lugre["x_n_m"]
                )) * 1.0e6
            ),
        }

        prefix = f"stepsize{n}"
        data_payload.update({
            f"{prefix}_edge_times_s": edge_times,
            f"{prefix}_edge_signs": signs,
            f"{prefix}_command_time_s": command_t,
            f"{prefix}_command_travel_m": command_m,
            f"{prefix}_ids_time_aligned_s": real_t,
            f"{prefix}_ids_displacement_aligned_um": real_disp_um,
            f"{prefix}_frictionless_time_s": frictionless["time_s"],
            f"{prefix}_frictionless_x_n_m": frictionless["x_n_m"],
            f"{prefix}_full_detent_time_s": full_detent["time_s"],
            f"{prefix}_full_detent_x_n_m": full_detent["x_n_m"],
            f"{prefix}_lugre_time_s": lugre["time_s"],
            f"{prefix}_lugre_x_n_m": lugre["x_n_m"],
        })

    fig.suptitle(
        "Real hardware stepping tests: linear frictionless / full detent / LuGre 4.2\n"
        "command sequence reproduced from each file's own detected encoder edges; "
        "real trace time- and sign-aligned to the model's first edge -- "
        "Microstepping Test Data/ (StepSize8 excluded, faulty)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.88])
    temporary_figure = FIGURE.with_name(
        f"{FIGURE.stem}.tmp{FIGURE.suffix}"
    )
    fig.savefig(temporary_figure, dpi=130)
    plt.close(fig)
    temporary_figure.replace(FIGURE)

    temporary_data = DATA.with_name(f"{DATA.stem}.tmp{DATA.suffix}")
    np.savez_compressed(temporary_data, **data_payload)
    temporary_data.replace(DATA)

    temporary_summary = SUMMARY.with_name(
        f"{SUMMARY.stem}.tmp{SUMMARY.suffix}"
    )
    temporary_summary.write_text(
        json.dumps(summary_per_size, indent=2) + "\n", encoding="utf-8"
    )
    temporary_summary.replace(SUMMARY)
    print(json.dumps(summary_per_size, indent=2), flush=True)
    print(f"Wrote {FIGURE}")
    print(f"Wrote {DATA}")
    print(f"Wrote {SUMMARY}")


if __name__ == "__main__":
    main()

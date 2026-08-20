#!/usr/bin/env python3
"""Broadband frequency-response identification of the nonlinear 15-state
LuGre model, via a linear chirp and a PRBS excitation of theta_cmd, at
several RMS amplitudes bracketing the guideway breakaway threshold.

See broadband_signals.py and broadband_estimators.py module docstrings for
the excitation design and estimator rationale; lugre_model_broadband.py for
why the electromagnetic/detent torques must be nonlinear (sin), not the
linearized k_EM/k_d, at these amplitudes; broadband_worker.py for the
per-run simulation pipeline itself.

The 5 amplitudes x {chirp, prbs} = 10 runs are fully independent (each is
its own solve_ivp call from t=0), so they're dispatched across a
ProcessPoolExecutor rather than run serially -- a single PRBS run measured
~5 minutes wall time, and the 60 s chirps are considerably more, so serial
execution would run into hours. Same single-threaded-BLAS-before-numpy
pattern as the earlier sinusoidal sweep's worker (see broadband_worker.py).

Solver: LSODA only (Radau/BDF diverge on this model's v=0 kink -- see
backlog.md / lugre_friction/README.md).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from broadband_signals import F_HI_HZ, F_LO_HZ, choose_rms_amplitudes
from broadband_worker import init_worker, run_one
from lugre_model_broadband import LuGreModelNonlinearDrive

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"

N_WORKERS = max(1, (os.cpu_count() or 4) - 2)


def build_work_items(amplitudes: np.ndarray) -> list[dict]:
    items = [{"A": float(A), "signal_type": st} for A in amplitudes for st in ("chirp", "prbs")]
    # Chirp (60 s/run) is far more expensive than PRBS (12 PRBS periods,
    # ~11 s/run) -- chirp first so the pool isn't left running one alone
    # after every PRBS job has already finished.
    items.sort(key=lambda it: (it["signal_type"] != "chirp", -it["A"]))
    return items


def main() -> None:
    model = LuGreModelNonlinearDrive()
    p = model.p

    amplitudes, a_thresh = choose_rms_amplitudes(p)
    a_kinematic_5um = 5e-6 * 2.0 * np.pi / p["L"]
    print(f"Guideway breakaway command angle: {a_thresh*1e3:.3f} mrad "
          f"(sanity check: 5 um travel = {a_kinematic_5um*1e3:.3f} mrad -- should be same order)")
    print(f"Chosen RMS amplitudes ({len(amplitudes)}, log-spaced, +-1 decade around threshold):")
    for A in amplitudes:
        print(f"  {A*1e3:9.4f} mrad RMS  (chirp peak {A*np.sqrt(2)*1e3:9.4f} mrad, "
              f"PRBS peak {A*1e3:9.4f} mrad)")

    work_items = build_work_items(amplitudes)
    print(f"\n{len(work_items)} work items, {N_WORKERS} worker processes")

    results_by_key: dict[tuple[float, str], dict] = {}
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=init_worker) as ex:
        futures = [ex.submit(run_one, item) for item in work_items]
        n_done = 0
        for fut in as_completed(futures):
            r = fut.result()
            results_by_key[(r["A"], r["signal_type"])] = r
            n_done += 1
            elapsed = time.time() - t_start
            sat = r["saturation"]
            print(f"  {n_done:2d}/{len(work_items)} done ({elapsed:7.1f}s elapsed): "
                  f"A={r['A']*1e3:9.4f} mrad {r['signal_type']:5s} "
                  f"{r['wall_s']:7.1f}s  sat sb={sat['sb']:.3f} nut={sat['nut']:.3f} way={sat['way']:.3f}",
                  flush=True)
    total_wall = time.time() - t_start
    print(f"\nSweep done in {total_wall:.1f}s")

    prbs_results = [results_by_key[(float(A), "prbs")] for A in amplitudes]
    chirp_results = [results_by_key[(float(A), "chirp")] for A in amplitudes]

    ASSET_DIR.mkdir(exist_ok=True)
    NPZ_DIR = ASSET_DIR / "npz"
    NPZ_DIR.mkdir(exist_ok=True)
    npz_payload = {"amplitudes": amplitudes, "breakaway_angle": a_thresh}
    for i, A in enumerate(amplitudes):
        rp, rc = prbs_results[i], chirp_results[i]
        npz_payload[f"prbs_f_amp{i}"] = rp["f"]
        npz_payload[f"prbs_G_amp{i}"] = rp["G"]
        npz_payload[f"prbs_gamma2_amp{i}"] = rp["gamma2"]
        npz_payload[f"chirp_f_amp{i}"] = rc["f"]
        npz_payload[f"chirp_mag_db_amp{i}"] = rc["mag_db"]
        npz_payload[f"chirp_phase_deg_amp{i}"] = rc["phase_deg"]
        for port in ("sb", "nut", "way"):
            npz_payload[f"prbs_sat_{port}_amp{i}"] = rp["saturation"][port]
            npz_payload[f"chirp_sat_{port}_amp{i}"] = rc["saturation"][port]
    np.savez(NPZ_DIR / "broadband_id_data.npz", **npz_payload)

    # ---- Plot: magnitude, phase, coherence -- solid=PRBS, dashed=chirp ----
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, len(amplitudes) - 1)) for i in range(len(amplitudes))]
    COH_GATE = 0.9

    fig, (ax_mag, ax_phase, ax_coh) = plt.subplots(3, 1, figsize=(10.0, 11.0), sharex=True)
    for i, A in enumerate(amplitudes):
        color = colors[i]
        rp, rc = prbs_results[i], chirp_results[i]

        # Welch/rfft return frequencies out to the full Nyquist (fs/2 = 10 kHz),
        # but nothing was ever excited above F_HI_HZ (2000 Hz) -- that region is
        # pure noise floor (PRBS's sinc envelope has rolled off, chirp never
        # swept there). Slicing to the actual excitation band BEFORE plotting
        # (not after) also keeps that noise out of the phase display; np.unwrap
        # is sequential, so noise-driven runaway beyond F_HI_HZ never contaminated
        # the in-band unwrapped values themselves, but plotting the full range
        # made the in-band data unreadable next to it.
        mp = (rp["f"] >= F_LO_HZ) & (rp["f"] <= F_HI_HZ)
        mc = (rc["f"] >= F_LO_HZ) & (rc["f"] <= F_HI_HZ)
        fp, fc_ = rp["f"][mp], rc["f"][mc]

        prbs_mag_db = 20.0 * np.log10(np.maximum(np.abs(rp["G"][mp]), 1e-300))
        # Wrapped to (-180, 180], not raw np.unwrap: this is a broadband/nonlinear
        # G(f), not a single smooth linear resonance chain -- unwrap's "phase
        # changes slowly between adjacent bins" assumption doesn't hold well
        # here, and raw unwrap accumulated into the thousands of degrees even
        # within-band, which is unreadable and not obviously more meaningful
        # than the wrapped value.
        prbs_phase_deg = ((np.angle(rp["G"][mp]) * 180.0 / np.pi + 180.0) % 360.0) - 180.0
        chirp_phase_deg = ((rc["phase_deg"][mc] + 180.0) % 360.0) - 180.0
        low_coh = rp["gamma2"][mp] < COH_GATE

        ax_mag.plot(fp, np.ma.masked_where(low_coh, prbs_mag_db), color=color,
                    linewidth=1.0, label=f"A={A*1e3:.3f} mrad RMS (PRBS)")
        ax_mag.plot(fp, np.ma.masked_where(~low_coh, prbs_mag_db), color="#bbbbbb",
                    linewidth=0.8, zorder=0)
        ax_mag.plot(fc_, rc["mag_db"][mc], color=color, linewidth=1.0, linestyle="--")

        ax_phase.plot(fp, np.ma.masked_where(low_coh, prbs_phase_deg), color=color,
                      linewidth=0.0, marker=".", markersize=1.5)
        ax_phase.plot(fp, np.ma.masked_where(~low_coh, prbs_phase_deg), color="#bbbbbb",
                      linewidth=0.0, marker=".", markersize=1.2, zorder=0)
        ax_phase.plot(fc_, chirp_phase_deg, color=color, linewidth=0.0,
                      marker=".", markersize=1.5)

        ax_coh.plot(fp, rp["gamma2"][mp], color=color, linewidth=1.0)

    ax_coh.axhline(COH_GATE, color="#333333", linestyle=":", linewidth=1.0)

    ax_mag.set_xscale("log")
    ax_mag.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title("x_n / theta_cmd magnitude -- solid=PRBS (grey where coherence<0.9), dashed=chirp")
    ax_mag.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=7, ncol=2)

    ax_phase.set_xscale("log")
    ax_phase.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_phase.set_ylim(-180.0, 180.0)
    ax_phase.set_yticks([-180, -90, 0, 90, 180])
    ax_phase.set_ylabel("Phase (deg, wrapped)")
    ax_phase.set_title("x_n / theta_cmd phase -- wrapped to (-180, 180] (dots: PRBS + chirp, "
                        "not a continuous unwrap -- see caption)")
    ax_phase.grid(True, which="both", linewidth=0.4, color="#cccccc")

    ax_coh.set_xscale("log")
    ax_coh.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_coh.set_xlabel("Frequency (Hz)")
    ax_coh.set_ylabel("Coherence gamma^2 (PRBS only)")
    ax_coh.set_title("Coherence -- <1 is nonlinear distortion, not measurement noise (noiseless sim)")
    ax_coh.set_ylim(0.0, 1.05)
    ax_coh.grid(True, which="both", linewidth=0.4, color="#cccccc")

    fig.suptitle("Rev 4 LuGre nonlinear-drive broadband ID: chirp + PRBS, "
                  f"{len(amplitudes)} RMS amplitudes")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.97])
    montage_path = ASSET_DIR / "broadband_id_montage.png"
    fig.savefig(montage_path, dpi=110)
    plt.close(fig)
    print(f"\nWrote {montage_path}")
    print(f"Wrote {NPZ_DIR / 'broadband_id_data.npz'}")

    # ---- Summary ----
    print("\nPer-port peak |z| / (Fs-or-Ts/sigma0) (static breakaway, all three ports):")
    for i, A in enumerate(amplitudes):
        rp, rc = prbs_results[i], chirp_results[i]
        print(f"  A={A*1e3:9.4f} mrad RMS:")
        print(f"    PRBS : sb={rp['saturation']['sb']:.3f} nut={rp['saturation']['nut']:.3f} "
              f"way={rp['saturation']['way']:.3f}  (n_segments={rp['n_segments']})")
        print(f"    chirp: sb={rc['saturation']['sb']:.3f} nut={rc['saturation']['nut']:.3f} "
              f"way={rc['saturation']['way']:.3f}")
    print(f"\nTotal wall time: {total_wall:.1f}s across {N_WORKERS} workers "
          f"({sum(r['wall_s'] for r in prbs_results):.1f}s PRBS + "
          f"{sum(r['wall_s'] for r in chirp_results):.1f}s chirp of serial work)")


if __name__ == "__main__":
    main()

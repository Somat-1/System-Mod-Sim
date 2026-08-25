#!/usr/bin/env python3
"""Sanity check requested 2026-08-20: is the frequency shift seen when going
from the frictionless baseline Bode (../build_bode_rev4.py) to the LuGre
sub-branch's linear baseline (run_local_linearization_bode.py) caused by
the friction implementation itself, or by a parameter mismatch between
../model_parameters.json and model_parameters.json (this directory)?

Two independent checks:

1. Direct parameter diff. The 19 keys shared between the two JSON files
   (inertias, rotational-chain k/c, k_brg/c_brg, c_EM, N_r/T_hold/T_d, L,
   masses) are compared value-for-value. k_nut/c_nut exist ONLY in the
   parent (the rigid nut coupling this sub-branch replaces); the eighteen
   sigma0/1/2, Fc/Fs/Tc/Ts, vs parameters exist ONLY here (the LuGre ports).
   Neither set being present in the other file is expected, not a mismatch.

2. Re-run build_bode_rev4's frictionless Bode using THIS file's parameters
   for every shared key, with k_nut/c_nut taken from the parent (since this
   file doesn't have them -- there is nothing else to use). If check 1 finds
   no mismatches, this run is mathematically guaranteed to reproduce
   ../rendered_assets/bode_rev4.png's data exactly; it's generated anyway
   because "re-run it" was the explicit ask, and doing so is what actually
   proves the match to a plot rather than a printed diff.

Both frictionless curves (parent-run and rerun-with-this-file's-params) are
then plotted against run_local_linearization_bode.py's saved LuGre linear
baseline, to show what actually changes when friction replaces the rigid
nut coupling. Structural note, found while wiring this up: the LuGre
baseline's K_eq_way term grounds x_n DIRECTLY (K_eq_way ~= sigma0_way, added
to x_n's own diagonal in build_linearized_matrices) -- a stiffness path
that does not exist at all in the frictionless baseline, where x_n has no
direct-to-ground term and is coupled to the rest of the system only through
the nut (k_nut) cross-terms. Adding a new grounding spring generally raises
the natural frequency of modes involving that DOF, which is a structural
(not parameter-mismatch) explanation for an upward shift.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parents[1]
sys.path.insert(0, str(REV4_DIR / "scripts"))

from build_bode_rev4 import build_matrices, build_state_space, laplace_transfer_function, INPUT_LABELS  # noqa: E402

PARENT_PARAM_FILE = REV4_DIR / "model_parameters.json"
LUGRE_PARAM_FILE = LUGRE_DIR / "model_parameters.json"
PARENT_BODE_DATA = REV4_DIR / "rendered_assets" / "npz" / "bode_rev4_data.npz"
LUGRE_LINEAR_BODE_DATA = LUGRE_DIR / "rendered_assets" / "npz" / "local_linearization_bode_data.npz"
OUT_DIR = LUGRE_DIR / "rendered_assets"


def load_params(path: Path) -> dict[str, float]:
    return json.loads(path.read_text(encoding="utf-8"))["parameters"]


def diff_shared_parameters(p_parent: dict, p_lugre: dict) -> int:
    shared = sorted(set(p_parent) & set(p_lugre))
    only_parent = sorted(set(p_parent) - set(p_lugre))
    only_lugre = sorted(set(p_lugre) - set(p_parent))
    print(f"Shared keys ({len(shared)}): {shared}")
    print(f"Only in parent (expected: rigid nut coupling): {only_parent}")
    print(f"Only in lugre_friction (expected: LuGre port params): {only_lugre}")
    n_mismatch = 0
    for k in shared:
        if p_parent[k] != p_lugre[k]:
            print(f"  MISMATCH {k}: parent={p_parent[k]!r}  lugre={p_lugre[k]!r}")
            n_mismatch += 1
    if n_mismatch:
        print(f"-> {n_mismatch}/{len(shared)} shared parameters DIFFER.")
    else:
        print(f"-> All {len(shared)} shared parameters are identical. No mismatch.")
    return n_mismatch


def run_frictionless_bode(p: dict[str, float], frequencies_hz: np.ndarray):
    M, C, K, B_u = build_matrices(p)
    A, B, C_y = build_state_space(M, C, K, B_u)
    b_col = B[:, INPUT_LABELS.index("theta_cmd")]
    response = np.array([
        laplace_transfer_function(A, b_col, C_y, 1j * 2.0 * np.pi * f) for f in frequencies_hz
    ])
    magnitude_db = 20.0 * np.log10(np.maximum(np.abs(response), 1e-15))
    phase_deg = np.unwrap(np.angle(response)) * 180.0 / np.pi
    return magnitude_db, phase_deg


def main() -> None:
    p_parent = load_params(PARENT_PARAM_FILE)
    p_lugre = load_params(LUGRE_PARAM_FILE)

    print("=== Check 1: direct parameter diff ===")
    n_mismatch = diff_shared_parameters(p_parent, p_lugre)

    print("\n=== Check 2: re-run frictionless Bode with this file's parameters ===")
    # p_lugre has no k_nut/c_nut (removed by design, Sec. "What changed" in
    # README.md) -- build_matrices() needs them, so they come from the
    # parent, which is the only place they exist.
    p_merged = dict(p_lugre)
    p_merged["k_nut"] = p_parent["k_nut"]
    p_merged["c_nut"] = p_parent["c_nut"]

    d_parent = np.load(PARENT_BODE_DATA)
    freq_hz = d_parent["frequencies_hz"]
    mag_parent, phase_parent = d_parent["magnitude_db"], d_parent["phase_deg"]

    mag_rerun, phase_rerun = run_frictionless_bode(p_merged, freq_hz)

    mag_diff = np.max(np.abs(mag_rerun - mag_parent))
    phase_diff = np.max(np.abs(phase_rerun - phase_parent))
    print(f"Max |magnitude_db diff| between parent run and rerun-with-lugre-params: {mag_diff:.3e} dB")
    print(f"Max |phase_deg diff|: {phase_diff:.3e} deg")
    print("-> Rerun reproduces the parent Bode exactly." if mag_diff < 1e-9 and phase_diff < 1e-9
          else "-> Rerun DIFFERS from the parent Bode -- investigate further.")

    d_lugre_lin = np.load(LUGRE_LINEAR_BODE_DATA)
    freq_lugre, mag_lugre, phase_lugre = (
        d_lugre_lin["freq_hz"], d_lugre_lin["magnitude_db"], d_lugre_lin["phase_deg"],
    )

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True)
    ax_mag.plot(freq_hz, mag_parent, color="#2b6cb0", linewidth=1.8, alpha=0.6,
                label="frictionless (parent params)")
    ax_mag.plot(freq_hz, mag_rerun, color="#2b6cb0", linewidth=1.0, linestyle=":",
                label="frictionless (rerun, lugre_friction shared params + parent k_nut/c_nut)")
    ax_mag.plot(freq_lugre, mag_lugre, color="#c05621", linewidth=1.2,
                label="LuGre linear baseline (frozen bristles, K_eq/C_eq)")
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title("Frictionless vs. LuGre-linearized command-to-stage Bode -- parameter-match sanity check")
    ax_mag.grid(True, linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=8)

    ax_phase.plot(freq_hz, phase_parent, color="#2b6cb0", linewidth=1.8, alpha=0.6)
    ax_phase.plot(freq_hz, phase_rerun, color="#2b6cb0", linewidth=1.0, linestyle=":")
    ax_phase.plot(freq_lugre, phase_lugre, color="#c05621", linewidth=1.2)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_xlim(0.0, 2000.0)
    ax_phase.grid(True, linewidth=0.4, color="#cccccc")

    fig.tight_layout()
    out_path = OUT_DIR / "frictionless_bode_match_check.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"\nWrote {out_path}")

    print("\n=== Verdict ===")
    if n_mismatch == 0 and mag_diff < 1e-9:
        print("No parameter mismatch: the frictionless Bode is byte-identical whether built from")
        print("../model_parameters.json or this directory's shared parameters + parent k_nut/c_nut.")
        print("Any frequency shift visible against the orange (LuGre linear baseline) curve is")
        print("therefore structural, not a parameter error. One concrete mechanism found while")
        print("building this check: the LuGre baseline's K_eq_way term (~sigma0_way) grounds x_n")
        print("DIRECTLY -- a stiffness path with no counterpart in the frictionless model, where x_n")
        print("has no direct-to-ground term at all and couples to the rest of the system only")
        print("through the (now-removed) k_nut cross-terms. Adding a new grounding spring raises the")
        print("natural frequency of modes that move x_n, which is a plausible source of an upward shift.")
    else:
        print("Mismatch(es) found -- see Check 1/2 output above before trusting any frequency-shift claim.")


if __name__ == "__main__":
    main()

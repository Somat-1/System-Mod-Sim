#!/usr/bin/env python3
"""Plot the raw hardware stepping-test measurements in
"Microstepping Test Data/" (project root, sibling to Rev 4) --
IDSdata.txt (IDS interferometer stage displacement, ~88 ms sample time, pm).
IDSdata.txt is not a single run -- it concatenates five separate runs, each
re-starting its own clock at t=0, one per step size (1, 2, 4, 8, 16), each
introduced by its own "Date: ... Step size N" header line. StepSize8 was
dropped (2026-08-27): its encoder CSV was flagged as faulty and removed from
the data folder. Only the remaining step sizes with a matching encoder CSV
(1, 2, 16) are plotted here.

2026-08-27: the encoder "Counter value" panels (Beckhoff EL5101, 10 ms
sample time, UINT32) were dropped from this overview at request -- only the
IDS interferometer displacement is shown now. `parse_encoder_csv` is kept
(not deleted) because generate_full_step_model_comparison.py in this same
folder still imports it for edge detection.

Displacement sign is flipped from the raw IDS convention (multiplied by -1)
so a commanded step reads as a positive-going move here, matching the sign
convention used in generate_full_step_model_comparison.py's model curves --
otherwise this plot and that one point opposite directions for the same
physical motion.

This data was previously unused by any script in the project -- nothing here
references or post-processes it. This script only visualizes it as measured;
it does not fit, filter, or compare it against any model (see
generate_full_step_model_comparison.py in this same folder for the
model-vs-model comparison at the matching full-step size).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # .../Sytem Mod & Sim
DATA_DIR = PROJECT_ROOT / "Microstepping Test Data"
OUT_DIR = DATA_DIR / "rendered_assets"

STEP_SIZES = [1, 2, 16]
COLORS = {1: "#2b6cb0", 2: "#c05621", 16: "#805ad5"}


def parse_encoder_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Returns (t_s, counts_delta) -- counts_delta is the UINT32 counter
    reinterpreted as a signed rollover and referenced to the first sample."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    data_start = None
    for i, line in enumerate(lines):
        if line.strip() == "" and i > 0:
            data_start = i + 1
    t_ms, counts = [], []
    for line in lines[data_start:]:
        parts = line.strip().split("\t")
        if len(parts) == 2 and parts[0].isdigit():
            t_ms.append(int(parts[0]))
            counts.append(int(parts[1]))
    t = np.array(t_ms, dtype=np.float64) / 1000.0
    c = np.array(counts, dtype=np.int64)
    c_signed = np.where(c > 2**31, c - 2**32, c)
    return t, c_signed - c_signed[0]


def parse_ids_displacement(path: Path) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """IDSdata.txt is not one run -- it concatenates five separate runs, each
    re-starting at t=0, one per step size (1, 2, 4, 8, 16), each introduced by
    its own "Date: ... Step size N" header line followed by its own column
    header line. Returns {step_size: (t_s, displacement_um)} for Axis1 (the
    only populated displacement axis), each block's time zeroed at its own
    first sample."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    blocks: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    current_step: int | None = None
    t_ms: list[float] = []
    disp_pm: list[float] = []

    def flush() -> None:
        if current_step is not None and t_ms:
            t = np.array(t_ms) / 1000.0
            disp_um = (np.array(disp_pm) - disp_pm[0]) * 1.0e-6
            blocks[current_step] = (t, disp_um)

    for line in lines:
        if line.startswith("Date:"):
            flush()
            step_token = line.rsplit(" ", 1)[-1]
            current_step = int(step_token)
            t_ms, disp_pm = [], []
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        try:
            t_val = float(parts[0])
            d_val = float(parts[1])
        except ValueError:
            continue  # column header line
        t_ms.append(t_val)
        disp_pm.append(d_val)
    flush()
    return blocks


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    ids_blocks = parse_ids_displacement(DATA_DIR / "IDSdata.txt")
    extra_ids = sorted(set(ids_blocks) - set(STEP_SIZES))
    if extra_ids:
        print(f"IDSdata.txt also has step size(s) {extra_ids} with no matching "
              f"encoder CSV -- not plotted here.")

    fig, axes = plt.subplots(1, len(STEP_SIZES), figsize=(10.5, 4.5))

    for i, n in enumerate(STEP_SIZES):
        ax = axes[i]
        t, disp_um = ids_blocks[n]
        ax.plot(t, -disp_um, color=COLORS[n], linewidth=0.8)
        ax.set_title(f"StepSize{n}", fontsize=10)
        ax.set_xlabel("Time (s)")
        if i == 0:
            ax.set_ylabel(r"IDS displacement ($\mu$m)"
                           "\n(relative to first sample, sign-flipped)")
        ax.grid(True, linewidth=0.4, color="#cccccc")

    fig.suptitle(
        "Rev 4 hardware stepping test data (as measured, not model output)\n"
        "IDS interferometer stage displacement, ~88 ms/sample -- "
        "Microstepping Test Data/",
        fontsize=11,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.85])

    out_path = OUT_DIR / "microstepping_test_data_overview.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()

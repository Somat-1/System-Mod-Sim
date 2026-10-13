#!/usr/bin/env python3
"""Render the generated motor-identification plan as PNG and SVG."""

from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sequence_design import build_sequence, timeline, summary

COLORS = {
    "warmup": "#8093a1", "alignment": "#6f5cc2", "static_map": "#2468a2",
    "current_sweep": "#d78526", "dynamics": "#b54545", "velocity": "#2b8a6e",
    "reversal": "#9b59a3", "close": "#44515a",
}


def render(out: Path, current_ma: int, stair_ms: int,
           include_reversal: bool) -> None:
    actions = build_sequence(current_ma, stair_ms, include_reversal)
    rows, info = timeline(actions), summary(actions)
    fig, axes = plt.subplots(3, 1, figsize=(16, 9), sharex=True,
                             gridspec_kw={"height_ratios": [1.35, 0.8, 1.0]})
    ax_pos, ax_current, ax_blocks = axes

    for row in rows:
        x0, x1 = row["t_start_s"] / 60, row["t_end_s"] / 60
        color = COLORS[row["block"]]
        for ax in axes:
            ax.axvspan(x0, x1, color=color, alpha=0.075, linewidth=0)
        y0, y1 = row["start_u256"] / 256, row["end_u256"] / 256
        ax_pos.plot([x0, x1], [y0, y1], color=color, linewidth=1.5)
        ax_current.plot([x0, x1], [row["current_ma"], row["current_ma"]],
                        color=color, linewidth=2.2)

    # One labelled bar per contiguous block.
    grouped = []
    for row in rows:
        if not grouped or grouped[-1][0] != row["block"]:
            grouped.append([row["block"], row["t_start_s"], row["t_end_s"]])
        else:
            grouped[-1][2] = row["t_end_s"]
    for block, start, end in grouped:
        x0, width = start / 60, (end - start) / 60
        ax_blocks.barh(0, width, left=x0, height=0.5,
                       color=COLORS[block], edgecolor="white")
        if width > 0.7:
            ax_blocks.text(x0 + width / 2, 0, block.replace("_", "\n"),
                           ha="center", va="center", color="white",
                           fontsize=9, fontweight="bold")

    ax_pos.set_ylabel("Command phase\n[full steps]")
    ax_current.set_ylabel("RMS current\n[mA]")
    ax_blocks.set_yticks([])
    ax_blocks.set_ylim(-0.55, 0.55)
    ax_blocks.set_xlabel("Planned elapsed time [min]")
    ax_pos.grid(axis="y", alpha=0.25)
    ax_current.grid(axis="y", alpha=0.25)
    ax_pos.set_title("Motor identification - native MRES=16, SpreadCycle, "
                     "MicroPlyer OFF, StealthChop OFF", loc="left",
                     fontsize=15, fontweight="bold")
    fig.text(0.01, 0.925,
             f"{info['total_minutes']:.1f} min | {info['actions']} actions | "
             "shaft encoder up to 2 kHz | full-revolution staircase out and back",
             fontsize=10)
    handles = [Patch(facecolor=color, label=name.replace("_", " "))
               for name, color in COLORS.items() if name in info["blocks"]]
    ax_pos.legend(handles=handles, ncol=4, frameon=False, loc="upper right")
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    out.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        target = out / f"motor_identification_sequence.{suffix}"
        fig.savefig(target, dpi=180 if suffix == "png" else None,
                    bbox_inches="tight")
        print(f"wrote {target}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nominal-current-ma", type=int, default=360)
    parser.add_argument("--stair-dwell-ms", type=int, default=150)
    parser.add_argument("--no-reversal", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    out = args.out or Path(__file__).resolve().parent.parent / "sequence"
    render(out, args.nominal_current_ma, args.stair_dwell_ms, not args.no_reversal)


if __name__ == "__main__":
    main()

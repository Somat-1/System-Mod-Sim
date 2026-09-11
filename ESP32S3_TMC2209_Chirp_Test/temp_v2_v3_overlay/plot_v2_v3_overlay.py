"""TEMPORARY: overlay the v2 and v3 chirp Bode magnitude for the primary axis.

Both runs used the same firmware, schedule and excitation parameters (cruise
3.0 full steps, notch 176.7 Hz +- 56.7 Hz, depth 0.10, MRES 1/16, SpreadCycle,
MicroPlyer off), so raw magnitude is directly comparable without normalising.

Channel: column 0 of each capture. v2 calls it "AI 1", v3 calls it "X". The
.npz keys are hardcoded to the column index in both scripts, so `up_mag_AI1`
is column 0 regardless of the display label.

CAVEAT: that the two column-0 channels are the same physical accelerometer
axis is an assumption from the rig being unchanged between runs, not something
either recording states.

Encoding: colour carries the run, line style carries the sweep direction.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
V2 = ROOT / "v2" / "analysis" / "plots" / "bode_data.npz"
V3 = ROOT / "v3" / "analysis" / "plots" / "bode_data.npz"

C_V2 = "#1769aa"
C_V3 = "#d94801"
C_INK = "#1a1a1a"
C_MUTED = "#6b7280"
C_GRID = "#d7dbe0"

F_N_HZ = 176.7
NOTCH_HALFWIDTH_HZ = 56.7

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": C_GRID, "axes.labelcolor": C_INK,
    "axes.titlesize": 11, "axes.titleweight": "semibold", "axes.labelsize": 9,
    "axes.grid": True, "grid.color": C_GRID, "grid.linewidth": 0.6,
    "grid.alpha": 0.7, "xtick.color": C_MUTED, "ytick.color": C_MUTED,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.frameon": False,
    "legend.fontsize": 8.5, "font.size": 9, "text.color": C_INK,
})

d2, d3 = np.load(V2), np.load(V3)
f = d2["freq_grid"]
assert np.allclose(f, d3["freq_grid"]), "frequency grids differ"

up2, dn2, fl2 = d2["up_mag_AI1"], d2["down_mag_AI1"], d2["floor_AI1"]
up3, dn3, fl3 = d3["up_mag_AI1"], d3["down_mag_AI1"], d3["floor_AI1"]

fig, (ax, axr) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                              gridspec_kw={"height_ratios": [2.2, 1]})

lo, hi = F_N_HZ - NOTCH_HALFWIDTH_HZ, F_N_HZ + NOTCH_HALFWIDTH_HZ
for a in (ax, axr):
    a.axvspan(lo, hi, color=C_MUTED, alpha=0.10, lw=0, zorder=0)

ax.plot(f, up2, color=C_V2, lw=1.5, label="v2  up-sweep")
ax.plot(f, dn2, color=C_V2, lw=1.2, ls="--", alpha=0.75, label="v2  down-sweep")
ax.plot(f, up3, color=C_V3, lw=1.5, label="v3  up-sweep")
ax.plot(f, dn3, color=C_V3, lw=1.2, ls="--", alpha=0.75, label="v3  down-sweep")
ax.plot(f, fl2, color=C_V2, lw=0.8, ls=":", alpha=0.5, label="v2  noise floor")
ax.plot(f, fl3, color=C_V3, lw=0.8, ls=":", alpha=0.5, label="v3  noise floor")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("Acceleration magnitude [m/s²]")
ax.set_title("v2 vs v3 chirp — primary axis (v2 “AI 1” / v3 “X”), raw magnitude",
             loc="left", pad=22)
ax.text(0.0, 1.012,
        "same firmware, schedule and excitation; shaded band = commanded notch "
        f"({lo:.0f}–{hi:.0f} Hz, depth 0.10)",
        transform=ax.transAxes, fontsize=8.5, color=C_MUTED, va="bottom")
ax.legend(loc="upper left", ncol=3)

with np.errstate(divide="ignore", invalid="ignore"):
    ratio_up = np.where(up2 > 0, up3 / up2, np.nan)
    ratio_dn = np.where(dn2 > 0, dn3 / dn2, np.nan)
axr.axhline(1.0, color=C_INK, lw=1.0, ls="-", alpha=0.5)
axr.plot(f, ratio_up, color=C_V3, lw=1.4, label="up-sweep")
axr.plot(f, ratio_dn, color=C_V3, lw=1.2, ls="--", alpha=0.75,
         label="down-sweep")
axr.set_xscale("log")
axr.set_yscale("log")
axr.set_ylim(0.2, 5.0)
axr.set_xlabel("Frequency [Hz]")
axr.set_ylabel("v3 / v2")
axr.set_title("Ratio — 1.0 means the two runs agree", loc="left")
axr.legend(loc="upper left", ncol=2)

fig.tight_layout()
out = HERE / "v2_vs_v3_primary_axis.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"wrote {out}")

band = (f >= 20) & (f <= 900) & np.isfinite(ratio_up)
print(f"\nratio v3/v2 over 20-900 Hz (up-sweep): median {np.nanmedian(ratio_up[band]):.3f}, "
      f"10th pct {np.nanpercentile(ratio_up[band], 10):.3f}, "
      f"90th pct {np.nanpercentile(ratio_up[band], 90):.3f}")
for lbl, up, dn in (("v2", up2, dn2), ("v3", up3, dn3)):
    pk = int(np.nanargmax(np.where((f > 100) & (f < 260), up, np.nan)))
    print(f"{lbl}: up-sweep peak in 100-260 Hz at {f[pk]:.1f} Hz "
          f"= {up[pk]:.3f} m/s²")

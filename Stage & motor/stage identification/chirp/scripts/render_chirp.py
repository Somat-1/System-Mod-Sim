#!/usr/bin/env python3
"""Draw both chirp profiles: frequency against time, and what each one covers.

Writes ../sequence/chirp_sweeps.png.

Three panels. The two upper ones are the profiles as they run — commanded
frequency against time, with the excitation amplitude shaded behind. The lower
one puts both on a frequency axis so the thing that actually matters is visible
at a glance: which band each profile touches, and the 150–250 Hz band that
neither of them drives.

Palette is Okabe-Ito; the resonance band is the one place a status colour is
used, because it means "do not drive", not "series 3".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chirp_design as cd  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

INK = {"full": "#0072B2", "high": "#009E73"}
DANGER = "#D55E00"
GRID = "#D7D7D7"
MUTED = "#6B6B6B"
TEXT = "#222222"


def curve(profile):
    """Commanded frequency and amplitude against time, as step/ramp segments."""
    t, f, a = [], [], []
    for t0, t1, s in profile.timeline():
        if s.kind == "silent":
            t += [t0, t1]; f += [np.nan, np.nan]; a += [0.0, 0.0]
        elif s.duration_s == 0.0:
            # the gap jump: a vertical break, drawn as a discontinuity
            t += [t0, t0]; f += [s.f_start, s.f_end]; a += [np.nan, np.nan]
        else:
            # the amplitude follows the law along the sweep, so sample it
            for k in range(41):
                x = k / 40
                fk = s.f_start + (s.f_end - s.f_start) * x
                t.append(t0 + (t1 - t0) * x)
                f.append(fk)
                a.append(cd.amplitude(fk) if s.amplitude_microsteps else 0.0)
    return np.array(t), np.array(f), np.array(a)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = Path(args.out) if args.out else Path(__file__).resolve().parent.parent / "sequence"
    out.mkdir(parents=True, exist_ok=True)

    profs = cd.profiles()
    summary = cd.summarise(profs)
    lo, hi = cd.RESONANCE_GAP

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#999999",
                         "axes.labelcolor": TEXT, "text.color": TEXT,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig = plt.figure(figsize=(14.5, 9.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.15, 0.85], hspace=0.62,
                          left=0.075, right=0.975, top=0.9, bottom=0.085)

    for i, p in enumerate(profs):
        ax = fig.add_subplot(gs[i])
        t, f, a = curve(p)
        # the band nobody drives
        ax.axhspan(lo, hi, color=DANGER, alpha=0.13, lw=0, zorder=0)
        ax.axhline(lo, color=DANGER, lw=0.8, ls="--", alpha=0.75)
        ax.axhline(hi, color=DANGER, lw=0.8, ls="--", alpha=0.75)
        ax.plot(t, f, lw=1.9, color=INK[p.key], zorder=3, solid_capstyle="round")

        # amplitude behind, on its own scale but as a shaded envelope only
        amp = np.nan_to_num(a)
        ax.fill_between(t, 0, amp / cd.CRUISE_MICROSTEPS * 90.0, step=None,
                        color=INK[p.key], alpha=0.10, lw=0, zorder=1)

        for t0, t1, s in p.timeline():
            if s.kind == "marker":
                ax.axvspan(t0, t1, color="#555555", alpha=0.35, lw=0, zorder=2)
            if s.kind == "dwell" and s.duration_s == 0.0:
                ax.annotate("", xy=(t0, s.f_end), xytext=(t0, s.f_start),
                            arrowprops=dict(arrowstyle="->", color=DANGER, lw=1.6),
                            zorder=4)
                ax.text(t0 + 1.5, (s.f_start + s.f_end) / 2,
                        "jump — band not visited", fontsize=7.6, color=DANGER,
                        va="center")
        info = summary["profiles"][p.key]
        ax.set_xlim(0, p.duration_s)
        ax.set_ylim(0, 1080)
        ax.set_ylabel("commanded frequency  [Hz]")
        ax.set_xlabel("time  [s]")
        ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
        ax.set_title(f"{p.title}   ·   {info['duration_min']:.2f} min",
                     fontsize=11.5, color=INK[p.key], weight="bold", loc="left")
        ax.text(0.995, 0.055, p.purpose, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=8, color=MUTED)

    # ---------------------------------------------------- coverage comparison
    ax = fig.add_subplot(gs[2])
    ax.axvspan(lo, hi, color=DANGER, alpha=0.16, lw=0)
    for j, p in enumerate(profs):
        y = len(profs) - 1 - j
        for a, b in summary["profiles"][p.key]["bands_hz"]:
            ax.barh(y, b - a, left=a, height=0.42, color=INK[p.key],
                    edgecolor="white", linewidth=2)
        ax.text(1.5, y + 0.33, p.key.upper(), fontsize=9.5, color=INK[p.key],
                weight="bold", va="bottom")
    ax.text((lo + hi) / 2, len(profs) - 0.5,
            f"{lo:g}–{hi:g} Hz\nnot driven by either",
            ha="center", va="center", fontsize=8.4, color=DANGER, weight="bold",
            linespacing=1.4)
    ax.set_xscale("log")
    ax.set_xlim(1, 1100)
    ax.set_ylim(-0.6, len(profs) - 0.25)
    ax.set_yticks([])
    ax.set_xlabel("frequency  [Hz]")
    ax.grid(True, axis="x", color=GRID, lw=0.6, alpha=0.8)
    ax.set_title("What each profile covers", fontsize=11.5, weight="bold", loc="left")

    fig.suptitle("Mechanical chirp sweeps — stage identification",
                 fontsize=13.5, weight="bold", x=0.075, ha="left", y=0.968)
    fig.text(0.075, 0.935,
             f"sweep rate {cd.SWEEP_RATE_HZ_PER_S:g} Hz/s · amplitude "
             f"{cd.CRUISE_FULL_STEPS:g} full step to {cd.F_CORNER_HZ:g} Hz, then 1/f² "
             f"({cd.CRUISE_MICROSTEPS:g} µsteps "
             f"at MRES {cd.MRES}) · {cd.CURRENT_RMS_MA} mA RMS · "
             f"{cd.SYNC_FREQ_HZ:g} Hz sync bursts at each end",
             fontsize=8.6, color=MUTED)

    handles = [Line2D([], [], color=INK[p.key], lw=2.4, label=p.title) for p in profs]
    handles.append(Patch(facecolor=DANGER, alpha=0.25,
                         label=f"{lo:g}–{hi:g} Hz — uncontrolled motion 2026-09-15; not driven"))
    handles.append(Patch(facecolor="#555555", alpha=0.35, label="sync marker burst"))
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.975, 0.975),
               frameon=False, fontsize=8.3)

    png = out / "chirp_sweeps.png"
    fig.savefig(png, dpi=170)
    plt.close(fig)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()

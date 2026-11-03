#!/usr/bin/env python3
"""Render the E5 sequence, the command correction and the expected result.

    F1  rendered/e5_sequence_windows.png    the session on one time axis, each
                                            window labelled with what it validates
    F2  rendered/e5_command_correction.png  what C runs command at each ISO target
    F3  rendered/e5_expected_result.png     U vs C as the frozen model predicts,
                                            against the pass lines of §6

The expected C result is an offline replay: E0 v4's measured deviations plus
the correction the table commands at each target and direction (quantisation
included). It is a prediction to hold E5 against, not a result.

    DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib \
      /opt/homebrew/bin/python3 render_e5_sequence.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered"
sys.path.insert(0, str(HERE))
from e5_geometry import FrozenModel  # noqa: E402
from e5_make_table import move_time_s  # noqa: E402

V4_RESULTS = HERE.parents[2] / "iso230-2" / "analysis" / "iso230_2_v4_results_blind.json"

INK, INK2, GRID = "#252b30", "#59636b", "#dfe3e6"
U_COL, C_COL, B3_COL, REF_COL = "#b5533a", "#2d6fa6", "#b8862b", "#7c5a9b"
PASS_COL = "#2f6b4f"
GAP_MIN = 3.0                      # drawn width of the operator pause between recordings

# Acceptance criteria, as fixed in the procedure (§6).
OFFLINE_M, OFFLINE_E, MARGIN = 1.335, 1.560, 0.4


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": "#7d858c", "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c", "ytick.color": "#41474c",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "legend.fontsize": 7.8, "axes.unicode_minus": True,
    })


def finish(ax) -> None:
    ax.grid(color=GRID, lw=.45)
    ax.spines[["top", "right"]].set_visible(False)


def save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)


# ---- data ------------------------------------------------------------------
def num(v: str):
    if v == "":
        return None
    try:
        return int(v)
    except ValueError:
        return float(v)


def load():
    config = json.loads((HERE / "e5_config.json").read_text())
    summary = json.loads((HERE / "sequence" / "e5_plan_summary.json").read_text())
    text = ("state", "action", "block", "mode", "note")
    with (HERE / "sequence" / "e5_measurement_plan.csv").open(newline="") as fh:
        plan = [{k: (v if k in text else num(v)) for k, v in r.items()}
                for r in csv.DictReader(fh)]
    with (HERE / "sequence" / "e5_command_table.csv").open(newline="") as fh:
        table = [{k: (v if k in ("block", "mode") else num(v)) for k, v in r.items()}
                 for r in csv.DictReader(fh)]
    model = FrozenModel((HERE / config["model_file"]).resolve())
    assert model.sha256 == summary["model_sha256"], "model changed since the table was made"
    return config, summary, plan, table, model


def walk(plan: list[dict], c: dict):
    """Trace points (t_min, x_mm) and validation windows on one continuous clock.

    Recording 2 starts GAP_MIN after recording 1 ends, so the operator pause
    is visible. A window is a contiguous stretch with one purpose: a reference
    group, one ISO run (with E0 v4's own x_ref visit at its start), or one
    B3 repetition.
    """
    trace, windows = [], []
    t, pos, rec = 0.0, c["x_ref_mm"], 1
    purpose, start = None, 0.0

    def close(now):
        if purpose is not None and now > start:
            windows.append((purpose, start, now))

    for r in plan:
        if r["recording"] != rec:
            close(t); purpose = None
            rec, t = r["recording"], t + GAP_MIN
            windows.append((("gap",), t - GAP_MIN, t))
        if r["action"] == "state_enter":
            state = r["state"]
            if state in ("DATUM", "MID_REF", "END_REF", "B3_REF", "B3_END_REF"):
                new = ("ref", r["block"])
            elif state == "ISO_REF":
                new = ("iso", c["b2_mode_order"][r["run"] - 1], r["run"])
            elif state == "ISO":
                new = purpose
            elif state == "B3":
                new = ("b3", r["run"])
            else:
                new = None
            if new != purpose:
                close(t); purpose, start = new, t
            continue
        if r["u_cmd"] is None:
            continue
        x = r["u_cmd"] / c["ustep_per_mm"]
        t1 = t + move_time_s(x - pos, c) / 60
        trace.append((t, pos)); trace.append((t1, x))
        t, pos = t1, x
        if r["action"] == "stop":
            t += (c["settle_s"] + c["record_s"]) / 60
            trace.append((t, x))
    close(t)
    return trace, windows, t


# ---- F1: sequence with validation windows -----------------------------------
def figure_sequence(config, summary, plan) -> None:
    trace, windows, end = walk(plan, config)
    rec_min = {int(k): v for k, v in summary["recording_min"].items()}
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    fig.subplots_adjust(left=.06, right=.93, top=.66, bottom=.17)

    def colour_of(p):
        if p[0] == "iso":
            return U_COL if p[1] == "U" else C_COL
        return {"b3": B3_COL, "ref": REF_COL}.get(p[0], INK2)

    ts = np.array([p[0] for p in trace]); xs = np.array([p[1] for p in trace])
    for p, a, b in windows:
        if p[0] == "gap":
            continue
        m = (ts >= a - 1e-9) & (ts <= b + 1e-9)
        ax.plot(ts[m], xs[m], color=colour_of(p), lw=.75, solid_joinstyle="round")

    for p, a, b in windows:
        mid = .5 * (a + b)
        if p[0] in ("iso", "b3"):
            label = f"{p[1]}{p[2]}" if p[0] == "iso" else f"B3·{p[1]}"
            ax.text(mid, 126, label, ha="center", va="center", fontsize=8, color="white",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=.25", fc=colour_of(p), ec="none"))
        elif p[0] == "ref":
            ax.axvspan(a, b, color=REF_COL, alpha=.22, lw=0)
            ax.annotate(p[1], xy=(mid, 0), xytext=(mid, -12), ha="center", va="top",
                        fontsize=7.6, color=REF_COL, fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color=REF_COL, lw=.8),
                        annotation_clip=False)
        elif p[0] == "gap":
            ax.axvspan(a, b, color="#eef0f2", lw=0)
            ax.text(mid, 62, "pause: save IDS file 1, start file 2 (no re-home)",
                    ha="center", va="center", fontsize=6.8, color=INK2, rotation=90)

    def bracket(a, b, y, colour, head, body=""):
        tr = ax.get_xaxis_transform()
        ax.plot([a, a, b, b], [y - .03, y, y, y - .03], color=colour, lw=1.2,
                transform=tr, clip_on=False)
        ax.text(.5 * (a + b), y + .015, head, transform=tr, ha="center", va="bottom",
                fontsize=8.6, color=colour, fontweight="bold", clip_on=False)
        if body:
            ax.text(.5 * (a + b), y + .075, body, transform=tr, ha="center", va="bottom",
                    fontsize=7.4, color=INK, clip_on=False, linespacing=1.3)

    b2 = [(a, b) for p, a, b in windows if p[0] == "iso"]
    b3 = [(a, b) for p, a, b in windows if p[0] == "b3"]
    bracket(b2[0][0], b2[-1][1], 1.04, INK,
            "B2   ISO 230-2, paired U vs C   (5 U + 5 C, order ABBA)",
            "§5.2  ISO metrics, U vs C     §5.3  achieved vs commanded correction (gain k)\n"
            "§5.4  C residual vs screw angle θ and x, reversal B")
    bracket(b3[0][0], b3[-1][1], 1.04, B3_COL, "B3   off-grid",
            "§5.5  residual vs\nelectrical phase")
    bracket(0, rec_min[1], 1.36, INK2, f"IDS recording 1   ({rec_min[1]:.1f} min)")
    bracket(rec_min[1] + GAP_MIN, end, 1.36, INK2, f"recording 2   ({rec_min[2]:.1f} min)")

    tr = ax.get_xaxis_transform()
    ax.text(end + .8, -.03, "B1a–e\ndrift anchors\n(§5.1, D2)", transform=tr, ha="left",
            va="top", fontsize=7.4, color=REF_COL, fontweight="bold", clip_on=False)
    ax.text(end + .8, .5, "§5.6 settling:\nread at every\nstop from the\ncontinuous\nIDS stream",
            transform=tr, ha="left", va="center", fontsize=7.2, color=INK2, clip_on=False)

    ax.axhline(config["x_ref_mm"], color=REF_COL, lw=.5, ls=(0, (2, 3)), alpha=.6)
    ax.set_xlim(-.2, end + .2); ax.set_ylim(0, 132)
    ax.set_yticks([0, 30, 61.35, 90, 120], ["0", "30", "x_ref", "90", "120"])
    ax.set_xlabel("session time [min]", labelpad=16)
    ax.set_ylabel("position [mm]")
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("E5 session: what each window validates", x=.06, ha="left",
                 fontsize=12, color=INK, y=.985)
    save(fig, "e5_sequence_windows")


# ---- F2: command correction -------------------------------------------------
def figure_correction(config, summary, table, model) -> None:
    rows = [r for r in table if r["block"] == "B2" and r["mode"] == "C" and r["run"] == 2]
    up = sorted((r for r in rows if r["dir"] == 1), key=lambda r: r["x_nom_mm"])
    dn = sorted((r for r in rows if r["dir"] == -1), key=lambda r: r["x_nom_mm"])
    x = np.array([r["x_nom_mm"] for r in up])
    cu = np.array([r["corr_um"] for r in up]); cd = np.array([r["corr_um"] for r in dn])
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    fig.subplots_adjust(left=.1, right=.97, top=.82, bottom=.14)
    ax.axhline(0, color=U_COL, lw=1.4, label="U runs: commanded = nominal")
    ax.plot(x, cu, "-o", color=C_COL, ms=4, lw=1.1, label="C, approached from below (up)")
    ax.plot(x, cd, "--o", color=C_COL, mfc="white", ms=4, lw=1.1,
            label="C, approached from above (down)")
    i = int(np.argmax(np.abs(np.r_[cu, cd])))
    xi, ci = np.r_[x, x][i], np.r_[cu, cd][i]
    ax.annotate(f"largest {ci:+.2f} µm", xy=(xi, ci), xytext=(xi - 28, ci + .3),
                fontsize=7.6, color=INK, va="center",
                arrowprops=dict(arrowstyle="-", color=INK2, lw=.7))
    j = 9
    ax.annotate("", xy=(x[j], cu[j]), xytext=(x[j], cd[j]),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=.9, shrinkA=3, shrinkB=3))
    ax.annotate(f"up − down = {-model.beta0_um:.2f} µm\n= reversal correction −β₀",
                xy=(x[j], .5 * (cu[j] + cd[j])), xytext=(x[j] + 4, -3.6),
                fontsize=7.4, va="center", color=INK,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=.7))
    ax.set_xlim(0, 122.8); ax.set_ylim(-5, 5.5)
    ax.set_xlabel("ISO target x* [mm]")
    ax.set_ylabel("commanded correction [µm]")
    ax.legend(loc="upper left", frameon=False)
    finish(ax)
    fig.suptitle("E5 command: C runs pre-correct each ISO target with the frozen model",
                 x=.1, ha="left", fontsize=11, color=INK, y=.97)
    fig.text(.1, .875, "correction = −(Ĝ(x*) − a₀) − (d/2)·β₀, rounded to whole µsteps "
             f"(0.625 µm)     model {summary['model_sha256'][:12]}", fontsize=7.3, color=INK2)
    save(fig, "e5_command_correction")


# ---- F3: expected result ----------------------------------------------------
def v4_baseline():
    d = json.loads(V4_RESULTS.read_text())["per_target"]
    rows = [d[k] for k in sorted(d, key=int)]
    get = lambda key: np.array([r[key] for r in rows])
    return {"x": get("cmd_mm"), "up": get("xbar_up"), "dn": get("xbar_dn"),
            "s_up": get("s_up"), "s_dn": get("s_dn")}


def iso_metrics(up, dn, s_up, s_dn) -> dict:
    mid = 0.5 * (up + dn)
    return {"E": float(max(up.max(), dn.max()) - min(up.min(), dn.min())),
            "M": float(np.ptp(mid)),
            "A": float(np.max(np.maximum(up + 2*s_up, dn + 2*s_dn)) -
                       np.min(np.minimum(up - 2*s_up, dn - 2*s_dn))),
            "B": float(np.max(np.abs(up - dn))), "B_mean": float(np.mean(up - dn)),
            "R_up": float(4*s_up.max()), "R_dn": float(4*s_dn.max()),
            "R": float(max(np.max(2*s_up + 2*s_dn + np.abs(up - dn)),
                           4*s_up.max(), 4*s_dn.max()))}


def replay(table: list[dict], base: dict) -> dict:
    """E0 v4 deviations plus the correction each C row of run 2 commands."""
    rows = [r for r in table if r["block"] == "B2" and r["mode"] == "C" and r["run"] == 2]
    pick = lambda d: np.array([r["corr_um"] for r in sorted(
        (r for r in rows if r["dir"] == d), key=lambda r: r["target_id"])])
    up, dn = base["up"] + pick(1), base["dn"] + pick(-1)
    return {"up": up, "dn": dn,
            "U": iso_metrics(base["up"], base["dn"], base["s_up"], base["s_dn"]),
            "C": iso_metrics(up, dn, base["s_up"], base["s_dn"])}


def figure_result(config, summary, table, model) -> dict:
    base = v4_baseline(); rp = replay(table, base)
    U, C = rp["U"], rp["C"]
    fig = plt.figure(figsize=(9.6, 4.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=.25,
                          left=.07, right=.98, top=.79, bottom=.14)
    fig.suptitle("E5 expected result: the frozen model replayed on the E0 v4 baseline",
                 x=.07, ha="left", fontsize=11, color=INK, y=.97)
    fig.text(.07, .875, "C = E0 v4 measured deviation + the correction this table commands. "
             "A prediction to test E5 against, not a measurement.", fontsize=7.4, color=INK2)

    ax = fig.add_subplot(gs[0])
    x = base["x"]
    for up, dn, col, lab in ((base["up"], base["dn"], U_COL, "U (E0 v4)"),
                             (rp["up"], rp["dn"], C_COL, "C (expected)")):
        ax.fill_between(x, np.minimum(up, dn), np.maximum(up, dn), color=col, alpha=.18, lw=0)
        ax.plot(x, .5 * (up + dn), color=col, lw=1.6, label=lab)
    ax.axhline(0, color=INK2, lw=.5)
    ax.text(121, 3.9, "line = mean of up and down;  band = up-to-down gap",
            ha="right", fontsize=7.2, color=INK2)
    ax.set_xlim(0, 122.8); ax.set_ylim(-5.5, 4.6)
    ax.set_xlabel("ISO target x* [mm]"); ax.set_ylabel("deviation from x* [µm]")
    ax.legend(loc="lower left", frameon=False)
    finish(ax)
    ax.set_title("Deviation along the travel", loc="left")

    ax = fig.add_subplot(gs[1])
    keys = [("E", "E"), ("M", "M"), ("|mean B|", "B_mean")]
    idx = np.arange(len(keys))
    uvals = [abs(U[k]) for _, k in keys]; cvals = [abs(C[k]) for _, k in keys]
    ax.bar(idx - .2, uvals, width=.38, color=U_COL, label="U (E0 v4)")
    ax.bar(idx + .2, cvals, width=.38, color=C_COL, label="C (expected)")
    passes = [OFFLINE_E + MARGIN, OFFLINE_M + MARGIN, .3]
    for i, (v, c) in enumerate(zip(passes, cvals)):
        ax.plot([i - .45, i + .45], [v, v], color=INK, lw=1.4,
                label="pass line (§6)" if i == 0 else None)
        ok = c <= v
        ax.text(i + .2, max(c, v) + .2, f"{c:.2f}\n{'pass' if ok else 'fail'}", ha="center",
                va="bottom", fontsize=7.4, color=PASS_COL if ok else U_COL, fontweight="bold")
    for i, v in enumerate(uvals):
        ax.text(i - .2, v + .2, f"{v:.2f}", ha="center", va="bottom", fontsize=7.2, color=INK)
    ax.set_xticks(idx, [k for k, _ in keys]); ax.set_ylim(0, 10.5)
    ax.set_ylabel("µm")
    ax.legend(loc="upper right", frameon=False)
    finish(ax)
    ax.set_title("ISO 230-2 headline metrics", loc="left")
    save(fig, "e5_expected_result")
    return rp


def main() -> int:
    style()
    config, summary, plan, table, model = load()
    for stale in ("e5_session_sequence", "e5_command_table_checks", "e5_validation_plan"):
        (OUT / f"{stale}.png").unlink(missing_ok=True)
    figure_sequence(config, summary, plan)
    figure_correction(config, summary, table, model)
    rp = figure_result(config, summary, table, model)
    out = {"U_E0_v4": rp["U"], "C_expected": rp["C"],
           "note": "offline replay: E0 v4 deviations + commanded correction; R unchanged by construction"}
    (HERE / "sequence" / "e5_expected_metrics.json").write_text(json.dumps(out, indent=2) + "\n")
    print("U (E0 v4): " + "  ".join(f"{k} {rp['U'][k]:.3f}" for k in ("E", "M", "A", "B_mean")))
    print("C expect : " + "  ".join(f"{k} {rp['C'][k]:.3f}" for k in ("E", "M", "A", "B_mean")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

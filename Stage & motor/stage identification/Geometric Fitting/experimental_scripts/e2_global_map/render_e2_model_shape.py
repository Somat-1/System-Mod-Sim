"""Two views of the fitted E2 geometric model, for MODEL.md.

1. An exploded component view: one strip per term, all on the same micrometre
   scale, so the size of each part along the travel is readable at a glance.
2. The (x, theta) map of the periodic part - the compensation table itself -
   with per-turn slices at three positions.

The model is rebuilt from the configuration render_e2_findings.py selected
(H=4, Ks=2, Kc=8, lambda=1.0), so the costly cross-validation sweep is not
repeated and the figures show the same fit as the E2 findings.

Writes rendered/model/e2_model_exploded_components.png
       rendered/model/e2_model_xtheta_table.png
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import render_e2_findings as m                                      # noqa: E402

OUT = HERE / "rendered" / "model"
INK, INK2, GRID = "#252b30", "#59636b", "#dfe3e6"
MINUS = "−"
ZOOM = (50.0, 54.0)          # two screw turns, for the detail column
SLICES = (10.0, 50.0, 78.0)  # where the per-turn shape is drawn

COLOURS = {"gain": "#99483e", "slow": "#356449", "k1": "#5f8bc5",
           "k2": "#c88e5c", "k34": "#7c6aa8", "sum": "#252b30"}
DIVERGING = LinearSegmentedColormap.from_list(
    "e2_table", ["#2f5f99", "#8fb0d6", "#eceeef", "#d8a877", "#99483e"])


def num(value, dp=2, sign=False):
    text = ("%+." + str(dp) + "f" if sign else "%." + str(dp) + "f") % value
    return text.replace("-", MINUS)


def tidy(ax):
    ax.grid(color=GRID, lw=.45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)


def model_fit():
    """The selected model, rebuilt from the configuration it was chosen with."""
    run = m.e2_run.load()
    data = m.preprocess(run, m.record_quality(run))
    return data, m.fit_selected(data, {"H": 4, "K_s": 2, "K_C": 8,
                                       "lambda": 1.0})


def harmonic_term(model, x, theta, k):
    """A_k(x) cos k0 + B_k(x) sin k0, exactly as the design matrix builds it."""
    alpha, beta = model["state"]["groups"][f"C{k}"]
    if model["state"]["Kc"] is not None and k <= 2:
        basis = m.raw_bspline(x, model["state"]["Kc"], model["state"]["domain"])
    else:
        basis = np.ones((x.size, 1))
    a = basis @ model["beta"][alpha]
    b = basis @ model["beta"][beta]
    return a * np.cos(k * theta) + b * np.sin(k * theta)


def amplitude(model, x, k):
    return m.harmonic_profile(model, x, k)["amplitude"]


def parts(model, x):
    """Every term of the model at the commanded angle that belongs to x."""
    theta = m.screw_theta(x)
    comp = m.components(model, x)
    k1 = harmonic_term(model, x, theta, 1)
    k2 = harmonic_term(model, x, theta, 2)
    k34 = sum(harmonic_term(model, x, theta, k)
              for k in range(3, model["state"]["H"] + 1))
    total = m.predict(model, x, theta) - model["beta"][0]
    return {"gain": comp["gain"], "slow": comp["slow"], "k1": k1, "k2": k2,
            "k34": k34, "sum": total, "theta": theta}


def figure_exploded(model):
    rows = [
        ("gain", "gain  $a_1(x-\\bar{x})$", None),
        ("slow", "slow trend  $s(x)$", None),
        ("k1", "once per turn  $k=1$", 1),
        ("k2", "twice per turn  $k=2$", 2),
        ("k34", "$k=3$ and $k=4$", None),
        ("sum", "model  $\\hat{G}(x,\\theta)-a_0$", None),
    ]
    lo, hi = model["state"]["domain"]
    x = np.linspace(lo, hi, 24000)
    z = np.linspace(*ZOOM, 4000)
    full, zoom = parts(model, x), parts(model, z)
    limit = 1.08 * np.max(np.abs(full["sum"]))

    fig = plt.figure(figsize=(9.9, 8.2))
    gs = fig.add_gridspec(len(rows), 2, width_ratios=[3.25, 1.0], left=.085,
                          right=.975, top=.862, bottom=.085, hspace=.46,
                          wspace=.07)
    for i, (key, label, k) in enumerate(rows):
        ax = fig.add_subplot(gs[i, 0])
        bx = fig.add_subplot(gs[i, 1], sharey=ax)
        colour = COLOURS[key]
        if k is not None:
            env = amplitude(model, x, k)
            ax.fill_between(x, -env, env, color=colour, alpha=.22, lw=0)
            bx.fill_between(z, -amplitude(model, z, k),
                            amplitude(model, z, k), color=colour, alpha=.22, lw=0)
        ax.axhline(0, color="#aeb4b9", lw=.6)
        bx.axhline(0, color="#aeb4b9", lw=.6)
        ax.plot(x, full[key], color=colour, lw=.55)
        bx.plot(z, zoom[key], color=colour, lw=1.1)
        ax.set_ylim(-limit, limit)
        ax.set_xlim(lo, hi)
        bx.set_xlim(*ZOOM)
        ax.set_ylabel("[µm]")
        # Both captions live outside the axes, so no strip can cover its own
        # label however large the term is.
        ax.set_title(label, loc="left", pad=3, fontsize=8.4, color=INK)
        span = np.ptp(full[key])
        bx.set_title(f"peak-to-peak {num(span)} µm", loc="right", pad=3,
                     fontsize=7.4, color=INK2)
        for a in (ax, bx):
            tidy(a)
        plt.setp(bx.get_yticklabels(), visible=False)
        if i < len(rows) - 1:
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(bx.get_xticklabels(), visible=False)
    ax.set_xlabel("commanded position x [mm]")
    bx.set_xlabel("x [mm]")

    fig.text(.085, .962, "What the E2 model is made of, and how big each part "
             "is along the travel", fontsize=12.5, color="#20262b", va="top")
    fig.text(.085, .930,
             "Left: the whole travel, every strip on the same µm scale - the "
             "shaded envelope is $\\pm A_k(x)$, the line is the term at the "
             "screw angle that actually belongs to x.\nRight: two screw turns "
             f"({ZOOM[0]:g}{MINUS}{ZOOM[1]:g} mm) at the same scale, where the "
             "shape of each term is visible. The strips add up to the bottom "
             "one.",
             fontsize=8, color=INK2, va="top", linespacing=1.4)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e2_model_exploded_components.png", dpi=300)
    plt.close(fig)
    return full


def figure_table(model):
    lo, hi = model["state"]["domain"]
    x = np.linspace(lo, hi, 900)
    theta = np.linspace(0, 2 * np.pi, 361)
    grid = np.zeros((theta.size, x.size))
    for k in range(1, model["state"]["H"] + 1):
        alpha, beta = model["state"]["groups"][f"C{k}"]
        if model["state"]["Kc"] is not None and k <= 2:
            basis = m.raw_bspline(x, model["state"]["Kc"],
                                  model["state"]["domain"])
        else:
            basis = np.ones((x.size, 1))
        a = basis @ model["beta"][alpha]
        b = basis @ model["beta"][beta]
        grid += (np.cos(k * theta)[:, None] * a[None, :]
                 + np.sin(k * theta)[:, None] * b[None, :])
    limit = float(np.max(np.abs(grid)))

    fig = plt.figure(figsize=(9.9, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1.0],
                          width_ratios=[1.0, .022], left=.075, right=.925,
                          top=.845, bottom=.095, hspace=.40, wspace=.03)
    ax = fig.add_subplot(gs[0, 0])
    image = ax.pcolormesh(x, np.degrees(theta), grid, cmap=DIVERGING,
                          norm=TwoSlopeNorm(0, -limit, limit), shading="gouraud")
    for position in SLICES:
        ax.axvline(position, color="#20262b", lw=.9, ls=(0, (4, 2)))
        # Inside the map, so the panel title above keeps the whole width.
        ax.text(position, 14, f"{position:g} mm", ha="center", va="bottom",
                fontsize=7.6, color="white").set_path_effects(
                    [pe.withStroke(linewidth=1.8, foreground="#20262b")])
    ax.set(xlabel="commanded position x [mm]",
           ylabel="screw angle $\\theta$ [deg]", ylim=(0, 360))
    ax.set_yticks(np.arange(0, 361, 90))
    ax.set_title("(a)  The periodic part of the model over position and screw "
                 "angle - the compensation table", loc="left", pad=8)

    bar = fig.add_subplot(gs[0, 1])
    fig.colorbar(image, cax=bar).set_label("deviation [µm]", fontsize=8.4)

    bx = fig.add_subplot(gs[1, 0])
    for position, dash in zip(SLICES, ("solid", (0, (5, 2)), (0, (1.4, 1.6)))):
        column = np.argmin(np.abs(x - position))
        bx.plot(np.degrees(theta), grid[:, column], color=INK, lw=1.2, ls=dash,
                label=f"x = {position:g} mm   "
                      f"(peak-to-peak {num(np.ptp(grid[:, column]))} µm)")
    bx.axhline(0, color="#aeb4b9", lw=.6)
    bx.set(xlabel="screw angle $\\theta$ [deg]", ylabel="deviation [µm]",
           xlim=(0, 360))
    bx.set_xticks(np.arange(0, 361, 45))
    bx.set_title("(b)  One screw turn, at three positions along the travel",
                 loc="left", pad=6)
    bx.legend(frameon=False, loc="lower right", ncol=3, fontsize=7.4,
              columnspacing=1.4)
    tidy(bx)
    bx.set_ylim(-1.35 * limit, 1.35 * limit)

    fig.text(.075, .975,
             "The table the firmware would apply: one screw revolution is "
             "2 mm of travel", fontsize=12.5, color="#20262b", va="top")
    fig.text(.075, .938,
             "Harmonics 1 and 2 carry position-dependent coefficients, 3 and 4 "
             "are constant along the travel; the slow trend and the gain are "
             "not shown here - they are a separate one-dimensional table.",
             fontsize=8, color=INK2, va="top")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e2_model_xtheta_table.png", dpi=300)
    plt.close(fig)
    return grid


def main():
    m.style()
    plt.rcParams.update({"axes.titlesize": 9.4, "font.size": 8})
    data, model = model_fit()
    print(f"model edf {model['edf']:.2f}, p {model['p']}, rms {model['rms']:.6f}")
    full = figure_exploded(model)
    grid = figure_table(model)
    theta = m.screw_theta(data["x"])
    check = np.abs(m.predict(model, data["x"], theta)
                   - (model["beta"][0] + sum(
                       parts(model, data["x"])[key]
                       for key in ("gain", "slow", "k1", "k2", "k34")))).max()
    print(f"decomposition closes to {check:.3e} um")
    for key in ("gain", "slow", "k1", "k2", "k34", "sum"):
        print(f"  {key:5s} peak-to-peak {np.ptp(full[key]):6.3f} um, "
              f"rms {np.sqrt(np.mean(full[key]**2)):6.3f} um")
    print(f"periodic table: max |value| {np.max(np.abs(grid)):.3f} um")


if __name__ == "__main__":
    main()

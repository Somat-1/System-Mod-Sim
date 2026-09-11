"""Decompose the TMC2209ExpRun 25 mm triangular ramps.

Answers four questions about the campaign capture:

1. Why the resting position drifts downward - by splitting each leg into its
   outbound and return travel, so the asymmetry that produces the drift is
   visible rather than inferred.
2. Settling error at the top and bottom endpoints - how far the axis creeps
   during the 1 s dwell when it is nominally stopped.
3. DIRECT vs INDIVIDUAL step generation at matched MRES and rate.
4. (documented, not computed) The chopper configuration of this run.

Campaign structure, from esp32_v4_mres_trajectory_campaign.ino: for each MRES
in (1, 4, 16, 32), six 25 mm out-and-back trajectories run in the order
DIRECT slow/moderate/fast then INDIVIDUAL slow/moderate/fast. So within each
group of six detected legs, legs 1-3 are DIRECT and legs 4-6 are INDIVIDUAL,
and leg k pairs with leg k+3 at the same commanded rate.

CONFIGURATION OF THIS RUN: MicroPlyer interpolation disabled, chopper left at
its power-on default -> StealthChop. Verified from git: the binary flashed at
16:59 and re-triggered at 17:20 contains intpol(false) and no en_spreadCycle
call; the SpreadCycle fix was committed at 17:51, after this run ended.
"""
import numpy as np
import matplotlib.pyplot as plt

from ids_common import (ANALYSIS_DIR, RUNS_DIR, C_MEASURED, C_COMMANDED,
                        C_INK, C_MUTED, apply_style, load_ids)
from plot_tmc2209_exp_run import find_trajectories, TRAJECTORY_MM

OUT_DIR = ANALYSIS_DIR / "tmc2209_exp_run"
CSV = RUNS_DIR / "TMC2209ExpRun.csv"

MRES_VALUES = [1, 4, 16, 32]
MODES = ["DIRECT", "DIRECT", "DIRECT", "INDIVIDUAL", "INDIVIDUAL", "INDIVIDUAL"]
RATE_NAMES = ["slow", "moderate", "fast"] * 2
RATES_FSPS = [27.5, 70.0, 200.0] * 2
STILL_UM_S = 20.0          # local slope below this counts as "stopped"
DWELL_S = 1.0              # commanded endpoint dwell


def dwell_after(still, i, dt, need_s=0.40, measure_s=1.0,
                limit_idx=None):
    """Measure the settling window that begins when motion actually stops.

    Walks forward from index i to the first point that is still and STAYS
    still for `need_s`, which stops the search latching onto a momentary
    quiet patch mid-ramp, then returns the following `measure_s`. Anchoring on
    the moment motion ceases is what makes the number a settling error rather
    than just the flattest stretch nearby.
    """
    n = still.size
    need = int(need_s / dt)
    span = int(measure_s / dt)
    j = max(0, int(i))
    # A fixed forward window cannot work here: the slow legs descend from
    # the 8 mm detection threshold to rest at 275 um/s, which takes 29 s,
    # while a fast leg takes 4 s. The caller supplies the real limit.
    hard = int(i) + int(45.0 / dt) if limit_idx is None else int(limit_idx)
    limit = min(n - need, hard)
    while j < limit:
        if still[j] and still[j:j + need].all():
            b = j
            while b < n - 1 and still[b + 1]:
                b += 1
            return j, min(b, j + span)
        j += 1
    return None


def main():
    apply_style()
    t, pos, meta = load_ids(CSV)
    pos_mm = pos / 1000.0
    dt = meta["dt_s"]
    # Sample-to-sample slope is useless here: 3 nm of encoder noise over a
    # 1 ms sample is 4.5 um/s of pure jitter, which fragments every dwell.
    # Difference over a 50 ms window instead - noise contribution drops to
    # ~0.1 um/s while the slowest commanded ramp is 275 um/s.
    wv = int(round(0.050 / dt))
    fwd = np.empty_like(pos)
    fwd[:-wv] = np.abs(pos[wv:] - pos[:-wv]) / (wv * dt)
    fwd[-wv:] = fwd[-wv - 1]
    bwd = np.empty_like(pos)
    bwd[wv:] = fwd[:-wv]
    bwd[:wv] = fwd[0]
    slope = np.minimum(fwd, bwd)      # still if either side is quiet
    still = slope < STILL_UM_S
    max_span = int(3.0 / dt)

    spans = find_trajectories(t, pos, dt)
    print(f"legs detected: {len(spans)}")

    rows = []
    for i, (s0, apex, e0) in enumerate(spans):
        g0 = spans[i - 1][2] if i else 0
        pre, post = pos_mm[g0:s0], pos_mm[e0:spans[i + 1][0] if i + 1 < len(spans) else pos.size]
        pre_s, post_s = still[g0:s0], still[e0:e0 + len(post)]
        rest_before = float(np.median(pre[pre_s])) if pre_s.sum() > 50 else float(np.median(pre))
        rest_after = float(np.median(post[post_s])) if post_s.sum() > 50 else float(np.median(post))

        # Top dwell: settling from the moment the outbound ramp stops.
        # Search from BEFORE the apex - creep during the dwell puts the
        # argmax at the dwell's end, and searching forward from there
        # finds only the descent.
        top = dwell_after(still, apex - int(3.0 / dt), dt,
                          limit_idx=apex + int(3.0 / dt))
        top_drift = (pos[top[1]] - pos[top[0]]) if top else np.nan
        top_dur = (top[1] - top[0]) * dt if top else np.nan

        # Bottom dwell: from where the return ramp stops, searched up to
        # the start of the next leg.
        nxt = spans[i + 1][0] if i + 1 < len(spans) else pos.size - 1
        bot = dwell_after(still, e0, dt, limit_idx=nxt)
        bot_drift = (pos[bot[1]] - pos[bot[0]]) if bot else np.nan
        bot_dur = (bot[1] - bot[0]) * dt if bot else np.nan

        out_mm = pos_mm[apex] - rest_before
        back_mm = pos_mm[apex] - rest_after
        rows.append(dict(
            leg=i + 1, mres=MRES_VALUES[i // 6], mode=MODES[i % 6],
            rate=RATE_NAMES[i % 6], fsps=RATES_FSPS[i % 6],
            rest_before=rest_before, apex=pos_mm[apex], rest_after=rest_after,
            out=out_mm, back=back_mm, net=rest_after - rest_before,
            top_drift_um=top_drift, top_dur=top_dur,
            bot_drift_um=bot_drift, bot_dur=bot_dur))

    # ---- 1. drift decomposition -----------------------------------------
    print("\n=== 1. OUTBOUND vs RETURN TRAVEL (why it drifts) ===")
    print(" leg MRES mode        rate      out[mm]   back[mm]   out-back[um]   rest_after[mm]")
    for r in rows:
        print(f"{r['leg']:4d} {r['mres']:4d} {r['mode']:11s} {r['rate']:9s} "
              f"{r['out']:9.4f} {r['back']:10.4f} {1000*(r['out']-r['back']):13.1f} "
              f"{r['rest_after']:14.4f}")
    net = np.array([r["net"] for r in rows])
    print(f"\n  total drift {1000*net.sum():+.0f} um; "
          f"legs with |net| > 50 um: {[r['leg'] for r in rows if abs(r['net'])>0.05]}")

    # ---- 2. endpoint settling -------------------------------------------
    print("\n=== 2. SETTLING DURING THE 1 s ENDPOINT DWELLS ===")
    print(" leg MRES mode        rate      top_drift[um] (dur)   bottom_drift[um] (dur)")
    for r in rows:
        print(f"{r['leg']:4d} {r['mres']:4d} {r['mode']:11s} {r['rate']:9s} "
              f"{r['top_drift_um']:13.3f} ({r['top_dur']:4.2f}s) "
              f"{r['bot_drift_um']:16.3f} ({r['bot_dur']:4.2f}s)")
    for lbl, key in (("top", "top_drift_um"), ("bottom", "bot_drift_um")):
        v = np.array([r[key] for r in rows], dtype=float)
        v = v[np.isfinite(v)]
        print(f"  {lbl:6s}: mean {v.mean():+7.3f} um, |mean| {np.abs(v).mean():6.3f}, "
              f"max |drift| {np.abs(v).max():6.3f} um")

    # ---- 3. DIRECT vs INDIVIDUAL ----------------------------------------
    print("\n=== 3. DIRECT vs INDIVIDUAL (same MRES, same rate) ===")
    print(" MRES rate       out_D[mm] out_I[mm]  d_out[um] | net_D[um] net_I[um] | topD topI [um]")
    pairs = []
    for gi, mres in enumerate(MRES_VALUES):
        g = rows[gi * 6:(gi + 1) * 6]
        for k in range(3):
            D, I = g[k], g[k + 3]
            pairs.append((mres, D["rate"], D, I))
            print(f"{mres:5d} {D['rate']:10s} {D['out']:9.4f} {I['out']:9.4f} "
                  f"{1000*(I['out']-D['out']):9.1f} | {1000*D['net']:8.1f} {1000*I['net']:9.1f} "
                  f"| {D['top_drift_um']:5.2f} {I['top_drift_um']:5.2f}")
    for lbl, sel in (("DIRECT", lambda r: r["mode"] == "DIRECT"),
                     ("INDIVIDUAL", lambda r: r["mode"] == "INDIVIDUAL")):
        sub = [r for r in rows if sel(r)]
        err = np.array([abs(r["out"] - TRAJECTORY_MM) for r in sub]) * 1000
        top = np.abs(np.array([r["top_drift_um"] for r in sub], dtype=float))
        nets = np.abs(np.array([r["net"] for r in sub])) * 1000
        print(f"  {lbl:11s} travel error {err.mean():7.1f} um mean | "
              f"top settling {np.nanmean(top):6.3f} um | |net drift| {nets.mean():7.1f} um")

    # ---- plot ------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11.5, 7.2), sharex=True)
    legs = np.arange(1, len(rows) + 1)
    isD = np.array([r["mode"] == "DIRECT" for r in rows])
    outs = np.array([r["out"] for r in rows])
    backs = np.array([r["back"] for r in rows])
    a1.plot(legs, 1000 * (outs - backs), color=C_MEASURED, lw=1.4, zorder=2)
    a1.scatter(legs[isD], 1000 * (outs - backs)[isD], s=52, color=C_MEASURED,
               zorder=3, label="DIRECT")
    a1.scatter(legs[~isD], 1000 * (outs - backs)[~isD], s=52, facecolors="white",
               edgecolors=C_MEASURED, linewidths=1.6, zorder=3, label="INDIVIDUAL")
    a1.axhline(0, color=C_COMMANDED, ls="--", lw=1.3, label="no asymmetry")
    a1.set_ylabel("outbound − return [µm]")
    a1.set_title("Per-leg travel asymmetry — positive means the return "
                 "overshoots, pulling the origin down", loc="left")
    a1.legend(loc="upper right", ncol=3)

    tops = np.abs(np.array([r["top_drift_um"] for r in rows], dtype=float))
    bots = np.abs(np.array([r["bot_drift_um"] for r in rows], dtype=float))
    a2.plot(legs, tops, color=C_MEASURED, marker="o", ms=5, lw=1.4,
            label="top endpoint")
    a2.plot(legs, bots, color=C_COMMANDED, marker="s", ms=5, lw=1.4, ls="--",
            label="bottom endpoint")
    a2.set_yscale("log")
    a2.set_xlabel("Trajectory leg number (chronological)")
    a2.set_ylabel("|creep during 1 s dwell| [µm]")
    a2.set_title("Settling error while nominally stopped", loc="left")
    a2.set_xticks(legs[::2])
    a2.legend(loc="upper right", ncol=2)
    for a in (a1, a2):
        for g in range(6, len(rows), 6):
            a.axvline(g + 0.5, color=C_MUTED, lw=0.8, ls=":")
        for g in range(3, len(rows), 6):
            a.axvline(g + 0.5, color=C_MUTED, lw=0.6, ls="-", alpha=0.35)
    fig.suptitle("TMC2209ExpRun — ramp asymmetry and endpoint settling   "
                 "(dotted = MRES group, faint solid = DIRECT|INDIVIDUAL split)",
                 x=0.005, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT_DIR / "ramp_asymmetry_and_settling.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("\nwrote ramp_asymmetry_and_settling.png")


if __name__ == "__main__":
    main()

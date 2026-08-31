#!/usr/bin/env python3
"""Build human-screened A2 montages with Rev 4 model overlays."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("preemptive", HERE / "analyze_preemptive_run.py")
p = importlib.util.module_from_spec(spec); spec.loader.exec_module(p)
MODELS = np.load(HERE / "curated_model_simulations.npz")
SIZES = (1, 2, 4, 8, 16, 32)


def a2_events(rows, run, n):
    return [r for r in rows if int(r["run_index"] or 0) == run and r["block"] == "A2"
            and r["label"] == f"N{n}_alternating" and r["event"] == "PULSE_BATCH"][:12]


def robust_levels(time, values, edges):
    levels=[]
    for i, edge in enumerate(edges):
        end = edges[i+1] if i+1 < len(edges) else edge + 0.4
        mask=(time >= edge+0.06) & (time <= end-0.04)
        if np.any(mask): levels.append(float(np.median(values[mask])))
    low=np.median(levels[1::2]); high=np.median(levels[0::2])
    return float(low), float(high)


def main():
    ids_t, ids_y=p.load_ids(); rows=p.load_esp(); esp0,_,cmd_t,cmd_y=p.command_trace(rows)
    offset,*_=p.find_alignment(ids_t,ids_y,cmd_t,cmd_y)
    payload={}; calibration=[]
    for n in SIZES:
        for run in (1,2):
            ev=a2_events(rows,run,n); first=int(ev[0]["timestamp_us"])/1e6-esp0
            edge=np.asarray([int(r["timestamp_us"])/1e6-esp0-first for r in ev])
            end=edge[-1]+0.45; local=np.arange(-0.1,end,0.001)
            measured=np.interp(local+first+offset,ids_t,ids_y)
            lo,hi=robust_levels(local,measured,edge)
            amplitude=hi-lo
            expected=n*5.0/16.0
            calibration.append((n,run,amplitude/expected))
            payload[(n,run)] = (local, (measured-lo)/amplitude, edge)

    # Primary montage: synchronized A2 family, normalized per panel because the
    # EL5101 export has no documented counts-to-length scale in this dataset.
    fig,axes=plt.subplots(3,2,figsize=(15,12),sharex=False,constrained_layout=True)
    for ax,n in zip(axes.flat,SIZES):
        for run,color,name in ((1,"#1769aa","measured I_lo"),(2,"#31a354","measured I_mid")):
            t,y,_=payload[(n,run)]
            ax.plot(t,y,color=color,lw=.75,alpha=.85,label=name)
        amp=n*5.0/16.0
        for kind,color,style,name in (
            ("tangent","#6a51a3","-","frictionless tangent"),
            ("detent","#d95f0e","--","exact-detent frictionless"),
            ("lugre","#cb181d","-.","Rev 4.2 parallel LuGre")):
            mt=MODELS[f"N{n}_{kind}_t"]; my=MODELS[f"N{n}_{kind}_um"]/amp
            ax.plot(mt,my,color=color,ls=style,lw=1.15,label=name)
        ct=MODELS[f"N{n}_command_t"]; cy=MODELS[f"N{n}_command_um"]/amp
        ax.step(ct,cy,where="post",color="#555",lw=.8,alpha=.75,label="command")
        ax.set_title(f"A2 alternating, N={n} pulses ({amp:.4g} µm ideal amplitude)")
        ax.set_xlabel("Time from first pulse [s]"); ax.set_ylabel("Normalized displacement")
        ax.grid(alpha=.22); ax.set_xlim(-.1,ct[-1])
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",ncol=3,fontsize=9)
    fig.suptitle("Best synchronized family: A2 alternating responses vs Rev 4 models\n"
                 "Shape/timing comparison only — per-panel normalization hides invalid absolute gain",fontsize=14)
    fig.savefig(HERE/"plots"/"02_curated_A2_model_overlay_montage.png",dpi=180)
    plt.close(fig)

    # Cross-current repeatability montage in raw counts, deliberately without
    # per-trace normalization so acquisition quality remains visible.
    fig,axes=plt.subplots(3,2,figsize=(15,11),constrained_layout=True)
    for ax,n in zip(axes.flat,SIZES):
        for run,color,name in ((1,"#1769aa","I_lo"),(2,"#31a354","I_mid")):
            t,y,_=payload[(n,run)]
            raw_amp=next(v for nn,rr,v in calibration if nn==n and rr==run)
            ax.plot(t,y*raw_amp*(n*5/16),color=color,lw=.75,label=name)
        ax.set_title(f"A2 N={n}: baseline-referenced encoder counts")
        ax.set(xlabel="Time from first pulse [s]",ylabel="Encoder-count change")
        ax.grid(alpha=.22); ax.legend()
    fig.suptitle("Accepted hardware evidence: A2 repeats coherently at both completed currents")
    fig.savefig(HERE/"plots"/"03_curated_A2_measured_repeatability_montage.png",dpi=180)
    plt.close(fig)

    # Gain-consistency check. Normalize hardware gain to N=16 separately for
    # each current; a displacement sensor/model fit should stay near one.
    cal={(n,r):v for n,r,v in calibration}
    fig,ax=plt.subplots(figsize=(9,5.5),constrained_layout=True)
    for run,color,name in ((1,"#1769aa","measured I_lo"),(2,"#31a354","measured I_mid")):
        relative=np.asarray([cal[(n,run)]/cal[(16,run)] for n in SIZES])
        ax.plot(SIZES,relative,"o-",color=color,label=name)
    model_gain=[]
    for n in SIZES:
        ideal=n*5/16; y=MODELS[f"N{n}_lugre_um"]
        model_gain.append((np.percentile(y,90)-np.percentile(y,10))/ideal)
    model_gain=np.asarray(model_gain); model_gain/=model_gain[SIZES.index(16)]
    ax.plot(SIZES,model_gain,"s--",color="#cb181d",label="Rev 4.2 LuGre")
    ax.axhline(1,color="#555",lw=.8); ax.set_xscale("log",base=2); ax.set_yscale("log")
    ax.set_xticks(SIZES,labels=[str(n) for n in SIZES]); ax.grid(alpha=.25,which="both")
    ax.set(xlabel="A2 step size N [driver pulses]",ylabel="Relative gain (normalized at N=16)",
           title="Gain-consistency rejection: measured plateau does not scale with commanded travel")
    ax.legend(); fig.savefig(HERE/"plots"/"04_A2_gain_consistency_rejection.png",dpi=180); plt.close(fig)

    # A single-scale overlay anchored only at N=16. Unlike the normalized
    # shape montage, this deliberately exposes the amplitude mismatch at all
    # other N values. It is a diagnostic, not a fitted calibration claim.
    anchor_scale=float(np.median([cal[(16,1)],cal[(16,2)]]))
    fig,axes=plt.subplots(3,2,figsize=(15,12),constrained_layout=True)
    for ax,n in zip(axes.flat,SIZES):
        for run,color,name in ((1,"#1769aa","measured I_lo"),(2,"#31a354","measured I_mid")):
            t,y,_=payload[(n,run)]; amp=cal[(n,run)]*(n*5/16)
            ax.plot(t,y*amp,color=color,lw=.75,alpha=.85,label=name)
        for kind,color,style,name in (
            ("tangent","#6a51a3","-","frictionless tangent"),
            ("detent","#d95f0e","--","exact-detent frictionless"),
            ("lugre","#cb181d","-.","Rev 4.2 parallel LuGre")):
            ax.plot(MODELS[f"N{n}_{kind}_t"],MODELS[f"N{n}_{kind}_um"]*anchor_scale,
                    color=color,ls=style,lw=1.15,label=name)
        ax.set_title(f"A2 N={n}; one scale anchored at N=16")
        ax.set(xlabel="Time from first pulse [s]",ylabel="Encoder counts / model-equivalent counts")
        ax.grid(alpha=.22); ax.set_xlim(-.1,MODELS[f"N{n}_command_t"][-1])
    handles,labels=axes.flat[0].get_legend_handles_labels(); fig.legend(handles,labels,loc="lower center",ncol=3,fontsize=9)
    fig.suptitle(f"Single-scale model overlay (anchor={anchor_scale:.3g} counts/µm at N=16)\n"
                 "Amplitude disagreement away from the anchor demonstrates that this is not a valid model fit")
    fig.savefig(HERE/"plots"/"05_A2_single_scale_model_overlay_rejection_montage.png",dpi=180); plt.close(fig)

    with (HERE/"a2_plateau_calibration_diagnostic.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["N_pulses","run_index","inferred_counts_per_um"]); w.writerows(calibration)

    with (HERE/"curated_selection.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["block","condition","decision","reason"])
        for n in SIZES: w.writerow(["A2",f"N{n}_alternating","timing/shape only","coherent reversals, but plateau amplitude fails N-scaling check"])
        w.writerow(["B","descending","secondary","repeatable but switching artifacts and weaker command-shape agreement"])
        w.writerow(["B","minor","secondary","repeatable but not clean enough for primary model-fit claims"])
        w.writerow(["A1","all","rejected","creep/drift dominates; weak monotonic command tracking"])
        w.writerow(["E","all","rejected","isolated spikes dominate; doublet duration approaches acquisition resolution"])
        w.writerow(["BLOCK_0","start/end","diagnostic only","repeatable fingerprint but not a clean gain comparison"])

    ratios=np.asarray([x[2] for x in calibration])
    report=f"""# Curated preemptive-run assessment

## Best synchronized evidence (not an accepted quantitative fit)

Block A2 at N = 1, 2, 4, 8, 16, and 32 pulses is the clearest synchronized
family. All conditions show the alternating timing and reproduce across the two
completed current settings. However, the raw plateau excursion remains roughly
similar while commanded travel changes by 32x. Therefore A2 is retained only as
timing/shape evidence, not as a valid quantitative model fit.

## Model overlay convention

The ESP commands drive the Rev 4 frictionless tangent, exact periodic-detent
frictionless, and Rev 4.2 parallel-LuGre models. Model output is stage coordinate
`x_n`. The EL5101 CSV contains counter values but no documented counts-to-length
calibration, so each measured panel is offset and amplitude-normalized to its own
two plateau medians. The overlays test timing and normalized response shape, not
absolute gain. A physical gain fit would be misleading until counter calibration
is supplied. The command-inferred ratios span {ratios.min():.3g} to {ratios.max():.3g}
counts/µm, confirming that one global inferred scale is not defensible here.
`04_A2_gain_consistency_rejection.png` exposes this rather than normalizing it away.
`05_A2_single_scale_model_overlay_rejection_montage.png` uses one scale anchored
at N=16 and shows the resulting amplitude mismatch at every other step size.

## Bottom line

No segment in this preemptive capture supports an absolute measured-versus-model
fit. The A2 montage is the best timing-aligned comparison available; its normalized
overlay must not be interpreted as parameter validation or gain agreement.

## Excluded from fit claims

- A1: dominated by drift/creep and inconsistent monotonic response.
- E: spike-dominated; doublets are too fast relative to this acquisition for a
  trustworthy displacement comparison.
- B descending/minor: retained only as secondary qualitative evidence.
- Block 0: useful as a repeatability fingerprint, not a clean model-gain test.
"""
    (HERE/"CURATED_REPORT.md").write_text(report,encoding="utf-8")
    print(report)


if __name__=="__main__": main()

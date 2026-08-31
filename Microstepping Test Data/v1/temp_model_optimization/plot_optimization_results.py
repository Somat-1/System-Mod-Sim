#!/usr/bin/env python3
"""Plot the multi-start Guyan/LuGre optimization output."""
import json
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent

def main():
    d=np.load(HERE/"optimization_results.npz",allow_pickle=False); summary=json.loads((HERE/"optimization_summary.json").read_text()); out=HERE/"plots"; out.mkdir(exist_ok=True)
    fig,axes=plt.subplots(1,3,figsize=(16,5),constrained_layout=True)
    for ax,n in zip(axes,(1,2,16)):
        t=d[f"step{n}_time_s"]; ax.plot(t,d[f"step{n}_measured_um"],color="#222",lw=1,label="IDS measured"); ax.plot(t,d[f"step{n}_simulated_um"],color="#c0392b",lw=1.2,label="optimized Guyan+LuGre");
        for edge in d[f"step{n}_edges_s"]: ax.axvline(edge,color="#999",lw=.35,alpha=.4)
        ax.set(title=f"StepSize{n}",xlabel="Time [s]",ylabel="Stage displacement [µm]"); ax.grid(alpha=.2)
    axes[0].legend(); fig.suptitle("Best multi-start reduced Guyan + Rev 4.2 LuGre calibration"); fig.savefig(out/"best_fit_overlay.png",dpi=180); plt.close(fig)
    solutions=d["solutions"]; fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True); ax.scatter(solutions[:,0],solutions[:,1],label="initial",s=55); ax.scatter(solutions[:,0],solutions[:,2],label="accepted terminal",s=55); ax.set_yscale("log"); ax.set(xlabel="Start index",ylabel="Weighted squared physical residual",title="Optimization outcomes"); ax.grid(alpha=.25); ax.legend(); fig.savefig(out/"multistart_costs.png",dpi=180); plt.close(fig)
    history_path=HERE/"refined_softl1_run_log.jsonl"
    if history_path.exists():
        records=[json.loads(line) for line in history_path.read_text().splitlines() if line.strip()]
        starts={r["start"]:r["initial_cost"] for r in records if r.get("event")=="START"}
        terminals={}
        for r in records:
            if r.get("event")=="GROUP_COMPLETE": terminals[r["start"]]=r["cost"]
        ids=sorted(starts)
        fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True)
        ax.scatter(ids,[starts[i] for i in ids],s=70,label="refined start (consistent SSE)",color="#2471a3")
        ax.scatter(ids,[terminals[i] for i in ids],s=80,marker="x",linewidths=2,label="soft-L1 endpoint (rejected)",color="#c0392b")
        selected=min(ids,key=lambda i: starts[i])
        ax.scatter([selected],[summary["optimized_cost"]],s=170,marker="*",label="selected + monotonic polish",color="#1e8449",zorder=5)
        ax.set_yscale("log"); ax.set_xticks(ids); ax.set(xlabel="Refined start index",ylabel="Weighted squared physical residual",title="Refined multi-start selection and corrected polish"); ax.grid(alpha=.25); ax.legend()
        fig.savefig(out/"refined_multistart_selection.png",dpi=180); plt.close(fig)
    names=[str(x) for x in d["variable_names"]]; x=d["best_x"]; fig,ax=plt.subplots(figsize=(12,6),constrained_layout=True); ax.bar(np.arange(len(x)),np.exp(x),color="#5e3c99"); ax.axhline(1,color="#555",lw=.8); ax.set_yscale("log"); ax.set_xticks(np.arange(len(x)),names,rotation=65,ha="right"); ax.set(ylabel="Multiplier relative to Rev 4.2",title="Optimized detent and friction parameter multipliers"); ax.grid(axis="y",alpha=.25); fig.savefig(out/"optimized_parameter_multipliers.png",dpi=180); plt.close(fig)
    print("Wrote",out)
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Generate optimized-parameter Bode and stepping montage."""
from __future__ import annotations
import importlib.util, json
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent; PROJECT=HERE.parents[2]
MODEL_PATH=PROJECT/"Rev 4"/"Guyan Model Reduction"/"Guyan w Friction"/"guyan_friction_model.py"
spec=importlib.util.spec_from_file_location("guyan_friction",MODEL_PATH); gm=importlib.util.module_from_spec(spec); spec.loader.exec_module(gm)

def response(model,f):
    A,B,C,D=model.analytical_linearization(); eye=np.eye(A.shape[0]); return np.asarray([(C@np.linalg.solve(2j*np.pi*ff*eye-A,B)+D)[0,0] for ff in f])

def main():
    summary=json.loads((HERE/"optimization_summary.json").read_text()); baseline=gm.load_parameters(); optimized=dict(baseline); optimized.update(summary["optimized_parameters"])
    f=np.logspace(-2,np.log10(8000),2200); h0=response(gm.GuyanFrictionModel(baseline),f); h1=response(gm.GuyanFrictionModel(optimized),f); out=HERE/"plots"; out.mkdir(exist_ok=True)
    fig,axes=plt.subplots(2,1,figsize=(10.5,8),sharex=True,constrained_layout=True)
    for h,color,label in ((h0,"#777","Rev 4.2 starting parameters"),(h1,"#b2182b","selected refined parameters")):
        axes[0].semilogx(f,20*np.log10(np.maximum(np.abs(h),1e-300)),color=color,lw=1.2,label=label); axes[1].semilogx(f,np.unwrap(np.angle(h))*180/np.pi,color=color,lw=1.2)
    axes[0].set_ylabel(r"Magnitude $|x_n/\theta_{cmd}|$ [dB re m/rad]"); axes[1].set(ylabel="Phase [deg]",xlabel="Frequency [Hz]")
    for ax in axes: ax.grid(True,which="both",alpha=.25)
    axes[0].legend(); fig.suptitle("Reduced Guyan + LuGre Bode: selected v1 calibration candidate\nSmall-signal rest-point tangent")
    fig.savefig(out/"optimized_bode.png",dpi=180); plt.close(fig)
    np.savez_compressed(HERE/"optimized_bode_data.npz",frequency_hz=f,baseline_response=h0,optimized_response=h1)

    d=np.load(HERE/"optimization_results.npz",allow_pickle=False); fig,axes=plt.subplots(3,1,figsize=(12,10),constrained_layout=True)
    for ax,n in zip(axes,(1,2,16)):
        t=d[f"step{n}_time_s"]; ax.plot(t,d[f"step{n}_measured_um"],color="#222",lw=1,label="IDS physical"); ax.plot(t,d[f"step{n}_simulated_um"],color="#b2182b",lw=1.2,label="selected Guyan+LuGre")
        for edge in d[f"step{n}_edges_s"]: ax.axvline(edge,color="#999",lw=.35,alpha=.35)
        rmse=summary["physical_rmse_um_by_step_size"][str(n)]; ax.set(title=f"StepSize{n} — RMSE {rmse:.4g} µm",ylabel="Stage displacement [µm]"); ax.grid(alpha=.22)
    axes[-1].set_xlabel("Time [s]"); axes[0].legend(); fig.suptitle("v1 stepping montage with selected refined parameters")
    fig.savefig(out/"optimized_stepping_montage.png",dpi=180); plt.close(fig)
    print("Wrote optimized_bode.png and optimized_stepping_montage.png")
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Recover the lowest consistent weighted-SSE candidate from the refined log."""
import importlib.util, json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("opt",HERE/"optimize_guyan_lugre.py"); opt=importlib.util.module_from_spec(spec); spec.loader.exec_module(opt)

def main():
    events=[json.loads(line) for line in (HERE/"optimization_run_log.jsonl").read_text().splitlines() if line.strip()]
    candidates=[e for e in events if e["event"] in ("START","GROUP_COMPLETE")]
    best=min(candidates,key=lambda e:e.get("initial_cost",e.get("cost",float("inf"))))
    best_cost=float(best.get("initial_cost",best.get("cost"))); x=np.asarray(best["x"],float)
    config=json.loads(opt.CONFIG_FILE.read_text()); base=json.loads(opt.PARAM_FILE.read_text())["parameters"]; data=opt.prepare_data(config); p=opt.decode(base,x); predictions=[opt.simulate(p,d) for d in data]
    starts=[e for e in events if e["event"]=="START"]; terminal=[]
    for s in starts:
        groups=[e for e in events if e["event"]=="GROUP_COMPLETE" and e["start"]==s["start"]]
        terminal_cost=groups[-1]["cost"] if groups else s["initial_cost"]
        terminal.append([s["start"],s["initial_cost"],terminal_cost,0.0]+s["x"])
    arrays={"variable_names":np.asarray(opt.VARIABLES),"solutions":np.asarray(terminal),"history":np.asarray([[e["start"],e["cycle"],e["group"],e["cost"],0] for e in events if e["event"]=="GROUP_COMPLETE"]),"best_x":x}
    for d,pred in zip(data,predictions): arrays["step%d_time_s"%d["n"]]=d["t"]; arrays["step%d_measured_um"%d["n"]]=d["y_um"]; arrays["step%d_simulated_um"%d["n"]]=pred; arrays["step%d_edges_s"%d["n"]]=d["edge_t"]
    np.savez_compressed(opt.OUTPUT,**arrays)
    rmse={str(d["n"]):float(np.sqrt(np.mean((pred-d["y_um"])**2))) for d,pred in zip(data,predictions)}
    summary={"selection":"lowest consistent weighted SSE among all logged starts/group endpoints","source_event":best["event"],"best_start":best["start"],"optimized_cost":best_cost,"physical_rmse_um_by_step_size":rmse,"variables":dict(zip(opt.VARIABLES,x.tolist())),"optimized_parameters":{k:p[k] for k in ["T_d"]+[v for v in p if any(v.startswith(s) for s in ("sigma","Fc_","Fs_","Tc_","Ts_","vs_"))]},"locked_structural_parameters":{k:base[k] for k in opt.STRUCTURAL},"config":config,"warning":"terminal soft-L1 points were rejected when their weighted squared physical residual increased"}
    opt.SUMMARY.write_text(json.dumps(summary,indent=2)+"\n"); print(json.dumps(summary,indent=2))
if __name__=="__main__":main()

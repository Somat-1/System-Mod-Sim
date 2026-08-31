#!/usr/bin/env python3
"""Multi-start, alternating-group calibration of reduced Guyan + LuGre.

Run with macOS system Python 3.9, which owns this machine's SciPy install:
    /usr/bin/python3 optimize_guyan_lugre.py
The script writes numerical/checkpoint files only. Plotting is separate.
"""

import csv
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.signal import find_peaks, medfilt


HERE=Path(__file__).resolve().parent; V1=HERE.parent; PROJECT=V1.parents[1]
PARAM_FILE=PROJECT/"Rev 4"/"lugre_friction"/"Rev 4.2"/"model_parameters.json"
CONFIG_FILE=HERE/"optimization_config.json"; OUTPUT=HERE/"optimization_results.npz"
SUMMARY=HERE/"optimization_summary.json"; CHECKPOINT=HERE/"optimization_checkpoint.json"
RUN_LOG=HERE/"optimization_run_log.jsonl"
FULL_STEP_M=5e-6; EDGE_THRESHOLDS={1:800,2:400,16:100}; PORTS=("way","nut","sb")
STRUCTURAL=("k_c","k_s1","k_s2","k_nut","k_brg","c_c","c_s1","c_s2","c_nut","c_brg","c_EM")
VARIABLES=("T_d",
 "sigma0_way","sigma1_way","sigma2_way","Fc_way","Fs_gap_way","vs_way",
 "sigma0_nut","sigma1_nut","sigma2_nut","Fc_nut","Fs_gap_nut","vs_nut",
 "sigma0_sb","sigma1_sb","sigma2_sb","Tc_sb","Ts_gap_sb","vs_sb")
GROUPS=((0,),tuple(range(1,7)),tuple(range(7,13)),tuple(range(13,19)))


def parse_encoder(path):
    t=[]; c=[]
    for line in path.read_text(errors="replace").splitlines():
        f=line.strip().split("\t")
        if len(f)==2 and f[0].isdigit() and f[1].isdigit(): t.append(int(f[0])/1000); c.append(int(f[1]))
    c=np.asarray(c,dtype=np.int64); c=np.where(c>2**31,c-2**32,c); return np.asarray(t),c-c[0]


def parse_ids(path):
    blocks={}; current=None; t=[]; y=[]
    def flush():
        if current is not None and t:
            tt=np.asarray(t)/1000; yy=(np.asarray(y)-y[0])*1e-6; blocks[current]=(tt,yy)
    for line in path.read_text(encoding="utf-8-sig",errors="replace").splitlines():
        if line.startswith("Date:"):
            flush(); current=int(line.rsplit(" ",1)[-1]); t=[]; y=[]; continue
        f=line.strip().split("\t")
        if len(f)>=2:
            try: t.append(float(f[0])); y.append(float(f[1]))
            except ValueError: pass
    flush(); return blocks


def edges(t,c,n):
    d=np.diff(medfilt(c.astype(float),kernel_size=5)); peaks,_=find_peaks(np.abs(d),height=EDGE_THRESHOLDS[n],distance=50)
    return t[peaks+1],np.sign(d[peaks])


def align_ids(t,y,edge_t,edge_sign):
    d=np.diff(y); amp=np.percentile(np.abs(y-np.median(y)),95); peaks,_=find_peaks(np.abs(d),height=.25*amp,distance=3)
    first=peaks[0]+1; sign=np.sign(d[peaks[0]]); return t+(edge_t-t[first]), y if sign==edge_sign else -y


def prepare_data(config):
    ids=parse_ids(V1/"IDSdata.txt"); datasets=[]
    for n in config["step_sizes"]:
        te,ce=parse_encoder(V1/("StepSize%d.csv"%n)); et,sg=edges(te,ce,n)
        ti,yi=align_ids(*ids[n],et[0],sg[0]); last_edge=min(len(et)-1,2*config["fit_cycles"]-1); end=et[last_edge]+1.0
        mask=(ti>=0)&(ti<=end); ti=ti[mask]; yi=yi[mask]; yi-=np.mean(yi[ti<et[0]]) if np.any(ti<et[0]) else yi[0]
        datasets.append({"n":n,"edge_t":et[:last_edge+1],"sign":sg[:last_edge+1],"t":ti,"y_um":yi,"end":end})
    return datasets


def decode(base,x):
    p=dict(base); values={name:np.exp(x[i]) for i,name in enumerate(VARIABLES)}
    p["T_d"]=base["T_d"]*values["T_d"]
    for port,prefix in (("way","F"),("nut","F"),("sb","T")):
        for stem in ("sigma0","sigma1","sigma2","vs"): p[stem+"_"+port]=base[stem+"_"+port]*values[stem+"_"+port]
        ckey=prefix+"c_"+port; skey=prefix+"s_"+port
        p[ckey]=base[ckey]*values[ckey]
        gap_name=prefix+"s_gap_"+port; p[skey]=p[ckey]+(base[skey]-base[ckey])*values[gap_name]
    return p


def structure(p):
    ell=p["L"]/(2*np.pi); kem=p["N_r"]*p["T_hold"]
    M=np.diag([p["I_m"],p["I_c"],p["I_s"],p["I_sb"],p["M_screw"],p["M_s"]])
    K=np.array([[p["k_c"]+kem,-p["k_c"],0,0,0,0],[-p["k_c"],p["k_c"]+p["k_s1"],-p["k_s1"],0,0,0],[0,-p["k_s1"],p["k_s1"]+p["k_s2"]+ell**2*p["k_nut"],-p["k_s2"],ell*p["k_nut"],-ell*p["k_nut"]],[0,0,-p["k_s2"],p["k_s2"],0,0],[0,0,ell*p["k_nut"],0,p["k_brg"]+p["k_nut"],-p["k_nut"]],[0,0,-ell*p["k_nut"],0,-p["k_nut"],p["k_nut"]]],float)
    C=np.array([[p["c_c"]+p["c_EM"],-p["c_c"],0,0,0,0],[-p["c_c"],p["c_c"]+p["c_s1"],-p["c_s1"],0,0,0],[0,-p["c_s1"],p["c_s1"]+p["c_s2"]+ell**2*p["c_nut"],-p["c_s2"],ell*p["c_nut"],-ell*p["c_nut"]],[0,0,-p["c_s2"],p["c_s2"],0,0],[0,0,ell*p["c_nut"],0,p["c_brg"]+p["c_nut"],-p["c_nut"]],[0,0,-ell*p["c_nut"],0,-p["c_nut"],p["c_nut"]]],float)
    beta=p["k_nut"]/(p["k_nut"]+p["k_brg"]); kap=1/(1/p["k_nut"]+1/p["k_brg"]); kch=1/(1/p["k_c"]+1/p["k_s1"]); nu=ell**2*kap/(kch+ell**2*kap); mu=p["k_s1"]/(p["k_c"]+p["k_s1"])
    T=np.array([[1,0],[1-mu*nu,mu*nu/ell],[1-nu,nu/ell],[1-nu,nu/ell],[-ell*beta*(1-nu),beta*(1-nu)],[0,1]],float)
    ports={"way":np.array([0,0,0,0,0,1.]),"nut":np.array([0,0,ell,0,1.,-1.]),"sb":np.array([0,0,0,1.,0,0])}
    return T.T@M@T,T.T@C@T,T.T@K@T,T.T@np.array([kem,0,0,0,0,0.]),{k:v@T for k,v in ports.items()}


def simulate(p,data):
    M,C,K,b,J=structure(p); Mi=np.linalg.inv(M); state=np.zeros(7); cursor=0.; travel=0.; out=np.empty(len(data["t"])); out[:]=np.nan
    boundaries=list(data["edge_t"])+[data["end"]]; signs=list(data["sign"])+[0]
    for end,sg in zip(boundaries,signs):
        sample=np.where((data["t"]>=cursor-1e-10)&(data["t"]<=end+1e-10))[0]
        eval_t=np.unique(np.clip(np.r_[cursor,data["t"][sample],end],cursor,end))
        theta=travel/(p["L"]/(2*np.pi))
        def rhs(_t,y):
            q=y[:2]; v=y[2:4]; reaction=np.zeros(2); zd=[]
            for i,port in enumerate(PORTS):
                vel=float(J[port]@v); pref="T" if port=="sb" else "F"; s0=p["sigma0_"+port]; s1=p["sigma1_"+port]; s2=p["sigma2_"+port]; fc=p[pref+"c_"+port]; fs=p[pref+"s_"+port]; vs=p["vs_"+port]; speed=np.sqrt(vel*vel+p["smooth_velocity_epsilon"]**2); g=fc+(fs-fc)*np.exp(-(vel/vs)**2); zdi=vel-s0*speed/g*y[4+i]; force=s0*y[4+i]+s1*zdi+s2*vel; reaction-=J[port]*force; zd.append(zdi)
            det=np.array([p["T_d"]*np.sin(4*p["N_r"]*q[0]),0]); acc=Mi@(b*theta+reaction-C@v-K@q-det); return np.r_[v,acc,zd]
        def jac(_t,y):
            q=y[:2]; v=y[2:4]; A=np.zeros((7,7)); A[:2,2:4]=np.eye(2); pos=-K.copy(); pos[0,0]-=4*p["N_r"]*p["T_d"]*np.cos(4*p["N_r"]*q[0]); vel_block=-C.copy()
            for i,port in enumerate(PORTS):
                jp=J[port]; vv=float(jp@v); pref="T" if port=="sb" else "F"; s0=p["sigma0_"+port]; s1=p["sigma1_"+port]; s2=p["sigma2_"+port]; fc=p[pref+"c_"+port]; fs=p[pref+"s_"+port]; vs=p["vs_"+port]; eps=p["smooth_velocity_epsilon"]; exp=np.exp(-(vv/vs)**2); g=fc+(fs-fc)*exp; gp=(fs-fc)*exp*(-2*vv/vs**2); speed=np.sqrt(vv*vv+eps*eps); sp=vv/speed; decay=s0*speed/g; dp=s0*(sp/g-speed*gp/g**2); dzdv=1-dp*y[4+i]; dzdz=-decay; dfdv=s1*dzdv+s2; dfdz=s0+s1*dzdz; vel_block-=dfdv*np.outer(jp,jp); A[2:4,4+i]=Mi@(-jp*dfdz); A[4+i,2:4]=dzdv*jp; A[4+i,4+i]=dzdz
            A[2:4,:2]=Mi@pos; A[2:4,2:4]=Mi@vel_block; return A
        if end>cursor+1e-10:
            sol=solve_ivp(rhs,(cursor,end),state,method="Radau",jac=jac,t_eval=eval_t,rtol=2e-5,atol=1e-8)
            if not sol.success: raise RuntimeError(sol.message)
            state=sol.y[:,-1]
            if len(sample): out[sample]=np.interp(data["t"][sample],sol.t,sol.y[1])*1e6
        travel+=sg*FULL_STEP_M/data["n"]; cursor=end
    if np.any(~np.isfinite(out)): out=np.interp(data["t"],data["t"][np.isfinite(out)],out[np.isfinite(out)])
    return out


def main():
    config=json.loads(CONFIG_FILE.read_text()); base=json.loads(PARAM_FILE.read_text())["parameters"]; locked={k:base[k] for k in STRUCTURAL}; data=prepare_data(config); rng=np.random.default_rng(config["random_seed"])
    lo=np.full(19,config["log_multiplier_bounds"][0]); hi=np.full(19,config["log_multiplier_bounds"][1]); lo[0],hi[0]=config["detent_log_multiplier_bounds"]
    history=[]; solutions=[]; cache={}; RUN_LOG.write_text("")
    scales={d["n"]:max(np.percentile(d["y_um"],95)-np.percentile(d["y_um"],5),FULL_STEP_M/d["n"]*1e6*.25) for d in data}
    def log_event(event,payload):
        record={"timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),"event":event}; record.update(payload)
        with RUN_LOG.open("a") as handle: handle.write(json.dumps(record)+"\n")
    log_event("SETUP",{"config":config,"variables":list(VARIABLES),"groups":[list(g) for g in GROUPS],"locked_structural":locked,"trajectory_scales_um":scales})
    def residual(x):
        key=tuple(np.round(x,10));
        if key in cache:return cache[key]
        p=decode(base,x); assert all(p[k]==locked[k] for k in STRUCTURAL)
        parts=[]
        for d in data: parts.append((simulate(p,d)-d["y_um"])/scales[d["n"]])
        parts.append(np.sqrt(config["regularization_weight"])*x); value=np.concatenate(parts); cache[key]=value; return value
    initial_file=HERE/config.get("initial_summary_file","pilot_optimization_summary.json")
    previous=np.asarray([json.loads(initial_file.read_text())["variables"][name] for name in VARIABLES]) if config.get("initialize_from_previous_best") else np.zeros(19)
    for start in range(config["number_of_starts"]):
        x=previous.copy() if start==0 else np.clip(previous+rng.normal(0,config["restart_log_perturbation"],19),lo,hi); start_cost=.5*np.dot(residual(x),residual(x)); began=time.time(); log_event("START",{"start":start,"initial_cost":start_cost,"x":x.tolist()})
        for cycle in range(config["alternating_passes"]):
            for group_index,group in enumerate(GROUPS):
                fixed=x.copy()
                def local(z): candidate=fixed.copy(); candidate[list(group)]=z; return residual(candidate)
                result=least_squares(local,x[list(group)],bounds=(lo[list(group)],hi[list(group)]),max_nfev=config["max_function_evaluations_per_group"],loss=config["robust_loss"],f_scale=config["robust_f_scale"],x_scale="jac",**config["optimizer_tolerances"])
                candidate=x.copy(); candidate[list(group)]=result.x
                old_cost=.5*np.dot(residual(x),residual(x)); candidate_cost=.5*np.dot(residual(candidate),residual(candidate))
                accepted=bool(candidate_cost <= old_cost)
                if accepted: x=candidate
                cost=min(old_cost,candidate_cost); history.append([start,cycle,group_index,cost,time.time()-began])
                CHECKPOINT.write_text(json.dumps({"start":start,"cycle":cycle,"group":group_index,"cost":cost,"x":x.tolist()},indent=2))
                print("start %d cycle %d group %d cost %.6g"%(start,cycle,group_index,cost),flush=True)
                log_event("GROUP_COMPLETE",{"start":start,"cycle":cycle,"group":group_index,"cost":cost,"candidate_cost":candidate_cost,"accepted":accepted,"nfev":result.nfev,"x":x.tolist()})
        solutions.append([start,start_cost,.5*np.dot(residual(x),residual(x),),time.time()-began]+list(x))
    best=min(solutions,key=lambda row:row[2]); xb=np.asarray(best[4:]); pb=decode(base,xb); predictions=[simulate(pb,d) for d in data]
    arrays={"variable_names":np.asarray(VARIABLES),"solutions":np.asarray(solutions),"history":np.asarray(history),"best_x":xb}
    for d,pred in zip(data,predictions): arrays["step%d_time_s"%d["n"]]=d["t"]; arrays["step%d_measured_um"%d["n"]]=d["y_um"]; arrays["step%d_simulated_um"%d["n"]]=pred; arrays["step%d_edges_s"%d["n"]]=d["edge_t"]
    np.savez_compressed(OUTPUT,**arrays)
    rmse={str(d["n"]):float(np.sqrt(np.mean((pred-d["y_um"])**2))) for d,pred in zip(data,predictions)}
    result={"best_start":int(best[0]),"initial_cost":best[1],"optimized_cost":best[2],"elapsed_s":best[3],"physical_rmse_um_by_step_size":rmse,"trajectory_scales_um":scales,"variables":dict(zip(VARIABLES,xb.tolist())),"optimized_parameters":{k:pb[k] for k in ["T_d"]+[v for v in pb if any(v.startswith(s) for s in ("sigma","Fc_","Fs_","Tc_","Ts_","vs_"))]},"locked_structural_parameters":locked,"config":config}
    SUMMARY.write_text(json.dumps(result,indent=2)+"\n"); print(json.dumps(result,indent=2))


if __name__=="__main__":main()

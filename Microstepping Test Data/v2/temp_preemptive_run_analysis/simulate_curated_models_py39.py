#!/usr/bin/env python3
"""Simulate all A2 commands with Rev 4 tangent, exact detent, and Rev 4.2 LuGre.

This file intentionally targets macOS system Python 3.9, where this machine's
working SciPy installation resides. It only writes an NPZ; plotting remains in
the main curation script.
"""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ESP = HERE.parent / "data" / "esp32_runs" / "full_campaign_marked_complete_20260828.txt"
P_BASE = ROOT / "Rev 4" / "model_parameters.json"
P_LUGRE = ROOT / "Rev 4" / "lugre_friction" / "Rev 4.2" / "model_parameters.json"
FULL_STEP = np.deg2rad(1.8)
N_Q = 6


def events():
    header = None; rows = []
    for line in ESP.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.strip().strip("\r")
        if text.startswith("timestamp_us,event,"):
            header = next(csv.reader([text])); continue
        if header and text and text[0].isdigit():
            fields = next(csv.reader([text]))
            if len(fields) == len(header): rows.append(dict(zip(header, fields)))
    return rows


def structural(p, tangent_detent=False):
    lead = p["L"] / (2 * np.pi); kem = p["N_r"] * p["T_hold"]
    M = np.diag([p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"]])
    K = np.array([
        [p["k_c"] + kem, -p["k_c"], 0, 0, 0, 0],
        [-p["k_c"], p["k_c"] + p["k_s1"], -p["k_s1"], 0, 0, 0],
        [0, -p["k_s1"], p["k_s1"] + p["k_s2"] + lead**2*p["k_nut"], -p["k_s2"], lead*p["k_nut"], -lead*p["k_nut"]],
        [0, 0, -p["k_s2"], p["k_s2"], 0, 0],
        [0, 0, lead*p["k_nut"], 0, p["k_brg"] + p["k_nut"], -p["k_nut"]],
        [0, 0, -lead*p["k_nut"], 0, -p["k_nut"], p["k_nut"]],
    ], dtype=float)
    if tangent_detent: K[0, 0] += 4*p["N_r"]*p["T_d"]
    C = np.array([
        [p["c_c"] + p["c_EM"], -p["c_c"], 0, 0, 0, 0],
        [-p["c_c"], p["c_c"] + p["c_s1"], -p["c_s1"], 0, 0, 0],
        [0, -p["c_s1"], p["c_s1"] + p["c_s2"] + lead**2*p["c_nut"], -p["c_s2"], lead*p["c_nut"], -lead*p["c_nut"]],
        [0, 0, -p["c_s2"], p["c_s2"], 0, 0],
        [0, 0, lead*p["c_nut"], 0, p["c_brg"] + p["c_nut"], -p["c_nut"]],
        [0, 0, -lead*p["c_nut"], 0, -p["c_nut"], p["c_nut"]],
    ], dtype=float)
    B = np.array([kem, 0, 0, 0, 0, 0], dtype=float)
    return M, C, K, B


def lugre_force(p, port, v, z):
    prefix = "T" if port == "sb" else "F"
    s0, s1, s2 = p["sigma0_"+port], p["sigma1_"+port], p["sigma2_"+port]
    fc, fs, vs = p[prefix+"c_"+port], p[prefix+"s_"+port], p["vs_"+port]
    speed = np.sqrt(v*v + p["smooth_velocity_epsilon"]**2)
    g = fc + (fs-fc)*np.exp(-(v/vs)**2)
    zdot = v - s0*speed/g*z
    return s0*z + s1*zdot + s2*v, zdot


def simulate(kind, edge_t, edge_pos, duration, p):
    M, C, K, B = structural(p, tangent_detent=(kind == "tangent"))
    Minv = np.diag(1/np.diag(M)); lead = p["L"]/(2*np.pi)
    jac = {"way": np.array([0,0,0,0,0,1.]),
           "nut": np.array([0,0,lead,0,1.,-1.]),
           "sb": np.array([0,0,0,1.,0,0])}
    state = np.zeros(15 if kind == "lugre" else 12)
    parts_t=[]; parts_y=[]; cursor=0.0; command=0.0
    boundaries = list(zip(edge_t, edge_pos)) + [(duration, edge_pos[-1] if len(edge_pos) else 0.0)]
    for end, next_command in boundaries:
        if end > cursor + 1e-10:
            count = max(2, int(np.ceil((end-cursor)/0.005))+1)
            tout = np.linspace(cursor, end, count)
            def rhs(_t, y):
                q=y[:6]; v=y[6:12]
                force = B*command - C.dot(v) - K.dot(q)
                zdots=[]
                if kind != "tangent": force[0] -= p["T_d"]*np.sin(4*p["N_r"]*q[0])
                if kind == "lugre":
                    for i, port in enumerate(("way","nut","sb")):
                        fp, zd = lugre_force(p, port, float(jac[port].dot(v)), y[12+i])
                        force -= jac[port]*fp; zdots.append(zd)
                acc=Minv.dot(force)
                return np.concatenate((v,acc,np.asarray(zdots)))
            sol=solve_ivp(rhs,(cursor,end),state,method="Radau",t_eval=tout,rtol=1e-5,atol=1e-8)
            if not sol.success: raise RuntimeError(sol.message)
            state=sol.y[:,-1]; parts_t.append(sol.t[:-1]); parts_y.append(sol.y[5,:-1])
        command = next_command
        cursor = end
    return np.concatenate(parts_t), np.concatenate(parts_y)*1e6


def main():
    rows=events(); pbase=json.loads(P_BASE.read_text())["parameters"]; plugre=json.loads(P_LUGRE.read_text())["parameters"]
    payload={}
    for n in (1,2,4,8,16,32):
        group=[r for r in rows if r["run_index"]=="1" and r["block"]=="A2" and r["label"]==("N%d_alternating"%n) and r["event"]=="PULSE_BATCH"][:12]
        tzero=int(group[0]["timestamp_us"])/1e6
        edge_t=np.array([int(r["timestamp_us"])/1e6-tzero for r in group])
        edge_pos=np.array([float(r["position_u16"])/16*FULL_STEP for r in group])
        duration=(int(group[-1]["timestamp_us"])/1e6-tzero)+0.55
        payload["N%d_command_t"%n]=np.r_[0,edge_t,duration]
        payload["N%d_command_um"%n]=np.r_[0,edge_pos*pbase["L"]/(2*np.pi)*1e6,edge_pos[-1]*pbase["L"]/(2*np.pi)*1e6]
        for kind,p in (("tangent",pbase),("detent",plugre),("lugre",plugre)):
            print("N%d %s"%(n,kind),flush=True)
            t,y=simulate(kind,edge_t,edge_pos,duration,p)
            payload["N%d_%s_t"%(n,kind)]=t; payload["N%d_%s_um"%(n,kind)]=y
    np.savez_compressed(HERE/"curated_model_simulations.npz",**payload)


if __name__ == "__main__": main()

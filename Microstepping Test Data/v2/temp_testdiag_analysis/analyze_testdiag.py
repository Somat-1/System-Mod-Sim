#!/usr/bin/env python3
"""Align and assess the QUICK12 IDS diagnostic recording."""
from pathlib import Path
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/"data"/"testDiag.csv"
PLOTS=HERE/"plots"

def load_ids():
    i=[]; y=[]
    for line in SOURCE.read_text(encoding="utf-8-sig",errors="replace").splitlines():
        f=line.split("\t")
        if len(f)>=2 and f[0].isdigit() and f[1].isdigit():
            i.append(int(f[0])); v=int(f[1]); y.append(v-2**32 if v>=2**31 else v)
    return np.asarray(i)*.001,np.asarray(y,dtype=float)

class Trace:
    def __init__(self): self.t=[0.0]; self.x=[0.0]; self.events=[]
    def dwell(self,d,label): self.t.append(self.t[-1]+d); self.x.append(self.x[-1]); self.events.append((self.t[-2],self.t[-1],label,"dwell"))
    def move(self,delta,d,label):
        a=self.t[-1]; self.t.append(a+d); self.x.append(self.x[-1]+delta); self.events.append((a,self.t[-1],label,"move"))
    def marker(self,amp,label,positive=True):
        s=1 if positive else -1; self.dwell(.5,label); self.move(s*amp,amp/50,label); self.dwell(.5,label); self.move(-s*amp,amp/50,label); self.dwell(.5,label)

def build_trace():
    q=Trace(); q.dwell(.075,"preflight")
    sections=[]
    pattern=[32,-16,8,-4,2,-1,1,-2,4,-8,16,-32]
    for name,scale,start_amp in (("half_native",1,20),("whole_endpoint",2,30)):
        q.dwell(.038,"configure"); q.marker(start_amp,f"{name}_start",True)
        a=q.t[-1]
        for n in (1,4,16):
            for sign in (1,-1):
                for _ in range(5):
                    p=sign*n*scale; q.move(p/2,abs(p)/500,f"{name}_A1_N{n}"); q.dwell(.4,f"{name}_A1_N{n}")
        sections.append((name,"A1",a,q.t[-1])); q.marker(10,f"{name}_A2_marker",False); a=q.t[-1]
        for n in (1,4,16):
            for _ in range(5):
                for sign in (1,-1):
                    p=sign*n*scale; q.move(p/2,abs(p)/500,f"{name}_A2_N{n}"); q.dwell(.4,f"{name}_A2_N{n}")
        sections.append((name,"A2",a,q.t[-1])); q.marker(15,f"{name}_B_marker",False); a=q.t[-1]
        for p0 in pattern:
            p=p0*scale; q.move(p/2,abs(p)/500,f"{name}_B"); q.dwell(.3,f"{name}_B")
        sections.append((name,"B",a,q.t[-1])); q.marker(20,f"{name}_E_marker",False); a=q.t[-1]
        for n in (1,4,16):
            for _ in range(5):
                p=n*scale; q.move(p/2,p/500,f"{name}_E_N{n}"); q.move(-p/2,p/500,f"{name}_E_N{n}"); q.dwell(1,f"{name}_E_N{n}")
        sections.append((name,"E",a,q.t[-1]))
    return q,sections

def smooth(y,n=21): return np.convolve(y,np.ones(n)/n,mode="same")

def main():
    PLOTS.mkdir(parents=True,exist_ok=True)
    t,y=load_ids(); q,sections=build_trace(); qt=np.asarray(q.t); qx=np.asarray(q.x)
    duration=qt[-1]; grid=np.arange(0,duration,.01); cmd=np.interp(grid,qt,qx)
    offsets=np.arange(0,max(.001,t[-1]-duration)+.0001,.01); scores=[]
    dx=np.diff(cmd); dx-=dx.mean(); nx=np.linalg.norm(dx)
    ys=smooth(y,31)
    for off in offsets:
        yy=np.interp(grid+off,t,ys); dy=np.diff(yy); dy-=dy.mean(); scores.append(abs(np.dot(dx,dy)/(nx*np.linalg.norm(dy)+1e-30)))
    offset=float(offsets[int(np.argmax(scores))]); measured=np.interp(grid+offset,t,ys)
    A=np.column_stack((cmd,np.ones_like(cmd),grid)); scale,intercept,drift=np.linalg.lstsq(A,measured,rcond=None)[0]
    fit=A@np.array([scale,intercept,drift]); resid=measured-fit
    pre=y[t<max(offset-.5,.5)]; noise_sigma=float(1.4826*np.median(np.abs(pre-np.median(pre))))
    metrics=[]
    for config,block,a,b in sections:
        m=(grid>=a)&(grid<=b); c=cmd[m]; z=measured[m]
        design=np.column_stack((c,np.ones(len(c)),grid[m])); s,i0,d=np.linalg.lstsq(design,z,rcond=None)[0]
        pred=design@np.array([s,i0,d]); corr=float(np.corrcoef(c,z)[0,1]) if np.std(c)>0 else 0
        metrics.append(dict(config=config,block=block,start_s=a,end_s=b,counts_per_full_step=float(s),correlation=corr,residual_rms_counts=float(np.sqrt(np.mean((z-pred)**2))),span_counts=float(np.percentile(z,95)-np.percentile(z,5))))
    result=dict(samples=len(y),recording_duration_s=float(t[-1]),command_duration_s=float(duration),aligned_start_s=offset,alignment_score=float(max(scores)),baseline_noise_sigma_counts=noise_sigma,global_counts_per_full_step=float(scale),global_residual_rms_counts=float(np.sqrt(np.mean(resid**2))),global_signal_span_counts=float(np.percentile(measured,99)-np.percentile(measured,1)),end_minus_start_counts=float(np.median(measured[-100:])-np.median(measured[:100])),sections=metrics)
    (HERE/"summary.json").write_text(json.dumps(result,indent=2))
    fig,(a0,a1)=plt.subplots(2,1,figsize=(16,9),sharex=True,constrained_layout=True)
    a0.plot(t,y,color="#222",lw=.35); a0.axvspan(offset,offset+duration,color="#f4a261",alpha=.18); a0.set(ylabel="EL5101 counts",title="testDiag complete recording and aligned QUICK12 interval"); a0.grid(alpha=.2)
    a1.plot(grid+offset,measured,color="#1769aa",lw=.8,label="measured (21 ms smooth)"); ax=a1.twinx(); ax.plot(grid+offset,cmd,color="#d94801",lw=1,label="commanded full-step position"); a1.set(xlabel="IDS time [s]",ylabel="EL5101 counts"); ax.set_ylabel("Commanded position [full steps]"); a1.grid(alpha=.2); a1.legend(a1.lines+ax.lines,[x.get_label() for x in a1.lines+ax.lines],loc="upper right")
    fig.savefig(PLOTS/"00_testdiag_alignment.png",dpi=180); plt.close(fig)
    fig,axs=plt.subplots(4,2,figsize=(17,13),constrained_layout=True)
    for ax0,(config,block,a,b),met in zip(axs.flat,sections,metrics):
        m=(grid>=a)&(grid<=b); tt=grid[m]-a; zz=measured[m]-np.median(measured[m][:max(1,int(.15*sum(m)))]); cc=cmd[m]-cmd[m][0]
        ax0.plot(tt,zz,color="#1769aa",lw=.8,label="measured"); axr=ax0.twinx(); axr.plot(tt,cc,color="#d94801",lw=1,label="command"); ax0.set_title(f"{config} {block}: r={met['correlation']:.3f}, gain={met['counts_per_full_step']:.2f} count/full-step"); ax0.set_ylabel("Relative counts"); axr.set_ylabel("Full steps"); ax0.grid(alpha=.2)
    fig.savefig(PLOTS/"01_testdiag_sections.png",dpi=180); plt.close(fig)
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()

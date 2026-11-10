#!/usr/bin/env python3
"""F1 companion illustrations: motion, observable, and counterfactual.
No measured friction data are used. All expected curves share one 30-element
quasi-static Maxwell-slip solver; grey alternatives are explicitly illustrative.
"""
from __future__ import annotations
import json
import textwrap
from pathlib import Path
from functools import lru_cache
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, brentq
from f1_sequence import CONFIG
from render_f1_sequence import style, INK, MUTED, GRID

HERE=Path(__file__).resolve().parent
OUT=HERE/"rendered"/"companion"
C=json.loads(CONFIG.read_text(encoding="utf-8"))
COL=C["colors"]; BETA=C["beta_0_um"]
S90_CFG=C["f0b_dependent"]["s90_um"]; DMIN=C["f0b_dependent"]["d_min_um"]
LIMITS=np.geomspace(.01,20,30)
GHOST="#9aa0a6"; SETUP="#b8bdc2"
# Fit the log-normal centre and width: fixed sigma=1 and centre=s90/3
# cannot also satisfy both saturation and initial stiffness. Preserve beta,
# sum(gamma)=20 and the configured illustrative s90 instead.
def weights(params):
    mu,sigma=params
    raw=np.exp(-.5*((np.log(LIMITS)-mu)/sigma)**2)
    return raw*(BETA/2)/(raw@LIMITS)
def calibration(params):
    g=weights(params)
    dx=brentq(lambda v: 2*np.sum(g*np.minimum(v,LIMITS))-.9*BETA,0,20)
    return [np.log(g.sum()/20),np.log((2*dx+.9*BETA)/S90_CFG)]
CAL=least_squares(calibration,[-4,2],bounds=([-12,.3],[4,6]),
                  xtol=1e-12,ftol=1e-12,gtol=1e-12)
G=weights(CAL.x)
assert np.linalg.norm(calibration(CAL.x))<1e-6

def advance(z,du):
    """Exact event-to-event equilibrium: du = dx + gamma . dz."""
    d=np.sign(du); remaining=du; dx_total=0.
    for _ in range(len(LIMITS)+2):
        if abs(remaining)<1e-13: break
        slipping=(np.abs(z)>=LIMITS-1e-12)&(np.sign(z)==d)
        k=1+G[~slipping].sum()
        room=np.where(slipping,np.inf,LIMITS-d*z)
        dx=d*min(abs(remaining)/k,float(room.min()))
        z=np.where(slipping,z,np.clip(z+dx,-LIMITS,LIMITS))
        remaining-=dx*k; dx_total+=dx
    assert abs(remaining)<1e-9
    return z,dx_total

def run(z0,u):
    z=np.array(z0,float); e=[-G@z]
    for du in np.diff(u):
        old=z.copy()
        z,dx=advance(z,float(du))
        assert abs(du-dx-G@(z-old))<1e-8
        e.append(-G@z)
    return np.asarray(u),np.array(e),z

def path(points,n=180):
    return np.r_[points[0],np.concatenate([np.linspace(a,b,n+1)[1:] for a,b in zip(points[:-1],points[1:])])]

def virgin(direction,travel):
    u=np.r_[0,direction*np.asarray(travel)]
    _,e,_=run(-direction*LIMITS,u)
    return np.abs(u),-direction*(e-e[0])

@lru_cache(maxsize=256)
def loop(stroke,n=180):
    h=stroke/2
    u,e,_=run(LIMITS,path([0,-h,h,-h,h,-h,h,-h],n))
    return ((u[5*n:6*n+1],e[5*n:6*n+1]),
            (u[6*n:7*n+1],e[6*n:7*n+1]))

def width_at_centre(stroke):
    (u,e),(v,f)=loop(float(stroke))
    return float(np.interp(0,v[::-1],f[::-1])-np.interp(0,u,e))

def finish(ax):
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(color=GRID,lw=.5); ax.set_axisbelow(True)
    ax.tick_params(length=3)

def badge(ax,x,y,n,colour=INK,offset=(0,0)):
    ax.annotate(str(n),xy=(x,y),xytext=offset,textcoords="offset points",
                ha="center",va="center",fontsize=6.5,color=colour,zorder=12,
                bbox=dict(boxstyle="circle,pad=.20",fc="white",ec=colour,lw=.8))

def icon(ax,points,colour,numbers=None,title="command",times=None):
    t=np.arange(len(points)) if times is None else np.asarray(times)
    ax.plot(t,points,color=SETUP,lw=2.4,zorder=1)
    ax.plot(t,points,color=colour,lw=1.2,zorder=2)
    if numbers:
        for i,num in numbers:
            badge(ax,t[i],points[i],num,colour,(0,7 if points[i]>=np.mean(points) else -7))
    ax.set(xlabel="time (schematic)",ylabel="$u$",xticks=[],yticks=[])
    ax.set_title(title,loc="left",fontsize=8,pad=10)
    ax.margins(x=.10,y=.35)
    ax.spines[["top","right"]].set_visible(False)

def title(ax,label,text):
    ax.set_title(f"({label})  {text}",loc="left",fontsize=9,pad=10)

def caption(fig,text):
    fig.text(.07,.025,textwrap.fill(text,145 if fig.get_figwidth()>9 else 114),
             fontsize=7.5,color=MUTED,va="bottom",linespacing=1.35)
    fig.text(.98,.985,"ILLUSTRATIVE | NOT MEASURED",ha="right",va="top",fontsize=7,color=MUTED)

def save(fig,name):
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f"{name}.png",dpi=220)   # PNG only; the PDFs went unused
    plt.close(fig)
    print(name,flush=True)

SERIF_RC={"font.family":"serif",
          "font.serif":["STIX Two Text","Times New Roman","DejaVu Serif"],
          "mathtext.fontset":"stix"}

def virgin_scaled(gain,direction,travel):
    """Virgin curve from the same solver with every element weight scaled.

    Scaling the weights moves the saturation level, and so the lost motion at
    large travel, while leaving the slip limits and therefore the shape of the
    pre-sliding bend where they are. That is what a direction asymmetry looks
    like: the same curve family with a different ceiling, not a different stage.
    """
    global G
    saved=G
    try:
        G=saved*gain
        return virgin(direction,travel)
    finally:
        G=saved

def figure_v1():
  with plt.rc_context(SERIF_RC):
    col=COL["V"]; s=np.geomspace(.01,100,1400)
    sd,ld=virgin(-1,s); su,lu=virgin(1,s)
    sa,la=virgin_scaled(.70,1,s)
    s90=float(np.interp(.9*BETA,ld,sd))
    dead=2*LIMITS.min()*(1+G.sum())
    ylo,yhi=.008,3.
    fig=plt.figure(figsize=(9,5.8))
    ic=fig.add_axes([.09,.72,.40,.18])
    # The wipe is compressed schematically, not drawn on the ramp's scale.
    icon(ic,[0,-1,0,0,-.25,-1.05],col,[(2,1),(3,2),(4,3),(5,4)],
         "(a)  wipe, then reverse")
    ic.lines[1].set_data([3,4,5],[0,-.25,-1.05])
    ic.text(.03,.02,"wipe: 3 mm",transform=ic.transAxes,fontsize=7,color=MUTED)
    # Each speed sits on the segment it describes.
    # Below and above respectively, to clear the badges at each segment end.
    ic.annotate(f"slow {C['v_slow_um_s']:g} µm/s",xy=(3.5,-.125),xytext=(0,-12),
                textcoords="offset points",fontsize=7,color=col,ha="center")
    ic.annotate(f"{C['v_fast_um_s']:g} µm/s",xy=(4.5,-.65),xytext=(9,6),
                textcoords="offset points",fontsize=7,color=col,ha="left")

    # The lower margin carries tick labels, the regime labels and the x-label
    # above a six-line caption, so the panel sits a little higher than (a).
    rect=[.09,.285,.86,.335]
    ax=fig.add_axes(rect)
    ax.axvspan(.01,dead,color=col,alpha=.09,lw=0)
    ax.axvspan(dead,s90,color=col,alpha=.045,lw=0)

    # l = s, everything commanded is lost. On log-log this is a straight
    # diagonal, so "stuck" reads as a line the curve lies along rather than as
    # a flat region sitting on the floor of a linear axis.
    ax.plot(s,s,color=GHOST,lw=.8,ls="--")
    dec_x=fig.get_figwidth()*rect[2]/np.log10(100/.01)
    dec_y=fig.get_figheight()*rect[3]/np.log10(yhi/ylo)
    ax.text(.023,.034,"carriage stuck:  $\\ell = s$",fontsize=7.5,color=MUTED,
            rotation=np.degrees(np.arctan2(dec_y,dec_x)),rotation_mode="anchor",
            ha="left",va="bottom")

    ax.plot(sd[1:],ld[1:],color=col,lw=3,alpha=.35)
    ax.plot(su[1:],lu[1:],color=col,lw=1.3,ls=(0,(5,3)))
    ax.plot(sa[1:],la[1:],color=GHOST,lw=1.2,ls=(0,(5,3)))

    # The two directions coincide in the model, so they are labelled on the
    # shared curve from opposite sides.
    ax.annotate("down (solid)",xy=(.16,float(np.interp(.16,sd,ld))),
                xytext=(14,-17),textcoords="offset points",fontsize=7.5,color=col,
                arrowprops=dict(arrowstyle="-",color=col,lw=.7))
    ax.annotate("up (dashed)",xy=(.40,float(np.interp(.40,su,lu))),
                xytext=(-8,21),textcoords="offset points",fontsize=7.5,color=col,
                ha="right",arrowprops=dict(arrowstyle="-",color=col,lw=.7))
    ax.annotate("if up $\\neq$ down: geometric\nleak or asymmetric friction",
                xy=(33,float(np.interp(33,sa,la))),xytext=(-6,-30),
                textcoords="offset points",fontsize=7,color=GHOST,ha="right",
                arrowprops=dict(arrowstyle="-",color=GHOST,lw=.7))

    ax.axhline(BETA,color=GHOST,lw=.8,ls=":")
    # Top left, the one empty corner: the right end carries badge 4.
    ax.text(.0115,BETA*1.13,f"sliding: fixed loss  β₀ = {BETA:.2f} µm",
            ha="left",fontsize=7.5)
    ax.axvline(s90,color=col,lw=.8,ls="--")
    ax.text(s90*1.18,.0105,"$s_{90}$\npre-sliding length",fontsize=7.5,
            color=MUTED,va="bottom")

    # 2 is the reversal itself, at zero travel, which a log axis cannot show;
    # it stays in the command sketch only.
    for x,n in [(C["v_slow_factor_of_s90"]*S90_CFG,3),(100,4)]:
        badge(ax,x,float(np.interp(x,sd,ld)),n,col,(0,9))
    for left,right,lab in [(.01,dead,"nearly stuck"),(dead,s90,"pre-sliding"),
                           (s90,100,"sliding")]:
        ax.text(np.sqrt(left*right),-.155,lab,ha="center",fontsize=7.5,
                color=MUTED,transform=ax.get_xaxis_transform())

    ax.set(xscale="log",yscale="log",xlim=(.01,100),ylim=(ylo,yhi),
           ylabel="lost motion $\\ell=s-\\Delta x$ [µm]")
    # Pushed clear of the regime labels under the axis.
    ax.set_xlabel("commanded travel since reversal $s$ [µm]",labelpad=15)
    ax.set_title("(b)",loc="left",fontsize=9,pad=10)
    finish(ax)
    caption(fig,"After reversal, almost all commanded travel is initially lost; the carriage then follows partly as contacts slip, and finally slides with fixed loss β₀. Both axes are logarithmic, so the three regimes are geometric: the curve runs along the ℓ = s diagonal while the carriage is stuck, bends away through pre-sliding, and flattens at β₀. Coincidence of the down and up curves is the expected symmetry check; the grey dashed curve shows what a failure would look like, an up curve saturating at a different level, and is illustrative rather than measured. Point 2 is the reversal itself, at zero travel, which a logarithmic axis cannot show, so it appears in the command sketch only. Motion timing and wipe depth are schematic. F0b s90 = 3 µm and d_min = 5 nm remain placeholders.")
    save(fig,"f1_V1_virgin_curve")

def main():
    style()
    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["DejaVu Sans"],
                         "axes.titlesize":9,"legend.frameon":False})
    for render in (figure_v1,): render()
    report=dict(elements=len(G),beta_um=BETA,initial_stiffness=float(G.sum()),
                initial_lost_motion_slope=float(G.sum()/(1+G.sum())),
                lognormal_centre_um=float(np.exp(CAL.x[0])),lognormal_sigma=float(CAL.x[1]),
                s90_um=S90_CFG,d_min_um=DMIN,illustrative=True,
                note="Both log-normal parameters fitted to saturation, stiffness and placeholder s90; not identified from data.")
    (OUT/"illustrative_model.json").write_text(json.dumps(report,indent=2)+"\n")
    print(report)
if __name__=="__main__":main()

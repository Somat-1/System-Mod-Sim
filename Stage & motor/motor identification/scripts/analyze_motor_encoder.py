#!/usr/bin/env python3
"""Render a measured-data rotary-encoder report.

The v1 radial runout is reported but not applied as angle calibration because
scale radius, phase registration, and a shared instrument clock are missing.
"""
from pathlib import Path
import csv, json, re, textwrap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parent.parent
MAIN=ROOT/'data'/'2026-10-13_motor_identification'
FINE=ROOT/'data'/'2026-10-13_microstep_resolution'
V1=ROOT/'v1_2026-10-02'
OUT=ROOT/'analysis'/'motor_encoder_report'
PDF=ROOT/'analysis'/'motor_encoder_analysis_report.pdf'
CPR=20000.; DPC=360/CPR; APC=DPC*3600; C2C=CPR/(200*256)
C={'navy':'#18324A','blue':'#2474A6','cyan':'#49A6B8','orange':'#D97706',
   'red':'#B83A3A','green':'#287A59','purple':'#7356A8','grey':'#65727E'}

def load_data():
    p=MAIN/'motor_identification_run.csv'
    a=np.loadtxt(p,delimiter=',',comments='#',skiprows=3,usecols=range(1,10),dtype=np.int64)
    names=('seq','t','action','event','command','encoder','index','current','mres')
    d={n:a[:,i] for i,n in enumerate(names)}
    with (MAIN/'plan.csv').open(newline='',encoding='utf-8') as f: plan=list(csv.DictReader(f))
    ends=[]
    pat=re.compile(r'^# ACTION_END,index=(\d+),tag=[^,]+,command_u256=(-?\d+),encoder=(-?\d+)')
    with p.open(encoding='utf-8') as f:
        for line in f:
            m=pat.match(line)
            if m:
                i,cmd,enc=map(int,m.groups()); ends.append((i,enc-cmd*C2C))
    return d,plan,np.asarray([x[1] for x in ends])

def settled(d,action):
    m=d['action']==action; cmd=d['command'][m]; enc=d['encoder'][m]
    br=np.r_[0,np.flatnonzero(np.diff(cmd)!=0)+1,len(cmd)]
    out=[]
    for a,b in zip(br[:-1],br[1:]):
        n=max(3,int(np.ceil((b-a)*.28))); z=enc[max(a,b-n):b]
        out.append((cmd[a],np.mean(z),np.std(z)))
    out=np.asarray(out); return out[np.argsort(out[:,0])].T

def pair(d,a,b):
    ca,ea,sa=settled(d,a); cb,eb,sb=settled(d,b)
    c,ia,ib=np.intersect1d(ca,cb,return_indices=True)
    return c,ea[ia],eb[ib],sa[ia],sb[ib]

def v1_result():
    p=V1/'dial_tracker'/'output'/'Video_Project_5_dial.csv'
    rows=[]
    with p.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(x for x in f if not x.startswith('#')):
            if r['flagged']=='0': rows.append((float(r['t_s']),float(r['delta_um'])))
    z=np.asarray(rows); ang=z[:,0]/60.13*360; th=np.radians(ang)
    A=np.c_[np.ones(len(th)),np.cos(th),np.sin(th)]; co=np.linalg.lstsq(A,z[:,1],rcond=None)[0]
    fit=A@co
    stats={'eccentricity_um':float(np.hypot(co[1],co[2])),
           'phase_deg':float(np.degrees(np.arctan2(co[2],co[1]))%360),
           'tir_um':float(np.ptp(z[:,1])),
           'residual_rms_um':float(np.sqrt(np.mean((z[:,1]-fit)**2)))}
    return ang,z[:,1],fit,stats

def static_result(d):
    cmd,ef,er,sf,sr=pair(d,15,16); ang=cmd/51200*360; ideal=cmd*C2C
    f,r=ef-ideal,er-ideal; mean=(f+r)/2; hys=(f-r)/2
    th=np.radians(ang); A=np.c_[np.ones(len(ang)),ang/360,np.cos(th),np.sin(th)]
    co=np.linalg.lstsq(A,mean,rcond=None)[0]; one=co[2]*np.cos(th)+co[3]*np.sin(th)
    amp=float(np.hypot(co[2],co[3])); phase=float(np.degrees(np.arctan2(co[3],co[2]))%360)
    y=mean[:-1].copy(); x=np.linspace(-1,1,len(y),endpoint=False); y-=np.polyval(np.polyfit(x,y,1),x)
    spec=2*np.abs(np.fft.rfft(y))/len(y)
    return {'cmd':cmd,'angle':ang,'f':f,'r':r,'mean':mean,'hys':hys,'one':one,
            'corrected':mean-one,'amp':amp,'phase':phase,'spec':spec,
            'noise':float(np.median(np.r_[sf,sr]))}

def current_result(d,a,b):
    cmd,ef,er,_,_=pair(d,a,b); ideal=cmd*C2C
    mean=((ef-ideal)+(er-ideal))/2; hys=(ef-er)/2; pc=np.mod(cmd,1024); ph=np.unique(pc)
    fm=np.array([np.mean(mean[pc==x]) for x in ph]); fh=np.array([np.mean(hys[pc==x]) for x in ph])
    return ph/1024*360,fm-np.mean(fm),fh

def responses(d,actions,duration=.3):
    grid=np.linspace(-.025,duration,326); traces=[]
    for action,sgn in actions.items():
        for edge in np.flatnonzero((d['action']==action)&(d['event']==sgn)):
            lo=np.searchsorted(d['t'],d['t'][edge]-30000); hi=np.searchsorted(d['t'],d['t'][edge]+duration*1e6)
            t=(d['t'][lo:hi]-d['t'][edge])*1e-6; y=d['encoder'][lo:hi].astype(float); pre=y[t<-.003]
            if len(pre)>2: traces.append(np.interp(grid,t,sgn*(y-np.median(pre))))
    return grid,np.asarray(traces)

def velocity_result(d,plan):
    out=[]
    for action in range(46,68):
        ix=np.flatnonzero(d['action']==action)
        # Velocity telemetry uses event=0; commanded-position changes mark the
        # actual step edges for these actions.
        edges=ix[np.r_[False,np.diff(d['command'][ix])!=0]]
        if len(edges)<8: continue
        span=edges[-1]-edges[0]; lo=edges[0]+int(.1*span); hi=edges[-1]-int(.1*span)
        t=(d['t'][lo:hi]-d['t'][lo])*1e-6; cmd=d['command'][lo:hi]*C2C; enc=d['encoder'][lo:hi]
        sc=np.polyfit(t,cmd,1)[0]; se=np.polyfit(t,enc,1)[0]; sgn=1 if sc>0 else -1
        base=np.median(d['encoder'][max(ix[0],edges[0]-600):edges[0]]-d['command'][max(ix[0],edges[0]-600):edges[0]]*C2C)
        out.append((abs(float(plan[action]['speed_fs_s'])),sgn,abs(se/sc),sgn*(np.median(enc-cmd)-base)))
    return np.asarray(out)

def reversal_result(d):
    out=[]
    for action in range(68,80):
        m=d['action']==action; cmd=d['command'][m]*C2C; enc=d['encoder'][m].astype(float)
        pp=abs(np.median(enc[cmd==cmd.max()])-np.median(enc[cmd==cmd.min()])); ratio=pp/np.ptp(cmd)
        out.append({'action':action,'amp':2**((action-68)%6),'freq':.2 if action<74 else 1.,
                    'ratio':ratio,'cmd':cmd-np.median(cmd),'enc':enc-np.median(enc)})
    return out

def fine_result():
    summary=[]
    with (FINE/'microstep_resolution_run_summary.csv').open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            summary.append({k:(int(v) if k=='mres' else float(v) if v else None) for k,v in r.items()})
    steps={}
    with (FINE/'microstep_resolution_run_steps.csv').open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            b=steps.setdefault(int(r['mres']),{'dir':[],'cmd':[],'pos':[]})
            b['dir'].append(r['direction']); b['cmd'].append(float(r['command_counts'])); b['pos'].append(float(r['settled_counts']))
    return summary,steps

def style():
    mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':9.2,'axes.titlesize':12,
        'axes.grid':True,'grid.color':'#DCE3E8','grid.linewidth':.65,
        'axes.edgecolor':'#50606D','legend.framealpha':.94,'legend.fontsize':8.2})

def page(title,kicker=''):
    f=plt.figure(figsize=(8.27,11.69)); title_size=16 if len(title)>45 else 19
    f.text(.07,.965,title,fontsize=title_size,weight='bold',color=C['navy'],va='top')
    f.text(.07,.935,kicker,fontsize=9.2,color=C['grey'],va='top')
    f.add_artist(plt.Line2D([.07,.93],[.918,.918],color=C['cyan'],lw=2.2)); return f

def para(f,x,y,text,width=92,size=9.3,color='#1F2933'):
    f.text(x,y,textwrap.fill(text,width),ha='left',va='top',fontsize=size,color=color,linespacing=1.35)

def box(f,x,y,w,h,head,text,color='blue',size=9):
    col=C[color]; f.add_artist(mpl.patches.FancyBboxPatch((x,y),w,h,transform=f.transFigure,
        boxstyle='round,pad=.012,rounding_size=.012',facecolor=mpl.colors.to_rgba(col,.08),edgecolor=col,lw=1))
    f.text(x+.018,y+h-.025,head,color=col,weight='bold',fontsize=10,va='top')
    para(f,x+.018,y+h-.062,text,max(26,int(w*112)),size)

def save(pdf,f,n,name):
    f.text(.07,.026,'Rotary encoder measurement report · generated from recorded data',color=C['grey'],fontsize=7.5)
    f.text(.93,.026,f'{n:02d}',color=C['grey'],fontsize=8,ha='right'); pdf.savefig(f,dpi=170)
    f.savefig(OUT/'plots'/name,dpi=170); plt.close(f)

def build():
    style(); (OUT/'plots').mkdir(parents=True,exist_ok=True); PDF.parent.mkdir(parents=True,exist_ok=True)
    d,plan,ends=load_data(); st=static_result(d); va,vd,vfit,vs=v1_result()
    currents={180:current_result(d,18,19),270:current_result(d,21,22),360:current_result(d,24,25)}
    tf,rf=responses(d,{26:1,27:-1}); tm,rm=responses(d,{28:1,29:-1})
    vel=velocity_result(d,plan); rev=reversal_result(d); fine,steps=fine_result()
    duration=(d['t'][-1]-d['t'][0])*1e-6; gaps=int(np.sum(np.clip(np.diff(d['seq'])-1,0,None)))
    jumps=np.diff(ends,prepend=ends[0]); ji=int(np.argmax(np.abs(jumps))); final=float(ends[-1])
    metrics={'samples':len(d['seq']),'duration_s':duration,'sequence_gaps':gaps,
      'final_residual_counts':final,'final_residual_deg':final*DPC,
      'largest_boundary_jump':{'action':ji,'label':plan[ji]['label'],'counts':float(jumps[ji]),'degrees':float(jumps[ji]*DPC)},
      'static_map':{'raw_mean_pp_counts':float(np.ptp(st['mean'])),'raw_mean_pp_deg':float(np.ptp(st['mean'])*DPC),
        'hysteresis_pp_deg':float(np.ptp(st['hys'])*DPC),'one_per_rev_amplitude_counts':st['amp'],
        'one_per_rev_amplitude_deg':st['amp']*DPC,'one_per_rev_phase_deg':st['phase'],
        'diagnostic_corrected_pp_deg':float(np.ptp(st['corrected'])*DPC),'settled_noise_median_counts':st['noise']},
      'v1_radial_runout':vs,'v1_applied_to_angle_data':False,
      'v1_nonapplication_reason':'Missing encoder scale radius, phase registration, and common clock.','fine_microstep_summary':fine}
    (OUT/'motor_encoder_metrics.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')

    with PdfPages(PDF) as pdf:
        n=1
        f=plt.figure(figsize=(8.27,11.69)); f.patches.append(mpl.patches.Rectangle((0,.72),1,.28,transform=f.transFigure,color=C['navy'],zorder=-1))
        f.text(.075,.925,'ROTARY ENCODER',color='#8FD3E1',fontsize=12,weight='bold')
        f.text(.075,.865,'Measured-data analysis',color='white',fontsize=29,weight='bold')
        f.text(.075,.815,'Motor identification campaign · 13 October 2026',color='white',fontsize=13)
        f.text(.075,.765,'A reader for students: evidence, calibration limits, and interpretation',color='#D6E6EE',fontsize=10.5)
        box(f,.075,.56,.39,.115,'WHAT WAS RECORDED',f"{len(d['seq']):,} samples over {duration/60:.1f} min, 83 actions, 20,000 counts/rev, {gaps} sequence gaps.",'blue')
        box(f,.535,.56,.39,.115,'RUN INTEGRITY',f"Final residual {final:+.0f} counts ({final*DPC:+.2f}°). Largest boundary change at action {ji}: {plan[ji]['label']}.",'red')
        box(f,.075,.385,.39,.125,'STATIC MAP',f"Direction-pair mean spans {np.ptp(st['mean']):.1f} counts ({np.ptp(st['mean'])*DPC:.2f}°); fitted 1/rev amplitude {st['amp']:.1f} counts ({st['amp']*DPC:.3f}°).",'green')
        box(f,.535,.385,.39,.125,'ECCENTRICITY EVIDENCE',f"V1 dial result: {vs['eccentricity_um']:.1f} µm eccentricity, {vs['tir_um']:.1f} µm TIR. It is not applied as angle calibration.",'orange')
        para(f,.075,.315,'Central lesson. Encoder resolution is not absolute accuracy. Quantisation, radial runout, torque harmonics, friction, driver state, and lost position can occupy the same error trace. This report separates observations from inferences.',105,10.5)
        f.text(.075,.14,'Reproduce with scripts/analyze_motor_encoder.py',color=C['grey'],fontsize=8.5)
        save(pdf,f,n,'01_cover.png'); n+=1

        f=page('What is measured—and what can be calibrated','The calibration decision is part of the result, not hidden preprocessing.')
        ax=f.add_axes([.08,.50,.84,.34]); ax.plot(va,vd,color=C['blue'],lw=1.2,label='v1 dial displacement'); ax.plot(va,vfit,color=C['orange'],lw=2,label='best 1/rev sinusoid')
        ax.set(xlim=(0,365),xlabel='rotation angle inferred from dial period [deg]',ylabel='radial displacement [µm]'); ax.legend(); ax.set_title('Independent v1 radial-runout measurement')
        ax.text(.02,.06,f"eccentricity {vs['eccentricity_um']:.1f} µm · TIR {vs['tir_um']:.1f} µm · phase {vs['phase_deg']:.0f}°",transform=ax.transAxes,color=C['navy'],weight='bold')
        box(f,.08,.31,.84,.13,'V1 RADIAL RESULT: NOT NUMERICALLY APPLIED TO ENCODER ANGLE','Centre offset e can create angular error proportional to e/R, where R is effective scale radius. R was not recorded, the high spot was not registered to encoder zero, and the instruments had no common clock. Applying 159.5 µm would manufacture an unknown amplitude and phase.','orange',9.1)
        box(f,.08,.13,.40,.12,'WHAT IS USED','Raw encoder counts, commanded phase, direction pairs, action boundaries, and harmonic decomposition of the measured static map.','green')
        box(f,.52,.13,.40,.12,'SENSITIVITY VIEW','Later pages subtract a fitted 1/rev term from the same static map, explicitly labeled diagnostic—not v1 calibration or ground truth.','purple')
        save(pdf,f,n,'02_calibration.png'); n+=1

        f=page('Complete campaign overview','Command and encoder use the same 20,000-count-per-revolution scale.')
        t=(d['t']-d['t'][0])*1e-6/60; take=max(1,len(t)//45000)
        a=f.add_axes([.08,.53,.84,.32]); a.plot(t[::take],d['command'][::take]*C2C,color='#AAB7C2',lw=.8,label='command'); a.plot(t[::take],d['encoder'][::take],color=C['blue'],lw=.65,label='encoder')
        a.set(ylabel='position [encoder counts]',xlim=(t[0],t[-1])); a.ticklabel_format(axis='y',style='sci',scilimits=(0,0)); a.legend(ncol=2); a.set_title('All recorded motion')
        b=f.add_axes([.08,.18,.84,.27],sharex=a); err=d['encoder']-d['command']*C2C; b.plot(t[::take],err[::take]*DPC,color=C['red'],lw=.7); b.axhline(0,color='#222',lw=.8)
        b.set(xlabel='elapsed time [min]',ylabel='encoder − command [deg]'); b.set_title('Raw residual: dynamics, direction effects, and lost position')
        for q in (261.36,1221.36,1417.68,1563.58,1857.58,2157.58): a.axvline(q/60,color=C['grey'],lw=.5); b.axvline(q/60,color=C['grey'],lw=.5)
        para(f,.08,.115,'Excursions during motion are expected lag and settling. Persistent changes after an action ends instead indicate a changed reference or lost position. The next page separates these cases.',112)
        save(pdf,f,n,'03_overview.png'); n+=1

        f=page('Position closure at every action boundary','Boundary residual = encoder minus command after the programmed action finishes.')
        a=f.add_axes([.08,.42,.84,.43]); actions=np.arange(len(ends)); cols=np.where(np.abs(jumps)>=100,C['red'],C['blue'])
        a.vlines(actions,0,ends*DPC,color=cols,alpha=.65,lw=1); a.scatter(actions,ends*DPC,c=cols,s=18); a.axhline(0,color='#222',lw=.9)
        a.set(xlabel='action index',ylabel='end residual [deg]',xlim=(-1,83)); a.set_title('One position reference is not maintained through all blocks')
        a.annotate(f"largest change\naction {ji}: {jumps[ji]:+.0f} counts",xy=(ji,ends[ji]*DPC),xytext=(ji-25,4.8),arrowprops={'arrowstyle':'->','color':C['red']},color=C['red'],weight='bold')
        box(f,.08,.24,.84,.11,'KEY INTEGRITY FINDING',f"Final residual {final:+.0f} counts = {final/100:+.2f} full steps = {final*DPC:+.2f}°. Largest change: action {ji} ({plan[ji]['label']}), {jumps[ji]:+.0f} counts. Later blocks require local interpretation.",'red')
        para(f,.08,.16,'The static map (actions 15–16) closes to 1 count and supports direction-paired analysis. The closing reference inherits offsets accumulated during later tests.',112)
        save(pdf,f,n,'04_closure.png'); n+=1

        f=page('Quasi-static angle map: raw and sensitivity-corrected','One shaft turn at native 1/16 microsteps; settled tail of every 150 ms dwell.')
        a=f.add_axes([.08,.55,.84,.29]); a.plot(st['angle'],st['f']*DPC,color=C['blue'],lw=.8,label='forward'); a.plot(st['angle'],st['r']*DPC,color=C['orange'],lw=.8,label='reverse')
        a.plot(st['angle'],st['mean']*DPC,color=C['navy'],lw=1.8,label='direction-pair mean'); a.plot(st['angle'],st['one']*DPC,color=C['purple'],lw=1.6,ls='--',label='fitted 1/rev')
        a.set(xlim=(0,360),ylabel='encoder angle error [deg]'); a.legend(ncol=2); a.set_title('Raw measured error')
        b=f.add_axes([.08,.18,.84,.28],sharex=a); b.plot(st['angle'],st['corrected']*DPC,color=C['purple'],lw=1,label='mean minus fitted 1/rev'); b.plot(st['angle'],st['hys']*DPC,color=C['green'],lw=.9,label='half forward−reverse')
        b.axhline(0,color='#222',lw=.8); b.set(xlabel='commanded mechanical angle [deg]',ylabel='angle [deg]'); b.legend(loc='lower right')
        b.text(.01,.98,'DIAGNOSTIC ONLY — FITTED FROM THIS MAP; V1 NOT APPLIED',transform=b.transAxes,va='top',color=C['purple'],weight='bold',fontsize=7.2,bbox={'fc':'white','ec':C['purple'],'boxstyle':'round'})
        para(f,.08,.115,f"Raw mean span {np.ptp(st['mean'])*DPC:.3f}°. Fitted 1/rev amplitude {st['amp']*DPC:.3f}° at phase {st['phase']:.1f}°. Subtraction leaves {np.ptp(st['corrected'])*DPC:.3f}° peak-to-peak, but also removes any real 1/rev motor error.",115)
        save(pdf,f,n,'05_static_map.png'); n+=1

        f=page('Spatial harmonic content of the static map','Cycles per mechanical revolution after removing offset and linear closure trend.')
        a=f.add_axes([.08,.40,.84,.44]); orders=np.arange(len(st['spec'])); amp=st['spec']*APC; a.semilogy(orders[1:251],amp[1:251],color=C['blue'],lw=.9)
        for order,col in ((1,'orange'),(2,'green'),(50,'purple'),(100,'red'),(200,'navy')):
            a.scatter(order,amp[order],color=C[col],s=34,zorder=4); a.annotate(f'{order}/rev\n{amp[order]:.0f} arcsec',(order,amp[order]),xytext=(5,8),textcoords='offset points',color=C[col],fontsize=8)
        a.set(xlim=(0,250),xlabel='spatial order [cycles/rev]',ylabel='amplitude [arcsec, log scale]'); a.set_title('Direction-pair mean')
        box(f,.08,.19,.84,.14,'HOW TO READ THIS','Order 1 is consistent with eccentric mounting but not unique to it. Orders 50 and 100 correspond to 50 electrical cycles/rev and its second harmonic; order 200 coincides with full-step positions and detent/step structure. Periodicity does not prove mechanism.','blue')
        save(pdf,f,n,'06_harmonics.png'); n+=1

        f=page('Current dependence over the electrical cycle','Four cycles folded and averaged; forward/reverse means suppress directional friction.')
        a=f.add_axes([.08,.51,.84,.33]); pal={180:'orange',270:'green',360:'blue'}
        for current,(ph,mn,hy) in currents.items(): a.plot(ph,mn*APC,color=C[pal[current]],lw=1.35,label=f'{current} mA')
        a.set(xlim=(0,360),ylabel='zero-mean angle error [arcsec]'); a.legend(title='requested current'); a.set_title('Repeatable electrical-cycle waveform')
        b=f.add_axes([.08,.18,.84,.23]); levels=sorted(currents); x=np.arange(3); rp=[np.ptp(currents[q][1])*APC for q in levels]; hp=[np.ptp(currents[q][2])*APC for q in levels]
        b.bar(x-.18,rp,.36,color=C['blue'],label='mean-error p-p'); b.bar(x+.18,hp,.36,color=C['green'],label='directional half-difference p-p'); b.set_xticks(x,[f'{q} mA' for q in levels]); b.set(ylabel='peak-to-peak [arcsec]'); b.legend(ncol=2)
        para(f,.08,.105,'Detent-torque angle error is expected to shrink roughly as current rises; phase-current waveform distortion is more nearly current-independent. This sweep compares trends but does not identify either cause uniquely.',112)
        save(pdf,f,n,'07_current.png'); n+=1

        f=page('Step response at the zero-position test point','Repeated excitation edges aligned at t = 0; shading is the central 80%.')
        a=f.add_axes([.08,.54,.84,.28]); lo,hi=np.percentile(rf,[10,90],axis=0); a.fill_between(tf*1000,lo,hi,color=C['blue'],alpha=.18); a.plot(tf*1000,np.mean(rf,axis=0),color=C['blue'],lw=1.6,label=f'mean of {len(rf)}')
        a.axhline(100,color='#222',ls='--',label='ideal full step'); a.axvline(0,color=C['grey'],lw=.8); a.set(xlim=(-25,300),ylabel='direction-normalised motion [counts]'); a.legend(); a.set_title('One full step')
        b=f.add_axes([.08,.19,.84,.25],sharex=a); lo,hi=np.percentile(rm,[10,90],axis=0); b.fill_between(tm*1000,lo,hi,color=C['purple'],alpha=.18); b.plot(tm*1000,np.mean(rm,axis=0),color=C['purple'],lw=1.6,label=f'mean of {len(rm)}')
        b.axhline(6.25,color='#222',ls='--',label='ideal 1/16 step'); b.axvline(0,color=C['grey'],lw=.8); b.set(xlabel='time from command edge [ms]',ylabel='motion [counts]'); b.legend(); b.set_title('One native 1/16 microstep')
        full_level=float(np.mean(rf[:,tf>.2])); micro_level=float(np.mean(rm[:,tm>.2]))
        para(f,.08,.11,f'One encoder count is 0.018°. The ensemble settles near {full_level:.0f} counts for the nominal 100-count full step and {micro_level:.1f} counts for the nominal 6.25-count microstep. These include the driver-state transition and should not be treated as a simple linear transfer function.',112)
        save(pdf,f,n,'08_step_response.png'); n+=1

        f=page('Constant-rate tracking','Linear fits use the central 80% between first and last step edges.')
        a=f.add_axes([.08,.53,.84,.29])
        for sgn,mark,col in ((1,'o','blue'),(-1,'s','orange')):
            q=vel[vel[:,1]==sgn]; a.semilogx(q[:,0],100*q[:,2],marker=mark,color=C[col],label=('+' if sgn>0 else '−')+' direction')
        a.axhline(100,color='#222',ls='--'); a.set(ylabel='measured / commanded speed [%]'); a.legend(); a.set_title('Mean speed gain')
        b=f.add_axes([.08,.20,.84,.24],sharex=a)
        for sgn,mark,col in ((1,'o','blue'),(-1,'s','orange')):
            q=vel[vel[:,1]==sgn]; b.semilogx(q[:,0],q[:,3]*DPC,marker=mark,color=C[col])
        b.axhline(0,color='#222',lw=.8); b.set(xlabel='commanded speed [full steps/s]',ylabel='direction-normalised lag [deg]'); b.set_title('Local lag relative to pre-motion error')
        para(f,.08,.115,'Speed gain tests count accumulation rate. A local baseline prevents earlier lost position from masquerading as velocity lag; it does not repair global reference changes.',112)
        save(pdf,f,n,'09_velocity.png'); n+=1

        f=page('Small triangular reversals and friction deadband','Response measured between samples at commanded positive and negative extrema.')
        a=f.add_axes([.08,.53,.84,.29])
        for freq,mark,col in ((.2,'o','blue'),(1.,'s','orange')):
            q=[x for x in rev if x['freq']==freq]; a.plot([x['amp'] for x in q],[100*x['ratio'] for x in q],marker=mark,color=C[col],label=f'{freq:g} Hz')
        a.axhline(100,color='#222',ls='--'); a.set_xscale('log',base=2); a.set_xticks([1,2,4,8,16,32],['1','2','4','8','16','32']); a.set(ylabel='measured p-p / command [%]'); a.legend(); a.set_title('Motion transmitted through reversal')
        b=f.add_axes([.08,.19,.84,.25])
        for action,col,label in ((73,'blue','32 µsteps · 0.2 Hz'),(79,'red','32 µsteps · 1 Hz')):
            q=next(x for x in rev if x['action']==action); take=max(1,len(q['cmd'])//5000); b.plot(q['cmd'][::take],q['enc'][::take],'.',ms=1.2,alpha=.38,color=C[col],label=label)
        b.set(xlabel='command relative to centre [counts]',ylabel='encoder relative to median [counts]'); b.legend(); b.set_title('Largest loops; action 79 includes a large within-block reference change')
        para(f,.08,.105,f"A closed narrow loop indicates repeatable following. Action 79 must not be reduced to one gain number: its {jumps[79]:+.0f}-count boundary change is qualitatively different from ordinary hysteresis.",112)
        save(pdf,f,n,'10_reversal.png'); n+=1

        f=page('Fine-microstep resolution campaign','Separate one-turn forward/back tests at 1/32, 1/64, and 1/128.')
        a=f.add_axes([.08,.47,.84,.37]); pal={32:'red',64:'green',128:'blue'}
        for mres,q in steps.items():
            dirs=np.asarray(q['dir']); cmds=np.asarray(q['cmd']); pos=np.asarray(q['pos'])
            for name,ls in (('fwd','-'),('back','--')):
                m=dirs==name; cmd=cmds[m]; p=pos[m]; order=np.argsort(cmd); cmd,p=cmd[order],p[order]; e=p-cmd; e-=e[0]; take=max(1,len(cmd)//2500)
                a.plot(cmd[::take]/CPR*360,e[::take]*DPC,ls=ls,lw=.8,color=C[pal[mres]],label=f'1/{mres} {name}')
        a.set(xlabel='commanded angle [deg]',ylabel='relative position error [deg]',xlim=(0,360)); a.legend(ncol=3,fontsize=7.5); a.set_title('Within-pass error shape; every pass re-zeroed at start')
        fs={x['mres']:x for x in fine}
        box(f,.08,.25,.40,.13,'1/64 AND 1/128',f"Mean steps {fs[64]['mean_step_counts']:.3f} and {fs[128]['mean_step_counts']:.3f} counts. Single-step responses {fs[64]['pm1_response_pct']:.0f}% and {fs[128]['pm1_response_pct']:.0f}% of ideal. Smooth 1/rev structure is consistent with eccentricity, not proof.",'green',8.7)
        box(f,.52,.25,.40,.13,'1/32 DATA QUALITY WARNING','The encoder accumulated +20,335 and −20,572 counts for nominal ±20,000 turns. Gradual over-count points to ringing or coupling slip. Repeat before using this trace.','red',8.7)
        para(f,.08,.16,'One encoder count is 64.8 arcsec. An ideal 1/128 microstep is only 0.781 count, so integer quantisation dominates individual increments; averaging is essential.',112)
        save(pdf,f,n,'11_fine_microsteps.png'); n+=1

        f=page('Conclusions, limitations, and a defensible next calibration','Concise interpretation for laboratory discussion.')
        box(f,.08,.73,.84,.12,'1 · WHAT THE DATA SUPPORT',f"The direction-paired static map resolves a {st['amp']*DPC:.3f}°-amplitude 1/rev term plus electrical/full-step harmonics. Average 1/64 and 1/128 motion is near nominal although individual microsteps are quantised and partly effective.",'green')
        box(f,.08,.55,.84,.12,'2 · WHAT THE DATA DO NOT SUPPORT','V1 radial runout cannot become absolute angular calibration with present metadata. A 1/rev harmonic is not uniquely eccentricity: the motor can contain once-per-rev structure. Later blocks do not share the opening reference.','orange')
        box(f,.08,.37,.84,.12,'3 · RECOMMENDED CALIBRATION EXPERIMENT','Record encoder angle and dial displacement from a common trigger over several slow turns. Register index to dial high spot; record scale radius and readhead geometry; repeat after remounting. Prefer an independent angular reference.','blue')
        para(f,.08,.27,'Reproducibility. Static values use the final 28% of each dwell. Direction samples are matched at identical commands. Spectra remove offset and linear closure trend. Step responses align marked edges. Velocity fits use the central 80%. No eccentricity correction is applied to source data.',112,9.2)
        para(f,.08,.16,'Sources: the 2026-10-13 motor-identification run, fine-microstep run, and v1 dial CSV. Machine-readable metrics and page images are stored in analysis/motor_encoder_report/.',112,8.6,C['grey'])
        save(pdf,f,n,'12_conclusions.png')

    print(f'Wrote {PDF}'); print(json.dumps(metrics,indent=2))

if __name__=='__main__': build()

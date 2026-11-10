"""Physical invariants for the illustrative companion-figure solver."""
import numpy as np
import render_f1_companions as r

assert len(r.G)==30 and np.all(r.G>0)
assert np.isclose(r.G.sum(),20)
assert np.isclose(2*r.G@r.LIMITS,r.BETA)
s=np.geomspace(1e-6,100,1000)
x,down=r.virgin(-1,s); _,up=r.virgin(1,s)
assert np.allclose(down,up,atol=1e-10)
assert np.all(np.diff(down)>=-1e-10)
assert np.all(down<=x+1e-10)
assert np.isclose(down[-1],r.BETA)
assert np.isclose(down[1]/s[0],20/21)
assert abs(np.interp(.9*r.BETA,down,x)-r.S90_CFG)<.003
widths=[r.width_at_centre(v) for v in (.1,.9,3.,10.,100.)]
assert np.all(np.diff(widths)>=-1e-10)
assert min(widths)>=0 and max(widths)<=r.BETA+1e-10
for stroke in (.9,3.,100.):
    (u,e),(v,f)=r.loop(stroke)
    d=u-u[0]
    _,initial,_=r.run(np.zeros_like(r.LIMITS),d/2)
    assert np.max(abs(e-(e[0]+2*initial)))<1e-8
    assert np.max(abs(f-(f[0]-2*initial)))<1e-8
outer=3*r.S90_CFG; inner=.4*outer;n=150
u,e,_=r.run(r.LIMITS,r.path([-outer/2,outer/2,0,inner,0,-outer/2],n))
assert abs(e[2*n]-e[4*n])<1e-8
_,e1,z1=r.run(r.LIMITS,np.array([0.,-7.]))
_,e2,z2=r.run(r.LIMITS,np.linspace(0.,-7.,1000))
assert np.allclose(z1,z2,atol=1e-10)
assert abs(e1[-1]-e2[-1])<1e-10
print("PASS: equilibrium, calibration, symmetry, saturation, monotonic width, Masing, memory, subdivision invariance.")

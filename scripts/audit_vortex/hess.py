import numpy as np
from vstruct import *
from scipy.optimize import minimize
p=B+'rb_c0p3_hxy/leg_0003/out/field.dat'
for i,c in [(2970,429.4)]:
    t,phi=read(p,i); pz=phi.mean(axis=2)
    def H(X,xi):
        e=np.exp(1j*kx*X)
        Psi=np.sum(pz[:,0]*e).real/Cxy
        v=2*np.real(np.sum(pz[:,1]*e)*np.exp(1j*ky*xi))/Cxy
        return Psi+v-c*X
    # search fixed points near x~660 (crit) : grad=0
    def g2(z):
        X,xi=z; h=1.0
        hx=(H(X+h,xi)-H(X-h,xi))/2/h; hy=(H(X,xi+h)-H(X,xi-h))/2/h
        return hx**2+(hy)**2*1e0
    Ly=2*np.pi/ky
    best=[]
    for X0 in np.linspace(450,900,10):
        for xi0 in np.linspace(0,Ly,12,endpoint=False):
            r=minimize(g2,[X0,xi0],method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-14,maxiter=4000))
            best.append((r.fun,r.x[0],r.x[1]%Ly))
    best=sorted(best)
    seen=[]
    for f,X,xi in best:
        if f>1e-6: continue
        if any(abs(X-a)<5 and abs(((xi-b+Ly/2)%Ly)-Ly/2)<20 for a,b in seen): continue
        seen.append((X,xi)); h=2.0
        Hxx=(H(X+h,xi)-2*H(X,xi)+H(X-h,xi))/h**2
        Hyy=(H(X,xi+h)-2*H(X,xi)+H(X,xi-h))/h**2
        Hxy=(H(X+h,xi+h)-H(X+h,xi-h)-H(X-h,xi+h)+H(X-h,xi-h))/4/h**2
        det=Hxx*Hyy-Hxy**2
        print(f"t={t:.1f} fixed pt X={X:.1f} xi={xi:.1f} |grad|^2={f:.1e} det={det:.3e} -> {'elliptic w_tr=%.4f'%np.sqrt(det) if det>0 else 'hyperbolic'}")

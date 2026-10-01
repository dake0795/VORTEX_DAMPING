import numpy as np
from vstruct import *
dx=lx/nx
def times_idx(path):
    return times(path)
runs={'leg3':B+'rb_c0p3_hxy/leg_0003/out/field.dat','leg4':B+'rb_c0p3_hxy/leg_0004/out/field.dat','leg5':B+'rb_c0p3_hxy/leg_0005/out/field.dat',
      'leg6':B+'rb_c0p3_hxy/leg_0006/out/field.dat','a03':B+'rb_c0p3_hxy_ampltest/amp0p3/out/field.dat','a003':B+'rb_c0p3_hxy_ampltest/amp0p03/out/field.dat'}
T={k:times(v) for k,v in runs.items()}
def phase_speed(path,i0,n=30,xsel=(0,lx/4)):
    ts=[];ph=[]
    for i in range(i0,i0+n):
        t,phi=read(path,i); c=phi[:,1,:].mean(axis=1)
        # vortex on positive jet: phase at x=0 region
        p1=np.sum(c*np.exp(1j*kx*xsel[1]/2)) # at x = lx/8
        ts.append(t);ph.append(np.angle(p1))
    w=-np.polyfit(ts,np.unwrap(ph),1)[0]; return w
def snap(key,tt,zsel='avg'):
    path=runs[key]; i=int(np.argmin(abs(T[key]-tt)))
    w=phase_speed(path,min(max(i-15,0),len(T[key])-31))
    c=w/ky
    o=analyse(path,i,zsel)
    U,S,Upp,A,p1=o['U'],o['S'],o['Upp'],o['A'],o['p1']
    # critical layer on positive jet: where U=c, first half box around jet at x~0
    res=[]
    for j in range(nxf):
        jn=(j+1)%nxf
        if (U[j]-c)*(U[jn]-c)<0:
            f=(c-U[j])/(U[jn]-U[j]); xc=xf[j]+f*(xf[1]-xf[0])
            Sc=S[j]+f*(S[jn]-S[j]); Ac=A[j]+f*(A[jn]-A[j])
            res.append((xc,Sc,Ac))
    jm=int(np.argmax(U))
    return o['t'],w,c,U[jm],Upp[jm],A[jm],res,o
print("hyp_x nu4=",2.4e-3*(dx/2)**4, 'dx',dx)
nu4=2.4e-3*(dx/2)**4
for key,tt in [('leg3',3000),('leg3',5900),('leg4',7500),('leg4',9000),('leg5',9800),('leg5',10300),('leg5',10600),('leg5',10900),('leg5',12000),('leg5',15000),('leg6',18500),
               ('a03',6000),('a03',6600),('a003',6000),('a003',6300),('a003',6600)]:
    t,w,c,Um,Uppm,Am,res,o=snap(key,tt)
    print(f"{key} t={t:8.1f} omega={w:+.4f} c={c:+7.1f} Umax={Um:7.1f} c/Umax={c/Um:+.3f} U''max={Uppm:+.2e} A(jetmax)={Am:9.1f}")
    for xc,Sc,Ac in res:
        wtr=ky*np.sqrt(abs(Sc)*Ac); half=2*np.sqrt(Ac/abs(Sc)); dl=(nu4/(ky*abs(Sc)))**0.2
        print(f"     crit x={xc:7.1f} S={Sc:+.3f} A={Ac:9.1f} w_tr={wtr:.4f} half={half:6.1f} rho (={half/dx:.2f} dx) delta_hyp={dl:.1f} lam=(d/w)^5={(dl/half)**5:.1e}")

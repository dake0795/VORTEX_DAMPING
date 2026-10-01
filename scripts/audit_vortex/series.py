import numpy as np
from vstruct import *
dx=lx/nx
runs={'orig':[B+f'rb_c0p3_hxy/leg_000{l}/out/field.dat' for l in [3,4,5,6]],
      'a03':[B+'rb_c0p3_hxy_ampltest/amp0p3/out/field.dat'],'a003':[B+'rb_c0p3_hxy_ampltest/amp0p03/out/field.dat']}
c=425.0
def crit(phi):
    p=phi.mean(axis=2)
    U=toreal(1j*kx*p[:,0]).real/Cxy; S=toreal(-kx**2*p[:,0]).real/Cxy; A=2*np.abs(toreal(p[:,1]))/Cxy
    out=[]
    for sgn in [1,-1]:
        cc=sgn*c
        for j in range(nxf):
            jn=(j+1)%nxf
            if (U[j]-cc)*(U[jn]-cc)<0 and abs(S[j])>1.0:
                out.append((abs(S[j]),A[j]))
    return out
res={}
for k,paths in runs.items():
    rows=[]
    for p in paths:
        ts=times(p); step=10 if k!='orig' else 40
        for i in range(0,len(ts),step):
            t,phi=read(p,i); cr=crit(phi)
            S=np.mean([a[0] for a in cr]); A=np.mean([a[1] for a in cr])
            rows.append((t,S,A,2*np.sqrt(A/S),ky*np.sqrt(S*A)))
    res[k]=np.array(rows); np.save(f'crit_{k}.npy',res[k]); print(k,len(rows))

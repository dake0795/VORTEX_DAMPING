import numpy as np
from vstruct import *
p=B+'rb_c0p3_hxy/leg_0004/out/field.dat'
xs=[32,627,1239,1834,2464,3091]
idx=[int(round(x/lx*nxf)) for x in xs]
T=[];P=[]
for i in range(0,40):
    t,phi=read(p,i)
    pp=phi[:,:,16]
    p1=toreal(pp[:,1])
    T.append(t);P.append(np.angle(p1[idx]))
T=np.array(T);P=np.unwrap(np.array(P),axis=0)
for k,x in enumerate(xs):
    w=-np.polyfit(T,P[:,k],1)[0]
    print(f"x={x:5d} omega={w:+.4f}  c=omega/ky={w/ky:+.1f}")
o=analyse(p,0)
U=o['U'];A=o['A'];ph=np.angle(o['p1'])
for j in range(0,nxf,16):
    print(f"{xf[j]:7.1f} U={U[j]:+8.1f} S={o['S'][j]:+.3f} A={A[j]:9.1f} ph={ph[j]:+.2f}")

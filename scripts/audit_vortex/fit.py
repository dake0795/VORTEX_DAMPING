import numpy as np
B='/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/'
dx=3689.25/64
o=np.load('orig.npy'); t,Ez,En=o[:,0],o[:,1],o[:,2]
cr={k:np.load(f'crit_{k}.npy') for k in ['orig','a03','a003']}
def win(t,y,a,b):
    m=(t>=a)&(t<b); p=np.polyfit(t[m],np.log(y[m]),1); return p[0]
def wmean(c,a,b,col):
    m=(c[:,0]>=a)&(c[:,0]<b); return np.exp(np.mean(np.log(c[m,col])))
rows=[]
print(" window       Gamma    A_c     w/dx   w_tr    E_nz")
for a in range(5900,11000,300):
    b=a+300; g=win(t,En,a,b); A=wmean(cr['orig'],a,b,2); w=wmean(cr['orig'],a,b,3); wt=wmean(cr['orig'],a,b,4)
    rows.append((a,b,-g,A,w/dx)); print(f"{a}-{b} {-g:.2e} {A:8.1f} {w/dx:5.2f} {wt:.3f} {np.exp(np.mean(np.log(En[(t>=a)&(t<b)]))):.2e}")
r=np.array(rows)
m=r[:,0]<10400
pA=np.polyfit(np.log(r[m,3]),np.log(r[m,2]),1); print("orig 5900-10400: Gamma ~ A^%.2f"%pA[0])
m2=r[:,0]<8300; print("orig 5900-8300: A^%.2f"%np.polyfit(np.log(r[m2,3]),np.log(r[m2,2]),1)[0])
m3=(r[:,0]>=8300)&(r[:,0]<10400); print("orig 8300-10400: A^%.2f"%np.polyfit(np.log(r[m3,3]),np.log(r[m3,2]),1)[0])
for nm,k in [('amp0p3','a03'),('amp0p03','a003')]:
    x=np.loadtxt(B+f'rb_c0p3_hxy_ampltest/{nm}/out/plunk_e_time.dat',comments='#'); tt,e=x[:,0],x[:,11]
    for a,b in [(5908,6707),(5908,6307),(6307,6707),(6008,6707),(6100,6400),(6400,6707)]:
        g=win(tt,e,a,b); A=wmean(cr[k],a,b,2); w=wmean(cr[k],a,b,3)
        pred=np.exp(np.polyval(pA,np.log(A)))
        print(f"{nm} {a}-{b} Gamma={-g:.2e}  A_c={A:7.1f} w/dx={w/dx:.2f}  orig-curve Gamma(A)={pred:.2e}")
g0=win(t,En,5908,6707); A0=wmean(cr['orig'],5908,6707,2); print("orig 5908-6707 Gamma",-g0,"A_c",A0)

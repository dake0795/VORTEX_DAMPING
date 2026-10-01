import numpy as np
B='/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/'
def load(p):
    d=np.loadtxt(p,comments='#'); return d
legs=[load(B+f'rb_c0p3_hxy/leg_000{l}/out/plunk_e_time.dat') for l in [1,3,4,5,6]]
d=np.concatenate(legs); d=d[np.argsort(d[:,0])]
t=d[:,0]; Ez=d[:,10]; En=d[:,11]; Zn=d[:,15]; sh=d[:,16]
np.save('orig.npy',np.c_[t,Ez,En,Zn,sh])
# windowed rates using linear fit of ln En over windows
def rate(t,y,a,b):
    m=(t>=a)&(t<b); 
    if m.sum()<5: return np.nan,np.nan,np.nan
    p=np.polyfit(t[m],np.log(y[m]),1); return p[0], np.exp(np.mean(np.log(y[m]))), m.sum()
print("win  Gam_nz   <En>   En/Ez   Gam_z  shear  Zn/En")
for a in list(range(800,9000,500))+list(range(9000,12000,250))+list(range(12000,19000,1000)):
    b=a+(500 if a<9000 else 250 if a<12000 else 1000)
    g,e,n=rate(t,En,a,b); gz,ez,_=rate(t,Ez,a,b)
    m=(t>=a)&(t<b)
    print(f"{a:6d}-{b:6d} {g:+.2e} {e:.3e} {e/ez:.2e} {gz:+.1e} {sh[m].mean():.3f} {np.mean(Zn[m]/En[m]):.2f}")
for nm in ['amp0p3','amp0p03']:
    x=load(B+f'rb_c0p3_hxy_ampltest/{nm}/out/plunk_e_time.dat'); tt=x[:,0]; e=x[:,11]; ez=x[:,10]
    print(nm, tt[0], tt[-1], e[0], e[0]/ez[0])
    for a in range(5900,6700,100):
        g,em,n=rate(tt,e,a,a+100); print(f"  {a}-{a+100} {g:+.2e} {em:.3e}")

import numpy as np
from field import *
lx=3689.25; kxmin=2*np.pi/lx; ky=0.002; Cxy=0.92640709665123333
kx=np.fft.fftfreq(nx,1/nx)*kxmin
nxf=1024; xf=np.arange(nxf)*lx/nxf
def toreal(c):  # c(kx) FFT order -> f(x) on fine grid
    return np.array([np.sum(c*np.exp(1j*kx*x)) for x in xf])
def prof(phi,zsel):
    if zsel=='avg': p=phi.mean(axis=2)
    else: p=phi[:,:,zsel]
    Phi=toreal(p[:,0]).real
    U=toreal(1j*kx*p[:,0]).real/Cxy
    S=toreal(-kx**2*p[:,0]).real/Cxy
    Upp=toreal(-1j*kx**3*p[:,0]).real/Cxy
    p1=toreal(p[:,1]); p1x=toreal(1j*kx*p[:,1]); p1xx=toreal(-kx**2*p[:,1])
    return Phi,U,S,Upp,p1,p1x,p1xx
def analyse(path,i,zsel=16,c=None,verbose=True):
    t,phi=read(path,i)
    Phi,U,S,Upp,p1,p1x,p1xx=prof(phi,zsel)
    out=dict(t=t)
    # zonal stats
    out['Urange']=(U.min(),U.max()); out['Srms']=np.sqrt(np.mean(S**2))
    A=2*np.abs(p1)/Cxy   # streamfunction amplitude of cos
    out['Amax']=A.max()
    # local maxima of |p1|
    pk=[j for j in range(nxf) if A[j]>=A[j-1] and A[j]>=A[(j+1)%nxf] and A[j]>0.3*A.max()]
    out['peaks']=[]
    for j in pk:
        # radial curvature of amplitude envelope: kappa^2 = -A''/A ; use |p1|
        h=xf[1]-xf[0]; a0,am,ap=A[j],A[j-1],A[(j+1)%nxf]
        kap2=-(ap-2*a0+am)/h**2/a0
        out['peaks'].append(dict(x=xf[j],A=a0,kap=np.sqrt(max(kap2,0)),U=U[j],S=S[j],Upp=Upp[j],
                                 phase=np.angle(p1[j])))
    out['U']=U; out['S']=S; out['Upp']=Upp; out['p1']=p1; out['A']=A
    # zonal jet extrema
    ext=[j for j in range(nxf) if (U[j]-U[j-1])*(U[(j+1)%nxf]-U[j])<0]
    out['jets']=[(xf[j],U[j],Upp[j]) for j in ext]
    return out
if __name__=='__main__':
    import sys
    p=B+'rb_c0p3_hxy/leg_0004/out/field.dat'
    o=analyse(p,0)
    print('t',o['t'],'U range',o['Urange'],'Srms',o['Srms'],'Amax',o['Amax'])
    for jx in o['jets']: print(' jet ext x=%.0f U=%.3f U\'\'=%.3e'%jx)
    for pk in o['peaks']: print(' vortex peak',{k:(round(v,6) if isinstance(v,float) else v) for k,v in pk.items()})
    o2=analyse(p,0,'avg'); print('zavg Amax',o2['Amax'],'Urange',o2['Urange'],'Srms',o2['Srms'])
    for pk in o2['peaks']: print(' vortex peak',pk)

import numpy as np
from scipy.integrate import quad
from scipy.special import wofz
from scipy.optimize import brentq
import mpmath as mp
sp=np.sqrt(np.pi)
Z=lambda zeta: 1j*sp*wofz(zeta)
def K0c(z):  # causal (from Im z>0), continued with principal sqrt + entire Z
    r=np.sqrt(complex(z)); return -2*(1+r*Z(r))
def K0int(w): # integral definition off the positive real axis
    w=complex(w)
    return K0c(complex(w.real,abs(w.imag))) if w.imag>=0 else np.conj(K0c(np.conj(w)))
p=lambda u: 2/sp*np.sqrt(u)*np.exp(-u)
print("--- K0 identity ---")
for z in [1+0.5j, 2+0.1j, -1+0.3j, 0.3+2j, -3+1e-9j]:
    f=lambda s: 4/sp*s**2*np.exp(-s**2)/(z-s**2)
    re=quad(lambda s: f(s).real,0,np.inf,limit=500)[0]; im=quad(lambda s: f(s).imag,0,np.inf,limit=500)[0]
    print(z, re+1j*im, K0c(z))
# real z>0 with +i0: PV - i pi p(x)
for x in [0.5,2.0,5.0]:
    pv=quad(lambda u: p(u)/(x-u)*(1) ,0,2*x,weight='cauchy',wvar=x)[0] if False else None
    g=lambda u: p(u)
    pvv=-quad(g,0,2*x,weight='cauchy',wvar=x)[0]+quad(lambda u:p(u)/(x-u),2*x,np.inf)[0]
    print('x',x, pvv-1j*np.pi*p(x), K0c(x+0j))
# K1 recurrence
z=1.3+0.4j
K1=quad(lambda u:(p(u)*u/(z-u)).real,0,np.inf)[0]+1j*quad(lambda u:(p(u)*u/(z-u)).imag,0,np.inf)[0]
print('K1',K1,z*K0c(z)-1)
print("--- Im J ---")
for a in [1.5,-1.5]:
  for om in [2.0,-2.0]:
    kU=0.3;b0,b1=0.7,-0.4;eps=1e-4
    f=lambda u: p(u)*(b0+b1*u)/(om-kU-a*u+1j*eps)
    im=quad(lambda u:f(u).imag,0,np.inf,points=[(om-kU)/a] if (om-kU)/a>0 else None,limit=1000)[0]
    ur=(om-kU)/a
    form=-2*sp/abs(a)*np.sqrt(max(ur,0))*np.exp(-ur)*(b0+b1*ur)*(ur>0)
    print(a,om,im,form)
print("--- passive ---")
a=0.7
for t in [0.5,3,20]:
    c=quad(lambda u:(p(u)*np.exp(-1j*a*u*t)).real,0,np.inf,limit=2000)[0]+1j*quad(lambda u:(p(u)*np.exp(-1j*a*u*t)).imag,0,np.inf,limit=2000)[0]
    print(t,c,(1+1j*a*t)**-1.5,abs(c)**2,(1+a*a*t*t)**-1.5)

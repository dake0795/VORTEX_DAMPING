import numpy as np, mpmath as mp
from scipy.special import wofz
from scipy.optimize import brentq
from scipy.integrate import quad
exec(open('audit.py').read().split('print("--- K0 identity ---")')[0])
# Closed model: flute, local Doppler, isotropic drift a*u (a>0 proton), Omega*=-omega*_n/a, eta
# D(z) = 1+k2 - 1/2 P(z) - 1/2 E(z), z = (omega-kU)/a
def D(z,k2,Os,eta):
    c=1-1.5*eta
    K0p=K0c(z); K1p=z*K0p-1
    w=-z; K0e=K0int(w) if w.real<0 else None
    if K0e is None: raise ValueError
    K1e=w*K0e-1
    P=(z-Os*c)*K0p-Os*eta*K1p
    E=-((z+Os*c)*K0e+Os*eta*K1e)
    return 1+k2-0.5*P-0.5*E
def DR(x,k2,Os,eta): return D(x+0j,k2,Os,eta).real
def DI_formula(x,Os,eta): return sp*np.sqrt(x)*np.exp(-x)*(x-Os*(1-1.5*eta+eta*x))
def solve(k,Os,eta):
    k2=k*k
    xs=np.linspace(0.05,60,6000); v=[DR(x,k2,Os,eta) for x in xs]
    roots=[brentq(DR,xs[i],xs[i+1],args=(k2,Os,eta)) for i in range(len(xs)-1) if v[i]*v[i+1]<0]
    out=[]
    for xr in roots:
        h=1e-6; dDR=(DR(xr+h,k2,Os,eta)-DR(xr-h,k2,Os,eta))/(2*h)
        DI=D(xr+0j,k2,Os,eta).imag
        g_wd=DI/dDR  # gamma/a (amplitude damping, positive=damped)
        # exact complex root
        f=lambda zz: D(complex(zz),k2,Os,eta)
        try:
            zr=complex(mp.findroot(lambda zz: mp.mpc(f(complex(zz))), mp.mpc(xr,-g_wd)))
        except Exception as e: zr=np.nan
        out.append((xr,dDR,DI,DI_formula(xr,Os,eta),g_wd,zr))
    return out
if __name__=="__main__":
    # sanity: Z-pinch fluid threshold
    xi=np.linspace(-1,1,20001)
    wd2=np.trapezoid(((1+xi**2)/2)**2,xi)/2*15/4; wdu=np.trapezoid((1+xi**2)/2,xi)/2*15/4
    print('Zpinch <wd^2>/wh^2',wd2,'<wd u>',wdu,'<wd w*>/(wh w*n (1+eta)) =',np.trapezoid((1+xi**2)/2,xi)/2*1.5)
    # pair symmetry D(-z)? check D_R even
    for eta in [0.0,1.0]:
      for Os in [0.0,0.5,1.0]:
        print(f'=== eta={eta} Omega*_eff={Os}')
        for k in [0.1,0.2,0.3,0.5,0.7,1.0,1.5]:
            for (xr,dDR,DI,DIf,g,zr) in solve(k,Os,eta):
                ufl=np.sqrt(max(15/4-1.5*Os*(1+eta),0))/k
                print(f' k={k:4} x_r={xr:7.3f} (fluid {ufl:6.2f}) DR\'={dDR:+.3e} DI={DI:+.3e} (formula {DIf:+.3e}) gam_wd/a={g:+.3e} exact z={zr.real:7.3f}{zr.imag:+.3e}j')

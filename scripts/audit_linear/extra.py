exec(open('model.py').read().split('if __name__')[0])
import warnings; warnings.filterwarnings('ignore')
from math import gamma as G
mom=lambda m: G(m+1.5)/G(1.5)
def closed(k,Os,eta):
    c=1-1.5*eta
    m0=mom(2)-Os*(c*mom(1)+eta*mom(2)); m1=mom(4)-Os*(c*mom(3)+eta*mom(4))
    u=np.sqrt(m0/k**2+m1/m0)
    return u, np.sqrt(np.pi)/2*u**1.5*np.exp(-u)*(u-Os*(c+eta*u))/k**2
for (Os,eta) in [(0,0),(0.5,0),(1,0),(0.5,1),(1,1),(1.2,1),(0.8,2)]:
  for k in [0.05,0.1,0.15,0.2]:
    res=[r for r in solve(k,Os,eta) if r[1]>0]
    if not res: print(Os,eta,k,'no PV root with DR>0'); continue
    xr,dDR,DI,DIf,g,zr=res[-1]
    u,gc=closed(k,Os,eta)
    print(f'Os={Os} eta={eta} k={k}: x_r={xr:.3f} closed u_r={u:.3f}  gam_wd/a={g:+.3e} closed={gc:+.3e} exact={-zr.imag:+.3e}')

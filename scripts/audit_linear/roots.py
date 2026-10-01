exec(open('model.py').read().split('if __name__')[0])
import warnings; warnings.filterwarnings('ignore')
def allroots(k,Os,eta,xr=(0.02,25),yr=(-3,3)):
    found=[]
    for x0 in np.linspace(*xr,26):
        for y0 in np.linspace(*yr,13):
            try:
                r=complex(mp.findroot(lambda zz: mp.mpc(D(complex(zz),k*k,Os,eta)), mp.mpc(x0,y0),tol=1e-20,maxsteps=60))
            except Exception: continue
            if r.real<=1e-6 or abs(D(r,k*k,Os,eta))>1e-8: continue
            if all(abs(r-f)>1e-5 for f in found): found.append(r)
    return sorted(found,key=lambda r:-r.imag)
for (Os,eta) in [(0,0),(1,0),(1,1),(0.5,0)]:
    for k in [0.1,0.3]:
        rs=allroots(k,Os,eta)
        print(Os,eta,k,[f'{r.real:.3f}{r.imag:+.3e}j' for r in rs[:6]])

"""Jets (k_y = 0) by radial harmonic m: GENE's nonlinear transfer of electrostatic energy E and of free energy W, and the
parallel-streaming term in each budget (species summed), time-averaged over a window of rb_c0p3_hxoff (hyp_x = 0).
Question: at the Debye-scale harmonics (k_x lambda_D >~ 0.5), is the stress exchanging free energy with the vortices,
or only converting the jets' E into entropy Z = W - E (which streaming then converts back)?"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import spec2d
D = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxoff/leg_0001/out"
LAMD, LX = np.sqrt(5000.0), 3689.25
def load(name, a, b):
    out = 0
    for sp in ("electrons", "positrons"):
        kx, ky, t, _, c = spec2d.read(f"{D}/Spectral_2D_{name}_{sp}.dat")
        m = (t >= a) & (t < b); out = out + c[m].mean(0)
    return kx, out
for a, b in ((2700, 3700), (3700, 4700)):
    print(f"\nwindow {a}-{b}:   m   kx lamD   NL_E        NL_W        NL_Z=W-E     par_E       par_W      dE/dt(fe)")
    kx, nE = load("E_nonlinear", a, b); _, nW = load("W_nonlinear", a, b); _, pE = load("E_parallel", a, b); _, pW = load("W_parallel", a, b)
    _, fE = load("E_fe", a, b)
    mm = np.rint(np.abs(kx) * LX / (2 * np.pi)).astype(int)
    for m in range(0, 16):
        s = mm == m
        if not s.any(): continue
        g = lambda arr: arr[s, 0].sum()
        print(f"           {m:2d}   {m*2*np.pi/LX*LAMD:5.2f}   {g(nE):+.3e}  {g(nW):+.3e}  {g(nW)-g(nE):+.3e}  {g(pE):+.3e}  {g(pW):+.3e}")
    print(f"     k_y = 1 (vortices): NL_E {nE[:,1].sum():+.3e}  NL_W {nW[:,1].sum():+.3e}   k_y >= 2: NL_E {nE[:,2:].sum():+.3e}  NL_W {nW[:,2:].sum():+.3e}")
    print(f"     totals: NL_E {nE.sum():+.3e}   NL_W {nW.sum():+.3e}")

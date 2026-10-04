"""Quasilinear check: the entropy produced by the Landau damping of all k_y >= 1 modes, -E_parallel(k_y >= 1), should be
carried by the nonlinearity into the zonal (y-averaged) distribution, NL_Z(k_y = 0) = NL_W - NL_E at k_y = 0.  Also the
free-energy neutrality of the jets, NL_W(k_y = 0), and NL conservation of E and W.  Windows of 600 time units."""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import spec2d
RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
def series(dirs, nx, q, name):
    spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    T, C = [], []
    for d in dirs:
        f = f"{d}/Spectral_2D_{q}_{name}_electrons.dat"
        if not os.path.exists(f): continue
        c = 0
        for sp in ("electrons", "positrons"):
            t, cc = spec2d.by_ky(f"{d}/Spectral_2D_{q}_{name}_{sp}.dat"); c = c + cc
        T.append(t); C.append(c)
    t = np.concatenate(T); c = np.concatenate(C); o = np.argsort(t, kind="stable"); t, c = t[o], c[o]
    keep = np.r_[True, np.diff(t) > 0]; return t[keep], c[keep]
for case, dirs, nx, t0 in (("hxoff nx64", sorted(glob.glob(RB + "/rb_c0p3_hxoff/leg_*/out")), 64, 1200),
                           ("nx128_hxoff", [RB + "/rb_c0p3_draintest/nx128_hxoff/out", RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128, 2600),
                           ("reference (hyp_x on)", sorted(glob.glob(RB + "/rb_c0p3_hxy/leg_*/out"))[:4], 64, 1200)):
    tE, nE = series(dirs, nx, "E", "nonlinear"); tW, nW = series(dirs, nx, "W", "nonlinear"); tp, pE = series(dirs, nx, "E", "parallel")
    print(f"\n{case}:   t      Landau loss k_y>=1   NL_Z(jets)/loss   NL_W(jets)/loss   NL_E sum/loss  NL_W sum/loss")
    for a in np.arange(t0, min(tE[-1], tW[-1], tp[-1]) - 599, 600):
        b = a + 600; m = lambda t: (t >= a) & (t < b)
        E = nE[m(tE)].mean(0); W = nW[m(tW)].mean(0); L = -pE[m(tp)].mean(0)[1:].sum()
        print(f"           {a+300:6.0f}   {L:.3e}           {(W[0]-E[0])/L:+.3f}            {W[0]/L:+.3f}           {E.sum()/L:+.4f}        {W.sum()/L:+.4f}")

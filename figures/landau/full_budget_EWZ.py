"""GENE's exact budgets of E, W and Z = W - E for the jets (k_y = 0), the vortices (k_y = 1) and the rest (k_y >= 2),
species summed, window means, rb_c0p3_hxoff (hyp_x = 0) and nx128_hxoff.  Terms: fe (dX/dt), nonlinear, parallel,
curvature, hyp_v, hyp_z, hyp_kperp, drive, sources, hypzcomp, coll.  Residual = fe - sum(terms)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import spec2d
RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
TERMS = ("nonlinear", "parallel", "curvature", "hyp_v", "hyp_z", "hyp_kperp", "drive", "sources", "hypzcomp")
def load(dirs, nx, q, name, a, b):
    spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    T, C = [], []
    for d in dirs:
        f = f"{d}/Spectral_2D_{q}_{name}_electrons.dat"
        if not os.path.exists(f): return None
        c = 0
        for sp in ("electrons", "positrons"):
            t, cc = spec2d.by_ky(f"{d}/Spectral_2D_{q}_{name}_{sp}.dat"); c = c + cc
        T.append(t); C.append(c)
    t = np.concatenate(T); c = np.concatenate(C); m = (t >= a) & (t < b)
    if name == "fe":
        return np.array([np.polyfit(t[m], c[m, j], 1)[0] for j in range(c.shape[1])])
    return c[m].mean(0)
for case, dirs, nx in (("hxoff nx64", [RB + "/rb_c0p3_hxoff/leg_0001/out"], 64),
                       ("nx128_hxoff", [RB + "/rb_c0p3_draintest/nx128_hxoff/out", RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    a, b = 3000, 4400
    print(f"\n{case}, t {a}-{b}")
    R = {}
    for q in ("E", "W"):
        R[q] = {n: load(dirs, nx, q, n, a, b) for n in ("fe",) + TERMS}
    P1 = -R["E"]["parallel"][1]
    print(f"   unit: Landau loss of the vortices P1 = -E_parallel(k_y=1) = {P1:.3e}")
    for q in ("E", "W", "Z"):
        print(f"   {q}:          jets(ky0)   vortices(ky1)  rest(ky>=2)")
        for n in ("fe",) + TERMS:
            if q == "Z":
                if R["E"][n] is None or R["W"][n] is None: continue
                arr = R["W"][n] - R["E"][n]
            else:
                arr = R[q][n]
                if arr is None: continue
            v = (arr[0], arr[1], arr[2:].sum())
            if max(abs(x) for x in v) / P1 < 0.005: continue
            print(f"      {n:10s}  {v[0]/P1:+.3f}      {v[1]/P1:+.3f}        {v[2]/P1:+.3f}")

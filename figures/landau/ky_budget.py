"""GENE's exact electrostatic-energy budget by k_y (species summed): k_y = 0 (jets), 1 (vortices), >= 2 (the rest);
terms fe (dE/dt from the record), nonlinear, parallel, curvature, hyp_v, hyp_z, hyp_kperp, drive, coll; rb_c0p3_hxoff
and nx128_hxoff windows.  Rates per unit of the vortices' parallel loss -E_parallel(k_y = 1)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import spec2d
RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
TERMS = ("fe", "nonlinear", "parallel", "curvature", "hyp_v", "hyp_z", "hyp_kperp", "drive")
def load(dirs, nx):
    spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    T = {}
    for term in TERMS:
        ts, cs = [], []
        for d in dirs:
            c = 0
            for sp in ("electrons", "positrons"):
                t, cc = spec2d.by_ky(f"{d}/Spectral_2D_E_{term}_{sp}.dat"); c = c + cc
            ts.append(t); cs.append(c)
        t = np.concatenate(ts); c = np.concatenate(cs); o = np.argsort(t); T[term] = (t[o], c[o])
    return T
for name, dirs, nx, wins in (("hxoff nx64", [RB + "/rb_c0p3_hxoff/leg_0001/out"], 64, ((2700, 3300), (3300, 3900), (3900, 4400))),
                             ("nx128_hxoff", [RB + "/rb_c0p3_draintest/nx128_hxoff/out", RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128, ((2700, 3300), (3300, 3900), (3900, 4400)))):
    T = load(dirs, nx)
    print(f"\n{name}  (per unit of the vortices' parallel loss Pk1 = -parallel(k_y=1))")
    for a, b in wins:
        t, fe = T["fe"]; m = (t >= a) & (t < b)
        dE = np.array([np.polyfit(t[m], fe[m, j], 1)[0] for j in range(fe.shape[1])])
        r = {term: T[term][1][(T[term][0] >= a) & (T[term][0] < b)].mean(0) for term in TERMS}
        P1 = -r["parallel"][1]
        g = lambda arr: (arr[0], arr[1], arr[2:].sum())
        print(f"  t {a}-{b}: P(k_y=1) = {P1:.3e}")
        for lab, arr in (("dE/dt", dE), ("nonlinear", r["nonlinear"]), ("parallel", r["parallel"]), ("curvature", r["curvature"]), ("hyp_v", r["hyp_v"]), ("hyp_z", r["hyp_z"]), ("hyp_kperp", r["hyp_kperp"]), ("drive", r["drive"])):
            z, o, rest = g(arr)
            print(f"     {lab:10s}  jets {z/P1:+.3f}   vortices {o/P1:+.3f}   k_y>=2 {rest/P1:+.3f}")

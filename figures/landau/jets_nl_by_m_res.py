"""Jets' nonlinear electrostatic-energy input by radial harmonic m (k_y = 0), species summed, per unit of the vortices'
Landau loss P1 = -E_parallel(k_y = 1), at the reference resolution (rb_c0p3_hxoff) and doubled radial resolution
(nx128_hxoff), t 3000-4400: where does the doubled-resolution run lose supply?"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import spec2d
RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
LX, LAMD = 3689.25, np.sqrt(5000.0)
def load(dirs, nx, name, a, b):
    spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    acc = 0; n = 0
    for d in dirs:
        for sp in ("electrons", "positrons"):
            kx, ky, t, _, c = spec2d.read(f"{d}/Spectral_2D_E_{name}_{sp}.dat"); m = (t >= a) & (t < b)
            if m.any(): acc = acc + c[m].sum(0); n += m.sum() / 2
    return kx, acc / n
for case, dirs, nx in (("nx64", [RB + "/rb_c0p3_hxoff/leg_0001/out"], 64), ("nx128", [RB + "/rb_c0p3_draintest/nx128_hxoff/out", RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    kx, nl = load(dirs, nx, "nonlinear", 3000, 4400); _, par = load(dirs, nx, "parallel", 3000, 4400)
    spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    P1 = 0
    for d in dirs:
        for sp in ("electrons", "positrons"):
            t, c = spec2d.by_ky(f"{d}/Spectral_2D_E_parallel_{sp}.dat"); m = (t >= 3000) & (t < 4400)
            if m.any(): P1 = P1 - c[m, 1].sum() / m.sum() / len(dirs) * (1 if len(dirs) == 1 else 1)
    mm = np.rint(np.abs(kx) * LX / (2 * np.pi)).astype(int)
    rowsNL = [nl[mm == m, 0].sum() / P1 for m in range(1, 16)]; rowsP = [par[mm == m, 0].sum() / P1 for m in range(1, 16)]
    tail = nl[mm >= 16, 0].sum() / P1
    print(f"{case}: P1 {P1:.3e}\n  NL_E(m)/P1, m=1..15: " + " ".join(f"{x:+.3f}" for x in rowsNL) + f"   m>=16: {tail:+.3f}   total {nl[:,0].sum()/P1:+.3f}")
    print("  par_E(m)/P1:         " + " ".join(f"{x:+.3f}" for x in rowsP) + f"   m>=16: {par[mm>=16,0].sum()/P1:+.3f}")

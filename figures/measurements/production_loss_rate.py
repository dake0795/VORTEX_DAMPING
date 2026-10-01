#!/usr/bin/env python3
"""Linear loss rate of the electrostatic energy of the k_y = 1 vortices in the eta = 1 PRODUCTION run (box 4 x larger
in x and y than the reduced box): -(E_parallel + E_curvature)/E_fe and E_nonlinear/E_fe at k_y = 1, from the last
records of Spectral_2D_E_<term>_<species>.dat (text, 24 bytes per line; see ../budget/spec2d.py for the layout)."""
import os, sys
import numpy as np
D = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/eta1_production_NEW_sep26/leg_0002/out"
NX, NKY, W = 256, 64, 24
NKX = NX - 1; HDR = NKX + NKY; FRL = 2 + NKX * NKY
NREC, STRIDE = int(sys.argv[1]) if len(sys.argv) > 1 else 240, int(sys.argv[2]) if len(sys.argv) > 2 else 10


def last(name):
    tot = None
    for sp in ("electrons", "positrons"):
        p = f"{D}/Spectral_2D_{name}_{sp}.dat"
        n = (os.path.getsize(p) // W - HDR) // FRL - 1
        ts, out = [], []
        with open(p, "rb") as f:
            for i in range(n - NREC * STRIDE, n, STRIDE):
                f.seek((HDR + i * FRL) * W)
                a = np.array(f.read(FRL * W).split(), float)
                ts.append(a[0]); d = a[2:].reshape(NKX, NKY)
                out.append(np.r_[d[:, 0].sum(), 2 * d[:, 1].sum(), 2 * d[:, 2:].sum()])
        tot = np.array(out) if tot is None else tot + np.array(out)
    return np.array(ts), tot


res = {}
for name in ("E_fe", "E_parallel", "E_curvature", "E_nonlinear", "W_hyp_v", "W_fe"):
    t, c = last(name); res[name] = c
    print(name, "t %.1f .. %.1f, %d records" % (t[0], t[-1], len(t)), "means ky=0 %.4e | ky=1 %.4e | ky>=2 %.4e" % tuple(c.mean(axis=0)), flush=True)
E1 = res["E_fe"][:, 1].mean()
print("k_y = 1:  -E_par/E = %.3e   -E_curv/E = %.3e   E_NL/E = %.3e   Z/E = %.3f   -W_hypv/(E+Z) = %.3e" % (
    -res["E_parallel"][:, 1].mean() / E1, -res["E_curvature"][:, 1].mean() / E1, res["E_nonlinear"][:, 1].mean() / E1,
    res["W_fe"][:, 1].mean() / E1, -res["W_hyp_v"][:, 1].mean() / (E1 + res["W_fe"][:, 1].mean())))
print("E_nz/E_zon = %.3e (ky=1 only: %.3e)" % ((res["E_fe"][:, 1].mean() + res["E_fe"][:, 2].mean()) / res["E_fe"][:, 0].mean(), E1 / res["E_fe"][:, 0].mean()))

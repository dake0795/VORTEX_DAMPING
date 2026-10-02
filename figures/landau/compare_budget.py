"""Jets' energy budget per odd radial harmonic: GENE's exact (E_nonlinear + E_parallel at k_y = 0) against the theory
(local stress of the two vortices + momentum deposited by their Landau damping), both as rates per unit energy of the
harmonic, averaged over a window."""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); LAND = HERE
sys.path.insert(0, "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh")
from common import *
import measure as M
sys.argv = [sys.argv[0]]
exec(open("/rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/figures/drain/zonal_budget_gene.py").read().split("for (t0, t1) in")[0])
G = np.load("/rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh/cache/gene_hxoff.npz"); tg = G["t"]; psi0 = G["phi_avg"][:, :, 0] / CXY
dep = np.load(os.path.join(LAND, "deposit_cache_hxoff.npy"), allow_pickle=True).item()
kxg = np.fft.fftfreq(NX, 1.0 / NX) * KXMIN
for (t0, t1, key) in ((1000, 2000, "1000"), (3400, 4400, "3400")):
    kx, nl, _ = rd("nonlinear", t0, t1); _, fe, _ = rd("fe", t0, t1); _, par, _ = rd("parallel", t0, t1)
    mG = np.abs(np.round(kx / KXMIN)).astype(int)
    files = [f for f in glob.glob(os.path.join(LAND, "local_stress_*.npz")) if t0 <= int(f.split("_")[-1].split(".")[0]) <= t1]
    loc = np.mean([np.load(f)["loc"] for f in files], axis=0); flu = np.mean([np.load(f)["flute"] for f in files], axis=0)
    dz = dep[key]["dzk"]
    p0 = psi0[(tg >= t0) & (tg <= t1)].mean(0)
    print(f"\nwindow {t0}-{t1} ({len(files)} stress snapshots): rates per unit energy of the harmonic")
    print("   m   GENE exact (NL+par)   theory: local+deposit   [local   deposit]   flute-Euler+deposit")
    for m in (1, 3, 5, 7, 9):
        g = (nl[mG == m].sum() + par[mG == m].sum()) / fe[mG == m].sum()
        idx = [m, NX - m]; E = sum(M.GXX * kxg[i] ** 2 * abs(p0[i]) ** 2 for i in idx)
        r = lambda f: sum(np.real(np.conj(p0[i]) * (-f[i])) for i in idx) / E
        print(f"  {m:2d}     {g:+.3e}          {r(loc)+r(dz):+.3e}        [{r(loc):+.2e} {r(dz):+.2e}]     {r(flu)+r(dz):+.3e}")

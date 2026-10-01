"""Rayleigh growth rate on the running-mean jets of the two amplitude-reduced restarts (as growth_history.py, which
does the original run): cache/growth_restarts.npz, rows = (run index 0: eps 0.3, 1: eps 0.03; t; g_plus; w_plus;
g_minus; w_minus; g_gene)."""
import numpy as np
from common import *
import linear as L, measure as M

rows = []
for k, run in enumerate(("a03", "a003")):
    t, psi = M.load(run)
    for i0 in range(0, len(t) - 51, 25):
        sl = slice(i0, i0 + 51)
        p0 = psi[sl, :, 0].mean(0)
        lam, _, _ = L.eig_gene(p0, 0.0, nmode=256)
        om = (1j * lam).real
        jp = max((j for j in range(len(lam)) if om[j] > 0), key=lambda j: lam[j].real)
        jm = max((j for j in range(len(lam)) if om[j] < 0), key=lambda j: lam[j].real)
        lg, _, _ = L.eig_gene(p0, HYP, nmode=64)
        rows.append((k, t[sl].mean(), lam[jp].real, om[jp], lam[jm].real, -om[jm], lg.real.max()))
    print(run, len(rows), flush=True)
rows = np.array(rows)
np.savez(f"{CACHE}/growth_restarts.npz", rows=rows)
np.set_printoptions(linewidth=200, precision=4, suppress=True)
print(rows)

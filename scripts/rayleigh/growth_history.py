"""Linear (Rayleigh) growth rate and frequency of the k_y = 1 modes on the measured jet profile through the run:
running-mean zonal profile over ~60 time units (51 frames); inviscid operator at 256 radial modes (converged), both
unstable modes (one per direction); and the operator at GENE's 64 modes with GENE's hyperviscosity.
cache/growth_history.npz, columns in `cols`."""
import numpy as np
from common import *
import linear as L, measure as M

rows = []
for run in ("leg1", "leg3", "leg4", "leg5"):
    t, psi = M.load(run)
    for i0 in range(0, len(t) - 51, 25):
        sl = slice(i0, i0 + 51)
        p0 = psi[sl, :, 0].mean(0)
        if np.abs(p0).max() == 0:
            continue
        lam, _, _ = L.eig_gene(p0, 0.0, nmode=256)
        om = (1j * lam).real
        jp = max((j for j in range(len(lam)) if om[j] > 0), key=lambda j: lam[j].real)
        jm = max((j for j in range(len(lam)) if om[j] < 0), key=lambda j: lam[j].real)
        lg, _, _ = L.eig_gene(p0, HYP, nmode=64)
        jg = np.argmax(lg.real)
        Ez, En = M.energies(psi[sl])
        w, ap, am, r = M.two_freq(t[sl], psi[sl, :, 1], w0=0.95, span=0.35, nw=141)
        U = M.fine(p0, 512, 1).real
        rows.append((t[sl].mean(), lam[jp].real, om[jp], lam[jm].real, -om[jm], lg[jg].real, abs((1j * lg[jg]).real), Ez.mean(), En.mean(), w, r, U.max(), -U.min()))
    print(run, len(rows), flush=True)
rows = np.array(rows)
np.savez(f"{CACHE}/growth_history.npz", rows=rows, cols="t g_plus w_plus g_minus w_minus g_gene w_gene Ezon Enz w_meas resid Umax Umin")
np.set_printoptions(linewidth=200, precision=4, suppress=True)
print(rows[::10][:, [0, 1, 2, 3, 4, 5, 9, 10]])

"""Rayleigh growth rate gamma_R of the k_y = 1 modes on the running-mean jet of rb_c0p3_hxoff (the reference run from
noise with hyp_x = 0) and, for comparison, of the reference legs 1 and 3, every 25 frames over 51-frame windows, with the
fluid energies of the window.  cache/growth_hxoff.npz, rows (t, g_plus, w_plus, g_minus, w_minus, Ezon, Enz)."""
import numpy as np
from common import *
import linear as L, measure as M
w = metric()["w"]
p = f"{RB}/rb_c0p3_hxoff/leg_0001/out/field.dat"
n = nframes(p); T, A = [], []
for i in range(n):
    t, phi = read(p, i); T.append(t); A.append(phi @ w)
t = np.array(T); psi = np.array(A) / CXY
np.savez(f"{CACHE}/gene_hxoff.npz", t=t, phi_avg=np.array(A))
rows = []
for i0 in range(0, len(t) - 51, 25):
    sl = slice(i0, i0 + 51)
    if t[sl].mean() < 900: continue
    p0 = psi[sl, :, 0].mean(0)
    lam, _, _ = L.eig_gene(p0, 0.0, nmode=256)
    om = (1j * lam).real
    jp = max((j for j in range(len(lam)) if om[j] > 0), key=lambda j: lam[j].real)
    jm = max((j for j in range(len(lam)) if om[j] < 0), key=lambda j: lam[j].real)
    Ez, En = M.energies(psi[sl])
    rows.append((t[sl].mean(), lam[jp].real, om[jp], lam[jm].real, -om[jm], Ez.mean(), En.mean()))
rows = np.array(rows)
np.savez(f"{CACHE}/growth_hxoff.npz", rows=rows)
print(rows.shape); np.set_printoptions(linewidth=200, precision=4, suppress=True); print(rows[::8])

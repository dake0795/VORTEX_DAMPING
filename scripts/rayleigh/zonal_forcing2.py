"""Breathing-averaged zonal forcing of the jets by the two vortices, from the two-frequency fit of the k_y = 1 flute
potential in windows (a_+ e^{-i w t} + a_- e^{+i w t}); the zonal Reynolds-stress forcing of each travelling mode is
steady, the cross term oscillates at 2w and averages out.  Compared with the observed slow change of the zonal
vorticity between window centres.  Usage: zonal_forcing2.py RUNCACHE TSTART TEND [WIN STEP]"""
import sys, numpy as np
from common import *
import measure as M
src = open("zonal_forcing.py").read()
d = np.load(sys.argv[1]); t = d["t"]; psi = d["phi_avg"] / CXY
GXX, GYY = M.GXX, M.GYY
nxg, nyg = 3 * NX // 2, 3 * NKY
m = np.fft.fftfreq(nxg, 1.0 / nxg); n = np.arange(nyg // 2 + 1)
kx, ky = m[:, None] * KXMIN, n[None, :] * KYMIN
K = GXX * kx ** 2 + GYY * ky ** 2
exec(src[src.index("def pad"):src.index("kxg = np.fft")])
kxg = np.fft.fftfreq(NX, 1.0 / NX) * KXMIN
T0, T1 = float(sys.argv[2]), float(sys.argv[3]); WIN = float(sys.argv[4]) if len(sys.argv) > 4 else 60.0; STEP = float(sys.argv[5]) if len(sys.argv) > 5 else 100.0
rows, F, Z = [], [], []
for c in np.arange(T0, T1 + 1e-6, STEP):
    w = (t >= c - WIN / 2) & (t < c + WIN / 2)
    if w.sum() < 20: continue
    ww, ap, am, r = M.two_freq(t[w], psi[w, :, 1], w0=1.0, span=0.4, nw=161)
    p0 = psi[w, :, 0].mean(0)
    NH = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    bp = np.zeros((NX, NKY), complex); bm = np.zeros((NX, NKY), complex)
    tt = t[w] - t[w][0]
    for nn in range(1, NH + 1):
        B = np.stack([np.exp(-1j * nn * ww * tt), np.exp(1j * nn * ww * tt)], 1)
        sol, *_ = np.linalg.lstsq(B, psi[w, :, nn], rcond=None)
        bp[:, nn], bm[:, nn] = sol[0], sol[1]
    f = zonal_euler_tendency(bp) + zonal_euler_tendency(bm)
    HXZ = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0       # radial hyperdiffusion of the run: its zonal damping is added to the "Euler" forcing
    if HXZ > 0:
        nuz = HXZ * (DX_G * kxg / 2) ** 4 * (1 + 1 / (5000.0 * np.where(kxg == 0, 1, kxg ** 2)))
        f = f - nuz * (-(GXX * kxg ** 2) * p0) * (kxg != 0)
    F.append(f); Z.append(-(GXX * kxg ** 2) * p0)
    rate = float(np.sum((np.conj(p0) * (-f)).real))
    rows.append((t[w].mean(), ww, r, rate, 0.5 * np.sum(GXX * kxg ** 2 * np.abs(p0) ** 2), M.energies(psi[w])[1].mean()))
rows = np.array(rows); F = np.array(F); Z = np.array(Z)
tc = rows[:, 0]
dz = np.gradient(Z, tc, axis=0)
for i in range(1, len(tc) - 1, max(1, len(tc) // 12)):
    obs = dz[i]; eul = F[i]
    p0 = -Z[i] / np.where(kxg == 0, 1, GXX * kxg ** 2) * (kxg != 0)
    eo = float(np.sum((np.conj(p0) * (-obs)).real))
    c = np.vdot(eul, obs).real / np.vdot(eul, eul).real
    print(f"t {tc[i]:7.0f}: w {rows[i,1]:.3f} fit resid {rows[i,2]:.2f}; dE_zon/dt observed {eo:+.3e}, Euler {rows[i,3]:+.3e}; |obs| {np.linalg.norm(obs):.2e} |eul| {np.linalg.norm(eul):.2e} proj {c:+.2f} resid {np.linalg.norm(obs-c*eul)/np.linalg.norm(obs):.2f}")
# long-window comparison (less noise): mean over all
obs = (Z[-3:].mean(0) - Z[:3].mean(0)) / (tc[-3:].mean() - tc[:3].mean()); eul = F.mean(0)
p0 = -Z.mean(0) / np.where(kxg == 0, 1, GXX * kxg ** 2) * (kxg != 0)
print(f"WHOLE {tc[0]:.0f}-{tc[-1]:.0f}: dE_zon/dt observed {np.sum((np.conj(p0)*(-obs)).real):+.3e}, Euler {np.sum((np.conj(p0)*(-eul)).real):+.3e}; "
      f"|obs| {np.linalg.norm(obs):.2e} |eul| {np.linalg.norm(eul):.2e} proj {np.vdot(eul,obs).real/np.vdot(eul,eul).real:+.2f} "
      f"resid {np.linalg.norm(obs-np.vdot(eul,obs).real/np.vdot(eul,eul).real*eul)/np.linalg.norm(obs):.2f}")
np.savez(f"{CACHE}/zforce2_{sys.argv[1].split('gene_')[1].split('.')[0]}{'_h' if len(sys.argv) > 7 else ''}.npz", rows=rows, F=F, Z=Z)

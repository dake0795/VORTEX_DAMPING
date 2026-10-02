"""Zonal-flow budget of the flute potential in a GENE run: observed d zeta_zon/dt over long windows against the 2D-Euler
Reynolds-stress forcing computed from the same flute potential (Jacobian-weighted field-line average, averaged metric).
The residual is the zonal forcing that 2D Euler does not contain (kinetic: momentum deposited by the Landau damping,
or a sink).  Usage: zonal_forcing.py RUNCACHE t0 t1 [t0 t1 ...]  (RUNCACHE = cache/gene_<run>.npz)"""
import sys, numpy as np
from common import *
import measure as M
d = np.load(sys.argv[1]); t = d["t"]; psi = d["phi_avg"] / CXY
GXX, GYY = M.GXX, M.GYY
nxg, nyg = 3 * NX // 2, 3 * NKY
m = np.fft.fftfreq(nxg, 1.0 / nxg); n = np.arange(nyg // 2 + 1)
kx, ky = m[:, None] * KXMIN, n[None, :] * KYMIN
K = GXX * kx ** 2 + GYY * ky ** 2
def pad(p):
    a = np.zeros((nxg, nyg // 2 + 1), complex); h = NX // 2
    a[:h, :NKY] = p[:h]; a[nxg - (h - 1):, :NKY] = p[NX - (h - 1):]
    return a * nxg * nyg
def zonal_euler_tendency(p):
    """d zeta/dt at k_y = 0 from -[psi_x zeta_y - psi_y zeta_x], zeta = -K psi; returned as GENE kx-ordered array."""
    P = pad(p); Z = -K * P
    px = np.fft.irfft2(1j * kx * P, s=(nxg, nyg)); py = np.fft.irfft2(1j * ky * P, s=(nxg, nyg))
    zx = np.fft.irfft2(1j * kx * Z, s=(nxg, nyg)); zy = np.fft.irfft2(1j * ky * Z, s=(nxg, nyg))
    N = -np.fft.rfft2(px * zy - py * zx)[:, 0] / (nxg * nyg)
    h = NX // 2; out = np.zeros(NX, complex); out[:h] = N[:h]; out[NX - (h - 1):] = N[nxg - (h - 1):]
    return out
kxg = np.fft.fftfreq(NX, 1.0 / NX) * KXMIN
args = [float(x) for x in sys.argv[2:]]
for t0, t1 in zip(args[::2], args[1::2]):
    a = (t > t0 - 25) & (t < t0 + 25); b = (t > t1 - 25) & (t < t1 + 25); w = (t >= t0) & (t <= t1)
    z0 = -(GXX * kxg ** 2) * psi[a, :, 0].mean(0); z1 = -(GXX * kxg ** 2) * psi[b, :, 0].mean(0)
    obs = (z1 - z0) / (t[b].mean() - t[a].mean())
    idx = np.where(w)[0][::2]
    eul = np.mean([zonal_euler_tendency(psi[i]) for i in idx], axis=0)
    res = obs - eul
    # energy rates of the jets: E_zon = 1/2 sum K |psi|^2 -> dE/dt = sum Re(conj(psi0) * (-dzeta/dt))  (zeta = -K psi)
    p0 = psi[w, :, 0].mean(0)
    rate = lambda f: float(np.sum((np.conj(p0) * (-f)).real))
    Ez = 0.5 * np.sum(GXX * kxg ** 2 * np.abs(p0) ** 2)
    U = M.fine(p0, 512, 1).real
    to_x = lambda f: M.fine(f / np.where(kxg == 0, 1, -GXX * kxg ** 2) * (kxg != 0), 512, 1).real   # velocity tendency dU/dt
    print(f"window {t0:.0f}-{t1:.0f} ({len(idx)} frames): E_zon {Ez:.4e};  dE_zon/dt observed {rate(obs):+.3e}, Euler {rate(eul):+.3e}, residual {rate(res):+.3e};"
          f"  |res|/|obs| (L2) {np.linalg.norm(res)/np.linalg.norm(obs):.2f}, |eul|/|obs| {np.linalg.norm(eul)/np.linalg.norm(obs):.2f}")
    np.savez(f"{CACHE}/zforce_{sys.argv[1].split('gene_')[1].split('.')[0]}_{int(t0)}_{int(t1)}.npz", obs=obs, eul=eul, res=res, p0=p0, U=U,
             dUobs=to_x(obs), dUeul=to_x(eul), dUres=to_x(res))

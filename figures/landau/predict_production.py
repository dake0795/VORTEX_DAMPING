#!/usr/bin/env python3
"""The prediction of eq:Hrate for the eta = 1 PRODUCTION run (box four times larger in x and y: nx0 256, nky0 64,
kymin 0.0005, lx 14757), from the last 60 time units of leg_0002 (every 10th field frame), with the dissipation
function D(wt) of cache.npz (same flux-tube geometry).  Compare with the measured parallel-channel loss rate of
../measurements/production_loss_rate.py."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf
import bounce
P = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/eta1_production_NEW_sep26/leg_0002/out"
NX, NKY, NZ, LX, KY, CXY = 256, 64, 32, 14757.0, 0.0005, 0.92640709665123333
FR = 16 + 4 + NX * NKY * NZ * 16 + 4
KX = np.fft.fftfreq(NX, 1.0 / NX) * 2 * np.pi / LX
NXF = 1024


def to_x(c, deriv=0):                          # c[..., nx] -> fine grid
    c = np.asarray(c) * (1j * KX) ** deriv
    pad = np.zeros(c.shape[:-1] + (NXF,), complex); h = NX // 2
    pad[..., :h] = c[..., :h]; pad[..., NXF - h + 1:] = c[..., h + 1:]
    return np.fft.ifft(pad, axis=-1) * NXF


geo = bounce.geometry(P + "/dipole_fix.dat")
wz = geo["J"] / geo["J"].sum(); Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
C = np.load(os.path.join(HERE, "cache.npz")); wg, Dth = C["wg"], C["Dth"]
n = os.path.getsize(P + "/field.dat") // FR
T, P0, P1 = [], [], []
with open(P + "/field.dat", "rb") as f:
    for i in range(n - 750, n, 10):
        f.seek(i * FR + 4); t = np.frombuffer(f.read(8), "<f8")[0]
        f.seek(i * FR + 20); phi = np.frombuffer(f.read(NX * NKY * NZ * 16), "<c16").reshape((NX, NKY, NZ), order="F")
        T.append(t); P0.append(phi[:, 0, :].copy()); P1.append(phi[:, 1, :].copy())
t = np.array(T); p0 = np.array(P0); p1 = np.array(P1)
print("frames", len(t), "t", t[0], "..", t[-1], flush=True)
# dominant frequencies of the mid-plane k_y = 1 potential
mid = p1[:, :, NZ // 2]
k0 = np.argmax(np.abs(mid).mean(axis=0))
dt = np.median(np.diff(t)); sp = np.fft.fft(mid[:, k0] * np.hanning(len(t))); fr = -2 * np.pi * np.fft.fftfreq(len(t), dt)
wpos = fr[np.argmax(np.where(fr > 0, np.abs(sp), 0))]; wneg = fr[np.argmax(np.where(fr < 0, np.abs(sp), 0))]
print("dominant k_x index", k0, "spectral peaks at", wpos, wneg, flush=True)
w, a_mid, res = vf.two_freq_fit(t, mid, (wpos, wneg))
M = np.exp(-1j * np.outer(t - t.mean(), w))
a, *_ = np.linalg.lstsq(M, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, NX, NZ)
resid = np.sum(np.abs(p1.reshape(len(t), -1) - M @ a.reshape(2, -1)) ** 2) / np.sum(np.abs(p1) ** 2)
U = to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), deriv=1).real / CXY
K2L2 = 5000.0 * KY ** 2; dx = LX / NXF
num = num_me = den = 0.0
s = Gyy * (geo["gxx"] / Gxx - geo["gyy"] / Gyy)
for k in range(2):
    ph = to_x(a[k].T).T; dph = to_x(a[k].T, deriv=1).T
    ph0 = (ph * wz[None, :]).sum(axis=1); dph0 = (dph * wz[None, :]).sum(axis=1); ph1 = ph - ph0[:, None]
    wt = w[k] - KY * U
    num += 2 * K2L2 * KY ** 2 * (np.abs(ph0) ** 2 * np.interp(np.abs(wt), wg, Dth)).sum() * dx
    Dloc = -wt * (np.imag(ph1 / ph0[:, None]) * s[None, :] * wz[None, :]).sum(axis=1) / K2L2
    num_me += 2 * K2L2 * KY ** 2 * (np.abs(ph0) ** 2 * Dloc).sum() * dx
    den += (Gxx * np.abs(dph0) ** 2 + Gyy * KY ** 2 * np.abs(ph0) ** 2).sum() * dx
    print(f"vortex {k}: w {w[k]:+.4f}; |wt| range {np.abs(wt).min():.2f}..{np.abs(wt).max():.2f}; U range {U.min():.0f}..{U.max():.0f}; |wt| weighted by |phi0|^2: {(np.abs(wt)*np.abs(ph0)**2).sum()/(np.abs(ph0)**2).sum():.2f}")
print(f"two-frequency residual {resid:.3f}")
np.savez(os.path.join(HERE, "cache_production.npz"), predicted=num / den, via_measured_phi1=num_me / den, measured_parallel=1.70e-4, t0=t[0], t1=t[-1], w=w)
print(f"production: predicted L/E = {num / den:.3e};  from the measured non-flute part {num_me / den:.3e};  measured parallel-channel rate 1.70e-04 (t 741-936)")

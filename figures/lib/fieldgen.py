"""Resolution-independent reader and vortex analysis for the reduced-box runs (4 Oct 2026).
Run(dirs, nx): field.dat frames of any radial resolution (nky0 16, nz0 32, same lx), flute (Jacobian-weighted)
average along the line, two-frequency fit of the k_y = 1 vortices, non-flute projection on g^yy - <g^yy>."""
import os
import numpy as np
from scipy.optimize import minimize

RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
NKY, NZ, LX, KYMIN, CXY = 16, 32, 3689.25, 0.002, 0.92640709665123333
LAM2 = 5000.0
ZMID = NZ // 2


def geometry(path=f"{RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"):
    L = open(path).read().split("\n"); i0 = [i for i, l in enumerate(L) if l.strip() == "/"][0]
    d = np.array([[float(x) for x in l.split()] for l in L[i0 + 1:] if l.strip()])
    w = d[:, 10] / d[:, 10].sum()
    return dict(gxx=d[:, 0], gyy=d[:, 3], B=d[:, 6], J=d[:, 10], w=w, Gxx=(w * d[:, 0]).sum(), Gyy=(w * d[:, 3]).sum())


class Run:
    def __init__(self, dirs, nx):
        self.dirs = [d for d in dirs if os.path.exists(os.path.join(d, "field.dat"))]
        self.nx = nx; self.N = nx * NKY * NZ; self.FR = 16 + 4 + self.N * 16 + 4
        self.kx = np.fft.fftfreq(nx, 1.0 / nx) * 2 * np.pi / LX
        self.index = []                                   # (t, path, i), later legs win on overlaps
        last = -np.inf
        for d in self.dirs:
            p = os.path.join(d, "field.dat"); n = os.path.getsize(p) // self.FR
            with open(p, "rb") as f:
                for i in range(n):
                    f.seek(i * self.FR + 4); t = np.frombuffer(f.read(8), "<f8")[0]
                    if t > last + 1e-9:
                        self.index.append((t, p, i)); last = t
        self.t = np.array([r[0] for r in self.index])

    def frame(self, j):
        t, p, i = self.index[j]
        with open(p, "rb") as f:
            f.seek(i * self.FR + 20); b = f.read(self.N * 16)
        return t, np.frombuffer(b, "<c16").reshape((self.nx, NKY, NZ), order="F")

    def window(self, tc, half=30.0, kys=(0, 1)):
        J = np.where((self.t >= tc - half) & (self.t <= tc + half))[0]
        T, P = [], []
        for j in J:
            t, phi = self.frame(j); T.append(t); P.append(phi[:, list(kys), :].copy())
        return np.array(T), np.array(P)                 # P[nt, nx, len(kys), nz]

    def to_x(self, c, nxf=512, deriv=0):
        """coefficients c[nx, ...] (FFT order along axis 0) -> real-space profile on nxf points."""
        c = np.asarray(c); kk = self.kx.reshape((-1,) + (1,) * (c.ndim - 1))
        c = c * (1j * kk) ** deriv
        pad = np.zeros((nxf,) + c.shape[1:], complex); h = self.nx // 2
        pad[:h] = c[:h]; pad[nxf - h + 1:] = c[h + 1:]
        return np.fft.ifft(pad, axis=0) * nxf


def two_freq(t, p, w0):
    tc = t - t.mean(); norm = np.sum(np.abs(p) ** 2)
    def solve(w):
        M = np.exp(-1j * np.outer(tc, w)); a, *_ = np.linalg.lstsq(M, p.reshape(len(t), -1), rcond=None)
        return a, np.sum(np.abs(p.reshape(len(t), -1) - M @ a) ** 2) / norm
    r = minimize(lambda w: solve(w)[1], w0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-12))
    a, res = solve(r.x)
    w = r.x
    if w[0] < w[1]:
        w, a = w[::-1], a[::-1]
    return w, a.reshape((2,) + p.shape[1:]), res


def freq_guess(t, s):
    j = np.argmax((np.abs(s) ** 2).sum(0)); x = s[:, j] - s[:, j].mean()
    ts = np.linspace(t[0], t[-1], 2048); xs = np.interp(ts, t, x.real) + 1j * np.interp(ts, t, x.imag)
    fr = np.fft.fftfreq(len(ts), ts[1] - ts[0]) * 2 * np.pi; Sp = np.abs(np.fft.fft(xs)) ** 2
    return -fr[np.argmax(np.where(fr < -0.05, Sp, 0))], -fr[np.argmax(np.where(fr > 0.05, Sp, 0))]

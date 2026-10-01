r"""Bounce-resonant (Landau) damping of the non-flute part of a potential in the flux tube of the simulations.

Unperturbed orbits along the field line (GENE normalisation: z in [-pi, pi), v_par in units of v_T = sqrt(2T/m),
mu B in units of T, time in L_ref/c_ref):
    dz/dt = C_xy v_T v_par / (J B),     v_par = +- sqrt(eps) sqrt(1 - lambda B(z)),     eps = v_par^2 + mu B,  lambda = mu/eps,
so the orbit frequency is omega_b = sqrt(eps) omhat(lambda).  Passing orbits (lambda < 1/B_max) go round the periodic
flux tube; trapped ones (1/B_max < lambda < 1/B_min) bounce about the minimum of B at z = 0.

For a potential phi(z) = sum_m c_m e^{i m z} the n-th orbit harmonic is phi_n(lambda) = sum_m K[lambda, n, m] c_m with
    K = (1/2pi) oint e^{i m z(theta)} e^{-i n theta} d theta        (theta = omega_b t, the orbit angle).
The rate at which particles take energy from a potential oscillating at the (Doppler-shifted) frequency wt is
    V <<wt^2 sum_n |phi_n|^2 delta(wt - n omega_b)>> = 2 pi^2 C_xy v_T pi^{-3/2} G,
    G = int d lambda sum_n exp(-eps_n) |wt|^4 / (|n|^3 omhat^4) |phi_n|^2,     eps_n = (wt / (n omhat))^2,
(V = int J dz; << >> the Maxwellian phase-space average), the sum being over n != 0 of the co-moving harmonics for
passing orbits (the two directions of motion give n and -n) and over n >= 1 for trapped ones.
"""
import numpy as np

CXY, VT = 0.92640709665123333, np.sqrt(2.0)


def geometry(path):
    lines = open(path).read().split("\n/\n", 1)[1]
    g = np.array([[float(v) for v in l.split()] for l in lines.strip().splitlines()])
    return dict(gxx=g[:, 0], gyy=g[:, 3], B=g[:, 6], J=g[:, 10], nz=len(g))


def _interp(f, z):
    """Fourier interpolation of a periodic function sampled at z_k = -pi + 2 pi k/nz."""
    nz = len(f)
    c = np.fft.fft(f) / nz
    m = np.fft.fftfreq(nz, 1.0 / nz)
    c = c * np.exp(1j * m * np.pi)                 # samples start at z = -pi
    if nz % 2 == 0:
        keep = np.abs(m) < nz // 2
        out = (c[keep][None, :] * np.exp(1j * np.outer(z, m[keep]))).sum(axis=1)
        out = out + (np.fft.fft(f)[nz // 2] / nz) * np.cos((nz // 2) * (z + np.pi))
        return out.real
    return (c[None, :] * np.exp(1j * np.outer(z, m))).sum(axis=1).real


class Orbits:
    def __init__(self, geo, nlam_p=48, nlam_t=64, nmax=8, npsi=1200):
        self.geo, self.nmax = geo, nmax
        nz = geo["nz"]
        self.m = np.fft.fftfreq(nz, 1.0 / nz)                     # z harmonics of the potential (FFT order)
        zf = -np.pi + 2 * np.pi * (np.arange(4096) + 0.5) / 4096
        Bf = _interp(geo["B"], zf)
        self.Bmax, self.Bmin = Bf.max(), Bf.min()
        zz = np.linspace(-0.2, 0.2, 40001); self.zmin = zz[np.argmin(_interp(geo["B"], zz))]; self.Bmin = _interp(geo["B"], np.array([self.zmin]))[0]
        self.JBf = lambda z: _interp(geo["J"] * geo["B"], z)
        self.Bfun = lambda z: _interp(geo["B"], z)
        lam, wl, kind, omhat, K = [], [], [], [], []
        ns_p = np.r_[-np.arange(nmax, 0, -1), np.arange(1, nmax + 1)]
        ns_t = np.arange(1, nmax + 1)
        # --- passing: Gauss-Legendre nodes in s, lambda = (1 - s^2)/B_max  (clusters nodes at the trapped boundary)
        s, ws = np.polynomial.legendre.leggauss(nlam_p); s = 0.5 * (s + 1); ws = 0.5 * ws
        for si, wi in zip(s, ws):
            l = (1 - si ** 2) / self.Bmax
            z = -np.pi + 2 * np.pi * (np.arange(npsi * 4) + 0.5) / (npsi * 4); dz = 2 * np.pi / (npsi * 4)
            dt = self.JBf(z) / np.sqrt(np.maximum(1 - l * self.Bfun(z), 1e-14)) * dz / (CXY * VT)
            tau = dt.sum(); th = 2 * np.pi * (np.cumsum(dt) - 0.5 * dt) / tau
            k = np.array([[(np.exp(1j * mm * z - 1j * n * th) * dt).sum() / tau for mm in self.m] for n in ns_p])
            lam.append(l); wl.append(wi * 2 * si / self.Bmax); kind.append(0); omhat.append(2 * np.pi / tau); K.append(k)
        # --- trapped: lambda = 1/B_max + (1/B_min - 1/B_max) (1 - s^2)... nodes cluster at both ends
        u, wu = np.polynomial.legendre.leggauss(nlam_t); u = 0.5 * (u + 1); wu = 0.5 * wu
        l0, l1 = 1 / self.Bmax, 1 / self.Bmin
        for ui, wi in zip(u, wu):
            q = np.sin(0.5 * np.pi * ui) ** 2                      # 0..1, clustered at both ends
            l = l0 + (l1 - l0) * q
            dq = np.pi * np.sin(0.5 * np.pi * ui) * np.cos(0.5 * np.pi * ui)
            zl, zr = self._bounce(l)
            zc, zb = 0.5 * (zl + zr), 0.5 * (zr - zl)
            psi = -0.5 * np.pi + np.pi * (np.arange(npsi) + 0.5) / npsi; dpsi = np.pi / npsi
            z = zc + zb * np.sin(psi)
            dt = self.JBf(z) * zb * np.cos(psi) / np.sqrt(np.maximum(1 - l * self.Bfun(z), 1e-14)) * dpsi / (CXY * VT)
            half = dt.sum(); tau = 2 * half
            th = np.pi * (np.cumsum(dt) - 0.5 * dt) / half        # 0..pi on the forward half
            k = np.zeros((2 * nmax, nz), complex)
            for j, n in enumerate(ns_t):
                row = np.array([(np.exp(1j * mm * z) * np.cos(n * th) * dt).sum() / half for mm in self.m])
                k[nmax + j] = row                                  # stored in the n > 0 slots; n < 0 slots stay empty
            lam.append(l); wl.append(wi * (l1 - l0) * dq); kind.append(1); omhat.append(2 * np.pi / tau); K.append(k)
        self.lam, self.wl, self.kind = np.array(lam), np.array(wl), np.array(kind)
        self.omhat, self.K, self.ns = np.array(omhat), np.array(K), ns_p

    def _bounce(self, l):
        """The two bounce points B(z) = 1/lambda on either side of the minimum of B (bisection)."""
        out = []
        for sgn in (-1.0, 1.0):
            a, b = self.zmin, self.zmin + sgn * np.pi
            for _ in range(60):
                c = 0.5 * (a + b)
                if l * self.Bfun(np.array([c]))[0] < 1:
                    a = c
                else:
                    b = c
            out.append(0.5 * (a + b))
        return out

    def G(self, cm, wt):
        """cm[nx, nz]: z-FFT coefficients (FFT order, c_m = (1/nz) sum_k phi(z_k) e^{-i m z_k}); wt[nx]: Doppler-shifted
        frequency.  Returns G[nx] (see module docstring) and its passing / trapped parts."""
        phin = np.einsum("lnm,xm->xln", self.K, cm)               # [nx, nlam, 2 nmax]
        n = self.ns[None, None, :]
        om = self.omhat[None, :, None]
        w = np.abs(wt)[:, None, None]
        eps = (w / (np.abs(n) * om)) ** 2
        W = np.exp(-eps) * w ** 4 / (np.abs(n) ** 3 * om ** 4)
        terms = (W * np.abs(phin) ** 2).sum(axis=2) * self.wl[None, :]
        return terms.sum(axis=1), terms[:, self.kind == 0].sum(axis=1), terms[:, self.kind == 1].sum(axis=1)


def zfft(phi_xz):
    """phi[..., nz] sampled at z_k = -pi + 2 pi k/nz -> coefficients c_m of e^{i m z} (FFT order)."""
    nz = phi_xz.shape[-1]
    m = np.fft.fftfreq(nz, 1.0 / nz)
    return np.fft.fft(phi_xz, axis=-1) / nz * np.exp(1j * m * np.pi)


if __name__ == "__main__":
    RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
    geo = geometry(f"{RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat")
    print("B max/min", geo["B"].max(), geo["B"].min(), " B(z=-pi)", geo["B"][0], "B(z=0)", geo["B"][16], " J(-pi), J(0)", geo["J"][0], geo["J"][16], " V = int J dz =", geo["J"].sum() * 2 * np.pi / 32)
    O = Orbits(geo)
    print("omhat passing: lambda=0 ->", O.omhat[O.kind == 0][np.argmin(O.lam[O.kind == 0])], " expected 2 pi C vT / int J B dz =", 2 * np.pi * CXY * VT / ((geo["J"] * geo["B"]).sum() * 2 * np.pi / 32))
    # deeply trapped limit: omhat -> (C vT/(J0 B0)) sqrt(B''/(2 B0))
    z = np.array([-0.02, 0.0, 0.02]); B3 = O.Bfun(z); Bpp = (B3[0] - 2 * B3[1] + B3[2]) / 0.02 ** 2
    J0B0 = O.JBf(np.array([0.0]))[0]
    print("omhat deeply trapped ->", O.omhat[O.kind == 1][np.argmax(O.lam[O.kind == 1])], " expected", CXY * VT / J0B0 * np.sqrt(Bpp / (2 * B3[1])))
    print("sum of lambda weights (should be 1/B_min = %.4f):" % (1 / O.Bmin), O.wl.sum())
    # flute potential: no harmonics
    c = np.zeros((1, 32), complex); c[0, 0] = 1.0
    print("G for a flute potential:", O.G(c, np.array([1.0]))[0])
    # a single harmonic e^{iz}: passing harmonic n = 1 dominated
    c = np.zeros((1, 32), complex); c[0, 1] = 1.0
    for w in (0.3, 1.0, 2.0, 4.0):
        print("phi = e^{iz}, wt = %.1f: G, passing, trapped =" % w, [float(v[0]) for v in O.G(c, np.array([w]))])

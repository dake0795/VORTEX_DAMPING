r"""Self-consistent non-flute part of a nearly flute-like mode in the flux tube, and its Landau damping.

Field equation along the field line at a given radius, to leading order in k^2 lambda_D^2 (see appendices/appH_landau.tex):
    phi_1(l) - rho[phi_1](l) = S(l) phi_0,      S(l) = k^2 lambda_D^2 <g^yy> [ g^xx(l)/<g^xx> - g^yy(l)/<g^yy> ],
where rho[phi] is the density response of particles that stream along the field line and see the potential at the
Doppler-shifted frequency wt:
    rho[phi](l) = (1/n_0) int d^3v f_0 sum_n  wt/(wt - n omega_b + i0) phi_n e^{i n theta}.
In the basis e^{i m z} with the volume weight J:  (Jm - R(wt)) c = s phi_0,  with
    R_{m m'} = (pi/2) C v_T pi^{-3/2} int d lambda (2 pi/omhat) sum_n I_n(zeta) conj(K_nm) K_nm',   zeta = |wt|/(|n| omhat),
    I_n = 2 zeta^2 sqrt(pi) (2 zeta D(zeta) - 1) - 2 pi i sign(wt) zeta^3 e^{-zeta^2}  (n >= 1; D = Dawson),   I_0 = sqrt(pi),
(passing: both directions; trapped: half of I_0 and the cosine kernel).
"""
import numpy as np
from scipy.special import dawsn
import bounce
from bounce import CXY, VT


class Response(bounce.Orbits):
    def __init__(self, geo, **kw):
        super().__init__(geo, **kw)
        nz = geo["nz"]; self.nz = nz
        zk = -np.pi + 2 * np.pi * np.arange(nz) / nz
        # volume-weighted overlap of the Fourier basis: Jm[m, m'] = int J e^{i (m' - m) z} dz
        e = np.exp(1j * np.outer(zk, self.m))                       # [z, m]
        self.E = e
        self.Jm = (np.conj(e).T * geo["J"][None, :]) @ e * (2 * np.pi / nz)
        self.V = geo["J"].sum() * 2 * np.pi / nz
        # n = 0 kernels (orbit average of e^{i m z}) for every lambda
        self.K0 = np.zeros((len(self.lam), nz), complex)
        for i, (l, kind) in enumerate(zip(self.lam, self.kind)):
            if kind == 0:
                z = -np.pi + 2 * np.pi * (np.arange(4800) + 0.5) / 4800
                dt = self.JBf(z) / np.sqrt(np.maximum(1 - l * self.Bfun(z), 1e-14))
            else:
                zl, zr = self._bounce(l); zc, zb = 0.5 * (zl + zr), 0.5 * (zr - zl)
                psi = -0.5 * np.pi + np.pi * (np.arange(1200) + 0.5) / 1200
                z = zc + zb * np.sin(psi)
                dt = self.JBf(z) * zb * np.cos(psi) / np.sqrt(np.maximum(1 - l * self.Bfun(z), 1e-14))
            self.K0[i] = (np.exp(1j * np.outer(self.m, z)) * dt[None, :]).sum(axis=1) / dt.sum()

    def Rmat(self, wt):
        pref = 0.5 * np.pi * CXY * VT * np.pi ** -1.5
        R = np.zeros((self.nz, self.nz), complex)
        nabs = np.abs(self.ns)
        for i in range(len(self.lam)):
            tau = 2 * np.pi / self.omhat[i]
            zeta = abs(wt) / (nabs * self.omhat[i])
            In = 2 * zeta ** 2 * np.sqrt(np.pi) * (2 * zeta * dawsn(zeta) - 1) - 2j * np.pi * np.sign(wt) * zeta ** 3 * np.exp(-zeta ** 2)
            K = self.K[i]                                               # [2 nmax, nz]
            if self.kind[i] == 0:
                A = np.sqrt(np.pi) * np.outer(np.conj(self.K0[i]), self.K0[i]) + (np.conj(K).T * In[None, :]) @ K
            else:
                A = 0.5 * np.sqrt(np.pi) * np.outer(np.conj(self.K0[i]), self.K0[i]) + (np.conj(K).T * In[None, :]) @ K
            R += pref * self.wl[i] * tau * A
        return R

    def source(self, k2lam2):
        g = self.geo; wz = g["J"] / g["J"].sum()
        Gxx, Gyy = (wz * g["gxx"]).sum(), (wz * g["gyy"]).sum()
        return k2lam2 * Gyy * (g["gxx"] / Gxx - g["gyy"] / Gyy)        # S(z_k)

    def solve(self, wt, k2lam2):
        """phi_1(z_k)/phi_0 (zero volume average) for Doppler-shifted frequency wt."""
        S = self.source(k2lam2)
        rhs = (np.conj(self.E).T * self.geo["J"][None, :]) @ S * (2 * np.pi / self.nz)
        A = self.Jm - self.Rmat(wt)
        c, *_ = np.linalg.lstsq(A, rhs, rcond=1e-10)
        phi1 = self.E @ c
        wz = self.geo["J"] / self.geo["J"].sum()
        return phi1 - (wz * phi1).sum()

    def loss_density(self, wt, phi1):
        """V <<wt^2 sum_n |phi_n|^2 delta(wt - n omega_b)>> for phi_1(z_k) (per |phi_0|^2 if phi1 is the ratio)."""
        G = self.G(bounce.zfft(phi1[None, :]), np.array([wt]))[0][0]
        return 2 * np.pi ** 2 * CXY * VT * np.pi ** -1.5 * G


if __name__ == "__main__":
    RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
    geo = bounce.geometry(f"{RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat")
    Rsp = Response(geo)
    R = Rsp.Rmat(1.0)
    print("sum rule R_00 / V (flute response = flute):", R[0, 0] / Rsp.V, " |R_m0| max for m != 0 relative:", np.abs(R[1:, 0]).max() / Rsp.V, np.abs(R[0, 1:] - Rsp.Jm[0, 1:]).max() / Rsp.V)
    for wt in (-3.0, -1.38, -0.65, 0.3, 1.5, 3.0):
        p = Rsp.solve(wt, 5000.0 * 0.002 ** 2)
        print(f"wt {wt:+.2f}: phi_1/phi_0 at z=-pi {p[0]:.4f}, z=0 {p[16]:.4f};  rms {np.sqrt((np.abs(p)**2*geo['J']).sum()/geo['J'].sum()):.4f};  loss density per |phi_0|^2 {Rsp.loss_density(wt, p):.3e}")

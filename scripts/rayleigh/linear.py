"""Linear fluid problem for the k_y = 1 vortex on a frozen jet (Rayleigh's equation with the flux-tube metric):
    (d_t + i k U)(psi'' - ke^2 psi) - i k U'' psi = -nu(kx) (psi'' - ke^2 psi),   ke^2 = k^2 <g^yy>/<g^xx>.
  eig_gene(psi0, hyp)        eigenvalues of the operator truncated to GENE's 63 radial modes, with GENE's sink
  ivp(psi0, a0, N, T, ...)   initial-value problem at N radial points (pseudo-spectral, RK4), any sink
"""
import numpy as np
from common import *
from measure import GXX, GYY, KX, fine

K = KYMIN
KE2 = K ** 2 * GYY / GXX


def pad(c, n):
    """FFT-ordered coefficients on NX-type grid (Nyquist empty) -> n-point FFT-ordered array."""
    m = len(c); h = m // 2
    o = np.zeros(n, complex)
    hh = min(h, n // 2)
    o[:hh] = c[:hh]; o[-(hh - 1):] = c[m - (hh - 1):]
    return o


def nu_gene(kx, hyp=HYP, dx=DX_G, ky=K):
    return hyp * ((dx * kx / 2) ** 4 + (DY_G * ky / 2) ** 4)


def eig_gene(psi0, hyp=HYP, nmode=NX, dx=None, k=K):
    """Dense operator on the retained modes |m| <= nmode/2 - 1; products dealiased exactly. Returns (lam, vecs, kx)
    with zeta ~ exp(lam t): omega = Im-part convention lam = -i omega + growth."""
    dx = LX / nmode if dx is None else dx
    h = nmode // 2
    m = np.r_[0:h, -(h - 1):0]
    kx = m * KXMIN
    n = 4 * nmode
    c0 = pad(psi0, n)
    kxn = np.fft.fftfreq(n, 1 / n) * KXMIN
    U = np.fft.ifft(1j * kxn * c0).real * n
    Upp = np.fft.ifft((1j * kxn) ** 3 * c0).real * n
    idx = m % n
    L = np.zeros((len(m), len(m)), complex)
    ke2 = k ** 2 * GYY / GXX
    for j in range(len(m)):
        e = np.zeros(n, complex); e[idx[j]] = 1.0          # psi = single mode
        psi = np.fft.ifft(e) * n
        zeta = np.fft.ifft(-(kxn ** 2 + ke2) * e) * n
        rhs = -1j * k * U * zeta + 1j * k * Upp * psi
        L[:, j] = (np.fft.fft(rhs) / n)[idx] - nu_gene(kx, hyp, dx, k) * (-(kx ** 2 + ke2)) * (np.arange(len(m)) == j)
    # zeta = -(kx^2+ke2) psi  ->  d_t psi = D^{-1} L psi
    D = -(kx ** 2 + ke2)
    A = L / D[:, None]
    lam, V = np.linalg.eig(A)
    return lam, V, kx


def ivp(psi0, a0, N=4096, T=400.0, dt=0.05, hyp=0.0, dx=DX_G, k=K, nout=4, proj=None):
    """a0: initial psi_1 coefficients (FFT order, any length <= N). Returns t, psi_1(t) at the NX GENE modes,
    and E(t) = sum (kx^2+ke^2)|psi|^2."""
    kx = np.fft.fftfreq(N, 1 / N) * KXMIN
    c0 = pad(psi0, N)
    U = np.fft.ifft(1j * kx * c0).real * N
    Upp = np.fft.ifft((1j * kx) ** 3 * c0).real * N
    D = -(kx ** 2 + k ** 2 * GYY / GXX)
    mask = np.abs(np.fft.fftfreq(N, 1 / N)) < N / 3          # 2/3 rule
    nu = nu_gene(kx, hyp, dx, k)
    z = D * pad(a0, N)

    def rhs(z):
        psi = z / D
        zr = np.fft.ifft(z); pr = np.fft.ifft(psi)
        return mask * np.fft.fft(-1j * k * U * zr + 1j * k * Upp * pr) - nu * z

    ns = int(round(T / dt))
    ts, out, E = [], [], []
    keep = np.r_[0:NX // 2, N - (NX // 2 - 1):N]
    for i in range(ns + 1):
        if i % nout == 0:
            psi = z / D
            ts.append(i * dt); out.append(psi[keep].copy()); E.append(float((-D * np.abs(psi) ** 2).sum()))
        k1 = rhs(z); k2 = rhs(z + 0.5 * dt * k1); k3 = rhs(z + 0.5 * dt * k2); k4 = rhs(z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.array(ts), np.array(out), np.array(E)


if __name__ == "__main__":
    np.set_printoptions(linewidth=220, precision=5, suppress=True)
    for tag in ("leg4_0", "leg4_2950"):
        d = np.load(f"{CACHE}/measured_{tag}.npz")
        psi0 = d["psi0"]; w = float(d["w"])
        xf = np.arange(4096) * LX / 4096
        U = fine(psi0, 4096, 1).real; S = fine(psi0, 4096, 2).real; Upp = fine(psi0, 4096, 3).real
        c = w / K
        print(f"== {tag}: measured omega {w:.4f} c {c:.1f}; U range {U.min():.0f} {U.max():.0f}; ke/q = {np.sqrt(KE2)/KXMIN:.3f}")
        for sgn in (1, -1):
            cr = [j for j in range(4096) if (U[j] - sgn * c) * (U[(j + 1) % 4096] - sgn * c) < 0]
            print("  layers c =", sgn * c, " x", xf[cr], " U'", S[cr], " U''", Upp[cr], " (U'' max %.2e)" % np.abs(Upp).max())
        for hyp in (HYP, 0.0):
            lam, V, kx = eig_gene(psi0, hyp)
            om = 1j * lam                     # psi ~ exp(-i om t): om = i lam
            o = np.argsort(-lam.real)
            print(f"  GENE-grid operator, hyp = {hyp}: least-damped eigenvalues (omega_r, growth):")
            print("   ", [(round(float(om[j].real), 4), float("%.3e" % lam[j].real)) for j in o[:12]])
            near = [j for j in range(len(lam)) if abs(abs(om[j].real) - w) < 0.12]
            print("    near measured omega:", [(round(float(om[j].real), 4), float("%.3e" % lam[j].real)) for j in sorted(near, key=lambda j: -lam[j].real)[:10]])

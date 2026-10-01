"""Measured quantities from the GENE flute potential: jet profile, the two counter-propagating k_y = 1 vortices
(two-frequency least-squares fit a_+(kx) e^{-i w t} + a_-(kx) e^{+i w t} in a window), fluid energies."""
import numpy as np
from common import *

M = metric(); GXX, GYY = M["gxx"], M["gyy"]
KX = np.fft.fftfreq(NX, 1 / NX) * KXMIN
KY = np.arange(NKY) * KYMIN
K2 = GXX * KX[:, None] ** 2 + GYY * KY[None, :] ** 2


def load(run, kind="phi_avg"):
    d = np.load(f"{CACHE}/gene_{run}.npz")
    return d["t"], d[kind] / CXY          # streamfunction psi = phi / C_xy


def energies(psi):
    """E_zon, E_nz of the fluid model, (1/2) sum <k_perp^2> |psi_k|^2 over all (kx, +-ky)."""
    e = K2 * np.abs(psi) ** 2
    return 0.5 * e[..., 0].sum(-1), e[..., 1:].sum((-1, -2))


def fine(c, nf=1024, d=0):
    """Fourier coefficients c(kx) (FFT order, NX) -> real-space profile on nf points, d-th x-derivative."""
    n = len(c)
    kx = np.fft.fftfreq(n, 1 / n) * KXMIN
    cc = np.zeros(nf, complex)
    h = n // 2
    cc[:h] = ((1j * kx) ** d * c)[:h]; cc[-(h - 1):] = ((1j * kx) ** d * c)[h + 1:]
    return np.fft.ifft(cc) * nf


def two_freq(t, s, w0=0.85, span=0.25, nw=401):
    """s(t, nx): fit a_p e^{-i w t} + a_m e^{+i w t}; scan w for the least residual. Returns w, a_p, a_m, resid."""
    best = None
    for w in np.linspace(w0 * (1 - span), w0 * (1 + span), nw):
        B = np.stack([np.exp(-1j * w * (t - t[0])), np.exp(1j * w * (t - t[0]))], 1)
        sol, res, *_ = np.linalg.lstsq(B, s, rcond=None)
        r = np.linalg.norm(s - B @ sol) / np.linalg.norm(s)
        if best is None or r < best[3]:
            best = (w, sol[0], sol[1], r)
    return best


if __name__ == "__main__":
    np.set_printoptions(linewidth=220, precision=3, suppress=True)
    for run, i0 in (("leg4", 0), ("leg4", 2950)):
        t, psi = load(run)
        sl = slice(i0, i0 + 51)
        w, ap, am, r = two_freq(t[sl], psi[sl, :, 1])
        xf = np.arange(1024) * LX / 1024
        U = fine(psi[sl, :, 0].mean(0), d=1).real; Upp = fine(psi[sl, :, 0].mean(0), d=3).real
        Ap, Am = np.abs(fine(ap)), np.abs(fine(am))
        c = w / KYMIN
        print(f"{run} t = {t[sl][0]:.1f}-{t[sl][-1]:.1f}: omega = {w:.4f}, c = {c:.1f}, residual {r:.3f}, U max/min {U.max():.0f}/{U.min():.0f}, c/Umax {c/U.max():.3f}")
        print("  |a+| peaks at x =", xf[np.argmax(Ap)], "max", Ap.max(), " U there", U[np.argmax(Ap)], " |a-| peaks at x =", xf[np.argmax(Am)], "max", Am.max(), "U there", U[np.argmax(Am)])
        for sgn, name in ((1, "+"), (-1, "-")):
            cr = [j for j in range(1024) if (U[j] - sgn * c) * (U[(j + 1) % 1024] - sgn * c) < 0]
            print(f"  critical layers of vortex {name}: x =", xf[cr], " U' =", fine(psi[sl, :, 0].mean(0), d=2).real[cr], " U'' =", Upp[cr], " U''rms", Upp.std())
        np.savez(f"{CACHE}/measured_{run}_{i0}.npz", t=t[sl], w=w, ap=ap, am=am, psi0=psi[sl, :, 0].mean(0), resid=r)

"""Pseudo-spectral 2D Euler model of the jet-vortex state in the GENE box (fluid limit of the notes, see common.py):
    d_t zeta + psi_x zeta_y - psi_y zeta_x = -nu_k S_k zeta,     zeta_k = -K_k psi_k,  K = <g^xx> kx^2 + <g^yy> ky^2,
psi = phi/C_xy.  Modes: |m| <= nx/2 - 1 in kx and 0 <= n <= nky - 1 in ky (GENE's set for nx = 64, nky = 16), products
dealiased by 3/2 padding, RK4.
Sinks (GENE's perpendicular hyperdiffusion, nu_k = hyp [(dx kx/2)^4 + (dy ky/2)^4], dx, dy GENE's grid spacings):
  sink = 'h'    : the sink acts on h (hyp_on_h = T).  Its charge moment damps (1 + lambda_D^2 k_perp^2) phi, so
                  the vorticity is damped at nu_k (1 + 1/(lambda_D^2 K)) - the Boltzmann part of h is damped too;
  sink = 'vort' : nu_k on the vorticity alone;
  sink = 'none'.
Separate radial / binormal coefficients: hypx, hypy (default: both = hyp).
Landau damping of the vortices (2 Oct 2026): landau = nu_L, the energy loss rate of every k_y >= 1 mode given by the
loss formula of the notes (linear, independent of amplitude); the vorticity of those modes is damped at nu_L / 2.
kin = 1 (2 Oct 2026, night): kinetic Landau damping of the k_y = 1 vortices instead of a friction.  Each evaluation
splits psi_1(x) into the two counter-propagating vortices, psi_+- = (psi_1 +- i dpsi_1/dt / w)/2 (w^2 from the Euler
tendency), takes the verified local loss density in each vortex's local frame, P(x) = 2 k^2 lambda_D^2 k^2 |psi_+-|^2
D(|wt|), wt = +-w - k U(x) [D from figures/landau/cache.npz], changes the vortex energy at the lab-frame rate
-(w/wt) P and deposits the absorbed momentum k P / wt in the jets (d zeta_0/dt = i k_x f).  No free constant.
zeps scales the zonal part of the initial condition as eps scales the rest.
Usage: python euler.py NAME [key=value ...]   keys: eps, T, dt, nx, nky, hyp, sink, dxfac (dx = DX_G*dxfac),
       dyfac, run (GENE run of the initial condition), frame, kind (phi_avg|phi_mid), tout
Output: cache/euler_NAME.npz  (t, Ezon, Enz, Eky[1..4], enstrophy, psi samples at GENE's 64 x 16 modes every tsamp)
"""
import sys
import time
import numpy as np
import scipy.fft as sf
from common import *
import measure as M

LAMD2 = 5000.0          # lambda_D^2 = debye2/2 in rho_ref^2


def run(name, eps=1.0, T=800.0, dt=0.02, nx=64, nky=16, hyp=HYP, sink="h", dxfac=1.0, dyfac=1.0, run="leg4", frame=0,
        kind="phi_avg", tout=0.5, tsamp=1.0, zonal_from=None, quiet=False, hypx=None, hypy=None, landau=0.0, zeps=1.0, kin=0, zcut=0.0):
    t0g, psig = M.load(run, kind)
    p_in = psig[frame].copy()
    if zonal_from is not None:                   # zonal part from another (run, frame)
        p_in[:, 0] = M.load(zonal_from[0], kind)[1][zonal_from[1]][:, 0]
    p_in[:, 1:] *= eps
    p_in[:, 0] *= zeps
    nxg, nyg = 3 * nx // 2, 3 * nky
    hx = nx // 2
    m = sf.fftfreq(nxg, 1.0 / nxg)
    n = np.arange(nyg // 2 + 1)
    kx, ky = m[:, None] * KXMIN, n[None, :] * KYMIN
    keep = (np.abs(m)[:, None] <= hx - 1) & (n[None, :] <= nky - 1)
    K = M.GXX * kx ** 2 + M.GYY * ky ** 2
    K[0, 0] = 1.0
    hypx = hyp if hypx is None else hypx
    hypy = hyp if hypy is None else hypy
    nu = hypx * (DX_G * dxfac * kx / 2) ** 4 + hypy * (DY_G * dyfac * ky / 2) ** 4
    if sink == "h":
        nu = nu * (1.0 + 1.0 / (LAMD2 * K))
    elif sink == "none":
        nu = 0.0 * nu
    nu = nu + 0.5 * landau * (n[None, :] >= 1)
    nu = np.where(keep, nu, 0.0)
    # initial condition: GENE's 64 x 16 coefficients into the padded array
    psi = np.zeros((nxg, nyg // 2 + 1), complex)
    h0 = NX // 2
    hh = min(h0, hx)
    nk = min(NKY, nky)
    psi[:hh, :nk] = p_in[:hh, :nk]
    psi[nxg - (hh - 1):, :nk] = p_in[NX - (hh - 1):, :nk]
    psi *= nxg * nyg
    z = np.where(keep, -K * psi, 0.0)
    z[0, 0] = 0.0
    ikx, iky = 1j * kx, 1j * ky
    LC = np.load(os.path.join(HERE, "..", "..", "figures", "landau", "cache.npz"))
    WG, DTH = LC["wg"], LC["Dth"]
    K2L2 = LAMD2 * KYMIN ** 2
    kx1 = kx[:, 0]; K1 = K[:, 1]; keep1 = keep[:, 1]; keep0 = keep[:, 0] & (m != 0)
    kstat = dict(n=0, P=0.0, Wlab=0.0, Uf=0.0, w=0.0)

    def kinetic(z, N):
        """kinetic damping tendency for zeta (k_y = 0 and 1 columns) given z and its Euler tendency N"""
        out = np.zeros_like(z)
        p1 = np.where(keep1, -z[:, 1] / K1, 0); dp1 = np.where(keep1, -N[:, 1] / K1, 0)
        s1 = np.sum(np.abs(p1) ** 2)
        if s1 == 0:
            return out
        w = np.sqrt(np.sum(np.abs(dp1) ** 2) / s1)
        x1 = np.fft.ifft(p1) / nyg; dx1 = np.fft.ifft(dp1) / nyg                     # psi_1(x), d psi_1/dt (x)
        U = np.fft.ifft(np.where(keep[:, 0], 1j * kx1 * (-z[:, 0] / np.where(K[:, 0] == 0, 1, K[:, 0])) * (m != 0), 0)).real / nyg
        g = np.zeros(nxg, complex); f = np.zeros(nxg)
        for sgn in (1, -1):
            ps = 0.5 * (x1 + sgn * 1j * dx1 / w)
            wt = sgn * w - KYMIN * U
            P = 2 * K2L2 * KYMIN ** 2 * np.abs(ps) ** 2 * np.interp(np.abs(wt), WG, DTH)
            dps = np.fft.ifft(1j * kx1 * np.fft.fft(ps))
            e = M.GXX * np.abs(dps) ** 2 + M.GYY * KYMIN ** 2 * np.abs(ps) ** 2
            beta = (sgn * w / wt) * P / (2 * np.maximum(e, 1e-30 * e.max()))
            g -= beta * ps
            f += KYMIN * P / wt
            kstat["P"] += P.mean(); kstat["Wlab"] += (sgn * w / wt * P).mean(); kstat["Uf"] += (U * KYMIN * P / wt).mean()
        kstat["n"] += 1; kstat["w"] += w
        out[:, 1] = np.where(keep1, -K1 * np.fft.fft(g) * nyg, 0)
        out[:, 0] = np.where(keep0, 1j * kx1 * np.fft.fft(f) * nyg, 0)
        return out

    def rhs(z):
        p = -z / K
        px = sf.irfft2(ikx * p, s=(nxg, nyg)); py = sf.irfft2(iky * p, s=(nxg, nyg))
        zx = sf.irfft2(ikx * z, s=(nxg, nyg)); zy = sf.irfft2(iky * z, s=(nxg, nyg))
        N = np.where(keep, -sf.rfft2(px * zy - py * zx), 0.0)
        if kin:
            # linear (Rayleigh) tendency of the k_y = 1 column: advection of the vortex by the jets and of the jets' vorticity by the vortex
            z01 = np.zeros_like(z); z01[:, 0] = z[:, 0]; z11 = np.zeros_like(z); z11[:, 1] = z[:, 1]
            p01, p11 = -z01 / K, -z11 / K
            a = lambda f, k_: sf.irfft2(k_ * f, s=(nxg, nyg))
            lin = -(a(p01, ikx) * a(z11, iky) - a(p01, iky) * a(z11, ikx) + a(p11, ikx) * a(z01, iky) - a(p11, iky) * a(z01, ikx))
            NL = np.where(keep, sf.rfft2(lin), 0.0)
            if zcut > 0:                         # diagnostic: no Euler response of zonal modes with |k_x| lambda_D > zcut
                N[:, 0] = np.where(np.abs(kx1) * np.sqrt(LAMD2) > zcut, 0.0, N[:, 0])
            return N - nu * z + kinetic(z, NL)
        return N - nu * z

    norm = 1.0 / (nxg * nyg) ** 2

    def diag(z):
        p2 = np.abs(z) ** 2 / K * norm          # K |psi|^2
        ez = 0.5 * p2[:, 0].sum()
        eky = p2[:, 1:].sum(0)
        ens = 0.5 * (np.abs(z[:, 0]) ** 2).sum() * norm + (np.abs(z[:, 1:]) ** 2).sum() * norm
        return ez, eky.sum(), eky[:4], ens

    ns = int(round(T / dt)); io = max(1, int(round(tout / dt))); isamp = max(1, int(round(tsamp / dt)))
    ts, EZ, EN, EK, ENS, tsm, PS, KS = [], [], [], [], [], [], [], []
    gi = np.r_[0:h0, nxg - (h0 - 1):nxg] if nx >= NX else None
    w0 = time.time()
    for i in range(ns + 1):
        if i % io == 0:
            ez, en, ek, ens = diag(z)
            ts.append(i * dt); EZ.append(ez); EN.append(en); EK.append(ek); ENS.append(ens)
            if kin and kstat["n"]:
                KS.append((i * dt, kstat["P"] / kstat["n"], kstat["Wlab"] / kstat["n"], kstat["Uf"] / kstat["n"], kstat["w"] / kstat["n"]))
                kstat.update(n=0, P=0.0, Wlab=0.0, Uf=0.0, w=0.0)
        if gi is not None and i % isamp == 0:
            ps = (-z / K)[gi][:, :min(NKY, nky)] / (nxg * nyg)
            tsm.append(i * dt); PS.append(np.concatenate([ps[:h0], np.zeros((1, ps.shape[1])), ps[h0:]]).astype(np.complex64))
        if not quiet and i % (50 * io * 20) == 0:
            print(f"  {name} t {i*dt:8.1f} Ezon {EZ[-1]:.5e} Enz {EN[-1]:.5e}  ({time.time()-w0:.0f} s)", flush=True)
        k1 = rhs(z); k2 = rhs(z + 0.5 * dt * k1); k3 = rhs(z + 0.5 * dt * k2); k4 = rhs(z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if not np.isfinite(z[1, 1]):
            print("blew up at", i * dt); break
    out = dict(t=np.array(ts), Ezon=np.array(EZ), Enz=np.array(EN), Eky=np.array(EK), ens=np.array(ENS),
               tsamp=np.array(tsm), psi=np.array(PS), t0=t0g[frame], kstat=np.array(KS),
               par=str(dict(kin=kin, zcut=zcut, eps=eps, zeps=zeps, landau=landau, hypx=hypx, hypy=hypy, T=T, dt=dt, nx=nx, nky=nky, hyp=hyp, sink=sink, dxfac=dxfac, dyfac=dyfac, run=run, frame=frame, kind=kind)))
    os.makedirs(CACHE, exist_ok=True)
    np.savez(f"{CACHE}/euler_{name}.npz", **out)
    return out


if __name__ == "__main__":
    kw = {}
    for a in sys.argv[2:]:
        k, v = a.split("=")
        kw[k] = v if k in ("sink", "run", "kind") else bool(float(v)) if k == "quiet" else (int(v) if k in ("nx", "nky", "frame", "kin") else float(v))
    o = run(sys.argv[1], **kw)
    print(sys.argv[1], "done: Enz(T)/Enz(0) =", o["Enz"][-1] / o["Enz"][0], " Ezon(T)/Ezon(0) =", o["Ezon"][-1] / o["Ezon"][0])

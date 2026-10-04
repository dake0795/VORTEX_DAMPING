r"""Rosenbluth-Hinton-type prediction of the final (phase-mixed) non-flute part of a jet harmonic from any checkpoint.

Linear k_y = 0 dynamics:  d_t h + v_par d_l h = (q F0/T) d_t phi  along orbits (no drift at k_y = 0, no drive).
Hence, on every orbit (eps, mu, and sigma = sign v_par for passing), I = <h>_orbit - F0 phibar is conserved, and the
final state has h = I + F0 phibar_final (orbit-constant).  Poisson (odd part, units q phi / T):
        (1 + k^2 lam^2 g^xx(l)) phi_f(l) = n[I](l) + rho0[phi_f](l).
Given a checkpoint (any time), compute I by orbit-averaging h at fixed mu (interpolating in v_par along the orbit),
solve for phi_f, and project its non-flute part on the static shape.  Usage: python rh_predict.py OUTDIR [KX_INDEX]
(OUTDIR holds checkpoint, field.dat, dipole_fix.dat, parameters.dat; the checkpoint's nx/nky are read from its header)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import bounce, response
OUT = sys.argv[1]; KX = int(sys.argv[2]) if len(sys.argv) > 2 else 1
geo = bounce.geometry(f"{OUT}/dipole_fix.dat"); B = geo["B"]; gxx = geo["gxx"]; J = geo["J"]; wz = J / J.sum(); nz = geo["nz"]
par = open(f"{OUT}/parameters.dat").read()
getp = lambda k: float(par.split(f"\n{k}")[1].split("=")[1].split()[0])
nv, nw, lv, lw, lx = int(getp("nv0")), int(getp("nw0")), getp("lv"), getp("lw"), getp("lx")
debye2 = getp("debye2") if "\ndebye2" in par else 5000.0
v = np.linspace(-lv, lv, nv); wv = np.full(nv, 2 * lv / (nv - 1))
for i, f in enumerate((17, 59, 43, 49)):
    wv[i] = wv[-1 - i] = f / 48 * 2 * lv / (nv - 1)
x, w = np.polynomial.laguerre.laggauss(nw); w = w * np.exp(x); fac = lw / w.sum(); mu = x * fac; wmu = w * fac
hd = open(f"{OUT}/checkpoint", "rb").read(46); tchk = np.frombuffer(hd[6:14], "<f8")[0]
nx, nky, nz_, nv_, nw_, ns = np.frombuffer(hd[22:46], "<i4")
G = np.memmap(f"{OUT}/checkpoint", dtype="<c16", mode="r", offset=46, shape=(ns, nw, nv, nz, nky, nx))
g = np.array(G[:, :, :, :, 0, KX])                                  # [s, mu, v, z]
kx = 2 * np.pi / lx * (KX if KX <= nx // 2 else KX - nx); kl2 = kx ** 2 * debye2 / 2   # lambda_D^2 = debye2/2 (pair plasma; rho = k^2 debye2 g^xx phi, checked)
# field at the checkpoint time: from the checkpoint itself via Poisson (odd charge / (k^2 lam^2 g^xx)); exact for k_y = 0
W3 = np.pi * B[None, None, :] * wmu[:, None, None] * wv[None, :, None]
F0 = np.pi ** -1.5 * np.exp(-v[None, :, None] ** 2 - mu[:, None, None] * B[None, None, :])
q = np.array([-1.0, 1.0]); rho = sum(q[s] * (W3 * g[s]).sum((0, 1)) for s in range(2))
phi = rho / (2 * kl2 * gxx)
ho = 0.5 * ((g[1] + phi * F0) - (g[0] - phi * F0))                 # odd part of h, units of phi
# ---- orbit average at fixed mu: for grid point (iw, iv, iz), energy eps = v^2 + mu B(z); along the orbit at z',
#      v' = sigma sqrt(eps - mu B(z')); weight J B dz / |v'|.  Passing: one direction sigma; trapped: both, in the well.
JB = J * B
def orbit_avg_field(F, iw):
    """F[v, z] at fixed mu index iw -> orbit average at every (v, z) grid point."""
    out = np.zeros((nv, nz), complex)
    for iz in range(nz):
        for iv in range(nv):
            e = v[iv] ** 2 + mu[iw] * B[iz]; vp2 = e - mu[iw] * B
            ok = vp2 > 1e-10
            if not ok.all():                                    # trapped: the well containing iz
                reg = np.zeros(nz, bool); j = iz
                while ok[j % nz] and not reg[j % nz]: reg[j % nz] = True; j += 1
                j = iz - 1
                while ok[j % nz] and not reg[j % nz]: reg[j % nz] = True; j -= 1
                ok = reg; sig = (1, -1)
            else:
                sig = (1,) if v[iv] >= 0 else (-1,)
            vp = np.sqrt(np.maximum(vp2, 1e-10)); wt = np.where(ok, JB / vp, 0.0)
            acc = 0
            for sgn in sig:
                vals = np.array([np.interp(sgn * vp[k], v, F[:, k].real) + 1j * np.interp(sgn * vp[k], v, F[:, k].imag) for k in range(nz)])
                acc = acc + (wt * vals).sum() / wt.sum()
            out[iv, iz] = acc / len(sig)
    return out
H = np.zeros_like(ho); P = np.zeros_like(ho)
phiF = phi[None, None, :] * F0                                   # F0 phi(z)
for iw in range(nw):
    H[iw] = orbit_avg_field(ho[iw], iw)
    # F0 phibar: F0 is constant on an orbit (function of eps), so <F0 phi>_orbit = F0 phibar
    P[iw] = orbit_avg_field(phiF[iw], iw)
I = H - P
nI = (W3 * I).sum((0, 1))
Rs = response.Response(geo); R0 = Rs.Rmat(1e-6); E = Rs.E; dz = 2 * np.pi / nz
proj = lambda f: (np.conj(E).T * J[None, :]) @ f * dz
# solve (1 + kl2 gxx) phi_f - rho0[phi_f] = n[I]  in weak form:  (Jm + kl2 M_gxx - R0) c = proj(nI)
Mg = (np.conj(E).T * (J * gxx)[None, :]) @ E * dz
c = np.linalg.solve(Rs.Jm + kl2 * Mg - R0, proj(nI)); phi_f = E @ c
Gxx = (wz * gxx).sum(); pf0 = (wz * phi_f).sum()
def nonflute_coef(p):
    c_s = np.linalg.lstsq(Rs.Jm - R0, -kl2 * proj(gxx - Gxx) * (wz * p).sum(), rcond=1e-10)[0]; ps = E @ c_s; ps -= (wz * ps).sum()
    p1 = p - (wz * p).sum(); return np.vdot(ps * wz, p1).real / np.vdot(ps * wz, ps).real
print(f"checkpoint t {tchk:.2f}, k_x index {KX} (k_x lam_D {np.sqrt(kl2):.3f})")
print(f"   non-flute factor of the field NOW (from the checkpoint):        {nonflute_coef(phi):+.3f}")
print(f"   predicted FINAL factor after phase mixing (I conserved):        {nonflute_coef(phi_f):+.3f}")
print(f"   flute level final / now: {abs(pf0)/abs((wz*phi).sum()):.4f};  orbit constancy now: rms(h - <h>)/rms(h) = {np.sqrt((W3*np.abs(ho-H)**2/F0).sum()/(W3*np.abs(ho)**2/F0).sum()):.4f}")

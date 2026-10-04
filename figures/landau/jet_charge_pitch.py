r"""Where in velocity space is the charge of a jet harmonic?  (4 Oct 2026)
Linear-run checkpoint (zonal_static/m1_chpt/out/checkpoint: k_y = 0, k_x = +-1 of the reduced box, stationary).
GENE grids: v_par uniform on [-lv, lv] (Simpson end weights), mu Gauss-Laguerre scaled to [0, lw] (gau_lag, asreg);
F0 = pi^-3/2 exp(-v^2 - mu B), int d^3v = pi B int dv dmu.  g = delta f;  h = g + q phi F0.
1. Normalisation check: rho(z) = sum_s q_s n[g_s](z) against phi(z) gxx(z) (exact zonal Poisson: k^2 lam^2 gxx phi = rho).
2. Orbit constancy of h/F0 (stationary state): spread of h/F0 within (eps, lambda) bins across z.
3. The orbit-averaged deposit I = <h> - F0 phibar as a function of pitch angle: uniform (Maxwellian deposit -> static
   theory, factor 1) or not?  4. Non-flute potential from the non-uniform part, with the orbit-averaged response."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce
OUT = f"{vf.RB}/rb_c0p3_draintest/zonal_static/m1_chpt/out"
geo = bounce.geometry(f"{OUT}/dipole_fix.dat"); B = geo["B"]; gxx = geo["gxx"]; J = geo["J"]; wz = J / J.sum(); nz = geo["nz"]
nv, nw, lv, lw = 84, 32, 4.5, 24.0
v = np.linspace(-lv, lv, nv); wv = np.full(nv, 2 * lv / (nv - 1))
for i, f in enumerate((17, 59, 43, 49)):
    wv[i] = wv[-1 - i] = f / 48 * 2 * lv / (nv - 1)
x, w = np.polynomial.laguerre.laggauss(nw); w = w * np.exp(x); fac = lw / w.sum(); mu = x * fac; wmu = w * fac
h = open(f"{OUT}/checkpoint", "rb").read(46); t, dt = np.frombuffer(h[6:22], "<f8"); dims = np.frombuffer(h[22:46], "<i4")
nx, nky, nz_, nv_, nw_, ns = dims; assert (nv_, nw_, nz_) == (nv, nw, nz)
G = np.memmap(f"{OUT}/checkpoint", dtype="<c16", mode="r", offset=46, shape=(ns, nw, nv, nz, nky, nx))
g = np.array(G[:, :, :, :, 0, 1])                                   # [s, mu, v, z], k_x = +1
# field at the checkpoint time from field.dat (last frame)
N = 3 * nz; FR = 16 + 4 + N * 16 + 4; b = open(f"{OUT}/field.dat", "rb").read()
phi = np.frombuffer(b[-FR:][20:20 + N * 16], "<c16").reshape((3, nz), order="F")[1]
tf = np.frombuffer(b[-FR:][4:12], "<f8")[0]; print(f"checkpoint t {t:.2f}, field t {tf:.2f}")
F0 = np.pi ** -1.5 * np.exp(-v[None, :, None] ** 2 - mu[:, None, None] * B[None, None, :])   # [mu, v, z]
W3 = np.pi * B[None, None, :] * wmu[:, None, None] * wv[None, :, None]
print("int F0 d^3v (should be 1):", np.round((W3 * F0).sum((0, 1))[[0, 8, 16, 24]], 4))
q = np.array([-1.0, 1.0])
rho = sum(q[s] * (W3 * g[s]).sum((0, 1)) for s in range(2))
ratio = rho / (phi * gxx)
print("rho/(phi gxx) along z (should be constant = k^2 lam^2 x normalisation):", np.round(ratio.real[::4], 5), " imag", np.round(np.abs(ratio.imag).max(), 6))
kl2 = (2 * np.pi / vf.LX) ** 2 * 5000.0; print("k^2 lam^2 =", kl2)
# h for the positrons (odd part: (h_p - h_e)/2)
hp = g[1] + phi[None, None, :] * F0; he = g[0] - phi[None, None, :] * F0; ho = 0.5 * (hp - he)   # odd (charge) part
X = ho / F0                                                         # h/F0, dimensionless (units of phi)
eps = v[None, :, None] ** 2 + mu[:, None, None] * B[None, None, :]; lam = mu[:, None, None] / np.maximum(eps, 1e-12)
lam = np.broadcast_to(lam, X.shape); eps = np.broadcast_to(eps, X.shape)
# orbit-time weight at each z for a particle: dt ~ dl/|v_par| ~ J B dz / |v_par|; in velocity space the phase-space weight
# W3*F0 already contains it (Liouville), so bin averages with weight W3*F0 are orbit averages over the population.
wt = np.broadcast_to(W3 * F0, X.shape)
lb = np.linspace(0, 1 / B.min(), 25); eb = np.array([0, 0.5, 1.0, 1.5, 2.5, 4.0, 25])
il = np.clip(np.digitize(lam, lb) - 1, 0, len(lb) - 2); ie = np.clip(np.digitize(eps, eb) - 1, 0, len(eb) - 2)
Xm = np.zeros((len(eb) - 1, len(lb) - 1), complex); Ws = np.zeros_like(Xm.real); spread = np.zeros_like(Ws)
for a in range(len(eb) - 1):
    for c in range(len(lb) - 1):
        m = (ie == a) & (il == c)
        if wt[m].sum() == 0: continue
        Xm[a, c] = (wt[m] * X[m]).sum() / wt[m].sum(); Ws[a, c] = wt[m].sum()
        spread[a, c] = np.sqrt((wt[m] * np.abs(X[m] - Xm[a, c]) ** 2).sum() / wt[m].sum()) / max(abs(Xm[a, c]), 1e-300)
phif = (wz * phi).sum()
print(f"\nflute phi {abs(phif):.4e}; h/F0 orbit-binned, in units of the flute phi (rows: energy bins {eb[:-1]}, cols: lambda bins)")
np.set_printoptions(linewidth=250, precision=3, suppress=True)
print("Re(h/F0)/phi_f:\n", (Xm / phif).real[:, ::2])
print("relative spread within bins (orbit constancy):\n", spread[:, ::2])
print("population weight per bin:\n", (Ws / Ws.sum())[:, ::2])
np.savez(os.path.join(HERE, "cache_jetcharge.npz"), Xm=Xm / phif, Ws=Ws, lb=lb, eb=eb, rho=rho, phi=phi, ratio=ratio)

# ---------------------------------------------------------------- 4. static part + part from the pitch-angle structure
import response
Rs = response.Response(geo)
R0 = Rs.Rmat(1e-6)                                  # orbit-averaged (Boltzmann in phibar) response, weak form
E = Rs.E; dz = 2 * np.pi / nz
proj = lambda f: (np.conj(E).T * J[None, :]) @ f * dz
nh = (W3 * ho).sum((0, 1))                          # n[h](z), units of phi
cPhi = np.linalg.lstsq(E, phi, rcond=None)[0]       # Fourier coefficients of the measured phi(z)
rhs_I = proj(nh) - R0 @ cPhi                        # weak form of n[I] = n[h] - rho0[phi]
Gxx = (wz * gxx).sum()
rhs_s = -kl2 * proj(gxx - Gxx) * phif               # static source
def solve(rhs):
    rhs = rhs - (Rs.Jm[:, 0] / Rs.Jm[0, 0]) * rhs[0]           # remove the constant (solvability)
    c = np.linalg.lstsq(Rs.Jm - R0, rhs, rcond=1e-10)[0]; p = E @ c; return p - (wz * p).sum()
p_s = solve(rhs_s); p_I = solve(rhs_I - (rhs_s + kl2 * proj(gxx) * phif))   # extra: n[I] minus its static-equivalent part
p_meas = phi - phif
def coef(p, ref): return np.vdot(ref * wz, p) / np.vdot(ref * wz, ref)
print("\nnon-flute potential, projected on the static shape:")
print(f"   static theory          1.000")
print(f"   from pitch-angle part  {coef(p_I, p_s).real:+.3f}   (shape match {abs(np.vdot(p_I*wz, p_s))**2/(np.vdot(p_I*wz,p_I).real*np.vdot(p_s*wz,p_s).real):.3f})")
print(f"   sum                    {coef(p_s + p_I, p_s).real:+.3f}")
print(f"   measured               {coef(p_meas, p_s).real:+.3f}   (residual of sum vs measured: {np.sqrt((wz*np.abs(p_s+p_I-p_meas)**2).sum()/(wz*np.abs(p_meas)**2).sum()):.3f})")

r"""As jet_charge_tperp.py, but fitting the deposit I/F0 = h/F0 - phibar(lambda) (in units of phi_f), with phibar the
orbit average of the measured phi(z) for each grid point's lambda (passing: whole line; trapped: the well that contains
the point, weight J B / sqrt(1 - lambda B))."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "jet_charge_pitch.py")).read().replace("print(", "(lambda *a, **k: None)("))
JB = J * B
def phibar(lmb, iz):
    ok = 1 - lmb * B > 0
    if not ok.all():                               # trapped: contiguous region around iz (periodic)
        reg = np.zeros(nz, bool); j = iz
        while ok[j % nz] and not reg[j % nz]: reg[j % nz] = True; j += 1
        j = iz - 1
        while ok[j % nz] and not reg[j % nz]: reg[j % nz] = True; j -= 1
        ok = reg
    wgt = np.where(ok, JB / np.sqrt(np.maximum(1 - lmb * B, 1e-6)), 0.0)
    return (wgt * phi).sum() / wgt.sum()
PB = np.zeros(X.shape, complex)
for iz in range(nz):
    for iw in range(nw):
        for iv in range(nv):
            PB[iw, iv, iz] = phibar(lam[iw, iv, iz], iz)
Yi = ((X - PB) / phif).real                         # I/F0 in units of phi_f
mu3 = np.broadcast_to(mu[:, None, None], X.shape); Wt = wt
for cols in (("1", "e", "e2"), ("1", "e", "e2", "mu")):
    M = np.array([{"1": np.ones_like(Yi), "e": eps, "e2": eps ** 2, "mu": mu3}[c].ravel() for c in cols]).T
    sw = np.sqrt(Wt.ravel()); coef, *_ = np.linalg.lstsq(M * sw[:, None], Yi.ravel() * sw, rcond=None)
    res = Yi.ravel() - M @ coef; yb = (Wt.ravel() * Yi.ravel()).sum() / Wt.sum()
    print(f"fit of I/F0 on {cols}: R^2 = {1 - (Wt.ravel()*res**2).sum()/(Wt.ravel()*(Yi.ravel()-yb)**2).sum():.3f}, coefficients {np.round(coef, 5)}")
dmu = coef[-1]
p_T = solve(proj(dmu * phif / B))
# check: the deposit computed this way reproduces the full pitch-angle part
nI = (W3 * (X - PB) * F0).sum((0, 1))
p_I2 = solve(proj(nI) - (rhs_s + kl2 * proj(gxx) * phif))
pr = lambda p: np.vdot(p_s * wz, p).real / np.vdot(p_s * wz, p_s).real
print(f"\nprojected on the static shape: T_perp model (d = {dmu:.5f}) {pr(p_T):+.3f};  deposit from I = h - F0 phibar {pr(p_I2):+.3f};  earlier full part {pr(p_I):+.3f}")

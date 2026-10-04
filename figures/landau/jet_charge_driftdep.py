r"""Model of the jets' deposit as proportional to the orbit-averaged magnetic drift:  I/F0 = c(eps) + kappa <v_d>_orbit,
v_d = (mu B + 2 v_par^2) K_y / B (GENE's drift).  <v_d>_orbit is evaluated per grid point with the orbit weights of
jet_charge_tperp2.py (passing: whole line; trapped: the well).  Fit, then the extra non-flute part it gives."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "jet_charge_tperp2.py")).read().split("mu3 = np.broadcast_to")[0].replace("print(", "(lambda *a, **k: None)("))
raw = [l for l in open(f"{OUT}/dipole_fix.dat").read().split("/")[1].split("\n") if l.strip()]
dg = np.array([[float(x) for x in l.split()] for l in raw]); Ky = dg[:, 7] / 0.92640709665123333
def orbit_avg(fz, lmb, iz):
    ok = 1 - lmb * B > 0
    if not ok.all():
        reg = np.zeros(nz, bool); j = iz
        while ok[j % nz] and not reg[j % nz]: reg[j % nz] = True; j += 1
        j = iz - 1
        while ok[j % nz] and not reg[j % nz]: reg[j % nz] = True; j -= 1
        ok = reg
    wgt = np.where(ok, JB / np.sqrt(np.maximum(1 - lmb * B, 1e-6)), 0.0)
    return (wgt * fz).sum() / wgt.sum()
VD = np.zeros(X.shape)
for iz in range(nz):
    for iw in range(nw):
        for iv in range(nv):
            e = eps[iw, iv, iz]; l = lam[iw, iv, iz]
            VD[iw, iv, iz] = orbit_avg(e * (l * B + 2 * (1 - l * B)) * Ky / B, l, iz)   # mu B + 2 v_par^2 = eps (lam B + 2(1 - lam B))
pr = lambda p: np.vdot(p_s * wz, p).real / np.vdot(p_s * wz, p_s).real
Wv = wt.ravel(); sw = np.sqrt(Wv)
for name, cols in (("energy only", [np.ones_like(eps), eps, eps ** 2]),
                   ("energy + <v_d>", [np.ones_like(eps), eps, eps ** 2, VD]),
                   ("energy + <v_d> + eps<v_d>", [np.ones_like(eps), eps, eps ** 2, VD, eps * VD])):
    M = np.array([c.ravel() for c in cols]).T
    coef, *_ = np.linalg.lstsq(M * sw[:, None], Yi.ravel() * sw, rcond=None)
    res = Yi.ravel() - M @ coef; yb = (Wv * Yi.ravel()).sum() / Wv.sum()
    r2 = 1 - (Wv * res ** 2).sum() / (Wv * (Yi.ravel() - yb) ** 2).sum()
    fitI = (M @ coef).reshape(X.shape) * phif
    p = solve(proj((W3 * fitI * F0).sum((0, 1))) - (rhs_s + kl2 * proj(gxx) * phif))
    print(f"{name:28s} R^2 {r2:.3f} -> extra non-flute {pr(p):+.3f}   coef {np.round(coef, 5)}")
print("measured extra (exact deposit): +0.661 / +0.689")

r"""Energy-dependent anisotropy models of the jets' deposit I/F0 (fit on the grid) and the extra non-flute part each gives."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
exec(open("jet_charge_tperp2.py").read().split("dmu = coef[-1]")[0].replace("print(", "(lambda *a, **k: None)("))
pr = lambda p: np.vdot(p_s * wz, p).real / np.vdot(p_s * wz, p_s).real
basis = {"1": np.ones_like(Yi), "e": eps, "e2": eps ** 2, "mu": mu3, "mue": mu3 * eps, "mue2": mu3 * eps ** 2, "mu2": mu3 ** 2}
# density of each basis function times F0 along z (exact quadrature on the grid)
for cols in (("1", "e", "e2", "mu"), ("1", "e", "e2", "mue"), ("1", "e", "e2", "mu", "mue"), ("1", "e", "e2", "mu", "mue", "mu2")):
    M = np.array([basis[c].ravel() for c in cols]).T; sw = np.sqrt(Wt.ravel())
    coef, *_ = np.linalg.lstsq(M * sw[:, None], Yi.ravel() * sw, rcond=None)
    res = Yi.ravel() - M @ coef; yb = (Wt.ravel() * Yi.ravel()).sum() / Wt.sum()
    r2 = 1 - (Wt.ravel() * res ** 2).sum() / (Wt.ravel() * (Yi.ravel() - yb) ** 2).sum()
    fitI = (M @ coef).reshape(X.shape) * phif
    nI = (W3 * fitI * F0).sum((0, 1)); p = solve(proj(nI) - (rhs_s + kl2 * proj(gxx) * phif))
    print(f"{str(cols):38s} R^2 {r2:.3f}  -> extra non-flute {pr(p):+.3f}   coef {np.round(coef, 5)}")

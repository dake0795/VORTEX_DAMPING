r"""One-parameter model of the jets' extra non-flute part: a perpendicular-temperature perturbation.
Regress h/F0 (odd part, units of phi_f) on [1, eps, eps^2, mu] over the whole velocity grid and z (phase-space weight
W3 F0); the mu term, I = d phi_f mu F0, has density n[I](z) = d phi_f / B(z), which varies along the line; solve its
non-flute response with the orbit-averaged operator and project on the static shape (to compare with 0.689)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
src = open(os.path.join(HERE, "jet_charge_pitch.py")).read()
exec(src.replace("print(", "(lambda *a, **k: None)("))
Y = (X / phif).real - 1; Wt = wt
mu3 = np.broadcast_to(mu[:, None, None], X.shape)
for cols in (("1", "e", "e2"), ("1", "e", "e2", "mu")):
    M = np.array([{"1": np.ones_like(Y), "e": eps, "e2": eps ** 2, "mu": mu3}[c].ravel() for c in cols]).T
    sw = np.sqrt(Wt.ravel())
    coef, *_ = np.linalg.lstsq(M * sw[:, None], Y.ravel() * sw, rcond=None)
    res = Y.ravel() - M @ coef; ybar = (Wt.ravel() * Y.ravel()).sum() / Wt.sum()
    r2 = 1 - (Wt.ravel() * res ** 2).sum() / (Wt.ravel() * (Y.ravel() - ybar) ** 2).sum()
    print(f"fit on {cols}: R^2 = {r2:.3f}, coefficients {np.round(coef, 5)}")
dmu = coef[-1]
src_T = proj(dmu * phif / B)                      # weak form of n[I_mu]
p_T = solve(src_T)
print(f"\nT_perp model (d = {dmu:.5f}): non-flute part projected on the static shape {coef_s if False else 0:.0f}", end="")
print(f"\r  T_perp model, one parameter d = {dmu:.5f}:  {np.vdot(p_s*wz, p_T).real/np.vdot(p_s*wz, p_s).real:+.3f}  (shape match "
      f"{abs(np.vdot(p_T*wz, p_s))**2/(np.vdot(p_T*wz,p_T).real*np.vdot(p_s*wz,p_s).real):.3f});  full pitch-angle part {np.vdot(p_s*wz, p_I).real/np.vdot(p_s*wz, p_s).real:+.3f}")

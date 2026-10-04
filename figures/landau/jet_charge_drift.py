r"""Is the pitch-angle structure of the jets' distribution set by the magnetic drift?  Bin the drift factor
(mu B + 2 v_par^2) K_y / B (GENE's drift, K_y = dB/dx / C_xy) with the same orbit weights as h/F0 (jet_charge_pitch.py),
and regress h/F0 - 1 on [1, eps, eps^2, d(eps,lambda)] with d the orbit-averaged drift factor."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "jet_charge_pitch.py")).read().split("# ---------------------------------------------------------------- 4.")[0].replace("print(", "(lambda *a, **k: None)("))
geof = np.loadtxt(f"{OUT}/dipole_fix.dat", skiprows=16) if False else None
raw = [l for l in open(f"{OUT}/dipole_fix.dat").read().split("/")[1].split("\n") if l.strip()]
d = np.array([[float(x) for x in l.split()] for l in raw])
Bz, dBdx = d[:, 6], d[:, 7]; Ky = dBdx / 0.92640709665123333
drift = (mu[:, None, None] * Bz[None, None, :] + 2 * v[None, :, None] ** 2) * Ky[None, None, :] / Bz[None, None, :]
drift = np.broadcast_to(drift, X.shape)
Dm = np.zeros_like(Ws)
for a in range(len(eb) - 1):
    for c in range(len(lb) - 1):
        m = (ie == a) & (il == c)
        if wt[m].sum() == 0: continue
        Dm[a, c] = (wt[m] * drift[m]).sum() / wt[m].sum()
Ec = np.zeros_like(Ws)
for a in range(len(eb) - 1):
    for c in range(len(lb) - 1):
        m = (ie == a) & (il == c)
        if wt[m].sum(): Ec[a, c] = (wt[m] * eps[m]).sum() / wt[m].sum()
Y = (Xm / phif).real - 1; ok = Ws > 0
np.set_printoptions(linewidth=250, precision=3, suppress=True)
print("orbit-averaged drift factor / eps (rows energy, cols lambda):\n", (Dm / np.maximum(Ec, 1e-9))[:, ::2])
for cols, name in ((["1", "e", "e2"], "energy only"), (["1", "e", "e2", "d"], "energy + drift"), (["1", "e", "e2", "e*lam"], "energy + eps*lambda")):
    M = []
    for cname in cols:
        M.append({"1": np.ones_like(Y), "e": Ec, "e2": Ec ** 2, "d": Dm, "e*lam": Ec * 0.5 * (lb[:-1] + lb[1:])[None, :]}[cname][ok])
    M = np.array(M).T; wgt = np.sqrt(Ws[ok])
    coef, *_ = np.linalg.lstsq(M * wgt[:, None], Y[ok] * wgt, rcond=None)
    res = Y[ok] - M @ coef; r2 = 1 - (Ws[ok] * res ** 2).sum() / (Ws[ok] * (Y[ok] - (Ws[ok] * Y[ok]).sum() / Ws[ok].sum()) ** 2).sum()
    print(f"{name:22s}: R^2 = {r2:.3f}   coefficients {np.round(coef, 4)}")

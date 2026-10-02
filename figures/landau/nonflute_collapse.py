import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gyy = (wz * geo["gyy"]).sum(); sh = geo["gyy"] - Gyy; sh -= (wz * sh).sum()
NXF = predict.NXF
Tr = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)["orig"]
R2 = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
bins = np.linspace(-1.5, 1.5, 16); out = {}
for run, tc in (("orig", 2000.0), ("orig", 4000.0), ("orig", 6100.0), ("orig", 8000.0), ("orig", 9800.0), ("hxoff", 1500.0), ("hxoff", 3900.0)):
    t, p0, p1 = predict.window_field(run, tc)
    if run == "orig": i = np.argmin(abs(Tr[:, 0] - tc)); w0 = (Tr[i, 1], Tr[i, 2])
    else: ww = R2[np.argmin(abs(R2[:, 0] - tc)), 1]; w0 = (ww, -ww)
    w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], w0)
    Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, vf.NX, vf.NZ)
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    row = []
    for kk in range(2):
        ph = predict.to_x(a[kk]); ph0 = (ph * wz[None, :]).sum(1); m = ph - ph0[:, None]
        r = (m * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum() / ph0
        wt = w[kk] - vf.KYMIN * U; s = np.abs(ph0) > 0.3 * np.abs(ph0).max()
        rr = np.conj(r) if w[kk] < 0 else r               # mirror the negative-frequency vortex
        wt2 = -wt if w[kk] < 0 else wt
        row.append([np.mean(rr[s & (wt2 >= a_) & (wt2 < b_)]) if (s & (wt2 >= a_) & (wt2 < b_)).sum() else np.nan for a_, b_ in zip(bins[:-1], bins[1:])])
    out[(run, tc)] = row
c = 0.5 * (bins[1:] + bins[:-1])
print("wt (mirrored for w<0):   " + " ".join(f"{x:+6.2f}" for x in c))
for key, row in out.items():
    for kk in range(2):
        print(f"{key[0]:5s} {key[1]:6.0f} v{kk}  real: " + " ".join(f"{v.real:+6.2f}" if np.isfinite(v) else "   -  " for v in row[kk]))

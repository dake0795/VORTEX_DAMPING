"""Does the growth of the jets' non-flute part (m = 1, 3) after the burst track the cumulative Landau loss of the vortices?
|phi1|(t) from jets_factor_history windows, against W(t) = int_{700}^t P_1 dt' with P_1 = -E_parallel(k_y = 1) (rb_c0p3_hxy)."""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import vortex_field as vf, predict, spec2d
wz = None
import bounce
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
ts, cs = [], []
for d in vf.RUNS["orig"]:
    c = 0
    for sp in ("electrons", "positrons"):
        t, cc = spec2d.by_ky(f"{d}/Spectral_2D_E_parallel_{sp}.dat"); c = c + cc
    ts.append(t); cs.append(c)
t = np.concatenate(ts); c = np.concatenate(cs); o = np.argsort(t); t, P = t[o], -c[o, 1]
keep = np.r_[True, np.diff(t) > 0]; t, P = t[keep], P[keep]
W = np.concatenate([[0], np.cumsum(0.5 * (P[1:] + P[:-1]) * np.diff(t))])
TC = np.arange(700, 3001, 100.0); A = {1: [], 3: []}; Wc = []
for tc in TC:
    tt, p0, p1 = predict.window_field("orig", float(tc), half=10.0)
    pz = p0.mean(0); pz1 = pz - (pz * wz[None, :]).sum(1)[:, None]
    for m in (1, 3):
        A[m].append(np.sqrt(np.sum(wz * np.abs(pz1[m]) ** 2)))
    Wc.append(np.interp(tc, t, W) - np.interp(700, t, W))
Wc = np.array(Wc)
print("    t     W (cum. Landau loss)   d|phi1| m=1   d|phi1| m=3     ratio m=1 dphi1/W")
for i, tc in enumerate(TC):
    d1 = A[1][i] - A[1][0]; d3 = A[3][i] - A[3][0]
    print(f"{tc:6.0f}   {Wc[i]:12.4e}        {d1:9.2f}     {d3:9.2f}      {d1/Wc[i] if Wc[i] else 0:.3e}")

"""Non-flute part of the vortex (projection on g^yy - <g^yy>, relative to phi_0) at its critical layers and in the jet,
for the reference run and the amplitude-reduced restarts (eps = 0.03: linear vortex on unchanged jets)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, predict
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gyy = (wz * geo["gyy"]).sum(); sh = geo["gyy"] - Gyy; sh -= (wz * sh).sum()
NXF = predict.NXF; xf = np.arange(NXF) * vf.LX / NXF
Tr = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)
def analyse(run, tc, w0):
    t, p0, p1 = predict.window_field(run, tc)
    if len(t) < 20: return
    w, _, res = vf.two_freq_fit(t, p1[:, :, vf.ZMID], w0)
    Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, vf.NX, vf.NZ)
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    k = 0; ph = predict.to_x(a[k]); ph0 = (ph * wz[None, :]).sum(1); meas = ph - ph0[:, None]
    A = (meas * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum(); r = A / ph0
    wt = w[k] - vf.KYMIN * U
    strong = np.abs(ph0) > 0.5 * np.abs(ph0).max()
    crit = strong & (np.abs(wt) < 0.15); mid = strong & (np.abs(wt) > 1.1)
    print(f"{run:5s} t {t.mean():7.1f}: w {w[0]:+.3f}, fit resid {res:.3f}, |phi0| max {np.abs(ph0).max():.3e};  A/phi0 at critical layers {np.mean(r[crit].real):+.3f}{np.mean(r[crit].imag):+.3f}i,"
          f"  in mid-jet {np.mean(r[mid].real):+.3f}{np.mean(r[mid].imag):+.3f}i")
for run, tcs in (("orig", (5950.0, 6100.0, 6300.0)), ("a003", (5950.0, 6000.0, 6100.0, 6300.0, 6600.0)), ("a03", (5950.0, 6100.0, 6300.0))):
    T = Tr[run] if run in Tr.files else Tr["orig"]
    for tc in tcs:
        i = np.argmin(abs(T[:, 0] - tc)); analyse(run, tc, (T[i, 1], T[i, 2]))

import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum(); sh = geo["gyy"] - Gyy; sh -= (wz * sh).sum()
L2 = 5000.0; k = vf.KYMIN; NXF = predict.NXF
Tr = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)["orig"]
R2 = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
print(" run    t     bump A/phi0   lam^2 zeta1/phi0 at crit   lam^2 d_x^2 phi0/phi0   ratio bump/(lam^2 zeta/phi0)")
for run, tc in (("orig", 2000.0), ("orig", 4000.0), ("orig", 6100.0), ("orig", 8000.0), ("orig", 9800.0), ("hxoff", 1500.0), ("hxoff", 3900.0), ("a003", 6100.0)):
    t, p0, p1 = predict.window_field(run, tc)
    if run != "hxoff": i = np.argmin(abs(Tr[:, 0] - tc)); w0 = (Tr[i, 1], Tr[i, 2])
    else: ww = R2[np.argmin(abs(R2[:, 0] - tc)), 1]; w0 = (ww, -ww)
    w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], w0)
    Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, vf.NX, vf.NZ)
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    for kk in range(1):
        a0 = (a[kk] * wz[None, :]).sum(1)
        ph0 = vf.to_x(a0, nxf=NXF); d2 = vf.to_x(a0, nxf=NXF, deriv=2)
        ph = predict.to_x(a[kk]); m = ph - (ph * wz[None, :]).sum(1)[:, None]
        r = (m * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum() / ph0
        wt = w[kk] - k * U; s = (np.abs(ph0) > 0.3 * np.abs(ph0).max()) & (np.abs(wt) < 0.1)
        zeta = L2 * (Gxx * d2 - Gyy * k ** 2 * ph0) / ph0
        print(f" {run:5s} {tc:6.0f}   {np.mean(r[s].real):+.3f}         {np.mean(zeta[s].real):+.3f}{np.mean(zeta[s].imag):+.3f}i            {np.mean((L2*d2/ph0)[s].real):+.3f}             {np.mean(r[s].real)/np.mean(zeta[s].real):+.3f}")

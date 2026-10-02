"""Predicted non-flute part of each vortex, phi_1(x, z) = lambda_D^2 [ d_x^2 phi_0 * chi_x(z; wt(x)) - k_y^2 phi_0 * chi_y(z; wt(x)) ],
chi_x / chi_y = response to the sources g^xx(l) - <g^xx>, g^yy(l) - <g^yy> (response.Response, Doppler-shifted frequency
wt = w - k_y U(x)), against the measured non-flute part of the fitted vortex patterns in rb_c0p3_hxoff."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat")
Rs = response.Response(geo); wz = geo["J"] / geo["J"].sum()
Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
WT = np.r_[-np.linspace(4.2, 0.02, 60), np.linspace(0.02, 4.2, 60)]
tab = {}
for nm, src in (("x", geo["gxx"] - Gxx), ("y", geo["gyy"] - Gyy), ("s", Gyy * (geo["gxx"] / Gxx - geo["gyy"] / Gyy))):
    Rs.source = lambda k2, s=src: k2 * s
    tab[nm] = np.array([Rs.solve(w, 1.0) for w in WT])        # [nw, nz], per unit lambda^2 (k2lam2 = 1 times source)
np.savez(os.path.join(HERE, "chi_tables.npz"), WT=WT, chix=tab["x"], chiy=tab["y"], chis=tab["s"])
def chi(nm, wt):          # wt[nx] -> [nx, nz]
    T = tab[nm]; i = np.clip(np.searchsorted(WT, wt), 1, len(WT) - 1); a = (wt - WT[i - 1]) / (WT[i] - WT[i - 1])
    return (1 - a)[:, None] * T[i - 1] + a[:, None] * T[i]
L2 = 5000.0; NXF = predict.NXF
R2 = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
for tc in (1500.0, 2700.0, 3900.0):
    t, p0, p1 = predict.window_field("hxoff", tc)
    w0 = R2[np.argmin(abs(R2[:, 0] - tc)), 1]
    w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], (w0, -w0))
    Mx = np.exp(-1j * np.outer(t - t.mean(), w))
    a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, vf.NX, vf.NZ)
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    for k in range(2):
        ph = predict.to_x(a[k]); ph0 = (ph * wz[None, :]).sum(axis=1)
        meas = ph - ph0[:, None]
        d2 = vf.to_x((a[k] * wz[None, :]).sum(axis=1), nxf=NXF, deriv=2)
        wt = w[k] - vf.KYMIN * U
        for sgn, form in ((1, "appendix: k^2 lamD^2 chi_s phi_0"), (1, "metric source on d_x^2 and d_y^2")):
            if form.startswith("appendix"):
                pred = L2 * vf.KYMIN ** 2 * ph0[:, None] * chi("s", wt)
            else:
                pred = sgn * L2 * (d2[:, None] * chi("x", wt) - vf.KYMIN ** 2 * ph0[:, None] * chi("y", wt))
            c = np.vdot(pred.ravel() * np.repeat(wz[None, :], NXF, 0).ravel(), meas.ravel()) / np.vdot(pred.ravel() * np.repeat(wz[None, :], NXF, 0).ravel(), pred.ravel())
            r = np.linalg.norm((meas - c * pred) * np.sqrt(wz)[None, :]) / np.linalg.norm(meas * np.sqrt(wz)[None, :])
            print(f"t {tc:.0f} vortex {k} {form}: |measured non-flute|/|flute| {np.linalg.norm(meas*np.sqrt(wz))/np.linalg.norm(ph0)/np.sqrt(1):.3e}; best scale of prediction {abs(c):.3f} at phase {np.degrees(np.angle(c)):+5.0f}; residual {r:.3f}")

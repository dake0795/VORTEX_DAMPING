"""Non-flute part of the jets (k_y = 0) per radial harmonic m, measured (window mean, Jacobian-weighted average removed)
against the static kinetic response phi_z1 - rho[phi_z1] = -q, q = -lambda_D^2 (g^xx(l) - <g^xx>) d_x^2 phi_z0
(response.Response at wt -> 0, conjugated per the convention), projected on its own leading shape along the line."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
vf.RUNS["c0p1"] = [f"{vf.RB}/rb_c0p1_hxy/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx = (wz * geo["gxx"]).sum(); L2 = 5000.0
R = response.Response(geo); R.source = lambda k2: k2 * (geo["gxx"] - Gxx)
chi = np.conj(R.solve(0.02, 1.0))                          # response per unit lambda^2 k_x^2 ... source = k2 * (gxx - <gxx>)
kxg = np.fft.fftfreq(vf.NX, 1.0 / vf.NX) * 2 * np.pi / vf.LX
for run, tc in [(a.split(":")[0], float(a.split(":")[1])) for a in sys.argv[1:]] or [("orig", 4000.0), ("orig", 8000.0), ("hxoff", 3900.0)]:
    t, p0, p1 = predict.window_field(run, tc)
    pz = p0.mean(0); pz0 = (pz * wz[None, :]).sum(1); pz1 = pz - pz0[:, None]
    print(f"\n{run} t {tc:.0f}:  m   kx lamD   |phi_z1|/|phi_z0|   measured/predicted (best complex scale)   shape match")
    for m in (1, 3, 5, 7, 9, 11, 13):
        for i in (m,):
            pred = (L2 * kxg[i] ** 2) * chi * pz0[i]            # d_x^2 -> -k_x^2 and the -q sign cancel: phi_1 = lam^2 k_x^2 chi phi_0
            meas = pz1[i]
            c = np.vdot(pred * wz, meas) / np.vdot(pred * wz, pred)
            res = np.sqrt(np.sum(wz * np.abs(meas - c * pred) ** 2) / np.sum(wz * np.abs(meas) ** 2))
            print(f"         {m:2d}   {kxg[i]*np.sqrt(L2):5.2f}      {np.sqrt(np.sum(wz*np.abs(meas)**2))/abs(pz0[i]):.3e}          {abs(c):.3f} at {np.degrees(np.angle(c)):+5.0f}               {1-res**2:.3f}")

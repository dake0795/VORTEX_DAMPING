r"""Hypothesis (4 Oct 2026): inside the cat's eye the vortex's charge is the jets' charge, stirred; it keeps the jets'
profile along the field line, c_J(x, l) = Q_z(x, l) / <Q_z>(x), Q_z = -lambda_D^2 g^xx(l) d_x^2 phi_z(x, l) (jets'
polarisation charge, non-flute part of the jets included).  Then the non-flute part of the vortex obeys
    phi_1 - rho[phi_1] = sigma(l),  sigma = [ zeta_bar c_J(l) - lambda_D^2 (g^xx phi_0'' - g^yy k^2 phi_0) ]_nonflute,
zeta_bar = lambda_D^2 (<g^xx> phi_0'' - <g^yy> k^2 phi_0)  (flute vorticity of the vortex).  With c_J = g^xx/<g^xx> this is
the appendix source.  Static response at the critical layer.  Tests:
 (1) the measured charge profile of the vortex, Q_v(l) = -lambda_D^2 (g^xx phi_v'' - g^yy k^2 phi_v) from the full 3D
     vortex, against c_J at the critical layers (shape match);
 (2) the predicted non-flute bump against the measured one."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import fieldgen as fg, bounce, response
G = fg.geometry(); w = G["w"]; gxx, gyy, Gxx, Gyy = G["gxx"], G["gyy"], G["Gxx"], G["Gyy"]
sh = gyy - Gyy; sh -= (w * sh).sum()
geo = bounce.geometry(f"{fg.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); Rs = response.Response(geo)
def static(src):
    Rs.source = lambda k2, s_=src: k2 * s_
    return np.conj(Rs.solve(0.02, 1.0))
L2, k, NXF = fg.LAM2, fg.KYMIN, 512
CASES = [("hxoff nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64, (3300.0, 4300.0)),
         ("ref nx64", [fg.RB + f"/rb_c0p3_hxy/leg_000{l}/out" for l in (3, 4, 5)], 64, (6100.0, 9600.0)),
         ("nx128_hxoff", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out"], 128, (3300.0,))]
for name, dirs, nx, tcs in CASES:
    R = fg.Run(dirs, nx)
    for tc in tcs:
        t, P = R.window(tc); p0, p1 = P[:, :, 0, :], P[:, :, 1, :]
        wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        pz = p0.mean(0)
        U = R.to_x((pz * w[None, :]).sum(1), NXF, 1).real / fg.CXY
        Qz = -L2 * gxx[None, :] * R.to_x(pz, NXF, 2).real                       # jets' charge [x, l]
        cJ = Qz / (Qz * w[None, :]).sum(1)[:, None]
        res_tab = []
        for kk in range(2):
            ph = R.to_x(a[kk], NXF); ph0 = (ph * w[None, :]).sum(1); nf = ph - ph0[:, None]
            d2 = R.to_x(a[kk], NXF, 2)
            Qv = -L2 * (gxx[None, :] * d2 - gyy[None, :] * k ** 2 * ph)                 # vortex charge [x, l] from the full 3D vortex
            cV = Qv / (Qv * w[None, :]).sum(1)[:, None]
            a0 = (a[kk] * w[None, :]).sum(1)
            p0d2 = R.to_x(a0, NXF, 2); p00 = R.to_x(a0, NXF)
            zb = L2 * (Gxx * p0d2 - Gyy * k ** 2 * p00)
            wt = wf[kk] - k * U
            sel = np.where((np.abs(ph0) > 0.3 * np.abs(ph0).max()) & (np.abs(wt) < 0.08))[0]
            for j in sel[:: max(1, len(sel) // 4)]:
                # (1) shape of the vortex charge vs the jets' charge, both minus 1 (non-flute part of the profile)
                dv, dj, dg = cV[j] - 1, cJ[j] - 1, gxx / Gxx - 1
                cos = lambda u, v: np.real(np.vdot(u * w, v)) / np.sqrt(np.real(np.vdot(u * w, u)) * np.real(np.vdot(v * w, v)))
                amp = lambda u, v: np.real(np.vdot(v * w, u)) / np.real(np.vdot(v * w, v))
                # (2) predicted non-flute part with the jets' profile (and with the plain metric, = appendix)
                loc = L2 * (gxx * p0d2[j] - gyy * k ** 2 * p00[j])
                out = {}
                for lab, c in (("jets", cJ[j]), ("metric", gxx / Gxx)):
                    sig = zb[j] * c - loc; sig = sig - (w * sig).sum()
                    Rs.source = lambda k2, s_=sig / p00[j]: k2 * s_
                    phi1 = np.conj(Rs.solve(0.02, 1.0))                                # per unit phi_0
                    out[lab] = (phi1 * w * sh).sum() / (w * sh ** 2).sum()
                Am = (nf[j] * w * sh).sum() / (w * sh ** 2).sum() / ph0[j]
                res_tab.append((wt[j], cos(dv, dj), amp(dv, dj), cos(dv, dg), Am, out["jets"], out["metric"]))
        print(f"\n{name} t {tc:.0f} (w {wf[0]:+.3f}): at critical-layer points")
        print("   wt      shape(vortex charge ~ jets' charge)  amp   shape(~ g^xx)   bump measured   predicted (jets' profile)   predicted (metric only)")
        for r in res_tab:
            print(f"  {r[0]:+.3f}          {r[1]:+.3f}                    {r[2]:5.2f}       {r[3]:+.3f}         {r[4].real:+.3f}            {r[5].real:+.3f}{r[5].imag:+.3f}i            {r[6].real:+.3f}")

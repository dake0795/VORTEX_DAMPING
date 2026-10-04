r"""Balance of the vortices with the absorbed power distributed by the MEASURED local dissipation function (from the lag
of their own non-flute part), instead of the theoretical D(wt) (4 Oct 2026).  Same windows as figure balance."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; gxx, gyy, Gxx, Gyy = G["gxx"], G["gyy"], G["Gxx"], G["Gyy"]
s_app = Gyy * (gxx / Gxx - gyy / Gyy)
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
L2, k = fg.LAM2, fg.KYMIN; K2L2 = L2 * k ** 2; nxf = 512
src = open(os.path.join(HERE, "drain_from_stress.py")).read(); exec(src[src.index("def tend"):src.index("def fdep")])
def terms(R, tc):
    t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
    wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
    psz = (p0.mean(0) * w[None, :]).sum(1); U = R.to_x(psz, nxf, 1).real
    out = dict(flu=0.0, Wth=0.0, Wme=0.0, Pth=0.0, Pme=0.0, E=0.0)
    for kk in range(2):
        a0 = (a[kk] * w[None, :]).sum(1)
        out["flu"] += float(np.sum(np.real(np.conj(psz) * tend(R, a0, Gxx, Gyy))))
        ph = R.to_x(a[kk], nxf); ph0 = (ph * w[None, :]).sum(1); ph1 = ph - ph0[:, None]
        wt = wf[kk] - k * U
        Dme = -wt * (np.imag(ph1 / ph0[:, None]) * s_app[None, :] * w[None, :]).sum(1) / K2L2
        Dth = np.interp(np.abs(wt), WG, DTH)
        strong = np.abs(ph0) > 0.05 * np.abs(ph0).max()
        Dme = np.where(strong, Dme, Dth)                     # where the vortex is negligible the ratio is noise
        for key, D in (("th", Dth), ("me", Dme)):
            Pd = 2 * K2L2 * k ** 2 * np.abs(ph0) ** 2 * D
            out["P" + key] += Pd.mean(); out["W" + key] += np.mean(wf[kk] / np.where(np.abs(wt) < 1e-6, 1e-6, wt) * Pd)
        dph0 = R.to_x(a0, nxf, 1); out["E"] += np.mean(Gxx * np.abs(dph0) ** 2 + Gyy * k ** 2 * np.abs(ph0) ** 2)
    return out
for key, dirs, nx, t0, t1 in (("nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64, 1200.0, 4500.0),
                              ("nx128", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128, 2600.0, 4300.0)):
    R = fg.Run(dirs, nx); C = np.load(os.path.join(HERE, "..", "balance", "cache.npz"))[key]
    print(f"\n{key}:   t     measured/P   stress/P  | damping/P theory D   measured D | sum theory  sum measured-D | P(meas D)/P(theory D)")
    for row in C:
        tc = row[0]
        if tc > t1: continue
        x = [terms(R, tc - 75.0), terms(R, tc + 75.0)]
        f = {kk: np.mean([y[kk] for y in x]) for kk in x[0]}
        Pt = f["Pth"]
        print(f"        {tc:6.0f}    {row[3]/Pt:+.2f}       {-f['flu']/Pt:+.2f}    |      {-f['Wth']/Pt:+.2f}            {-f['Wme']/Pt:+.2f}     |   {(-f['flu']-f['Wth'])/Pt:+.2f}        {(-f['flu']-f['Wme'])/Pt:+.2f}       |   {f['Pme']/Pt:.2f}")

r"""Drain rate of the vortices from flute stress + kinetic damping (4 Oct 2026):
dE_nz/dt (pred) = -sum_m (flute Reynolds-stress energy rate into the jets) - sum_k < (w_k / wt_k) P_k >,
the second term the lab-frame energy change of vortex k under its Landau damping (absorbed power P at local wt).
Measured: slope of E_nz (flute, fluid units, k_y = 1, both vortices) in the window.  nx 64 and nx 128, hyp_x = 0."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; gxx, gyy, Gxx, Gyy = G["gxx"], G["gyy"], G["Gxx"], G["Gyy"]
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
L2, k = fg.LAM2, fg.KYMIN; nxf = 512
src = open(os.path.join(HERE, "drain_from_stress.py")).read(); exec(src[src.index("def tend"):src.index("def fdep")])
def window_terms(R, tc):
    t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
    wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
    psz = (p0.mean(0) * w[None, :]).sum(1); U = R.to_x(psz, nxf, 1).real
    flu = 0; Wlab = 0; Pabs = 0; Enz = 0
    for kk in range(2):
        a0 = (a[kk] * w[None, :]).sum(1)
        flu += float(np.sum(np.real(np.conj(psz) * tend(R, a0, Gxx, Gyy))))
        ph0 = R.to_x(a0, nxf); dph0 = R.to_x(a0, nxf, 1); wt = wf[kk] - k * U
        Pd = 2 * (L2 * k ** 2) * k ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH)
        Wlab += np.mean(wf[kk] / np.where(np.abs(wt) < 1e-6, 1e-6, wt) * Pd); Pabs += Pd.mean()
        Enz += np.mean(Gxx * np.abs(dph0) ** 2 + Gyy * k ** 2 * np.abs(ph0) ** 2)
    return dict(flu=flu, Wlab=Wlab, P=Pabs, Enz=Enz)
CASES = [("hxoff nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64),
         ("nx128_hxoff", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)]
for name, dirs, nx in CASES:
    R = fg.Run(dirs, nx)
    print(f"\n{name}:   window      measured dE_nz/dt   pred = -flute - Wlab      [-flute    -Wlab ]    absorbed P   E_nz    Gamma meas / pred")
    for (t0, t1) in ((2700.0, 3300.0), (3300.0, 3900.0), (3900.0, 4400.0)):
        J = np.where((R.t >= t0) & (R.t <= t1))[0][::5]
        En = []
        for j in J:
            tt, phi = R.frame(j); pk = phi[:, 1, :] / fg.CXY; a0 = (pk * w[None, :]).sum(1)
            En.append((tt, np.mean(Gxx * np.abs(R.to_x(a0, nxf, 1)) ** 2 + Gyy * k ** 2 * np.abs(R.to_x(a0, nxf)) ** 2)))
        En = np.array(En); meas = np.polyfit(En[:, 0], En[:, 1], 1)[0]
        tw = [window_terms(R, tc) for tc in np.linspace(t0 + 60, t1 - 60, 4)]
        fl = np.mean([x["flu"] for x in tw]); Wl = np.mean([x["Wlab"] for x in tw]); Pa = np.mean([x["P"] for x in tw]); E = np.mean([x["Enz"] for x in tw])
        pred = -fl - Wl
        print(f"            {t0:.0f}-{t1:.0f}    {meas:+.3e}          {pred:+.3e}          [{-fl:+.3e} {-Wl:+.3e}]   {Pa:.3e}  {E:.3e}   {-meas/E:.2e} / {-pred/E:.2e}")

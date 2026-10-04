"""Fraction of the absorbed power P (loss formula, theory D) where wt and w have opposite signs (jet cores, outrun by the
vortex), and the lab-frame damping term -<(w/wt)P>/<P> (positive = the damping gives the vortices energy). hyp_x = 0 runs."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
L2, k = fg.LAM2, fg.KYMIN; nxf = 512
for name, dirs, nx in (("hxoff nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64),
                       ("nx128_hxoff", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    R = fg.Run(dirs, nx)
    for tc in (2000.0, 3000.0, 4200.0):
        if tc > R.t[-1] - 40 or tc < R.t[0] + 40: continue
        t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
        wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        U = R.to_x((p0.mean(0) * w[None, :]).sum(1), nxf, 1).real
        Pc = Pt = Wl = 0.0
        for kk in range(2):
            ps = R.to_x((a[kk] * w[None, :]).sum(1), nxf); wt = wf[kk] - k * U
            Pd = np.abs(ps) ** 2 * np.interp(np.abs(wt), WG, DTH)
            Pt += Pd.sum(); Pc += Pd[np.sign(wt) != np.sign(wf[kk])].sum(); Wl += (wf[kk] / wt * Pd).sum()
        print(f"{name} t {tc:.0f}: fraction of P in the jet cores {Pc/Pt:.2f};  lab-frame damping term -<(w/wt)P>/<P> = {-Wl/Pt:+.2f}")

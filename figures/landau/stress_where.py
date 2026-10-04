"""Where does the vortices' flute Reynolds stress on the jets' m = 1 come from, and where does it differ between the
reference and doubled radial resolution?  Mask the vortex stream functions within +-d of their critical layers
(smooth mask, d in Debye lengths) and recompute the m = 1 energy input (fluid units), t 3000-4000."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; Gxx, Gyy = G["Gxx"], G["Gyy"]; k = fg.KYMIN; LD = np.sqrt(fg.LAM2)
src = open(os.path.join(HERE, "drain_from_stress.py")).read(); exec(src[src.index("def tend"):src.index("def fdep")])
for case, dirs, NX in (("nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64), ("nx128", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    R = fg.Run(dirs, NX); nxf = 1024; X = np.arange(nxf) * fg.LX / nxf
    mm = np.abs(np.round(R.kx / (2 * np.pi / fg.LX))).astype(int)
    res = {}
    for tc in np.linspace(3060, 3940, 5):
        t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
        wf, a, _ = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        psz = (p0.mean(0) * w[None, :]).sum(1); U = R.to_x(psz, nxf, 1).real
        for d in (0.0, 1.0, 2.0, 4.0):
            tot = np.zeros(NX, complex)
            for kk in range(2):
                a0 = (a[kk] * w[None, :]).sum(1); ps = R.to_x(a0, nxf); c = wf[kk] / k
                s = np.sign(U - c); idx = np.where(s != np.roll(s, -1))[0]
                mask = np.ones(nxf)
                for i in idx:
                    dist = np.abs((X - X[i] + fg.LX / 2) % fg.LX - fg.LX / 2) / LD
                    mask *= 1 - np.exp(-(dist / max(d, 1e-9)) ** 8) if d > 0 else 1
                cf = np.fft.fft(ps * mask) / nxf; cm = np.zeros(NX, complex); h = NX // 2; cm[:h] = cf[:h]; cm[NX - h + 1:] = cf[nxf - h + 1:]
                tot += tend(R, cm, Gxx, Gyy)
            res.setdefault(d, []).append(np.real(np.conj(psz) * tot)[mm == 1].sum())
    print(case, "  m = 1 flute-stress energy input, critical layers masked within +-d lambda_D:  " + "   ".join(f"d={d:.0f}: {np.mean(v):+.3f}" for d, v in res.items()))

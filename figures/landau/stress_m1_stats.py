"""m = 1 jet forcing (fluid units): flute Reynolds stress of the two vortices and their deposition f = kP/wt, over many
windows (every 60 time units, t 2800-4300), at both radial resolutions, with standard errors - is the resolution
difference of the stress larger than its window-to-window scatter?"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; Gxx, Gyy = G["Gxx"], G["Gyy"]; k = fg.KYMIN; L2 = fg.LAM2
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
src = open(os.path.join(HERE, "drain_from_stress.py")).read(); exec(src[src.index("def tend"):src.index("def fdep")])
for case, dirs, NX in (("nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64), ("nx128", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    R = fg.Run(dirs, NX); nxf = 512; mm = np.abs(np.round(R.kx / (2 * np.pi / fg.LX))).astype(int)
    ST, DP = [], []
    for tc in np.arange(2830, 4290, 60.0):
        t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
        wf, a, _ = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        psz = (p0.mean(0) * w[None, :]).sum(1); U = R.to_x(psz, nxf, 1).real
        flu = np.zeros(NX, complex); f = np.zeros(nxf)
        for kk in range(2):
            a0 = (a[kk] * w[None, :]).sum(1); flu += tend(R, a0, Gxx, Gyy)
            ph0 = R.to_x(a0, nxf); wt = wf[kk] - k * U
            Pd = 2 * (L2 * k ** 2) * k ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH)
            f += k * Pd / np.where(np.abs(wt) < 1e-6, 1e-6, wt)
        c = np.fft.fft(f) / nxf; fc = np.zeros(NX, complex); h = NX // 2; fc[:h] = c[:h]; fc[NX - h + 1:] = c[nxf - h + 1:]
        dep = -1j * R.kx * fc
        ST.append(np.real(np.conj(psz) * flu)[mm == 1].sum()); DP.append(np.real(np.conj(psz) * dep)[mm == 1].sum())
    ST, DP = np.array(ST), np.array(DP); n = len(ST)
    print(f"{case}: {n} windows   stress {ST.mean():+.3f} +- {ST.std()/np.sqrt(n):.3f} (scatter {ST.std():.3f})   deposition {DP.mean():+.3f} +- {DP.std()/np.sqrt(n):.3f}   sum {np.mean(ST+DP):+.3f} +- {np.std(ST+DP)/np.sqrt(n):.3f}")

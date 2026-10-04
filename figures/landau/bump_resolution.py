"""Critical-layer non-flute bump of the vortices against the vortex vorticity there, at nx 64 and nx 128 (4 Oct 2026).
ratio = (projection of the non-flute part on g^yy - <g^yy>)/phi_0  divided by  lambda_D^2 zeta_1/phi_0, at |wt| < 0.1."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; sh = G["gyy"] - G["Gyy"]; sh -= (w * sh).sum()
RB = fg.RB; D = RB + "/rb_c0p3_draintest/"; S = RB + "/rb_c0p3_hxy_sinktest/"
CASES = [("hxoff nx64", [RB + "/rb_c0p3_hxoff/leg_0001/out"], 64, (2800.0, 3300.0, 3800.0, 4300.0)),
         ("nx128_hxoff", [D + "nx128_hxoff/out", D + "nx128_hxoff_leg2/out"], 128, (2800.0, 3300.0, 3800.0, 4300.0)),
         ("ref nx64", [RB + f"/rb_c0p3_hxy/leg_000{l}/out" for l in (4, 5)], 64, (9600.0, 10000.0)),
         ("nx128 (sinktest)", [S + "nx128/out"], 128, (9600.0, 10000.0)),
         ("nx128_hx16 (control)", [S + "nx128_hx16/out"], 128, (9600.0, 10000.0))]
NXF = 512
for name, dirs, nx, tcs in CASES:
    R = fg.Run(dirs, nx)
    for tc in tcs:
        if not ((R.t > tc - 30).any() and (R.t < tc + 30).any()): continue
        t, P = R.window(tc)
        p0, p1 = P[:, :, 0, :], P[:, :, 1, :]
        w0 = fg.freq_guess(t, p1[:, :, fg.ZMID])
        wf, a, res = fg.two_freq(t, p1, w0)
        U = R.to_x((p0.mean(0) * w[None, :]).sum(1), NXF, 1).real / fg.CXY
        out = []
        for k in range(2):
            ph = R.to_x(a[k], NXF); ph0 = (ph * w[None, :]).sum(1); nf = ph - ph0[:, None]
            A = (nf * w[None, :] * sh[None, :]).sum(1) / (w * sh ** 2).sum()
            a0 = (a[k] * w[None, :]).sum(1)
            zeta = fg.LAM2 * (G["Gxx"] * R.to_x(a0, NXF, 2) - G["Gyy"] * fg.KYMIN ** 2 * R.to_x(a0, NXF))
            wt = wf[k] - fg.KYMIN * U
            s = (np.abs(ph0) > 0.3 * np.abs(ph0).max()) & (np.abs(wt) < 0.1)
            bump = np.mean((A / ph0)[s].real); zr = np.mean((zeta / ph0)[s].real)
            mid = (np.abs(ph0) > 0.3 * np.abs(ph0).max()) & (np.abs(wt) > 1.0)
            out.append((bump, zr, bump / zr, np.mean((A / ph0)[mid])))
        b = np.mean([o[0] for o in out]); z = np.mean([o[1] for o in out]); r = np.mean([o[2] for o in out]); m = np.mean([o[3] for o in out])
        print(f"{name:22s} t {tc:6.0f}: w {wf[0]:+.3f} (fit resid {res:.3f});  bump {b:+.3f}  lam^2 zeta/phi0 {z:+.4f}  ratio {r:5.2f}   mid-jet {m.real:+.3f}{m.imag:+.3f}i")

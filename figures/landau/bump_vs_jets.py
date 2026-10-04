"""Critical-layer bump ratio of the vortices (bump / lambda^2 zeta_1, static theory 4.85) and the jets' m = 1 non-flute
factor, on one time grid through the reference run (rb_c0p3_hxy legs 1-6): do they track each other?"""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import fieldgen as fg, bounce, response
G = fg.geometry(); w = G["w"]; sh = G["gyy"] - G["Gyy"]; sh -= (w * sh).sum(); NXF = 512
geo = bounce.geometry(f"{fg.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Rs = response.Response(geo); Rs.source = lambda k2: k2 * (geo["gxx"] - (wz * geo["gxx"]).sum()); chi = np.conj(Rs.solve(0.02, 1.0))
kx1 = 2 * np.pi / fg.LX
R = fg.Run(sorted(glob.glob(fg.RB + "/rb_c0p3_hxy/leg_000[1-6]/out")), 64)
print("   t     bump ratio (static 4.85)   /4.85    jets m=1 factor")
for tc in (1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 9800):
    if not ((R.t > tc - 30).any() and (R.t < tc + 30).any()): continue
    t, P = R.window(float(tc)); p0, p1 = P[:, :, 0, :], P[:, :, 1, :]
    wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
    U = R.to_x((p0.mean(0) * w[None, :]).sum(1), NXF, 1).real / fg.CXY
    rs = []
    for k in range(2):
        ph = R.to_x(a[k], NXF); ph0 = (ph * w[None, :]).sum(1); nf = ph - ph0[:, None]
        A = (nf * w[None, :] * sh[None, :]).sum(1) / (w * sh ** 2).sum(); a0 = (a[k] * w[None, :]).sum(1)
        zeta = fg.LAM2 * (G["Gxx"] * R.to_x(a0, NXF, 2) - G["Gyy"] * fg.KYMIN ** 2 * R.to_x(a0, NXF))
        wt = wf[k] - fg.KYMIN * U; s = (np.abs(ph0) > 0.3 * np.abs(ph0).max()) & (np.abs(wt) < 0.1)
        rs.append(np.mean((A / ph0)[s].real) / np.mean((zeta / ph0)[s].real))
    pz = p0.mean(0); pz0 = (pz * wz[None, :]).sum(1); pz1 = pz[1] - pz0[1]
    pred = fg.LAM2 * kx1 ** 2 * chi * pz0[1]; c = np.vdot(pred * wz, pz1) / np.vdot(pred * wz, pred)
    r = np.mean(rs); print(f"{tc:6d}       {r:5.2f}                 {r/4.85:5.2f}      {abs(c):5.3f}")

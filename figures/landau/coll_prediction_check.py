"""Pre-registered test (written 4 Oct 2026, 21:05 [corrected from 21:50], before coll_hxoff ran): if the jets' extra non-flute factor (~1.6 x the
static response, STATUS (w)) is the pitch-angle structure of their charge, pitch-angle scattering (Landau collisions,
nu = 0.01, 1/nu = 100) should relax it to the static value 1 within a few hundred time units of the restart at t 2500,
while the collisionless control (rb_c0p3_hxoff) stays at ~1.6.  Also prints the critical-layer bump ratio / 4.85."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import fieldgen as fg, bounce, response
G = fg.geometry(); w = G["w"]; sh = G["gyy"] - G["Gyy"]; sh -= (w * sh).sum(); NXF = 512
geo = bounce.geometry(f"{fg.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Rs = response.Response(geo); Rs.source = lambda k2: k2 * (geo["gxx"] - (wz * geo["gxx"]).sum()); chi = np.conj(Rs.solve(0.02, 1.0))
kx1 = 2 * np.pi / fg.LX
D = fg.RB + "/rb_c0p3_draintest/"
for name, dirs in (("coll_hxoff", [D + "coll_hxoff/out", D + "coll_hxoff_leg2/out"]), ("hxoff (control)", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"])):
    dirs = [d for d in dirs if os.path.exists(d + "/field.dat")]
    if not dirs: print(name, ": no output yet"); continue
    R = fg.Run(dirs, 64)
    print(f"\n{name}:   t     jets m=1 factor    bump/4.85")
    for tc in (2600, 2700, 2800, 3000, 3300, 3600, 4000, 5000, 6000, 7000):
        if not ((R.t > tc - 30).any() and (R.t < tc + 30).any()): continue
        t, P = R.window(float(tc)); p0, p1 = P[:, :, 0, :], P[:, :, 1, :]
        pz = p0.mean(0); pz0 = (pz * wz[None, :]).sum(1); pz1 = pz[1] - pz0[1]
        pred = fg.LAM2 * kx1 ** 2 * chi * pz0[1]; c = np.vdot(pred * wz, pz1) / np.vdot(pred * wz, pred)
        wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        U = R.to_x((p0.mean(0) * w[None, :]).sum(1), NXF, 1).real / fg.CXY; rs = []
        for k in range(2):
            ph = R.to_x(a[k], NXF); ph0 = (ph * w[None, :]).sum(1); nf = ph - ph0[:, None]
            A = (nf * w[None, :] * sh[None, :]).sum(1) / (w * sh ** 2).sum(); a0 = (a[k] * w[None, :]).sum(1)
            zeta = fg.LAM2 * (G["Gxx"] * R.to_x(a0, NXF, 2) - G["Gyy"] * fg.KYMIN ** 2 * R.to_x(a0, NXF))
            wt = wf[k] - fg.KYMIN * U; s = (np.abs(ph0) > 0.3 * np.abs(ph0).max()) & (np.abs(wt) < 0.1)
            rs.append(np.mean((A / ph0)[s].real) / np.mean((zeta / ph0)[s].real))
        print(f"          {tc:6d}     {abs(c):5.3f}            {np.mean(rs)/4.85:5.2f}")

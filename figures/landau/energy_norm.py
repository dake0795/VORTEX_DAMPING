"""Zonal electrostatic energy per radial harmonic: GENE's Spectral_2D_E_fe (k_y = 0, species summed) against the
field-based fluid energy 1/2 <g^xx> k_x^2 |psi|^2 (flute average, psi = phi/C_xy), and against field-based candidates
with the local metric and with GENE's own form lambda^2 k_perp^2 |phi|^2 integrated with the Jacobian."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import fieldgen as fg, spec2d
G = fg.geometry(); w = G["w"]; gxx, Gxx, J = G["gxx"], G["Gxx"], G["J"]
D = fg.RB + "/rb_c0p3_hxoff/leg_0001/out"
R = fg.Run([D], 64)
kx_s, ky_s, ts, tot, d = None, None, None, None, 0
for sp in ("electrons", "positrons"):
    a = spec2d.read(f"{D}/Spectral_2D_E_fe_{sp}.dat"); kx_s, ky_s, ts = a[0], a[1], a[2]; d = d + a[4]
for tc in (2000.0, 4000.0):
    i = np.argmin(abs(ts - tc)); j = np.argmin(abs(R.t - tc))
    tt, phi = R.frame(j)
    print(f"\nt GENE {ts[i]:.1f}, field {tt:.1f}")
    print("  m   GENE E_fe(ky=0)   fluid 1/2<gxx>kx^2|psi0|^2   ratio   local lam^2 int J gxx kx^2|phi|^2 / int J   ratio")
    for m in (1, 2, 3, 4, 5, 7, 9):
        kx = m * 2 * np.pi / fg.LX
        ig = [np.argmin(abs(kx_s - kx)), np.argmin(abs(kx_s + kx))]
        Eg = d[i, ig, 0].sum()
        ph = phi[[m, 64 - m], 0, :]                                          # +-m, all z
        ps0 = (ph * w[None, :]).sum(1) / fg.CXY
        Ef = 0.5 * np.sum(Gxx * kx ** 2 * np.abs(ps0) ** 2)
        El = fg.LAM2 * np.sum((J[None, :] * gxx[None, :] * kx ** 2 * np.abs(ph) ** 2).sum(1)) / J.sum()
        print(f"  {m:2d}   {Eg:.4e}          {Ef:.4e}                {Eg/Ef:8.3f}   {El:.4e}   {Eg/El:8.4f}")

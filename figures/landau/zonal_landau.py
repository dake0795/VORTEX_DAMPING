r"""Landau damping of the jets' radial harmonics that oscillate at the breathing frequency (2 Oct 2026).

A zonal harmonic of radial wavenumber k_x has a non-flute part forced by the variation of g^xx along the field line,
S(l) = k_x^2 lambda_D^2 <g^xx> [g^xx(l)/<g^xx> - 1] (the k_y = 0 counterpart of the vortices' source); oscillating at
frequency w (2 w_vortex, the breathing) it is damped by the particles that stream along the line, at the rate
nu_z = 2 k_x^2 lambda_D^2 D_x(w) / <g^xx>,  D_x(w) = -w <s_x Im chi>   (the loss formula of the notes with k_y -> k_x).
Test: GENE's parallel-streaming term of the k_y = 0 electrostatic energy, per |m|, against nu_z(m) times the energy in
the oscillating part of that harmonic (from the flute potential: time variance about the slow mean).
"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat")
R = response.Response(geo)
wz = geo["J"] / geo["J"].sum(); Gxx = (wz * geo["gxx"]).sum()
sx = Gxx * (geo["gxx"] / Gxx - 1.0)
R.source = lambda k2lam2: k2lam2 * sx
def Dx(w):
    return -w * (wz * np.imag(R.solve(w, 1.0)) * sx).sum()
L2 = 5000.0; KXMIN = 2 * np.pi / vf.LX
# oscillating fraction of each zonal harmonic in GENE (flute potential, Jacobian-weighted average)
G = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "gene_hxoff.npz")); t = G["t"]; z = G["phi_avg"][:, :, 0]
for (t0, t1, w) in ((1000, 2000, 2 * 1.144), (3000, 4400, 2 * 0.915)):
    s = (t >= t0) & (t <= t1); zt = z[s]
    # slow mean = running mean over 20 t.u. (several breathing periods); oscillating part = the rest
    nwin = max(3, int(round(20.0 / np.median(np.diff(t[s])))))
    k = np.ones(nwin) / nwin
    slow = np.array([np.convolve(zt[:, i], k, mode="same") for i in range(zt.shape[1])]).T
    osc = zt - slow
    frac = (np.abs(osc[nwin:-nwin]) ** 2).mean(0) / (np.abs(zt[nwin:-nwin]) ** 2).mean(0)
    D = Dx(w)
    print(f"\nwindow {t0}-{t1}: breathing frequency 2w = {w:.3f}; D_x(2w) = {D:.4e}")
    print("   m  kx lamD   osc. fraction   nu_z predicted")
    out = []
    for m in range(1, 16):
        kx = m * KXMIN; nu = 2 * kx ** 2 * L2 * D / Gxx
        f = 0.5 * (frac[m] + frac[-m])
        out.append((m, kx * np.sqrt(L2), f, nu)); print(f"  {m:2d}   {kx*np.sqrt(L2):5.2f}     {f:.3e}       {nu:.3e}")
    np.save(os.path.join(HERE, f"zonal_landau_{t0}.npy"), np.array(out))

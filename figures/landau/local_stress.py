r"""Zonal forcing of the jets by the two vortices computed LOCALLY along the field line and then field-line averaged,
against the same forcing computed from the field-line-averaged potential (the flute/Euler estimate).  At each z:
charge rho_z = -lambda_D^2 (g^xx(z) d_x^2 + g^yy(z) d_y^2) phi_z (polarisation, local metric), advected by the local
E x B flow of phi_z (C_xy factor as in the notes).  Each vortex = travelling pattern from the two-frequency fit at every
z (harmonics 1..NH of k_y, frequencies n w); breathing-averaged stress = sum of the two self terms.  Compared with
GENE's exact budget per |m| (E_nonlinear + E_parallel at k_y = 0)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat")
J = geo["J"]; wz = J / J.sum(); gxx, gyy = geo["gxx"], geo["gyy"]
NX, NKY, NZ = vf.NX, vf.NKY, vf.NZ
kx = np.fft.fftfreq(NX, 1.0 / NX) * 2 * np.pi / vf.LX; ky = np.arange(NKY) * vf.KYMIN
nxg, nyg = 3 * NX // 2, 3 * NKY
mg = np.fft.fftfreq(nxg, 1.0 / nxg); KXg = (mg * 2 * np.pi / vf.LX)[:, None]; KYg = (np.arange(nyg // 2 + 1) * vf.KYMIN)[None, :]
def pad(p):
    a = np.zeros((nxg, nyg // 2 + 1), complex); h = NX // 2
    a[:h, :NKY] = p[:h]; a[nxg - (h - 1):, :NKY] = p[NX - (h - 1):]
    return a * nxg * nyg
def zonal_tend(p, Gxx, Gyy):
    """k_y = 0 tendency of rho = (Gxx k_x^2 + Gyy k_y^2) psi under advection by psi (psi = phi/C_xy) -> GENE kx order"""
    P = pad(p); Z = -(Gxx * KXg ** 2 + Gyy * KYg ** 2) * P
    f = lambda A, k: np.fft.irfft2(1j * k * A, s=(nxg, nyg))
    Nn = -np.fft.rfft2(f(P, KXg) * f(Z, KYg) - f(P, KYg) * f(Z, KXg))[:, 0] / (nxg * nyg)
    h = NX // 2; out = np.zeros(NX, complex); out[:h] = Nn[:h]; out[NX - (h - 1):] = Nn[nxg - (h - 1):]
    return out
NH = 4
TCS = [float(x) for x in sys.argv[1:]] or [1500.0, 2100.0, 2700.0, 3300.0, 3900.0]
for tc in TCS:
    T, P0, P1 = predict.window_field("hxoff", tc)          # only k_y = 0, 1; read harmonics separately below
    # all k_y: re-read frames
    frames = []
    ts = vf.frame_times(vf.RUNS["hxoff"][0] + "/field.dat")
    idx = np.where((ts >= tc - 30) & (ts <= tc + 30))[0]
    for i in idx:
        t, phi = vf.read_frame(vf.RUNS["hxoff"][0] + "/field.dat", i); frames.append(phi)
    F = np.array(frames) / vf.CXY; t = ts[idx]; tt = t - t.mean()
    Rw = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
    w0 = Rw[np.argmin(abs(Rw[:, 0] - tc)), 1]
    w, _, _ = vf.two_freq_fit(t, F[:, :, 1, vf.ZMID], (w0, -w0))
    loc = np.zeros(NX, complex); avg_p = [np.zeros((NX, NKY), complex) for _ in range(2)]; allp = []
    for iz in range(NZ):
        pats = [np.zeros((NX, NKY), complex) for _ in range(2)]
        for n in range(1, NH + 1):
            B = np.stack([np.exp(-1j * n * w[0] * tt), np.exp(-1j * n * w[1] * tt)], 1)
            sol, *_ = np.linalg.lstsq(B, F[:, :, n, iz], rcond=None)
            pats[0][:, n], pats[1][:, n] = sol[0], sol[1]
        for k in range(2):
            loc += wz[iz] * zonal_tend(pats[k], gxx[iz], gyy[iz])
            avg_p[k] += wz[iz] * pats[k]
        allp.append(pats)
    Gxx, Gyy = (wz * gxx).sum(), (wz * gyy).sum()
    flute = zonal_tend(avg_p[0], Gxx, Gyy) + zonal_tend(avg_p[1], Gxx, Gyy)
    metric_only = sum(wz[iz] * (zonal_tend(avg_p[0], gxx[iz], gyy[iz]) + zonal_tend(avg_p[1], gxx[iz], gyy[iz])) for iz in range(NZ))
    print(f"\nt {tc:.0f}: w {w[0]:+.3f} {w[1]:+.3f};   |m|   local-then-averaged    flute (averaged potential)   ratio")
    for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13):
        a = abs(loc[m]) + abs(loc[-m]); b = abs(flute[m]) + abs(flute[-m])
        c = abs(metric_only[m]) + abs(metric_only[-m])
        print(f"        {m:2d}      {a:.3e}              {b:.3e}             {a/b:.3f}   phase {np.degrees(np.angle(loc[m]/flute[m])):+5.0f}   metric-only/flute {c/b:.3f}")
    np.savez(os.path.join(HERE, f"local_stress_{int(tc)}.npz"), loc=loc, flute=flute, metric_only=metric_only)

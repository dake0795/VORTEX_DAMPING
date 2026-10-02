r"""Non-flute part of a vortex on jets that are themselves not flute-like (2-3 Oct 2026).
phi_1 - rho[phi_1] = -q_1,  q_1 = non-flute part of the vortex's polarisation charge generated following the flow:
  (i)  the appendix term (metric variation acting on the flute vortex, jets' flute vorticity included): -k^2 lam^2 s(l) phi_0;
  (ii) NEW: the jets' own non-flute polarisation charge Q_z1 = -lam^2 [g^xx(l) d_x^2 phi_z1]_nf advected radially by the
       vortex:  q_1 = -(k / wt) psi_0 d_x Q_z1   (psi_0 = phi_0 / C_xy, wt = w - k U).
Both solved with the kinetic response (response.Response) at the local wt, x by x; compared with the measured
non-flute part (projection on g^yy - <g^yy>, relative to phi_0) of the fitted vortex in the reference run (t 6100)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response, predict
RUN = sys.argv[1] if len(sys.argv) > 1 else "orig"; TC = float(sys.argv[2]) if len(sys.argv) > 2 else 6100.0
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
s_app = Gyy * (geo["gxx"] / Gxx - geo["gyy"] / Gyy)
sh = geo["gyy"] - Gyy; sh -= (wz * sh).sum()
Rs = response.Response(geo); L2 = 5000.0; k = vf.KYMIN; NXF = predict.NXF; xf = np.arange(NXF) * vf.LX / NXF
Tr = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)["orig"]
t, p0, p1 = predict.window_field(RUN, TC)
i = np.argmin(abs(Tr[:, 0] - TC)); w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], (Tr[i, 1], Tr[i, 2]))
Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, vf.NX, vf.NZ)
pz = p0.mean(axis=0)                                            # zonal potential [nx, nz]
U = vf.to_x((pz * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
pz0 = (pz * wz[None, :]).sum(axis=1); pz1 = pz - pz0[:, None]  # jets' non-flute part (GENE coefficients)
d2z1 = predict.to_x(pz1, deriv=2).real                          # d_x^2 phi_z1 [NXF, nz]
Qz1 = -L2 * geo["gxx"][None, :] * d2z1; Qz1 -= (Qz1 * wz[None, :]).sum(1)[:, None]
kf = np.fft.fftfreq(NXF, 1.0 / NXF) * 2 * np.pi / vf.LX; dQ = np.fft.ifft(1j * kf[:, None] * np.fft.fft(Qz1, axis=0), axis=0).real
print(f"jets: |phi_z1|/|phi_z0| = {np.sqrt((np.abs(predict.to_x(pz1))**2*wz).sum()/ (np.abs(predict.to_x(pz0[:,None]))**2).mean()):.3e}")
def solve(src, wt):
    Rs.source = lambda k2, s_=src: k2 * s_
    return Rs.solve(wt, 1.0)
SIG = float(os.environ.get("SIG", "-1"))          # convention: v_x = SIG * i k psi (GENE's Fourier sign); response conjugated if SIG = -1
U1 = predict.to_x(pz1, deriv=1).real / vf.CXY      # jets' non-flute velocity [NXF, nz]
def solve_c(src, wt):
    r = solve(src, wt); return np.conj(r) if SIG < 0 else r
for kk in range(2):
    ph = predict.to_x(a[kk]); ph0 = (ph * wz[None, :]).sum(1); meas = ph - ph0[:, None]
    d2v = vf.to_x((a[kk] * wz[None, :]).sum(axis=1), nxf=NXF, deriv=2)
    qv0 = -L2 * (Gxx * d2v - Gyy * k ** 2 * ph0)
    wt = w[kk] - k * U
    print(f"\nvortex {kk}, t {TC:.0f} ({RUN}), SIG {SIG:+.0f}:   x      wt     measured A/phi0      (i) appendix     (ii) jets' nf charge   (iii) nf jet flow     total")
    rows = []
    for j in range(0, NXF, 6):
        if abs(ph0[j]) < 0.3 * abs(ph0).max() or abs(wt[j]) < 0.02: continue
        Am = (meas[j] * wz * sh).sum() / (wz * sh ** 2).sum() / ph0[j]
        p_i = -k ** 2 * L2 * solve_c(s_app, wt[j])
        p_ii = solve_c(SIG * (k / wt[j]) * (1.0 / vf.CXY) * dQ[j], wt[j])
        q3 = -SIG * (k * U1[j] / wt[j]) * qv0[j] / ph0[j]; q3 = q3 - (wz * q3).sum()
        p_iii = solve_c(-q3, wt[j])
        proj = lambda p: (p * wz * sh).sum() / (wz * sh ** 2).sum()
        A1, A2, A3 = proj(p_i), proj(p_ii), proj(p_iii); At = A1 + A2 + A3
        rows.append((xf[j], wt[j], Am, At))
        print(f"            {xf[j]:7.1f} {wt[j]:+.3f}   {Am.real:+.3f}{Am.imag:+.3f}i    {A1.real:+.3f}{A1.imag:+.3f}i   {A2.real:+.3f}{A2.imag:+.3f}i      {A3.real:+.3f}{A3.imag:+.3f}i    {At.real:+.3f}{At.imag:+.3f}i")
    R = np.array(rows, dtype=object); Am = np.array([r[2] for r in rows]); At = np.array([r[3] for r in rows])
    print(f"   overall: |measured - theory| / |measured| = {np.linalg.norm(Am - At) / np.linalg.norm(Am):.3f};  best scale {np.vdot(At, Am) / np.vdot(At, At):.3f}")

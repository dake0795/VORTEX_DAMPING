r"""Local stress of the vortices on the jets with the vortices' non-flute part PREDICTED (terms i-iii of nonflute_jets.py)
instead of measured, against the stress with the measured non-flute part and with none (flute), per odd jet harmonic,
in rb_c0p3_hxoff (the run of compare_budget.py)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
os.environ["SIG"] = "-1"
import vortex_field as vf, bounce, response, predict
NXF_ = int(os.environ.get('NXF', 64)); predict.NXF = NXF_
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
gxx, gyy = geo["gxx"], geo["gyy"]; Gxx, Gyy = (wz * gxx).sum(), (wz * gyy).sum()
s_app = Gyy * (gxx / Gxx - gyy / Gyy)
Rs = response.Response(geo); L2 = 5000.0; k = vf.KYMIN; NX, NKY, NZ = vf.NX, vf.NKY, vf.NZ; NXF = int(os.environ.get("NXF", NX))
def solve(src, wt):
    Rs.source = lambda k2, s_=src: k2 * s_
    return np.conj(Rs.solve(wt, 1.0))
src = open(os.path.join(HERE, "local_stress.py")).read()
exec(src[src.index("kx = np.fft.fftfreq"):src.index("NH = 4")])
R2 = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
def to_coef(f):       # f[NXF, nz] real-space (vf.to_x convention) -> GENE coefficients [NX, nz]
    c = np.fft.fft(f, axis=0) / NXF; out = np.zeros((NX, f.shape[1]), complex); h = NX // 2
    out[:h] = c[:h]; out[NX - h + 1:] = c[NXF - h + 1:]; return out
for tc in [float(x) for x in sys.argv[1:]] or [1500.0, 2700.0, 3900.0]:
    t, p0, p1 = predict.window_field("hxoff", tc)
    w0 = R2[np.argmin(abs(R2[:, 0] - tc)), 1]; w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], (w0, -w0))
    Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, NX, NZ)
    pz = p0.mean(axis=0); pz0 = (pz * wz[None, :]).sum(1); pz1 = pz - pz0[:, None]
    U = vf.to_x(pz0, nxf=NXF, deriv=1).real / vf.CXY
    Qz1 = -L2 * gxx[None, :] * predict.to_x(pz1, deriv=2).real; Qz1 -= (Qz1 * wz[None, :]).sum(1)[:, None]
    kf = np.fft.fftfreq(NXF, 1.0 / NXF) * 2 * np.pi / vf.LX; dQ = np.fft.ifft(1j * kf[:, None] * np.fft.fft(Qz1, axis=0), axis=0).real; U1 = predict.to_x(pz1, deriv=1).real / vf.CXY
    loc_m = np.zeros(NX, complex); loc_p = np.zeros(NX, complex); flu = np.zeros(NX, complex)
    for kk in range(2):
        ph = predict.to_x(a[kk]); ph0 = (ph * wz[None, :]).sum(1)
        d2v = vf.to_x((a[kk] * wz[None, :]).sum(axis=1), nxf=NXF, deriv=2); qv0 = -L2 * (Gxx * d2v - Gyy * k ** 2 * ph0)
        wt = w[kk] - k * U
        dU = np.gradient(U, vf.LX / NXF); DEL = float(os.environ.get("DELFAC", "1")) * k * np.abs(dU) * (vf.LX / vf.NX)
        inv = wt / (wt ** 2 + DEL ** 2)                                   # 1/wt regularised over one GENE grid cell
        php = np.zeros((NXF, NZ), complex)
        for j in range(NXF):
            if abs(wt[j]) < 1e-3: continue
            q3 = (k * U1[j] * inv[j]) * qv0[j]; q3 = q3 - (wz * q3).sum()
            php[j] = (-k ** 2 * L2 * solve(s_app, wt[j]) * ph0[j] + solve(-(k * inv[j]) / vf.CXY * dQ[j] * ph0[j], wt[j]) + solve(-q3, wt[j]))
        pat_m = a[kk] / vf.CXY; pat_p = (to_coef(ph0[:, None] + php)) / vf.CXY; pat_f = np.repeat((a[kk] * wz[None, :]).sum(1)[:, None], NZ, 1) / vf.CXY
        for iz in range(NZ):
            for pat, acc in ((pat_m, loc_m), (pat_p, loc_p)):
                P = np.zeros((NX, NKY), complex); P[:, 1] = pat[:, iz]; acc += wz[iz] * zonal_tend(P, gxx[iz], gyy[iz])
        P = np.zeros((NX, NKY), complex); P[:, 1] = pat_f[:, 0]; flu += zonal_tend(P, Gxx, Gyy)
    sh = gyy - Gyy; sh = sh - (wz * sh).sum()
    pr = lambda f: (f * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum()
    kk = 0; ph = predict.to_x(a[kk]); ph0 = (ph * wz[None, :]).sum(1)
    Am = pr(ph - ph0[:, None]) / ph0; Ap = pr(php) / ph0          # php of the last vortex (kk = 1) -> recompute for kk = 0 not stored; use kk = 1
    ph1 = predict.to_x(a[1]); ph10 = (ph1 * wz[None, :]).sum(1); Am1 = pr(ph1 - ph10[:, None]) / ph10; Ap1 = pr(php) / ph10
    strong = np.abs(ph10) > 0.3 * np.abs(ph10).max()
    print(f"  t {tc:.0f}: non-flute structure of vortex 1, predicted vs measured: |diff|/|meas| = {np.linalg.norm((Ap1-Am1)[strong])/np.linalg.norm(Am1[strong]):.3f}, scale {np.vdot(Ap1[strong], Am1[strong])/np.vdot(Ap1[strong], Ap1[strong]):.3f}")
    if os.environ.get("SHOW"):
        Ur = U; wt1 = w[1] - k * Ur
        for j in range(NXF):
            if strong[j]: print(f"     x {j*vf.LX/NXF:7.1f} wt {wt1[j]:+.3f}  measured {Am1[j].real:+.3f}{Am1[j].imag:+.3f}i  predicted {Ap1[j].real:+.3f}{Ap1[j].imag:+.3f}i")
    print(f"\nt {tc:.0f}:   m    stress with measured non-flute      with PREDICTED non-flute      flute only      pred/meas")
    for m in (1, 3, 5, 7, 9, 11):
        A = loc_m[m] + 0; B = loc_p[m]; C = flu[m]
        print(f"        {m:2d}    {abs(A):.3e} ({np.degrees(np.angle(A)):+4.0f})          {abs(B):.3e} ({np.degrees(np.angle(B)):+4.0f})          {abs(C):.3e} ({np.degrees(np.angle(C)):+4.0f})   {abs(B/A):.2f} {np.degrees(np.angle(B/A)):+4.0f}")
    np.savez(os.path.join(HERE, f"stress_pred_{int(tc)}.npz"), loc_m=loc_m, loc_p=loc_p, flu=flu)

r"""Stress of the vortices on the jets with the non-flute part in rank-one form phi_1 = A(x) sh(z), sh = g^yy - <g^yy>:
(a) A measured (projection of GENE's non-flute part); (b) A predicted = C lam^2 zeta_1(x) + lam^2 k^2 chi_s-term (appH),
zeta_1 = (<g^xx> d_x^2 - <g^yy> k^2) phi_0 the vortex's flute vorticity; C = 4.85 (static response, predicted) or C_fit.
Compared with the stress from GENE's full non-flute part, per odd jet harmonic, rb_c0p3_hxoff."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, predict
predict.NXF = 64
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
gxx, gyy = geo["gxx"], geo["gyy"]; Gxx, Gyy = (wz * gxx).sum(), (wz * gyy).sum()
sh = gyy - Gyy; sh -= (wz * sh).sum()
L2 = 5000.0; k = vf.KYMIN; NX, NKY, NZ = vf.NX, vf.NKY, vf.NZ; NXF = 64
T = np.load(os.path.join(HERE, "chi_tables.npz")); WT = T["WT"]
pcs = np.array([np.conj((c * wz * sh).sum() / (wz * sh ** 2).sum()) for c in T["chis"]])
src = open(os.path.join(HERE, "local_stress.py")).read()
exec(src[src.index("kx = np.fft.fftfreq"):src.index("NH = 4")])
R2 = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
def to_coef(f):
    c = np.fft.fft(f, axis=0) / NXF; out = np.zeros((NX,) + f.shape[1:], complex); h = NX // 2
    out[:h] = c[:h]; out[NX - h + 1:] = c[NXF - h + 1:]; return out
CF = float(os.environ.get("CFIT", "7.0"))
for tc in [float(x) for x in sys.argv[1:]] or [1500.0, 2700.0, 3900.0]:
    t, p0, p1 = predict.window_field("hxoff", tc)
    w0 = R2[np.argmin(abs(R2[:, 0] - tc)), 1]; w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], (w0, -w0))
    Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, NX, NZ)
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    acc = {key: np.zeros(NX, complex) for key in ("full", "rank1", "predC", "predfit", "flute")}
    for kk in range(2):
        ph = predict.to_x(a[kk]); ph0 = (ph * wz[None, :]).sum(1); m = ph - ph0[:, None]
        Am = (m * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum()
        a0 = (a[kk] * wz[None, :]).sum(1); d2 = vf.to_x(a0, nxf=NXF, deriv=2)
        zeta = L2 * (Gxx * d2 - Gyy * k ** 2 * ph0)
        wt = w[kk] - k * U
        app = -L2 * k ** 2 * np.interp(wt, WT, pcs.real) * ph0 - 1j * L2 * k ** 2 * np.interp(wt, WT, pcs.imag) * ph0
        if w[kk] < 0: app = np.conj(app) if False else app
        crit = np.exp(-(wt / 0.25) ** 2)                                         # static form where |wt| small
        Ap = lambda C: C * zeta * crit + app * (1 - crit)
        pats = {"full": ph, "rank1": ph0[:, None] + Am[:, None] * sh[None, :], "predC": ph0[:, None] + Ap(4.85)[:, None] * sh[None, :],
                "predfit": ph0[:, None] + Ap(CF)[:, None] * sh[None, :], "flute": np.repeat(ph0[:, None], NZ, 1)}
        for key, pat in pats.items():
            pc = to_coef(pat) / vf.CXY
            for iz in range(NZ):
                P = np.zeros((NX, NKY), complex); P[:, 1] = pc[:, iz]; acc[key] += wz[iz] * zonal_tend(P, gxx[iz], gyy[iz])
        if kk == 0:
            s = np.abs(ph0) > 0.3 * np.abs(ph0).max()
            print(f"\nt {tc:.0f}: structure of vortex 0, rank-one predicted vs measured A: residual C=4.85 {np.linalg.norm((Ap(4.85)-Am)[s])/np.linalg.norm(Am[s]):.2f}, C={CF} {np.linalg.norm((Ap(CF)-Am)[s])/np.linalg.norm(Am[s]):.2f}")
    print("   m     full measured        rank-one measured       predicted C=4.85      predicted C=fit        flute")
    for mm in (1, 3, 5, 7, 9):
        row = "  ".join(f"{abs(acc[k_][mm]):.2e}({np.degrees(np.angle(acc[k_][mm])):+4.0f})" for k_ in ("full", "rank1", "predC", "predfit", "flute"))
        print(f"  {mm:2d}   {row}")

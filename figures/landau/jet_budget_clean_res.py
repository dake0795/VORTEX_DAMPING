r"""Jets' energy budget per radial harmonic, one normalisation (4 Oct 2026).  All rates in fluid energy units:
E_m = 1/2 * sum_{+-m} <g^xx> k_x^2 |psi_0|^2 (psi_0 flute zonal stream function); GENE's E_fe(k_y = 0) = 2 lambda^2 C^2 E
(checked to 0.05 %, energy_norm.py), so GENE's E_nonlinear and E_parallel at k_y = 0 are divided by the same factor
(computed per harmonic from E_fe / E).  Theory pieces from the two fitted vortices: flute Reynolds stress, local
(non-flute-resolved) stress, and the momentum deposited by their Landau damping, f = k P / wt."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import fieldgen as fg, spec2d
G = fg.geometry(); w = G["w"]; gxx, gyy, Gxx, Gyy = G["gxx"], G["gyy"], G["Gxx"], G["Gyy"]
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
L2, k = fg.LAM2, fg.KYMIN
src = open(os.path.join(HERE, "drain_from_stress.py")).read(); exec(src[src.index("def tend"):src.index("def fdep")])
CASE = sys.argv[1] if len(sys.argv) > 1 else "nx64"
DIRS = {"nx64": [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], "nx128": [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"]}[CASE]
NX = 64 if CASE == "nx64" else 128; R = fg.Run(DIRS, NX)
spec2d.NX = NX; spec2d.NKX = NX - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
S = {}
for term in ("fe", "nonlinear", "parallel"):
    accs, tss = [], []
    for D in DIRS:
        acc = 0
        for sp in ("electrons", "positrons"):
            a = spec2d.read(f"{D}/Spectral_2D_E_{term}_{sp}.dat"); acc = acc + a[4][:, :, 0]; kxs, tt = a[0], a[2]
        accs.append(acc); tss.append(tt)
    S[term] = np.concatenate(accs); ts = np.concatenate(tss)
mG = np.abs(np.round(kxs / (2 * np.pi / fg.LX))).astype(int)
mm = np.abs(np.round(R.kx / (2 * np.pi / fg.LX))).astype(int)
nxf = 512
for (t0, t1) in ((3000.0, 4000.0),):
    sel = (ts >= t0) & (ts <= t1)
    tcs = np.linspace(t0 + 60, t1 - 60, 6)
    acc = {key: np.zeros(NX) for key in ("flu", "loc", "dep", "E")}
    for tc in tcs:
        t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
        wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        psz = (p0.mean(0) * w[None, :]).sum(1); U = R.to_x(psz, nxf, 1).real
        loc = np.zeros(NX, complex); flu = np.zeros(NX, complex); f = np.zeros(nxf)
        for kk in range(2):
            for iz in range(fg.NZ):
                loc += w[iz] * tend(R, a[kk][:, iz], gxx[iz], gyy[iz])
            a0 = (a[kk] * w[None, :]).sum(1); flu += tend(R, a0, Gxx, Gyy)
            ph0 = R.to_x(a0, nxf); wt = wf[kk] - k * U
            Pd = 2 * (L2 * k ** 2) * k ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH)
            f += k * Pd / np.where(np.abs(wt) < 1e-6, 1e-6, wt)
        c = np.fft.fft(f) / nxf; fc = np.zeros(NX, complex); h = NX // 2; fc[:h] = c[:h]; fc[NX - h + 1:] = c[nxf - h + 1:]
        dep = -1j * R.kx * fc                       # d rho_0/dt, rho = -lap psi = K psi  (dU/dt = f, U = psi')
        acc["flu"] += np.real(np.conj(psz) * flu) / len(tcs); acc["loc"] += np.real(np.conj(psz) * loc) / len(tcs)
        acc["dep"] += np.real(np.conj(psz) * dep) / len(tcs); acc["E"] += 0.5 * Gxx * R.kx ** 2 * np.abs(psz) ** 2 / len(tcs)
    print(f"\nwindow {t0:.0f}-{t1:.0f} (rates in fluid energy units per unit time; per harmonic, +-m summed)")
    print("  m     E_m        GENE: nonlinear   parallel  |  theory: flute stress   local stress   deposition   | flute+dep   local+dep")
    for m in (1, 2, 3, 4, 5, 7, 9):
        q = mm == m; g = mG == m
        fac = S["fe"][sel][:, g].sum(1).mean() / acc["E"][q].sum()          # GENE units per fluid unit (2 lam^2 C^2 ~ 8500)
        nl = S["nonlinear"][sel][:, g].sum(1).mean() / fac; pa = S["parallel"][sel][:, g].sum(1).mean() / fac
        fl, lo, de = acc["flu"][q].sum(), acc["loc"][q].sum(), acc["dep"][q].sum()
        print(f"  {m:2d}  {acc['E'][q].sum():.3e}    {nl:+.3e}   {pa:+.3e}  |   {fl:+.3e}        {lo:+.3e}     {de:+.3e}   | {fl+de:+.3e}  {lo+de:+.3e}")

r"""Landau loss of the harmonics k_y = 2..4 of the vortices: loss formula (theory D, extended table) on each frequency
component j w_+ + (n - j) w_- of the k_y = n potential, against GENE's -E_parallel(k_y = n), relative to k_y = 1."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import fieldgen as fg, spec2d
G = fg.geometry(); w = G["w"]; Gxx, Gyy = G["Gxx"], G["Gyy"]
C = np.load(os.path.join(HERE, "cache_Dext.npz")); WG, DTH = C["wg"], C["Dth"]
L2, k = fg.LAM2, fg.KYMIN; nxf = 512; NMAX = 4
def gene_par(dirs, nx, a, b):
    spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    out = 0
    for sp in ("electrons", "positrons"):
        ts, cs = [], []
        for d in dirs:
            t, c = spec2d.by_ky(f"{d}/Spectral_2D_E_parallel_{sp}.dat"); ts.append(t); cs.append(c)
        t = np.concatenate(ts); c = np.concatenate(cs); m = (t >= a) & (t < b); out = out + c[m].mean(0)
    return -out
for name, dirs, nx in (("hxoff nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64),
                       ("nx128_hxoff", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    R = fg.Run(dirs, nx)
    print(f"\n{name}:   window     n   predicted P_n/P_1    GENE P_n/P_1    (fit residual of k_y = n)")
    for tc in (3000.0, 3600.0, 4200.0):
        t, P = R.window(tc, kys=tuple(range(NMAX + 1)))
        p0 = (P[:, :, 0, :].mean(0) * w[None, :]).sum(1) / fg.CXY; U = R.to_x(p0, nxf, 1).real
        p1 = P[:, :, 1, :]; wf, _, _ = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
        tt = t - t.mean(); Pn = {}
        for n in range(1, NMAX + 1):
            Om = np.array([j * wf[0] + (n - j) * wf[1] for j in range(n + 1)])
            pn = (P[:, :, n, :] * w[None, None, :]).sum(2)                     # flute part of k_y = n [nt, nx]
            M = np.exp(-1j * np.outer(tt, Om)); A, *_ = np.linalg.lstsq(M, pn, rcond=None)
            resid = np.linalg.norm(pn - M @ A) / np.linalg.norm(pn)
            kn = n * k; tot = 0.0
            for j in range(n + 1):
                ph0 = R.to_x(A[j], nxf); wt = Om[j] - kn * U
                tot += np.mean(2 * (L2 * kn ** 2) * kn ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH, right=0.0))
            Pn[n] = (tot, resid)
        g = gene_par(dirs, nx, tc - 300, tc + 300)
        for n in range(2, NMAX + 1):
            print(f"            {tc:6.0f}     {n}       {Pn[n][0]/Pn[1][0]:.3f}            {g[n]/g[1]:.3f}          ({Pn[n][1]:.2f})")
        print(f"            {tc:6.0f}   sum 2-{NMAX}   {sum(Pn[n][0] for n in range(2, NMAX+1))/Pn[1][0]:.3f}            {g[2:NMAX+1].sum()/g[1]:.3f}    all k_y>=2 in GENE: {g[2:].sum()/g[1]:.3f}")

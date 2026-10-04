r"""Debye-length scan (4 Oct 2026): loss rate of the vortices to entropy through parallel streaming, measured
(-Spectral_2D_E_parallel / Spectral_2D_E_fe, k_y >= 1, species summed) against the loss formula (eq:loss of the notes;
D(wt) from cache.npz, lambda_D^2 = debye2/2 of each run), in windows after the condensation of each run.
Runs: rb_c0p3_hxy (debye2 1e4), rb_c0p3_debye_x4 (4e4), rb_c0p3_debye_d4 (2.5e3)."""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import vortex_field as vf, bounce, predict, spec2d
RB = vf.RB
RUNS = {"ref (lamD^2 5000)": ([f"{RB}/rb_c0p3_hxy/leg_000{l}/out" for l in (1, 3)], 5000.0, (1500.0, 2500.0, 3500.0)),
        "debye_x4 (lamD^2 20000)": (sorted(glob.glob(f"{RB}/rb_c0p3_debye_x4/leg_*/out")), 20000.0, None),
        "debye_d4 (lamD^2 1250)": (sorted(glob.glob(f"{RB}/rb_c0p3_debye_d4/leg_*/out")), 1250.0, None)}
for key in RUNS:
    d_, L_, t_ = RUNS[key]; RUNS[key] = ([x for x in d_ if os.path.exists(x + "/plunk_e_time.dat") and os.path.getsize(x + "/field.dat") > 0], L_, t_)
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
geo = bounce.geometry(f"{RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
NXF = predict.NXF; dx = vf.LX / NXF; k = vf.KYMIN
for name, (dirs, L2, tcs) in RUNS.items():
    # measured parallel-channel rate, k_y >= 1
    T, F, P = [], [], []
    for d in dirs:
        try:
            for sp in ("electrons", "positrons"):
                t, f = spec2d.by_ky(f"{d}/Spectral_2D_E_fe_{sp}.dat"); _, p = spec2d.by_ky(f"{d}/Spectral_2D_E_parallel_{sp}.dat")
                if sp == "electrons": fe, pa = f[:, 1:].sum(1), p[:, 1:].sum(1)
                else: fe, pa = fe + f[:, 1:].sum(1), pa + p[:, 1:].sum(1)
            T.append(t); F.append(fe); P.append(pa)
        except Exception as e:
            print("  skip", d, e)
    t = np.concatenate(T); o = np.argsort(t); t, fe, pa = t[o], np.concatenate(F)[o], np.concatenate(P)[o]
    vf.RUNS[name] = dirs
    if tcs is None:
        # windows: from 300 t.u. after the first time the run's zonal fraction exceeds 0.9, three windows
        zz = []
        for d in dirs:
            a = np.loadtxt(d + "/plunk_e_time.dat", comments="#"); zz.append(a)
        a = np.vstack(zz); a = a[np.argsort(a[:, 0])]; frac = a[:, 10] / (a[:, 10] + a[:, 11])
        t0 = a[np.argmax(frac > 0.9), 0] + 300.0; t1 = a[-1, 0] - 50.0
        tcs = tuple(np.linspace(t0, t1, 3)) if t1 > t0 + 200 else (0.5 * (t0 + t1),)
    print(f"\n{name}:")
    for tc in tcs:
        m = (t >= tc - 60) & (t < tc + 60)
        meas = -pa[m].mean() / fe[m].mean() if m.sum() else np.nan
        tt, p0, p1 = predict.window_field(name, tc)
        if len(tt) < 20: print(f"   t {tc:.0f}: too few frames"); continue
        # frequency guess from the spectrum of the k_y = 1 mid-plane potential (strongest radial mode)
        sig = p1[:, :, vf.ZMID]; j = np.argmax((np.abs(sig) ** 2).sum(0)); x = sig[:, j] - sig[:, j].mean()
        ts = np.linspace(tt[0], tt[-1], 2048); xs = np.interp(ts, tt, x.real) + 1j * np.interp(ts, tt, x.imag)
        fr = np.fft.fftfreq(len(ts), ts[1] - ts[0]) * 2 * np.pi; Sp = np.abs(np.fft.fft(xs)) ** 2
        wpos = fr[np.argmax(np.where(fr > 0.05, Sp, 0))]; wneg = fr[np.argmax(np.where(fr < -0.05, Sp, 0))]
        w, _, res = vf.two_freq_fit(tt, p1[:, :, vf.ZMID], (-wneg, -wpos))
        Mx = np.exp(-1j * np.outer(tt - tt.mean(), w)); aa, *_ = np.linalg.lstsq(Mx, p1.reshape(len(tt), -1), rcond=None); aa = aa.reshape(2, vf.NX, vf.NZ)
        U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
        num = den = 0.0
        for kk in range(2):
            ph, dph = predict.to_x(aa[kk]), predict.to_x(aa[kk], deriv=1)
            ph0 = (ph * wz[None, :]).sum(1); dph0 = (dph * wz[None, :]).sum(1)
            wt = w[kk] - k * U
            num += 2 * (L2 * k ** 2) * k ** 2 * (np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH)).sum() * dx
            den += (Gxx * np.abs(dph0) ** 2 + Gyy * k ** 2 * np.abs(ph0) ** 2).sum() * dx
        pred = num / den
        print(f"   t {tc:7.0f}: w {w[0]:+.3f} {w[1]:+.3f} (fit resid {res:.3f}); loss rate measured {meas:.3e}, predicted {pred:.3e}, ratio {meas/pred:.2f}")

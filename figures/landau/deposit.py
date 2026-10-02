r"""Momentum deposited in the jets by the Landau damping of the vortices (2 Oct 2026).

Particles that absorb power P(x) from a wave of Doppler-shifted frequency wt = w - k U(x) absorb its momentum
k P / wt with it.  In the magnetised pair plasma a binormal force f drives a radial current f/B, whose divergence
changes the zonal charge: the jets feel the force density f(x) = k P(x) / wt(x), dE_zon/dt = <U f>.
P(x) is the loss density of the loss formula (bounce.Orbits.G), normalised to the predicted rate of each vortex.
Compared with the zonal forcing in GENE that 2D Euler does not contain (observed minus Reynolds stress,
scripts/rayleigh/zonal_forcing2.py, cache zforce2_hxoff.npz), on the jet harmonics.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf
import bounce, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
RAY = os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache")
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat")
O = bounce.Orbits(geo)
wz = geo["J"] / geo["J"].sum(); GXX = (wz * geo["gxx"]).sum(); GYY = (wz * geo["gyy"]).sum()
NXF = predict.NXF; LAMD2 = predict.LAMD2; K2L2 = LAMD2 * vf.KYMIN ** 2
LC = np.load(os.path.join(HERE, "cache.npz"))
Zc = np.load(os.path.join(RAY, "zforce2_hxoff.npz")); rows, F, Z = Zc["rows"], Zc["F"], Zc["Z"]
tcz = rows[:, 0]; kxg = np.fft.fftfreq(vf.NX, 1.0 / vf.NX) * 2 * np.pi / vf.LX
out = {}
for (a, b) in ((1000, 2000), (2200, 3200), (3400, 4400)):
    tc = 0.5 * (a + b)
    t, p0, p1 = predict.window_field("hxoff", tc)
    j = np.argmin(abs(rows[:, 0] - tc)); w0 = (rows[j, 1], -rows[j, 1])
    w, amid, res = vf.two_freq_fit(t, p1[:, :, vf.ZMID], w0)
    M = np.exp(-1j * np.outer(t - t.mean(), w))
    aa, *_ = np.linalg.lstsq(M, p1.reshape(len(t), -1), rcond=None); aa = aa.reshape(2, vf.NX, vf.NZ)
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    f = np.zeros(NXF); Ptot = 0.0; Elab = 0.0
    dx = vf.LX / NXF
    for k in range(2):
        ph, dph = predict.to_x(aa[k]), predict.to_x(aa[k], deriv=1)
        ph0 = (ph * wz[None, :]).sum(axis=1); dph0 = (dph * wz[None, :]).sum(axis=1)
        wt = w[k] - vf.KYMIN * U
        dens = 2 * K2L2 * vf.KYMIN ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), LC["wg"], LC["Dth"])     # loss density (notes, eq. loss)
        den = ((GXX * np.abs(dph0) ** 2 + GYY * vf.KYMIN ** 2 * np.abs(ph0) ** 2)).sum() * dx
        rate = dens.sum() * dx / den
        a1 = ph0 / vf.CXY
        Ek = np.mean(GXX * np.abs(np.gradient(a1, dx)) ** 2 + GYY * vf.KYMIN ** 2 * np.abs(a1) ** 2)       # fluid energy of the vortex (per x average)
        P = rate * Ek * dens / dens.mean()
        print(f"   vortex {k}: predicted rate {rate:.3e}, fluid energy {Ek:.3e}")
        f += vf.KYMIN * P / wt
        Ptot += P.mean(); Elab += (w[k] / wt * P).mean()
    fk = np.fft.fft(f) / NXF
    dzeta_kin = np.zeros(vf.NX, complex); h = vf.NX // 2
    dzeta_kin[:h] = (1j * kxg * np.r_[fk[:h], np.zeros(vf.NX - h)])[:h]
    dzeta_kin[vf.NX - h + 1:] = (1j * kxg * np.r_[np.zeros(h + 1), fk[NXF - h + 1:]])[vf.NX - h + 1:]
    # GENE: observed and Euler zonal vorticity tendency over the window
    ia, ib = np.argmin(abs(tcz - a)), np.argmin(abs(tcz - b))
    obs = (Z[ib] - Z[ia]) / (tcz[ib] - tcz[ia]); eul = F[ia:ib + 1].mean(0); resid = obs - eul
    ps0 = (p0.mean(axis=0) * wz[None, :]).sum(axis=1) / vf.CXY
    er = lambda dz_: float(np.sum((np.conj(ps0) * (-dz_ / np.where(kxg == 0, 1, -GXX * kxg ** 2) * (kxg != 0)) * GXX * kxg ** 2).real))
    print(f"\nwindow {a}-{b}: w {w[0]:+.3f} {w[1]:+.3f}, fit resid {res:.2f}; absorbed power <P> {Ptot:.3e} (lab-frame wave loss {Elab:.3e}); "
          f"jets' energy rate: kinetic force <U f> {np.mean(U*f):+.3e}; GENE observed {er(obs):+.3e}, Euler {er(eul):+.3e}, non-Euler residual {er(resid):+.3e}")
    for mm in (1, 3, 5, 7):
        r = dzeta_kin[mm] / resid[mm]
        print(f"   m {mm}: kinetic / non-Euler residual = {abs(r):.2f} at phase {np.degrees(np.angle(r)):+5.0f} deg   (|resid| {abs(resid[mm]):.2e}, |Euler| {abs(eul[mm]):.2e})")
    out[f"{a}"] = dict(f=f, U=U, dzk=dzeta_kin, resid=resid, eul=eul, obs=obs)
np.save(os.path.join(HERE, "deposit_cache.npy"), out, allow_pickle=True)

r"""Drain of the vortices from the pseudomomentum (wave-activity) budget of the Rayleigh mode (4 Oct 2026).

Linearised Euler about the jets, with the damping as a vorticity source S on the vortex:
    ik[(U - c) zeta - Z_x psi] = S,   Z_x = G^xx U''   (fluid units: psi = phi/C_xy, U = d_x psi_zon).
Kelvin's wave activity A = <zeta'^2>/(2 Z_x) obeys dA/dt + d_x(flux) = <zeta' S>/Z_x, so the Reynolds stress only
redistributes it and   d/dt int A = int Re(psi* S)/(2(U - c)).   The damping takes from the vortex the lab-frame power
(w/wt) P(x) = Re(psi* S)/2 (Sec. 3.5), so   d/dt int A = -int k w P / wt^2,   and with the mode relation zeta = Z_x psi/(U-c)
    Gamma_A = -(d/dt int A)/int A = int k w P/wt^2 dx  /  int k^2 Z_x |psi|^2 / wt^2 dx      (per vortex; P in the
normalisation of the loss formula, E = <G^xx|psi'|^2 + G^yy k^2|psi|^2>).  Inside the cat's eyes the fluid is mixed and
carries no wave activity: |x - x_c| < f_w w is excluded, f_w = 0.5, 1, 2 to test the sensitivity.
Compared with the measured Gamma_E of the vortices (both, flute, slope of E_nz) in the runs without radial hyperdiffusion."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; Gxx, Gyy = G["Gxx"], G["Gyy"]
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
L2, k = fg.LAM2, fg.KYMIN; nxf = 512


def terms(R, tc, fws=(0.5, 1.0, 2.0)):
    t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
    wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
    psz = (p0.mean(0) * w[None, :]).sum(1)
    U = R.to_x(psz, nxf, 1).real; U1 = R.to_x(psz, nxf, 2).real; U2 = R.to_x(psz, nxf, 3).real
    X = np.arange(nxf) * fg.LX / nxf; dx = X[1]
    out = {"E": 0.0, "P": 0.0, "Wlab": 0.0}
    for fw in fws:
        out[fw] = 0.0
    for kk in range(2):
        a0 = (a[kk] * w[None, :]).sum(1); ps = R.to_x(a0, nxf); dps = R.to_x(a0, nxf, 1)
        om = wf[kk]; c = om / k; wt = om - k * U
        Pd = 2 * (L2 * k ** 2) * k ** 2 * np.abs(ps) ** 2 * np.interp(np.abs(wt), WG, DTH)
        Ev = np.mean(Gxx * np.abs(dps) ** 2 + Gyy * k ** 2 * np.abs(ps) ** 2)
        out["E"] += Ev; out["P"] += Pd.mean(); out["Wlab"] += np.mean(om / wt * Pd)
        # critical layers of this vortex and their cat's-eye half-widths
        s = np.sign(U - c); idx = np.where(s != np.roll(s, -1))[0]
        mask = {fw: np.ones(nxf, bool) for fw in fws}
        for i in idx:
            j = (i + 1) % nxf; xc = X[i] + dx * (c - U[i]) / (U[j] - U[i] + 1e-300)
            Sh = abs(U1[i]); A = 2 * abs(ps[i]); wc = 2 * np.sqrt(A / max(Sh, 1e-300))
            dist = np.abs((X - xc + fg.LX / 2) % fg.LX - fg.LX / 2)
            for fw in fws:
                mask[fw] &= dist > fw * wc
        num_den = {}
        for fw in fws:
            m = mask[fw]
            num = np.mean(np.where(m, k * om * Pd / wt ** 2, 0.0))
            den = np.mean(np.where(m, k ** 2 * Gxx * U2 * np.abs(ps) ** 2 / wt ** 2, 0.0))
            out[fw] += (num / den) * Ev            # Gamma_A of this vortex times its energy = its predicted drain
        out.setdefault("nlayers", []).append(len(idx))
    return out


CASES = [("hxoff nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64, ((2700.0, 3300.0), (3300.0, 3900.0), (3900.0, 4400.0), (5000, 5600), (6500, 7100))),
         ("nx128_hxoff", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128, ((2700.0, 3300.0), (3300.0, 3900.0), (3900.0, 4400.0)))]
for name, dirs, nx, wins in CASES:
    R = fg.Run(dirs, nx)
    print(f"\n{name}:  window      Gamma_E meas    pred (cut 0.5w, 1w, 2w)              L/E      meas/(L/E)  pred(1w)/(L/E)   layers")
    for (t0, t1) in wins:
        if t1 > R.t[-1]:
            continue
        J = np.where((R.t >= t0) & (R.t <= t1))[0][::5]
        En = []
        for j in J:
            tt, phi = R.frame(j); pk = phi[:, 1, :] / fg.CXY; a0 = (pk * w[None, :]).sum(1)
            En.append((tt, np.mean(Gxx * np.abs(R.to_x(a0, nxf, 1)) ** 2 + Gyy * k ** 2 * np.abs(R.to_x(a0, nxf)) ** 2)))
        En = np.array(En); q = np.polyfit(En[:, 0], np.log(En[:, 1]), 1); Gm = -q[0]
        tw = [terms(R, tc) for tc in np.linspace(t0 + 60, t1 - 60, 4)]
        E = np.mean([x["E"] for x in tw]); LE = np.mean([x["P"] for x in tw]) / E
        pr = [np.mean([x[fw] for x in tw]) / E for fw in (0.5, 1.0, 2.0)]
        print(f"           {t0:.0f}-{t1:.0f}   {Gm:+.3e}     {pr[0]:+.3e} {pr[1]:+.3e} {pr[2]:+.3e}    {LE:.3e}   {Gm/LE:+.3f}       {pr[1]/LE:+.3f}        {tw[0]['nlayers']}")

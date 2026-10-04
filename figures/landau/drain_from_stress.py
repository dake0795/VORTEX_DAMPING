r"""Energy budget of the jets, all harmonics (4 Oct 2026): measured dE_zon/dt (flute zonal energy, fluid units) against
the theory = local (non-flute-resolved) Reynolds stress of the two fitted vortices + momentum deposited by their
Landau damping (loss formula), + hyp_y on the jets (zero for k_y = 0); at nx 64 and nx 128 without radial
hyperdiffusion.  If both resolutions are captured, the resolution dependence of the drain lies in the measured
structure, which the stress then converts correctly."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import fieldgen as fg
G = fg.geometry(); w = G["w"]; gxx, gyy, Gxx, Gyy = G["gxx"], G["gyy"], G["Gxx"], G["Gyy"]
LC = np.load(os.path.join(HERE, "cache.npz")); WG, DTH = LC["wg"], LC["Dth"]
L2, k = fg.LAM2, fg.KYMIN
def tend(R, pat, gx, gy):
    """zonal tendency of rho = (gx kx^2 + gy k^2) psi from the k_y = 1 pattern coefficients pat[nx] (psi = pat e^{iky} + cc)"""
    nx = R.nx; nxf = 3 * nx // 2 * 2
    f = R.to_x(pat, nxf); fx = R.to_x(pat, nxf, 1)
    r = -(gx * R.to_x(pat, nxf, 2) - gy * k ** 2 * f); rx = -(gx * R.to_x(pat, nxf, 3) - gy * k ** 2 * fx)
    T = -2 * np.real(fx * np.conj(1j * k * r) - 1j * k * f * np.conj(rx))           # -[psi_x rho_y - psi_y rho_x]_0, real space
    c = np.fft.fft(T) / nxf; out = np.zeros(nx, complex); h = nx // 2
    out[:h] = c[:h]; out[nx - h + 1:] = c[nxf - h + 1:]
    return out                                                                         # d rho_0/dt coefficients
def fdep(R, a, wf, psz, U, nxf):
    f = np.zeros(nxf)
    for kk in range(2):
        a0 = (a[kk] * w[None, :]).sum(1); ph0 = R.to_x(a0, nxf); wt = wf[kk] - k * U
        P_ = 2 * (L2 * k ** 2) * k ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH)
        f += k * P_ / np.where(np.abs(wt) < 1e-6, 1e-6, wt)
    c = np.fft.fft(f) / nxf; nx = R.nx; h = nx // 2; fc = np.zeros(nx, complex); fc[:h] = c[:h]; fc[nx - h + 1:] = c[nxf - h + 1:]
    dr = 1j * R.kx * fc                                     # d rho_0/dt from the force (rho = K psi convention, sign as <U f>)
    mm = np.abs(np.round(R.kx / (2 * np.pi / fg.LX))).astype(int)
    tot = float(np.sum(np.real(np.conj(psz) * dr))); lo = float(np.sum(np.real(np.conj(psz) * dr)[mm <= 4]))
    return (tot, lo)
def analyse(R, tc):
    t, P = R.window(tc); p0, p1 = P[:, :, 0, :] / fg.CXY, P[:, :, 1, :] / fg.CXY
    wf, a, res = fg.two_freq(t, p1, fg.freq_guess(t, p1[:, :, fg.ZMID]))
    psz = (p0.mean(0) * w[None, :]).sum(1)
    K = Gxx * R.kx ** 2
    nxf = 512; U = R.to_x(psz, nxf, 1).real; dx = fg.LX / nxf
    loc = np.zeros(R.nx, complex); flu = np.zeros(R.nx, complex); dep = 0.0; Lsum = 0.0; Enz = 0.0
    for kk in range(2):
        for iz in range(fg.NZ):
            loc += w[iz] * tend(R, a[kk][:, iz], gxx[iz], gyy[iz])
        a0 = (a[kk] * w[None, :]).sum(1); flu += tend(R, a0, Gxx, Gyy)
        ph0 = R.to_x(a0, nxf); dph0 = R.to_x(a0, nxf, 1); wt = wf[kk] - k * U
        P_ = 2 * (L2 * k ** 2) * k ** 2 * np.abs(ph0) ** 2 * np.interp(np.abs(wt), WG, DTH)       # loss density, flute-fluid units (same as eq:loss numerator)
        e_ = Gxx * np.abs(dph0) ** 2 + Gyy * k ** 2 * np.abs(ph0) ** 2
        Lsum += P_.mean(); Enz += e_.mean()
        dep += np.mean(U * k * P_ / np.where(np.abs(wt) < 1e-6, 1e-6, wt))                     # <U f>, f = k P / wt
    # jets' energy rate from a zonal charge tendency: E_zon = 1/2 sum K |psi|^2 per mode... use the same normalisation as e_ (per-x mean of |psi'|^2)
    mm = np.abs(np.round(R.kx / (2 * np.pi / fg.LX))).astype(int)
    rate = lambda dr, sel=None: float(np.sum((np.real(np.conj(psz) * dr))[(mm <= 4) if sel == "lo" else (mm >= 5) if sel == "hi" else slice(None)]))
    Ez = 0.5 * float(np.sum(K * np.abs(psz) ** 2))
    # deposit by harmonic: force density f(x) -> d rho_0/dt = i k_x f  (coefficients)
    return dict(t=t.mean(), w=wf[0], res=res, loc=rate(loc), flu=rate(flu), dep=dep, L=Lsum, Enz=Enz, Ez=Ez,
                loc_lo=rate(loc, "lo"), loc_hi=rate(loc, "hi"), flu_lo=rate(flu, "lo"), fdep=fdep(R, a, wf, psz, U, nxf))
CASES = [("hxoff nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64), ("nx128_hxoff", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)]
for name, dirs, nx in CASES:
    R = fg.Run(dirs, nx)
    # measured jets' energy (same normalisation) through the window, from frames every ~10 t.u.
    sel = np.arange(0, len(R.t), max(1, int(len(R.t) / max(1, (R.t[-1] - R.t[0]) / 10.0))))
    print(f"\n{name}:")
    for (t0, t1) in ((3000.0, 3600.0), (3600.0, 4300.0)):
        J = [j for j in sel if t0 <= R.t[j] <= t1]
        Ezs = []
        for j in J:
            tt, phi = R.frame(j); ps = (phi[:, 0, :] * w[None, :]).sum(1) / fg.CXY
            Ezs.append((tt, 0.5 * np.sum(Gxx * R.kx ** 2 * np.abs(ps) ** 2)))
        Ezs = np.array(Ezs); meas = np.polyfit(Ezs[:, 0], Ezs[:, 1], 1)[0]
        an = [analyse(R, tc) for tc in np.linspace(t0 + 60, t1 - 60, 4)]
        loc = np.mean([x["loc"] for x in an]); flu = np.mean([x["flu"] for x in an]); dep = np.mean([x["dep"] for x in an])
        L = np.mean([x["L"] for x in an]); Enz = np.mean([x["Enz"] for x in an])
        llo = np.mean([x["loc_lo"] for x in an]); lhi = np.mean([x["loc_hi"] for x in an]); flo = np.mean([x["flu_lo"] for x in an])
        dtot = np.mean([x["fdep"][0] for x in an]); dlo = np.mean([x["fdep"][1] for x in an])
        print(f"  t {t0:.0f}-{t1:.0f}: m<=4: local {llo:+.3e} + deposit {dlo:+.3e} = {llo+dlo:+.3e}  (flute {flo:+.3e} + deposit = {flo+dlo:+.3e});  m>=5 local {lhi:+.3e};  deposit spectral total {dtot:+.3e} vs <U f> {dep:+.3e}")
        print(f"  t {t0:.0f}-{t1:.0f}: jets dE/dt measured {meas:+.3e};  theory local+deposit {loc+dep:+.3e} [local {loc:+.3e}, deposit {dep:+.3e}];  flute+deposit {flu+dep:+.3e};  vortex loss L {L:.3e} (L/E_nz {L/Enz:.2e})")

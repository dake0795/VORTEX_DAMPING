"""Reynolds-stress energy input into the jets' m = 1 by binormal harmonic k_y = 1..5 (flute, field-line averaged
potential) and for k_y = 1 also computed locally along the line (local potential, local metric, then averaged),
from the raw field frames (no fit), t 2800-4300, every 5th frame, both radial resolutions, mean +- s.e.
Compare with GENE's E_nonlinear at k_y = 0, m = 1 (fluid units via E_fe / E, as jet_budget_clean.py)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import fieldgen as fg, spec2d
G = fg.geometry(); w = G["w"]; gxx, gyy, Gxx, Gyy = G["gxx"], G["gyy"], G["Gxx"], G["Gyy"]; k1 = fg.KYMIN
def tend(R, pat, gx, gy, k):
    nx = R.nx; nxf = 3 * nx // 2 * 2
    f = R.to_x(pat, nxf); fx = R.to_x(pat, nxf, 1)
    r = -(gx * R.to_x(pat, nxf, 2) - gy * k ** 2 * f); rx = -(gx * R.to_x(pat, nxf, 3) - gy * k ** 2 * fx)
    T = -2 * np.real(fx * np.conj(1j * k * r) - 1j * k * f * np.conj(rx))
    c = np.fft.fft(T) / nxf; out = np.zeros(nx, complex); h = nx // 2; out[:h] = c[:h]; out[nx - h + 1:] = c[nxf - h + 1:]
    return out
for case, dirs, NX in (("nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64), ("nx128", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128)):
    R = fg.Run(dirs, NX); mm = np.abs(np.round(R.kx / (2 * np.pi / fg.LX))).astype(int)
    J = np.where((R.t >= 2800) & (R.t <= 4300))[0][::5]
    rows = []
    for j in J:
        t, phi = R.frame(j); phi = phi / fg.CXY
        psz = (phi[:, 0, :] * w[None, :]).sum(1)
        row = []
        for ky in range(1, 6):
            a0 = (phi[:, ky, :] * w[None, :]).sum(1)
            row.append(np.real(np.conj(psz) * tend(R, a0, Gxx, Gyy, ky * k1))[mm == 1].sum())
        loc = sum(w[iz] * tend(R, phi[:, 1, iz], gxx[iz], gyy[iz], k1) for iz in range(fg.NZ))
        row.append(np.real(np.conj(psz) * loc)[mm == 1].sum())
        Em = 0.5 * Gxx * (R.kx ** 2 * np.abs(psz) ** 2)[mm == 1].sum(); row.append(Em)
        rows.append(row)
    A = np.array(rows); n = len(A); m_, s_ = A.mean(0), A.std(0) / np.sqrt(n)
    # GENE m = 1 nonlinear in fluid units
    spec2d.NX = NX; spec2d.NKX = NX - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
    nl = fe = 0; cnt = 0
    for d in dirs:
        for term in ("nonlinear", "fe"):
            for sp in ("electrons", "positrons"):
                kx, ky_, ts, _, c = spec2d.read(f"{d}/Spectral_2D_E_{term}_{sp}.dat"); sel = (ts >= 2800) & (ts <= 4300)
                mG = np.abs(np.round(kx / (2 * np.pi / fg.LX))).astype(int) == 1
                v = c[sel][:, mG, 0].sum(1)
                if term == "nonlinear": nl = nl + v.sum()
                else: fe = fe + v.sum()
            cnt += sel.sum() if term == "fe" else 0
    fac = (fe / cnt) / m_[-1]; gene = (nl / cnt) / fac
    print(f"{case}: {n} frames; m = 1 energy input from the flute stress of k_y = 1..5: " + "  ".join(f"{m_[i]:+.3f}+-{s_[i]:.3f}" for i in range(5)))
    print(f"        k_y = 1 local stress {m_[5]:+.3f}+-{s_[5]:.3f};  GENE nonlinear m = 1: {gene:+.3f}")

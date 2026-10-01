r"""Prediction of the loss rate L/E_nz of the electrostatic energy of the vortices by bounce-resonant (Landau) damping of
their non-flute part, from the measured potential (no free parameter):

    L/E = 4 pi^{3/2} C_xy v_T  sum_pm int dx G_pm(x)  /  ( lambda_D^2 sum_pm int dx int J dz [gxx |d_x phi_pm|^2 + gyy k^2 |phi_pm|^2] ),

G from bounce.Orbits.G with the Doppler-shifted frequency wt_pm(x) = omega_pm - k U(x) and the z-harmonics of the
potential phi_pm(x, z) of each of the two counter-propagating vortices (two-frequency fit of the k_y = 1 potential at
every z in a 60-time-unit window, frequencies from the mid-plane fit).
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf
import bounce

NXF = 256
LAMD2 = 5000.0


def window_field(run, tc, half=30.0):
    """All frames of `run` within tc +- half: t[nt], phi0[nt, nx, nz], phi1[nt, nx, nz]."""
    T, P0, P1 = [], [], []
    last = -np.inf
    for out in vf.RUNS[run]:
        path = out + "/field.dat"
        ts = vf.frame_times(path)
        for i in np.where((ts >= tc - half) & (ts <= tc + half) & (ts > last + 1e-9))[0]:
            t, phi = vf.read_frame(path, i)
            T.append(t); P0.append(phi[:, 0, :].copy()); P1.append(phi[:, 1, :].copy())
        if len(ts):
            last = max(last, ts[-1])
    return np.array(T), np.array(P0), np.array(P1)


def to_x(c, deriv=0):
    """c[nx(kx), nz] -> f[NXF, nz]."""
    return vf.to_x(c.T, nxf=NXF, deriv=deriv).T


def predict(run, tc, O, geo, w0):
    t, p0, p1 = window_field(run, tc)
    if len(t) < 20:
        return None
    w, amid, res = vf.two_freq_fit(t, p1[:, :, vf.ZMID], w0)
    M = np.exp(-1j * np.outer(t - t.mean(), w))
    a, *_ = np.linalg.lstsq(M, p1.reshape(len(t), -1), rcond=None)
    a = a.reshape(2, vf.NX, vf.NZ)
    resid = np.sum(np.abs(p1.reshape(len(t), -1) - M @ a.reshape(2, -1)) ** 2) / np.sum(np.abs(p1) ** 2)
    wz = geo["J"] / geo["J"].sum()
    U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
    dz, dx = 2 * np.pi / vf.NZ, vf.LX / NXF
    num = den = 0.0; parts = []
    nf_num = nf_den = 0.0
    for k in range(2):
        ph, dph = to_x(a[k]), to_x(a[k], deriv=1)
        wt = w[k] - vf.KYMIN * U
        G, Gp, Gt = O.G(bounce.zfft(ph), wt)
        n_k = 4 * np.pi ** 1.5 * bounce.CXY * bounce.VT * G.sum() * dx
        d_k = LAMD2 * ((geo["gxx"][None, :] * np.abs(dph) ** 2 + geo["gyy"][None, :] * vf.KYMIN ** 2 * np.abs(ph) ** 2) * geo["J"][None, :]).sum() * dz * dx
        num += n_k; den += d_k
        parts.append((n_k / d_k, Gp.sum() / G.sum()))
        nf_num += (np.abs(ph - ph.mean(axis=1, keepdims=True)) ** 2).sum(); nf_den += (np.abs(ph) ** 2).sum()
    return dict(t=t.mean(), rate=num / den, w=w, resid=resid, rate_p=parts[0][0], rate_m=parts[1][0],
                fpass=parts[0][1], nonflute=nf_num / nf_den, Umax=U.max())


if __name__ == "__main__":
    geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat")
    O = bounce.Orbits(geo)
    C = np.load(os.path.join(HERE, "..", "channels", "cache.npz"))
    T = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)["orig"]
    for tc in (2000.0, 4000.0, 6000.0, 8000.0, 9500.0, 10300.0):
        i = np.argmin(np.abs(T[:, 0] - tc)); w0 = (T[i, 1], T[i, 2])
        r = predict("orig", tc, O, geo, w0)
        j = np.argmin(np.abs(C["tc_orig"] - tc))
        print(f"t {r['t']:8.1f}  w {r['w'][0]:+.3f} {r['w'][1]:+.3f}  resid {r['resid']:.3f}  predicted L/E {r['rate']:.3e} (+ {r['rate_p']:.3e}, - {r['rate_m']:.3e}; passing share {r['fpass']:.2f})"
              f"   measured -E_par/E {C['par_orig'][j]:.3e}   ratio {r['rate'] / C['par_orig'][j]:.2f}   nonflute power {r['nonflute']:.2e}", flush=True)

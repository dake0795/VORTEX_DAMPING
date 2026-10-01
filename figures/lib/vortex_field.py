r"""Field-data helpers shared by figures/geometry and figures/trapping (vortex-damping notes, 1 Oct 2026).

Definitions are those of scripts/audit_vortex (field.py, vstruct.py, island.py, series.py), collected in one place:

  * field.dat frame = 16 + 4 + nx*nky*nz*16 + 4 bytes; time (f8) at offset 4; phi complex128 in Fortran order
    (nx, nky, nz) at offset 20; k_x in FFT order.  Real field: phi = phi_0 + 2 Re sum_{ky>0} phi_ky e^{i ky y}.
  * profiles in x through f(x) = sum_k c_k e^{+i k_x x} (the convention of vstruct.toreal), here by a zero-padded FFT;
  * E x B streamfunction Psi = phi / C_xy, so the zonal flow is U = d_x Phi / C_xy and its shear S = d_x U;
  * cosine amplitude of a single k_y: A = 2 |phi_ky| / C_xy;
  * the two counter-propagating vortices are separated by the two-frequency least-squares fit
        phi_{ky=1}(x, t) = a_+(x) e^{-i w_+ t} + a_-(x) e^{-i w_- t},   w_+ > 0 > w_-,
    in 60-time-unit windows (the frequencies by variable projection, the profiles linearly);
  * critical layers: U(x_c) = c = w / k_y; pendulum estimates w_tr = k_y sqrt(|S| A), half-width w = 2 sqrt(A / |S|).

All measurements at the outboard mid-plane (z index NZ/2), as in the notes.
"""
import os
import numpy as np

RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
NX, NKY, NZ = 64, 16, 32
N = NX * NKY * NZ
FR = 16 + 4 + N * 16 + 4
LX, KYMIN, CXY = 3689.25, 0.002, 0.92640709665123333
LY = 2.0 * np.pi / KYMIN
KX = np.fft.fftfreq(NX, 1.0 / NX) * 2.0 * np.pi / LX
DX = LX / NX
LAMBDA_D = np.sqrt(1.0e4 / 2.0)             # 70.71 rho_ref
ZMID = NZ // 2
NXF = 1024                                   # fine x grid for profiles
XF = np.arange(NXF) * LX / NXF
T_RESTART = 5906.678243865176

RUNS = {"orig": [f"{RB}/rb_c0p3_hxy/leg_000{l}/out" for l in (1, 3, 4, 5, 6)],
        "a03": [f"{RB}/rb_c0p3_hxy_ampltest/amp0p3/out"],
        "a003": [f"{RB}/rb_c0p3_hxy_ampltest/amp0p03/out"]}


def frame_times(path):
    n = os.path.getsize(path) // FR
    ts = np.empty(n)
    with open(path, "rb") as f:
        for i in range(n):
            f.seek(i * FR + 4)
            ts[i] = np.frombuffer(f.read(8), "<f8")[0]
    return ts


def read_frame(path, i):
    """Full potential of frame i: (t, phi[nx, nky, nz])."""
    with open(path, "rb") as f:
        f.seek(i * FR)
        b = f.read(FR)
    return np.frombuffer(b[4:12], "<f8")[0], np.frombuffer(b[20:20 + N * 16], "<c16").reshape((NX, NKY, NZ), order="F")


def midplane_series(run, tmin=-np.inf, tmax=np.inf, iz=ZMID):
    """Mid-plane potential of every output frame of a run in [tmin, tmax]: t[nt], phi[nt, nx, nky] (legs spliced,
    the repeated first frame of each leg dropped)."""
    T, P = [], []
    last = -np.inf
    for out in RUNS[run]:
        path = out + "/field.dat"
        ts = frame_times(path)
        with open(path, "rb") as f:
            for i in np.where((ts >= tmin) & (ts <= tmax) & (ts > last + 1e-9))[0]:
                f.seek(i * FR + 20 + 16 * NX * NKY * iz)
                P.append(np.frombuffer(f.read(16 * NX * NKY), "<c16").reshape((NX, NKY), order="F"))
                T.append(ts[i])
        if len(ts):
            last = max(last, ts[-1])
    return np.array(T), np.array(P)


def to_x(c, nxf=NXF, deriv=0):
    """c[..., nx] in FFT order -> d^deriv f / dx^deriv on the fine grid XF (last axis)."""
    c = np.asarray(c) * (1j * KX) ** deriv
    pad = np.zeros(c.shape[:-1] + (nxf,), complex)
    h = NX // 2
    pad[..., :h] = c[..., :h]
    pad[..., nxf - h + 1:] = c[..., h + 1:]
    return np.fft.ifft(pad, axis=-1) * nxf


def two_freq_fit(t, p, w0):
    """Fit p[nt, nx] = a_+[nx] e^{-i w_+ t} + a_-[nx] e^{-i w_- t} (variable projection over the two frequencies).
    Returns (w_+, w_-), a[2, nx], relative residual.  Times are measured from the window centre."""
    from scipy.optimize import minimize
    tc = t - t.mean()
    norm = np.sum(np.abs(p) ** 2)

    def solve(w):
        M = np.exp(-1j * np.outer(tc, w))                          # nt x 2
        a, *_ = np.linalg.lstsq(M, p, rcond=None)
        return a, np.sum(np.abs(p - M @ a) ** 2) / norm

    r = minimize(lambda w: solve(w)[1], w0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-12))
    w = r.x
    a, res = solve(w)
    if w[0] < w[1]:
        w, a = w[::-1], a[::-1]
    return w, a, res


def crossings(U, c):
    """Fine-grid positions x_c (linear interpolation) and indices where U = c."""
    out = []
    for j in range(NXF):
        jn = (j + 1) % NXF
        if (U[j] - c) * (U[jn] - c) < 0:
            f = (c - U[j]) / (U[jn] - U[j])
            out.append((XF[j] + f * LX / NXF, j, f))
    return out


def interp(F, j, f):
    return F[j] + f * (F[(j + 1) % NXF] - F[j])


def window_analysis(t, phi, w0):
    """One window of mid-plane data: frequencies, per-vortex profiles and critical-layer quantities.
    t[nt], phi[nt, nx, nky].  Returns a dict."""
    w, a, res = two_freq_fit(t, phi[:, :, 1], w0)
    p0 = phi[:, :, 0].mean(axis=0)                                  # window-mean zonal potential
    U = to_x(p0, deriv=1).real / CXY
    S = to_x(p0, deriv=2).real / CXY
    out = dict(t=t.mean(), w=w, res=res, a=a, p0=p0, U=U, S=S, Umax=U.max(), Umin=U.min(),
               E=np.sum(np.abs(a) ** 2, axis=1))
    for k, name in enumerate(("p", "m")):
        A = 2.0 * np.abs(to_x(a[k])) / CXY
        c = w[k] / KYMIN
        cr = [(xc, abs(interp(S, j, f)), interp(A, j, f)) for xc, j, f in crossings(U, c)]
        out["A_" + name] = A
        out["c_" + name] = c
        out["crit_" + name] = cr
        if cr:
            Sc = np.mean([q[1] for q in cr])
            Ac = np.mean([q[2] for q in cr])
            out["S_" + name], out["Ac_" + name] = Sc, Ac
            out["wtr_" + name] = KYMIN * np.sqrt(Sc * Ac)
            out["half_" + name] = 2.0 * np.sqrt(Ac / Sc)
        else:
            out["S_" + name] = out["Ac_" + name] = out["wtr_" + name] = out["half_" + name] = np.nan
    return out

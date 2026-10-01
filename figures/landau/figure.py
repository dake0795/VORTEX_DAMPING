#!/usr/bin/env python3
r"""Figure `landau`: Landau damping of the non-flute part of the vortices - prediction against simulation (1 Oct 2026).

  python figure.py --rebuild   solve the kinetic field equation along the field line for a grid of Doppler-shifted
                               frequencies (response.py), read the potential of the reduced-box run and of the two
                               amplitude-reduced restarts in 60-time-unit windows (predict.py), write cache.npz
  python figure.py             draw ../landau.pdf from cache.npz

Theory (appendices/appH_landau.tex).  The non-flute part of a nearly flute-like vortex is phi_1 = k^2 lambda_D^2 chi phi_0,
with (1 - R_wt) chi = s(l), s = <g^yy> [g^xx/<g^xx> - g^yy/<g^yy>], R_wt the response of particles streaming along the
field line at the Doppler-shifted frequency wt.  The electrostatic energy is lost at the rate
    L/E_nz = 2 (k lambda_D)^2 k^2 sum_pm int dx |phi_0|^2 D(wt) / sum_pm int dx (<g^xx> |phi_0'|^2 + <g^yy> k^2 |phi_0|^2),
    D(wt) = - wt <Im chi(wt) s>     (the dissipation function; depends only on the geometry and on wt).
Measured counterpart of D at each radius: - wt <Im(phi_1/phi_0) s> / (k lambda_D)^2, from the potential of each vortex.
(a) D(wt): theory (line) and the measured values at all radii where the vortex is strong, both vortices, several times.
(b) L/E_nz predicted from the flute-like part of the vortices and the jets alone against the rate measured by the code's
    energy diagnostic in the parallel-streaming channel (../channels/cache.npz): the reduced-box run at every window, the
    two amplitude-reduced restarts, and the production run (box four times larger; predict_production.py).
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import notes_style as ns
import traces
import vortex_field as vf
import bounce, response, predict as P

CACHE = os.path.join(HERE, "cache.npz")
K2L2 = P.LAMD2 * vf.KYMIN ** 2
TWIN = {"orig": np.arange(1500.0, 10801.0, 620.0), "a03": np.array([6100.0, 6350.0, 6600.0]), "a003": np.array([6100.0, 6350.0, 6600.0])}


def rebuild():
    geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat")
    R = response.Response(geo)
    wz = geo["J"] / geo["J"].sum()
    Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
    s = Gyy * (geo["gxx"] / Gxx - geo["gyy"] / Gyy)
    wg = np.r_[np.linspace(0.02, 1.0, 26), np.linspace(1.08, 4.2, 40)]
    Dth = np.array([-w * (wz * np.imag(R.solve(w, 1.0)) * s).sum() for w in wg])
    print("theory D(wt):", np.round(Dth[::6], 4), flush=True)
    out = dict(wg=wg, Dth=Dth)
    T = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)
    dx = vf.LX / P.NXF
    for run, tcs in TWIN.items():
        rows, pts = [], []
        tr = T[run]
        for tc in tcs:
            i = np.argmin(np.abs(tr[:, 0] - tc)); w0 = (tr[i, 1], tr[i, 2])
            t, p0, p1 = P.window_field(run, tc)
            if len(t) < 20:
                continue
            w, amid, res = vf.two_freq_fit(t, p1[:, :, vf.ZMID], w0)
            M = np.exp(-1j * np.outer(t - t.mean(), w))
            a, *_ = np.linalg.lstsq(M, p1.reshape(len(t), -1), rcond=None)
            a = a.reshape(2, vf.NX, vf.NZ)
            U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=P.NXF, deriv=1).real / vf.CXY
            num_th = num_me = den = 0.0
            for k in range(2):
                ph, dph = P.to_x(a[k]), P.to_x(a[k], deriv=1)
                ph0 = (ph * wz[None, :]).sum(axis=1); dph0 = (dph * wz[None, :]).sum(axis=1); ph1 = ph - ph0[:, None]
                wt = w[k] - vf.KYMIN * U
                Dloc = -wt * (np.imag(ph1 / ph0[:, None]) * s[None, :] * wz[None, :]).sum(axis=1) / K2L2
                D_th = np.interp(np.abs(wt), wg, Dth)
                num_th += 2 * K2L2 * vf.KYMIN ** 2 * (np.abs(ph0) ** 2 * D_th).sum() * dx
                num_me += 2 * K2L2 * vf.KYMIN ** 2 * (np.abs(ph0) ** 2 * Dloc).sum() * dx
                den += (Gxx * np.abs(dph0) ** 2 + Gyy * vf.KYMIN ** 2 * np.abs(ph0) ** 2).sum() * dx
                strong = np.abs(ph0) > 0.12 * np.abs(ph0).max()
                for j in np.where(strong)[0][::3]:
                    pts.append((t.mean(), k, wt[j], Dloc[j], abs(ph0[j]) / np.abs(ph0).max()))
            rows.append((t.mean(), num_th / den, num_me / den, w[0], w[1], res))
            print(run, f"t {t.mean():8.1f}: predicted L/E {num_th / den:.3e}   from the measured non-flute part {num_me / den:.3e}", flush=True)
        out[f"rows_{run}"] = np.array(rows); out[f"pts_{run}"] = np.array(pts)
    np.savez(CACHE, **out)


def draw():
    C = np.load(CACHE)
    M = np.load(os.path.join(HERE, "..", "channels", "cache.npz"))
    ns.apply_style()
    fig, (a, b) = ns.two_panel(wspace=0.80)
    # (a) the dissipation function
    pts = C["pts_orig"]
    sel = (pts[:, 0] > 3000) & (pts[:, 0] < 10000)
    a.plot(np.abs(pts[sel, 2]), pts[sel, 3], "o", ms=2.0, c="0.55", mec="none", zorder=2)
    a.plot(C["wg"], C["Dth"], c=ns.PAL["theory"], lw=1.6, zorder=3)
    a.set_xlim(0, 3.3); a.set_ylim(-0.04, 0.38)
    a.axhline(0, c="0.7", lw=0.6, zorder=1)
    a.set_xlabel(r"$|\tilde\omega|\;(c_{\rm ref}/L_{\rm ref})$")
    a.set_ylabel(r"$\mathcal{D}\;(c_{\rm ref}/L_{\rm ref})$")
    traces.label(a, 0.2, 0.22, "theory", c=ns.PAL["theory"]); traces.label(a, 0.95, 0.335, "simulation", c="0.35")
    ns.style_axes(a)
    # (b) predicted against measured loss rate, every run
    r = C["rows_orig"]
    pm = np.interp(r[:, 0], M["tc_orig"], M["par_orig"])
    lim = np.array([1.0e-4, 3.5e-3])
    b.loglog(lim, lim, c="0.6", lw=0.8, zorder=1)
    b.loglog(pm, r[:, 1], "o", c=ns.PAL["orig"], ms=3.4, zorder=3)
    for run, mk in (("a03", "s"), ("a003", "^")):
        rr = C[f"rows_{run}"]
        b.loglog(np.interp(rr[:, 0], M[f"tc_{run}"], M[f"par_{run}"]), rr[:, 1], mk, c=ns.PAL[run], mfc="w", ms=5, mew=1.0, zorder=4)
    Pr = np.load(os.path.join(HERE, "cache_production.npz"))
    b.loglog([float(Pr["measured_parallel"])], [float(Pr["predicted"])], "*", c=ns.PAL["aux"], ms=10, zorder=5)
    b.set_xlim(*lim); b.set_ylim(*lim)
    b.set_xlabel(r"measured $\mathcal{L}/E_{\rm nz}\;(c_{\rm ref}/L_{\rm ref})$")
    b.set_ylabel(r"predicted $\mathcal{L}/E_{\rm nz}\;(c_{\rm ref}/L_{\rm ref})$")
    traces.label(b, 2.2e-4, 1.35e-4, "larger box", c=ns.PAL["aux"]); traces.label(b, 3.3e-4, 2.2e-3, "reduced box", c="k")
    ns.tag(a, "(a)", x=0.86, y=0.14, va="bottom"); ns.tag(b, "(b)", x=0.86, y=0.14, va="bottom")
    ns.save(fig, "landau")
    for run in ("orig", "a03", "a003"):
        r = C[f"rows_{run}"]; tcm = M[f"tc_{run}"]; pm = M[f"par_{run}"]
        for row in r:
            j = np.argmin(np.abs(tcm - row[0]))
            print(f"{run:5s} t {row[0]:8.1f}  predicted {row[1]:.3e}  via measured phi_1 {row[2]:.3e}  measured (GENE) {pm[j]:.3e}  ratio pred/meas {row[1] / pm[j]:.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rebuild", action="store_true"); A = ap.parse_args()
    if A.rebuild or not os.path.exists(CACHE):
        rebuild()
    draw()

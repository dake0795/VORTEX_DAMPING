#!/usr/bin/env python3
r"""Figure `geometry` of the vortex-damping notes: where the vortices sit on the jets, and the cat's eyes.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../geometry.pdf (+ .png)
  python3 figure.py --rebuild                         read the field data and rewrite cache.npz first

Data (read only): mid-plane potential of the reduced-box run rb_c0p3_hxy over the 60-time-unit window centred on
T0 (leg 3 / leg 4).  Definitions in ../lib/vortex_field.py: the zonal flow U = d_x Phi / C_xy of the window-mean
k_y = 0 potential; the two vortices from the two-frequency fit a_+ e^{-i w_+ t} + a_- e^{-i w_- t} of the k_y = 1
potential; phase speeds c_pm = w_pm / k_y.  Panel (b): the streamfunction in the frame of the + vortex,
    H(x, xi) = Psi_zon(x) + 2 Re[a_+(x) e^{i k_y xi}] / C_xy - c_+ x,     xi = y - c_+ t,
with Psi_zon = Phi / C_xy; the - vortex is left out (in this frame it sweeps past at the breathing frequency).
The x axis is shifted by a quarter of the box so that neither jet is cut by the periodic boundary.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
CACHE = os.path.join(HERE, "cache.npz")
T0, WIN = 5910.0, 60.0


def rebuild():
    import vortex_field as vf
    t, p = vf.midplane_series("orig", T0 - WIN / 2, T0 + WIN / 2)
    o = vf.window_analysis(t, p, np.array([0.86, -0.86]))
    print(f"t = {o['t']:.1f}, w = {o['w']}, residual = {o['res']:.2e}, c/Ujet = {o['c_p'] / o['Umax']:.3f}, {o['c_m'] / o['Umin']:.3f}")
    np.savez_compressed(CACHE, t=o["t"], w=o["w"], res=o["res"], p0=o["p0"], a=o["a"],
                        crit_p=np.array(o["crit_p"]), crit_m=np.array(o["crit_m"]))


def saddles_and_centres(H0x, ap, apx, c_shift, xs, vf):
    """Fixed points of H near the critical layers: d_xi H = 0 fixes k_y xi = -arg a_+(x) + j pi; then solve d_x H = 0."""
    from scipy.optimize import brentq
    out = []
    for xc in xs:
        for j in (0, 1):
            def f(x):
                a, ax = ap(x), apx(x)
                ph = -np.angle(a) + j * np.pi
                return H0x(x) + 2.0 * np.real(ax * np.exp(1j * ph)) / vf.CXY
            lo, hi = xc - 4 * vf.DX, xc + 4 * vf.DX
            try:
                x = brentq(f, lo, hi)
            except ValueError:
                continue
            out.append((x, ((-np.angle(ap(x)) + j * np.pi) / vf.KYMIN) % vf.LY, j))
    return out


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--rebuild", action="store_true")
    if ap_.parse_args().rebuild:
        rebuild()
    import notes_style as ns
    import vortex_field as vf
    ns.apply_style()
    PAL, FS = ns.PAL, ns.FS
    C = np.load(CACHE)
    p0, a, w = C["p0"], C["a"], C["w"]
    lam = vf.LAMBDA_D
    nxf = vf.NXF
    sh = nxf // 4                                           # shift: x runs over [-L/4, 3L/4)
    x = (np.arange(nxf) - sh) * vf.LX / nxf
    roll = lambda f: np.roll(f, sh, axis=-1)
    U = roll(vf.to_x(p0, deriv=1).real / vf.CXY)
    Psi = roll(vf.to_x(p0).real / vf.CXY)
    Ap, Am = (roll(2.0 * np.abs(vf.to_x(a[k])) / vf.CXY) for k in (0, 1))
    cp, cm = w / vf.KYMIN
    U0 = np.abs(U).max()
    wrap = lambda xc: (xc + vf.LX / 4) % vf.LX - vf.LX / 4

    fig, (ax, bx) = ns.two_panel(height=2.27, left=0.62, wspace=0.70)      # (b) has equal scales in x and y

    # (a) jets, vortex profiles and phase speeds
    ax.axhline(0, c="0.8", lw=0.6)
    ax.plot(x / lam, U / U0, c="k")
    ax.plot(x / lam, Ap / Ap.max(), c=PAL["plus"])
    ax.plot(x / lam, -Am / Am.max(), c=PAL["minus"])
    for c, col, cr in ((cp, PAL["plus"], C["crit_p"]), (cm, PAL["minus"], C["crit_m"])):
        ax.axhline(c / U0, c=col, ls=":", lw=0.9)
        ax.plot(wrap(cr[:, 0]) / lam, np.full(len(cr), c / U0), "o", ms=4.0, mfc="w", mec=col, mew=1.0, zorder=5)
    ax.set_xlim(x[0] / lam, (x[-1] + vf.LX / nxf) / lam)
    ax.set_ylim(-1.12, 1.12)
    ax.set_xlabel(r"$x\;(\lambda_D)$")
    ax.set_ylabel(r"$U/U_{\max}$")
    ax.text(0.45, 0.50, r"$U$", ha="center", va="bottom")
    ax.text(14.2, 0.66, r"$|\psi_+|$", ha="left", va="center", color=PAL["plus"])
    ax.text(-7.3, -0.74, r"$-|\psi_-|$", ha="left", va="center", color=PAL["minus"])
    ax.text(26.0, cp / U0 + 0.04, r"$c_+$", ha="center", va="bottom", color=PAL["plus"])
    ax.text(4.9, cm / U0 - 0.05, r"$c_-$", ha="center", va="top", color=PAL["minus"])
    ns.style_axes(ax)

    # (b) streamlines in the frame of the + vortex
    ny = 400
    xi = np.arange(ny + 1) * vf.LY / ny
    apx_f = roll(vf.to_x(a[0]))
    H0 = Psi - cp * x
    H = H0[None, :] + 2.0 * np.real(apx_f[None, :] * np.exp(1j * vf.KYMIN * xi[:, None])) / vf.CXY
    # fixed points: continuous evaluation of the Fourier series (x measured in the unshifted coordinate)
    ev = lambda cf, d=0: (lambda xx: np.sum(cf * (1j * vf.KX) ** d * np.exp(1j * vf.KX * xx)))
    H0x = lambda xx: ev(p0, 1)(xx).real / vf.CXY - cp
    fp = saddles_and_centres(H0x, ev(a[0]), ev(a[0], 1), cp, C["crit_p"][:, 0], vf)
    Hval = lambda xx, y: ev(p0)(xx).real / vf.CXY - cp * wrap(xx) + 2.0 * np.real(ev(a[0])(xx) * np.exp(1j * vf.KYMIN * y)) / vf.CXY
    # background streamlines, evenly spaced in H; then, for each critical layer, the separatrix and the trapped orbits
    bx.contour(x / lam, xi / lam, H, levels=np.linspace(H.min(), H.max(), 26), colors="0.65", linewidths=0.5,
               linestyles="solid")
    for xc in C["crit_p"][:, 0]:
        pts = [(Hval(q[0], q[1]), q) for q in fp if abs(q[0] - xc) < 4 * vf.DX]
        if len(pts) < 2:
            continue
        m = np.abs(wrap(xc) - x) < 6 * vf.DX
        hx = np.array([p[0] for p in pts])
        # saddle (X-point) or centre (O-point) by the sign of the Hessian determinant at the fixed point
        dets = []
        for _, (xx, y, j) in pts:
            h = 2.0
            Hxx = (Hval(xx + h, y) - 2 * Hval(xx, y) + Hval(xx - h, y)) / h ** 2
            Hyy = (Hval(xx, y + h) - 2 * Hval(xx, y) + Hval(xx, y - h)) / h ** 2
            Hxy = (Hval(xx + h, y + h) - Hval(xx + h, y - h) - Hval(xx - h, y + h) + Hval(xx - h, y - h)) / (4 * h ** 2)
            dets.append(Hxx * Hyy - Hxy ** 2)
        isad, icen = int(np.argmin(dets)), int(np.argmax(dets))
        Hs, Ho = hx[isad], hx[icen]
        print(f"critical layer x = {xc:.0f}: saddle H = {Hs:.4e} (det {dets[isad]:.2e}), centre H = {Ho:.4e} "
              f"(w_tr = {np.sqrt(max(dets[icen], 0)):.3f}); X-point at x = {pts[isad][1][0]:.0f}, O-point at x = {pts[icen][1][0]:.0f}")
        sub = np.where(m)[0]
        xs_, Hs_ = x[sub] / lam, H[:, sub]
        lev = np.sort(Hs + (Ho - Hs) * np.array([0.25, 0.5, 0.75, 0.95]))
        bx.contourf(xs_, xi / lam, Hs_, levels=sorted([Hs, Ho + 0.01 * (Ho - Hs)]), colors=[PAL["plus"]], alpha=0.12)
        bx.contour(xs_, xi / lam, Hs_, levels=lev, colors=PAL["plus"], linewidths=0.7, linestyles="solid")
        bx.contour(xs_, xi / lam, Hs_, levels=[Hs], colors="k", linewidths=1.3, linestyles="solid")
        bx.axvline(wrap(xc) / lam, c="k", ls=":", lw=0.8, zorder=1)
        ax.axvline(wrap(xc) / lam, c="k", ls=":", lw=0.8, zorder=1)
    bx.set_xlim(*ax.get_xlim())
    bx.set_ylim(0, vf.LY / lam)
    bx.set_xlabel(r"$x\;(\lambda_D)$")
    bx.set_ylabel(r"$y-c_+t\;(\lambda_D)$")
    ns.style_axes(bx)
    ns.tag(ax, "(a)", x=0.955, ha="right")
    ns.tag(bx, "(b)", x=0.955, ha="right")
    ns.save(fig, "geometry")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""Figure `corrugation`: the vorticity strips of the jets through the run, and the growth rate they give (1 Oct 2026).

  python figure.py --rebuild   read field.dat of the reduced-box run (legs 1, 3-5) and the Rayleigh growth history
                               (../../scripts/rayleigh/cache/growth_history.npz); write cache.npz
  python figure.py             draw ../corrugation.pdf from cache.npz

Corrugation of the jets = the contrast of the vorticity strips: in each shear layer the zonal vorticity U' has two
extrema, on the critical layers of the two vortices, and a dip between them; the contrast is the mean |U'| at the two
strips minus |U'| at the middle of the layer, divided by the largest |U'|, averaged over the two shear layers.  U is
the mid-plane zonal flow averaged over 60 time units.
(a) the contrast against time, through the decay and after the collapse; (b) the Rayleigh growth rate gamma_R of the
k_y = 1 mode on the same jets against the contrast.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns
import traces
import vortex_field as vf

CACHE = os.path.join(HERE, "cache.npz")


def rebuild():
    legs = {l: f"{vf.RB}/rb_c0p3_hxy/leg_000{l}/out/field.dat" for l in (1, 3, 4, 5)}
    times = {l: vf.frame_times(p) for l, p in legs.items()}
    x = vf.XF
    tc = np.arange(900.0, 15901.0, 200.0); con = []
    for c in tc:
        acc, n = 0, 0
        for l, p in legs.items():
            for i in np.where((times[l] > c - 30) & (times[l] < c + 30))[0][::3]:
                acc = acc + vf.read_frame(p, i)[1][:, 0, vf.ZMID]; n += 1
        p0 = acc / n
        U = vf.to_x(p0, deriv=1).real / vf.CXY; S = vf.to_x(p0, deriv=2).real / vf.CXY
        # shear layers = between consecutive zeros of U'' ... located through the zeros of U (jet edges): the two layers
        # are centred on the zeros of U; strips = extrema of |U'| on either side of each centre
        vals = []
        for j in np.where(np.diff(np.sign(U)) != 0)[0]:
            w = int(0.14 * vf.NXF)                                   # search +- 0.14 box on either side of the centre
            idx = np.arange(j - w, j + w) % vf.NXF
            a = np.abs(S[idx]); mid = a[w]
            vals.append((0.5 * (a[:w].max() + a[w:].max()) - mid) / np.abs(S).max())
        con.append(np.mean(vals))
    G = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "growth_history.npz"), allow_pickle=True)["rows"]
    tg, g = G[:, 0], 0.5 * (G[:, 1] + G[:, 3])
    gR = np.array([g[(tg > c - 100) & (tg < c + 100)].mean() if ((tg > c - 100) & (tg < c + 100)).any() else np.nan for c in tc])
    np.savez(CACHE, tc=tc, contrast=np.array(con), gR=gR)


def draw():
    C = np.load(CACHE)
    tc, con, g = C["tc"], C["contrast"], C["gR"]
    ns.apply_style()
    fig, (a, b) = ns.two_panel()
    traces.collapse_band(a)
    a.plot(tc / 1e3, con, c="k")
    a.set_xlim(0, 16); a.set_ylim(0, 0.36)
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$"); a.set_ylabel(r"contrast of the strips")
    m = np.isfinite(g) & (tc <= traces.COLLAPSE[1] + 100)
    b.plot(con[m], 1e2 * np.maximum(g[m], 0.0), "o", c="k", ms=3.2)
    b.set_xlim(0.0, 0.36); b.set_ylim(-0.3, 9.5)
    b.set_xlabel(r"contrast of the strips"); b.set_ylabel(r"$\gamma_{\rm R}\;(10^{-2}c_{\rm ref}/L_{\rm ref})$")
    ns.style_axes(a); ns.style_axes(b)
    ns.tag(a, "(a)"); ns.tag(b, "(b)")
    ns.save(fig, "corrugation")
    for i in range(0, len(tc), 6):
        print(f"t {tc[i]:6.0f} contrast {con[i]:.3f} gamma_R {g[i]:.4f}")


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        rebuild()
    draw()

#!/usr/bin/env python3
r"""Figure `balance`: the energy of the vortices is the small difference of two large exchanges with the jets (4 Oct 2026).

  python figure.py --rebuild   fit the two vortices in windows of the runs without radial hyperdiffusion (nx 64:
                               rb_c0p3_hxoff; nx 128: rb_c0p3_draintest/nx128_hxoff + leg 2) and evaluate, per window,
                               (i) the Reynolds-stress exchange (energy the flute stress of the vortices passes to the jets),
                               (ii) the lab-frame effect of the Landau damping on the vortices, -<(w/wt) P>, and
                               (iii) the measured rate of change of the energy of the vortices; write cache.npz
  python figure.py             draw ../balance.pdf
(a) the three rates divided by the absorbed power P, against time, nx 64; (b) the same at nx 128.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "landau"))
CACHE = os.path.join(HERE, "cache.npz")


def rebuild():
    import fieldgen as fg
    globals()["fg"] = fg; globals()["np"] = np
    src = open(os.path.join(HERE, "..", "landau", "vortex_drain_pred.py")).read()
    exec(src[src.index("G = fg.geometry()"):src.index("CASES = [")].replace('os.path.join(HERE, "cache.npz")', 'os.path.join(HERE, "..", "landau", "cache.npz")').replace('os.path.join(HERE, "drain_from_stress.py")', 'os.path.join(HERE, "..", "landau", "drain_from_stress.py")'), globals())
    out = {}
    for key, dirs, nx, t0, t1 in (("nx64", [fg.RB + "/rb_c0p3_hxoff/leg_0001/out"], 64, 1200.0, 4500.0),
                                  ("nx128", [fg.RB + "/rb_c0p3_draintest/nx128_hxoff/out", fg.RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128, 2600.0, 5400.0)):
        R = fg.Run(dirs, nx); t1 = min(t1, R.t[-1] - 50)
        rows = []
        for a in np.arange(t0, t1 - 299, 300.0):
            b = a + 300.0
            J = np.where((R.t >= a) & (R.t <= b))[0][::4]
            En = []
            for j in J:
                tt, phi = R.frame(j); pk = phi[:, 1, :] / fg.CXY; a0 = (pk * w[None, :]).sum(1)
                En.append((tt, np.mean(Gxx * np.abs(R.to_x(a0, nxf, 1)) ** 2 + Gyy * k ** 2 * np.abs(R.to_x(a0, nxf)) ** 2)))
            En = np.array(En); meas = np.polyfit(En[:, 0], En[:, 1], 1)[0]
            tw = [window_terms(R, tc) for tc in (a + 75.0, a + 225.0)]
            fl = np.mean([x["flu"] for x in tw]); Wl = np.mean([x["Wlab"] for x in tw]); Pa = np.mean([x["P"] for x in tw]); E = np.mean([x["Enz"] for x in tw])
            rows.append((0.5 * (a + b), -fl, -Wl, meas, Pa, E))
            print(key, f"t {0.5*(a+b):.0f}: stress {-fl/Pa:+.2f} P, damping {-Wl/Pa:+.2f} P, sum {(-fl-Wl)/Pa:+.2f} P, measured {meas/Pa:+.2f} P", flush=True)
        out[key] = np.array(rows)
    np.savez(CACHE, **out)


def draw():
    import notes_style as ns
    C = np.load(CACHE)
    ns.apply_style()
    fig, (a, b) = ns.two_panel(sharey=True, height=2.95, bottom=0.92)
    for ax, key, tag in ((a, "nx64", "(a)"), (b, "nx128", "(b)")):
        r = C[key]; t = r[:, 0] / 1e3; P = r[:, 4]
        ax.axhline(0, c="0.75", lw=0.6)
        ax.plot(t, r[:, 1] / P, "-", c=ns.PAL["W"], lw=1.2)
        ax.plot(t, r[:, 2] / P, "-", c=ns.PAL["Z"], lw=1.2)
        ax.plot(t, (r[:, 1] + r[:, 2]) / P, "--", c="0.45", lw=1.1)
        ax.plot(t, r[:, 3] / P, "o", c="k", ms=3.5)
        ax.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
        ns.style_axes(ax); ns.tag(ax, tag)
    a.set_ylabel(r"$\dot E_{\rm nz}/P$")
    a.set_ylim(-1.6, 1.6)
    from matplotlib.lines import Line2D
    ns.legend_below(fig, [Line2D([], [], c=ns.PAL["W"], lw=1.2), Line2D([], [], c=ns.PAL["Z"], lw=1.2), Line2D([], [], c="0.45", ls="--", lw=1.1), Line2D([], [], c="k", marker="o", ls="none", ms=3.5)],
                    ["Reynolds stress", "damping (frame of the box)", "sum", "measured"], ncol=2, fontsize=ns.FS)
    ns.save(fig, "balance")


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        rebuild()
    draw()

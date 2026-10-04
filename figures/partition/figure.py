#!/usr/bin/env python3
r"""Figure `partition`: where the energy goes, from GENE's exact electrostatic-energy budget by k_y (4 Oct 2026).

  python figure.py --rebuild   read Spectral_2D_E_{fe,nonlinear,parallel,curvature}_<species> of the runs without radial
                               hyperdiffusion (rb_c0p3_hxoff, all legs; rb_c0p3_draintest/nx128_hxoff + leg 2) in windows of
                               400 time units; write cache.npz
  python figure.py             draw ../partition.pdf
All rates are divided by the Landau loss of the vortices, P_1 = -E_parallel(k_y = 1): the nonlinear supply from the jets,
-NL(k_y = 0); the part of it received by the vortices, NL(k_y = 1); the part received by their harmonics, NL(k_y >= 2);
the Landau loss of the harmonics; and the rate of change of the energy of the vortices.  (a) reference resolution,
(b) doubled radial resolution.
"""
import glob, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
CACHE = os.path.join(HERE, "cache.npz")
RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
TERMS = ("fe", "nonlinear", "parallel", "curvature")


def rebuild():
    import spec2d
    out = {}
    for key, dirs, nx, t0 in (("nx64", sorted(glob.glob(RB + "/rb_c0p3_hxoff/leg_*/out")), 64, 1000.0),
                              ("nx128", [RB + "/rb_c0p3_draintest/nx128_hxoff/out", RB + "/rb_c0p3_draintest/nx128_hxoff_leg2/out"], 128, 2600.0)):
        spec2d.NX = nx; spec2d.NKX = nx - 1; spec2d.REC = 2 + spec2d.NKX * spec2d.NKY
        T = {}
        for term in TERMS:
            ts, cs = [], []
            for d in dirs:
                if not os.path.exists(f"{d}/Spectral_2D_E_{term}_electrons.dat"):
                    continue
                c = 0
                for sp in ("electrons", "positrons"):
                    t, cc = spec2d.by_ky(f"{d}/Spectral_2D_E_{term}_{sp}.dat"); c = c + cc
                ts.append(t); cs.append(c)
            t = np.concatenate(ts); c = np.concatenate(cs); o = np.argsort(t, kind="stable"); t, c = t[o], c[o]
            keep = np.r_[True, np.diff(t) > 0]; T[term] = (t[keep], c[keep])
        rows = []
        t_all = T["fe"][0]
        for a in np.arange(t0, t_all[-1] - 799, 400.0):
            b = a + 800.0
            r = {}
            for term in TERMS:
                t, c = T[term]; m = (t >= a) & (t < b)
                r[term] = c[m].mean(0) if term != "fe" else np.array([np.polyfit(t[m], c[m, j], 1)[0] for j in range(c.shape[1])])
            P1 = -r["parallel"][1]
            rows.append((0.5 * (a + b), -r["nonlinear"][0] / P1, r["nonlinear"][1] / P1, r["nonlinear"][2:].sum() / P1,
                         -r["parallel"][2:].sum() / P1, r["fe"][1] / P1, -r["curvature"][1] / P1, P1))
        out[key] = np.array(rows); print(key, len(rows), "windows", t_all[0], t_all[-1])
    np.savez(CACHE, **out)


def draw():
    import notes_style as ns
    C = np.load(CACHE)
    ns.apply_style()
    fig, (a, b) = ns.two_panel(sharey=True, height=3.2, bottom=1.25)
    for ax, key, tag in ((a, "nx64", "(a)"), (b, "nx128", "(b)")):
        r = C[key]; t = r[:, 0] / 1e3
        ax.axhline(0, c="0.75", lw=0.6)
        ax.plot(t, r[:, 1], "-", c="k", lw=1.4)
        ax.plot(t, r[:, 2], "-", c=ns.PAL["W"], lw=1.2)
        ax.plot(t, r[:, 3], "-", c=ns.PAL["aux"], lw=1.2)
        ax.plot(t, -r[:, 4], "--", c=ns.PAL["aux"], lw=1.0)
        ax.plot(t, r[:, 5], "o", c=ns.PAL["Z"], ms=3.0)
        ax.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$"); ns.style_axes(ax); ns.tag(ax, tag)
    a.set_ylabel(r"rate $/\,P_1$"); a.set_ylim(-0.7, 2.7)
    from matplotlib.lines import Line2D
    ns.legend_below(fig, [Line2D([], [], c="k", lw=1.4), Line2D([], [], c=ns.PAL["W"], lw=1.2), Line2D([], [], c=ns.PAL["aux"], lw=1.2),
                          Line2D([], [], c=ns.PAL["aux"], ls="--", lw=1.0), Line2D([], [], c=ns.PAL["Z"], marker="o", ls="none", ms=3.0)],
                    ["supply from the jets", "received by the vortices", "received by their harmonics", "Landau loss of the harmonics", r"$\dot E_{\rm nz}$ of the vortices"],
                    ncol=2, fontsize=ns.FS)
    ns.save(fig, "partition")
    for key in ("nx64", "nx128"):
        r = C[key]
        print(key, f"supply {r[:,1].mean():.2f}  vortices {r[:,2].mean():.2f}  harmonics {r[:,3].mean():.2f}  harm. loss {r[:,4].mean():.2f}  dE_v {r[:,5].mean():.2f}  drift loss {r[:,6].mean():.3f}")


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        rebuild()
    draw()

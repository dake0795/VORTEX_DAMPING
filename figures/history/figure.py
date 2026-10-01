#!/usr/bin/env python3
r"""Figure `history`: the vortices decay and collapse while the zonal flow stays.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../history.pdf (+ .png)
  env -u LD_LIBRARY_PATH python3 figure.py --rebuild  read the run (legs 1, 3-6) and rewrite cache.npz first

(a) E_zon and E_nz as running means over W_MEAN time units (about ten breathing periods), in units of E_0, the
    largest running-mean E_zon.  ONE POINT: the vortices lose many decades and collapse; the zonal flow does not move.
(b) Gamma_E = -d ln E_nz/dt and -d ln E_zon/dt from least-squares fits of ln E against t in sliding windows of
    W_FIT time units every S_FIT, from the end of the post-burst relaxation to the end of the collapse.  ONE POINT:
    the vortex decay rate rises without a plateau as the vortices weaken; the zonal rate is orders smaller and falls.
The grey band is the collapse (traces.COLLAPSE), the same in every time-axis panel of the notes.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
CACHE = os.path.join(HERE, "cache.npz")
W_MEAN = 37.0
W_FIT, S_FIT = 400.0, 100.0
T_FIT0 = 800.0


def rebuild():
    import traces as T
    t, Ez, En, Zn = T.plunk("orig")
    Ezm, Enm = T.runmean(t, Ez, W_MEAN), T.runmean(t, En, W_MEAN)
    print(f"E_zon first exceeds E_nz at t = {t[np.argmax(Ez > En)]:.1f}; E_0 = {np.nanmax(Ezm):.4e}")
    tc, g, eg, _ = T.sliding_rate(t, En, W_FIT, S_FIT, tmin=T_FIT0, tmax=T.COLLAPSE[1])
    tz, gz, egz, _ = T.sliding_rate(t, Ez, W_FIT, S_FIT, tmin=T_FIT0, tmax=T.COLLAPSE[1])
    print("min Gamma_E %.2e at t = %.0f; last window %.2e; zonal %.2e -> %.2e" % (g.min(), tc[np.argmin(g)], g[-1], gz[1], gz[-1]))
    assert (g > 0).all() and (gz[1:] > 0).all()
    k = slice(None, None, 4)
    np.savez(CACHE, t=t[k], Ezm=Ezm[k], Enm=Enm[k], E0=np.nanmax(Ezm), tc=tc, g=g, eg=eg, tz=tz, gz=gz, egz=egz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    if ap.parse_args().rebuild or not os.path.exists(CACHE):
        rebuild()
    import notes_style as ns
    import traces as T
    from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator
    ns.apply_style()
    C = np.load(CACHE)
    fig, (a, b) = ns.two_panel()
    x = C["t"] / 1e3
    a.semilogy(x, C["Ezm"] / C["E0"], c=ns.PAL["zon"], lw=1.6)
    a.semilogy(x, C["Enm"] / C["E0"], c=ns.PAL["E"])
    T.label(a, 5.5, 0.13, r"zonal flow", c="0.35", ha="center")
    T.label(a, 5.5, 1.5e-5, r"vortices", ha="center")
    a.set_xlim(0, 20); a.set_xticks([0, 5, 10, 15, 20])
    a.set_ylim(1e-8, 1e2); a.set_yticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e2])
    a.set_ylabel(r"$E_{\rm zon}/E_0,\ E_{\rm nz}/E_0$")
    m = (C["gz"] > 0) & (C["tz"] > T_FIT0 + W_FIT)       # the first windows still see the zonal flow forming
    b.semilogy(C["tz"][m] / 1e3, C["gz"][m], c=ns.PAL["zon"], lw=1.6)
    b.semilogy(C["tc"] / 1e3, C["g"], c=ns.PAL["E"])
    T.label(b, 6.0, 2.0e-3, r"vortices", ha="center")
    T.label(b, 7.5, 2.2e-5, r"zonal flow", c="0.35", ha="center")
    b.set_xlim(0, 12); b.set_xticks([0, 4, 8, 12])
    b.set_ylim(1e-6, 1e-1)
    b.set_ylabel(r"decay rate $(c_{\rm ref}/L_{\rm ref})$")
    for ax, s in ((a, "(a)"), (b, "(b)")):
        T.collapse_band(ax)
        ax.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=12))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ns.tag(ax, s)
    a.set_yticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e2])
    ns.save(fig, "history")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""Figure `saturation`: the trapping-saturation law  omega_tr = alpha gamma_R  (1 Oct 2026).

  python figure.py      draw ../saturation.pdf

Inputs (both made from the run data by other scripts of this repository; nothing is read from the runs here):
  ../trapping/cache.npz                          omega_tr of the two vortices against time (pendulum estimate at the
                                                 critical layers; figures/trapping/figure.py --rebuild)
  ../../scripts/rayleigh/cache/growth_history.npz growth rate gamma_R of the unstable Rayleigh mode on the measured
                                                 running-mean jets against time (scripts/rayleigh/growth_history.py)
The copy of the two series used here is written to cache.npz so the figure can be redrawn without them.
(a) omega_tr (black) and alpha gamma_R (green dashed) against time, alpha = the geometric mean of their ratio;
(b) omega_tr against gamma_R with the line omega_tr = alpha gamma_R.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns
import traces

CACHE = os.path.join(HERE, "cache.npz")
T0, T1, WIN = 1000.0, 10850.0, 300.0


def rebuild():
    G = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "growth_history.npz"), allow_pickle=True)["rows"]
    T = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)["orig"]
    tg, gR = G[:, 0], 0.5 * (G[:, 1] + G[:, 3])
    tt, wtr = T[:, 0], np.sqrt(T[:, 6] * T[:, 7])
    tc = np.arange(T0, T1 + 1, 100.0)
    g = np.array([gR[(tg > c - WIN / 2) & (tg < c + WIN / 2)].mean() for c in tc])
    w = np.array([np.exp(np.log(wtr[(tt > c - WIN / 2) & (tt < c + WIN / 2)]).mean()) for c in tc])
    # the two amplitude-reduced restarts: gamma_R on their own jets (scripts/rayleigh/growth_restarts.py)
    R = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "growth_restarts.npz"))["rows"]
    Tc = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"), allow_pickle=True)
    extra = {}
    for k, run in enumerate(("a03", "a003")):
        gk = R[R[:, 0] == k]; tr = Tc[run]
        cc = np.arange(6050.0, 6651.0, 150.0)
        extra[f"g_{run}"] = np.array([0.5 * (gk[(gk[:, 1] > c - 75) & (gk[:, 1] < c + 75), 2] + gk[(gk[:, 1] > c - 75) & (gk[:, 1] < c + 75), 4]).mean() for c in cc])
        wt = np.sqrt(tr[:, 6] * tr[:, 7])
        extra[f"w_{run}"] = np.array([np.exp(np.log(wt[(tr[:, 0] > c - 75) & (tr[:, 0] < c + 75)]).mean()) for c in cc])
    np.savez(CACHE, tc=tc, gR=g, wtr=w, **extra)


def draw():
    C = np.load(CACHE)
    tc, g, w = C["tc"], C["gR"], C["wtr"]
    alpha = np.exp(np.mean(np.log(w / g)))
    p = np.polyfit(np.log(g), np.log(w), 1)
    print(f"alpha = {alpha:.2f}; ratio range {np.min(w / g):.2f} .. {np.max(w / g):.2f}; log-log slope {p[0]:.2f}")
    ns.apply_style()
    fig, (a, b) = ns.two_panel()
    traces.collapse_band(a)
    a.semilogy(tc / 1e3, w, c="k")
    a.semilogy(tc / 1e3, alpha * g, c=ns.PAL["theory"], ls="--")
    a.set_xlim(0, 11.1); a.set_ylim(0.03, 1.5)
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$"); a.set_ylabel(r"frequency $(c_{\rm ref}/L_{\rm ref})$")
    traces.label(a, 5.2, 0.62, r"$\omega_{\rm tr}$", c="k"); traces.label(a, 2.0, 0.27, r"$C\gamma_{\rm R}$", c=ns.PAL["theory"])
    gg = np.array([4e-3, 0.12])
    b.loglog(gg, alpha * gg, c=ns.PAL["theory"], ls="--")
    b.loglog(g, w, "o", c="k", ms=3.2)
    for run, mk in (("a03", "s"), ("a003", "^")):
        b.loglog(C[f"g_{run}"], C[f"w_{run}"], mk, c=ns.PAL[run], mfc="w", ms=5, mew=1.0, zorder=5)
        print(run, "ratio", C[f"w_{run}"] / C[f"g_{run}"])
    b.set_xlim(4e-3, 0.12); b.set_ylim(0.03, 1.5)
    b.set_xlabel(r"$\gamma_{\rm R}\;(c_{\rm ref}/L_{\rm ref})$"); b.set_ylabel(r"$\omega_{\rm tr}\;(c_{\rm ref}/L_{\rm ref})$")
    traces.label(b, 0.0075, 0.2, r"$\omega_{\rm tr}=C\gamma_{\rm R}$", c=ns.PAL["theory"])
    ns.tag(a, "(a)", y=0.14, va="bottom"); ns.tag(b, "(b)")
    ns.save(fig, "saturation")


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        rebuild()
    draw()

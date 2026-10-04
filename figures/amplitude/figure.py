#!/usr/bin/env python3
r"""Figure `amplitude`: the amplitude-reduced restarts against the original run.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../amplitude.pdf (+ .png)
  env -u LD_LIBRARY_PATH python3 figure.py --rebuild  read the three runs and rewrite cache.npz first

(a) E_nz/(eps^2 E_nz^orig(t0)) against t - t0, t0 the restart time, as running means over W_MEAN time units (about
    ten breathing periods; for a restart the mean is taken across t0 on the composite series whose part before t0
    is the original, of which the restart is the rescaled copy).  E_nz^orig(t0) is the running mean of the original run at t0, so the original starts at
    one.  ONE POINT: a linear process on a frozen zonal state would keep the three curves together; they separate.
(b) Gamma_E = -d ln E_nz/dt from least-squares fits of ln E_nz against the amplitude measure
    a = (E_nz/E_zon)^{1/2} (geometric mean over the window): original run in windows of W_ORIG time units every S_ORIG
    from T_FIT0 (the end of the post-burst relaxation, where Gamma_E is smallest; see figure `history`) to the end of
    the collapse; restarts in windows of W_FIT every S_FIT after their initial transient (t > T_TRANS).  Guide
    lines: Gamma_E proportional to a^{-5/8}, a^{-1}, a^{-5/2}, all three through the original run at T_ANCHOR.  ONE POINT: the weaker the
    vortices, the faster they decay, and the three runs lie close to one curve.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
CACHE = os.path.join(HERE, "cache.npz")
W_MEAN = 37.0
W_FIT, S_FIT = 200.0, 100.0
W_ORIG, S_ORIG = 400.0, 200.0
T_FIT0, T_COLLAPSE = 3600.0, 11100.0
T_TRANS = 5980.0
RUNS = ("orig", "a03", "a003")


def rebuild():
    import traces as T
    out = {}
    to, Ezo, Eno, _ = T.plunk("orig")
    i0 = np.searchsorted(to, T.T_RESTART)            # the last record of leg 3 = the checkpoint
    assert abs(to[i0] - T.T_RESTART) < 1e-6
    E00 = T.runmean(to, Eno, W_MEAN)[i0]             # breathing-averaged E_nz of the original at the restart
    print(f"E_nz^orig(t0): checkpoint record {Eno[i0]:.4e}, running mean {E00:.4e}")
    for r in RUNS:
        t, Ez, En, Zn = T.plunk(r)
        m = (t >= T.T_RESTART - 1e-6) & (t <= 6707.5)
        y = En / (T.EPS[r] ** 2 * E00)
        # the state of a restart before t0 is the original with its k_y >= 1 part scaled by eps, so its rescaled
        # energy there is that of the original: the running mean is taken across t0 on this composite series
        pre = to < T.T_RESTART - 1e-6
        tcmp, ycmp = (t, y) if r == "orig" else (np.concatenate([to[pre], t]), np.concatenate([Eno[pre] / E00, y]))
        ymean = T.runmean(tcmp, ycmp, W_MEAN)[len(tcmp) - len(t):] if r != "orig" else T.runmean(t, y, W_MEAN)
        out[f"t_{r}"], out[f"y_{r}"], out[f"ym_{r}"] = t[m] - T.T_RESTART, y[m], ymean[m]
        t0, t1 = (T_FIT0, T_COLLAPSE) if r == "orig" else (T_TRANS, t[-1])
        w, st = (W_ORIG, S_ORIG) if r == "orig" else (W_FIT, S_FIT)
        tc, g, eg, en = T.sliding_rate(t, En, w, st, tmin=t0, tmax=t1)
        _, _, _, ez = T.sliding_rate(t, Ez, w, st, tmin=t0, tmax=t1)
        out[f"tc_{r}"], out[f"g_{r}"], out[f"eg_{r}"], out[f"a_{r}"] = tc, g, eg, np.sqrt(en / ez)
        if r != "orig":
            k = np.argmax(y[m]); print(f"{r}: first sample {y[m][0]:.3f}; peak of rescaled E_nz {y[m][k]:.2f} at t = {t[m][k]:.0f}")
            p = np.polyfit(t[t > T_TRANS], np.log(En[t > T_TRANS]), 1); print(f"   Gamma_E over t > {T_TRANS:.0f}: {-p[0]:.3e}")
        p = np.polyfit(t[(t >= 5908) & (t <= 6707)], np.log(En[(t >= 5908) & (t <= 6707)]), 1)
        print(f"{r}: Gamma_E over 5908-6707: {-p[0]:.3e}")
        for a_, b_, c_, d_ in zip(tc, g, eg, out[f"a_{r}"]):
            if True:
                print(f"   {a_:7.0f}  Gamma {b_:+.2e} +- {c_:.1e}   a {d_:.3e}")
    mo = (out["tc_orig"] > 5900) & (out["tc_orig"] < 10400)
    for lab, mm in (("5900-10400", mo), ("5900-8300", mo & (out["tc_orig"] < 8300)), ("8300-10400", mo & (out["tc_orig"] >= 8300))):
        print(f"orig {lab}: Gamma ~ a^{np.polyfit(np.log(out['a_orig'][mm]), np.log(out['g_orig'][mm]), 1)[0]:.2f}")
    np.savez(CACHE, **out)


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
    lab = {"orig": r"original", "a03": r"$\epsilon=0.3$", "a003": r"$\epsilon=0.03$"}
    ls = {"orig": "-", "a03": "--", "a003": "-."}
    mk = {"orig": "o", "a03": "s", "a003": "^"}
    fig, (a, b) = ns.two_panel()
    for r in RUNS:
        a.semilogy(C[f"t_{r}"], C[f"ym_{r}"], c=ns.PAL[r], ls=ls[r], lw=1.4)
    T.label(a, 330, 9.5, lab["a003"], c=ns.PAL["a003"], ha="left", va="bottom")
    T.label(a, 330, 1.45, lab["a03"], c=ns.PAL["a03"], ha="left", va="bottom")
    T.label(a, 60, 0.42, lab["orig"], c=ns.PAL["orig"], ha="left", va="center")
    a.set_xlim(0, 800); a.set_ylim(0.15, 40)
    a.xaxis.set_major_locator(MultipleLocator(200)); a.xaxis.set_minor_locator(MultipleLocator(50))
    a.set_xlabel(r"$t-t_0\;(L_{\rm ref}/c_{\rm ref})$")
    a.set_ylabel(r"$E_{\rm nz}/\epsilon^2E_{\rm nz}^{\rm orig}(t_0)$")
    b.loglog(C["a_orig"], C["g_orig"], c=ns.PAL["orig"], lw=1.4, zorder=3)
    for r in ("a03", "a003"):
        m = C[f"g_{r}"] > 0
        b.loglog(C[f"a_{r}"][m], C[f"g_{r}"][m], ls="none", c=ns.PAL[r], marker=mk[r], ms=4.2, mfc="w", mew=0.9, zorder=4)
    b.set_xlim(0.2, 5e-4); b.set_ylim(1.5e-4, 3e-2)
    b.set_xlabel(r"$a=(E_{\rm nz}/E_{\rm zon})^{1/2}$")
    b.set_ylabel(r"$\Gamma_E\;(c_{\rm ref}/L_{\rm ref})$")
    for ax in (a, b):
        ax.yaxis.set_minor_formatter(NullFormatter())
    b.xaxis.set_major_locator(LogLocator(base=10, numticks=6)); b.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    fig.delaxes(b); pa = a.get_position(); a.set_position([0.5 - pa.width / 2, pa.y0, pa.width, pa.height])  # panel (b) removed 4 Oct 2026: it read the numerical collapse as physics
    ns.save(fig, "amplitude")


if __name__ == "__main__":
    main()

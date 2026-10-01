#!/usr/bin/env python3
r"""Figure `remnant`: the non-zonal electrostatic energy after the collapse of the vortices.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../remnant.pdf (+ .png)
  env -u LD_LIBRARY_PATH python3 figure.py --rebuild  read the run (legs 4-6) and rewrite cache.npz first

(a) E_nz/E_0 for t >= T_SHOW (E_0 as in figure `history`): block means over W_BLOCK time units (dots; geometric
    means); an exponential fitted to the block means of the
    remnant, t >= T_REM; and the passive drift-mixing law (t - t_c)^{-3} measured from the collapse (t_c = T_C) and
    anchored at the first block of the remnant.  ONE POINT: what is left after the collapse decays slowly and
    irregularly, far more slowly than the algebraic tail of passive phase mixing started at the collapse.
(b) the block means of the remnant against t - t0 on logarithmic axes, for t0 = t_c (the collapse) and t0 = t_b
    (the burst, T_BURST), with a slope of -3 through the first block of each.  ONE POINT: a t^{-3} law fits only if
    time is counted from the burst, where over this window it cannot be told from a slow exponential.
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
W_BLOCK = 250.0
T_SHOW, T_C, T_REM, T_BURST = 9000.0, 10900.0, 11100.0, 665.0


def rebuild():
    import traces as T
    t, Ez, En, Zn = T.plunk("orig")
    E0 = np.nanmax(T.runmean(t, Ez, W_MEAN))
    tb, eb = [], []
    a = T_SHOW + (T_REM - T_SHOW) % W_BLOCK          # block edges aligned on T_REM
    while a + W_BLOCK <= t[-1]:
        m = (t >= a) & (t < a + W_BLOCK)
        tb.append(a + W_BLOCK / 2); eb.append(np.exp(np.mean(np.log(En[m]))))
        a += W_BLOCK
    tb, eb = np.array(tb), np.array(eb)
    r = tb > T_REM
    x, y = tb[r], np.log(eb[r])
    pe = np.polyfit(x, y, 1)
    rms = lambda f: np.sqrt(np.mean((y - f) ** 2))
    print(f"remnant blocks: {r.sum()}, t = {x[0]:.0f} .. {x[-1]:.0f}; E_nz first {eb[r][0]:.0f}, max {eb[r].max():.0f} at "
          f"{x[np.argmax(eb[r])]:.0f}, last {eb[r][-1]:.0f}; drop in ln E_nz {y[0] - y[-1]:.2f}")
    print(f"exponential: Gamma_E = {-pe[0]:.3e}, rms log residual {rms(np.polyval(pe, x)):.3f}")
    h = len(x) // 2
    print("   halves:", " ".join(f"{-np.polyfit(x[s], y[s], 1)[0]:.2e}" for s in (slice(0, h), slice(h, None))))
    for t0 in (0.0, T_BURST, T_C):
        pp = np.polyfit(np.log(x - t0), y, 1)
        print(f"power law from t0 = {t0:.0f}: exponent {pp[0]:+.2f}, rms log residual {rms(np.polyval(pp, np.log(x - t0))):.3f}")
    print(f"(t - T_C)^-3 anchored at the first block predicts a drop in ln E_nz of {3 * np.log((x[-1] - T_C) / (x[0] - T_C)):.1f}")
    np.savez(CACHE, tend=t[-1], E0=E0, tb=tb, eb=eb, pe=pe)


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
    E0, tb, eb = C["E0"], C["tb"], C["eb"] / C["E0"]
    r = tb > T_REM
    fig, (a, b) = ns.two_panel()
    a.semilogy(tb / 1e3, eb, ls="none", marker="o", ms=3.0, c="k", zorder=4)
    tt = np.linspace(T_REM, float(C["tend"]), 200)
    a.semilogy(tt / 1e3, np.exp(np.polyval(C["pe"], tt)) / E0, c=ns.PAL["W"], lw=1.4, zorder=3)
    tp = np.linspace(tb[r][0], float(C["tend"]), 300)
    a.semilogy(tp / 1e3, eb[r][0] * ((tp - T_C) / (tb[r][0] - T_C)) ** -3.0, c=ns.PAL["Z"], ls="--", lw=1.4, zorder=3)
    T.label(a, 15.6, 6.0e-6, r"exponential", c=ns.PAL["W"], ha="center")
    T.label(a, 13.6, 4.0e-9, r"$\propto(t-t_{\rm c})^{-3}$", c=ns.PAL["Z"], ha="left")
    T.collapse_band(a)
    a.set_xlim(9, 20); a.xaxis.set_major_locator(MultipleLocator(5)); a.xaxis.set_minor_locator(MultipleLocator(1))
    a.set_ylim(1e-9, 1e-2)
    a.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=12)); a.yaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
    a.set_ylabel(r"$E_{\rm nz}/E_0$")
    for t0, c, mk in ((T_C, "k", "o"), (T_BURST, ns.PAL["aux"], "s")):
        x = tb[r] - t0
        b.loglog(x, eb[r], ls="none", marker=mk, ms=3.4, c=c, mfc=c if mk == "o" else "w", mew=0.8)
        xg = np.array([x[0], x[0] * 2.2])
        b.loglog(xg, eb[r][0] * (xg / x[0]) ** -3.0, c="0.45", ls="--", lw=0.8, zorder=1)
    T.label(b, 1.9e3, 2.3e-6, r"$t_0=t_{\rm c}$", c="k", ha="center")
    T.label(b, 1.5e4, 2.3e-6, r"$t_0=t_{\rm b}$", c=ns.PAL["aux"], ha="center")
    b.text(1.0e3, 3.0e-8, r"$\propto(t-t_0)^{-3}$", color="0.3", ha="left", va="center")
    b.set_xlim(1.5e2, 4e4); b.set_ylim(1e-8, 1e-5)
    b.set_xlabel(r"$t-t_0\;(L_{\rm ref}/c_{\rm ref})$")
    b.set_ylabel(r"$E_{\rm nz}/E_0$")
    b.yaxis.set_minor_formatter(NullFormatter())
    ns.tag(a, "(a)", x=0.96, ha="right"); ns.tag(b, "(b)")
    ns.save(fig, "remnant")


if __name__ == "__main__":
    main()

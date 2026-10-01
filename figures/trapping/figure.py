#!/usr/bin/env python3
r"""Figure `trapping` of the vortex-damping notes: trapping frequency against damping rate, and damping rate against
the width of the cat's eye.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../trapping.pdf (+ .png)
  python3 figure.py --rebuild                         read the field data and rewrite cache.npz first

Data (read only): mid-plane potential of the reduced-box run rb_c0p3_hxy (legs 1, 3-6) and of its two
amplitude-reduced restarts (rb_c0p3_hxy_ampltest/amp0p3, amp0p03).  Definitions in ../lib/vortex_field.py:
two-frequency fit of the k_y = 1 potential in 60-time-unit windows (stepped by 20), which separates the two
counter-propagating vortices; critical layers U = c_pm of the window-mean zonal flow; pendulum trapping frequency
w_tr = k_y sqrt(|S| A) and separatrix half-width w = 2 sqrt(A/|S|) there (mean of the two layers of each vortex);
per-vortex amplitude damping rate gamma_amp = -(1/2) d ln E_pm / dt from linear fits of ln E_pm over 200 time units
[points of panel (b)] or 400 time units [sliding, the curve of panel (a)],
E_pm = sum_kx |a_pm(k_x)|^2; breathing frequency w_b = w_+ - w_-.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
CACHE = os.path.join(HERE, "cache.npz")
WIN, STEP, FIT = 60.0, 20.0, 200.0
FIT_A = 400.0                            # the curve of panel (a) uses longer fits (smoother); the points of (b) FIT
T_BURST, T_END = 700.0, 11000.0          # windows of the original run
T_TRANS = 5980.0                         # end of the initial transient of the restarts
RES_MAX = 0.10                           # two-vortex fit accepted below this relative residual
T_EARLY = 5000.0                         # panel (b): the original run from here on (earlier, the rate still fluctuates)
T_COLLAPSE = (10600.0, 11100.0)          # shaded in every time-axis panel of the notes


def series(run, t0, t1):
    import vortex_field as vf
    t, p = vf.midplane_series(run, t0 - WIN, t1 + WIN)
    rows, w0 = [], None
    for tc in np.arange(t0 + WIN / 2, t1 - WIN / 2 + 1e-6, STEP):
        m = (t >= tc - WIN / 2) & (t <= tc + WIN / 2)
        best = None
        for g in ([w0] if w0 is not None else []) + [np.array([q, -q]) for q in (0.8, 0.86, 1.0, 1.2, 1.35)]:
            o = vf.window_analysis(t[m], p[m], g)
            if best is None or o["res"] < best["res"]:
                best = o
            if w0 is not None and best["res"] < 0.02:
                break
        o = best
        w0 = o["w"] if o["res"] < 0.05 else None
        rows.append([o["t"], o["w"][0], o["w"][1], o["res"], o["E"][0], o["E"][1], o["wtr_p"], o["wtr_m"],
                     o["half_p"], o["half_m"], o["Ac_p"], o["Ac_m"], o["S_p"], o["S_m"], o["Umax"], o["Umin"]])
    return np.array(rows)


def rates(r, centres, fit=FIT):
    """gamma_amp of the two vortices and the mean half-widths from fits of length `fit` centred on `centres`."""
    out = []
    for tc in centres:
        m = (r[:, 0] >= tc - fit / 2) & (r[:, 0] <= tc + fit / 2)
        if m.sum() < 5:
            continue
        g = [-0.5 * np.polyfit(r[m, 0], np.log(r[m, 4 + k]), 1)[0] for k in (0, 1)]
        out.append([tc, g[0], g[1], np.exp(np.log(r[m, 8]).mean()), np.exp(np.log(r[m, 9]).mean()), r[m, 3].max()])
    return np.array(out)


def rebuild():
    c = {}
    for run, (t0, t1) in (("orig", (T_BURST, T_END)), ("a03", (5906.7, 6707.0)), ("a003", (5906.7, 6707.0))):
        r = series(run, t0, t1)
        lo = r[0, 0] + FIT / 2 if run == "orig" else T_TRANS + FIT / 2
        c[run] = r
        c[run + "_sl"] = rates(r, np.arange(lo + (FIT_A - FIT) / 2, r[-1, 0] - FIT_A / 2 + 1e-6, STEP), FIT_A)   # sliding, for the curve of (a)
        c[run + "_bl"] = rates(r, np.arange(lo, r[-1, 0] - FIT / 2 + 1e-6, FIT if run == "orig" else FIT / 2))
        print(run, r.shape, c[run + "_sl"].shape, c[run + "_bl"].shape)
    np.savez_compressed(CACHE, cols="t w+ w- res E+ E- wtr+ wtr- half+ half- A+ A- S+ S- Umax Umin", **c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    if ap.parse_args().rebuild:
        rebuild()
    import notes_style as ns
    import vortex_field as vf
    ns.apply_style()
    PAL, FS = ns.PAL, ns.FS
    C = np.load(CACHE)
    o, sl = C["orig"], C["orig_sl"]
    ok = o[:, 3] < RES_MAX
    fig, (ax, bx) = ns.two_panel(left=0.66, wspace=0.80)

    # (a) the one point: the trapping frequency exceeds the damping rate throughout.  One line per quantity (the
    # geometric mean of the two vortices) on a band spanning the two.
    ax.axvspan(*T_COLLAPSE, color="0.9", lw=0)
    ok = o[:, 3] < RES_MAX
    ax.fill_between(o[ok, 0], o[ok, 6:8].min(axis=1), o[ok, 6:8].max(axis=1), color="k", alpha=0.2, lw=0)
    ax.semilogy(o[ok, 0], np.sqrt(o[ok, 6] * o[ok, 7]), c="k")
    g = np.where((sl[:, 1:3] > 0) & (sl[:, 5:6] < RES_MAX), sl[:, 1:3], np.nan)
    ax.fill_between(sl[:, 0], g.min(axis=1), g.max(axis=1), color=PAL["theory"], alpha=0.25, lw=0)
    ax.semilogy(sl[:, 0], np.sqrt(g[:, 0] * g[:, 1]), c=PAL["theory"])
    ax.text(7300, 0.50, r"$\omega_{\rm tr}$", ha="center", va="bottom")
    ax.text(7300, 4.6e-4, r"$\gamma_{\rm amp}$", ha="center", va="bottom", color=PAL["theory"])
    ax.set_xlabel(r"$t\;(L_{\rm ref}/c_{\rm ref})$")
    ax.set_ylabel(r"frequency, rate $(c_{\rm ref}/L_{\rm ref})$")
    ax.set_xlim(0, 11500)
    ax.set_xticks([0, 5000, 10000])
    ax.set_ylim(5e-5, 8.0)
    ns.tag(ax, "(a)")

    # (b) the one point: the damping rate rises steeply as the cat's eye shrinks towards one grid cell, which is
    # also about a Debye length.  One point per vortex and per fit window; the original run from T_EARLY on.
    for run, mk, ms, lab in (("orig", "o", 3.0, "original"), ("a03", "s", 3.6, r"$\epsilon=0.3$"),
                             ("a003", "^", 3.8, r"$\epsilon=0.03$")):
        b = C[run + "_bl"]
        b = b[(b[:, 5] < RES_MAX) & (b[:, 0] >= T_EARLY)]
        w_, g_ = b[:, 3:5].ravel() / vf.DX, b[:, 1:3].ravel()
        bx.loglog(w_[g_ > 0], g_[g_ > 0], mk, ms=ms, c=PAL[run], mew=0, label=lab)
    bx.axvline(1.0, c="k", ls="--", lw=0.8)
    bx.axvline(vf.LAMBDA_D / vf.DX, c="k", ls=":", lw=0.9)
    bx.text(0.93, 1.25e-4, r"$\Delta x$", ha="right", va="bottom")
    bx.text(vf.LAMBDA_D / vf.DX * 1.07, 1.25e-4, r"$\lambda_D$", ha="left", va="bottom")
    bx.set_xlabel(r"$w/\Delta x$")
    bx.set_ylabel(r"$\gamma_{\rm amp}\;(c_{\rm ref}/L_{\rm ref})$")
    bx.set_xlim(0.4, 5.0)
    bx.set_xticks([0.5, 1, 2, 4])
    bx.set_xticks([], minor=True)
    bx.set_xticklabels(["0.5", "1", "2", "4"])
    bx.set_ylim(1e-4, 4e-2)
    bx.legend(loc="upper right", handletextpad=0.1, borderaxespad=0.2)
    ns.tag(bx, "(b)")
    ns.save(fig, "trapping")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""Figure `entropy`: electrostatic energy, entropy and free energy of the non-zonal part.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../entropy.pdf (+ .png)
  env -u LD_LIBRARY_PATH python3 figure.py --rebuild  read the three runs and rewrite cache.npz first

(a) E_nz and Z_nz of the original run (their sum is the non-zonal free energy W_nz) as running means over W_MEAN
    time units (about ten breathing periods), in units of E_0, the largest running-mean E_zon (as in figure
    `history`).  ONE POINT: the entropy of the vortices falls with their electrostatic energy, so W_nz is not
    conserved; only in the collapse (grey band, traces.COLLAPSE) does the entropy stay behind.
(b) Z_nz/E_nz (ratio of the running means) against the amplitude measure a = (E_nz/E_zon)^{1/2} (running means):
    original run from T0 to the end of the collapse, and the two restarts after their initial transient
    (t > T_TRANS).  ONE POINT: at the same amplitude the restarts carry less entropy than the original - the state
    is not a function of the amplitude alone.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
CACHE = os.path.join(HERE, "cache.npz")
W_MEAN = 37.0
T0 = 1000.0
T_TRANS = 5980.0
T_REMNANT = 11100.0                       # end of the collapse: the remnant after it is not drawn in (b)
RUNS = ("orig", "a03", "a003")


def rebuild():
    import traces as T
    out = {}
    for r in RUNS:
        t, Ez, En, Zn = T.plunk(r)
        Ezm, Enm, Znm = (T.runmean(t, y, W_MEAN) for y in (Ez, En, Zn))
        if r == "orig":
            k = slice(None, None, 4)
            out.update(t=t[k], Enm=Enm[k], Znm=Znm[k], E0=np.nanmax(Ezm))
            for a_, b_ in ((1000, 3000), (5000, 7000), (9000, 9500), (10000, 10300), (10300, 10600), (10600, 10900),
                           (10900, 11100), (11100, 19700)):
                m = (t >= a_) & (t < b_)
                rr = [-np.polyfit(t[m], np.log(y[m]), 1)[0] for y in (En, Zn, En + Zn)]
                print(f"orig {a_:6d}-{b_:6d}: Z/E {np.mean(Zn[m]) / np.mean(En[m]):6.2f}   rates E, Z, W (1e-3): "
                      + " ".join(f"{1e3 * x:6.2f}" for x in rr))
        m = np.isfinite(Enm) & (t > (T0 if r == "orig" else T_TRANS))
        k = slice(None, None, 4 if r == "orig" else 1)
        out[f"a_{r}"], out[f"q_{r}"], out[f"tq_{r}"] = np.sqrt(Enm[m] / Ezm[m])[k], (Znm[m] / Enm[m])[k], t[m][k]
        if r != "orig":
            mm = t > T_TRANS
            rr = [-np.polyfit(t[mm], np.log(y[mm]), 1)[0] for y in (En, Zn, En + Zn)]
            print(f"{r}: Z/E from {out[f'q_{r}'][0]:.2f} to {out[f'q_{r}'][-1]:.2f}; rates E, Z, W (1e-3): "
                  + " ".join(f"{1e3 * x:6.2f}" for x in rr))
    for r in ("a03", "a003"):                        # the original at the amplitudes the restart passes through
        ao, qo = out["a_orig"], out["q_orig"]
        sel = out["tq_orig"] > 5000
        for a_ in (out[f"a_{r}"][0], out[f"a_{r}"][-1]):
            j = np.argmin(np.abs(np.log(ao[sel] / a_)))
            print(f"   {r} at a = {a_:.2e}: original has Z/E = {qo[sel][j]:.2f} (t = {out['tq_orig'][sel][j]:.0f})")
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
    fig, (a, b) = ns.two_panel()
    x = C["t"] / 1e3
    a.semilogy(x, C["Znm"] / C["E0"], c=ns.PAL["Z"], ls="--", lw=1.4)
    a.semilogy(x, C["Enm"] / C["E0"], c=ns.PAL["E"], lw=1.2)
    T.label(a, 15.8, 1.6e-4, r"entropy, $Z_{\rm nz}$", c=ns.PAL["Z"], ha="center")
    T.label(a, 15.2, 2.2e-8, r"energy, $E_{\rm nz}$", c=ns.PAL["E"], ha="center")
    T.collapse_band(a)
    a.set_xlim(0, 20); a.set_xticks([0, 5, 10, 15, 20]); a.xaxis.set_minor_locator(MultipleLocator(1))
    a.set_ylim(1e-9, 1e2); a.set_yticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=12)); a.yaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
    a.set_ylabel(r"$E_{\rm nz}/E_0,\ Z_{\rm nz}/E_0$")
    ls = {"orig": "-", "a03": "--", "a003": "-."}
    for r in RUNS:
        m = C[f"tq_{r}"] <= T_REMNANT
        b.loglog(C[f"a_{r}"][m], C[f"q_{r}"][m], c=ns.PAL[r], ls=ls[r], lw=1.2 if r == "orig" else 1.6,
                 zorder=3 if r != "orig" else 2)
    T.label(b, 9.0e-2, 1.5, r"original", c=ns.PAL["orig"], ha="left", va="bottom")
    T.label(b, 2.6e-2, 0.36, r"$\epsilon=0.3$", c=ns.PAL["a03"], ha="left", va="center")
    T.label(b, 2.6e-3, 0.36, r"$\epsilon=0.03$", c=ns.PAL["a003"], ha="left", va="center")
    b.set_xlim(0.3, 3e-4); b.set_ylim(0.2, 20)
    b.set_xlabel(r"$a=(E_{\rm nz}/E_{\rm zon})^{1/2}$")
    b.set_ylabel(r"$Z_{\rm nz}/E_{\rm nz}$")
    b.xaxis.set_major_locator(LogLocator(base=10, numticks=6)); b.yaxis.set_minor_formatter(NullFormatter())
    ns.tag(a, "(a)"); ns.tag(b, "(b)")
    ns.save(fig, "entropy")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""Figure `channels`: the channels through which the vortices lose electrostatic energy, against the drift-resonance
prediction (1 Oct 2026).

  python figure.py --rebuild   read Spectral_2D_E_{fe,parallel,curvature}_<species>.dat of the reduced-box run
                               (legs 1, 3-5) and of the two amplitude-reduced restarts; write cache.npz
  python figure.py             draw ../channels.pdf from cache.npz

For the non-zonal part (k_y >= 1, species summed) the script stores, in sliding windows, the mean native electrostatic
energy E and the mean rates of change of E through GENE's parallel-streaming and curvature (magnetic-drift) terms.
PROD holds the same two rates for the k_y = 1
vortices of the production run, from ../measurements/production_loss_rate.py.
(a) the two rates of the original run against time; (b) against the amplitude a, with the
restarts and the production run.
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, os.path.join(HERE, "..", "budget"))
import notes_style as ns
import traces
import spec2d

CACHE = os.path.join(HERE, "cache.npz")
PROD = dict(a=0.25, par=1.70e-4, curv=2.04e-4)            # production_loss_rate.py (t = 741-936 of leg 2)
SPEC = {"orig": [f"{ns.RB}/rb_c0p3_hxy/leg_000{l}/out" for l in (1, 3, 4, 5)],
        "a03": [f"{ns.RB}/rb_c0p3_hxy_ampltest/amp0p3/out"], "a003": [f"{ns.RB}/rb_c0p3_hxy_ampltest/amp0p03/out"]}
WINDOWS = {"orig": (400.0, 100.0, 1200.0, 10900.0), "a03": (240.0, 120.0, 6100.0, 6587.0), "a003": (240.0, 120.0, 6100.0, 6587.0)}


def rebuild():
    out = {}
    for run, dirs in SPEC.items():
        S = {k: [] for k in ("fe", "parallel", "curvature")}; T = []
        for d in dirs:
            for k in S:
                c = 0
                for sp in ("electrons", "positrons"):
                    t, cc = spec2d.by_ky(f"{d}/Spectral_2D_E_{k}_{sp}.dat"); c = c + cc[:, 1:].sum(axis=1)
                S[k].append(c)
            T.append(t)
        t = np.concatenate(T); o = np.argsort(t); t = t[o]
        S = {k: np.concatenate(v)[o] for k, v in S.items()}
        tp, ez, en, zn = traces.plunk(run)
        wid, step, c0, c1 = WINDOWS[run]
        tc = np.arange(c0, c1 + 1, step); r = {k: [] for k in ("E", "par", "curv", "a")}
        for c in tc:
            m = (t >= c - wid / 2) & (t < c + wid / 2); mp = (tp >= c - wid / 2) & (tp < c + wid / 2)
            E = S["fe"][m].mean()
            r["E"].append(E); r["par"].append(-S["parallel"][m].mean() / E); r["curv"].append(-S["curvature"][m].mean() / E)
            r["a"].append(np.sqrt(en[mp].mean() / ez[mp].mean()))
        out[f"tc_{run}"] = tc
        for k, v in r.items():
            out[f"{k}_{run}"] = np.array(v)
        print(run, len(tc), "windows; parallel", np.array(r["par"])[::max(1, len(tc) // 6)], "curvature", np.array(r["curv"])[::max(1, len(tc) // 6)])
    np.savez(CACHE, **out)


def draw():
    C = np.load(CACHE)
    ns.apply_style()
    fig, (a, b) = ns.two_panel(sharey=True)
    tc = C["tc_orig"] / 1e3
    traces.collapse_band(a)
    pos = lambda y: np.where(y > 0, y, np.nan)
    a.semilogy(tc, pos(C["par_orig"]), c="k", ls="-")
    a.semilogy(tc, pos(C["curv_orig"]), c="k", ls="--")
    a.set_xlim(1, 11.1); a.set_ylim(2e-5, 1e-2)
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$"); a.set_ylabel(r"loss rate $(c_{\rm ref}/L_{\rm ref})$")
    traces.label(a, 4.2, 3.1e-3, "parallel streaming", c="k"); traces.label(a, 4.6, 2.6e-4, "magnetic drift", c="k")
    am = C["a_orig"]
    b.loglog(am, pos(C["par_orig"]), c="k", ls="-"); b.loglog(am, pos(C["curv_orig"]), c="k", ls="--")
    for run, mk in (("a03", "s"), ("a003", "^")):
        b.loglog(C[f"a_{run}"], pos(C[f"par_{run}"]), mk, c=ns.PAL[run], mfc=ns.PAL[run], ms=4.2)
        b.loglog(C[f"a_{run}"], pos(C[f"curv_{run}"]), mk, c=ns.PAL[run], mfc="w", ms=4.2, mew=0.9)
    b.loglog([PROD["a"]], [PROD["par"]], "*", c=ns.PAL["aux"], ms=9); b.loglog([PROD["a"]], [PROD["curv"]], "*", c=ns.PAL["aux"], mfc="w", ms=9, mew=0.9)
    b.set_xlim(0.5, 1.5e-3)
    b.set_xlabel(r"$a=(E_{\rm nz}/E_{\rm zon})^{1/2}$")
    __import__("matplotlib.pyplot").pyplot.setp(b.get_yticklabels(), visible=False)
    ns.tag(a, "(a)", y=0.14, va="bottom"); ns.tag(b, "(b)", y=0.14, va="bottom")
    ns.save(fig, "channels")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rebuild", action="store_true"); A = ap.parse_args()
    if A.rebuild or not os.path.exists(CACHE):
        rebuild()
    draw()

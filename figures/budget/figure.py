#!/usr/bin/env python3
r"""Figures `budget` and `dissipation`: the free-energy budget of the vortices, mode by mode (1 Oct 2026).

  python figure.py --rebuild   read Spectral_2D_<term>.dat of the reduced-box run (legs 1, 3-5) and of the two
                               amplitude-reduced restarts, and plunk_e_time.dat for the amplitude; write cache.npz
  python figure.py             draw ../budget.pdf and ../dissipation.pdf from cache.npz

GENE's free-energy diagnostic (diag_fespec2d) resolves every term of the free-energy balance in (k_x, k_y).  Summed
over k_x and over k_y >= 1 it gives the budget of the non-zonal part alone,
    dW_nz/dt = T - S,     S = -(hyp_v + hyp_z + hyp_kperp)_nz >= 0,
where T is the nonlinear transfer from k_y = 0 (the gradient drive of the non-zonal part is negligible and is counted
in T).  In sliding windows the script stores the window means of W_nz and of each sink, the slope of W_nz(t), and
the transfer as the residual T = dW_nz/dt + S.  (The direct time-mean of the sampled nonlinear term is also stored,
as Tdir; it has the same sign and size but is aliased by the breathing, which the output samples about three times
per period.)

budget:      (a) S/W_nz, T/W_nz and -dln W_nz/dt of the original run against time;
             (b) S/W_nz and T/W_nz against the amplitude a = (E_nz/E_zon)^{1/2} for the original run and the restarts.
dissipation: (a) the perpendicular sink, (b) the velocity-space sink, each split into the part acting on the zonal
             (k_y = 0) and on the non-zonal (k_y >= 1) part of the distribution, against time.
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import notes_style as ns
import traces
import spec2d

CACHE = os.path.join(HERE, "cache.npz")
TERMS = ("fe", "hyp_v", "hyp_z", "hyp_kperp", "nonlinear")
SPEC = {"orig": [f"{ns.RB}/rb_c0p3_hxy/leg_000{l}/out" for l in (1, 3, 4, 5)],
        "a03": [f"{ns.RB}/rb_c0p3_hxy_ampltest/amp0p3/out"], "a003": [f"{ns.RB}/rb_c0p3_hxy_ampltest/amp0p03/out"]}
# windows: (width, step, first centre, last centre)
WINDOWS = {"orig": (400.0, 100.0, 1200.0, 15900.0), "a03": (240.0, 120.0, 6100.0, 6587.0), "a003": (240.0, 120.0, 6100.0, 6587.0)}


def rebuild():
    out = {}
    for run, dirs in SPEC.items():
        S = {k: [] for k in TERMS}; T = []
        for d in dirs:
            for k in TERMS:
                t, c = spec2d.by_ky(f"{d}/Spectral_2D_{k}.dat")
                S[k].append(np.c_[c[:, 0], c[:, 1:].sum(axis=1)])            # zonal, non-zonal
            T.append(t)
        t = np.concatenate(T); o = np.argsort(t); t = t[o]
        S = {k: np.concatenate(v)[o] for k, v in S.items()}
        tp, ez, en, zn = traces.plunk(run)
        wid, step, c0, c1 = WINDOWS[run]
        tc = np.arange(c0, c1 + 1, step)
        r = {k: [] for k in ("W", "dW", "Sv", "Sz", "Sp", "Tdir", "Svz", "Spz", "Wz", "a")}
        for c in tc:
            m = (t >= c - wid / 2) & (t < c + wid / 2); mp = (tp >= c - wid / 2) & (tp < c + wid / 2)
            r["W"].append(S["fe"][m, 1].mean()); r["dW"].append(np.polyfit(t[m], S["fe"][m, 1], 1)[0])
            r["Wz"].append(S["fe"][m, 0].mean())
            r["Sv"].append(-S["hyp_v"][m, 1].mean()); r["Sz"].append(-S["hyp_z"][m, 1].mean())
            r["Sp"].append(-S["hyp_kperp"][m, 1].mean()); r["Tdir"].append(S["nonlinear"][m, 1].mean())
            r["Svz"].append(-S["hyp_v"][m, 0].mean()); r["Spz"].append(-S["hyp_kperp"][m, 0].mean())
            r["a"].append(np.sqrt(en[mp].mean() / ez[mp].mean()))
        out[f"tc_{run}"] = tc
        for k, v in r.items():
            out[f"{k}_{run}"] = np.array(v)
        print(run, "windows", len(tc))
    np.savez(CACHE, **out)


def rates(C, run):
    W = C[f"W_{run}"]; S = C[f"Sv_{run}"] + C[f"Sz_{run}"] + C[f"Sp_{run}"]
    return S / W, (C[f"dW_{run}"] + S) / W, -C[f"dW_{run}"] / W


def draw_budget(C):
    tc = C["tc_orig"]; pre = tc <= traces.COLLAPSE[1] - 200.0
    S, T, G = rates(C, "orig")
    fig, (a, b) = ns.two_panel()
    traces.collapse_band(a)
    for ax in (a, b):
        ax.axhline(0.0, c="0.6", lw=0.6, zorder=1)
    a.plot(tc[pre] / 1e3, 1e3 * S[pre], c=ns.PAL["Z"], ls="--", label="removal")
    a.plot(tc[pre] / 1e3, 1e3 * T[pre], c=ns.PAL["W"], ls="-.", label="supply")
    a.plot(tc[pre] / 1e3, 1e3 * G[pre], c="k", ls="-", label="net decay")
    a.set_xlim(1, 11.1); a.set_ylim(-1.4, 4.6)
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
    a.set_ylabel(r"rate $(10^{-3}c_{\rm ref}/L_{\rm ref})$")
    a.legend(loc="upper left", bbox_to_anchor=(0.16, 1.0), handlelength=2.0)
    # (b) against amplitude
    am = C["a_orig"]
    b.semilogx(am[pre], 1e3 * S[pre], c=ns.PAL["Z"], ls="--")
    b.semilogx(am[pre], 1e3 * T[pre], c=ns.PAL["W"], ls="-.")
    for run, mk in (("a03", "s"), ("a003", "^")):
        s, t_, g = rates(C, run)
        b.semilogx(C[f"a_{run}"], 1e3 * s, mk, c=ns.PAL["Z"], mfc="w", ms=4.5, mew=0.9)
        b.semilogx(C[f"a_{run}"], 1e3 * t_, mk, c=ns.PAL["W"], mfc="w", ms=4.5, mew=0.9)
    b.set_xlim(0.2, 1.4e-3); b.set_ylim(-7.4, 4.6)
    b.set_xlabel(r"$a=(E_{\rm nz}/E_{\rm zon})^{1/2}$")
    b.set_ylabel(r"rate $(10^{-3}c_{\rm ref}/L_{\rm ref})$")
    traces.label(b, 0.13, 2.6, "removal", c=ns.PAL["Z"]); traces.label(b, 0.13, -1.2, "supply", c=ns.PAL["W"])
    ns.style_axes(a)
    ns.tag(a, "(a)"); ns.tag(b, "(b)")
    ns.save(fig, "budget")


def draw_dissipation(C):
    tc = C["tc_orig"] / 1e3
    fig, (a, b) = ns.two_panel(sharey=True)
    n = C["Spz_orig"][np.argmin(np.abs(C["tc_orig"] - 6000.0))]              # unit: the zonal perpendicular sink
    for ax, kz, kn in ((a, "Spz_orig", "Sp_orig"), (b, "Svz_orig", "Sv_orig")):
        traces.collapse_band(ax)
        ax.semilogy(tc, C[kz] / n, c=ns.PAL["zon"], ls="-", lw=1.6)
        ax.semilogy(tc, C[kn] / n, c="k", ls="-")
        ax.set_xlim(1, 16); ax.set_ylim(2e-4, 20)
        ax.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
    a.set_ylabel(r"removal rate (normalised)")
    traces.label(a, 3.2, 2.4, "on the zonal flow", c="0.4"); traces.label(a, 2.0, 0.012, "on the vortices", c="k")
    traces.label(b, 6.3, 0.9, "on the vortices", c="k"); traces.label(b, 1.6, 0.012, "on the zonal flow", c="0.4")
    ns.tag(a, "(a) perpendicular", x=0.30, y=0.96); ns.tag(b, "(b) velocity space", x=0.30, y=0.96)
    plt = __import__("matplotlib.pyplot").pyplot
    plt.setp(b.get_yticklabels(), visible=False)
    ns.save(fig, "dissipation")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--rebuild", action="store_true"); A = ap.parse_args()
    if A.rebuild or not os.path.exists(CACHE):
        rebuild()
    C = np.load(CACHE)
    ns.apply_style()
    draw_dissipation(C)                       # draw_budget(C) is kept as a cross-check of figures/exchange; not in the notes
    for run in ("orig", "a03", "a003"):
        S, T, G = rates(C, run)
        st = 5 if run == "orig" else 1
        for i in range(0, len(S), st):
            sv = C[f"Sv_{run}"][i] / (C[f"Sv_{run}"][i] + C[f"Sz_{run}"][i] + C[f"Sp_{run}"][i])
            print(f"{run:5s} t {C['tc_'+run][i]:6.0f} a {C['a_'+run][i]:.2e}  S/W {1e3*S[i]:6.3f}  T/W {1e3*T[i]:7.3f}  net {1e3*G[i]:6.3f}  v-share {sv:.2f}")


if __name__ == "__main__":
    main()

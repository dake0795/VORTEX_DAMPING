#!/usr/bin/env python3
r"""Figure `euler`: the two-dimensional fluid model of the jet-vortex state against GENE (1 Oct 2026).

  python figure.py --rebuild   cache.npz from scripts/rayleigh/cache (GENE flute potentials gene_*.npz from extract.py,
                               fluid-model runs euler_*.npz from euler.py / batch1.sh)
  python figure.py             figures/euler.pdf from cache.npz alone

(a) E_nz(t)/E_nz(t_0), GENE against the fluid model started from GENE's potential at t_0;   point: the model keeps the
    vortices but they decay several times more slowly
(b) the eps = 0.03 restart, E_nz/E_nz(t_0), GENE, fluid model, and growth at the Rayleigh rate;   point: the regrowth
    is the shear-flow instability, at the full linear rate in the model and at about half of it in GENE
(c) fluid model with the perpendicular sink times 4, 1, 1/4;   point: the sink sets the decay and brings on the collapse
(d) decay rate against the sink strength; resolution check; GENE;   point: rate proportional to the sink, converged in
    resolution, below GENE's
All energies are those of the fluid model, (1/2) sum <k_perp^2> |psi_k|^2, evaluated identically for GENE's
field-line-averaged potential, as running means over ten breathing periods [two in (b)].
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns
import matplotlib.pyplot as plt

RAY = os.path.join(HERE, "..", "..", "scripts", "rayleigh")
CACHE = os.path.join(HERE, "cache.npz")
T0 = 5906.678
W_LONG, W_SHORT = 37.0, 7.4


def smooth(t, E, w):
    ti = np.arange(t[0], t[-1], 0.5)
    Ei = np.interp(ti, t, np.log(np.maximum(E, 1e-300)))
    n = max(1, int(w / 0.5)); k = np.ones(n) / n
    return ti[n // 2:n // 2 + len(Ei) - n + 1], np.exp(np.convolve(Ei, k, "valid"))


def rate(t, E, a, b):
    m = (t >= a) & (t <= b)
    return -np.polyfit(t[m], np.log(E[m]), 1)[0]


def rebuild():
    sys.path.insert(0, RAY)
    import measure as M
    out = {}

    def gene(runs):
        T, E = [], []
        for r in runs:
            t, psi = M.load(r); T.append(t); E.append(M.energies(psi)[1])
        T, E = np.concatenate(T), np.concatenate(E)
        u = np.r_[True, np.diff(T) > 0]
        return T[u] - T0, E[u]

    for key, runs in (("orig", ("leg4", "leg5")), ("a003", ("a003",))):
        t, E = gene(runs)
        for w, tagw in ((W_LONG, "L"), (W_SHORT, "S")):
            ts, Es = smooth(t, E, w)
            out[f"g_{key}_{tagw}_t"], out[f"g_{key}_{tagw}_E"] = ts[::4], Es[::4] / E[0]
    for name in ("base", "a003_long", "hx4", "hd4", "nx128p", "nx256p"):
        d = np.load(f"{RAY}/cache/euler_{name}.npz")
        for w, tagw in ((W_LONG, "L"), (W_SHORT, "S")):
            ts, Es = smooth(d["t"], d["Enz"], w)
            out[f"m_{name}_{tagw}_t"], out[f"m_{name}_{tagw}_E"] = ts[::4], Es[::4] / d["Enz"][0]
    g = np.load(f"{RAY}/cache/growth_history.npz")["rows"]
    out["gamma_R"] = float(np.interp(T0 + 30.0, g[:, 0], 0.5 * (g[:, 1] + g[:, 3])))
    np.savez(CACHE, **out)
    print("wrote", CACHE)


def main():
    if "--rebuild" in sys.argv:
        rebuild()
    C = np.load(CACHE)
    ns.apply_style()
    GEN, MOD = "k", ns.PAL["theory"]
    fig, ax = ns.grid_2xN(2, panel_h=1.62, hspace=0.62, left=0.66, wspace=0.80, top=0.12)
    (a, b), (c, d) = ax
    tl = r"$(t-t_0)\;(10^3L_{\rm ref}/c_{\rm ref})$"
    ye = r"$E_{\rm nz}(t)/E_{\rm nz}(t_0)$"

    # (a) decay: GENE against the model
    a.semilogy(C["g_orig_L_t"] / 1e3, C["g_orig_L_E"], c=GEN)
    a.semilogy(C["m_base_L_t"] / 1e3, C["m_base_L_E"], c=MOD)
    a.set_xlim(0, 5.6); a.set_ylim(1.5e-4, 3)
    a.set_xlabel(tl); a.set_ylabel(ye)
    a.text(2.9, 0.035, "GENE", color=GEN, ha="right"); a.text(3.3, 0.85, "fluid model", color=MOD, ha="left", va="bottom")
    ns.tag(a, "(a)", y=0.16)

    # (b) regrowth of the scaled-down restart
    tt = np.linspace(0, 55, 50)
    b.semilogy(tt, 1.25 * np.exp(2 * float(C["gamma_R"]) * tt), c="0.45", ls=":", lw=1.2, label="Rayleigh growth")
    b.semilogy(C["g_a003_S_t"], C["g_a003_S_E"], c=GEN)
    b.semilogy(C["m_a003_long_S_t"], C["m_a003_long_S_E"], c=MOD)
    b.set_xlim(0, 200); b.set_ylim(0.7, 80)
    b.set_xlabel(r"$(t-t_0)\;(L_{\rm ref}/c_{\rm ref})$"); b.set_ylabel(ye)
    b.text(150, 5.2, "GENE", color=GEN, ha="center", va="top"); b.text(150, 40, "fluid model", color=MOD, ha="center", va="bottom")
    b.legend(loc="lower right", handlelength=1.4, labelcolor="0.35")
    ns.tag(b, "(b)")

    # (c) the sink in the model
    for name, lab, (xt, yt) in (("hx4", r"$4\nu_\perp$", (1.62, 0.36)), ("base", r"$\nu_\perp$", (2.25, 0.60)),
                                ("hd4", r"$\nu_\perp/4$", (2.25, 0.89))):
        c.plot(C[f"m_{name}_L_t"] / 1e3, C[f"m_{name}_L_E"], c=MOD)
        c.text(xt, yt, lab, color=MOD, ha="left", va="top")
    c.set_xlim(0, 3.0); c.set_ylim(0, 1.25)
    c.set_xlabel(tl); c.set_ylabel(ye)
    ns.tag(c, "(c)", y=0.16)

    # (d) rate against sink strength
    W = (100.0, 800.0)
    r = {n: rate(C[f"m_{n}_L_t"], C[f"m_{n}_L_E"], *W) for n in ("hd4", "base", "hx4", "nx128p", "nx256p")}
    rg = rate(C["g_orig_L_t"], C["g_orig_L_E"], *W)
    xs = np.array([0.18, 5.5])
    d.loglog(xs, r["base"] * xs, c="0.6", ls="--", lw=0.9)
    d.loglog([0.25, 1, 4], [r["hd4"], r["base"], r["hx4"]], "o", c=MOD, ms=5.5)
    d.loglog([1, 1], [r["nx128p"], r["nx256p"]], "s", mfc="none", mec=MOD, ms=7.5, mew=1.0)
    d.loglog([1], [rg], "*", c=GEN, ms=10)
    d.set_xlim(0.15, 6.5); d.set_ylim(1.5e-5, 2.5e-3)
    d.set_xticks([0.25, 1, 4]); d.set_xticklabels([r"$1/4$", r"$1$", r"$4$"]); d.minorticks_off()
    d.set_xlabel(r"sink strength $/\,\nu_\perp$"); d.set_ylabel(r"$\Gamma_E\;(c_{\rm ref}/L_{\rm ref})$")
    d.text(1.25, rg, "GENE", color=GEN, ha="left", va="center")
    d.text(1.35, r["base"] * 0.62, "fluid model", color=MOD, ha="left", va="center")
    ns.tag(d, "(d)")
    for axx in (a, b, c):
        axx.minorticks_on()
    ns.save(fig, "euler")


if __name__ == "__main__":
    main()

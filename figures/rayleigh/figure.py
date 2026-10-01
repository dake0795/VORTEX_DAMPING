#!/usr/bin/env python3
r"""Figure `rayleigh`: the vortices are the shear-flow (Rayleigh) instability of the jets (1 Oct 2026).

  python figure.py --rebuild   cache.npz from scripts/rayleigh/cache (growth_history.npz, eigmode_leg4_0.npz,
                               measured_leg4_0.npz; made by scripts/rayleigh/{extract,measure,linear,growth_history}.py)
  python figure.py             figures/rayleigh.pdf from cache.npz alone

(a) zonal vorticity U'(x) of the measured jets against the single harmonic of the same amplitude; the critical layers
    U = +-c of the two vortices marked;           point: each shear layer is split into two vorticity strips, one per layer
(b) |psi(x)| of the unstable Rayleigh mode against the measured vortex;   point: same radial structure
(c) frequency of the unstable Rayleigh mode on the running-mean jet against the measured one, through the run
(d) its growth rate through the run;              point: positive throughout, gone at the collapse
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
LAMD = 50.0 * np.sqrt(2.0)
T_COLLAPSE = (10600.0, 11000.0)


def rebuild():
    sys.path.insert(0, RAY)
    import measure as M
    from common import KYMIN
    g = np.load(f"{RAY}/cache/growth_history.npz")["rows"]
    e = np.load(f"{RAY}/cache/eigmode_leg4_0.npz")
    m = np.load(f"{RAY}/cache/measured_leg4_0.npz")
    nf = len(e["x"])
    p1 = np.zeros(64, complex); p1[1] = m["psi0"][1]; p1[-1] = m["psi0"][-1]
    U = e["U"]; c = float(m["w"]) / KYMIN
    lay = {s: [float(e["x"][j]) for j in range(nf) if (U[j] - s * c) * (U[(j + 1) % nf] - s * c) < 0] for s in (1, -1)}
    np.savez(CACHE, x=e["x"], U=U, S=M.fine(m["psi0"], nf, 2).real, S1=M.fine(p1, nf, 2).real, eig=np.abs(e["eig"]),
             meas=np.abs(e["meas"]), xplus=lay[1], xminus=lay[-1], hist=g, t_profile=float(m["t"].mean()))
    print("wrote", CACHE)


def main():
    if "--rebuild" in sys.argv:
        rebuild()
    C = np.load(CACHE)
    ns.apply_style()
    fig, ax = ns.grid_2xN(2, panel_h=1.62, hspace=0.62, left=0.60, wspace=0.78, top=0.12)
    (a, b), (c, d) = ax
    x = C["x"] / LAMD
    lay = [(xx / LAMD, ns.PAL["plus"]) for xx in C["xplus"]] + [(xx / LAMD, ns.PAL["minus"]) for xx in C["xminus"]]

    # (a) zonal vorticity
    a.plot(x, C["S1"], c=ns.PAL["zon"], ls="--", lw=1.0)
    a.plot(x, C["S"], c="k")
    for xx, col in lay:
        a.axvline(xx, c=col, lw=0.8, ls=":")
    a.set_xlim(x[0], x[-1] + x[1]); a.set_ylim(-2.75, 2.75)
    a.set_xlabel(r"$x/\lambda_D$"); a.set_ylabel(r"$U'\;(c_{\rm ref}/L_{\rm ref})$")
    wb = dict(facecolor="white", edgecolor="none", pad=1.5)
    a.text(0.965, 0.04, "single harmonic", transform=a.transAxes, ha="right", va="bottom", color="0.4", bbox=wb, zorder=7)
    a.text(0.965, 0.16, "jets", transform=a.transAxes, ha="right", va="bottom", color="k", bbox=wb, zorder=7)
    ns.tag(a, "(a)")

    # (b) mode structure
    b.plot(x, C["meas"] / C["meas"].max(), c=ns.PAL["plus"], lw=1.6, label="vortex")
    b.plot(x, C["eig"] / C["eig"].max(), c="k", ls="--", lw=1.0, label="Rayleigh mode")
    for xx, col in lay[:2]:
        b.axvline(xx, c=col, lw=0.8, ls=":")
    b.set_xlim(x[0], x[-1] + x[1]); b.set_ylim(0, 1.55)
    b.set_xlabel(r"$x/\lambda_D$"); b.set_ylabel(r"$|\psi_+|/\max|\psi_+|$")
    b.legend(loc="upper center", bbox_to_anchor=(0.51, 1.0), handlelength=1.3)
    b.set_yticks([0, 0.5, 1.0])
    ns.tag(b, "(b)")

    # (c), (d) history
    h = C["hist"]
    t = h[:, 0] / 1e3
    ok = (h[:, 1] > 2e-3) & (h[:, 3] > 2e-3) & (h[:, 0] > 800)   # both unstable modes exist
    meas = (h[:, 10] < 0.3) & (h[:, 0] > 700)                # the two-vortex fit holds
    for axx in (c, d):
        axx.axvspan(T_COLLAPSE[0] / 1e3, T_COLLAPSE[1] / 1e3, color="0.88", lw=0)
        axx.set_xlim(0, 11.6); axx.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$")
    c.plot(t[meas], h[meas, 9], c=ns.PAL["plus"], lw=1.6, label="vortex")
    c.plot(t[ok], 0.5 * (h[ok, 2] + h[ok, 4]), c="k", ls="--", lw=1.0, label="Rayleigh mode")
    c.set_ylim(0.7, 1.3); c.set_ylabel(r"$|\omega|\;(c_{\rm ref}/L_{\rm ref})$")
    c.legend(loc="upper center", bbox_to_anchor=(0.56, 1.0))
    ns.tag(c, "(c)", x=0.04, y=0.16)
    g = np.where((h[:, 1] > 2e-3) & (h[:, 3] > 2e-3), 0.5 * (h[:, 1] + h[:, 3]), np.maximum(np.maximum(h[:, 1], h[:, 3]), 0.0))
    g[h[:, 0] < 700] = np.nan
    d.plot(t, np.where(g > 1e-3, g, 0.0) * 1e2, c="k")
    d.set_ylim(0, 9.5); d.set_ylabel(r"$\gamma_{\rm R}\;(10^{-2}c_{\rm ref}/L_{\rm ref})$")
    d.text(T_COLLAPSE[0] / 1e3 - 0.25, 7.3, "collapse", rotation=90, ha="right", va="center", color="0.4")
    ns.tag(d, "(d)", x=0.12, y=0.16)
    for axx in (a, b, c, d):
        ns.style_axes(axx)
    ns.save(fig, "rayleigh")


if __name__ == "__main__":
    main()

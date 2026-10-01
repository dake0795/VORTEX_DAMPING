#!/usr/bin/env python3
r"""Figure `sinktests`: what sets the decay of the vortices - restarts with one sink changed (1 Oct 2026).

  python figure.py --rebuild   read plunk_e_time.dat of leg 5 of the reference run (baseline) and of the restarts
                               rb_c0p3_hxy_sinktest/{hv_d4,hxy_x4,hxy_d4,nx128}; write cache.npz
  python figure.py             draw ../sinktests.pdf from cache.npz

All runs start from the same state (t0 = 9423.17, the accelerating phase of the decay of the reference run):
  hv_d4   velocity-space hyperdiffusion / 4
  hxy_x4  perpendicular (x and y) hyperdiffusion x 4
  hxy_d4  perpendicular hyperdiffusion / 4
  nx128   radial resolution doubled at the same grid-normalised coefficient: radial hyperdiffusion / 16 at a given k_x
(a) E_nz(t)/E_nz(t0) (running means over ten breathing periods);
(b) the mean decay rate Gamma_E over the common window against the strength of the radial hyperdiffusion relative to
    the baseline, with the loss rate to entropy L/E_nz (figures/channels) as a horizontal line.
The runs are still advancing: rebuild to extend the curves.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns
import traces

T0, WIN = 9423.17, 37.0
CACHE = os.path.join(HERE, "cache.npz")
RUNS = {"base": f"{ns.RB}/rb_c0p3_hxy/leg_0005/out",
        **{r: f"{ns.RB}/rb_c0p3_hxy_sinktest/{r}/out" for r in ("hv_d4", "hxy_x4", "hxy_d4", "nx128")}}
STRENGTH = {"hxy_x4": 4.0, "base": 1.0, "hxy_d4": 0.25, "nx128": 1.0 / 16.0}
LOSS = 1.8e-3                         # L/E_nz, the loss to entropy through parallel streaming (figures/channels)


def load(d):
    rows = []
    for l in open(os.path.join(d, "plunk_e_time.dat")):
        if l.startswith("#"):
            continue
        v = l.split()
        if len(v) >= 17:
            try:
                rows.append([float(v[0]), float(v[11])])
            except ValueError:
                pass
    return np.array(rows)


def rebuild():
    out = {}
    for r, d in RUNS.items():
        a = load(d); out[r] = a[a[:, 0] <= 11500.0]
        print(r, f"t {a[0, 0]:.1f} .. {a[-1, 0]:.1f}")
    np.savez(CACHE, **out)


def draw():
    C = np.load(CACHE)
    tend = min(C[r][-1, 0] for r in C.files if r != "base")
    ns.apply_style()
    fig, (a, b) = ns.two_panel()
    sty = {"base": ("k", "-", 1.5), "hv_d4": (ns.PAL["Z"], "--", 1.1), "hxy_x4": (ns.PAL["W"], "-", 1.1),
           "hxy_d4": (ns.PAL["W"], "--", 1.1), "nx128": (ns.PAL["aux"], "-", 1.1)}
    rate = {}
    for r in ("base", "hxy_x4", "hxy_d4", "nx128", "hv_d4"):
        t, e = C[r][:, 0], C[r][:, 1]
        es = traces.runmean(t, e, WIN, log=True); ok = np.isfinite(es); ts, es = t[ok], es[ok]
        m = ts <= (tend if r == "base" else ts[-1])
        eb = traces.runmean(C["base"][:, 0], C["base"][:, 1], WIN, log=True); kb = np.isfinite(eb)
        e0 = eb[kb][0]
        c, ls, lw = sty[r]
        a.semilogy((ts[m] - T0), es[m] / e0, ls, c=c, lw=lw)
        mm = (ts >= T0 + 60.0) & (ts <= tend - WIN / 2)
        rate[r] = -np.polyfit(ts[mm], np.log(es[mm]), 1)[0]
        print(r, f"mean Gamma_E over t0+60 .. {tend:.0f}: {rate[r]:.3e}")
    a.set_xlim(0, tend - T0 + 20); a.set_ylim(0.03, 1.6)
    a.set_xlabel(r"$t-t_0\;(L_{\rm ref}/c_{\rm ref})$"); a.set_ylabel(r"$E_{\rm nz}(t)/E_{\rm nz}(t_0)$")
    xe = tend - T0
    traces.label(a, 0.40 * xe, 0.11, r"$\nu_\perp\times4$", c=ns.PAL["W"])
    traces.label(a, 0.04 * xe, 0.40, r"reference and $\nu_v/4$", c="k")
    traces.label(a, 0.60 * xe, 0.56, r"$\nu_\perp/4$", c=ns.PAL["W"])
    traces.label(a, 0.62 * xe, 1.25, r"$\Delta x/2$", c=ns.PAL["aux"])
    xs = np.array([STRENGTH[r] for r in ("hxy_x4", "base", "hxy_d4", "nx128")])
    ys = np.array([rate[r] for r in ("hxy_x4", "base", "hxy_d4", "nx128")])
    b.loglog([1.0], [rate["hv_d4"]], "s", c=ns.PAL["Z"], mfc="w", ms=10, mew=1.1)
    b.loglog(xs, ys, "o", c="k", ms=5)
    b.axhline(LOSS, c=ns.PAL["theory"], ls=":", lw=1.4)
    b.set_xlim(0.03, 8); b.set_ylim(8e-5, 2e-2)
    b.set_xlabel(r"radial hyperdiffusion / reference"); b.set_ylabel(r"$\Gamma_E\;(c_{\rm ref}/L_{\rm ref})$")
    from matplotlib.ticker import FixedLocator, NullFormatter
    b.xaxis.set_major_locator(FixedLocator([1 / 16, 0.25, 1, 4])); b.set_xticklabels([r"$1/16$", r"$1/4$", r"$1$", r"$4$"])
    b.xaxis.set_minor_formatter(NullFormatter())
    traces.label(b, 0.04, 2.6e-3, r"loss to entropy", c=ns.PAL["theory"])
    ns.tag(a, "(a)", y=0.14, va="bottom"); ns.tag(b, "(b)")
    ns.save(fig, "sinktests")


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        rebuild()
    draw()

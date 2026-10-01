#!/usr/bin/env python3
r"""Figure `sinktests`: the sink scan and the grid test (restarts of 1 Oct 2026 from t = 9423).

  python figure.py --rebuild   read plunk_e_time.dat of leg 5 (baseline) and of the five restarts, write cache.npz
  python figure.py             draw ../sinktests.pdf from cache.npz

(a) E_nz(t)/E_nz(t0) of the baseline and of the four restarts with one sink changed by a factor of four.
(b) the same for the baseline and the restart with twice the radial resolution.
All curves are running means over ten breathing periods (37 time units).
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns

T0, TEND, WIN = 9423.17, 12500.0, 37.0
CACHE = os.path.join(HERE, "cache.npz")
RUNS = {"base": f"{ns.RB}/rb_c0p3_hxy/leg_0005/out", **{r: f"{ns.RB}/rb_c0p3_hxy_sinktest/{r}/out"
                                                        for r in ("hv_x4", "hv_d4", "hxy_x4", "hxy_d4", "nx128")}}


def load(d):
    p = os.path.join(d, "plunk_e_time.dat")
    if not os.path.exists(p):
        return None
    rows = []
    for l in open(p):
        if l.startswith("#"):
            continue
        v = l.split()
        if len(v) >= 17:
            try:
                rows.append([float(v[0]), float(v[10]), float(v[11]), float(v[15])])
            except ValueError:
                pass
    a = np.array(rows)
    return a[(a[:, 0] >= T0 - 1) & (a[:, 0] <= TEND)] if len(a) else None


def smooth(t, y):
    """Running geometric mean over WIN time units on the (nearly uniform) output grid."""
    n = max(1, int(round(WIN / np.median(np.diff(t)))))
    if len(y) < 2 * n:
        return t[:0], y[:0]
    k = np.ones(n) / n
    return np.convolve(t, k, "valid"), np.exp(np.convolve(np.log(y), k, "valid"))


def rebuild():
    out = {}
    for r, d in RUNS.items():
        a = load(d)
        if a is None or len(a) < 10:
            print(r, "no data yet"); continue
        out[r] = a; print(r, f"t {a[0, 0]:.1f} .. {a[-1, 0]:.1f}, {len(a)} records")
    np.savez(CACHE, **out)


def draw():
    C = np.load(CACHE)
    ns.apply_style()
    fig, (a, b) = ns.two_panel(sharey=True)
    e0 = None
    sty = {"base": ("k", "-", 1.4, "reference"), "hv_x4": (ns.PAL["Z"], "-", 1.1, r"$\nu_v\times4$"),
           "hv_d4": (ns.PAL["Z"], "--", 1.1, r"$\nu_v/4$"), "hxy_x4": (ns.PAL["W"], "-", 1.1, r"$\nu_\perp\times4$"),
           "hxy_d4": (ns.PAL["W"], "--", 1.1, r"$\nu_\perp/4$"), "nx128": (ns.PAL["aux"], "-", 1.1, r"$\Delta x/2$")}
    for r in ("base", "hv_x4", "hv_d4", "hxy_x4", "hxy_d4", "nx128"):
        if r not in C.files:
            continue
        t, en = C[r][:, 0], C[r][:, 2]
        ts, es = smooth(t, en)
        if e0 is None:
            e0 = es[0]
        c, ls, lw, lab = sty[r]
        for ax in ((a, b) if r == "base" else (b,) if r == "nx128" else (a,)):
            ax.semilogy(ts - T0, es / e0, ls, c=c, lw=lw, label=lab)
    for ax in (a, b):
        ax.set_xlabel(r"$t-t_0\;(L_{\rm ref}/c_{\rm ref})$"); ax.legend(loc="lower left")
    a.set_ylabel(r"$E_{\rm nz}(t)/E_{\rm nz}(t_0)$")
    ns.tag(a, "(a)", x=0.86); ns.tag(b, "(b)", x=0.86)
    ns.save(fig, "sinktests")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rebuild", action="store_true"); A = ap.parse_args()
    if A.rebuild or not os.path.exists(CACHE):
        rebuild()
    draw()

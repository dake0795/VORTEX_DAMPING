#!/usr/bin/env python3
r"""Figure `exchange` of the vortex-damping notes: where the electrostatic energy of the vortices goes.

  env -u LD_LIBRARY_PATH python3 figure.py            plot from cache.npz -> ../exchange.pdf (+ .png)
  python3 figure.py --rebuild                         read the GENE rate diagnostics and rewrite cache.npz first

Data (read only): GENE's mode-resolved budget of the Debye electrostatic energy E_D, Spectral_2D_Eplunk_<species>.dat
and Spectral_2D_EDrate<channel>_<species>.dat (every 100 steps; text, one value per line: header = k_x axis (NX0 - 1
values) and k_y axis, then per frame: time, total, data[k_x, k_y]), of the reduced-box run rb_c0p3_hxy (legs 1, 3-5)
and its two amplitude-reduced restarts.  Conventions of ep_turbulence_paper/scripts_letter/lfig09_breathing: species
summed, k_y > 0 counted twice, E_D = Eplunk / 2, rates = EDrate<ch> / 2, parallel = parallel_raw - hyp_v - hyp_z.

For the non-zonal part (k_y >= 1) each channel is fitted, in consecutive windows of WIN time units, to
    const + slope tau + two harmonics of the breathing frequency w_b
(w_b from ../trapping/cache.npz, so run ../trapping/figure.py --rebuild first) and the constant, divided by the
window-mean non-zonal E_D, is kept: the cycle-averaged rate per unit vortex energy.  Channels:
  NL      nonlinear exchange with the zonal flow (exactly minus the zonal NL rate; the same number follows from the
          potential alone, z by z, as the Reynolds-stress work: ../measurements/reynolds_transfer_z.py);
  linear  parallel + curvature, i.e. the transit-averaged magnetic-drift term, which moves electrostatic energy into
          entropy;
  hyp     the perpendicular hyperdiffusion acting directly on E_D (negligible).
The companion figure `budget` is the same decomposition for the free energy W.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
CACHE = os.path.join(HERE, "cache.npz")
WIN = 100.0
W = 24                                           # bytes per value line
STEMS = ("Eplunk", "EDrateNL", "EDratecurvature", "EDrateparallel_raw", "EDratehyp_v", "EDratehyp_z", "EDratehyp_kperp")
T_COLLAPSE = (10600.0, 11100.0)


def read_ky(path, nkx, nky):
    """All whole frames of one Spectral_2D file, summed over k_x: t[n], data[n, nky]."""
    hdr, frl = nkx + nky, 2 + nkx * nky
    n = (os.path.getsize(path) // W - hdr) // frl
    T, D = [], []
    with open(path, "rb") as f:
        f.seek(hdr * W)
        left = n
        while left > 0:
            m = min(left, 400)
            a = np.array(f.read(m * frl * W).split(), float).reshape(m, frl)
            T.append(a[:, 0]); D.append(a[:, 2:].reshape(m, nkx, nky).sum(axis=1))
            left -= m
    return np.concatenate(T), np.concatenate(D)


def series(run, vf):
    nkx, nky = vf.NX - 1, vf.NKY
    wgt = np.full(nky, 2.0); wgt[0] = 1.0
    out, last = None, -np.inf
    for o in vf.RUNS[run][:-1] if run == "orig" else vf.RUNS[run]:      # leg 6 (t > 16,300) is not needed here
        leg = {}
        for stem in STEMS:
            xs = [read_ky(f"{o}/Spectral_2D_{stem}_{sp}.dat", nkx, nky) for sp in ("electrons", "positrons")]
            n = min(len(x[0]) for x in xs)
            leg[stem] = 0.5 * (xs[0][1][:n] + xs[1][1][:n]) * wgt
            t = xs[0][0][:n]
        n = min(len(v) for v in leg.values())
        leg = {k: v[:n] for k, v in leg.items()}
        t = t[:n]
        m = t > last + 1e-9
        last = t[-1]
        leg = {k: v[m] for k, v in leg.items()}
        leg["t"] = t[m]
        out = leg if out is None else {k: np.concatenate([out[k], leg[k]]) for k in leg}
        print(run, o.split("/")[-2], m.sum(), "frames", flush=True)
    return out


def fitc(tau, y, wb):
    M = np.stack([np.ones_like(tau), tau, np.cos(wb * tau), np.sin(wb * tau), np.cos(2 * wb * tau), np.sin(2 * wb * tau)], axis=1)
    return np.linalg.lstsq(M, y, rcond=None)[0]


def rebuild():
    import vortex_field as vf
    tr = np.load(os.path.join(HERE, "..", "trapping", "cache.npz"))
    c = {}
    for run, t0, t1 in (("orig", 1000.0, 11100.0), ("a03", 5980.0, 6705.0), ("a003", 5980.0, 6705.0)):
        s = series(run, vf)
        t = s["t"]
        nz = lambda k: s[k][:, 1:].sum(axis=1)
        lin = nz("EDrateparallel_raw") - nz("EDratehyp_v") - nz("EDratehyp_z") + nz("EDratecurvature")
        rows = []
        for a in np.arange(t0, t1 - WIN + 1e-6, WIN):
            m = (t >= a) & (t < a + WIN)
            tc = a + WIN / 2
            wb = np.interp(tc, tr[run][:, 0], tr[run][:, 1] - tr[run][:, 2])
            tau = t[m] - tc
            En, Ez = nz("Eplunk")[m].mean(), s["Eplunk"][m, 0].mean()
            rows.append([tc, En / Ez, fitc(tau, nz("EDrateNL")[m], wb)[0] / En, fitc(tau, lin[m], wb)[0] / En,
                         fitc(tau, nz("EDratehyp_kperp")[m], wb)[0] / En, fitc(tau, nz("Eplunk")[m], wb)[1] / En,
                         fitc(tau, s["EDrateNL"][m, 0], wb)[0] / En])
        c[run] = np.array(rows)
    np.savez_compressed(CACHE, cols="t_c  E_nz/E_zon  NL  linear(parallel+curvature)  hyp_kperp  dlnE_nz/dt  NL_zonal  (rates per unit non-zonal E_D)", **c)


def block(r, edges):
    """means and standard errors of the rows over the groups of windows [i0, i1) in `edges`."""
    mu = np.array([r[a:b].mean(axis=0) for a, b in edges])
    se = np.array([r[a:b].std(axis=0, ddof=1) / np.sqrt(b - a) for a, b in edges])
    return mu, se


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--table", action="store_true", help="print the cached rates (units 1e-3)")
    a = ap.parse_args()
    if a.rebuild:
        rebuild()
    import notes_style as ns
    ns.apply_style()
    PAL = ns.PAL
    C = np.load(CACHE)
    o = C["orig"]
    if a.table:
        for run in ("orig", "a03", "a003"):
            print(run, str(C["cols"]))
            print(np.array2string(np.c_[C[run][:, 0], C[run][:, 1], C[run][:, 2:] * 1e3], precision=3, suppress_small=True, max_line_width=160))
    fig, (ax, bx) = ns.two_panel(left=0.62, wspace=0.80)
    u = 1e3
    c_nl, c_lin = PAL["aux"], PAL["theory"]
    e4 = [(i, i + 4) for i in range(0, len(o) - 3, 4)]

    # (a) the one point: the loss to entropy is steady; what changes is the nonlinear exchange with the jet
    ax.axvspan(*T_COLLAPSE, color="0.9", lw=0)
    ax.axhline(0, c="0.6", lw=0.6)
    mu, se = block(o, e4)
    for col, c_ in ((3, c_lin), (2, c_nl), (5, "k")):
        ax.fill_between(mu[:, 0], u * (mu[:, col] - se[:, col]), u * (mu[:, col] + se[:, col]), color=c_, alpha=0.22, lw=0)
        ax.plot(mu[:, 0], u * mu[:, col], c=c_)
    ax.text(5000, 2.05, "from the jet", color=c_nl, ha="center", va="bottom")
    ax.text(5000, -2.35, "to entropy", color=c_lin, ha="center", va="top")
    ax.text(6500, -0.72, "net", color="k", ha="center", va="top")
    ax.set_xlim(0, 11500)
    ax.set_xticks([0, 5000, 10000])
    ax.set_ylim(-9.5, 4.8)
    ax.set_xlabel(r"$t\;(L_{\rm ref}/c_{\rm ref})$")
    ax.set_ylabel(r"$\dot E_{\rm nz}/E_{\rm nz}\;(10^{-3}\,c_{\rm ref}/L_{\rm ref})$")
    ns.style_axes(ax)
    ns.tag(ax, "(a)", y=0.06, va="bottom")

    # (b) the one point: the exchange with the jet is a function of the vortex amplitude, the same in the restarts,
    # and it changes sign at small amplitude
    lin = -u * np.median(o[:, 3])
    bx.axhline(0, c="0.6", lw=0.6)
    bx.axhline(lin, c=c_lin, ls=":", lw=1.0)
    bx.text(1.8e-1, lin + 0.15, "loss to entropy", color=c_lin, ha="left", va="bottom")
    for run, mk, ms, edges, lab in (("orig", "o", 3.2, e4, "original"), ("a03", "s", 3.8, [(0, 7)], r"$\epsilon=0.3$"),
                                    ("a003", "^", 4.0, [(0, 2), (2, 4), (4, 7)], r"$\epsilon=0.03$")):
        mu, se = block(C[run], edges)
        bx.errorbar(np.sqrt(mu[:, 1]), u * mu[:, 2], yerr=u * se[:, 2], fmt=mk, ms=ms, c=PAL[run], mew=0, elinewidth=0.7, label=lab)
    bx.set_xscale("log")
    bx.set_xlim(2.2e-1, 9e-4)
    bx.set_ylim(-9.5, 4.8)
    bx.set_xlabel(r"$(E_{\rm nz}/E_{\rm zon})^{1/2}$")
    bx.set_ylabel(r"from the jet $(10^{-3}\,c_{\rm ref}/L_{\rm ref})$")
    bx.legend(loc="lower left", handletextpad=0.1, borderaxespad=0.2)
    ns.tag(bx, "(b)", x=0.955, ha="right")
    ns.save(fig, "exchange")


if __name__ == "__main__":
    main()

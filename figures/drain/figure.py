#!/usr/bin/env python3
r"""Figures `drain` and `supply`: the slow drain of the vortices and who pays for it (2 Oct 2026).

  python figure.py --rebuild   read plunk_e_time.dat of the reference run rb_c0p3_hxy (legs 1, 3-6), of the same run
                               from noise without radial hyperdiffusion rb_c0p3_hxoff (all legs), and of the restarts
                               rb_c0p3_hxy_sinktest/{nx128, hxy_d4} with their continuation legs; write cache.npz
  python figure.py             draw ../drain.pdf and ../supply.pdf from cache.npz

rb_c0p3_hxoff is the reference deck with hyp_x = 0 (hyp_y kept), started from the same noise.

drain   (a) E_nz(t)/E_1 (running means over ten breathing periods), E_1 = E_nz at T1, for the two runs, with the law
            E_nz = E_1 [1 + (3/4) K E_1^{3/4} (t - T1)]^{-4/3};
        (b) (E_nz/E_1)^{-3/4} against t: a straight line under dE_nz/dt = -K E_nz^{7/4}.
        K is fitted to the run without radial hyperdiffusion over T1 <= t <= its end (one number, printed).
supply  (a) the rate -dE_zon/dt at which the jets lose energy, measured in windows of WBIN time units, against
            (nu_L - Gamma_E) E_nz + nu_h E_zon: the supply T = L - Gamma_E E_nz of eq. (budget), with nu_L the loss rate to
            entropy of figures/channels (parallel + drift channel), plus the radial hyperdiffusion of the box-scale
            jet, nu_h = 3.4e-6 times the strength of the radial hyperdiffusion relative to the reference run (0 for
            rb_c0p3_hxoff, 1/16 for nx128, 1/4 for hxy_d4); no adjustable parameter;
        (b) [removed 4 Oct 2026: the E_nz^{1/4} drain law failed after t ~4500] E_nz^{1/4} against the electrostatic energy E_zon + E_nz of the
            run without radial hyperdiffusion: a straight line if the trapping frequency is proportional to the energy
            in excess of that of the stable state.
"""
import glob, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns
import traces

CACHE = os.path.join(HERE, "cache.npz")
T1, WIN, WBIN = 1000.0, 37.0, 400.0
NU_H = 3.4e-6                                  # energy decay rate of the box-scale jet under the radial hyperdiffusion (scripts/rayleigh/README.md)
STRENGTH = {"ref": 1.0, "hxoff": 0.0, "nx128": 1.0 / 16.0, "hxy_d4": 0.25}
T_RANGE = {"ref": (T1, 10500.0), "hxoff": (T1, 1e9), "nx128": (9500.0, 1e9), "hxy_d4": (9500.0, 1e9)}


def load(dirs):
    rows = []
    for d in dirs:
        p = os.path.join(d, "plunk_e_time.dat")
        if not os.path.exists(p):
            continue
        for l in open(p):
            if l.startswith("#"):
                continue
            v = l.split()
            if len(v) >= 17:
                try:
                    rows.append([float(v[0]), float(v[10]), float(v[11])])
                except ValueError:
                    pass
    a = np.array(rows); a = a[np.argsort(a[:, 0], kind="stable")]
    return a[np.concatenate([[True], np.diff(a[:, 0]) > 0])]


def rebuild():
    S = f"{ns.RB}/rb_c0p3_hxy_sinktest"
    runs = {"ref": traces.RUNS["orig"], "hxoff": sorted(glob.glob(f"{ns.RB}/rb_c0p3_hxoff/leg_*/out")),
            "nx128": [f"{S}/nx128/out", f"{S}/nx128_leg2/out"], "hxy_d4": [f"{S}/hxy_d4/out", f"{S}/hxy_d4_leg2/out"]}
    out = {}
    for r, dirs in runs.items():
        a = load(dirs); a = a[a[:, 0] <= 13000.0][::2]
        out[r] = a.astype(np.float64); print(r, f"t {a[0, 0]:.1f} .. {a[-1, 0]:.1f}, {len(a)} records")
    ch = np.load(os.path.join(HERE, "..", "channels", "cache.npz"))
    m = ch["tc_orig"] < 9000.0
    out["nuL"] = float((ch["par_orig"] + ch["curv_orig"])[m].mean())
    np.savez(CACHE, **out)


def smooth(a):
    t = a[:, 0]; e = traces.runmean(t, a[:, 2], WIN, log=True); ok = np.isfinite(e)
    return t[ok], e[ok]


def bins(a, t0, t1):
    """Windows of WBIN: centre, log-mean E_nz, mean E_zon, Gamma_E, -dE_zon/dt."""
    out = []
    t1 = min(t1, a[-1, 0])
    for c in np.arange(t0 + WBIN / 2, t1 - WBIN / 2 + 1e-6, WBIN):
        m = (a[:, 0] >= c - WBIN / 2) & (a[:, 0] < c + WBIN / 2)
        if m.sum() < 20:
            continue
        out.append([c, np.exp(np.log(a[m, 2]).mean()), a[m, 1].mean(), -np.polyfit(a[m, 0], np.log(a[m, 2]), 1)[0],
                    -np.polyfit(a[m, 0], a[m, 1], 1)[0]])
    return np.array(out)


def draw():
    C = np.load(CACHE)
    nuL = float(C["nuL"])
    ns.apply_style()
    col = {"ref": "k", "hxoff": ns.PAL["W"], "nx128": ns.PAL["aux"], "hxy_d4": ns.PAL["W"]}
    # ------------------------------------------------------------------ the law, fitted to the run without radial hyperdiffusion
    th, eh = smooth(C["hxoff"]); tr, er = smooth(C["ref"])
    E1 = float(np.interp(T1, th, eh)); E1r = float(np.interp(T1, tr, er))
    m = th >= T1
    slope = np.polyfit(th[m] - T1, (eh[m] / E1) ** -0.75 - 1.0, 1)
    s = float(np.sum((th[m] - T1) * ((eh[m] / E1) ** -0.75 - 1.0)) / np.sum((th[m] - T1) ** 2))     # line through (T1, 1)
    K = s / 0.75 / E1 ** 0.75
    res = (eh[m] / E1) ** -0.75 - (1.0 + s * (th[m] - T1))
    print(f"nu_L = {nuL:.3e};  run without radial hyperdiffusion: t to {th[-1]:.0f};  (E/E_1)^(-3/4) = 1 + {s:.3e} (t - T1), "
          f"rms residual {res.std():.3f} of a range {(eh[m][-1]/E1)**-0.75 - 1:.2f};  K = {K:.3e} (free-intercept slope {slope[0]:.3e}, intercept {1+slope[1]:.3f})")
    mr = (tr >= T1) & (tr <= th[-1])
    print(f"reference over the same interval: E_1 ratio ref/hxoff {E1r/E1:.3f}; mean |ln(E_ref/E_hxoff)| = "
          f"{np.abs(np.log(np.interp(tr[mr], th, eh) / er[mr])).mean():.3f}; (E/E_1)^(-3/4) at the end: ref {(er[mr][-1]/E1r)**-0.75:.2f}, hxoff {(eh[-1]/E1)**-0.75:.2f}")
    # ------------------------------------------------------------------ figure drain
    fig, (a, b) = ns.two_panel()
    TMAX = 12.0
    tl = np.linspace(T1, TMAX * 1e3, 400); law = (1.0 + s * (tl - T1)) ** (-4.0 / 3.0)
    a.semilogy(tr / 1e3, er / E1r, c="k", lw=1.5)
    a.semilogy(th / 1e3, eh / E1, c=col["hxoff"], lw=1.3)
    a.semilogy(tl / 1e3, law, "--", c=ns.PAL["theory"], lw=1.3)
    a.set_xlim(0, TMAX); a.set_ylim(1e-3, 30)
    a.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$"); a.set_ylabel(r"$E_{\rm nz}/E_1$")
    b.plot(tr / 1e3, (er / E1r) ** -0.75, c="k", lw=1.5)
    b.plot(th / 1e3, (eh / E1) ** -0.75, c=col["hxoff"], lw=1.3)
    b.plot(tl / 1e3, 1.0 + s * (tl - T1), "--", c=ns.PAL["theory"], lw=1.3)
    b.set_xlim(0, TMAX); b.set_ylim(0, 16)
    b.set_xlabel(r"$t\;(10^3L_{\rm ref}/c_{\rm ref})$"); b.set_ylabel(r"$(E_{\rm nz}/E_1)^{-3/4}$")
    from matplotlib.lines import Line2D
    b.legend([Line2D([], [], c="k", lw=1.5), Line2D([], [], c=col["hxoff"], lw=1.3), Line2D([], [], c=ns.PAL["theory"], ls="--", lw=1.3)],
             ["reference", r"$\nu_x=0$", r"equation~(\ref{eq:drainlaw})" if False else "the drain law"], loc="lower right", labelspacing=0.2)
    ns.style_axes(b)
    ns.tag(a, "(a)", x=0.86); ns.tag(b, "(b)")
    ns.save(fig, "drain")
    # ------------------------------------------------------------------ figure supply
    fig, (a, b) = ns.two_panel()
    E0 = C["ref"][:, 1].max()
    mk = {"ref": dict(marker="o", c="k", ms=4.2, mfc="k"), "hxoff": dict(marker="o", c=col["hxoff"], ms=4.2, mfc=col["hxoff"]),
          "hxy_d4": dict(marker="s", c=col["hxy_d4"], ms=4.6, mfc="w"), "nx128": dict(marker="^", c=col["nx128"], ms=5.0, mfc="w")}
    allr = []
    for r in ("ref", "hxy_d4", "nx128", "hxoff"):
        B = bins(C[r], *T_RANGE[r])
        pred = (nuL - B[:, 3]) * B[:, 1] + NU_H * STRENGTH[r] * B[:, 2]
        a.loglog(pred / E0, B[:, 4] / E0, ls="none", mew=1.0, **mk[r])
        ratio = B[:, 4] / pred; allr.append(ratio)
        print(f"supply {r:7s}: {len(B)} windows, measured/predicted jet loss: median {np.median(ratio):.3f}, range {ratio.min():.2f}-{ratio.max():.2f}; "
              f"hyperdiffusion share of the prediction {np.median(NU_H*STRENGTH[r]*B[:,2]/pred):.2f}")
    lim = (2e-7, 4e-5)
    a.loglog(lim, lim, "-", c=ns.PAL["theory"], lw=1.0, zorder=0)
    a.set_xlim(*lim); a.set_ylim(*lim)
    a.set_xlabel(r"$[(\nu_{\rm L}-\Gamma_E)E_{\rm nz}+\nu_{\rm h}E_{\rm zon}]/E_0$"); a.set_ylabel(r"$-\dot E_{\rm zon}/E_0\;(c_{\rm ref}/L_{\rm ref})$")
    a.legend([Line2D([], [], ls="none", mew=1.0, **mk[r]) for r in ("ref", "hxy_d4", "nx128", "hxoff")],
             ["reference", r"$\nu_\perp/4$", r"$\Delta x/2$", r"$\nu_x=0$"], loc="lower right", labelspacing=0.15, handletextpad=0.2)
    fig.delaxes(b); pa = a.get_position(); a.set_position([0.5 - pa.width / 2, pa.y0, pa.width, pa.height])
    ns.style_axes(a)
    ns.save(fig, "supply")


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        rebuild()
    draw()

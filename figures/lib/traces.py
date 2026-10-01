r"""Time traces of the reduced-box run and its amplitude-reduced restarts, shared by the time-trace figures
(history, amplitude, entropy, dissipation, remnant).  Read only; used by the figures' --rebuild step.

plunk_e_time.dat (0-based columns): 0 t, 10 E_zon, 11 E_nz (reduced electrostatic energy, zonal / non-zonal),
15 Z_nz (non-zonal entropy).  fe_time.dat: 0 t, 1 W, 2 drive, 3 velocity-space, 4 parallel, 5 perpendicular
hyperdiffusion (all <= 0), 6 collisions (0, collision_op = 'none').
"""
import os
import numpy as np

RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
LEGS = (1, 3, 4, 5, 6)                    # leg 2 ran no steps
T_RESTART = 5906.678243865176             # end of leg 3 = start of the amplitude tests
RUNS = {"orig": [f"{RB}/rb_c0p3_hxy/leg_{l:04d}/out" for l in LEGS],
        "a03": [f"{RB}/rb_c0p3_hxy_ampltest/amp0p3/out"],
        "a003": [f"{RB}/rb_c0p3_hxy_ampltest/amp0p03/out"]}
EPS = {"orig": 1.0, "a03": 0.3, "a003": 0.03}


def _splice(paths, name):
    """Concatenate the legs in time order, dropping any record that does not advance the time."""
    d = np.concatenate([np.loadtxt(os.path.join(p, name), comments="#") for p in paths])
    keep = np.concatenate([[True], np.maximum.accumulate(d[:-1, 0]) < d[1:, 0]])
    return d[keep]


def plunk(run):
    """t, E_zon, E_nz, Z_nz of a run ('orig', 'a03', 'a003')."""
    d = _splice(RUNS[run], "plunk_e_time.dat")
    return d[:, 0], d[:, 10], d[:, 11], d[:, 15]


def sinks(run="orig"):
    """t, W, drive, and the three sinks as positive rates (velocity-space, parallel, perpendicular)."""
    d = _splice(RUNS[run], "fe_time.dat")
    return d[:, 0], d[:, 1], d[:, 2], -d[:, 3], -d[:, 4], -d[:, 5], d


def runmean(t, y, width, log=False):
    """Running mean of y over a window of `width` time units centred on each sample (time based, so that the
    slightly different output cadence of the legs does not matter); NaN where the window leaves the record.
    log=True averages ln y (a geometric mean), which is the natural mean for a quantity decaying by decades."""
    z = np.log(y) if log else np.asarray(y, float)
    c = np.concatenate([[0.0], np.cumsum(z)])
    i0 = np.searchsorted(t, t - width / 2.0, side="left")
    i1 = np.searchsorted(t, t + width / 2.0, side="right")
    m = (c[i1] - c[i0]) / np.maximum(i1 - i0, 1)
    m[(t - width / 2.0 < t[0]) | (t + width / 2.0 > t[-1])] = np.nan
    return np.exp(m) if log else m


def sliding_rate(t, y, width, step, tmin=None, tmax=None):
    """Least-squares slope of ln y against t in windows of `width` every `step` time units.
    Returns the window centres, the decay rate -d ln y/dt, its standard error, and the geometric mean of y."""
    tmin = t[0] if tmin is None else tmin
    tmax = t[-1] if tmax is None else tmax
    tc, g, eg, ym = [], [], [], []
    a = tmin
    while a + width <= tmax + 1e-9:
        m = (t >= a) & (t < a + width)
        if m.sum() >= 10 and np.all(y[m] > 0):
            x, ly = t[m] - t[m].mean(), np.log(y[m])
            s = (x * ly).sum() / (x * x).sum()
            res = ly - ly.mean() - s * x
            # the breathing makes neighbouring residuals correlated; the error is the white-noise one and is a
            # lower bound
            err = np.sqrt((res ** 2).sum() / max(m.sum() - 2, 1) / (x * x).sum())
            tc.append(a + width / 2.0); g.append(-s); eg.append(err); ym.append(np.exp(ly.mean()))
        a += step
    return np.array(tc), np.array(g), np.array(eg), np.array(ym)


# ---- shared furniture of the time-axis panels (same in every figure) ------------------------------------------------
COLLAPSE = (10600.0, 11100.0)             # the collapse of the vortices: the last, fastest part of their decay
BAND = dict(color="0.87", lw=0, zorder=0)


def collapse_band(ax, scale=1e3):
    """Shade the collapse interval (the same pale grey band in every panel whose abscissa is t / scale)."""
    return ax.axvspan(COLLAPSE[0] / scale, COLLAPSE[1] / scale, **BAND)


def label(ax, x, y, s, c="k", **kw):
    """Direct label next to a curve, in the curve's colour, on a white patch (data coordinates)."""
    opts = dict(color=c, ha="left", va="center", zorder=7, bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    opts.update(kw)
    return ax.text(x, y, s, **opts)

#!/usr/bin/env python3
r"""Measurement B, exact version (no figure): the budget of the electrostatic (Debye) energy E_D by k_y, from GENE's own
mode-resolved rate diagnostics (Spectral_2D_EDrate<channel>_<species>.dat, written every 100 steps; text, one value
per line: header = k_x axis (NX0 - 1 values) + k_y axis, then per frame time, total, data[k_x, k_y]).

Conventions as ep_turbulence_paper/scripts_letter/lfig09_breathing: species summed, k_y > 0 counted twice,
E_D = Eplunk / 2 and its rates = EDrate<ch> / 2; parallel = parallel_raw - hyp_v - hyp_z (hyp_on_h = T).
Each channel is fitted in sub-windows to  const + slope tau + two harmonics of the breathing frequency, and the
constants are averaged: the cycle-averaged rate.  Everything is divided by the mean non-zonal E_D of the sub-window,
so the numbers are rates to be compared with Gamma_E = -d ln E_nz / dt.
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import vortex_field as vf

NKX, NKY, W = vf.NX - 1, vf.NKY, 24
HDR, FRL = NKX + NKY, 2 + NKX * NKY
CH = ("NL", "curvature", "drive", "parallel_raw", "hyp_v", "hyp_z", "hyp_kperp", "hypzcomp", "coll", "sources")


def frames(path, t0, t1):
    n = (os.path.getsize(path) // W - HDR) // FRL
    with open(path, "rb") as f:
        ts = np.empty(n)
        for i in range(n):
            f.seek((HDR + i * FRL) * W)
            ts[i] = float(f.read(W))
        idx = np.where((ts >= t0) & (ts <= t1))[0]
        if len(idx) == 0:
            return np.empty(0), np.empty((0, NKX, NKY))
        f.seek((HDR + idx[0] * FRL) * W)
        a = np.array(f.read((idx[-1] - idx[0] + 1) * FRL * W).split(), float).reshape(-1, FRL)
    return a[:, 0], a[:, 2:].reshape(-1, NKX, NKY)


def load(stem, t0, t1, run="orig"):
    """species-summed, k_x-summed, k_y-weighted series [nt, nky] over the legs of a run."""
    T, D, last = [], [], -np.inf
    for out in vf.RUNS[run]:
        tt, dd = None, 0.0
        for sp in ("electrons", "positrons"):
            p = f"{out}/Spectral_2D_{stem}_{sp}.dat"
            if not os.path.exists(p):
                return None, None
            tt, d = frames(p, t0, t1)
            dd = dd + d.sum(axis=1)
        m = tt > last + 1e-9
        if m.any():
            T.append(tt[m]); D.append(dd[m]); last = tt[m][-1]
    wgt = np.full(NKY, 2.0); wgt[0] = 1.0
    return np.concatenate(T), np.concatenate(D) * wgt


def fitmean(tau, y, wb, nh=2):
    M = [np.ones_like(tau), tau]
    for h in range(1, nh + 1):
        M += [np.cos(h * wb * tau), np.sin(h * wb * tau)]
    c, *_ = np.linalg.lstsq(np.stack(M, axis=1), y, rcond=None)
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--win", type=float, default=100.0)
    ap.add_argument("--run", default="orig")
    ap.add_argument("--bounce", action="store_true", help="the orbit-variance (bounce) part of E instead of the Debye part")
    ap.add_argument("--windows", type=float, nargs="*", default=[5910.0, 6710.0, 9430.0, 10430.0, 10430.0, 10830.0])
    a = ap.parse_args()
    tr = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trapping", "cache.npz"))[a.run]
    for (t0, t1) in zip(a.windows[0::2], a.windows[1::2]):
        t, E = load("Eplunk", t0, t1, a.run); E = 0.5 * E
        R = {}
        for ch in CH:
            tt, r = load(("EbounceRate" if a.bounce else "EDrate") + ch, t0, t1, a.run)
            if r is not None:
                assert np.allclose(tt, t)
                R[ch] = r if a.bounce else 0.5 * r
        if a.bounce:
            E = load("Ebounce", t0, t1, a.run)[1]
            En_ref = 0.5 * load("Eplunk", t0, t1, a.run)[1][:, 1:].sum(axis=1).mean()
        R["parallel"] = R.pop("parallel_raw") - R["hyp_v"] - R["hyp_z"]
        groups = {"zonal (ky=0)": [0], "vortex (ky=1)": [1], "ky>=2": list(range(2, NKY)), "non-zonal": list(range(1, NKY))}
        acc = {g: {ch: [] for ch in R} for g in groups}
        dE = {g: [] for g in groups}
        for tc in np.arange(t0 + a.win / 2, t1 - a.win / 2 + 1e-6, a.win):
            m = (t >= tc - a.win / 2) & (t < tc + a.win / 2)
            wb = np.interp(tc, tr[:, 0], tr[:, 1] - tr[:, 2])
            tau = t[m] - tc
            En = En_ref if a.bounce else E[m][:, 1:].sum(axis=1).mean()
            for g, ks in groups.items():
                for ch in R:
                    acc[g][ch].append(fitmean(tau, R[ch][m][:, ks].sum(axis=1), wb)[0] / En)
                dE[g].append(fitmean(tau, E[m][:, ks].sum(axis=1), wb)[1] / En)
        print(f"\n=== {t0:.0f}-{t1:.0f}: cycle-averaged rates of E_D by k_y group, divided by the non-zonal E_D (units 1e-4 c_ref/L_ref);"
              f"  E_nz/E_zon = {E[:, 1:].sum(axis=1).mean() / E[:, 0].mean():.3e},  E(ky=1)/E_nz = {E[:, 1].mean() / E[:, 1:].sum(axis=1).mean():.3f}")
        chs = list(R)
        print(f"{'group':15s}" + "".join(f"{c:>11s}" for c in chs) + f"{'sum':>11s}{'dE/dt':>11s}")
        for g in groups:
            v = [np.mean(acc[g][c]) * 1e4 for c in chs]
            e = [np.std(acc[g][c], ddof=1) / np.sqrt(len(acc[g][c])) * 1e4 for c in chs]
            print(f"{g:15s}" + "".join(f"{x:+11.3f}" for x in v) + f"{sum(v):+11.3f}{np.mean(dE[g]) * 1e4:+11.3f}")
            print(f"{'  +-':15s}" + "".join(f"{x:11.3f}" for x in e))


if __name__ == "__main__":
    main()

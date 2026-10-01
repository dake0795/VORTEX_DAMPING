#!/usr/bin/env python3
r"""Measurement B (no figure): is the electrostatic energy lost by the vortices handed to the jet?

Fluid (Euler) limit of the notes: d_t U = -d_x <v_x v_y>_y with v = (-d_y phi, d_x phi) / C_xy, so the zonal energy
E_U = (1/2) int U^2 dx changes through the Reynolds stress at the rate
    T = -int U d_x <v_x v_y> dx = int S <v_x v_y> dx,     <v_x v_y> = -(2 / C_xy^2) sum_{ky>0} k_y Im[phi_k^* d_x phi_k],
and r(t) = T / E_U is the rate of change of ln E_zon by the nonlinear exchange.  It is compared with the GENE energy
diagnostic y(t) = E_nz / E_zon (plunk_e_time.dat, same output times), for which the exchange alone gives
dy/dt = -r.  Both series are fitted in windows to  const + slope tau + harmonics of the breathing frequency w_b:
  * check of the method (normalisation, sign, flute assumption): the w_b component of r against that of -dy/dt;
  * the answer: the window mean of r against the secular -dy/dt (the decay of the vortices).
The potential is taken at the outboard mid-plane, or (--zavg) Jacobian-averaged along the field line.
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import vortex_field as vf


def series(t0, t1, zavg):
    if not zavg:
        return vf.midplane_series("orig", t0, t1)
    T, P, last = [], [], -np.inf
    for out in vf.RUNS["orig"]:
        geo = np.array(open(out + "/dipole_fix.dat").read().split("/")[1].split(), float).reshape(vf.NZ, -1)
        J = geo[:, 10] / geo[:, 10].sum()
        ts = vf.frame_times(out + "/field.dat")
        for i in np.where((ts >= t0) & (ts <= t1) & (ts > last + 1e-9))[0]:
            t, phi = vf.read_frame(out + "/field.dat", i)
            T.append(t); P.append(phi @ J)
        if len(ts):
            last = max(last, ts[-1])
    return np.array(T), np.array(P)


def transfer(p):
    """r = T / E_U for each frame; p[nt, nx, nky]."""
    U = vf.to_x(p[:, :, 0], nxf=256, deriv=1).real / vf.CXY
    S = vf.to_x(p[:, :, 0], nxf=256, deriv=2).real / vf.CXY
    ky = vf.KYMIN * np.arange(vf.NKY)
    f = vf.to_x(np.moveaxis(p[:, :, 1:], 1, 2), nxf=256)                 # [nt, nky-1, nxf]
    fx = vf.to_x(np.moveaxis(p[:, :, 1:], 1, 2), nxf=256, deriv=1)
    stress = -(2.0 / vf.CXY ** 2) * np.sum(ky[None, 1:, None] * np.imag(np.conj(f) * fx), axis=1)
    return (S * stress).mean(axis=1) / (0.5 * (U ** 2).mean(axis=1))


def fit(tau, y, wb, nh=2):
    M = [np.ones_like(tau), tau]
    for h in range(1, nh + 1):
        M += [np.cos(h * wb * tau), np.sin(h * wb * tau)]
    M = np.stack(M, axis=1)
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return c, y - M @ c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zavg", action="store_true")
    ap.add_argument("--win", type=float, default=100.0)
    a = ap.parse_args()
    d = np.concatenate([np.loadtxt(o + "/plunk_e_time.dat", comments="#") for o in vf.RUNS["orig"]])
    d = d[np.argsort(d[:, 0])]
    for (t0, t1) in ((5910.0, 6710.0), (9430.0, 10430.0)):
        t, p = series(t0, t1, a.zavg)
        r = transfer(p)
        y = np.interp(t, d[:, 0], d[:, 11] / d[:, 10])
        assert np.abs(np.interp(t, d[:, 0], d[:, 0]) - t).max() < 1e-6
        rows = []
        for tc in np.arange(t0 + a.win / 2, t1 - a.win / 2 + 1e-6, a.win):
            m = (t >= tc - a.win / 2) & (t < tc + a.win / 2)
            o = vf.window_analysis(t[m], p[m], np.array([0.85, -0.85]))
            wb = o["w"][0] - o["w"][1]
            tau = t[m] - tc
            cr, er = fit(tau, r[m], wb)
            cy, ey = fit(tau, y[m], wb)
            rb = cr[2] - 1j * cr[3]                    # r = Re[rb e^{+i wb tau}] with this convention
            yb = cy[2] - 1j * cy[3]
            rho = rb / (-1j * wb * yb)                 # r_b against the w_b component of -dy/dt
            rows.append([tc, wb, cr[0], -cy[1], cy[0], abs(rb), rho.real, rho.imag, er.std(), ey.std()])
        R = np.array(rows)
        print(f"\nwindow {t0:.0f}-{t1:.0f} ({'field-line average' if a.zavg else 'mid-plane'}), sub-windows of {a.win:.0f}:")
        print("   t_c     w_b     <r>        -dy/dt      y        |r_b|     rho=r_b/(-dy/dt)_b    <r>/(-dy/dt)   <r>/|r_b|")
        for q in R:
            print(f"{q[0]:7.0f} {q[1]:6.3f} {q[2]:+.3e} {q[3]:+.3e} {q[4]:.3e} {q[5]:.3e}   {q[6]:+.3f}{q[7]:+.3f}i        {q[2] / q[3]:+7.3f}       {q[2] / q[5]:+.1e}")
        F = R[:, 2] / R[:, 3]
        print(f"mean <r> = {R[:, 2].mean():+.3e} +- {R[:, 2].std(ddof=1) / np.sqrt(len(R)):.1e};  mean -dy/dt = {R[:, 3].mean():+.3e};"
              f"  ratio of means = {R[:, 2].mean() / R[:, 3].mean():+.3f};  mean ratio = {F.mean():+.3f} +- {F.std(ddof=1) / np.sqrt(len(F)):.3f} (scatter {F.std(ddof=1):.3f});"
              f"  |rho| = {np.hypot(R[:, 6], R[:, 7]).mean():.3f}, arg = {np.degrees(np.angle(R[:, 6] + 1j * R[:, 7])).mean():+.1f} deg")
        # whole-interval secular rate for comparison
        m = (d[:, 0] >= t0) & (d[:, 0] <= t1)
        g = -np.polyfit(d[m, 0], np.log(d[m, 11]), 1)[0]
        print(f"Gamma_E over the interval = {g:.3e};  Gamma_E * mean y = {g * (d[m, 11] / d[m, 10]).mean():.3e}")


if __name__ == "__main__":
    main()

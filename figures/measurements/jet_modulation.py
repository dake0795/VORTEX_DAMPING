#!/usr/bin/env python3
r"""Measurement A (no figure): how much does the breathing modulate the jet?

For a window of mid-plane data the k_y = 0 potential is fitted, mode by mode in k_x, to
    Phi(k_x, t) = p0 + p1 tau + b e^{-i w_b tau} + b' e^{+i w_b tau} (+ the same at 2 w_b),   tau = t - t_c,
with w_b = w_+ - w_- from the two-frequency fit of the k_y = 1 potential in the same window.  Reported: the amplitude
of the oscillation of U = d_x Phi / C_xy at w_b relative to the jet speed (maximum over x and at the critical layers),
of the shear S at the critical layers, the displacement of the critical layers dU / |S| against the separatrix
half-width w and the grid spacing, and E_nz / E_zon for comparison.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import vortex_field as vf

WIN = 120.0


def plunk():
    d = np.concatenate([np.loadtxt(o + "/plunk_e_time.dat", comments="#") for o in vf.RUNS["orig"]])
    d = d[np.argsort(d[:, 0])]
    return d[:, 0], d[:, 10], d[:, 11]


def main():
    tp, Ez, En = plunk()
    print("   t     w_b    E_nz/E_zon  max dU/Umax  rms dU/Umax  dU/(U-c scale) at crit  dS/S at crit  dx_c/rho  dx_c/w   dx_c/dx   2w_b: max dU/Umax   resid")
    for tc in (3000.0, 4500.0, 5850.0, 7500.0, 9000.0, 10000.0, 10500.0):
        t, p = vf.midplane_series("orig", tc - WIN / 2, tc + WIN / 2)
        o = vf.window_analysis(t, p, np.array([0.86, -0.86]) * (1.15 if tc < 4000 else 1.0))
        wb = o["w"][0] - o["w"][1]
        tau = t - t.mean()
        M = np.stack([np.ones_like(tau), tau, np.exp(-1j * wb * tau), np.exp(1j * wb * tau),
                      np.exp(-2j * wb * tau), np.exp(2j * wb * tau)], axis=1)
        z = p[:, :, 0]
        coef, *_ = np.linalg.lstsq(M, z, rcond=None)
        resid = np.sum(np.abs(z - M @ coef) ** 2) / np.sum(np.abs(z - M[:, :2] @ coef[:2]) ** 2)
        amp = lambda k, d: np.abs(vf.to_x(coef[k], deriv=d)) / vf.CXY + np.abs(vf.to_x(coef[k + 1], deriv=d)) / vf.CXY
        dU, dS, dU2 = amp(2, 1), amp(2, 2), amp(4, 1)
        U, S = o["U"], o["S"]
        Um = np.abs(U).max()
        cr = o["crit_p"] + o["crit_m"]
        jj = [int(round(q[0] / vf.LX * vf.NXF)) % vf.NXF for q in cr]
        dUc = np.mean([dU[j] for j in jj]); dSc = np.mean([dS[j] / abs(S[j]) for j in jj])
        dxc = np.mean([dU[j] / abs(S[j]) for j in jj])
        half = 0.5 * (o["half_p"] + o["half_m"])
        m = (tp >= tc - WIN / 2) & (tp <= tc + WIN / 2)
        print(f"{tc:7.0f} {wb:6.3f}  {np.mean(En[m] / Ez[m]):.3e}   {dU.max() / Um:.3e}   {np.sqrt(np.mean(dU ** 2)) / Um:.3e}"
              f"      {dUc / Um:.3e}           {dSc:.3e}    {dxc:7.2f}  {dxc / half:.3f}   {dxc / vf.DX:.3f}       {dU2.max() / Um:.2e}     {resid:.2f}")


if __name__ == "__main__":
    main()

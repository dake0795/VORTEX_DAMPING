#!/usr/bin/env python3
r"""Cross-check of mode_budget.py: the nonlinear zonal E_D rate rebuilt from the potential alone, z by z.

With the charge density proportional to k_perp^2 lambda_D^2 phi at each point along the field line, the E x B
nonlinearity gives, at each z, d_t U = -d_x <v_x v_y>_y (the metric drops out when g^xy = 0), so
    dE_D,zon/dt |_NL  =  sum_z J(z) g^xx(z) int S <v_x v_y> dx  /  normalisation,
compared here with E_D,zon = sum_z J g^xx (1/2) int U^2 dx.  Printed: the cycle-averaged r = T / E_zon from (i) the
mid-plane alone, (ii) the Jacobian-averaged (flute) potential, (iii) the z-resolved sum, per sub-window.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import vortex_field as vf
from reynolds_transfer import fit

NXF = 128
ky = vf.KYMIN * np.arange(vf.NKY)


def TU(p):
    """p[nx, nky] -> (int S tau dx, (1/2) int U^2 dx) per unit length."""
    U = vf.to_x(p[:, 0], nxf=NXF, deriv=1).real / vf.CXY
    S = vf.to_x(p[:, 0], nxf=NXF, deriv=2).real / vf.CXY
    f = vf.to_x(p[:, 1:].T, nxf=NXF); fx = vf.to_x(p[:, 1:].T, nxf=NXF, deriv=1)
    tau = -(2.0 / vf.CXY ** 2) * np.sum(ky[1:, None] * np.imag(np.conj(f) * fx), axis=0)
    return (S * tau).mean(), 0.5 * (U ** 2).mean()


def main():
    t0, t1, win = float(sys.argv[1]), float(sys.argv[2]), 100.0
    tr = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trapping", "cache.npz"))["orig"]
    T, R = [], []
    last = -np.inf
    for out in vf.RUNS["orig"]:
        geo = np.array(open(out + "/dipole_fix.dat").read().split("/")[1].split(), float).reshape(vf.NZ, -1)
        J, gxx = geo[:, 10] / geo[:, 10].sum(), geo[:, 0]
        ts = vf.frame_times(out + "/field.dat")
        for i in np.where((ts >= t0) & (ts <= t1) & (ts > last + 1e-9))[0]:
            t, phi = vf.read_frame(out + "/field.dat", i)
            a = TU(phi[:, :, vf.ZMID]); b = TU(phi @ J)
            c = np.array([TU(phi[:, :, z]) for z in range(vf.NZ)])
            w = J * gxx
            T.append(t); R.append([a[0] / a[1], b[0] / b[1], (w * c[:, 0]).sum() / (w * c[:, 1]).sum()])
        if len(ts):
            last = max(last, ts[-1])
    T, R = np.array(T), np.array(R)
    rows = []
    for tc in np.arange(t0 + win / 2, t1 - win / 2 + 1e-6, win):
        m = (T >= tc - win / 2) & (T < tc + win / 2)
        wb = np.interp(tc, tr[:, 0], tr[:, 1] - tr[:, 2])
        rows.append([fit(T[m] - tc, R[m, k], wb)[0][0] for k in range(3)] + [np.hypot(*fit(T[m] - tc, R[m, 2], wb)[0][2:4])])
    rows = np.array(rows)
    print(f"{t0:.0f}-{t1:.0f}: cycle-averaged T/E_zon:  mid-plane {rows[:, 0].mean():+.3e} +- {rows[:, 0].std(ddof=1) / np.sqrt(len(rows)):.1e};"
          f"  flute average {rows[:, 1].mean():+.3e} +- {rows[:, 1].std(ddof=1) / np.sqrt(len(rows)):.1e};"
          f"  z-resolved {rows[:, 2].mean():+.3e} +- {rows[:, 2].std(ddof=1) / np.sqrt(len(rows)):.1e};  oscillation amplitude {rows[:, 3].mean():.3e}")


if __name__ == "__main__":
    main()

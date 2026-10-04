"""Burst-rerun test (run after rb_c0p3_draintest/burst_rerun legA, lin_m1, lin_m3, legB):
 1. phase-mixing prediction of the final jets' non-flute factor from legA's checkpoint (t 700), m = 1 and m = 3
    (rh_predict.py: I = <h>_orbit - F0 phibar conserved by the linear k_y = 0 dynamics);
 2. the factor's evolution in the LINEAR runs (lin_m1, lin_m3: k_y = 0 only, from that checkpoint);
 3. the factor's evolution in the NONLINEAR continuation (legB) and in the original reference run at the same times.
If 1 = end of 2 ~ 3, the jets' non-flute factor is the phase-mixed remnant of the burst."""
import os, sys, subprocess, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
BR = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_draintest/burst_rerun"
for m in (1, 3):
    print(f"--- rh_predict, legA checkpoint, m = {m}")
    print(subprocess.run([sys.executable, os.path.join(HERE, "rh_predict.py"), f"{BR}/legA/out", str(m)], capture_output=True, text=True).stdout)
import bounce, response, vortex_field as vf
geo = bounce.geometry(f"{BR}/legA/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum(); Gxx = (wz * geo["gxx"]).sum(); nz = geo["nz"]
Rs = response.Response(geo); Rs.source = lambda k2: k2 * (geo["gxx"] - Gxx); chi = np.conj(Rs.solve(0.02, 1.0))
def factor(pk, kx):
    p0 = (pk * wz).sum(); pred = 5000.0 * kx ** 2 * chi * p0; return abs(np.vdot(pred * wz, pk - p0) / np.vdot(pred * wz, pred))
for m in (1, 3):
    out = f"{BR}/lin_m{m}/out"
    if not os.path.exists(out + "/field.dat"): print(f"lin_m{m}: no output"); continue
    N = 3 * nz; FR = 16 + 4 + N * 16 + 4; b = open(out + "/field.dat", "rb").read(); n = len(b) // FR
    kx = 2 * np.pi * m / vf.LX; rows = []
    for j in np.linspace(0, n - 1, 9).astype(int):
        f = b[j * FR:(j + 1) * FR]; t = np.frombuffer(f[4:12], "<f8")[0]
        pk = np.frombuffer(f[20:20 + N * 16], "<c16").reshape((3, nz), order="F")[1]; rows.append((t, factor(pk, kx)))
    print(f"lin_m{m} (linear):  " + "  ".join(f"t {t:.0f}: {x:.3f}" for t, x in rows))
import fieldgen as fg
for name, dirs in (("legB (nonlinear rerun)", [f"{BR}/legA/out", f"{BR}/legB/out"]), ("reference", [fg.RB + "/rb_c0p3_hxy/leg_0001/out"])):
    dirs = [d for d in dirs if os.path.exists(d + "/field.dat")]
    R = fg.Run(dirs, 64); rows = []
    for tc in (650, 700, 750, 800, 900, 1000, 1200, 1500):
        j = np.argmin(np.abs(R.t - tc))
        if abs(R.t[j] - tc) > 15: continue
        t, phi = R.frame(j); rows.append((t, factor(phi[1, 0, :], 2 * np.pi / vf.LX), factor(phi[3, 0, :], 6 * np.pi / vf.LX)))
    print(f"{name}: " + "  ".join(f"t {t:.0f}: m1 {a:.3f} m3 {c:.3f}" for t, a, c in rows))

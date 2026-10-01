#!/usr/bin/env python3
r"""Do the spectra and nonlinear times of the simulations follow the scalings of two-dimensional hydrodynamics?

For a run and a time window: the mid-plane potential phi(k_x, k_y) from field.dat (every STRIDE-th frame) gives the
E x B kinetic-energy spectrum per mode, e = (g^xx k_x^2 + g^yy k_y^2)|phi|^2 / (2 C_xy^2)  (k_y > 0 counted twice),
binned in shells of k = (g^xx k_x^2 + g^yy k_y^2)^{1/2} to E(k) (energy per unit k), with the field-line-averaged
metric.  From it: the local eddy rate [k^3 E(k)]^{1/2} and the strain of all larger scales [int_0^k k'^2 E dk']^{1/2}.
2D hydrodynamics: E ~ k^-5/3 and rate ~ k^{2/3} in the inverse energy range; E ~ k^-3 and a rate independent of k
in the forward enstrophy range.
Usage: hydro_scalings.py OUTDIR NX NKY LX KYMIN T0 T1 [STRIDE]
"""
import os, sys
import numpy as np
out, NX, NKY, LX, KYMIN, T0, T1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
STRIDE = int(sys.argv[8]) if len(sys.argv) > 8 else 5
NZ, CXY, GXX, GYY, LAMD = 32, 0.92640709665123333, 0.573, 0.956, 70.71
N = NX * NKY * NZ; FR = 16 + 4 + N * 16 + 4
p = os.path.join(out, "field.dat"); n = os.path.getsize(p) // FR
with open(p, "rb") as f:
    ts = np.empty(n)
    for i in range(n):
        f.seek(i * FR + 4); ts[i] = np.frombuffer(f.read(8), "<f8")[0]
    idx = np.where((ts >= T0) & (ts <= T1))[0][::STRIDE]
    P = np.zeros((NX, NKY)); Pz = np.zeros((NX, NKY))
    for i in idx:
        f.seek(i * FR + 20)
        phi = np.frombuffer(f.read(N * 16), "<c16").reshape((NX, NKY, NZ), order="F")
        P += np.abs(phi[:, :, NZ // 2]) ** 2 / len(idx)
kx = np.fft.fftfreq(NX, 1.0 / NX) * 2 * np.pi / LX; ky = np.arange(NKY) * KYMIN
KX, KY = np.meshgrid(kx, ky, indexing="ij")
k2 = GXX * KX ** 2 + GYY * KY ** 2; k = np.sqrt(k2)
w = np.where(KY > 0, 2.0, 1.0)
e = 0.5 * k2 * P * w / CXY ** 2                      # energy per mode (rho_ref c_ref / L_ref)^2
ez = e[:, 0].sum(); enz = e[:, 1:].sum()
print(f"{out}\n frames {len(idx)}, t {ts[idx[0]]:.0f}-{ts[idx[-1]]:.0f};  u_rms: zonal {np.sqrt(2*ez):.1f}, non-zonal {np.sqrt(2*enz):.1f} (thermal drift speed ~ 3);  E_nz/E_zon {enz/ez:.3e}")
for name, mask in (("non-zonal (k_y >= 1)", KY > 0), ("zonal (k_y = 0)", KY == 0)):
    kk = k[mask]; ee = e[mask]
    edges = np.geomspace(kk[kk > 0].min() * 0.999, kk.max() * 1.001, 15)
    kc = np.sqrt(edges[1:] * edges[:-1]); E = np.array([ee[(kk >= a) & (kk < b)].sum() / (b - a) for a, b in zip(edges[:-1], edges[1:])])
    ok = E > 0
    strain = np.sqrt(np.cumsum(np.where(ok, kc ** 2 * E * np.diff(edges), 0.0)))
    local = np.sqrt(kc ** 3 * E)
    sl = np.gradient(np.log(np.where(ok, E, np.nan)), np.log(kc))
    print(f" {name}:   k*lambda_D    E(k)        local slope   [k^3E]^1/2   strain of larger scales")
    for j in range(len(kc)):
        if ok[j]:
            print(f"            {kc[j]*LAMD:8.3f}   {E[j]:10.3e}   {sl[j]:7.2f}      {local[j]:8.3f}     {strain[j]:8.3f}")
    for lo, hi in ((0.03, 0.2), (0.2, 1.0), (1.0, 4.0)):
        m = ok & (kc * LAMD >= lo) & (kc * LAMD < hi)
        if m.sum() >= 3:
            print(f"   fitted slope of E(k) for {lo} < k lambda_D < {hi}: {np.polyfit(np.log(kc[m]), np.log(E[m]), 1)[0]:.2f}")

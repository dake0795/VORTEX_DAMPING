"""Euler invariants of the flute potential through each run: energy and enstrophy of the zonal and non-zonal parts, and
the enstrophy removed by the perpendicular hyperdiffusion, from field.dat (Jacobian-weighted field-line average, metric
<g^xx>, <g^yy>).  psi = phi / C_xy; per mode e = k'^2 |psi|^2, w = k'^4 |psi|^2, k'^2 = <g^xx> k_x^2 + <g^yy> k_y^2
(k_y > 0 twice); hyperdiffusive rate of mode k [GENE, on h: charge moment (1 + lambda_D^2 k^2) phi]:
nu_k = hx (dx k_x/2)^4 + hy (dy k_y/2)^4, damping of the vorticity at nu_k (1 + 1/(lambda_D^2 k_perp^2)).
Usage: invariants.py NAME NX HX HY STRIDE DIR [DIR ...]  ->  inv_NAME.npz"""
import sys, os, numpy as np
name, NX, HX, HY, ST = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
dirs = sys.argv[6:]
NKY, NZ, LX, KYMIN, CXY = 16, 32, 3689.25, 0.002, 0.92640709665123333
LAM2 = 1.0e4 / 2.0
g = np.loadtxt(os.path.join(dirs[0], "dipole_fix.dat"), skiprows=19)
w = g[:, 10] / g[:, 10].sum(); gxx, gyy = (w * g[:, 0]).sum(), (w * g[:, 3]).sum()
kx = np.fft.fftfreq(NX, 1.0 / NX) * 2 * np.pi / LX; ky = np.arange(NKY) * KYMIN
KX, KY = np.meshgrid(kx, ky, indexing="ij")
k2 = gxx * KX ** 2 + gyy * KY ** 2; kp2 = KX ** 2 + KY ** 2
dx, dy = LX / NX, (2 * np.pi / KYMIN) / (2 * NKY)
nu = HX * (dx * KX / 2) ** 4 + HY * (dy * KY / 2) ** 4
fac = np.where(kp2 > 0, 1 + 1 / (LAM2 * np.where(kp2 > 0, kp2, 1)), 0)
wt = np.where(KY > 0, 2.0, 1.0)
N = NX * NKY * NZ; FR = 16 + 4 + N * 16 + 4
rows = []
for d in dirs:
    p = os.path.join(d, "field.dat"); n = os.path.getsize(p) // FR
    with open(p, "rb") as f:
        for i in range(0, n, ST):
            f.seek(i * FR + 4); t = np.frombuffer(f.read(8), "<f8")[0]
            f.seek(i * FR + 20)
            phi = np.frombuffer(f.read(N * 16), "<c16").reshape((NX, NKY, NZ), order="F")
            P = np.abs((phi * w).sum(-1) / CXY) ** 2 * wt
            e, ens = k2 * P, k2 ** 2 * P
            dens = 2 * nu * fac * ens                       # enstrophy removed per unit time
            z = KY == 0
            # fine zonal structure: |k_x| lambda_D > 0.5
            fz = z & (np.abs(KX) * np.sqrt(LAM2) > 0.5)
            rows.append([t, e[z].sum(), e[~z].sum(), ens[z].sum(), ens[~z].sum(), dens[z].sum(), dens[~z].sum(), ens[fz].sum(), (2 * nu * fac * e)[z].sum()])
    print(d, n, flush=True)
a = np.array(rows); a = a[np.argsort(a[:, 0])]
np.savez(f"inv_{name}.npz", a=a, cols="t Ezon Enz Wzon Wnz DWzon_h DWnz_h Wzon_fine DEzon_h")
print("wrote", name, a.shape)

"""Shared constants and readers for the fluid-limit (Rayleigh / 2D Euler) calculations of the notes (1 Oct 2026).

Fluid limit used here (notes, limit (b)): for a flute-like potential the charge moment of the transit-averaged
equations gives, at leading order in k_perp lambda_D,
    d_t zeta + (1/C_xy) [d_x phi d_y zeta - d_y phi d_x zeta] = sinks,   zeta = (<g^xx> d_x^2 + <g^yy> d_y^2) phi,
where <.> is the Jacobian-weighted average along the field line (the weight that turns the velocity integral at fixed
l into the orbit-space integral), x, y are GENE's flux-tube coordinates in rho_ref and t is in L_ref/c_ref.
"""
import os
import numpy as np

RB = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
NX, NKY, NZ = 64, 16, 32
LX, KYMIN = 3689.25, 0.002
LY = 2 * np.pi / KYMIN
KXMIN = 2 * np.pi / LX
CXY = 0.92640709665123333
HYP = 0.0024                       # hyp_x = hyp_y of the run
DX_G, DY_G = LX / NX, LY / (2 * NKY)   # GENE's deli, delj (dy = ly/2/nky0)
FR = 4 + 8 + 4 + 4 + NX * NKY * NZ * 16 + 4


def geometry(leg=5):
    """Columns of GENE's out/dipole_fix.dat: gxx gxy gxz gyy gyz gzz B dBdx dBdy dBdz jacobian ..."""
    L = open(f"{RB}/rb_c0p3_hxy/leg_{leg:04d}/out/dipole_fix.dat").read().split("\n")
    i0 = [i for i, l in enumerate(L) if l.strip() == "/"][0]
    return np.array([[float(x) for x in l.split()] for l in L[i0 + 1:] if l.strip()])


def metric():
    d = geometry()
    w = d[:, 10] / d[:, 10].sum()
    return dict(w=w, gxx=float((w * d[:, 0]).sum()), gyy=float((w * d[:, 3]).sum()), gxx_mid=d[16, 0], gyy_mid=d[16, 3])


def nframes(path):
    return os.path.getsize(path) // FR


def read(path, i):
    with open(path, "rb") as f:
        f.seek(i * FR)
        b = f.read(FR)
    t = np.frombuffer(b[4:12], "<f8")[0]
    phi = np.frombuffer(b[20:20 + NX * NKY * NZ * 16], "<c16").reshape((NX, NKY, NZ), order="F")
    return t, phi


def times(path):
    n = nframes(path)
    ts = np.empty(n)
    with open(path, "rb") as f:
        for i in range(n):
            f.seek(i * FR + 4)
            ts[i] = np.frombuffer(f.read(8), "<f8")[0]
    return ts

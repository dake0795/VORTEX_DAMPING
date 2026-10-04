"""Extend the theoretical dissipation function D(wt) beyond |wt| = 4.2 (for the harmonics k_y = 2, 3 of the vortices)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat")
R = response.Response(geo, nmax=24)
wz = geo["J"] / geo["J"].sum(); Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
s = Gyy * (geo["gxx"] / Gxx - geo["gyy"] / Gyy)
C = np.load(os.path.join(HERE, "cache.npz"))
wg = np.r_[C["wg"], np.arange(4.4, 12.01, 0.2)]
D = np.array([-w * (wz * np.imag(R.solve(w, 1.0)) * s).sum() for w in wg])
print("check overlap with the stored table (nmax 24 vs default):", np.round(D[:66:8], 4), np.round(C["Dth"][::8], 4))
print("extended:", np.round(D[66::5], 4), "at", wg[66::5])
np.savez(os.path.join(HERE, "cache_Dext.npz"), wg=wg, Dth=D)

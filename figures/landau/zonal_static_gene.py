"""Linear GENE test (rb_c0p3_draintest/zonal_static/m1_chpt: nx0 3, seeded with the m = 1 jet of the rb_c0p3_hxoff checkpoint at t 9360.8, k_y = 0 only, reduced-box geometry and velocity grid):
non-flute part of one zonal harmonic (k_x lambda_D = 0.12) against the static kinetic response used in jets_nonflute.py.
In the nonlinear runs the jets' non-flute part is the predicted shape times ~1.6; is that factor already in linear GENE?"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response
OUT = f"{vf.RB}/rb_c0p3_draintest/zonal_static/m1_chpt/out"
NX, NZ = 3, 32; N = NX * NZ; FR = 16 + 4 + N * 16 + 4
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx = (wz * geo["gxx"]).sum(); L2 = 5000.0
R = response.Response(geo); R.source = lambda k2: k2 * (geo["gxx"] - Gxx)
chi = np.conj(R.solve(0.02, 1.0))
kx = 2 * np.pi / vf.LX
b = open(f"{OUT}/field.dat", "rb").read(); nfr = len(b) // FR
print(f"{nfr} frames; frame size check {len(b) % FR}")
for j in list(range(0, nfr, max(1, nfr // 12))) + [nfr - 1]:
    f = b[j * FR:(j + 1) * FR]; t = np.frombuffer(f[4:12], "<f8")[0]
    phi = np.frombuffer(f[20:20 + N * 16], "<c16").reshape((NX, NZ), order="F")[1]
    p0 = (phi * wz).sum(); p1 = phi - p0
    pred = L2 * kx ** 2 * chi * p0
    c = np.vdot(pred * wz, p1) / np.vdot(pred * wz, pred)
    res = np.sqrt(np.sum(wz * np.abs(p1 - c * pred) ** 2) / np.sum(wz * np.abs(p1) ** 2))
    print(f"t {t:7.2f}  |phi0| {abs(p0):.3e}  |phi1|/|phi0| {np.sqrt(np.sum(wz*np.abs(p1)**2))/abs(p0):.3e}  meas/pred {abs(c):.3f} at {np.degrees(np.angle(c)):+5.0f}  shape {1-res**2:.3f}")

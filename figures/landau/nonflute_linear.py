import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce
P = f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/field.dat"
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gyy = (wz * geo["gyy"]).sum(); Gxx = (wz * geo["gxx"]).sum()
shapes = {"gyy": geo["gyy"] - Gyy, "gxx": geo["gxx"] - Gxx, "B": geo["B"] - (wz * geo["B"]).sum(), "J": geo["J"] - (wz * geo["J"]).sum()}
for k in shapes: shapes[k] = shapes[k] - (wz * shapes[k]).sum()
ts = vf.frame_times(P)
for tsel in (250.0, 350.0, 450.0, 550.0):
    i = np.argmin(abs(ts - tsel)); t, phi = vf.read_frame(P, i)
    print(f"\nt {t:.1f} (linear phase from noise, no jets yet)")
    print("  ky  ky lamD   |nonflute|/|flute|   proj on gyy: real  imag    (gxx)     (B)    best-shape share")
    for ky in (1, 2, 3, 4, 6, 8, 12):
        f = phi[:, ky, :]; f0 = (f * wz[None, :]).sum(1); f1 = f - f0[:, None]
        tot = np.sum(np.abs(f1) ** 2 * wz[None, :]); nf = np.sqrt(tot / np.sum(np.abs(f0) ** 2))
        out = []
        for nm in ("gyy", "gxx", "B"):
            sh = shapes[nm]; A = (f1 * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum()
            c = np.vdot(f0, A) / np.vdot(f0, f0)        # A = c * phi0 best fit
            out.append(c)
        sh = shapes["gyy"]; A = (f1 * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum()
        share = 1 - np.sum(np.abs(f1 - A[:, None] * sh[None, :]) ** 2 * wz[None, :]) / tot
        print(f"  {ky:2d}  {ky*vf.KYMIN*np.sqrt(5000):5.2f}     {nf:.3e}          {out[0].real:+.4f} {out[0].imag:+.4f}   {abs(out[1]):.4f}  {abs(out[2]):.4f}   {share:.3f}")

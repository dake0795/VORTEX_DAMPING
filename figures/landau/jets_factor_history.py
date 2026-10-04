"""History of the jets' m = 1, 3 harmonics through and after the burst (rb_c0p3_hxy): flute amplitude |phi_0|,
non-flute amplitude |phi_1|, and the factor |phi_1| / (static response to phi_0) - is the rise of the factor at t 600-800
a change in phi_1 or in phi_0?  Also the vortex energy (k_y = 1) for timing."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, response, predict
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx = (wz * geo["gxx"]).sum(); L2 = 5000.0
R = response.Response(geo); R.source = lambda k2: k2 * (geo["gxx"] - Gxx)
chi = np.conj(R.solve(0.02, 1.0)); nchi = np.sqrt(np.sum(wz * np.abs(chi) ** 2))
kxg = np.fft.fftfreq(vf.NX, 1.0 / vf.NX) * 2 * np.pi / vf.LX
print("   t     m=1: |phi0|     |phi1|    factor  |   m=3: |phi0|     |phi1|    factor  |  |phi(ky=1)| flute rms")
for tc in np.arange(550, 1301, 50):
    t, p0, p1 = predict.window_field("orig", float(tc), half=10.0) if "half" in predict.window_field.__code__.co_varnames else predict.window_field("orig", float(tc))
    pz = p0.mean(0); pz0 = (pz * wz[None, :]).sum(1); pz1 = pz - pz0[:, None]
    row = f"{tc:5.0f}"
    for m in (1, 3):
        a1 = np.sqrt(np.sum(wz * np.abs(pz1[m]) ** 2)); a0 = abs(pz0[m])
        row += f"   {a0:10.3e} {a1:10.3e}  {a1/(L2*kxg[m]**2*nchi*a0):6.3f}  |"
    v = np.sqrt(np.mean(np.abs((p1 * wz).sum(-1)) ** 2))
    print(row + f"   {v:10.3e}")

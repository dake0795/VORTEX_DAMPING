import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib")); sys.path.insert(0, HERE)
import vortex_field as vf, bounce, predict
vf.RUNS["hxoff"] = [f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out"]
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxoff/leg_0001/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx, Gyy = (wz * geo["gxx"]).sum(), (wz * geo["gyy"]).sum()
T = np.load(os.path.join(HERE, "chi_tables.npz")); WT = T["WT"]
R2 = np.load(os.path.join(HERE, "..", "..", "scripts", "rayleigh", "cache", "zforce2_hxoff.npz"))["rows"]
NXF = predict.NXF
tc = 2700.0
t, p0, p1 = predict.window_field("hxoff", tc)
w0 = R2[np.argmin(abs(R2[:, 0] - tc)), 1]; w, _, _ = vf.two_freq_fit(t, p1[:, :, vf.ZMID], (w0, -w0))
Mx = np.exp(-1j * np.outer(t - t.mean(), w)); a, *_ = np.linalg.lstsq(Mx, p1.reshape(len(t), -1), rcond=None); a = a.reshape(2, vf.NX, vf.NZ)
U = vf.to_x((p0.mean(axis=0) * wz[None, :]).sum(axis=1), nxf=NXF, deriv=1).real / vf.CXY
# zonal non-flute part of the jets themselves
z0 = predict.to_x(p0.mean(axis=0)); z0f = (z0 * wz[None, :]).sum(1)
print("jets: |non-flute|/|flute| =", np.linalg.norm((z0 - z0f[:, None]) * np.sqrt(wz)) / np.linalg.norm(z0f))
shapes = {"gxx-<gxx>": geo["gxx"] - Gxx, "gyy-<gyy>": geo["gyy"] - Gyy, "s(l)": Gyy * (geo["gxx"] / Gxx - geo["gyy"] / Gyy), "B-<B>": geo["B"] - (wz * geo["B"]).sum() if "B" in geo else None}
k = 0
ph = predict.to_x(a[k]); ph0 = (ph * wz[None, :]).sum(1); meas = ph - ph0[:, None]
tot = np.sum(np.abs(meas) ** 2 * wz[None, :])
for nm, sh in shapes.items():
    if sh is None: continue
    sh = sh - (wz * sh).sum()
    A = (meas * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum()
    res = np.sum(np.abs(meas - A[:, None] * sh[None, :]) ** 2 * wz[None, :]) / tot
    print(f"z-shape {nm:10s}: fraction of non-flute variance explained by A(x) * shape(z): {1 - res:.3f}")
# z-structure of the measured non-flute part: SVD
Uu, S, Vh = np.linalg.svd(meas * np.sqrt(wz)[None, :], full_matrices=False)
print("SVD energy fractions:", np.round(S[:5] ** 2 / np.sum(S ** 2), 3))
v1 = Vh[0] / np.sqrt(wz); print("leading z-shape (normalised) vs gxx-<gxx>, gyy-<gyy>, s:",
    [round(abs(np.vdot(v1 * wz, sh - (wz * sh).sum())) / np.sqrt(np.vdot(v1 * wz, v1).real * (wz * (sh - (wz * sh).sum()) ** 2).sum()), 3) for sh in (geo["gxx"], geo["gyy"], shapes["s(l)"])])
xf = np.arange(NXF) * vf.LX / NXF
amp = np.sqrt((np.abs(meas) ** 2 * wz[None, :]).sum(1)); wt = w[k] - vf.KYMIN * U
print("x where measured non-flute amplitude peaks (top 6) and wt there, |phi0| there / max:")
for j in np.argsort(amp)[::-1][:6]: print(f"   x {xf[j]:7.1f}  amp {amp[j]:.3e}  wt {wt[j]:+.3f}  |phi0| {abs(ph0[j])/abs(ph0).max():.2f}")
print("critical layers (wt = 0):", xf[np.where(np.diff(np.sign(wt)) != 0)[0]])
sh = geo["gyy"] - Gyy; sh = sh - (wz * sh).sum()
A = (meas * wz[None, :] * sh[None, :]).sum(1) / (wz * sh ** 2).sum()          # measured amplitude of the gyy-shape, per x
Ts = T["chiy"]                                                                  # response to source gyy - <gyy>, per unit lambda^2
def proj(c): return (c * wz * sh).sum() / (wz * sh ** 2).sum()
pc = np.array([proj(c) for c in Ts])
L2 = 5000.0
print("\n  x       wt      |phi0|/max   measured A/phi0      predicted -k^2 lam^2 chi_y(wt) [proj]     ratio meas/pred")
d2 = vf.to_x((a[k] * wz[None, :]).sum(axis=1), nxf=NXF, deriv=2)
Tx = T["chix"]; pcx = np.array([proj(c) for c in Tx])
for j in range(0, NXF, 8):
    if abs(ph0[j]) < 0.2 * abs(ph0).max(): continue
    py = -L2 * vf.KYMIN ** 2 * np.interp(wt[j], WT, pc.real) - 1j * L2 * vf.KYMIN ** 2 * np.interp(wt[j], WT, pc.imag)
    px = L2 * (d2[j] / ph0[j]) * (np.interp(wt[j], WT, pcx.real) + 1j * np.interp(wt[j], WT, pcx.imag))
    r = A[j] / ph0[j]
    print(f"{xf[j]:7.1f}  {wt[j]:+.3f}   {abs(ph0[j])/abs(ph0).max():.2f}      {abs(r):.3e} {np.degrees(np.angle(r)):+5.0f}     y-term {abs(py):.3e}  x-term {abs(px):.3e}  sum {abs(py+px):.3e} {np.degrees(np.angle(py+px)):+5.0f}   {abs(r)/abs(py+px):.2f}")
Tsn = T["chis"] if "chis" in T.files else None

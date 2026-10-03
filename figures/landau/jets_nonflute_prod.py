import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bounce, response
D = sys.argv[1]; NX, NKY, NZ = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); LX = float(sys.argv[5]); T0, T1 = float(sys.argv[6]), float(sys.argv[7])
geo = bounce.geometry(D + "/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum(); Gxx = (wz * geo["gxx"]).sum(); L2 = 5000.0
R = response.Response(geo); R.source = lambda k2: k2 * (geo["gxx"] - Gxx); chi = np.conj(R.solve(0.02, 1.0))
N = NX * NKY * NZ; FR = 16 + 4 + N * 16 + 4; p = D + "/field.dat"; n = os.path.getsize(p) // FR
acc = 0; cnt = 0
with open(p, "rb") as f:
    for i in range(0, n, 10):
        f.seek(i * FR + 4); t = np.frombuffer(f.read(8), "<f8")[0]
        if t < T0 or t > T1: continue
        f.seek(i * FR + 20)
        # read only the k_y = 0 column for all kx, z: Fortran order (NX, NKY, NZ) -> stride
        phi = np.frombuffer(f.read(N * 16), "<c16").reshape((NX, NKY, NZ), order="F")[:, 0, :]
        acc = acc + phi; cnt += 1
pz = acc / cnt; pz0 = (pz * wz[None, :]).sum(1); pz1 = pz - pz0[:, None]
kx = np.fft.fftfreq(NX, 1.0 / NX) * 2 * np.pi / LX
print(f"{D}: {cnt} frames, nz {NZ}")
print("  m  kx lamD  |phi_z1|/|phi_z0|  measured/predicted   shape")
order = np.argsort(-np.abs(pz0[1:NX//2]))[:6] + 1
for m in sorted(order):
    pred = L2 * kx[m] ** 2 * chi * pz0[m]; meas = pz1[m]
    c = np.vdot(pred * wz, meas) / np.vdot(pred * wz, pred); res = np.sqrt(np.sum(wz * np.abs(meas - c * pred) ** 2) / np.sum(wz * np.abs(meas) ** 2))
    print(f" {m:3d}  {kx[m]*np.sqrt(L2):5.2f}   {np.sqrt(np.sum(wz*np.abs(meas)**2))/abs(pz0[m]):.3e}       {abs(c):.3f} at {np.degrees(np.angle(c)):+4.0f}     {1-res**2:.3f}")

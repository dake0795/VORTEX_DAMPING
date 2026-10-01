#!/usr/bin/env python3
r"""Test of the drift-resonance calculation at the critical layers against the measured loss of E_nz to entropy.

Prediction (notes, eq:gamma_layer): the rate at which the magnetic-drift term converts the electrostatic energy of a
flute-like vortex into entropy, per unit of that energy,
    Gamma_d = (2 pi <w_d (w_d - w_*)> / (lambda_D^2 k)) * sum_c (|phi_c|^2 / |S_c|) / int dx (gxx |phi'|^2 + gyy k^2 |phi|^2),
with w_d = k D (mu B + 2 v_par^2), D = <|K_y|/B> along the field line, w_* = (k/C_xy) [omn_eff + omt_eff (v^2 - 3/2)], so
    <w_d (w_d - w_*)> = k^2 D [7 D - 2 (omn_eff + omt_eff)/C_xy]
for a Maxwellian, where omn_eff, omt_eff are the gradients AT THE CRITICAL LAYER, as modified by the zonal state:
    omn_eff = omn - d_x <dn>_zonal,   omt_eff = omt - d_x <dT>_zonal   (GENE units; dn is the moment of delta f, which is the
non-Boltzmann zonal distribution delta F of the notes).
Inputs measured here: the vortex profile and the shear (geometry cache), the zonal density and temperature (mom files).
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import vortex_field as vf

OMN = OMT = 2.7816
NX, NKY, NZ = vf.NX, vf.NKY, vf.NZ
PL = NX * NKY * NZ * 16
FM = 16 + 6 * (4 + PL + 4)
out = f"{vf.RB}/rb_c0p3_hxy/leg_0004/out"
lines = open(f"{out}/dipole_fix.dat").read().split("\n/\n", 1)[1]
geo = np.array([[float(v) for v in l.split()] for l in lines.strip().splitlines()])
gxx, gyy, B, dBdx, J = geo[:, 0], geo[:, 3], geo[:, 6], geo[:, 7], geo[:, 10]
wz = J / J.sum()
D = (wz * np.abs(dBdx) / vf.CXY / B).sum()
GXX, GYY = (wz * gxx).sum(), (wz * gyy).sum()
print(f"D = <|K_y|/B> = {D:.3f};  <gxx> = {GXX:.3f}, <gyy> = {GYY:.3f}; mid-plane gxx {gxx[16]:.3f} gyy {gyy[16]:.3f}")


def zonal_moments(nframes=10):
    """Field-line (Jacobian) averaged zonal density and temperature of each species, averaged over the first frames."""
    res = {}
    for sp in ("positrons", "electrons"):
        p = f"{out}/mom_{sp}.dat"
        n = min(nframes, os.path.getsize(p) // FM)
        acc = np.zeros((3, NX), complex); tt = []
        with open(p, "rb") as f:
            for i in range(n):
                f.seek(i * FM + 4); tt.append(np.frombuffer(f.read(8), "<f8")[0])
                for k in range(3):
                    f.seek(i * FM + 16 + k * (4 + PL + 4) + 4)
                    a = np.frombuffer(f.read(PL), "<c16").reshape((NX, NKY, NZ), order="F")[:, 0, :]
                    acc[k] += (a * wz[None, :]).sum(axis=1) / n
        res[sp] = (acc[0], (acc[1] + 2 * acc[2]) / 3.0)
        print(sp, "frames", n, "t", tt[0], "..", tt[-1])
    return res


C = np.load(os.path.join(HERE, "..", "geometry", "cache.npz"))
k = vf.KYMIN
U = vf.to_x(C["p0"], deriv=1).real / vf.CXY
S = vf.to_x(C["p0"], deriv=2).real / vf.CXY
mom = zonal_moments()
dn = {sp: vf.to_x(m[0], deriv=1).real for sp, m in mom.items()}
dT = {sp: vf.to_x(m[1], deriv=1).real for sp, m in mom.items()}
x = vf.XF
for iv, name in ((0, "+"), (1, "-")):
    c = C["w"][iv] / k
    ph = vf.to_x(C["a"][iv]); dph = vf.to_x(C["a"][iv], deriv=1)
    den = np.trapz(GXX * np.abs(dph) ** 2 + GYY * k ** 2 * np.abs(ph) ** 2, x)
    idx = np.where(np.diff(np.sign(U - c)) != 0)[0]
    tot_flat = tot_eff = tot_unfl = 0.0
    for j in idx:
        geomf = np.abs(ph[j]) ** 2 / abs(S[j]) / den
        gn = np.mean([dn[sp][j] for sp in dn]); gT = np.mean([dT[sp][j] for sp in dT])
        omn_e, omt_e = OMN - gn, OMT - gT
        m_eff = k ** 2 * D * (7 * D - 2 * (omn_e + omt_e) / vf.CXY)
        m_flat = k ** 2 * D * 7 * D
        m_unfl = k ** 2 * D * (7 * D - 2 * (OMN + OMT) / vf.CXY)
        pref = 2 * np.pi / (vf.LAMBDA_D ** 2 * k)
        tot_eff += pref * m_eff * geomf; tot_flat += pref * m_flat * geomf; tot_unfl += pref * m_unfl * geomf
        print(f"vortex {name}: layer x = {x[j]:7.1f}  S = {S[j]:+.2f}  |phi_c|^2/(|S| int) = {geomf:.3e}   d_x dn = {gn:+.3f} (p {dn['positrons'][j]:+.3f}, e {dn['electrons'][j]:+.3f})"
              f"  d_x dT = {gT:+.3f}  -> omn_eff {omn_e:.2f} omt_eff {omt_e:.2f}  <wd(wd-w*)>/k^2 = {m_eff / k**2:+.2f}")
    print(f"vortex {name}: predicted Gamma_d = {tot_eff:.2e} (measured gradients at the layers), {tot_flat:.2e} (fully flattened), {tot_unfl:.2e} (unflattened);  measured loss to entropy ~ 2.0e-3")

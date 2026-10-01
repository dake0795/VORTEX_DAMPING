"""Check of the magnetic-drift frequency in GENE units (the notes owe this: footnote of 'which limit applies').
GENE (gene-stress/src: dgdxy_terms.F90, prefactors.F90 get_curv, geometry.F90 set_curvature) advances
    d_t g = ... - (T_a/q_a) (mu B + 2 v_par^2)/B * K_y * d_y g,      K_y = (dB/dx)/C_xy   (dipole_fix: g^xy = g^yz = 0),
with v_par in units of sqrt(2T/m), mu B in units of T, lengths in rho_ref, time in L_ref/c_ref.  A thermal particle has
<v_par^2> = 1/2 and <mu B> = 1, i.e. mu B + 2 v_par^2 = 2 on average and = 2 u for an isotropic particle of energy u T.
"""
import numpy as np
from common import *

d = geometry()
B, dBdx, J = d[:, 6], d[:, 7], d[:, 10]
w = J / J.sum()
Ky = dBdx / CXY
vd_loc = 2.0 * np.abs(Ky) / B                     # thermal drift speed at each z (rho_ref c_ref/L_ref)
print("thermal drift speed |v_d|: outboard mid-plane %.3f, inboard %.3f, Jacobian-weighted mean %.3f" % (vd_loc[16], vd_loc[0], (w * vd_loc).sum()))
# orbit averages: deeply trapped (v_par = 0 at the mid-plane, energy u): v_d = u K_y/B there; a passing particle with
# lambda = 0 (all parallel): v_d = 2 u <K_y/B> weighted by dl/v_par ~ dz J B
print("deeply trapped, energy u T: v_d = %.3f u;  lambda = 0 passing: v_d = %.3f u" % (abs(Ky[16]) / B[16], 2 * (w * B * np.abs(Ky) / B).sum() / (w * B).sum()))
for name, v in (("mid-plane thermal", vd_loc[16]), ("field-line mean thermal", (w * vd_loc).sum())):
    print(f"k_y v_d ({name}) = {KYMIN * v:.2e}")
omn = 2.7816
print("k_y v_* = k_y omn (T/q) / C_xy = %.2e  [gradient drive, per unit of omn = omt]" % (KYMIN * omn / CXY))
m = np.load(f"{CACHE}/measured_leg4_0.npz")
from measure import fine
U = fine(m["psi0"], 2048, 1).real; S = fine(m["psi0"], 2048, 2).real
c = float(m["w"]) / KYMIN
print("k_y (U_jet - c) = %.3f;  k_y U_jet = %.3f;  |S| at the layers ~ %.2f" % (KYMIN * (U.max() - c), KYMIN * U.max(), 1.75))
print("ratio (U_jet - c)/v_d(thermal, mean) = %.0f;  shift of the resonant layer v_d/|S| = %.1f rho_ref per thermal energy" % ((U.max() - c) / (w * vd_loc).sum(), (w * vd_loc).sum() / 1.75))

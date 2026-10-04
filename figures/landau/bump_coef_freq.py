"""Critical-layer response coefficient (projection of chi on g^yy - <g^yy>, source g^xx/<g^xx> - 1) vs frequency: does a finite trapping frequency raise the bump above 4.85?"""
import sys, numpy as np
sys.path.insert(0, "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/figures/landau"); sys.path.insert(0, "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/figures/lib")
import vortex_field as vf, bounce, response
geo = bounce.geometry(f"{vf.RB}/rb_c0p3_hxy/leg_0004/out/dipole_fix.dat"); wz = geo["J"] / geo["J"].sum()
Gxx = (wz * geo["gxx"]).sum(); Gyy = (wz * geo["gyy"]).sum(); sh = geo["gyy"] - Gyy; sh -= (wz * sh).sum()
R = response.Response(geo); R.source = lambda k2: k2 * (geo["gxx"] / Gxx - 1)
C0 = None
for wt in (0.002, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5):
    chi = R.solve(wt, 1.0); c = (wz * chi * sh).sum() / (wz * sh ** 2).sum()
    C0 = C0 or c
    print(f"wt {wt:5.3f}: projection {c.real:+.4f}{c.imag:+.4f}i   ratio to static {abs(c)/abs(C0):.3f}")

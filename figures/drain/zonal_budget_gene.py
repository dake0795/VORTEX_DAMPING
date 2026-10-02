import numpy as np, os, sys
D = "/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxoff/leg_0001/out"
NKX, NKY = 63, 16; REC = 2 + NKX*NKY; LB = 24
def rd(name, t0, t1, stride=2):
    out = 0; n = 0
    for sp in ("electrons", "positrons"):
        p = f"{D}/Spectral_2D_E_{name}_{sp}.dat"
        nrec = (os.path.getsize(p) - (NKX+NKY)*LB)//(REC*LB)
        acc = 0; cnt = 0
        with open(p, "rb") as f:
            head = np.array(f.read((NKX+NKY)*LB).split(), float)
            for i in range(0, nrec, stride):
                f.seek((NKX+NKY)*LB + i*REC*LB); t = float(f.read(LB))
                if t < t0 or t > t1: continue
                f.seek((NKX+NKY)*LB + i*REC*LB)
                a = np.array(f.read(REC*LB).split(), float)
                acc = acc + a[2:].reshape(NKX, NKY)[:, 0]; cnt += 1
        out = out + acc/cnt
    return head[:NKX], out, cnt
terms = ("fe", "nonlinear", "parallel", "curvature", "hyp_z", "hyp_v", "hyp_kperp", "drive", "hypzcomp", "sources", "coll")
for (t0, t1) in ((1000, 2000), (3000, 4400)):
    R = {}
    for nm in terms:
        kx, R[nm], cnt = rd(nm, t0, t1)
    L = np.sqrt(5000.0); m = np.abs(np.round(kx / (2*np.pi/3689.25))).astype(int)
    print(f"\nwindow {t0}-{t1} ({cnt} records/species); k_y = 0 electrostatic-energy budget by |m| (both signs summed)")
    print("  m  kxlamD    E_fe       nonlin      parallel    curvature   hyp_z      hyp_v     hypzcomp   sources")
    for mm in range(1, 32):
        s = m == mm
        print(f" {mm:2d}  {mm*2*np.pi/3689.25*L:5.2f}  " + " ".join(f"{R[nm][s].sum():+.3e}" for nm in ("fe", "nonlinear", "parallel", "curvature", "hyp_z", "hyp_v", "hypzcomp", "sources")))
    print("  totals: " + " ".join(f"{nm} {R[nm].sum():+.3e}" for nm in terms))

"""Extract the flute part of the potential from the GENE field files: Jacobian-weighted average along the field line
and the outboard mid-plane value, as (nt, nx, nky) complex arrays in GENE's (kx FFT order, ky = 0..15) layout.
  cache/gene_<run>.npz : t, phi_avg, phi_mid
"""
import sys
import numpy as np
from common import *

RUNS = {"leg1": f"{RB}/rb_c0p3_hxy/leg_0001/out/field.dat", "leg3": f"{RB}/rb_c0p3_hxy/leg_0003/out/field.dat", "leg4": f"{RB}/rb_c0p3_hxy/leg_0004/out/field.dat",
        "leg5": f"{RB}/rb_c0p3_hxy/leg_0005/out/field.dat",
        "a03": f"{RB}/rb_c0p3_hxy_ampltest/amp0p3/out/field.dat", "a003": f"{RB}/rb_c0p3_hxy_ampltest/amp0p03/out/field.dat"}
os.makedirs(CACHE, exist_ok=True)
w = metric()["w"]
for run in sys.argv[1:]:
    p = RUNS[run]
    n = nframes(p)
    sel = range(n)
    if run == "leg5":
        ts = times(p); sel = [i for i in range(n) if ts[i] < 11600.0]
    T, A, M = [], [], []
    for i in sel:
        t, phi = read(p, i)
        T.append(t); A.append(phi @ w); M.append(phi[:, :, 16].copy())
    np.savez(f"{CACHE}/gene_{run}.npz", t=np.array(T), phi_avg=np.array(A), phi_mid=np.array(M))
    print(run, len(T), T[0], T[-1], flush=True)

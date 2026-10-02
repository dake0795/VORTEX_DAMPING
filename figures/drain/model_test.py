"""Closed model of the decay: E_nz^{1/4} = A (E - E_c) [saturation omega_tr = C gamma_R with gamma_R linear in the
distance of the jets from marginal stability], dE/dt = -nu_L E_nz - s nu_h E_zon [Euler conserves E; Landau loss;
radial hyperdiffusion of the jets].  A, E_c from the hyp_x-off run (t < 4400) only; predictions for every other run."""
import numpy as np, sys, os
sys.path.insert(0, "../lib")
C = np.load("cache.npz"); nuL = float(C["nuL"]); NU_H = 3.4e-6
def bins(a, t0, t1, w=200.0):
    out = []
    for c in np.arange(t0 + w/2, min(t1, a[-1,0]) - w/2 + 1, w):
        m = (a[:,0] >= c-w/2) & (a[:,0] < c+w/2)
        if m.sum() > 10: out.append([c, np.exp(np.log(a[m,2]).mean()), a[m,1].mean()])
    return np.array(out)
h = bins(C["hxoff"], 1000, 1e9)
E = h[:,1] + h[:,2]
q = np.polyfit(E, h[:,1]**0.25, 1); A, Ec = q[0], -q[1]/q[0]
print(f"fit on nu_x = 0: A = {A:.4e}, E_c = {Ec:.6e}  (E_zon at t=1000: {h[0,2]:.6e})")
def run(t0, Enz0, Ezon0, s, t1, dt=1.0):
    # initial E from the slaving relation is not imposed: integrate E, slave E_nz, but start from measured E_nz via offset
    Etot = Enz0 + Ezon0
    off = Enz0**0.25 - A*(Etot - Ec)          # offset of this run from the universal line (reported)
    t, out = t0, []
    while t < t1:
        y = A*(Etot - Ec) + off
        Enz = max(y, 0)**4
        out.append((t, Enz, Etot - Enz))
        if y <= 0: break
        Etot += -(nuL*Enz + s*NU_H*(Etot - Enz))*dt; t += dt
    return np.array(out), off
S = {"ref": 1.0, "hxoff": 0.0, "nx128": 1/16, "hxy_d4": 0.25}
for r, t0 in (("ref", 1000.0), ("ref", 4400.0), ("ref", 9423.2), ("hxoff", 1000.0), ("nx128", 9500.0), ("hxy_d4", 9500.0)):
    a = C[r]; b = bins(a, t0 - 100, 1e9)
    m0 = np.argmin(abs(b[:,0] - t0 - 0)); 
    M, off = run(b[0,0], b[0,1], b[0,2], S[r], 13000)
    tc = M[-1,0] if M[-1,1] == 0 else np.nan
    print(f"\n{r} from t={b[0,0]:.0f}: offset {off/ b[0,1]**0.25:+.3f} of E_nz^1/4; predicted collapse (E_nz -> 0) at t = {tc:.0f}")
    for tt in np.arange(b[0,0], min(b[-1,0], 12500)+1, 600):
        i = np.argmin(abs(b[:,0]-tt)); j = np.argmin(abs(M[:,0]-tt))
        if abs(M[j,0]-tt) < 2: print(f"   t {tt:7.0f}: E_nz measured {b[i,1]:.3e} predicted {M[j,1]:.3e}  ratio {b[i,1]/max(M[j,1],1e-30):.2f}")

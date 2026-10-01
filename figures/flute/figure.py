#!/usr/bin/env python3
r"""Figure `flute`: how flute-like the vortex is, scale by scale (1 Oct 2026).

  python figure.py --rebuild   read field.dat of the reduced-box run (legs 3-5), write cache.npz
  python figure.py             draw ../flute.pdf from cache.npz

For the k_y = 1 potential phi(k_x, z) of the original run the script accumulates, over all output frames in two
windows (the slow phase of the decay and the approach to the collapse), the power spectrum sum_z |phi|^2 and the
non-flute part sum_z |phi - <phi>_z|^2 (uniform weights along the field line), per radial wavenumber.
(a) the spectrum against k_x lambda_D; (b) the non-flute fraction against k_x lambda_D.
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import notes_style as ns

NX, NKY, NZ = 64, 16, 32
N = NX * NKY * NZ
FR = 16 + 4 + N * 16 + 4
LX, LAMD = 3689.25, 70.71
WINDOWS = {"slow": (5000.0, 7000.0), "late": (10300.0, 10800.0)}
CACHE = os.path.join(HERE, "cache.npz")


def rebuild():
    acc = {k: [np.zeros(NX), np.zeros(NX), 0] for k in WINDOWS}
    for leg in (3, 4, 5):
        p = f"{ns.RB}/rb_c0p3_hxy/leg_000{leg}/out/field.dat"
        n = os.path.getsize(p) // FR
        with open(p, "rb") as f:
            for i in range(n):
                f.seek(i * FR + 4); t = np.frombuffer(f.read(8), "<f8")[0]
                w = [k for k, (a, b) in WINDOWS.items() if a <= t < b]
                if not w:
                    continue
                f.seek(i * FR + 20)
                phi = np.frombuffer(f.read(N * 16), "<c16").reshape((NX, NKY, NZ), order="F")[:, 1, :]
                tot = (np.abs(phi) ** 2).sum(axis=1)
                nf = (np.abs(phi - phi.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
                acc[w[0]][0] += tot; acc[w[0]][1] += nf; acc[w[0]][2] += 1
    kx = np.fft.fftfreq(NX, 1.0 / NX) * 2 * np.pi / LX
    np.savez(CACHE, kx=kx, **{f"{k}_tot": v[0] / v[2] for k, v in acc.items()},
             **{f"{k}_nf": v[1] / v[2] for k, v in acc.items()}, nframes=[v[2] for v in acc.values()])
    print("frames per window:", {k: v[2] for k, v in acc.items()})


def draw():
    C = np.load(CACHE)
    kx = C["kx"]
    kk = np.arange(1, NX // 2)                                   # |k_x| index 1..31
    x = kk * 2 * np.pi / LX * LAMD
    ns.apply_style()
    fig, (a, b) = ns.two_panel()
    for key, ls, lab in (("slow", "-", "slow decay"), ("late", "--", "before the collapse")):
        tot = C[f"{key}_tot"]; nf = C[f"{key}_nf"]
        T = np.array([tot[k] + tot[NX - k] for k in kk]); F = np.array([nf[k] + nf[NX - k] for k in kk])
        a.loglog(x, T / T.max(), ls, c="k", label=lab)
        b.loglog(x, F / T, ls, c="k")
    from matplotlib.ticker import FixedLocator, NullFormatter
    for ax in (a, b):
        ax.axvline(1.0, c="0.5", lw=0.8, ls=":")
        ax.axvspan(1.0, 4.2, color="0.93", lw=0, zorder=0)       # scales below the Debye length
        ax.set_xlim(0.105, 4.2)
        ax.set_xlabel(r"$k_x\lambda_D$")
        ax.xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1.0, 2.0]))
        ax.set_xticklabels([r"$0.2$", r"$0.5$", r"$1$", r"$2$"])
        ax.xaxis.set_minor_formatter(NullFormatter())
    a.set_ylabel(r"power of the vortex (normalised)")
    b.set_ylabel(r"non-flute fraction")
    b.set_ylim(top=1.5)
    a.legend(loc="lower left")
    ns.tag(a, "(a)", x=0.86); ns.tag(b, "(b)")
    ns.save(fig, "flute")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rebuild", action="store_true"); A = ap.parse_args()
    if A.rebuild or not os.path.exists(CACHE):
        rebuild()
    draw()

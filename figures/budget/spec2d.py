"""Reader for GENE's Spectral_2D_<term>.dat (diag_fespec2d.F90, formatted): header = (nx0-1) k_x values (k_x < 0 first,
Nyquist omitted) + nky0 k_y values; then per record: time, total, and the (nx0-1) x nky0 array (k_x outer, k_y inner).
The total is sum(k_y = 0) + 2 sum(k_y >= 1).  by_ky() returns t and the per-k_y contributions to that total."""
import numpy as np

NX, NKY = 64, 16
NKX = NX - 1
REC = 2 + NKX * NKY


def read(path):
    a = np.fromstring(open(path).read(), sep="\n")
    head = NKX + NKY
    n = (len(a) - head) // REC
    kx, ky = a[:NKX], a[NKX:head]
    r = a[head:head + n * REC].reshape(n, REC)
    return kx, ky, r[:, 0], r[:, 1], r[:, 2:].reshape(n, NKX, NKY)


def by_ky(path):
    kx, ky, t, tot, d = read(path)
    c = d.sum(axis=1)                      # (n, nky)
    c[:, 1:] *= 2.0
    assert np.allclose(c.sum(axis=1), tot, rtol=2e-3, atol=2e-3 * np.abs(c).sum(axis=1).max()), "total does not match"
    return t, c

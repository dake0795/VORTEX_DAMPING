# scripts/rayleigh: the fluid limit on the measured jets (1 Oct 2026)

Fluid limit of the notes (limit (b)): for a flute-like potential the charge moment of the transit-averaged equations
is 2D Euler with the field-line-averaged (Jacobian-weighted) metric, `zeta = (<g^xx> d_x^2 + <g^yy> d_y^2) psi`,
`psi = phi/C_xy`; `<g^xx> = 0.573`, `<g^yy> = 0.956`, so the Rayleigh problem has `k_eff = k_y sqrt(<g^yy>/<g^xx>)`.

| File | What it does |
|---|---|
| `common.py` | constants, GENE field reader, geometry (`out/dipole_fix.dat`), averaged metric |
| `extract.py RUN...` | flute potential (Jacobian-weighted field-line average, and mid-plane) of GENE runs -> `cache/gene_<run>.npz` |
| `measure.py` | fluid energies, two-frequency fit of the two vortices, critical layers -> `cache/measured_*.npz` |
| `linear.py` | Rayleigh operator on a frozen jet: dense eigenproblem at any number of radial modes (with or without GENE's hyperviscosity), linear initial-value problem |
| `growth_history.py` | growth rate and frequency of the unstable Rayleigh modes on the running-mean jet through the run -> `cache/growth_history.npz` |
| `euler.py NAME key=value...` | pseudo-spectral nonlinear 2D Euler model in the GENE box, started from GENE's potential -> `cache/euler_NAME.npz` |
| `batch1.sh` | the set of model runs behind `figures/euler` (run as one Slurm job, 4 cores, ~15 min) |
| `drift.py` | magnetic-drift frequency in GENE units from the geometry file and GENE's drift prefactor |

Checks: inviscid model conserves energy to 1e-10 and enstrophy to 2e-8 over 100 time units; eigenvalues converged
between 256 and 1024 radial modes; linear initial-value problem (1024 and 4096 points) reproduces the eigenvalue.
Sink in the model: GENE's hyperdiffusion acts on h (`hyp_on_h = T`), whose charge moment is `(1 + lambda_D^2 k^2) phi`,
so the vorticity is damped at `nu_k (1 + 1/(lambda_D^2 <k_perp^2>))`; at `k_x = 1`, `k_y = 0` this gives an energy decay
rate of 3.4e-6, the measured late-time decay of the zonal flow (3.45e-6).
The `cache/` directory holds the extracted potentials and model outputs (about 470 MB, not for git; regenerate with the scripts).

# VORTEX_DAMPING

Working notes by D. Kennedy (UKAEA) and G. G. Plunk (IPP Greifswald) on why the two breathing vortices of the
strongly driven (eta = 1) dipole pair-plasma condensate decay, and how fast. The notes take the earlier analytical
note (30 Sep 2026), correct and complete it in line with two independent audits, and compare it with a reduced-box
GENE simulation and two amplitude-reduced restarts.

## Build

```
cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING
env -u LD_LIBRARY_PATH pdflatex notes
env -u LD_LIBRARY_PATH bibtex notes
env -u LD_LIBRARY_PATH pdflatex notes
env -u LD_LIBRARY_PATH pdflatex notes
```

`env -u LD_LIBRARY_PATH` is needed on the CSD3 login nodes so that TeX finds its fonts. The build has no errors,
no undefined references and no overfull boxes. The output is 21 pages.

## Contents

| Path | What it is |
|---|---|
| `notes.tex` | The working notes. JPP class; the preamble and macros are copied from `ep_turbulence_paper/main.tex` (lines 1-160 and 190-196). The voice follows `STYLE_GUIDE.md`. |
| `notes.bib` | Bibliography. Entries from the Letter's `references.bib` and `letter-references.bib` are copied verbatim. The others (O'Neil 1965, Landau 1946, Fried & Conte 1961, Briggs et al. 1970, Rayleigh 1880, Drazin & Reid 2004, Haberman 1972, Stewartson 1978, Warn & Warn 1978, Ivanov et al. 2025 = arXiv:2405.00854, and the Letter as an unpublished working draft) are added. |
| `notes.pdf`, `notes.aux`, `notes.bbl`, `notes.blg`, `notes.log`, `notes.out` | Build products. |
| `jpp.cls`, `jpp.bst` | The JPP class and bibliography style, as in the long paper. |
| `STYLE_GUIDE.md` | The voice of Ivanov et al. (2025), taken from its LaTeX source. |
| `STATUS.md` | Status snapshot from 1 Oct 2026, 05:35, written before the notes were drafted. Now out of date. |
| `reference_style/ivanov2025_source/` | LaTeX source of Ivanov et al. (2025), on which the style guide is based. |
| `reference_style/vortex_damping_analytical_note.{tex,pdf}` | The earlier analytical note. |
| `source_note/vortex_damping_analytical_note.{tex,pdf}` | The same note, copied from `DIPOLE_TEST/VORTEX_CALCULATION/`. The notes start from it. |
| `figures/test1_amplitude_dependence.{pdf,png}` | Figure 2 of the notes: amplitude-reduced restarts and decay rate against amplitude. |
| `figures/test2_late_remnant.{pdf,png}` | Figure 4: the remnant after the collapse, with exponential and power-law fits. |
| `figures/test3_trapping_vs_damping.{pdf,png}` | Figure 3: jet and vortex profiles, trapping frequency against damping rate, and decay rate against cat's-eye width. |
| `figures/test4_entropy_energy.{pdf,png}` | Figure 1: E_nz, Z_nz, W_nz and Z_nz/E_nz. |
| `scripts/audit_linear/` | Audit of Secs 1-3 of the earlier note: `audit.py` and `k0.py` (K_0 identity, Im J, C(t)), `model.py` and `roots.py` (closed uniform-flow model, exact second-sheet roots), `extra.py`. |
| `scripts/audit_vortex/` | Audit of Secs 4-5 against the simulation data: `field.py` (field.dat reader), `energy.py` (windowed rates), `vstruct.py` (radial profiles), `phase.py` (phase speeds), `island.py` (resonant layer and trapped-region table), `series.py`, `fit.py` (rate against A_c), `hess.py` (Hessian trapping frequency). |

The scripts were copied from the session scratchpad that produced them. Some read intermediate caches (`orig.npy`,
`crit_*.npy`, `cache/*.npz`) that were not copied, because they are large and can be regenerated from the run
data. The scripts that drew the four figures are not in this repository.

## Data

The notes use only CSD3 data (read only):
- the reduced-box eta = 1 run `DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy`, legs 1 and 3-6 (t up to about 19 000);
- its two amplitude-reduced restarts in `DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_ampltest`, from the t = 5907 checkpoint, with every k_y >= 1 component scaled by
  epsilon = 0.3 and 0.03.

The run parameters are in Table 1 of the notes and in `DIPOLE_TEST/nl_reducedbox_20260925/README.md`.

## Main conclusions

The full statement is in Sec. 7 of the notes.
- The simulated vortices are in the drift-free (critical-layer) limit, not the uniform-flow drift-resonance limit.
  - Their critical layers lie on the jet flanks.
  - The magnetic drift shifts the resonance by about 1 rho_ref, against a trapped half-width of more than 100 rho_ref.
- The decay depends on amplitude, so it is not linear. It is also not an O'Neil crossover: the trapping frequency
  exceeds the damping rate by 10-1800 throughout.
- The data are consistent with a sink-limited nonlinear critical layer.
- The collapse coincides with the trapped half-width reaching one grid cell, about one Debye length.

## Proposed next tests (none run)

- Restart from t of about 9000 with nx0 = 128. If the collapse is set by the grid, it should move to an E_nz about
  16 times lower.
- Compute the Rayleigh quasi-mode on the measured jet profile.
- Measure the mode-resolved k_y = 1 -> 0 transfer to see whether the vortex energy is dissipated or absorbed by the jet.
- Scan the perpendicular hyperdiffusion at fixed resolution.

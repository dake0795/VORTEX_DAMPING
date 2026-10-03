# VORTEX_DAMPING

Working notes by D. Kennedy (UKAEA) and G. G. Plunk (IPP Greifswald) on why the two breathing vortices of the
strongly driven (eta = 1) dipole pair-plasma condensate decay. `notes.pdf` is the compiled document.

## Where the notes stand (3 Oct 2026)

Full record, with every retraction: `STATUS.md`, items (a)-(j). The notes open with a list of explanations ruled out.

- **The vortices are the shear-flow instability of their own jets** (Rayleigh's equation on the measured jets gives
  their frequency and structure); the instability saturates by trapping.
- **Their loss is linear and kinetic**: Landau damping, along the field line, of their small non-flute part, at a
  calculated rate with no adjustable parameter that matches the measured one in every run.
- **The slow drain is physical** (the same without radial hyperdiffusion). **The jets pay** for most of the loss and
  their energy budget closes; the exchange is carried by the momentum the damped vortices deposit in the jets and by
  the stress of the vortices computed locally along the field line, where their non-flute part matters.
- **The late acceleration and the collapse** of the reference run are numerical (radial hyperdiffusion and grid).
- **Open:** the size of the vortices' non-flute part at their critical layers (the kinetic response gives it to within
  a factor of two), hence the drain rate; a fluid model local on both sides; the end state. Drain tests in
  `DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_draintest` are running.

## Build

```
env -u LD_LIBRARY_PATH pdflatex notes && env -u LD_LIBRARY_PATH bibtex notes
env -u LD_LIBRARY_PATH pdflatex notes && env -u LD_LIBRARY_PATH pdflatex notes
```

(`env -u LD_LIBRARY_PATH` is needed on the CSD3 login nodes so that TeX finds its fonts.)

## Contents

| Path | What it is |
|---|---|
| `notes.tex`, `notes.bib`, `notes.pdf` | The working notes (JPP class). Voice and house rules: `STYLE_GUIDE.md`. |
| `STYLE_GUIDE.md` | The voice of Ivanov et al. (2025, arXiv:2405.00854), from its source in `reference_style/`; section 7a holds the rules of these notes: figures in the long dipole paper's style (two panels per row, all text at body size), and no simulation numbers unless the argument needs them. |
| `figures/lib/notes_style.py` | Figure style (wraps `ep_turbulence_paper/scripts/lib/paper_style.py`). `traces.py`, `vortex_field.py`: shared readers. |
| `figures/<name>/figure.py` | One directory per figure: `--rebuild` reads the run data and writes `cache.npz`; the default draws `figures/<name>.pdf` from the cache alone. |
| `figures/cartoon.tex` | The TikZ cartoon (figure 1). |
| `figures/measurements/` | Measurements used in the text without a figure of their own: modulation of the jets by the breathing, mode-resolved budgets, the loss rate in the production run. |
| `figures/budget/` | Free-energy (W) budget of the vortices; cross-check of `figures/exchange` (not shown in the notes); also draws `dissipation`. `spec2d.py` reads GENE's `Spectral_2D_*.dat`. |
| `scripts/rayleigh/` | Rayleigh problem and two-dimensional Euler model on the measured jets (own README). |
| `scripts/audit_vortex/` | First field analysis (superseded by `figures/lib/vortex_field.py`). |
| `source_note/` | The earlier analytical note (30 Sep 2026) that these notes start from. |
| `jpp.cls`, `jpp.bst` | JPP class and bibliography style, as in the long paper. |

## Figures of the notes

| Figure | File | What to learn from it |
|---|---|---|
| 1 | `figures/cartoon.tex` | The mechanism: the vortices live on the instability of their jets and die of Landau damping when it fails. |
| 2 | `history` | The vortices decay ever faster and collapse; the zonal flow barely changes. |
| 3 | `amplitude` | A vortex made smaller on unchanged jets regrows, then decays faster. |
| 4 | `geometry` | Critical layers on the jet flanks, with a cat's eye at each. |
| 5 | `rayleigh` | The vortices are the unstable shear mode of their corrugated jets. |
| 6 | `corrugation` | The vortices maintain the corrugation; the jets are only just unstable. |
| 7 | `saturation` | Trapping frequency and Rayleigh growth rate decay together: omega_tr = C gamma_R. |
| 8 | `exchange` | The loss to entropy is steady and linear; the supply from the jets fades and reverses. |
| 9 | `channels` | The loss is carried by the streaming of particles along the field line. |
| 10 | `landau` | The loss is Landau damping of the non-flute part, at the calculated rate, with no adjustable parameter, in every run. |
| 11 | `dissipation` | The perpendicular sink acts on the zonal flow, the velocity-space sink on the vortices. |
| 12 | `euler` | A fluid with the same jets has the same vortices but not their kinetic loss. |
| 13 | `sinktests` | The decay follows the perpendicular sink, not the velocity-space one; with a weak perpendicular sink the vortices persist. |
| 15 | `supply` | The jets pay for the part of the loss that the vortices do not drain (energy budget of the jets, no adjustable constant). |
| 14 | `trapping` | The cat's eyes turn over many times during the decay; the collapse comes as they narrow to a grid cell. |

Appendices (`appendices/*.tex`): step-by-step derivations at undergraduate level - A: Euler limit and Rayleigh's
equation; B: trapping and saturation; C: Landau damping of the non-flute part (the loss formula).
`figures/flute` (non-flute power of the vortices by radial scale) is kept as a measurement but is not in the notes.

## Data

CSD3 only (read only): the reduced-box eta = 1 run `DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy` (legs 1, 3-6), its
amplitude-reduced restarts `rb_c0p3_hxy_ampltest/{amp0p3,amp0p03}`, the restarts `rb_c0p3_hxy_sinktest/*` (running),
and, for one comparison, the production run `DIPOLE_TEST/eta1_production_NEW_sep26`. Still to be read in when they
have run: `rb_c0p3_hxoff` and `rb_c0p3_debye_{x4,d4}` (see `STATUS.md`).

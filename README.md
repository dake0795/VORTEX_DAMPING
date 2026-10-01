# VORTEX_DAMPING

Working notes by D. Kennedy (UKAEA) and G. G. Plunk (IPP Greifswald) on why the two breathing vortices of the
strongly driven (eta = 1) dipole pair-plasma condensate decay. `notes.pdf` is the compiled document.

## Where the notes stand (1 Oct 2026)

The notes give the theory of the vortices and test it; they present only the calculation that describes the
simulations (house rule, `STYLE_GUIDE.md` section 7a).

- **The vortices are the saturated shear-flow instability of their own jets.** In the fluid limit the jets and
  vortices obey the two-dimensional Euler equation; jets corrugated at the Debye scale are Rayleigh-unstable; the
  instability saturates by trapping. Rayleigh's equation on the measured jets gives the vortex frequency and
  structure, and omega_tr = alpha gamma_R holds through the decay.
- **Their damping is linear and kinetic**: Landau damping, along the field line, of the part of the vortex potential
  that is not flute-like, at a rate independent of amplitude.
- **They decay because the jets slowly lose their instability**, and collapse when the jets become stable.
- **Open:** the kinetic loss from first principles (in progress); what erodes the jets' corrugation; grid scale
  against Debye scale. Five restarts that bear on these are running
  (`DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_sinktest`).

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
| 1 | `figures/cartoon.tex` | The mechanism: jets feed the vortices; phase mixing along the field line drains them. |
| 2 | `history` | The vortices decay ever faster and collapse; the zonal flow barely changes. |
| 3 | `amplitude` | A vortex made smaller on unchanged jets regrows, then decays faster. |
| 4 | `geometry` | Critical layers on the jet flanks, with a cat's eye at each. |
| 5 | `rayleigh` | The vortices are the unstable shear mode of their corrugated jets. |
| 6 | `saturation` | The vortex amplitude follows the jets' growth rate: omega_tr = alpha gamma_R. |
| 7 | `exchange` | The loss to entropy is steady and linear; the supply from the jets fades and reverses. |
| 8 | `channels` | The loss is carried by the streaming of particles along the field line. |
| 9 | `flute` | The vortices are flute-like at the box scale, not at the Debye scale. |
| 10 | `dissipation` | The perpendicular sink belongs to the zonal flow, the velocity-space sink to the vortices. |
| 11 | `euler` | A fluid with the same jets has the same vortices but not their kinetic loss. |
| 12 | `trapping` | The cat's eyes turn over many times during the decay; the collapse comes as they narrow to a grid cell. |

Appendices (`appendices/*.tex`): step-by-step derivations at undergraduate level (Euler limit and Rayleigh's equation;
trapping and saturation; Landau damping of the non-flute part - in preparation).

## Data

CSD3 only (read only): the reduced-box eta = 1 run `DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy` (legs 1, 3-6), its
amplitude-reduced restarts `rb_c0p3_hxy_ampltest/{amp0p3,amp0p03}`, the restarts `rb_c0p3_hxy_sinktest/*` (running),
and, for one comparison, the production run `DIPOLE_TEST/eta1_production_NEW_sep26`.

# Status (1 Oct 2026, 05:35)

## In this repository
- `STYLE_GUIDE.md`: the voice of Ivanov et al. (2025), arXiv:2405.00854, taken from its LaTeX source (`reference_style/ivanov2025_source/`).
- `jpp.cls`, `jpp.bst`: the JPP class, as in the long paper.
- `source_note/`: Dan's analytical note (30 Sep), which these working notes start from.

## Coming next
`notes.tex`, the working notes, is being drafted. Before drafting, the note is audited and tested against the simulations.

## Audit so far: Sections 1–3 of the note (linear response, Landau resonance, energy integral)
They are correct. Every equation was re-derived and checked numerically: K_0, Im J, and |C|² = (1+a²t²)^(−3/2).

There are two substantive gaps:
1. **Zonal Poisson fixes part of the background.** It forces H_a = e_a f_0 Φ/T + δF_z,a, so the drive is
   N_a = (e_a f_0/T)(ω − kU + k v*_a) − k δF′_z,a.
   This is the Doppler-invariant form, and b_0 must be built from ω − kU.
2. **The weak-damping formula describes a quasi-mode only.** In the pair plasma the drift resonances of the two species cover the whole real frequency axis. So the weak-damping formula describes a second-sheet quasi-mode, valid only when |γ| ≪ |ω − kU|.

## Data tests
The amplitude tests support amplitude-dependent damping. A vortex scaled to 0.3 decays 2× faster than the original, and one scaled to 0.03 about 10× faster.

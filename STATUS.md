# Status (2 Oct 2026, 23:40) - overnight work, see the block at the top

## 2 Oct 2026, night: what the slow drain is (work in progress; notes.tex not yet rewritten)

1. The slow drain is physical: `rb_c0p3_hxoff` (hyp_x = 0) drains like the reference to t 4,400; only the late
   acceleration and the collapse need the radial hyperdiffusion (and the coarse grid).
2. The supply budget closes in four runs with no free constant (`figures/supply.pdf`).
3. 2D Euler + the predicted Landau damping applied as a friction on the vortices does NOT reproduce GENE: the vortices
   decay at nu_L and the jets keep their energy (`scripts/rayleigh/batch_drain.sh`, cache/euler_dr_*.npz).
4. Why: in GENE the jets are forced by more than the Euler Reynolds stress. Breathing-averaged, the Euler stress of the
   vortices would ACCELERATE the jets; the jets decelerate (`scripts/rayleigh/zonal_forcing2.py`).
5. The missing force is the momentum the resonant particles absorb from the Landau-damped vortices, f(x) = k P(x)/wt(x),
   deposited in the zonal charge by the radial current it drives (`figures/landau/deposit.py`; P from the verified loss
   formula, no free constant). It accounts for the jets' non-Euler energy loss (8.8 vs 9.4, 3.65 vs 3.59, 2.0 vs 2.5 in
   three windows of 1000 t.u.) and for the forcing of the dominant jet harmonic m = 1 (amplitude 1.1, 1.1, 0.7; phase
   -3, -10, -32 deg). Higher harmonics (m >= 5) carry a large Euler stress cancelled by something not yet identified.
6. Where wt and w have opposite signs (resonance beyond the critical layer), absorption in the local frame is a GAIN of
   vortex energy in the lab frame: the damping moves energy from the jets to the vortices and to heat.
   Same test on the reference (t 6000-9400, its zonal hyperdiffusion added to the Euler forcing): phase again locked
   (0 to -4 deg on m = 1, 3), energy 72 %, 55 %, 34 % of the non-Euler loss as the collapse approaches; the remainder
   grows where the hyperdiffusion and the Debye-scale grid act on the fine structure of the jets.
7. Closed fluid model (`scripts/rayleigh/euler.py kin=1`): 2D Euler + the kinetic damping in each vortex's local frame
   (vortices split with the linear Rayleigh tendency) + the deposited momentum. Its kinetic terms reproduce GENE's
   ratios (lab-frame vortex gain / absorbed power -0.73 vs -0.72; jets' loss / absorbed power -1.73 vs -1.72) but the
   model vortices GROW (6.5e3 -> 2e4) while GENE's decay: the Euler stress acts on the jets' radial harmonics at
   k_x lambda_D >~ 0.5, where in GENE it is cancelled (measured: the observed tendency there is ~0, residual = -Euler).
   Diagnostic runs with the jets' Euler response removed above k_x lambda_D = 0.5 / 1.0 bracket GENE (vortices die too
   fast / stay steady). So the drain is decided by the zonal response at the Debye scale, where the flute limit
   fails; the missing piece is the kinetic zonal response (and kinetic stress) at k_x lambda_D ~ 0.5-1.
   Later the same evening: GENE's own k_y = 0 budget (`figures/drain/zonal_budget_gene.py`) shows that at m >= 5
   (k_x lambda_D >= 0.6) the exact nonlinear input to the jets' electrostatic energy is cancelled by the
   parallel-streaming term (m = 5: -737 / +747; m = 7: -686 / +714), up to m = 15; at m = 3 streaming is the main loss.
   Tested and REJECTED: Landau damping of jet harmonics oscillating at the breathing frequency (`figures/landau/
   zonal_landau.py`): wrong sign at m >= 5, 14-100x too small at m = 1, 3. The per-harmonic comparison of the flute-Euler
   stress with GENE's exact transfer is too noisy above m = 3 to support any reduction factor (earlier 'Euler cancelled at
   m >= 5' from the flute analysis is therefore not established; GENE's internal nonlinear/streaming cancellation is).
   2 Oct, ~22:00-23:45 (`figures/landau/local_stress.py`, `compare_budget.py`, `nonflute_check.py`, `nonflute_shape.py`):
   (a) The stress of the vortices on the jets must be computed LOCALLY along the field line (local potential, local
       metric) and then averaged: the stress of the averaged potential is wrong on the jets' (odd) harmonics, by sign at
       m = 3 and by factors 3-40 at m >= 5; even harmonics agree. The local metric alone changes nothing: the whole
       difference comes from the ~2 % non-flute part of the vortices.
   (b) THEORY = local stress + momentum deposited by the Landau damping reproduces GENE's exact jet budget
       (E_nonlinear + E_parallel at k_y = 0) on the dominant harmonics: m = 1: -2.9e-5 vs -2.0e-5 (t 1000-2000),
       -1.08e-5 vs -1.21e-5 (t 3400-4400); m = 3: -2.30e-4 vs -2.71e-4, -1.45e-4 vs -1.49e-4. (Rates per unit energy.)
   (c) The measured non-flute part of each vortex is a single shape along the line (98.9 %, ~ g^yy - <g^yy>), real and
       in antiphase with phi_0, largest at the critical layers (|phi_1/phi_0| ~ 0.28 per unit shape). The linear response
       of appH (chi_s, used by the loss formula) predicts +0.06 there (wrong sign, 4.5x too small); a static local
       metric source gives the right sign at the critical layers but half the size and fails inside the jets.
       So the reactive non-flute structure of the vortices near their (trapped, nonlinear) critical layers is not
       given by the linear theory, although its dissipative part (the loss formula) is.
   (d) Models: closed fluid model with deposition + flute stress fails (vortices grow); freezing the jets' response
       above k_x lambda_D ~ 0.5 also fails (jets gain energy): the local stress is needed at all harmonics.
   OPEN, the next calculation: the non-flute structure of a vortex at a trapped critical layer (static limit inside the
   cat's eye, with the metric source acting on the eye's fine radial structure), to put the local stress into the model.
8. Next: derive that response (zonal mode with its non-flute part s(l) k_x^2 lambda_D^2 phi_0; stress with the
   orbit-averaged potentials), put it in the model, then test against rb_c0p3_hxoff and the drain tests.
9. (superseded) 7. Next: a fluid model with this kinetic damping in each vortex's local frame and the momentum deposited in the jets
   (closed, no free constant), against `rb_c0p3_hxoff` and the drain tests in `rb_c0p3_draintest`.
10. Note: `figures/landau/predict.py`'s bounce.G estimate is 10x the verified rate; the verified formula is the one in
   `figures/landau/figure.py` (response.Response, D(wt)).


Take-over document for the whole project: `PROJECT_DIPOLE/audit/reports/HANDOVER_20261001.md` (rules, every job, what
to do when each test reports). This file is the state of the notes.

## Notes

- `notes.pdf` (31 pp) is current with `notes.tex` at the head of `main`. Structure: introduction with the cartoon ->
  the setting -> theory of the vortices (Euler limit, Rayleigh instability of the corrugated jets, trapping saturation
  omega_tr = C gamma_R, the loss formula, the decay and its predictions) -> comparison with the simulations -> summary;
  appendices A (Rayleigh), B (trapping), C (Landau damping of the non-flute part) at undergraduate level.
- House rules (`STYLE_GUIDE.md` Sec. 7a): only the calculation that describes the simulations, written as though the
  answer was always known; no simulation numbers unless essential (scalings in figures); at most two panels per row,
  body-size text, one message per panel. Discarded material is in the git history up to commit 35058e6.
- The loss formula is parameter-free and matches the measured parallel-channel loss to about 10 % in the reduced box,
  the amplitude-reduced restarts and the production box (`figures/landau`).

## Tests running (DIPOLE_TEST/nl_reducedbox_20260925; all restart from t = 9423 unless said)

Vortex energy relative to the restart, 100-time-unit log-means, read at 18:55 on 1 Oct:

| run | change | t reached | at t = 10,000 | at t = 10,250 |
|---|---|---|---|---|
| reference `rb_c0p3_hxy` | - | 21,710 | 0.41 | 0.22 (collapse at ~11,000) |
| `hv_d4` | velocity-space sink / 4 | 10,071 | 0.42 | - |
| `hxy_d4` (+ `hxy_d4_leg2` queued) | perpendicular sink / 4 | 10,079 | 0.68 | - |
| `hxy_x4` (ended) | perpendicular sink x 4 | 9,850 | collapsed (0.18 at 9,700) | - |
| `nx128` (+ `nx128_leg2` queued) | nx 128, same coefficient: radial damping 16x weaker at fixed k_x | 10,333 | 0.88 | 0.83 |
| `nx128_hx16` (+ `nx128_hx16_leg2` queued) | control: nx 128 with the reference's radial damping | 9,536 | too early | - |
| `rb_c0p3_hxoff` (from noise, 3-day legs, self-chaining) | hyp_x = 0 throughout | 1,328 | tracks the reference so far | decisive at t > 11,000, about Sun 4 Oct |
| `rb_c0p3_debye_x4`, `rb_c0p3_debye_d4` (from noise) | Debye length x 2 and / 2 | queued | start Fri 2 Oct ~19:00 and Sat 3 Oct ~01:25 | test of the k^2 lambda_D^2 scaling of the loss formula |

Conclusion so far (in the notes as the subsection "What sets the rate of the decay" and figure `sinktests`): the decay
follows the radial hyperdiffusion and is indifferent to the velocity-space sink, so the vortex lifetime of the
reference run is numerical. Not yet established: that the control decays like the reference; that the vortices of
`rb_c0p3_hxoff` outlive t ~ 11,000.

## To do as the tests report

1. `figures/sinktests/figure.py --rebuild` (last rebuilt with records to t ~ 10,000), adding the control and the
   continuation legs to its run list.
2. Notes: state the control's result; add the `rb_c0p3_hxoff` outcome; add the Debye scan as a test of the loss
   formula (`figures/landau/predict.py` is the template: predicted against measured loss rate).
3. Recompile, update this file and the "Where the notes stand" block of `README.md`, commit with `notes.pdf`, push.

## Open physics

The constant C of the saturation law; the corrugation at which the jets become stable; what erodes the corrugation
in a real (collisional) plasma; the long-time state when the vortices persist (expected: the Landau loss slowly
drains the jets). A field-based zonal energy came out about 3.3 times GENE's `E_fe(k_y = 0)` in a side calculation
and is unexplained; nothing in the notes uses it.

## Related, elsewhere

- `figures/measurements/hydro_scalings.py` (two-dimensional-hydrodynamics scalings in a run) was the first pass of what
  is now Appendix A.4 of the long paper (`ep_turbulence_paper`, branch `v3.1-rollback`, `figures/fig43_hydro_spectra`
  and `fig44_hydro_rates`), which uses the long paper's own runs.

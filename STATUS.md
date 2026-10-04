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
   (e) 3 Oct, night (`figures/landau/nonflute_amp.py`, `nonflute_linear.py`, `nonflute_jets.py`, `stress_predicted.py`):
       the vortices' non-flute structure is LINEAR (identical in the eps = 0.03 and 0.3 restarts and the reference:
       -0.39 to -0.44 at the critical layers, -0.137+0.058i mid-jet). In the linear phase from noise (no jets) it is the
       appH response in magnitude and k_y^2 scaling (-0.067 vs 0.063 at k_y = 1; the sign is a Fourier-convention
       flip, fixed consistently: response conjugated). With jets, a NEW linear term is needed: (ii) the jets' own
       non-flute polarisation charge Q_z1 = -lambda_D^2 [g^xx d_x^2 phi_z1]_nf advected radially by the vortex,
       phi_1 - rho[phi_1] = -(k/wt) psi_0 d_x Q_z1 (sign per the convention), plus (iii) the vortex's charge advected by
       the jets' non-flute flow U_1 (small). Terms (i)+(ii)+(iii) reproduce the measured structure with best scale
       0.84-1.09 for both vortices in both runs (no free constant; residual ~0.5, from a missing real part mid-jet).
       At the critical layer itself (ii) is singular (1/wt): the kinetic analogue of the Rayleigh critical-layer
       singularity (the vortex advecting the jets' pitch-angle-dependent distribution); on GENE's 64-point grid it
       spoils the predicted stress at m >= 3. Regularisation to implement next: flattening of the jets' distribution,
       every pitch-angle class, inside the cat's eye (half-width w = 2 (|psi_1|/|U'|)^(1/2), measured, no free constant).
   (f) 3 Oct, later in the night (`figures/landau/nonflute_collapse.py`, `bump_vs_vorticity.py`): CORRECTION to (e):
       term (ii) is antisymmetric across each critical layer, while GENE's non-flute part is a symmetric bump there; the
       earlier agreement came from points on one side. What the data show instead: the non-flute ratio phi_1/phi_0
       (projection on g^yy - <g^yy>) is ONE function of the Doppler-shifted frequency for both vortices, both runs and
       all amplitudes: universal away from the critical layers (-0.17 at wt -1.4, -0.07 at -1.0, +0.06 at +1.4) and a
       bump at the critical layer, about +-0.4 wide in wt (one Debye length in x), deepening through the decay
       (-0.25, -0.36, -0.43, -0.50, -0.54 at t 2000-9800). The bump is proportional to the vortex's own vorticity at its
       critical layer: bump / (lambda_D^2 zeta_1/phi_0) = 6.0-8.1 in every run, stage and amplitude.
       Physics: at the critical layer the fluid co-moves with the vortex, so the vortex's vorticity is a charge spread
       along the field line by the local metric; the static kinetic response to (g^xx/<g^xx> - 1) gives coefficient 4.85
       (prediction -0.23 vs measured -0.43 at t 6100): right sign and form, a factor ~2 low (in the linear phase,
       without critical layers, the same response matches to 10 %).
   (g) 3 Oct, ~00:00 (`figures/landau/stress_rank1.py`, `scripts/rayleigh/euler.py lst=C`): the non-flute part in rank-one
       form A(x) sh(z) reproduces the stress of GENE's full non-flute part at every harmonic (m = 1-9), so the problem
       reduces to A(x). A predicted (C lam^2 zeta_1 at the critical layer, C = 4.85 from the static response, appH
       response elsewhere) gives the stress on the jets' dominant harmonics m = 1, 3 to 0-40 % with the right phase
       (flute: wrong phase); m >= 5 still flute-like. Closed fluid model with this local stress on the jets only: the
       jets lose 13 % (GENE 4 %), the vortices grow (C = 4.85) or the run goes unstable (C = 7): the jet-vortex
       interaction must be local along the field line on BOTH sides (the vortex's exchange term too), which needs the
       jets' non-flute structure as well.
   (h) 3 Oct, ~01:30 (`figures/landau/jets_nonflute.py`, `jets_nonflute_prod.py`): the JETS' non-flute part is the
       static kinetic response to (g^xx(l) - <g^xx>) d_x^2 phi_z0 in shape (0.98-0.999 for m = 1-5) and sign, with
       magnitude 1.5-1.7x the prediction in the reduced box (every time, both runs) but 0.96-1.38x in the newer
       production run (same hyp_z 0.05; nz0 32): not geometry (g^xy = 0; GENE's k_perp^2 is the file's metric) and
       not hyp_z. Reading: in static or co-moving structures (jets; the vortices' critical layers) a collisionless plasma
       does not relax to the equilibrium response, so the non-flute structure keeps memory of how the structure formed;
       the equilibrium response gives the shape and 0.6-1x the size. The vortices' critical-layer coefficient (6-8
       measured vs 4.85 predicted) is the same effect.
   (i) 3 Oct, ~02:00, correction to (h): the reduced-box factor (1.5-1.65) is constant from t 800 (just after the burst)
       to t 8000 while the vortex energy falls ~100x, and both the reduced box and eta1_production_NEW_sep26 (factor
       0.96-1.38) are collisionless with the same hzmask binary, hyp_z 0.05 and identical geometry. So the 'collisionless
       memory' reading of (h) is not supported by anything specific, and the box dependence of the factor (1.63 vs 1.22
       at the same k_x lambda_D = 0.12) is UNEXPLAINED. Differences left between the two: box size (lx 3689 vs 14757),
       k_y,min (2e-3 vs 5e-4), nky0, stage of the decay. Treat the equilibrium static response as giving shape and sign,
       and the size to within a factor 1-1.7.
   (j) Drain tests at 3 Oct 04:20 (rb_c0p3_draintest, restarts from t 2500, hyp_x = 0):
       - rb_c0p3_hxoff vs reference: E_nz ratio 1.10-1.11 at t 3000-4000, then 1.39 (t 5000) and 1.25 (5500): the radial
         hyperdiffusion matters from t ~ 4500, not only at the collapse (correction to item 1 of the 2 Oct block).
       - vort0p5_hxoff (vortices halved, jets unchanged): NO regrowth over 500 t.u. (4.1e6 -> 3.8e6 vs 1.1-1.3e7 unscaled):
         the vortex amplitude is NOT slaved to the jets (against the drain-law reading and the early 'regrowing' remark).
       - all0p5_hxoff (everything halved): drains 5x slower (Gamma_E 1.0e-4 vs 5.3e-4); Euler scaling alone gives 2x.
   (k) 4 Oct (`figures/landau/debye_scan.py`): the LOSS FORMULA holds across the Debye scan with no free constant:
       debye_x4 (lambda_D^2 x4, vortex frequency 0.43) measured/predicted 0.80-0.99; reference 1.01-1.08; debye_d4
       (lambda_D^2 /4, frequency 3.4) 1.01-1.14. The measured rate changes only ~3x over 16x in lambda_D^2 because the
       frequency (8x) enters D(wt); the formula captures both.
   (l) 4 Oct: rb_c0p3_hxoff has NOT collapsed by t 10,216 (E_nz 8.4x the reference at t 10,000); its decay rate is
       roughly constant, 2-4e-4 (about 0.15 nu_L), from t ~4500: near-exponential drain; the E^(7/4) drain law fitted
       to t < 4400 fails after. hxyoff (both perpendicular sinks off) = hxoff (ratio 1.01 at t 3000). nx128_hxoff drains
       ~2x faster than hxoff at t 3500-4400 while its jets lose ~25 % less: the split of the Landau loss between jets
       and vortices is NOT converged in radial resolution (Debye scale ~1.2 dx at nx 64). nx256_hxoff queued
       (37266546, + leg 2 37266547; starts ~5 Oct 15:00 on cpu-r-hbm-[2,4,5,35]). Control nx128_hx16: 1.23x the
       reference at t 10,000, 1.98x at 10,500 (decays like the reference, collapse delayed).
   (m) 4 Oct, evening (`figures/lib/fieldgen.py`; `figures/landau/bump_resolution.py`, `eye_charge.py`, `energy_norm.py`,
       `jet_budget_clean.py`, `vortex_drain_pred.py`) - CORRECTIONS and new results:
       - The critical-layer bump/vorticity ratio is resolution-independent (6.0-6.2 at nx 64 and 128, t 2800-4300;
         8.0 both at t 9600-10000): the factor vs the static theory is physics, not numerics.
       - At the critical layers the vortex CHARGE is nearly flute (its profile along the line has 0.2-0.4 of the
         metric variation): parallel streaming equilibrates charge along the line; the potential carries the non-flute
         part. Flat charge + static response gives ~80 % (early) to ~60 % (late) of the measured bump.
       - GENE's E_fe(k_y = 0) = lambda^2 int J g^xx k_x^2 |phi|^2 / int J to 0.05 %: = 2 lambda^2 C^2 x the field-based fluid
         energy (~8540, falling ~12 % to fine harmonics). The old '3.3x' side issue is not real.
       - RETRACTION of (b), (g) and of the notes' 'local stress': (b)'s per-harmonic agreement had a factor-2 error
         (energy of both +-m without the 1/2). With one normalisation, GENE's nonlinear forcing of the jets at m = 1 is
         flute Reynolds stress + deposition (-1.87 vs -1.65 at t 3000-4000; -3.15 vs -2.05 at 1500-2500); the 'local'
         stress (local phi and local metric) overshoots 2-3.5x and is the wrong object (particles are advected by their
         orbit-averaged potential; the advected charge includes the non-adiabatic response). At m >= 5 GENE's
         nonlinear forcing IS the flute stress (m 5: -0.032 vs -0.041; m 7: -0.090 vs -0.073), returned by parallel
         streaming. m = 3 unexplained.
       - The vortex drain = small difference of two large exchanges with the jets: flute stress (vortex -> jets) and the
         lab-frame effect of the damping (jets -> vortex through the deposited momentum), each ~ the absorbed power,
         nearly cancelling; measured drain 15-35 % of the absorbed power; with each term uncertain by 10-50 % the
         predicted sum has the wrong sign in 3 of 6 windows. Predicting the drain needs both to a few per cent; it also
         explains the resolution sensitivity of the drain.
   (n) 4 Oct, late afternoon: Debye scan added to figure `landau`; new appendix `app:momentum` (appI_momentum.tex:
       energy/momentum ratio omega/k of a travelling wave, radial displacement dx = dp_y/(q B), current f/B, force on
       the jets f = kP/wt; the damping conserves the momentum of jets + vortices); new figure `balance` (the two
       exchanges and the measured drain, nx 64 and 128) with the momentum-budget explanation of the near-cancellation.
       At nx 128 the Reynolds stress alone matches the measured drain; the estimate of the damping exchange is too large
       by ~0.5 P there. Collisional restart queued (rb_c0p3_draintest/coll_hxoff 37267574 + leg 2 37267575, Landau
       1e-2, hyp_x = 0, cpu-r-hbm-[78-79] after debye_x4 leg 2). Literature search delegated (LITERATURE.md, pending).
       Not done: pitch-angle-resolved (g1.dat) stress diagnostic (judged non-predictive); the missing static factor;
       the m = 3 harmonic.
   (o) 4 Oct, evening: GENE's exact budget by k_y (`figures/landau/ky_budget.py`, figure `partition`): per unit of the
       vortices' Landau loss P1, the jets supply ~1.1-1.4 P1 nonlinearly; the vortices receive ~0.9-1.1, their harmonics
       (k_y >= 2) ~0.2-0.25, which are Landau damped as fast; drift loss 0.075; the vortices drain by 0.26 (nx 64) and 0.31
       (nx 128) on average. NL sums to zero over k_y in every window. The missing channel in (m)/(n) was the harmonics.
       Loss formula on the harmonics (`harmonics_loss.py`, extended D table `cache_Dext.npz`): k_y = 3 about right, k_y = 2
       at 55-65 %; total 60-85 % of GENE's harmonic loss, no free constant. Measured local D instead of theory D changes
       the damping exchange by only 2-6 % (`drain_measuredD.py`). Literature (`LITERATURE.md`, sub-agent, key entries
       re-verified on Crossref): closest prior art Crooks & O'Neil 1995 (rotational pumping of the diocotron mode);
       related-work paragraph added to the introduction. NOTE: github.com/dake0795/VORTEX_DAMPING is PUBLIC (checked via
       the GitHub API, 4 Oct) - for Dan.
   (p) 4 Oct, late: the jets' missing static factor (1.5-1.65) is NOT numerical in the kinetic solver (static response
       converged in pitch-angle grid x4, bounce harmonics x2, npsi x2: 1.4588 unchanged) and NOT the GENE time step
       (rb_c0p1_hxy, Courant 0.1: m = 1 factor 1.60/1.63 vs 1.58/1.61 at Courant 0.3, t 1200/1800). Decisive test queued:
       zonal_static/m1 (linear, nx0 3, ky0 only, reduced-box grid; job 37268642, 1 core pinned cpu-r-hbm-37, waits on
       Dan's cs450p ~6 h) - GENE's own linear static non-flute part of one zonal harmonic vs theory.
   (q) 4 Oct, night: linear GENE test of the jets' factor (zonal_static/m1_chpt, job 37270554, 18 min on 1 core): the
       m = 1 jet (k_y = 0, k_x = +-1) cut from the rb_c0p3_hxoff checkpoint at t 9360.8, evolved LINEARLY for 200 time
       units: |phi_1|/|phi_0| = 7.765e-3 and meas/pred = 1.667 (shape 0.997) at every frame, unchanged to 1e-4 - the same
       as the nonlinear run at that time (1.667). (`figures/landau/zonal_static_gene.py`.) So a collisionless k_y = 0
       structure is linearly stationary with the non-flute part it has: the size is a property of the distribution
       (any function of the invariants on each line is stationary), not fixed by phi_0, and the 'equilibrium static
       response' is one member of that family. This is specific support for the memory reading of (h), withdrawn in (i)
       for lack of a test; with it the constancy from t 800 and the box dependence (different bursts) are consistent.
       Open: what the burst selects (1.6 reduced box, 1.0-1.4 production). First try m1 (ppj init) gave phi = 0: equal
       species init carries no charge; the first job also failed on an srun cpus-per-task clash (fixed in the sbatch).
   (r) 4 Oct, 17:00-18:00, when the jets' factor is set (`jets_factor_history.py`, `jets_factor_vs_loss.py`, rb_c0p3_hxy):
       the jets form at t 650-700 (m = 1 flute |phi_0| 1e3 -> 2.9e5, then constant to 0.3 %); right at formation the
       factor is ~1 (t 600: m = 1 1.11, m = 5-9 1.05-1.17), at t 700 1.46, then the non-flute part keeps growing on its
       own (1982 -> 2229, factor -> 1.59) and stops at t ~2100. Not the cumulative Landau loss (W doubles after t 2100
       with no change; m = 3 decreases). The approach is ~1/(t - 700): (t - 700) x (gap to final) ~ 1.5e4 for t 900-1300,
       the signature of collisionless phase mixing of a transient left by the abrupt burst. Test submitted (18:00):
       burst_rerun/legA (rb_c0p3_hxy from t 0 to 700, 4 nodes after ref leg8 on cpu-r-hbm-[29,40-42]), legB nonlinear
       700 -> 1500, lin_m1 / lin_m3 linear k_y = 0 from legA's checkpoint (1 core each, cpu-r-hbm-37). If linear
       reproduces the rise, the factor is the phase-mixed residual of the burst (a Rosenbluth-Hinton-type initial-value
       problem), computable from the state at formation.
   (s) 4 Oct, evening, paper reframed (Dan: succinct, final answer, 30 pages max, important material in the main text,
       'make this as strong as you can'): 36 -> 20 pages. Landau derivation and momentum in Sec. 3 (3.4, 3.5); new 3.6
       energy budget with dE_v/dt = c dM_v/dt (balance-figure damping term is -<(w/wt)P> = -c int f, consistent);
       predictions rewritten; Sec. 4 reordered (4.1 slow drain physical / collapse numerical first, then geometry,
       Rayleigh, saturation, Landau, jets pay); intro, abstract, cartoon panel (b) and summary tell the momentum story;
       ruled-out list, fluid-model-alone, channels, corrugation, dissipation, trapping figures and the old Rayleigh and
       trapping appendices removed (files kept). Panels (b) of amplitude and exchange removed (they read the numerical
       collapse as physics). Injection renamed I (P is the absorbed power).
   (t) 4 Oct, correction to (s): the relation dE_v/dt = c dM_v/dt for EVERY exchange (put in Sec. 3.6 in (s)) is
       wrong for the Reynolds stress - the stress divergence integrates to zero over the box, so it cannot change the
       vortices' total momentum, yet it exchanges energy (fig. balance). Delta E = c Delta p holds for particles in a
       steady pattern, not fluid-fluid exchange. Removed; 3.6 now states only what 3.5 derives (lab-frame energy of the
       damping and its sign) and that the stress returns energy, measured.
   (u) 4 Oct, late: pseudomomentum (wave-activity) route to the drain (`figures/landau/drain_pseudomomentum.py`).
       The Rayleigh mode relation zeta = Z_x psi/(U - c) holds on the measured fields (coefficient 0.84-0.92, fit
       0.93-0.97; /tmp check). Kelvin's wave activity outside the cat's eyes has the sign of -c. With the damping as a
       local vorticity source of lab-frame power (w/wt)P, the outside-eye budget predicts GROWTH at ~2.2 L/E (cuts 0.5-2
       w: 1.5-2.3 L/E); measured is decay at 0.13-0.38 L/E. So the wave-activity exchange is dominated by the cat's
       eyes (critical-layer flux), which this excludes: the cancellation moves into the critical layers, no prediction.
       NOT for the paper. Check of a paper claim (`core_fraction.py`): 76-80 % of the absorbed power is in the jet cores
       (wt, w opposite signs); lab-frame damping term +0.45-0.64 P (gain for the vortices). Claim in Sec. 3.6 holds.
   (v) 4 Oct, night: E/W split of the exchange by k_y (`full_budget_EWZ.py`, `zonal_W_vs_E.py`, `quasilinear_check.py`).
       Exact: the nonlinearity conserves E and W to < 1e-4 in every window. E side (closes to +-0.08 P1 at k_y = 1):
       jets -1.25 / -0.96 P1, vortices +0.96 / +0.69, harmonics +0.29 / +0.27 (nx 64 / nx 128, t 3000-4400).
       CORRECTION within the same hour: the Spectral_2D_W_* files are the NATIVE ENTROPY (diagnostic_schema.dat: 'W_fe =
       native entropy; E_fe = native electrostatic; FE_fe = total W'; W_fe slope + E_fe slope = fe_time total), so my
       'Z = W - E' and the 'entropy carried into the jets / quasilinear alpha ~ 1' reading were WRONG. Read correctly,
       the nonlinearity moves almost no entropy into the jets (NL of the W file at k_y = 0: -0.02 / +0.17 P1); the
       vortices' entropy goes to k_y >= 2. The entropy budget does not close (gap 0.77 P1 at k_y = 1) because 'spectral
       parallel is RAW: embedded hyp_v/hyp_z must be subtracted conditionally' (schema) - entropy-side statements need
       that subtraction and are NOT for the paper. Note also 'reduced dipole E = Eplunk/2 + Ebounce; distinct from
       native electrostatic E': all budget statements in the paper use the native electrostatic E.
   (w) 4 Oct, night: THE JETS' EXTRA NON-FLUTE FACTOR IS EXPLAINED (`figures/landau/jet_charge_pitch.py`, linear-run
       checkpoint t 9560.8, k_y = 0, m = 1). GENE grids reproduced (int F0 d^3v = 1); exact zonal Poisson
       k^2 lam^2 g^xx phi = rho holds (rho/(phi g^xx) = 0.02901 = 2 k^2 lam^2 at every z). h/F0 is orbit-constant to
       < 2 % (stationary), ~0.99 phi_f at low energy, but at high energy it rises with pitch angle (0.97 -> 1.14 phi_f
       from passing to deeply trapped). Exact decomposition of the field equation: phi_1 - rho0[phi_1] =
       -k^2 lam^2 (g^xx - <g^xx>) phi_f  [static]  +  (n[I] - <n[I]>)  [pitch-angle structure of the deposit
       I = <h> - F0 phibar]. Projected on the static shape: static 1.000 + pitch-angle part 0.689 (shape match 0.985)
       = 1.689 vs measured 1.665, residual 1.5 %. So the factor is the pitch-angle distribution of the jets' charge,
       which a collisionless zonal state keeps (linear test (q)); the 'equilibrium' response assumes a Maxwellian deposit.
       Next: what made the pitch-angle structure (candidate: the interchange burst, h ~ F0 (w - w*T)/(w - wd(eps,lam)),
       through the bounce-averaged drift; check on burst_rerun legA at t 700), and whether the vortices' critical-layer
       excess (6-8 vs 4.85) is the same thing.
   (x) 4 Oct, night, follow-up to (w): the deposit I = <h> - F0 phibar (`jet_charge_tperp2.py`) reproduces the extra
       non-flute part directly (+0.661 vs +0.689 from n[h] - rho0[phi]). Profile in pitch angle (high energy): I/F0
       rises across the passing range (lambda 0 -> 1 = 1/B_max) and is flat over the trapped range, amplitude ~ eps
       (low energy: flat). Smooth models fit the shape but not the integral: linear in mu (T_perp) R^2 0.88 -> +1.58;
       mu, mu*eps, mu^2 R^2 0.92 -> +1.24; deposit proportional to the orbit-averaged magnetic drift R^2 0.89 -> +0.23
       (`jet_charge_aniso.py`, `jet_charge_driftdep.py`). So the exact statement stands (the factor is the
       pitch-angle structure of the jets' charge, closed to 1.5 %), but its origin is not pinned by a one-parameter
       model; the burst_rerun legA checkpoint at t 700 (and lin_m1/lin_m3) is the test of where it is made.
   (y) 4 Oct, 21:30 (`figures/landau/bump_vs_jets.py`, reference run): the vortices' critical-layer bump ratio / 4.85
       is 1.27-1.28 for t 1000-3000, then 1.32, 1.39, 1.49, 1.60, 1.66-1.71 at t 4000-9800 (approaching the numerical
       collapse), while the jets' m = 1 factor is 1.56 at t 1000 and 1.62-1.66 thereafter. The bump does NOT track the
       jets' pitch-angle factor in the slow phase; it rises to the jets' value only as the cat's eyes narrow towards the
       grid / Debye scale. Separate, eye-width-dependent effect; the (w) mechanism does not explain the slow-phase bump.
   (z) 4 Oct, 21:45 (`bump_coef_freq.py`): the fluid in a cat's eye circulates at omega_tr, so its charge is forced
       at ~omega_tr, not at wt = 0. The response coefficient falls with frequency (4.86 at 0.002, 4.73 at 0.1, 4.08 at
       0.3, 3.41 at 0.5), so a finite trapping frequency LOWERS the bump below 4.85; measured is above. Ruled out.
       The slow-phase bump excess (1.27) remains unexplained; candidates: pitch-angle structure of the eye's charge
       (as (w) for the jets) - needs the distribution function at the critical layer, which a single checkpoint
       mixes with the other vortex.
   (y') correction to (y): in rb_c0p3_hxoff (no radial hyperdiffusion, no collapse) the bump ratio / 4.85 also rises,
       1.24-1.28 (t 2600-4000), 1.31 (5000), 1.44 (6000), 1.57 (7000), with the jets' factor flat at 1.66: the rise
       follows the weakening of the vortices (narrowing eyes), not the numerical collapse.
   PRE-REGISTERED PREDICTION (4 Oct 2026, 21:50, before coll_hxoff has run): with Landau collisions nu = 0.01
       (pitch-angle scattering time ~100), the jets' m = 1 non-flute factor relaxes from ~1.6 to ~1.0 (static response)
       within a few hundred time units of the restart at t 2500, because it is the pitch-angle structure of the jets'
       charge (STATUS (w)); the collisionless control stays at ~1.6. Check with `figures/landau/coll_prediction_check.py`.
   (aa) 4 Oct, 23:00: phase-mixing (Rosenbluth-Hinton-type) predictor `figures/landau/rh_predict.py`: linear k_y = 0
       dynamics conserve I = <h>_orbit - F0 phibar on every orbit (eps, mu, sigma), so the final non-flute part follows
       from any checkpoint by solving (1 + k^2 lam^2 g^xx) phi_f - rho0[phi_f] = n[I]. Validated on the stationary m = 1
       checkpoint: predicted final 1.641 vs present 1.711 (from the checkpoint's own Poisson field; field.dat gives
       1.665) - good to a few %. NOTE GENE debye2 = 2 lambda_D^2 for the pair plasma (rho = k^2 debye2 g^xx phi).
       `burst_analysis.py` runs the whole burst test (prediction from legA t 700, linear lin_m1/3, nonlinear legB,
       reference) once the jobs finish.
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

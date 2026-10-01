# Status (1 Oct 2026, 10:15)

- Notes rewritten in the house style (figures two per row at body size; simulation numbers removed) and around three
  tested results: the vortices are the saturated Rayleigh instability of their jets (frequency, structure,
  omega_tr = alpha gamma_R); their damping is linear Landau damping of their non-flute part (drift resonance
  suppressed by trapping, far below the prediction of eq. gamma_layer); they decay as the jets lose their
  instability.
- Running (submitted 1 Oct 09:00-09:15, each needs ~28 h to reach the collapse): `rb_c0p3_hxy_sinktest/` hv_x4, hv_d4,
  hxy_x4, hxy_d4 (2 nodes each, -flat) and nx128 (4 nodes, hbm). Predictions they test are stated at the end of
  Sec. 6.8 and in the summary of the notes; `figures/sinktests/figure.py` will draw them.
- Not yet done: first-principles calculation of the kinetic loss (non-flute fraction f and its Landau damping);
  theory of alpha; model for the erosion of the vorticity strips (lifetime); the same analysis for the production run.
- Next: pedagogical pass on the mathematics (every derivation from a stated starting point, undergraduate level), Dan 1 Oct.

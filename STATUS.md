# Status (1 Oct 2026, 11:15)

- Notes restructured to the house rule "only the calculation that describes the simulations, written as though the
  answer was always known": theory (Euler limit, Rayleigh instability of the jets, trapping saturation, kinetic loss,
  the decay) -> six predictions -> tests. Discarded material (drift-resonance limits, sink models, remnant power laws)
  is in the git history up to commit 35058e6.
- Appendices A (Rayleigh), B (trapping) and C (Landau damping of the non-flute part) written and verified. The loss
  formula (eq. loss in the main text) is parameter-free and matches the measured parallel-channel loss to ~10 % in the
  reduced box, the amplitude-reduced restarts and the production box (figures/landau).
- Cold read by an independent agent done; its points are applied (commit 056d407 and later).
- New run: rb_c0p3_hxoff (reduced box from noise with hyp_x = 0, 3-day legs) started 1 Oct 10:23; first readings of the
  sink tests: loss rate independent of hyp_v, decay rate follows hyp_x/hyp_y.
- Running (submitted 1 Oct 09:00-09:15, ~28 h to reach the collapse): `rb_c0p3_hxy_sinktest/` hv_x4, hv_d4, hxy_x4,
  hxy_d4 (2 nodes each, -flat) and nx128 (4 nodes, hbm); `figures/sinktests/figure.py` will draw them.
- Open: theory of alpha; rate of erosion of the jets' corrugation (lifetime); same analysis for the production run.

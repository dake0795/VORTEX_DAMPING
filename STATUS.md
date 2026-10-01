# Status (1 Oct 2026, 11:15)

- Notes restructured to the house rule "only the calculation that describes the simulations, written as though the
  answer was always known": theory (Euler limit, Rayleigh instability of the jets, trapping saturation, kinetic loss,
  the decay) -> six predictions -> tests. Discarded material (drift-resonance limits, sink models, remnant power laws)
  is in the git history up to commit 35058e6.
- Appendices A (Rayleigh) and B (trapping) written and verified; appendix C (Landau damping of the non-flute part:
  formula, evaluation against the measured loss with no free parameter) in progress.
- Running (submitted 1 Oct 09:00-09:15, ~28 h to reach the collapse): `rb_c0p3_hxy_sinktest/` hv_x4, hv_d4, hxy_x4,
  hxy_d4 (2 nodes each, -flat) and nx128 (4 nodes, hbm); `figures/sinktests/figure.py` will draw them.
- Open: theory of alpha; rate of erosion of the jets' corrugation (lifetime); same analysis for the production run.

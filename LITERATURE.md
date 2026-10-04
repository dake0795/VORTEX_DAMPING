# Literature search: Landau damping of near-flute vortices in dipole pair plasma

Search date 4 Oct 2026, web search plus the Crossref API.
VERIFIED means the bibliographic record (authors, title, journal, volume, pages, year, DOI) was checked on Crossref, or on the publisher or arXiv page, or the first page of the PDF was read.
UNVERIFIED means the reference was seen only in a search snippet.
Where a paper's content was read, this is stated ("read").

Caveat: a web search returned github.com/dake0795/VORTEX_DAMPING. Its summary reproduces this paper's claims, which means the unpublished paper is publicly indexed. It is not prior art.

---------------------------------------------------------------------------------------------------
## A. Closest prior art for claim 2: 2D modes given z-structure by geometry and damped by bounce resonance (non-neutral plasma)

The non-neutral (Penning-Malmberg) plasma community has the nearest analogues. There, an E x B mode that is 2D to leading order (diocotron mode = Kelvin wave on a vortex) acquires a weak dependence along B from the trap geometry. Particles bouncing along B at n*omega_b, resonant with the Doppler-shifted mode frequency, then damp the mode.

1. S. M. Crooks & T. M. O'Neil (1995), "Rotational pumping and damping of the m=1 diocotron mode", Phys. Plasmas 2, 355-364. doi:10.1063/1.870962. VERIFIED, read.
   **Closest single analogue.** The end potentials are non-axisymmetric in the mode frame, so a flux tube's length varies as it E x B rotates, and the 2D mode gets a z-dependent part.
   - Sec. II (adiabatic limit, omega_b >> omega_R): collisional damping with rate proportional to lambda_D^2.
   - Sec. III: a collisionless (nu -> 0) "resonant particle transport" rate at bounce-rotation resonances l*omega_R = 2n*omega_B.
   - The damping follows from conservation of canonical angular momentum. Resonant particles receive a theta-force at reflection, take same-sign radial steps and carry the mode's angular momentum outward. This is also a direct analogue of claim 3.
   - Differences from our case: the z-structure comes from external end electrodes, not from Poisson's equation with a field-line-varying metric; there is a single m=1 mode on a rigid rotor; there are no zonal jets; and the rate is not tested against simulation.
2. B. P. Cluggish & C. F. Driscoll (1995), "Transport and damping from rotational pumping in magnetized electron plasmas", Phys. Rev. Lett. 74, 4213-4216. doi:10.1103/PhysRevLett.74.4213. VERIFIED. Experimental confirmation of item 1 (collisional regime).
3. F. Anderegg, M. Affolter, A. A. Kabantsev, D. H. E. Dubin, A. Ashourvan & C. F. Driscoll (2016), "Bounce harmonic Landau damping of plasma waves", Phys. Plasmas 23, 055706. doi:10.1063/1.4946021. VERIFIED, read abstract and introduction.
   A z-variation of the potential ("squeeze") creates spatial harmonics. Particles with bounce-averaged speed v_ph/n then Landau-damp the wave at a rate proportional to V_s^2 (the square of the geometric perturbation), measured quantitatively.
   This is the same "geometry-made harmonics, then n*omega_b resonance" logic. However, the wave is a Trivelpiece-Gould wave with k_par != 0, not a flute vortex.
4. H. L. Berk & D. L. Book (1969), "Plasma wave regeneration in inhomogeneous media", Phys. Fluids 12, 649-661. doi:10.1063/1.1692529. VERIFIED. The original bounce-harmonic damping theory, as cited by item 3.
5. T. J. Hilsabeck & T. M. O'Neil (2001), "Finite length diocotron modes", Phys. Plasmas 8, 407-422. doi:10.1063/1.1340856. VERIFIED (record only). Finite-length (z-dependent) corrections to the 2D diocotron mode; contents not read.
6. A. A. Kabantsev, C. F. Driscoll, T. J. Hilsabeck, T. M. O'Neil & J. H. Yu (2001), "Trapped-particle asymmetry modes in single-species plasmas", Phys. Rev. Lett. 87, 225002. doi:10.1103/PhysRevLett.87.225002. VERIFIED.
   Companion theory: Hilsabeck & O'Neil (2003), "Trapped-particle diocotron modes", Phys. Plasmas 10, 3492-3506. doi:10.1063/1.1599356. VERIFIED, read abstract.
   Here the modes have z-structure, and passing particles "run back and forth along the field lines attempting to Debye shield" the perturbation. The damping is collisional (separatrix scattering, rate ~ sqrt(nu)), not Landau.

## B. Geometry-induced parallel structure and transit/bounce damping in toroidal fusion plasmas

7. H. Sugama & T.-H. Watanabe (2006), "Collisionless damping of geodesic acoustic modes", J. Plasma Phys. 72, 825-828. doi:10.1017/S0022377806004958. VERIFIED.
   Also: Z. Gao, K. Itoh, H. Sanuki & J. Q. Dong (2008), "Eigenmode analysis of geodesic acoustic modes", Phys. Plasmas 15, 072511. doi:10.1063/1.2956993. VERIFIED.
   The GAM is the tokamak analogue: a flute-like zonal potential picks up poloidal sidebands from geometry and is transit-resonance (Landau) damped. The difference is that GAM coupling is via geodesic curvature and compressibility, whereas in our case it is via Poisson's equation with g^xx, g^yy varying in different proportions along l. That coupling matters only at finite k^2 lambda_D^2, so it is specific to non-quasineutral (pair or non-neutral) plasma.
8. M. N. Rosenbluth & F. L. Hinton (1998), "Poloidal flow driven by ion-temperature-gradient turbulence in tokamaks", Phys. Rev. Lett. 80, 724-727. doi:10.1103/PhysRevLett.80.724. VERIFIED.
   Also: F. L. Hinton & M. N. Rosenbluth (1999), "Dynamics of axisymmetric E x B and poloidal flows in tokamaks", Plasma Phys. Control. Fusion 41, A653-A662. doi:10.1088/0741-3335/41/3A/059. VERIFIED.
   The collisionless response of a flux-surface potential, partly lost through parallel (transit/bounce) dynamics in the geometry.
9. H. Sugama & T.-H. Watanabe (2006), "Collisionless damping of zonal flows in helical systems", Phys. Plasmas 13, 012501. doi:10.1063/1.2149311. VERIFIED.
   Also: A. Mishchenko, P. Helander & A. Könies (2008), "Collisionless dynamics of zonal flows in stellarator geometry", Phys. Plasmas 15, 072309. doi:10.1063/1.2963085. VERIFIED.
   Also: P. Helander, A. Mishchenko, R. Kleiber & P. Xanthopoulos (2011), "Oscillations of zonal flows in stellarators", Plasma Phys. Control. Fusion 53, 054006. doi:10.1088/0741-3335/53/5/054006. VERIFIED.
   These treat 3D geometry setting the collisionless damping or oscillation of zonal (k_y = 0) potentials; non-zonal vortices are not treated.
10. T. H. Stix (1973), "Decay of poloidal rotation in a tokamak plasma", Phys. Fluids 16, 1260-1267. doi:10.1063/1.1694506. VERIFIED.
    Also: A. B. Hassam & R. M. Kulsrud (1978), "Time evolution of mass flows in a collisional tokamak", Phys. Fluids 21, 2271-2279. doi:10.1063/1.862166. VERIFIED.
    Magnetic pumping: a flow forced through field-line-varying geometry is damped. This is the classical "geometric damping" (the namesake of Crooks & O'Neil's "rotational pumping").
11. H. Okuda & J. M. Dawson (1973), "Theory and numerical simulation on plasma diffusion across a magnetic field", Phys. Fluids 16, 408-426. doi:10.1063/1.1694356. VERIFIED.
    A precedent for nearly-2D convective cells being damped by small-k_par parallel kinetics (in a uniform field, not by geometry).

## C. Pair plasma and dipole (closed-field-line) gyrokinetics

12. P. Helander (2014), "Microstability of magnetically confined electron-positron plasmas", Phys. Rev. Lett. 113, 135003. doi:10.1103/PhysRevLett.113.135003. VERIFIED.
13. A. Mishchenko, G. G. Plunk & P. Helander (2018), "Electrostatic stability of electron-positron plasmas in dipole geometry", J. Plasma Phys. 84, 905840201. doi:10.1017/S0022377818000193. VERIFIED.
    A finite Debye length produces a fluid mode that removes the lambda_D = 0 singularity, and the paper treats Landau damping with drift resonances. This is the closest pair-dipole kinetic precedent; it is linear only.
14. D. Kennedy, A. Mishchenko, P. Xanthopoulos, P. Helander, A. Bañón Navarro & T. Görler (2020), "Linear gyrokinetics of electron-positron plasmas in closed field-line systems", J. Plasma Phys. 86, 905860208. doi:10.1017/S0022377820000276. VERIFIED.
    Also: D. Kennedy, A. Mishchenko, P. Xanthopoulos & P. Helander (2018), "Linear electrostatic gyrokinetics for electron-positron plasmas", J. Plasma Phys. 84, 905840606. doi:10.1017/S0022377818001150. VERIFIED.
15. A. Mishchenko, D. Kennedy, P. Helander, A. Könies, G. G. Plunk, P. Xanthopoulos et al. (2023), "Gyrokinetic applications in electron-positron and non-neutral plasmas", J. Plasma Phys. 89, 935890403. doi:10.1017/S0022377823000764. VERIFIED. A review, including nonlinear pair-plasma runs.
16. M. R. Stoneking, T. Sunn Pedersen, P. Helander, H. Chen, U. Hergenhahn, E. V. Stenson et al. (2020), "A new frontier in laboratory physics: magnetized electron-positron plasmas", J. Plasma Phys. 86, 155860601. doi:10.1017/S0022377820001385. VERIFIED. The APEX motivation.
17. Closed-field-line gyrokinetics (electron-ion), all VERIFIED:
    - P. Ricci, B. N. Rogers & W. Dorland (2006), "Small-scale turbulence in a closed-field-line geometry", Phys. Rev. Lett. 97, 245001. doi:10.1103/PhysRevLett.97.245001.
    - S. Kobayashi, B. N. Rogers & W. Dorland (2009), "Gyrokinetic simulations of turbulent transport in a ring dipole plasma", Phys. Rev. Lett. 103, 055003. doi:10.1103/PhysRevLett.103.055003.
    - S. Kobayashi, B. N. Rogers & W. Dorland (2010), "Particle pinch in gyrokinetic simulations of closed field-line systems", Phys. Rev. Lett. 105, 235004. doi:10.1103/PhysRevLett.105.235004.
    These papers report entropy-mode streamers, Kelvin-Helmholtz breakup and zonal flows. No vortex-damping theory was found in them.
18. Dipole experiments, both VERIFIED:
    - D. T. Garnier, A. Hansen, M. E. Mauel, E. Ortiz, A. Boxer, J. Ellsworth et al. (2006), "Production and study of high-beta plasma confined by a superconducting dipole magnet", Phys. Plasmas 13, 056111. doi:10.1063/1.2186616. (LDX)
    - Z. Yoshida, H. Saitoh, J. Morikawa, Y. Yano, S. Watanabe, Y. Ogawa et al. (2010), "Magnetospheric vortex formation: self-organized confinement of charged particles", Phys. Rev. Lett. 104, 235004. doi:10.1103/PhysRevLett.104.235004. (RT-1; a non-neutral vortex in a dipole)

## D. Inviscid (2D, spatial-Landau) damping and critical layers: claim 1 territory, distinct from claim 2

Do not conflate these with claim 2. This is 2D damping at a radius where the shear flow resonates with the mode (the critical radius), with no parallel kinetics.

19. R. J. Briggs, J. D. Daugherty & R. H. Levy (1970), "Role of Landau damping in crossed-field electron beams and inviscid shear flow", Phys. Fluids 13, 421-432. doi:10.1063/1.1692936. VERIFIED.
20. D. A. Schecter, D. H. E. Dubin, A. C. Cass, C. F. Driscoll, I. M. Lansky & T. M. O'Neil (2000), "Inviscid damping of asymmetries on a two-dimensional vortex", Phys. Fluids 12, 2397-2412. doi:10.1063/1.1289505. VERIFIED.
21. A. A. Kabantsev, C. Y. Chim, T. M. O'Neil & C. F. Driscoll (2014), "Diocotron and Kelvin mode damping from a flux through the critical layer", Phys. Rev. Lett. 112, 115003. doi:10.1103/PhysRevLett.112.115003. VERIFIED, read abstract.
    Also: C. Y. Chim & T. M. O'Neil (2016), "Flux-driven algebraic damping of m = 1 diocotron mode", Phys. Plasmas 23, 072113. doi:10.1063/1.4958317. VERIFIED.
    Also: C. Y. Chim & T. M. O'Neil (2021), "Flux-driven algebraic damping of m = 2 diocotron mode", Phys. Plasmas 28, 092105. doi:10.1063/5.0060022. VERIFIED, from the PDF title page.
    These papers include angular-momentum transfer to particles swept around the cat's eye.

## E. Critical-layer and cat's-eye saturation (claim 1)

22. D. J. Benney & R. F. Bergeron (1969), "A new class of nonlinear waves in parallel flows", Stud. Appl. Math. 48, 181-204. doi:10.1002/sapm1969483181. VERIFIED.
23. R. Haberman (1972), "Critical layers in parallel flows", Stud. Appl. Math. 51, 139-161. doi:10.1002/sapm1972512139. VERIFIED.
24. K. Stewartson (1977/78), "The evolution of the critical layer of a Rossby wave", Geophys. Astrophys. Fluid Dyn. 9, 185-200. doi:10.1080/03091927708242326. VERIFIED.
    Note: Crossref gives the issue year as 1977; notes.bib has 1978. Check this.
25. T. Warn & H. Warn (1978), "The evolution of a nonlinear critical level", Stud. Appl. Math. 59, 37-71. doi:10.1002/sapm197859137. VERIFIED.
26. S. A. Maslowe (1986), "Critical layers in shear flows", Annu. Rev. Fluid Mech. 18, 405-432. doi:10.1146/annurev.fl.18.010186.002201. VERIFIED. (Review)
27. S. M. Churilov & I. G. Shukhman (1987), "The nonlinear development of disturbances in a zonal shear flow", Geophys. Astrophys. Fluid Dyn. 38, 145-175. doi:10.1080/03091928708219202. VERIFIED.
    This is the weakly unstable shear flow saturating by critical-layer trapping, the fluid counterpart of omega_tr ~ gamma.
28. T. M. O'Neil, J. H. Winfrey & J. H. Malmberg (1971), "Nonlinear interaction of a small cold beam and a plasma", Phys. Fluids 14, 1204-1212. doi:10.1063/1.1693587. VERIFIED. The plasma statement of trapping saturation, omega_bounce ~ gamma.
29. N. J. Balmforth, S. G. Llewellyn Smith & W. R. Young (2001), "Disturbing vortices", J. Fluid Mech. 426, 95-133. doi:10.1017/S0022112000002159. VERIFIED. Nonlinear critical layers on vortices.

## F. Wave-mean-flow interaction and momentum deposition by absorption (claim 3)

30. Fluid and atmospheric theory, all VERIFIED:
    - J. G. Charney & P. G. Drazin (1961), "Propagation of planetary-scale disturbances from the lower into the upper atmosphere", J. Geophys. Res. 66, 83-109. doi:10.1029/JZ066i001p00083. (Non-acceleration theorem)
    - D. G. Andrews & M. E. McIntyre (1976), "Planetary waves in horizontal and vertical shear: the generalized Eliassen-Palm relation and the mean zonal acceleration", J. Atmos. Sci. 33, 2031-2048. doi:10.1175/1520-0469(1976)033<2031:PWIHAV>2.0.CO;2.
    - D. G. Andrews & M. E. McIntyre (1978), "An exact theory of nonlinear waves on a Lagrangian-mean flow", J. Fluid Mech. 89, 609-646. doi:10.1017/S0022112078002773.
    - O. Bühler (2014), Waves and Mean Flows, 2nd edn, Cambridge University Press. doi:10.1017/CBO9781107478701. (DOI prefix verified via a chapter record.)
    The result used here: the mean flow is accelerated only by wave transience or dissipation. Absorbed pseudomomentum is deposited in the mean flow, so dissipation of a wave is not friction on the mean flow.
31. P. D. Killworth & M. E. McIntyre (1985), "Do Rossby-wave critical layers absorb, reflect, or over-reflect?", J. Fluid Mech. 161, 449-492. doi:10.1017/S0022112085003019. VERIFIED. Momentum exchange between a wave and the flow at an absorbing critical layer.
32. Astrophysical analogue of the sign question in claim 3 (whether resonant absorption takes energy from or gives energy to the wave, depending on the frame and the resonance), both VERIFIED:
    - D. Lynden-Bell & A. J. Kalnajs (1972), "On the generating mechanism of spiral structure", Mon. Not. R. Astron. Soc. 157, 1-30. doi:10.1093/mnras/157.1.1.
    - P. Goldreich & S. Tremaine (1979), "The excitation of density waves at the Lindblad and corotation resonances by an external potential", Astrophys. J. 233, 857. doi:10.1086/157448.
    Resonant stars at Lindblad or corotation resonances absorb a wave's angular momentum, and the sign of the wave's energy in the rotating frame decides whether the wave grows or decays.
33. Magnetised plasma: absorbed wave momentum becomes a radial (cross-field) particle displacement or current, because canonical momentum is p_y = m v_y + qBx. All VERIFIED:
    - N. J. Fisch & J.-M. Rax (1992), "Interaction of energetic alpha particles with intense lower hybrid waves", Phys. Rev. Lett. 69, 612-615. doi:10.1103/PhysRevLett.69.612. (Alpha channeling: each quantum absorbed moves the particle across the field by a fixed amount.)
    - J.-M. Rax, R. Gueroult & N. J. Fisch (2017), "Efficiency of wave-driven rigid body rotation toroidal confinement", Phys. Plasmas 24, 032504. doi:10.1063/1.4977919. (Resonant absorption gives a radial current, which gives a radial electric field and E x B rotation.)
    - I. E. Ochs & N. J. Fisch (2021), "Wave-driven torques to drive current and rotation", Phys. Plasmas 28, 102506. doi:10.1063/5.0062034.
    - I. E. Ochs & N. J. Fisch (2022), "Momentum conservation in current drive and alpha-channeling-mediated rotation drive", Phys. Plasmas 29, 062106. doi:10.1063/5.0085821.
34. Zonal-flow momentum theorems, both VERIFIED:
    - P. H. Diamond, S.-I. Itoh, K. Itoh & T. S. Hahm (2005), "Zonal flows in plasma: a review", Plasma Phys. Control. Fusion 47, R35-R161. doi:10.1088/0741-3335/47/5/R01.
    - P. H. Diamond, Ö. D. Gürcan, T. S. Hahm, K. Miki, Y. Kosuga & X. Garbet (2008), "Momentum theorems and the structure of atmospheric jets and zonal flows in plasmas", Plasma Phys. Control. Fusion 50, 124018. doi:10.1088/0741-3335/50/12/124018.
    These give the Taylor identity, Charney-Drazin non-acceleration for drift waves and zonal flows, and the role of dissipation and resonant absorption in breaking non-acceleration. This is the closest plasma-turbulence framing of claims 3 and 4.

## G. Staircases, vortices in jets, 2D Euler and tertiary instability (claim 1 context)

All VERIFIED:

35. D. G. Dritschel & M. E. McIntyre (2008), "Multiple jets as PV staircases: the Phillips effect and the resilience of eddy-transport barriers", J. Atmos. Sci. 65, 855-874. doi:10.1175/2007JAS2227.1.
36. G. Dif-Pradalier, P. H. Diamond, V. Grandgirard, Y. Sarazin, J. Abiteboul, X. Garbet et al. (2010), "On the validity of the local diffusive paradigm in turbulent plasma transport", Phys. Rev. E 82, 025401. doi:10.1103/PhysRevE.82.025401.
    Also: G. Dif-Pradalier, G. Hornung, P. Ghendrih, Y. Sarazin, F. Clairet, L. Vermare et al. (2015), "Finding the elusive E x B staircase in magnetized plasmas", Phys. Rev. Lett. 114, 085004. doi:10.1103/PhysRevLett.114.085004.
37. A. Frishman, J. Laurie & G. Falkovich (2017), "Jets or vortices: what flows are generated by an inverse turbulent cascade?", Phys. Rev. Fluids 2, 032602. doi:10.1103/PhysRevFluids.2.032602.
    Also: F. Bouchet & E. Simonnet (2009), "Random changes of flow topology in two-dimensional and geophysical turbulence", Phys. Rev. Lett. 102, 094504. doi:10.1103/PhysRevLett.102.094504.
38. A. Hasegawa & K. Mima (1978), "Pseudo-three-dimensional turbulence in magnetized nonuniform plasma", Phys. Fluids 21, 87-92. doi:10.1063/1.862083.
39. Tertiary (Kelvin-Helmholtz) instability of zonal flows:
    - B. N. Rogers, W. Dorland & M. Kotschenreuther (2000), "Generation and stability of zonal flows in ion-temperature-gradient mode turbulence", Phys. Rev. Lett. 85, 5336-5339. doi:10.1103/PhysRevLett.85.5336.
    - E. Kim & P. H. Diamond (2002), "Dynamics of zonal flow saturation in strong collisionless drift wave turbulence", Phys. Plasmas 9, 4530-4539. doi:10.1063/1.1514641.
    - H. Zhu, Y. Zhou & I. Y. Dodin (2020), "Theory of the tertiary instability and the Dimits shift from reduced drift-wave models", Phys. Rev. Lett. 124, 055002. doi:10.1103/PhysRevLett.124.055002.
    - The same authors' JPP companion (2020): J. Plasma Phys. 86, 905860405. doi:10.1017/S0022377820000823.
    A search snippet (source not identified, UNVERIFIED) states that gyrokinetic Kelvin-Helmholtz modes are "strongly Landau damped" at finite k_par/k_perp. If it is wanted, trace it to its source before citing.

---------------------------------------------------------------------------------------------------
## NOVELTY VERDICTS

**Claim 1 (Rayleigh-unstable corrugated jets, then vortices saturated by cat's-eye trapping).**
Each ingredient is textbook (sections D, E, G). Applying them to pair-plasma gyrokinetics with Debye-scale corrugation was not found in prior work. It is novel as an application, not as theory.

**Claim 2 (metric-forced non-flute part of relative size k^2 lambda_D^2, bounce/transit Landau damping of a near-2D vortex, parameter-free rate checked against simulation).**
No prior paper doing this was found. The mechanism class is known, though, and the paper should cite it:
- Crooks & O'Neil (1995) computed collisionless bounce-rotation-resonant damping of a 2D E x B mode whose z-structure comes from geometry (trap ends), with damping set by angular-momentum conservation.
- Anderegg et al. (2016) and Berk & Book (1969) give bounce-harmonic Landau damping from geometry-induced harmonics.
- The GAM (Sugama & Watanabe 2006) is the tokamak analogue of a flute potential coupled to parallel structure by geometry and then transit-damped.

What appears new:
- (i) The forcing is through Poisson's equation with g^xx(l) and g^yy(l) varying in different proportions, so it exists only at finite k lambda_D. This needs non-quasineutral plasma and has no quasineutral fusion counterpart.
- (ii) It is applied to self-organised box-scale vortices embedded in zonal jets, with the Doppler shift by U(x) entering the resonance.
- (iii) The rate is parameter-free and verified against gyrokinetic simulations across amplitude, box size and lambda_D.

Present it as a new instance of a known class ("rotational pumping / bounce-harmonic damping"), not as a new kind of damping.

**Claim 3 (resonant particles absorb k P / omega_tilde, a binormal force drives a radial current, momentum goes into the jets, and damping is not friction).**
The general principle is well known:
- non-acceleration and dissipation-driven mean-flow change (Charney & Drazin; Andrews & McIntyre; Killworth & McIntyre; Bühler);
- in magnetised plasma, absorbed wave momentum gives a radial current, then E_r and rotation (Fisch & Rax 1992; Rax, Gueroult & Fisch 2017; Ochs & Fisch 2021, 2022);
- Crooks & O'Neil (1995) and Kabantsev et al. (2014) do this explicitly for diocotron and Kelvin modes, with canonical angular momentum passed to particles.

Applying it to vortex-to-zonal-jet transfer via Landau damping, and the consequence that damping can feed energy from the jets to the vortex in the lab frame, was not found. It is novel as an application and consequence; the principle is not new.

**Claim 4 (vortex drain = small difference of a Reynolds-stress exchange and a damping-induced exchange).**
Not found in prior work. Conceptually it is close to Diamond et al. (2008) momentum theorems and to Lynden-Bell & Kalnajs or Goldreich & Tremaine (wave energy sign in a rotating frame), which should be cited for framing.

Not found despite searching: any paper on Landau damping of vortices or convective cells in dipole or pair plasma, or any metric-anisotropy (g^xx/g^yy) mechanism for non-flute forcing.

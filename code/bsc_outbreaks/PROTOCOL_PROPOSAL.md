# Validation protocol proposed with the outbreak dataset

Section 5 of the working notes of the outbreak dataset (1 October 2026), called "the original protocol" in `../bsc_validation/PREREGISTRATION.md`, which is based on it and lists every amendment (Supplementary Table S24 of the article). It is released so that those amendments can be checked. This is a redacted copy: every change is listed in `../REDACTIONS.md`. References to other sections and tables of the working notes point to material that is not released; the records and statistics they refer to are in `../../data/bsc_outbreaks/`.

## 5. Validation protocol (fixed before any model run)

### 5.1 What the model has to output

**Single-exposure events** (all A-grade records; exposure of minutes to 10 h). The incubation period is much
longer than the event, so there is exactly one infectious person and no second generation. The quantity to
compute is the first-generation infection probability of each susceptible j at seat r_j, given the index at
seat r_0 and the event duration T:

    p_j(T) = 1 - E[ exp( -beta * Integral_0^T m(t) q(X_j(t)) 1{X_j(t) = X_0(t)} / N(X_j(t), t) dt ) ]

where X_0, X_j are the lattice walks of the index and of j. This is the per-cell hazard
lambda(r) = beta q(r) I_r/N_r with I = 1. To first order in the hazard, p_j is a function of the
co-occupation integral Integral_0^T Sum_r q(r) P(r,t|r_0) P(r,t|r_j) dt, which the spectral
decomposition [...] delivers in closed form. Two implementation requirements follow:

1. The simulation code must not let people infected during the event infect others during the same event
   (the SIR model makes new cases infectious immediately). Run with onward transmission switched off, or
   add a latent state.
2. Seats, index seat and T are inputs taken from the record, not fitted.

**Multi-generation records** (O1 Seoul call centre, C2 Jerusalem school). The day-scale SIR dynamics are
appropriate, with a working-hours schedule. A latent period is needed for any comparison with onset dates;
without it only final attack rates and spatial contrasts can be compared.

### 5.2 Models to be compared

- **M0, well-mixed Wells-Riley (null).** p = 1 - exp(-Lambda_k T_j), one infectiousness parameter per
  event, no spatial structure. [...] The spatial model has to
  beat it on the same data.
- **M1, same-site model.** Walkers are people; transmission only between people on the
  same site; relative profiles D(r)/D0 and q(r) fixed at the default values for the mapped scene. Free
  parameters: beta and D0 [...].
- **M2, framework with a non-local term.** The index emits quanta
  that perform the same confined lattice walk with diffusion D_air and are removed at rate kappa
  (ventilation + deposition + decay, taken from the measured or estimated air-change rate of the record);
  seated susceptibles absorb in proportion to the local quanta occupation. Limits: D_air -> infinity
  recovers M0, D_air -> 0 recovers same-site transmission. The propagator and spectral results [...]
  apply unchanged because the walker dynamics are the same.

### 5.3 Parameters

| Parameter | Status | Source |
|---|---|---|
| Room geometry, seat positions, index seat, T, occupancy | fixed | record |
| Removal rate kappa | fixed where measured (R1 0.56-0.77 per h; T2 1.72 and T3 3.22 L/s per person; H1 0.3-1.0 per h estimated), otherwise a scene prior | record |
| Infectiousness multiplier epsilon_k | one per index case, fitted or integrated out | per-hour hazards differ by more than an order of magnitude between events (section 3.3) |
| Mixing coefficient (D0 in M1, D_air in M2) | shared within a setting class (vehicle cabin; seated room) | calibration events |
| Near-field enhancement (same-site q relative to neighbours) | shared | calibration events |
| Mask factor | fixed from literature if used (C1, C3) | not fitted |

T2 and T3 share one index case within the same afternoon, so they must be fitted with the same epsilon:
the pair tests the predicted dependence on duration (150-200 min against 60 min) and ventilation
(1.72 against 3.22 L/s per person) with infectiousness held fixed.

### 5.4 Calibration and hold-out split

Primary analysis, fixed split:

| Role | Records | What is used |
|---|---|---|
| Calibration, vehicle cabins | T4 high-speed trains | 23-cell offset matrix and per-hour slopes (pooled over 2,334 index cases, so epsilon averages out) |
| Calibration, seated rooms | R1 Guangzhou restaurant | per-table counts and exposure times, measured ventilation |
| Held-out shape tests (epsilon_k profiled out) | T1 bus (two published splits), T2 coach (front/rear, left/right), T5 flight (near/far), C1 classroom (front/back rows), O1 call centre (wing contrast), O4 meat plant (after digitization) | stratified counts |
| Held-out level tests (epsilon from the predictive distribution of calibration events) | H1 choir, H2 church, R2 Jeonju, T3 minibus (epsilon tied to T2), O2 Zurich office, C5 fitness classes | overall counts |
| Negative controls | S1 supermarket customers 0/8,224; C3 masked classrooms 5/728; T5 adjacent cabin 0/35; O1 other call-centre floors 1/595; T7 no rail cluster | upper bounds |
| Multi-generation | O1 (floor attack rate, onset growth rate), C2 (grade profile, cases per class) | day-scale model |

Secondary analysis: leave-one-event-out cross-validation over the eight A-grade records, shared parameters
refitted each time.

### 5.5 Likelihood and metrics

- **Likelihood.** Binomial per stratum, k_s ~ Bin(n_s, p_s). For T4, where only percentages and CIs are
  available, a normal likelihood on the logit scale with standard error (logit(hi) - logit(lo))/3.92;
  replace with binomial once cell denominators are obtained.
- **Shape.** For each stratified held-out event, z = (log RR_pred - log RR_obs)/SE(log RR_obs), with the
  observed risk ratios and CIs of Table 3. Pass if |z| < 1.96.
- **Level.** Observed count inside the 90% posterior-predictive interval; coverage over held-out events;
  slope of log predicted against log observed cumulative hazard.
- **Against the null.** Difference in expected log predictive density between M1 or M2 and M0, summed over
  seat-resolved strata, with its standard error.
- **Clustering (O1).** Distribution, over index seats drawn uniformly in the north wing, of the wing
  contrast and of the join-count statistic of section 3.4; compare with the observed values.
- **Negative controls.** Posterior-predictive probability of exceeding the one-sided 95% upper bound of
  Table 5; any probability above 0.5 is a failure.

### 5.6 Decision rule for the wording of the abstract

A statement of agreement with documented outbreaks is supported only if all four hold:

1. at most one held-out shape test fails;
2. at least 80% of held-out level events fall inside their 90% predictive interval and the calibration
   slope lies in [0.67, 1.5];
3. no negative control fails;
4. the spatial model improves on M0 by more than two standard errors of the log-predictive-density
   difference.

If 1 and 3 hold but 2 or 4 fails, the supported statement is that the model "reproduces the direction and
size of documented near/far contrasts after calibration of two parameters". If 1 fails, the claim should be
limited to qualitative consistency. In every case the abstract can name only scenes for which a comparison
was run: with the present data that means office (O1, O2), classroom (C1, C2, C3) and vehicle cabins as a
stand-in for the metro carriage; for the supermarket only an upper-bound check (S1) is possible.

### 5.7 What the data can and cannot discriminate

The stratified events are not uniform in kind (Table 3). The flight, the trains and the classroom show
steep proximity effects (risk ratios 2.8 to 7.3; adjacent-seat train risk 20 to 40 times the 2-3 row
level). The Zhejiang bus, the Hunan coach and the choir are close to well mixed (risk ratios 1.0 to 1.8,
with CIs that include 1). The restaurant is organised by air-conditioning zone, not by distance. A model
with one isotropic kernel therefore has to explain both regimes through its inputs, and the only input
that differs systematically is air removal: cabins with high air-change rates show near-field patterns,
recirculating or unventilated rooms show room-wide spread. M2 makes that a quantitative prediction (the
contrast grows with kappa L^2 / D_air); M1 has no parameter that can produce it. Two records (R2 Jeonju,
H2 Sydney) show directed transport along an air jet or within a sector; an isotropic walk cannot
reproduce them and they are designated stress tests, not pass/fail tests.

### 5.8 Data gaps to close before the run

1. Cell denominators of Hu et al. Table 1 (supplement or authors).
2. Seat-level digitization of the published seat maps: Shen 2020 (bus 2), Ou 2022 Fig. 2 (coach and
   minibus), Khanh 2020 Fig. 1 (cabin), Lam-Hine 2021 Fig. 1 (classroom), Li 2021 Fig. 1 (table
   coordinates), Guenther 2020 Fig. 3A and Appendix Table S2.
3. Full texts for the six abstract-only sources.
4. A physical scale for the Seoul floor plan (building footprint), and the working hours of the call centre.
5. [...]

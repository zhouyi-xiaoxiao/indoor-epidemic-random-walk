# Pre-registration: validation of the lattice random-walk SIR model against documented outbreaks

[...] Written 2026-10-01, 05:00-05:40 BST, before any
transmission model was run against any outbreak count. Version 1.0. The file `PREREG_FREEZE.json` records the
SHA-256 of this file and of the event-geometry module at the moment of freezing.

## 0. Status and limits of this pre-registration

1. **Data-aware, model-blind.** The outbreak counts (Table 1 and Table 3 of
   `../bsc_outbreaks/memo.md`) were assembled, and read by the analyst, before this document was written.
   Blinding to the data is therefore impossible. What is fixed here, before any model output exists, is the
   observation model for every event, the three models, the calibration/hold-out split, the estimation
   procedure, the metrics, the decision rule and the list of sensitivity analyses. No model has been fitted,
   simulated or evaluated on any record at the time of writing.
2. **Protocol source.** Section 5 of `../bsc_outbreaks/memo.md` (the "original protocol") is the basis. The
   re-check of the dataset raised binding objections. Section 2 below lists every amendment made in
   response and every further deviation forced by missing data, with the reason.
3. **Post-hoc log.** Any choice made after a model has been run against data is recorded in `POSTHOC_LOG.md`
   with a timestamp, and is never used to upgrade the primary conclusion.
4. **Power addendum.** The power of the decision rule depends on the calibrated parameters. It is computed
   after calibration (calibration records only) and before any hold-out comparison, and recorded in
   `PREREG_ADDENDUM_A1_power.md`. The script that computes it reads hold-out geometry, stratum sizes and event
   totals, and does not read the hold-out stratum counts.

## 1. Question

[...] The question here: after calibration on designated records, does model M1,
or model M2 with the non-local term, predict the spatial pattern and the
level of infection in held-out documented outbreaks better than a well-mixed model, and which sentence about
agreement do the results license?

## 2. Amendments to the original protocol (all made before any model run)

### 2.1 Required by the verifier

| # | Verifier objection | Amendment |
|---|---|---|
| V1 | Guangzhou restaurant counts are an upper bound: Lu et al. state that at least one member of each of families B and C was infected at the lunch; the others may be within-family. | Table B contributes the interval-censored likelihood P(1 <= K_B <= 3 of 4), table C P(1 <= K_C <= 2 of 7). Sensitivity: upper bound (3, 2) and lower bound (1, 1). Table A (the index's own family) is excluded, as in the source. |
| V2 | Train cohort: companions (family, friends) inflate risk, most of all in the adjacent seat. | The adjacent-seat cell (same row, 1 seat apart) is excluded from the primary calibration likelihood. Sensitivity S-adj includes it; S-row excludes the whole index row. No "near-field enhancement" fitted on the adjacent cell is transferred to hold-out events in the primary analysis. |
| V3 | Three grade-A records do not meet the grade-A criteria (C1: no hours, no room size; T4: pooled index cases, no printed denominators; H1: no positions). | C1: room size, hours and ventilation enter as explicit priors (section 4) and the shape test conditions on the class total, so hours enter only through the people-walk model M1. T4: used as a pooled cohort with reconstructed cell counts (V5). H1: used for level only; no shape claim. |
| V4 | The shape criterion (|z| < 1.96 on the log risk ratio) is passed by any constant near/far ratio between 1.16 and 3.19; the level criterion is passed by a constant-hazard model with no spatial structure; only "beat the null by 2 SE" has power. | Decision rule replaced (section 8): (A1) the spatial model must beat the well-mixed null on a pooled, exact conditional likelihood-ratio test at one-sided alpha = 0.05; (A2) every primary shape event must lie inside the model's own 95% conditional predictive region; each event is classified as discriminating or non-discriminating; level claims are restricted to tests with defined power (section 7). Power is computed and stated before unblinding (addendum A1). |
| V5 | Train table has no cell denominators. | The verifier's reconstructed cell counts (`../bsc_outbreaks/verify/results/train_matrix_reconstruction.csv`, 230/72,692 against the reported 234/72,093) are used with a binomial likelihood. Reason: three cells have a printed lower limit of 0.00, for which the logit-normal likelihood of the original protocol is undefined. Sensitivity S-logit: logit-normal on the 20 cells with a positive lower limit. |
| V6 | "Three regimes" is not statistically established (heterogeneity p = 0.29). | No analysis below assumes the three-regime reading. |
| V7 | The call-centre wing contrast depends on 10 unmapped cases. | The wing contrast is not used as a pass/fail test for the spatial models (section 5.6); the within-wing join-count test is used. |

### 2.2 Forced by data gaps (section 5.8 of the original protocol was not closed)

| # | Gap | Consequence |
|---|---|---|
| D1 | Seat maps of Shen 2020 (bus), Ou 2022 (coach, minibus), Khanh 2020 (cabin), Lam-Hine 2021 (classroom) and Li 2021 (tables) are not digitized. | Zonal observation models: the published stratum sizes and the stated index seat are used; the unknown occupied seats inside each zone are drawn at random and marginalised (300 configurations). Nothing is tuned. |
| D2 | O4 (German meat plant) appendix counts not transcribed. | O4 is not evaluated. |
| D3 | C2 (Jerusalem school) needs the multi-room scene-coupling extension, which is not part of the verified simulator, and has no per-class denominators. | C2 is not evaluated. |
| D4 | H2 (Sydney church) has no room dimensions. | H2 is excluded from quantitative level tests. It stays a qualitative stress case together with R2 (directed air jet), as in the original protocol. |
| D5 | Square lattices with a = 0.5 m (vehicles) cannot place one seat per site for 3+2 seating at 0.75-1.0 m row pitch. | Vehicles use a rectangular lattice with one site per seat position and physical pitches (a_x, a_y); hop rates are D/a_x^2 and D/a_y^2, so the continuum limit is isotropic diffusion with coefficient D in m^2/h. Seated rooms use the square lattice a = 1 m of the original mapping. |
| D6 | The index seat on flight VN54 is not in the locally stored text. | Primary: index seat uniform over the 28 business seats. Sensitivity S-5K: seat 5K (row 5, right window), as recalled from the source figure; not verified. |

### 2.3 Other choices fixed now

- O1 (Seoul call centre) is treated by epidemic percolation (final size of a fixed-dose multi-generation chain),
  not by the day-scale SIR with latent period. Reason: the test statistic is a final spatial pattern; a
  generation-based final-size model needs no latent-period or working-hours parameters. The onset-curve
  comparison of the original protocol is dropped (it would need those parameters and they are not reported).
- Infection probabilities of M1 are computed from the pair-level mean co-occupation time (closed form), then
  checked against the verified exact individual-based simulator at the calibrated values (section 3.2).

## 3. Models

Notation: lattice sites r of a rectangular room with reflecting walls, G(r, s | r0) the continuous-time
lattice-walk propagator with diffusion coefficient D (probability that a walker started at r0 is at r after
time s), index case at r0 for the whole event, susceptible j at seat r_j for a duration T_j, one infectious
person and no second generation within an event (original protocol 5.1).

All three models have the form p_j = 1 - exp(-eps_k * X_j), with one infectiousness parameter eps_k per event
(per index case) and a model-specific exposure X_j:

### 3.1 M0, well-mixed Wells-Riley (null)

X_j = (1/M) [ T_j/kappa - (1 - exp(-kappa T_j))/kappa^2 ], M the number of lattice sites (V = M a_x a_y H).
With eps = E b / (a_x a_y H), E the quanta emission rate and b = 0.5 m^3/h a reference breathing rate, this
is the transient Wells-Riley equation. For a shape test with equal durations X_j is constant, so the
conditional distribution of cases over strata is hypergeometric.

### 3.2 M1, same-site model

People are the walkers; transmission only between people on the same site (hazard
lambda = beta q I_r / N_r); all people in these events are seated, so D(r) = D_p and q are uniform over the
seating area. Exposure is the expected co-occupation time of the pair (index, j):

X_j = Integral_0^T Sum_r G(r, t | r0) G(r, t | r_j) dt = Integral_0^T G(r_j, 2t | r0) dt,

the co-occupation integral of the original protocol (5.1), in closed form by the cosine-mode expansion [...]. eps_k = beta q E[1/N_r], the pair hazard while co-located. Free parameters: eps_k (per event) and
D_p (shared within a setting class). M1 has no ventilation parameter.

The exponentiated mean dose 1 - exp(-eps X) is exact to first order in eps and an upper bound otherwise. Check,
fixed now: at the posterior-mode parameters each hold-out event is simulated with the verified exact simulator
(`../bsc_sim`, per-cell rule, only the index infectious, no recovery, 4,000 replicates). If any
stratum attack rate differs from the closed form by more than 0.03 in absolute value, the simulated
probabilities replace the closed form for M1 in the primary tables.

### 3.3 M2, framework with the non-local term

The index emits quanta at rate E; quanta perform the same lattice walk with diffusion coefficient D_air and are
removed at rate kappa (ventilation + deposition + inactivation); a susceptible inhales in proportion to the
quanta occupation of its site:

X_j = Integral_0^T (T - s) exp(-kappa s) G(r_j, s | r0) ds, eps = E b / (a_x a_y H).

D_air -> infinity gives M0 exactly (M0 is nested in M2). Free parameters: eps_k (per event), D_air (shared
within a setting class). kappa is not fitted: measured where available, otherwise a prior that is integrated
out. In the primary analysis there is no separate near-field factor (amendment V2).

### 3.4 Fixed constants and priors (not fitted)

- Deposition 0.3 /h + inactivation 0.63 /h = 0.93 /h added to every ventilation rate (values used by Ou et al.
  2022 after Miller et al. 2021). kappa = ACH + 0.93.
- Ventilation: T2 4.8 ACH and T3 9.6 ACH (measured, Ou 2022; lognormal SD 0.2); R1 uniform 0.56-0.77 ACH
  (measured); H1 uniform 0.3-1.0 ACH (estimated by Miller 2021). Unmeasured, log-uniform priors: T4 trains 5-20;
  T1 bus (recirculation mode) 0.5-5; T5 aircraft cabin 10-30; C1 classroom (HEPA unit, doors and windows open)
  2-12; O1 office 1-6; O2 office 1 ACH (source statement, lognormal SD 0.3); R2 0.1-0.5; C5 0.5-5; C3 1-6;
  S1 1-4; T7 metro 5-20. These ranges are the analyst's engineering judgement, not measurements.
- Flat prior on log10 D over [0.01, 10^4] m^2/h for D_p and D_air (61 grid points).
- Train co-travel time: gamma distribution with mean 2.1 h and SD 1.8 h (source), 16 quadrature nodes.

## 4. Observation model per event

Coordinates: x across the cabin or room, y along it. "Marginalised" means averaged over 300 random
configurations (seed fixed).

| Event | Role | Lattice and pitches | Index | Susceptibles and strata | T | Notes |
|---|---|---|---|---|---|---|
| T4 trains | calibration, vehicle class | 17 rows x 6 columns (A B C aisle D F); a_x = 0.5 m, a_y = 1.0 m | every seat in turn | 23 offset cells (rows apart 0-3, columns apart 0-5); model value = mean over all seat pairs with that offset | gamma(2.1 h, 1.8 h) | row pitch 1.0 m assumed (source's 0.4 m is not a seat pitch); sensitivity S-pitch uses 0.4 m |
| T1 Zhejiang bus | hold-out shape | 15 rows x 6 columns (3 seats, aisle, 2 seats); a_x = 0.42 m, a_y = 0.75 m; H = 2.0 m | row 8, middle seat of the 3-seat side | rows 6-10: 23 of 24 other seats occupied; rows 5 and 11: all 10; other rows: 34 of 40; primary split rows 5-11 (33) vs rest (34); secondary split rows 6-10 (23) vs rest (44) | two rides of 50 min | bus width 2.5 m assumed |
| T2 Hunan coach | hold-out shape | 13 rows x 5 columns (A B aisle C D; seat 13E on the aisle site of row 13); a_x = 0.5 m, a_y = 0.877 m; V = 60.42 m^3 | 12D | driver side (C, D) fully occupied except 8C; opposite side: 2 empty seats in rows 1-7, 4 in rows 8-13; primary split rows 8-13 (19) vs rows 1-7 (26); secondary split driver side (24) vs opposite side (20), 13E ignored | 200 min (Ou 2022) | sensitivity S-Luo: 150 min |
| T3 minibus | tied-index level | 6 rows x 4 columns; V = 21.69 m^3 | row 4 | 17 persons (Ou 2022), seats marginalised | 60 min | eps tied to T2 |
| T5 flight VN54 | hold-out shape | business cabin 7 rows x 6 lateral sites (A, aisle, D, G, aisle, K), a_x = 0.9 m, a_y = 1.1 m, embedded in a 46-row cabin lattice (2 rows forward, 3 rows service area, 5 rows premium economy, 2 rows service, 27 rows economy) | uniform over 28 seats (D6) | 20 others on random seats; near = the 12 occupied seats closest to the index, far = the other 8 | 10 h | seat pitch assumed |
| R1 Guangzhou restaurant | calibration, seated-room class | 17 x 8 sites, a = 1 m, H = 3.14 m | table A | 18 tables on a 6 x 3 grid of table positions (pitch 2.83 m x 2.7 m); A, B, C adjacent along the window wall, T18 across from A; A's column and the positions of the 14 remote tables marginalised; patrons at their table's site; per-table overlap time with table A | per table (23-82 min) | table coordinates not digitized (D1) |
| C1 Marin classroom | hold-out shape | 11 x 13 sites, a = 1 m; desks at 2-site pitch, 5 rows x 5 columns | teacher at the front centre, one site in front of the first row's gap | front two rows 10 of 10 desks; back three rows 14 of 15 (empty desk marginalised) | 13 h (two symptomatic school days x 6.5 h); matters for M1 only | sensitivity S-corner: teacher's desk in a front corner |
| O1 Seoul call centre | hold-out shape, multi-generation | north wing: 137 digitized desks, coordinates in desk pitches, pitch uniform 1.0-1.6 m; rectangular hall = bounding box + 1 pitch | uniform over north-wing desks | epidemic percolation with p_ij = 1 - exp(-eps X_ij); eps set so that the mean final size of outbreaks reaching at least 40 cases is 79; outbreaks with 69-89 cases retained; statistic: join-count z-score (adjacent desks within 1.25 pitch, exact permutation moments) | steady state (8-h days) | 10 of 94 cases unmapped |
| H1, O2, R2, C5 | level | volume only: 810; 700 m^2 x 2.7 m; 96.6 m^2 x 3.2 m; 60 m^2 x 3 m | n/a | 52/60 (alt. 32/60); 8/12 (alt. 10/12); 2/13; 57/217 (8 instructors, two 50-min classes per student assumed) | 2.5 h; 11 h; 21 min; 2 x 50 min | ceiling heights 2.7, 3.2 and 3 m assumed |

Hold-out stratum counts are stored in `data/heldout_outcomes.json` and are read only by the evaluation
scripts (`03_*` onward).

## 5. Calibration and hold-out split (frozen)

### 5.1 Calibration

- Vehicle-cabin class: T4. Parameters: D_p and eps (M1); D_air and eps (M2); eps (M0). Likelihood: binomial on
  22 cells (adjacent-seat cell excluded). kappa_T4 integrated over its prior.
- Seated-room class: R1. Parameters as above. Likelihood: section 2.1 V1, 17 tables.
- Estimation: profile likelihood over eps on the log10 D grid; posterior weights proportional to the profile
  likelihood (flat prior on log10 D); 400 posterior draws carry parameter uncertainty into every prediction.
- Reported for the calibration fits: maximum log-likelihood, AIC and BIC of M0, M1, M2 (in-sample, with the
  caveat that kappa priors are not counted as parameters), deviance goodness of fit, and the fitted kernels.

### 5.2 Hold-out, primary shape tests (eps_k profiled from the event total; never from the strata)

T1 (rows 5-11 vs rest), T2 (rear vs front), T5 (near vs far), C1 (front two rows vs back three) and O1
(within-wing join count). Secondary splits (T1 rows 6-10; T2 driver side vs opposite side) are reported and not
counted in the decision rule.

### 5.3 Hold-out, level checks

- Tied index: eps fitted on T2 (7/46, 200 min), prediction for T3 (2/17, 60 min) with the measured volumes and
  ventilation rates.
- Unselected cohort transfer: eps distribution anchored on T4 (mean over 2,334 index cases), prediction for
  C3 (masked Utah classrooms, 5/728).
- Selected single-index outbreaks H1, O2, R2, C5: implied emission rate per event, and the dispersion test of
  section 7.

### 5.4 Negative controls

S1 supermarket customers (0/8,224; bound 0.0364%), C3 masked classrooms (5/728; bound 1.44%), T5 premium
economy (0/35; bound 8.2%), T7 no rail cluster (probability of at least 5 secondary cases from one ride). O1
other floors and the T1 unexposed bus contain no modelled exposure and are not counted. R1 remote tables are
part of the calibration likelihood and are not counted.

### 5.5 Secondary analysis

Leave-one-event-out over the stratified grade-A events of each class (vehicles: T4, T1, T2, T5; rooms: R1,
C1): D refitted on the remaining events (joint profile likelihood, event-specific eps), conditional log score
of the left-out event.

### 5.6 What is not tested

O1 wing contrast (79/137 against 4/62): the two wings are separate rooms on the published plan. A lattice with
a wall reproduces any contrast if the door coupling is free, so the contrast is not counted for M1 or M2. It
is reported as evidence against a single well-mixed floor only.

## 6. Shape metrics

For an event with strata "near" (n1 people) and "far" (n2), total K cases:

1. **Conditional predictive distribution.** For each posterior draw (D, kappa, seat configuration, index seat),
   eps is set so that the expected total equals K; the exact distribution of k_near given the total K follows
   from the two Poisson-binomial laws. The predictive pmf is the average over draws. This is the spatial
   pattern score: S = log pmf(k_near observed).
2. **Two-sided predictive p-value**: sum of pmf over outcomes no more probable than the observed one.
3. **Predicted risk ratio** with 90% predictive interval (parameter uncertainty and demographic noise), and
   the z statistic of the original protocol, z = (log RR_pred - log RR_obs)/SE(log RR_obs), for continuity.
4. **Against the null**: Delta = S(model) - S(M0) per event and pooled over the four single-event shape tests.

Event labels: *reproduced, discriminating* (model p >= 0.05, M0 p < 0.05); *reproduced, non-discriminating*
(both >= 0.05); *not reproduced* (model p < 0.05).

## 7. Level metrics

Levels of published outbreaks cannot be predicted from a typical-index distribution: investigated outbreaks are
selected for high infectiousness. The following is therefore fixed:

1. Implied emission rate E_k per event under M0/M2 (identical for a group spread over the room, because the
   spatial mean of X is the M0 value) and implied pair hazard under M1, with exact binomial 95% limits.
2. Dispersion test: SD of log implied infectiousness across H1, O2, R2, C5 under (a) the Wells-Riley scaling
   T/(V kappa) with transient, (b) the M1 scaling T/M, (c) a constant hazard per hour. A model's physical
   scaling is supported only if its SD is smaller than that of (c) with one-sided p < 0.05 (permutation of the
   structural factors across events). The power of this test is stated in addendum A1.
3. Tied-index check: observed T3 count inside the 90% predictive interval of the T2-fitted model.
4. Unselected transfer: observed C3 count inside the 90% predictive interval of the T4-anchored model.
   Infectiousness between index cases: lognormal with sigma = 0.94 (the 10-400 quanta/h range [...]
   read as a 95% range), mean fixed by T4. Masks: factor 0.35 on the dose (assumed, range 0.2-0.6;
   also reported without masks).
5. Negative controls: probability, over the predictive distribution, that the predicted attack rate exceeds
   the one-sided 95% upper bound; above 0.5 is a failure. S1 assumes that every screened customer made one
   27-minute visit while one infectious employee was in the store (maximal exposure); a failure under this
   assumption is reported as inconclusive, not as a model failure.

The original "at least 80% inside the 90% interval and slope in [0.67, 1.5]" is reported for continuity but
carries no weight in the decision rule (objection V4).

## 8. Decision rule for the wording (replaces original 5.6)

Evaluated separately for M1 and M2.

- **A1 (discrimination).** Pooled Delta over T1, T2, T5, C1 is positive and its one-sided Monte Carlo p-value
  under M0 (hypergeometric resampling of every event, 20,000 replicates) is below 0.05.
- **A2 (adequacy).** Two-sided predictive p >= 0.05 for each of T1, T2, T5, C1, and the observed O1 join-count
  z inside the central 95% predictive interval.
- **A3 (level checks).** Tied-index and unselected-transfer checks inside their 90% intervals.
- **A4 (negative controls).** No failure.

| Outcome | Sentence licensed |
|---|---|
| W1: A1, A2 with no failure, A3, A4, and power of A1 at least 0.8 | "good out-of-sample agreement with the spatial patterns of documented outbreaks" (never "excellent": five events) |
| W2: A1 and A4, at most one failure in A2 | "after calibration of [number] shared parameters, with one infectiousness parameter per event, the model reproduces out-of-sample the near/far contrasts of [events]; it does not capture [event]" |
| W3: A2 with at most one failure, A1 fails | "consistent with the documented spatial contrasts, but not distinguishable from a well-mixed model on these data" |
| W4: two or more failures in A2 | "reproduces [events] and fails [events]"; no general agreement claim |
| W5: A4 fails, or every discriminating event is not reproduced | no agreement claim |

In every case: the sentence names only settings for which a comparison was run; vehicle cabins are named as
stand-ins for the metro carriage; the supermarket is covered by an upper-bound check only; level agreement is
never claimed beyond A3. If only M2 passes, the non-local term must be part of the model
before the sentence may be used.

## 9. Sensitivity analyses (fixed list; reported, never used to upgrade the primary outcome)

S-adj, S-row (train cells); S-logit (train likelihood); S-pitch (train row pitch 0.4 m); S-aniso (vehicle
class: hop rate across rows multiplied by a fitted factor theta <= 1, a seat-back barrier in the sense of a
heterogeneous D); S-R1-upper, S-R1-lower (restaurant counts); S-5K (flight index seat); S-Luo (coach
150 min, 7/48); S-corner (classroom index); S-H1-confirmed (choir 32/60); S-nu (adjacent-seat near-field
factor from the train adjacent cell applied to hold-out adjacent seats).

## 10. Seeds and budget

Master seed 20261001. Each script checkpoints to `data/` and runs in under 20 minutes and under 4 GB.

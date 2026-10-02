# Pre-registration A1: second out-of-sample test of the distance kernel on newly collected seat-/position-level outbreaks

[...] Written 2026-10-01 between 09:20 and 09:50 BST, after the
new outbreak records had been collected and transcribed and before any transmission model was evaluated on
any of them. `PREREG_FREEZE.json` holds the SHA-256 of this file, of the transcribed seat maps
(`src/a1/digitized.py`) and the time of freezing. Later choices go to `POSTHOC_LOG.md`.

## 0. Status and limits

1. **Data-aware, model-blind.** The analyst collected and transcribed every count used below, so all
   outcomes, including those of the hold-out events, were seen before this document was written. The
   published near/far attack rates of several events were read. What is fixed here, before any model output
   exists, is the event list, the outcome and source definitions, the geometry, the strata, the models, the
   calibration/hold-out split, the metrics and the decision rule. This is weaker than a blind
   pre-registration and is stated as such in every summary.
2. **Split chosen by a rule that does not use outcomes**: publication era and pathogen (section 4).
3. **Power.** The power of the decision rule depends on the calibrated parameters. It is computed after
   calibration and before any hold-out outcome is scored, by a script that reads hold-out geometry, stratum
   sizes and event totals but not the stratum counts, and is written to `PREREG_ADDENDUM_power.md`.
4. The result of the first validation round (`../../bsc_validation`) is known: M1 failed; M2 beat the
   well-mixed model on one flight and one classroom, failed on a coach and (conservatively) a call centre; a
   constant near/far risk ratio could have passed the old rule. This round is designed around that criticism.

## 1. Question and hypotheses

Does the non-local airborne term (model M2: quanta emitted by the infectious
person random-walk on the lattice and are removed at rate kappa), with one mixing coefficient per setting
class calibrated on one set of outbreaks, predict the distance gradient of infection risk in other
outbreaks, and does it do so better than (a) a well-mixed model and (b) a generic rule "risk near the index
is a constant multiple of risk elsewhere"?

- H1: M2 out-predicts the well-mixed model M0 on the hold-out events.
- H2: M2 out-predicts the constant-risk-ratio model CRR on the hold-out events.
- H3: every hold-out event lies inside M2's own predictive distribution.

A failure of any of these is reported as a failure.

## 2. Data: new events

None of these events was used in round 1 (the German meat plant was catalogued there as O4 but not
evaluated; Khanh et al. flight VN54 was a round-1 hold-out and is not reused). Sources, tables/figures and
access dates are in `sources.md`; raw copies in `data/raw/`.

### 2.1 Events that are modelled

| ID | Event | Pathogen | Setting, duration | Risk set and outcome (primary) | Position data | Grade |
|---|---|---|---|---|---|---|
| F1 | Flight CA112 Hong Kong-Beijing, 15 Mar 2003 (Olsen 2003) | SARS-CoV-1 | B737-300, 3 h | 111 other passengers; 18 with SARS (laboratory-confirmed or probable; includes those not interviewed) | full seat map | A |
| F2 | Chicago-Honolulu, May 1994 (Kenyon 1996) | M. tuberculosis | wide-body aircraft (type not in the abstract available), 8.75 h | 68 contacts in the index's cabin section without other risk factors; 6 TST-positive | two strata (within 2 rows 4/13; rest of section 2/55); no seat map | A |
| F3 | Los Angeles-Auckland, 25 Apr 2009 (Baker 2010) | influenza A(H1N1)pdm09 | B747-400 rear section, 13 h | 107 susceptible passengers; 3 laboratory-confirmed post-flight cases | full seat map; 9 sources | A |
| F4 | Cancun-Birmingham, Apr 2009 (Young 2014) | influenza A(H1N1)pdm09 | B767, 9.5 h | passengers with data and known seat; 9 infected in flight (clinical definition) | full seat map; 6 sources | A |
| F5 | Domestic flight to Naha, 23 Mar 2020 (Toyokawa 2022) | SARS-CoV-2 | B737-800, 2 h | 140 other passengers on the transcribed map; 14 PCR-confirmed | full seat map; index 23C | A |
| F6 | Sydney-Perth, 19 Mar 2020 (Speake 2020) | SARS-CoV-2 | A330-200 economy mid cabin, 5 h | passengers of the mid cabin who were not primary cases; 11 secondary cases (8 flight-associated + 3 possibly) | full seat map; 4 sources | A |
| F7 | EK448 Dubai-Auckland, 29 Sep 2020 (Swadi 2021) | SARS-CoV-2 | B777-300ER, 18 h | 84 other passengers, all PCR-tested twice; 4 infected in flight (C, D, E, F) | seats of everyone within rows 23-30; the rest known only to sit elsewhere | A |
| F8 | Tel Aviv-Frankfurt, 9 Mar 2020 (Hoehl 2020) | SARS-CoV-2 | B737-900, 4.67 h | 71 interviewed passengers outside the tourist group; 2 likely transmissions | full seat map; 7 sources | A |
| W1 | Ward 8A, Prince of Wales Hospital, bedside assessment 6-7 Mar 2003 (Wong 2004) | SARS-CoV-1 | hospital ward, 40 min | 19 medical students; 7 cases | three strata by assigned bed (3/3, 4/8, 0/8) and floor plan | A |
| W2 | Same ward, inpatients 4-12 Mar 2003 (Yu 2005) | SARS-CoV-1 | hospital ward, days | 74 inpatients; 30 cases | three strata by bay (13/20, 11/21, 6/33) and floor plan | B (more than one generation; same index as W1) |
| P1 | Beef-processing line, May 2020 (Günther 2020) | SARS-CoV-2 | 32 m x 8.5 m line in a cooled hall, 3 shifts | 78 workers with fixed stations; 20 positive | cumulative counts by 1-m distance from the index (Appendix Table S2) | A |

Ten new grade-A events, one grade-B event.

### 2.2 Events catalogued but not modelled (reasons fixed now)

Murphy 2020 (Ireland flight; no identified source), Foxwell 2011 (survey response 42%, no denominators by
seat), Marsden 2003 and Desenclos 2004 (only the secondary tabulation of Hertzberg 2016 available; 2 cases on
AF171), Chen 2020, Choi 2020, Bae 2020, Eldin 2020 (one or two cases, or pre-flight exposure dominant),
Cho 2016 MERS emergency room (zone-level, index moved between zones; belongs to the zone-level avenue),
Moser 1979 and Riley 1978 (no positions), Katelaris 2021 (already H2 in round 1), Hertzberg 2016 (review; used
for cross-checks of transcribed counts only).

### 2.3 Transcription rules (fixed)

- Seat maps are transcribed by eye from the published figures; every map is checked against the totals
  printed in the source (`00_check_data.py`), and each discrepancy is listed in `sources.md`.
- F1: all other passengers are in the risk set, as in Olsen's own risk ratio (8/23 vs 10/88).
- F3: the three laboratory-confirmed post-flight cases are the outcome (the source counts two of them as
  in-flight infections and one as possible; the seat map does not say which symbol is which).
- F4: passengers without data ("X") and the immune passenger are outside the risk set, as in the source.
- F5: flight attendants are not on the map and are excluded. Probable cases count as non-cases in the
  primary analysis.
- F6: sources are the four infectious Ruby Princess passengers of the mid cabin (the reference set of the
  source's own "within 2 rows" statement). The two infectious passengers with B.1 lineage and the
  non-infectious primary cases are neither sources nor susceptibles.
- F7: sources are passengers A (26G) and B (26D). Passenger G (infected in quarantine according to the
  source) is a susceptible non-case. The four family members in row 24 D-G occupy those four seats.
- F8: members of the tourist group who tested negative at the airport are outside the risk set (exposed for
  seven days before the flight); the seven not interviewed are outside the risk set.

## 3. Models

All models give each susceptible j an infection probability p_j = 1 - exp(-eps X_j), with one free
infectiousness parameter eps per event. With several sources, X_j is the sum over sources (equal emission).

- **M0, well mixed**: X_j constant.
- **CRR, constant risk ratio**: p_j = rho q for susceptibles in the near zone, q elsewhere (q fixed by the
  event total), with one constant rho per setting class calibrated on the calibration events. Near zone:
  aircraft, within two rows of a source (rows |d| <= 2); rooms, the index's own sub-space (W1: the index's
  cubicle; W2: the index's bay; P1: within 8 m, the zone singled out by the source, a definition that favours
  CRR). Reported in addition with the published constant rho = 2.4 (Hertzberg and Weiss 2016).
- **M2, framework with the non-local term**: X_j = Integral_0^T (T - s) exp(-kappa s) G(r_j, s | r0) ds
  on a rectangular lattice with reflecting walls, G the continuous-time lattice-walk propagator with
  coefficient D_air (`src/a1/lattice.py`, the verified round-1 implementation, copied unchanged). One D_air per
  setting class (aircraft cabins; large rooms). kappa is not fitted (section 3.1).
- **M1cf, same-site model** (people walk, same-site transmission), closed-form
  co-occupation kernel X_j = Integral_0^T G(r_j, 2t | r0) dt with one D_p per class. Reported for
  completeness; round 1 showed with the exact simulator that the closed form overstates M1's reach, so a
  pass of M1cf would not rehabilitate M1, and a failure confirms round 1.
- **M3, exploratory two-scale extension (not model M1, not part of the decision rule for M2)**:
  X_j = (1 - w) X_j^{M2} + w <X^{M2}>, <.> the lattice average (the M0 exposure); w in
  {0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95, 1} and D_air calibrated jointly. It stands for cabin
  recirculation and people moving about. It is included because round 1 showed infections far from the index
  that a single diffusion length cannot produce; it is evaluated by the same rule and reported separately.

### 3.1 Fixed constants and priors (none fitted)

- Loss rate kappa = ACH + 0.93 /h (deposition 0.3 + inactivation 0.63, as in round 1, used for every
  pathogen). Aircraft cabins: ACH log-uniform 10-30 (round-1 prior). Ward 8A: 7.79 ACH measured (Wong 2004),
  lognormal SD 0.2. Meat plant: ACH uniform 0.3-1.0 ("< 1", Günther 2020). Seven kappa nodes per event.
- Flat prior on log10 D over [0.01, 10^4] m^2/h, 61 grid points, for D_air and D_p.
- Aircraft lattice: one site per seat plus one site per aisle; a_x = 0.5 m, a_y = 0.8 m (row pitch), the
  whole mapped cabin (or cabin section) as one rectangle; rows missing on the map (galleys, exits) are kept
  as empty lattice rows. F2: generic 3-4-3 section of 14 rows, index seat uniform over rows 3-12, 13 near
  and 55 far susceptibles on random seats (marginalised over 200 configurations). F7: 3-4-3, lattice rows
  17-50; the 71 passengers not shown sit on random seats of rows 17-22 and 31-50 (200 configurations).
  F6: the mid cabin is one closed rectangle (17 rows, 2-4-2).
- Ward 8A: 21 m x 16 m open rectangle, a = 1 m; bed coordinates from Wong 2004 Fig. 4 (cubicle 9.2 m x 6.1 m,
  five beds per row); partitions ignored (the verified kernel has uniform D). W1: students stand at the
  assigned bed; the three "adjacent" students at bed 12 (as marked on the plan); the eight same-cubicle and
  eight other-cubicle students on random beds of the respective bed lists (200 configurations); T = 40 min.
  W2: subjects spread evenly over the beds of each bay; steady exposure (T = 100 h).
- P1: lattice 12 (lateral) x 32 (longitudinal) m, a = 1 m; index on the lateral edge (x in {0, 1}), longitudinal
  position uniform over the proximal half; workers on random sites of their distance bin (300
  configurations); three 8-h shifts.

## 4. Strata

- Aircraft with seat maps (F1, F3, F4, F5, F6, F8): row distance to the nearest source, four bins
  {0-1}, {2}, {3-5}, {6+}. F7: three bins {0-1}, {2}, {3+}. F2: two bins {0-2}, {3+}.
- W1, W2: the three published strata. P1: distance (0,4], (4,8], (8,12], >12 m.

The observable of every event is the vector of case counts over its bins, conditional on the event total K.

## 5. Calibration / hold-out split

Rule: events caused by pathogens other than SARS-CoV-2 (all published before 2020) calibrate; all
SARS-CoV-2 events are held out. The split is independent of the round-1 calibration (trains and a
restaurant, both SARS-CoV-2), which is not used in the primary analysis.

- Calibration, aircraft class: F1, F2, F3, F4. Calibration, room class: W1, W2.
- Hold-out: F5, F6, F7, F8 (aircraft kernel) and P1 (room kernel).
- Estimation: for each class and model, the likelihood of a parameter value is the product over calibration
  events of the conditional probability of the observed bin vector given K (eps profiled so that the expected
  total equals K; kappa prior and seat configurations averaged inside each event). Posterior weights on the
  grid are proportional to that likelihood. Predictions average over the posterior (400 draws).

## 6. Metrics

1. S_e(model) = log predictive probability of the observed bin vector of event e given K_e.
2. Delta0 = sum over hold-out events of [S(M2) - S(M0)]; DeltaC = sum of [S(M2) - S(CRR)].
3. Monte Carlo reference distributions (exact enumeration where feasible, otherwise 200,000 draws): Delta0
   under M0; DeltaC under CRR (posterior predictive) and under M2.
4. Adequacy: two-sided predictive p-value of the observed bin vector under each model (sum of the
   probabilities of all vectors no more probable than the observed one).
5. Descriptive: observed and predicted near (rows <= 2 / near zone) vs far risk ratios per event, pooled by
   inverse-variance weights on the log scale, with Cochran's Q, by pathogen and by setting.

## 7. Decision rule (for M2; applied identically to M3 and M1cf in separate columns)

- **B1**: Delta0 > 0 with one-sided p < 0.05 under M0.
- **B2**: DeltaC > 0 with one-sided p < 0.05 under CRR. Mirror statement: DeltaC < 0 with one-sided
  p < 0.05 under M2 means the constant-ratio rule predicts better than M2.
- **B3**: adequacy p >= 0.05 for every hold-out event.

| Outcome | Condition | Sentence licensed |
|---|---|---|
| V1 | B1, B2, B3, power of (B1 and B2) at least 0.8 | "the calibrated kernel predicts held-out distance gradients better than both a well-mixed model and a constant near/far risk ratio" |
| V2 | B1 and B3 (at most one failure); B2 not met and mirror not met | "better than well mixed; not distinguishable from a constant near/far risk ratio" |
| V3 | B1; mirror of B2 met | "better than well mixed, but a constant near/far risk ratio predicts better: the kernel shape is not supported" |
| V4 | B1 not met | "no out-of-sample support for the kernel" |
| V5 | two or more B3 failures | list-type sentence naming the events reproduced and those not, whatever B1/B2 give |

The aircraft class (F5-F8) and the room class (P1) are also reported separately; the pooled rule decides.
Levels (attack rates) are not tested: eps is fitted per event.

## 8. Secondary analyses (fixed list; reported, never used to upgrade the primary outcome)

- S-T4 (frozen transfer): M2 with the round-1 train posterior for D_air (no parameter fitted on any new
  event) predicts all eight flights; M2 with the round-1 restaurant posterior predicts W1, W2, P1.
- S-swap: calibrate on the SARS-CoV-2 events, predict the earlier ones.
- S-LOO: leave-one-event-out within each class; per-event preferred D and likelihood-ratio test of a
  common D; results by pathogen and by setting.
- S-seat: seat-level conditional log score instead of bins (seat-map events).
- S-out: wider outcomes (F5 confirmed + probable; F6 flight-associated only, 8; F3 with the suspected
  post-flight case; F7 with passenger G).
- S-src: F6 with the three symptomatic sources only and with all six infectious passengers; F3 with the
  suspected in-flight cases as sources; F7 with passenger A only.
- S-hh: F5 with each household (14B/F/G/H; 25F/G) collapsed to one unit at its first-listed seat.
- S-pitch: a_y = 0.75 and 0.86 m. S-kappa: cabin ACH 5-15. S-noW2: room class calibrated on W1 only.
- S-int: F1 restricted to interviewed passengers and cases.

## 9. Seeds, budget

Master seed 20261002. Every script checkpoints to `out/`, runs under 20 minutes and 4 GB.

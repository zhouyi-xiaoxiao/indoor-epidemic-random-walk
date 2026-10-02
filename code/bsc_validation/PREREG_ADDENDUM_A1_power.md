# Addendum A1 to the pre-registration: power of the decision rule, and late specifications

Written 2026-10-01 05:17 BST, after calibration on T4 and R1 (`data/01_cal_*_primary.json`) and after the
conditional predictive distributions of the hold-out events had been computed (`data/02_predictive.json`),
before any script had read `data/heldout_outcomes.json`. Numbers are from `data/02_power.json`
(exact enumeration over the joint outcome space of the four single-event shape tests, event totals fixed).

Note on the header of `PREREGISTRATION.md`: it gives a nominal writing window of 05:00-05:40; the file was in
fact frozen at 05:12:03 (`PREREG_FREEZE.json`), before the first calibration run (05:13:30, `logs/01_primary.log`).

## A1.1 Power of A1 and A2 (T1, T2, T5, C1; O1 not included)

Rows: which model is being judged. Columns: which model generated the data ("truth"), at its calibrated
posterior.

| Judged model | Quantity | truth = M0 (well mixed) | truth = M1 | truth = M2 |
|---|---|---|---|---|
| M1 | P(A1 passes) | 0.00004 | 1.000 | 0.903 |
| M1 | P(A2: no failure) | 0.0000 | 0.906 | 0.043 |
| M1 | P(A2: at most one failure) | 0.008 | 0.997 | 0.778 |
| M1 | P(A1 and no A2 failure) | 0.0000 | 0.906 | 0.043 |
| M2 | P(A1 passes) | 0.002 | 1.000 | 0.997 |
| M2 | P(A2: no failure) | 0.0003 | 0.069 | 0.874 |
| M2 | P(A2: at most one failure) | 0.022 | 0.758 | 0.994 |
| M2 | P(A1 and no A2 failure) | 0.0002 | 0.069 | 0.874 |

Reading:

- If a calibrated spatial model is true, it passes A1 with probability at least 0.997 and passes A1 and A2
  together (no failure) with probability 0.87 (M2) or 0.91 (M1). The pre-registered requirement "power of A1
  at least 0.8" for sentence W1 is met for both models.
- If the outbreaks are in fact well mixed, a spatial model passes A1 with probability at most 0.002 and
  passes A2 with at most one failure with probability 0.008 (M1) or 0.022 (M2). The amended rule cannot be
  passed by a model whose gradients are not in the data. This answers verifier objection V4.
- The rule also separates the two spatial models from each other: if M2 is true, M1 passes A2 without failure
  with probability 0.04; if M1 is true, M2 does so with probability 0.07.
- Per event, probability that the event rejects the well-mixed null at the 5% level if the judged model is
  true: T1 1.00 (M1) / 0.61 (M2); T2 0.97 / 0.81; T5 0.29 / 0.76; C1 0.83 / 0.96. Every event carries
  information under at least one model.

## A1.2 What the calibrated models predict (stated before unblinding)

Expected near/far risk ratio, median and 90% range over parameter and seating uncertainty:

| Event | M0 | M1 | M2 |
|---|---|---|---|
| T1 bus, rows 5-11 vs rest | 1 | 21 (11-36) | 2.3 (1.6-3.9) |
| T2 coach, rear vs front | 1 | 35 (15-84) | 11 (4.9-30) |
| T5 flight, near vs far | 1 | 1.8 (1.4-2.6) | 3.4 (2.2-5.9) |
| C1 classroom, front vs back | 1 | 5.6 (1.1-6.9) | 6.0 (2.8-7.0) |

## A1.3 Power of the level dispersion test (section 7.2)

If the Wells-Riley scaling is true and infectiousness varies between index cases with log-SD 0.94, the
permutation test over H1, O2, R2, C5 reaches p < 0.05 with probability 0.29 (the smallest attainable p-value
with four events is 1/24 = 0.042). The test is therefore reported but cannot support a level claim; a
non-significant result is uninformative. The plain comparison "SD under Wells-Riley scaling smaller than
under a constant hazard per hour" would hold with probability 0.85 if Wells-Riley is true.

## A1.4 O1 join-count test

The re-check of the outbreak dataset computed a power of at least 0.94 for the join-count
statistic against neighbour-driven spread with local weight at least equal to background. The predictive
distributions of the statistic under M0, M1 and M2 are computed in `scripts/06_callcentre.py`; the observed
value is compared only after they are stored.

## A1.5 Late specifications (omitted from version 1.0, fixed now, before the corresponding scripts exist)

- C3 (Utah classrooms): room volume log-uniform 180-350 m^3; exposure uniform 6.5-13 h (one to two school
  days); ventilation prior as listed in 3.4; 51 index cases, 728 tested contacts, about 14 tested contacts per
  index.
- S1 (supermarket customers): scene default floor 50 m x 30 m, ceiling 4 m assumed (6,000 m^3); dwell 27 min;
  one infectious employee present at steady state.
- T7 (metro ride): metro scene 22.5 m x 3 m, height 2.2 m (148.5 m^3); 20 min; 310 riders; criterion:
  probability of at least 5 secondary cases from one ride with one infectious rider must not exceed 0.5.
- T5 premium economy: exposure averaged over the 30 lattice sites of the five premium-economy rows, 35
  passengers; evaluated at the eps profiled from the business-class total.
- T4-anchored infectiousness uses the M2 calibration (eps of each posterior draw times the cell volume divided
  by the breathing rate).
- Tied-index check: posterior of eps given T2 (7/46) with a flat prior on log eps; T3 has 17 exposed.

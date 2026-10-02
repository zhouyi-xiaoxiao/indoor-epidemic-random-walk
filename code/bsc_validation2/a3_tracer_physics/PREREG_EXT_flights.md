# Pre-registration, extension EXT-F: tracer-measured cabin kernel on four further SARS-CoV-2 flights

Written 2026-10-01 (BST), after the main A3 analysis was unblinded and before any model of this directory was
evaluated on these four flights. Hash and time in `PREREG_EXT_FREEZE.json`.

## 0. Status: exploratory, outcome-aware, does not enter the A3 decision

- The main pre-registration fixed "no further outbreak is used" and excluded two of these flights by name
  (Speake 2020: several sources; Swadi 2021: no seat-level denominators). Avenue 1 of this round
  (`../a1_more_outbreaks`) has since built seat maps for them with its own conventions. Using them here is a
  post-hoc extension. It cannot change the A3 outcome (S4) and is reported separately.
- I have read avenue 1's interim summary: the near/far counts of the four flights (7/23 vs 7/117; 8/46 vs
  3/52; 4/13 vs 0/71; 2/20 vs 0/51), and that its outbreak-calibrated M2 (D = 126 m2/h, 90% 52-438, cabin air
  change 10-30 /h) beat M0 by +10.2 and lost to its calibrated constant ratio by 0.77 on these flights, with
  one adequacy failure (Naha, p = 0.047). The physics D of A3 (112 and 708 m2/h) brackets that value, so the
  result is largely foreseeable. What is new is only that no outbreak enters the kernel.
- Avenue 1 is unverified; its seat maps were transcribed by eye from low-resolution figures (its
  `POSTHOC_LOG.md` P3 lists denominator discrepancies). Its code is imported read-only; hashes of the files
  used are frozen with this document.

## 1. Events, geometry, statistic (all taken from avenue 1 unchanged)

F5 Naha (Toyokawa 2022), F6 Sydney-Perth (Speake 2020), F7 EK448 (Swadi 2021), F8 Tel Aviv-Frankfurt (Hoehl
2020): `a1.events.build(ev)` gives lattice (0.5 m x 0.8 m), source seats, susceptible seats, row-distance bins
(F5, F6, F8: rows 0-1, 2, 3-5, >5; F7: 3 bins), durations and totals K. Statistic: the vector of case counts
per bin given K (conditional distribution over all compositions, `a1.engine.cond_pmf`); near = bins within
two rows.

## 2. Models

- **P (physics kernel, primary)**: M2 exposure (`a1.engine.x_m2`), D_air from the Kinahan forward-cabin
  fits already frozen in `results/01_kernels.json` (707.9 and 112.2 m2/h, equal weights), removal rate
  = air change rate uniform 30-36 /h (7 quantile nodes) + 0.93 /h, exponential dose-response, level profiled.
  Random-configuration events (F7): avenue 1's 200 configurations, each paired with one node in turn.
- **P-econ (sensitivity, physics only)**: D fitted, with the estimator of the main pre-registration (log least
  squares, free constant, positive counts, release seat excluded, grid 0.1-1e5), to each unmasked breathing
  release of the economy sections of Kinahan et al. (777 FWD-MID, MID-AFT, AFT; 767 FWD, FWD-MID, AFT) on a
  0.5 m x 0.8 m lattice with seat columns A B C | D E F G | J K L (3-4-3), A B | D E F G | K L (2-4-2), A B |
  D E F | K L (2-3-2), ten padding rows fore and aft, kappa 35 /h (777) or 32 /h (767); equal-weight mixture
  over all releases. Removal rate as for P.
- **M0**: equal risk. **CRR(2), CRR(3)**: risk ratio 2 or 3 for seats within two rows.
- For reference only (copied from avenue 1, not recomputed): its calibrated M2 and calibrated CRR.

## 3. Metrics and rule

- Adequacy per flight: two-sided predictive p of the observed bin vector (sum of probabilities of vectors no
  more probable than the observed one), 5% level.
- Log score per flight and pooled over the four flights: P - M0, P - CRR(2), P - CRR(3); one-sided p of the
  pooled difference under the competitor by Monte Carlo (200,000 joint draws, seed 20261001 + 40).
- Rule (same logic as the main rule, four events): E1 = no adequacy failure, beats M0 (difference > 0, p < 0.05)
  and both CRRs (> 0, p < 0.05); E2 = no failure, beats M0, not both CRRs; E3 = no failure, does not beat M0;
  E4 = one failure; E5 = two or more failures.
- Power, computed by `src/e01_ext_predict.py` from predictive tables and K only, before the observed vectors
  are read: probability of E1-E5 when the truth is P, M0, CRR(2), CRR(3).
- Also reported: near/far risk ratio predicted and observed; the comparison of the tracer D with avenue 1's
  outbreak-calibrated D (descriptive).

## 4. Forecast

P beats M0 clearly (the far strata are almost empty of cases), does not beat CRR(2)/CRR(3), and Naha is at the
edge of adequacy. Expected outcome E2 or E4.

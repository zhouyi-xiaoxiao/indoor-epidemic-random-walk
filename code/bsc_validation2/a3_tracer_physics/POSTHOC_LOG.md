# Post-hoc log (A3 tracer physics)

Every choice made after `PREREG_FREEZE.json` was written, with time (BST, 2026-10-01) and reason.
None of these is used to upgrade the primary conclusion.

- **09:22 P1 (data handling).** Ou et al. Fig. S3B/S4B (minibus): automatic pairing of the EMF text records
  failed where a CFD-only value of an A seat is printed directly above the measured/CFD pair of the B seat, so
  the 47 records were assigned by hand from the dump (`data/raw/ou2022_supp/emf_text_dump.txt`) and checked
  against the rendered page. The coach (Fig. S3A/S4A) was assigned automatically without conflict; 45
  susceptible seats carry a Fig. S4A value (47 seats minus the empty 8C and the index seat 12D). Fig. S4B
  carries 15 passenger values (two at 4A, mother and child), so the seat-level minibus analysis uses 15
  susceptibles, not the 17 (or 12) of the first validation's level check. T3 is secondary.
- **09:33 P2 (forced by the data).** The pre-registration lists a release at seat 7A in the forward cabin of
  the 777. There is none: 7A is a sensor seat (the sheet lists sensors 5A, 5G, 5L, 7A after row 10). The
  unmasked breathing releases in that cabin are 5A, 5G and 5L. K-A for the flight therefore mixes two measured
  kernels (5L as is, 5A mirrored), not three, and D_cabin is the geometric mean of two fits. The centre-seat
  release 5G is added as a post-hoc sensitivity variant (S-5G: release 5G treated as if it were the index seat,
  mapped G -> K column by a lateral shift is not possible, so it is used unshifted with the index in 5G; shown
  for information only).
- **09:33 P3 (digitisation detail).** Overlapping markers of Woodward Fig. 9 merge into one blob; each
  detected blob (37 for the middle release, 23 for the end release, after excluding points within 0.02 L of
  the source) has equal weight. Lateral position: red markers (y < 0) on lattice column 1, blue on column 4;
  the source is put on column 1. The plot frame was located at pixel rows 18/278 and 356/616 and columns
  40/700.5.
- **09:33 P1.** kappa_tracer for the coach point estimate is 6.5 /h; bootstrap draws use uniform 5.79-7.25.
- **09:45 unblinding.** `src/s04_evaluate.py` is the first script that reads stratum outcomes. Hashes of the
  pre-registration, the power addendum, the kernels and the predictions immediately before it was written
  are in `PREREG_FREEZE_2_before_unblinding.json`. A formatting bug (a `None` risk ratio for the post-hoc 5G
  variant) stopped the first run after the results file had been written; the print statement was fixed and
  the script re-run with identical numbers.
- **09:50 P4 (diagnostic, post hoc).** The flight passes the primary test (p = 0.14) only through the
  mixture of two measured kernels that disagree. Looking at why: the mirrored 5A release has integrated
  counts of exactly zero in all three replicates at three sensors that map onto far seats of VN54 (1G, 3A,
  4A), and those zeros give these seats zero risk. With the same zeros treated as missing (pre-registered
  variant S-zero) the flight is rejected (p = 0.003). The primary result is kept as registered; the memo
  reports the flight as a fragile pass and recommends the conservative reading.
- **09:50 P5 (post hoc, descriptive).** Added: pooled score without the flight, without the classroom and
  without the restaurant (`verify/v01_checks.py`); the near/far dose ratio each outbreak would need under
  the exponential law; implied emission for the K-B events; adequacy failures of M0, CRR(2), CRR(3) under
  the same rule. None enters the decision.
- **09:50 P6.** T4 uses 20 quantiles of the bootstrap distribution of D_train times 7 ventilation nodes (140
  kernels) instead of all 200 bootstrap draws; O1 uses 60 parameter draws for P and 20 for M0, 1,500
  simulated outbreaks each. Dose-response variants were run for the five stratum events only (in T4 the
  attack rates are below 2%, where the three laws coincide).
- **09:50 P7 (limit of the pre-registered rule).** T4 fails the deviance test for every model considered
  (M0, CRR(2), CRR(3), the outbreak-calibrated kernel and the physics kernel). The pre-registration counted
  it as a D1 failure of P without asking whether any competitor passes; the memo says so.
- **10:21 resumed stretch of work.** `memo.md` did not exist (that stretch of work stopped after writing `src/s07_memo.py`);
  generated it from the stored results. Corrected three hard-coded sentences of the memo script against the
  result files: the flight is rejected in four of the six pre-registered sensitivity variants (not "three",
  not "5 of 10"); the first validation's classroom emission is quoted as the re-check's median 1.4e7 quanta/h;
  the flight emission of about 1,060 quanta/h is called high, not plausible. No number changed.
- **10:22 P8 (verification).** `verify/v02_repro.py` re-ran s01-s05 into `verify/repro/` (frozen result files
  untouched): zero numerical differences. Raw records of Ou Fig. S3A and of the Kinahan 5A/5L releases
  re-read; one-dimensional analytic fit of the carriage profile added as a cross-check.
- **10:24 P9 (post-hoc extension EXT-F).** Four further SARS-CoV-2 flights from the second outbreak test, own protocol
  `PREREG_EXT_flights.md` (hash in `PREREG_EXT_FREEZE.json`), written after the main unblinding and with
  knowledge of the near/far counts and of the second outbreak test's interim result. Contradicts the main pre-registration's
  "no further outbreak" and its exclusion of Speake 2020 and Swadi 2021, so it is reported separately and does
  not enter the decision. Tables and power frozen 10:26 (`PREREG_EXT_FREEZE_2_before_unblinding.json`) before
  `src/e02_ext_evaluate.py` read the observed vectors.
- **10:27 P10.** The first version of the extension self-check integrated the matrix exponential by Simpson's
  rule on 400 steps and disagreed with the mode sum by up to 37%; the quadrature was too coarse for the fast
  modes. Replaced by the closed form B^-2 (exp(BT) - I - BT) v; agreement 3e-12. The analysis code was not
  changed. Single-coefficient variants P_D708 and P_D112 were computed in e01 for information (not listed in
  the extension protocol).

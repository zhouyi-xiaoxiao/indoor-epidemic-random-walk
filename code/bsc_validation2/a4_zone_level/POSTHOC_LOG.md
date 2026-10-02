# Post-hoc log — the zone-level analysis (zone level)

Choices made **before** the freeze are part of PREREGISTRATION.md. Design-stage history, for transparency:

* D1 (before freeze, synthetic data only): the external-hazard term was first a 7-day centred incidence proxy from other
  schools. On a synthetic epidemic this biased the free-fit shares towards the class level and made the F-vs-Z likelihood
  ratio non-χ² (the proxy absorbs diffuse school-level transmission). It was replaced by the causal city-level term
  `α0 + α1 P_city\s(t)` so that the fitted model and the simulator coincide. No real outcome was involved.
* D2 (before freeze): h5py was installed into `vendor/` (inside this directory) to read the JLD2 file; Julia is not available.

Entries below are made **after** the freeze.

* P0 (2026-10-01, after freeze): `src/04_matsumoto.py` was run once, unchanged (hash equal to the frozen one), 08:49–08:56 UTC;
  output `results/matsumoto_primary.json`. An earlier stretch of work was interrupted after this step. No re-run, no change of seed.
* P1 (resumed stretch of work, same day): `src/05_seoul_floors.py` executed for the first time inside this analysis, unchanged (hash
  equal to the frozen one); output `results/seoul_floors.json`, log `logs/05_seoul.log`. The same computation had
  been made separately before, for the overview of the four analyses (not distributed); the numbers agree (median 4, 95 % interval 0–19).
* P2 (post hoc, descriptive, non-decisional): `src/06_posthoc_descriptive.py` → `results/posthoc_descriptive.json`:
  (a) recount of students/cases/classes and recomputation of the heterogeneity statistic directly from the CSV with pandas
  (agrees: 3.383); (b) number of schools in which Z is ahead of H, X, C and leave-one-school-out range of the totals;
  (c) own-class/rest-of-school per-pair ratio implied by the Lyon shares in the Matsumoto class structure, for comparison
  with the ratio fitted by model C; (d) AIC of the full-data fits; (e) Seoul 11th-floor wings: coupling that would reproduce
  79/137 vs 4/62 — fitted quantity, no test (as registered in section 3.3). First version of (a) counted a school with zero
  cases in the degrees of freedom (3.343); corrected to the registered definition, which excludes it.
* P3: figures `src/07_figures.py` (reads stored results only). Panel (d) shows the binomial 95 % interval of the well-mixed
  competitor (62–96), which is not in the registered script; display only.
* P4 (wording): the registered verdict rule gives "pass". Because (i) the bootstrap interval of TV (0.047–0.168) crosses the
  tolerance 0.15, (ii) one pre-declared sensitivity variant (no winter-break switch) gives TV 0.162, and (iii) the generic
  model C is not beaten out of sample, the memo reports the outcome as **pass_weak**, not as a robust pass. This is a
  downgrade relative to the mechanical application of the rule; no upgrade of any tier was made.
* Not done: no additional dataset was analysed. Candidates with zone denominators but without independent mixing
  measurements are listed in the memo (section 6); fitting a coupling to them would not be a test.

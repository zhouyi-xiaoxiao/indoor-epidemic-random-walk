# Pre-registration A3: non-local airborne kernel calibrated on physics measurements, then tested on outbreaks

[...] Written 2026-10-01 (BST) before any transmission
model was evaluated on any outbreak count in this directory, and before any mixing coefficient was estimated
from any tracer data set. `PREREG_FREEZE.json` holds the SHA-256 of this file and of the raw physics files at
the moment of freezing. Later choices go to `POSTHOC_LOG.md`.

## 0. Status and limits

1. **Data-aware on both sides.** The outbreak counts are those of the first validation
   (`../../bsc_validation`), which I have read, including which events the calibrated model failed (Hunan
   coach, Seoul call centre) and which D each event prefers when fitted alone (flight 0.1, trains 25, bus 63,
   coach about 6,300 m2/h). The tracer sources were opened before this file was written, to find out what
   they contain: I have seen the seat-level values of Ou et al. Fig. S3/S4, the table of Li et al., the summary
   statistics and one summary sheet of Kinahan et al., Fig. 9 of Woodward et al. and the regression of Cheng
   et al. No kernel has been fitted, no D has been estimated, and no predicted distribution has been computed.
2. What is fixed here: the estimator that turns each physics data set into a kernel, which kernel is used
   for which outbreak, the observation model, the competitors, the metrics, the decision rule and the
   sensitivity list.
3. Power is computed after the kernels are estimated (physics data only) and before the comparison with
   outbreak counts; it is written to `PREREG_ADDENDUM_power.md`. The power script reads strata sizes and event
   totals, not stratum counts.

## 1. Hypothesis

H: the non-local airborne term of model M2, with its spatial kernel fixed by tracer-gas or
tracer-aerosol measurements made independently of who was infected, and with one infectiousness parameter per
event, predicts the spatial pattern of infection in documented outbreaks in the matching enclosure: it is
adequate for every event, and it predicts better than (a) the well-mixed model M0 and (b) a generic model
with a constant near/far risk ratio.

Two forms of "physically measured kernel" are tested.

- **K-A (direct kernel).** Exposure of seat j is proportional to the tracer dose measured, or computed by a
  tracer-validated CFD model of the same enclosure, at seat j for a source at the index seat. No lattice
  walk. This tests the non-local term and the dose-response law with the kernel taken from measurement.
- **K-B (M2 with measured D).** The M2 model of the first validation (lattice-walk propagator [...]
  for the quanta, removal rate kappa, `../../bsc_validation/src/valmod/lattice.py`) with D_air
  estimated from tracer data of the same enclosure type. This tests the diffusion form of the kernel.

Neither form contains mobile people or heterogeneous D(r); everyone is seated in every event. A pass would
support the non-local term with a physical kernel, not the same-site rule (M1), which is not re-tested.

## 2. Outbreak data (unchanged from the first validation)

Counts, strata and sources as in `../../bsc_validation/data/heldout_outcomes.json`, `../../bsc_outbreaks`:

| Event | Strata (near vs far) | Outcome used | Kernel (primary) |
|---|---|---|---|
| T2 Hunan coach (Ou 2022, Luo 2020) | rear rows 8-13 vs front rows 1-7 | 3/19 vs 4/26 | K-A, same vehicle |
| T1 Zhejiang bus (Shen 2020) | rows 5-11 vs other rows | 14/33 vs 9/34 | K-B, D from the Hunan coach |
| T5 flight VN54 (Khanh 2020) | <=2 seats vs >2 seats from index, business class | 11/12 vs 1/8 | K-A, Boeing 777 forward cabin |
| T4 high-speed trains (Hu 2021) | 22 seat-offset cells, adjacent seat excluded | reconstructed counts, 138 cases | K-B, D from a rail carriage |
| R1 Guangzhou restaurant (Li 2021, Lu 2020) | tables TB, TC, T18 vs 14 remote tables | 2 to 5 of 16 vs 0 of 63 | K-A, same room |
| C1 Marin classroom (Lam-Hine 2021) | front two rows vs back three | 8/10 vs 4/14 | K-B, D from room relation |
| O1 Seoul call centre (Park 2020) | join-count z in the north wing | z = -0.64 | K-B, D from room relation |
| T3 Hunan minibus (Ou 2022) | seat level, 2 cases | secondary only | K-A, same vehicle |

New relative to the first validation: R1 and T4 were calibration sets there and are test sets here, because
nothing is calibrated on outbreaks; the VN54 seat map is digitised from Khanh et al. Fig. 1
(`data/raw/khanh2020_F1.jpg`): index 5K; occupied K2-K5, G1, G3-G7, D1, D3-D7, A1-A5; cases 3K, 4K, 4G-7G,
3D-7D, 5A; near = {3K, 4K, 3G, 4G-7G, 3D-7D} (12 seats, 11 cases), far = {2K, 1G, 1D, 1A-5A} (8 seats, 1 case,
2A lost to follow-up and counted as a non-case as in the source's table). The Hunan coach seat map is taken
from Ou et al. Fig. S4A (occupied seats are those with a CFD value) with cases at 1D, 5C, 6A, 6D, 9C, 9D, 13D.

Searched for and not used (fixed now): Speake 2020 (Sydney-Perth flight: eleven infectious sources of unknown
relative infectiousness in two cabins; the cabin contrast tests infectiousness, not the kernel); Olsen 2003
(SARS flight: seat map not openly accessible, only one printed stratum); Swadi 2021 (no denominators by
seat). No further outbreak with an independently measured kernel of the same enclosure was found.

## 3. Physics data and estimators (no outbreak count enters)

All raw files are under `data/raw/`; sources, access dates and hashes in `data/SOURCES.md`.

**P1 Hunan coach and minibus.** Ou et al. 2022, Build. Environ. 207:108414, supplementary Fig. S3 (ethane
tracer measured at 9 points on the coach and 14 on the minibus, and CFD at the same conditions) and Fig. S4
(tracer-validated CFD at all seats under the conditions of the infection trip, index at 12D / 4C). Values are
read from the text records of the EMF figure files.
- K-A(T2), K-A(T3): the Fig. S4 value of each occupied seat. Steady state is assumed (as the source does);
  the four passengers with shorter rides get no correction in the primary analysis.
- Measurement uncertainty: sigma_log = SD of log(measured/CFD) over the Fig. S3 points of that bus (reference
  seat excluded). Predictive distributions are mixtures over 300 kernels with independent lognormal
  seat-level noise of that SD.
- D_coach for K-B: least-squares fit of log steady-state M2 concentration to log measured tracer (Fig. S3A,
  measured values, source 12C, reference seat excluded, free additive constant), on the coach lattice of the
  first validation, kappa_tracer = air change rate of the test condition (heating and fans on: 5.79 and 7.25
  ACH, mean 6.5; no deposition or inactivation for a gas). D on a log10 grid from 0.1 to 1e5 m2/h, step 0.05.
  Uncertainty: bootstrap of the tracer points (200 resamples) combined with kappa_tracer uniform 5.79-7.25.
  If the optimum is at the upper grid edge the kernel is well mixed and is used as such.

**P2 Aircraft cabin.** Kinahan et al. 2021, PLoS ONE 16(12):e0246916, raw workbooks on figshare
(10.6084/m9.figshare.13537349, .13537319): in-flight releases of 1 um fluorescent tracer particles from a
breathing mannequin, time-integrated counts at 40+ breathing-zone sensors. Forward (1-2-1) cabin of the
777-200, rows 1-10, columns A, D, G, L, releases at 5A, 5G, 5L, 7A; measured air change rate 35 /h.
- Per release seat: mean over the unmasked breathing replicates of the integrated count at each sensor seat;
  blanks are missing, zeros are kept as zeros (as in the source's statistics).
- K-A(T5): the VN54 index sat in window seat 5K. Three measured kernels are used: release 5L as is
  (L = K column), release 5A mirrored (A<->L, D<->G), release 7A mirrored and shifted two rows forward. A VN54
  seat with no value takes the value of the mirror seat of the same row in the same map, else the mean of
  its row, else the mean of the adjacent rows. The primary predictive distribution is the equal-weight
  mixture over the three kernels (release-to-release variation is the measurement uncertainty).
- D_cabin for K-B: for each of the three releases, least-squares fit of log M2 steady-state concentration to
  log integrated count (positive values only, release seat excluded, free constant) on the VN54 cabin lattice
  of the first validation with kappa = 35 /h; D grid as in P1. D_cabin is the geometric mean of the three
  fits; the uncertainty distribution is the three values with equal weight.
- Removal rate for VN54: air change rate uniform 30-36 /h (the 32 and 35 /h reported for the two wide-body
  airframes tested) plus 0.93 /h. The first validation's prior (10-30 /h) is a sensitivity analysis.

**P3 Rail carriage.** Woodward et al. 2022, Indoor Air 32:e13121, Fig. 9: steady-state concentrations of
nebulised salt aerosol (PM2.5) at about 18 monitors along a UK inter-city carriage (saloon 19.85 m x 2.69 m x
2.02 m, design 11-15 ACH), release in the middle (B) and at the end (A) of the saloon.
- Points are digitised from the figure image by colour segmentation (marker centroids; axes calibrated on
  the plot frame). Points within 0.02 L of the release position are excluded (the source flags them as "near
  source"). The solid 1D-model curve is digitised as a cross-check.
- D_train: least-squares fit as in P1 on a 6 x 20 lattice (pitches 0.45 m x 0.99 m), reflecting ends, source
  at the release position, kappa_tracer uniform 11-15 /h; primary: release B (middle), sensitivity: release A.
- Applied to the Chinese second-class coach of T4 with the lattice and the ventilation prior of the first
  validation (5-20 ACH + 0.93). Transfer across rolling stock is an assumption.

**P4 Restaurant.** Li et al. 2021, Build. Environ. 196:107788, Table 3: time-averaged ethane concentration at
tables (measured at 8 tables, tracer-validated CFD at all), normalised to table A.
- K-A(R1): tracer concentration of the table (measured where available, otherwise CFD) times the overlap
  time with table A. Secondary kernel: the source's "normalised predicted exposure to droplet nuclei" (CFD
  with deposition, filtration and inactivation), which is a model output, not a measurement.

**P5 Rooms without their own tracer test (classroom, office).** Cheng et al. 2011, Environ. Sci. Technol.
45:4016: K/L^2 = 0.52 ACH + 0.31 per hour, L the cube root of the room volume (11 CO-tracer experiments in two
naturally ventilated rooms, R2 = 0.92, K from 0.001 to 0.013 m2/s).
- D_room = (0.52 ACH + 0.31) V^(2/3) m2/h with ACH drawn from the event's ventilation prior of the first
  validation and a lognormal factor with 90% of its mass within a factor 2 (sigma = 0.42) for the scatter
  of the regression and the transfer to mechanically ventilated rooms. V is the volume of the event lattice.

## 4. Models compared

For every model the infection probability of susceptible j is p_j = 1 - exp(-eps X_j), with eps profiled
so that the expected number of cases equals the event total K, and the test statistic is the number of cases
in the near stratum given K (conditional Poisson-binomial distribution, `valmod.stats.cond_pmf`).

- **P (physics kernel)**: X_j from K-A or K-B as assigned in section 2 (the "primary assignment"). The
  other form, where both exist (T2, T5), is a sensitivity analysis.
- **M0 (well mixed)**: X_j constant; hypergeometric.
- **CRR(rho)**: constant near/far risk ratio, p_near = rho q, p_far = q, with q such that the expected total
  is K; rho = 2 and rho = 3 (the verifier's benchmark of the first validation). CRR-best (rho maximising the
  pooled likelihood) is reported for information only.
- **M2-cal** (the calibrated kernel of the first validation, D_air 25 m2/h vehicles, 0.5 m2/h rooms, at the
  posterior mode) is shown for reference where available; it is not re-fitted.

Event-specific points:
- T4: 22 cells; conditional multinomial given the total 138. Exposure of a cell = mean over all seat pairs
  with that offset and over the travel-time distribution (as in the first validation). CRR near = same row.
- R1: outcome interval-censored. Scores use K = 2 (one case at TB, one at TC: the number the source asserts
  was acquired at the lunch); adequacy is evaluated at K = 2, 3, 4, 5.
- O1: epidemic percolation with the kernel between desks, eps tuned to a mean final size of 79, outbreaks of
  69-89 cases retained (script logic of the first validation), probability by pooling over 60 draws of
  (ACH, pitch, D-factor); no stratum likelihood, so O1 enters adequacy only.

## 5. Metrics

1. Adequacy per event: two-sided predictive p-value of the observed near count (sum of probabilities of
   outcomes no more probable than the observed one); T4: deviance of the 22 cells against the model's
   conditional probabilities, p by parametric bootstrap (2,000); O1: P(z <= z_obs), failing outside
   [0.025, 0.975]; R1: fails only if p < 0.05 for every K from 2 to 5.
2. Log score per event: log P(observed near count | K) under each model; pooled over the five stratum
   events T2, T1, T5, C1, R1(K=2). T4 is scored separately (it would dominate any pooled sum).
3. Seat-level conditional log-likelihood of the observed case set against M0 for T2, T3, T5 (secondary).
4. Descriptive only: predicted risk ratios; D_phys against the D each event prefers when fitted alone
   (first validation, Table V16); implied emission in quanta/h.

## 6. Decision rule

Primary assignment, exponential dose-response, seven events (T2, T1, T5, T4, R1, C1, O1).

- **D1 adequacy**: number of events failing metric 1 at the 5% level.
- **D2 beats well-mixed**: pooled log-score difference P minus M0 is positive and its exact one-sided p-value
  under M0 (enumeration over the joint distribution of the five stratum events) is below 0.05; and for T4 the
  log-score difference is positive.
- **D3 beats the generic gradient**: pooled difference P minus CRR(2) and P minus CRR(3) both positive, with
  exact one-sided p < 0.05 when CRR(rho) is taken as the truth.

Outcomes:
- **S1** no D1 failure, D2 and D3 hold: "with the kernel fixed by physical measurement the framework predicts
  the spatial patterns, and better than a generic gradient".
- **S2** no D1 failure, D2 holds, D3 does not: "compatible with all events and better than well mixed; not
  shown to be better than a generic gradient".
- **S3** no D1 failure, D2 fails: "compatible, not discriminating".
- **S4** exactly one D1 failure: list-type sentence naming it.
- **S5** two or more D1 failures: the hypothesis fails; no agreement claim; the failures are reported.

Power (addendum, before unblinding): probability of each outcome when the truth is P, M0, CRR(2), CRR(3).
For CRR truths O1 is generated under M0 and T4 with near = same row. Per event, the probability that M0 and
each CRR are rejected when P is true, and that P is rejected when each of them is true, to label events as
discriminating or not. An event at which P and M0 predict nearly the same thing is declared non-discriminating
and a pass there is not counted as evidence for the kernel.

## 7. Sensitivity analyses (fixed list)

- S-lin: linear dose-response p proportional to X (low-dose limit). S-gam: p = eps X / (1 + eps X)
  (exponentially distributed susceptibility).
- S-form: K-B in place of K-A for T2 and T5; K-A with measured points only for T2 (each seat takes the value
  of the nearest measured point of Fig. S3A).
- S-zero: Kinahan zeros treated as missing. S-vent: VN54 ventilation prior of the first validation.
  S-2A: VN54 far stratum 1/7.
- S-R1: K = 5; the source's droplet-nuclei exposure column.
- S-end: D_train from the end release. S-Dx2, S-D/2: D_room doubled and halved. S-cfdD: D_coach fitted
  to the CFD values of Fig. S4A instead of the measured points.
- S-side: T2 driver side vs opposite side (6/24 vs 1/20), secondary split of the first validation.

## 8. Forecast (written before any computation)

My expectation, recorded so that it can be checked: the coach kernel is close to well mixed, so T2 passes
and is non-discriminating against M0; the cabin kernel has a near/far dose ratio of only a few, so with 11 of
12 near passengers infected the flight is probably not reproduced; the carriage kernel is nearly flat over
three rows, so the steep same-row excess of T4 is not reproduced; the restaurant tracer contrast (remote
tables at 0.4-0.9 of the index table) is too weak for 0 of 63 at K = 5 and marginal at K = 2; the room
relation gives a mild gradient, so the call centre probably passes and the classroom is uncertain. If this
forecast is right the outcome is S5 for the patterns that need a steep gradient, with the physics kernel
repairing the two earlier failures (coach, call centre).

## 9. Computing

Python environment `../../env/venv`, numpy/scipy; seeds fixed (20261001 + offsets); every run under 20 min
and 4 GB; outputs under this directory only. Code of the first validation is imported read-only.

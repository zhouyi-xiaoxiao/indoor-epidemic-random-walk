# Avenue 2 — PRE-REGISTRATION: the movement-and-contact component against real human contact data

Written before any model (generative or epidemic) was run on any real dataset. Frozen by SHA-256 in
`PREREG_FREEZE.json` together with the analysis code (`src/a2lib.py`, `src/pipeline.py`, `src/sir.c`,
`src/02_run.py`). Everything decided later goes to `POSTHOC_LOG.md`.

## 0. What had been looked at before this freeze (disclosure)

* Raw files were downloaded and two purely descriptive scripts were run in an earlier, interrupted stretch of work:
  `src/00_describe.py` (N, number of days, raw contact counts per day, number of people present, first/last
  interval, group labels) and `src/00a_align_check.py` (time-offset check between the face-to-face files and the
  co-presence files). No contact-structure statistic of section 4 and no model output was computed on real data.
* A power analysis was run on SYNTHETIC data only (`src/01_power.py`, `out/01_power*.json`), section 7.
* Prior knowledge of the analyst (cannot be un-known): the SocioPatterns literature reports broad contact-duration
  and inter-contact distributions, strong group (class/department) structure and day-to-day repetition of
  contacts. The expectations in section 6 are therefore informed guesses, and are written down as such.

## 1. Question and hypothesis

The model's defining mechanism [...]:
individuals perform a Markov nearest-neighbour random walk on a square lattice with a possibly heterogeneous
hopping rate, harmonic-mean rate at interfaces, reflecting walls; transmission/contact = two individuals on the
same site. Question: **are the contacts produced by this mechanism statistically and epidemiologically
equivalent to real indoor human close-proximity contacts?**

H_RW (the claim under test): a lattice random-walk co-location model, calibrated on part of the days of a
deployment, reproduces the contact structure of the remaining days and the SIR outcomes on the real temporal
network of those days.

## 2. Data (real, open, citable)

SocioPatterns face-to-face proximity (RFID badges, detection range ≈1–1.5 m, 20 s resolution), downloaded from
http://www.sociopatterns.org/datasets/ (index page saved as `data/datasets_index.html`, files under `data/raw/`,
SHA-256 in `data/raw_SHA256.txt`, accessed 2026-10-01):

| key | setting | file | reference |
|---|---|---|---|
| InVS13 | office building, 2013 | workplace_InVS_tij.dat.zip + metadata | Génois et al., Network Science 3:326 (2015), doi:10.1017/nws.2015.10 |
| InVS15 | same office building, 2015 | workplace_InVS15_tij.dat.gz + metadata | Génois & Barrat, EPJ Data Sci. 7:11 (2018), doi:10.1140/epjds/s13688-018-0140-1 |
| LyonSchool | primary school, 2 days | primaryschool.csv.gz + metadata | Stehlé et al., PLoS ONE 6:e23176 (2011), doi:10.1371/journal.pone.0023176 |
| Thiers13 | high school, 5 days | HighSchool2013_proximity_net.csv.gz + metadata | Mastrandrea, Fournet & Barrat, PLoS ONE 10:e0136497 (2015), doi:10.1371/journal.pone.0136497 |
| LH10 | hospital ward, 5 days | hospital_lyon_contacts.dat.gz (roles inline) | Vanhems et al., PLoS ONE 8:e73970 (2013), doi:10.1371/journal.pone.0073970 |
| SFHH | scientific conference, 2 days | SFHH_tij.dat.gz | Génois & Barrat 2018 (above); Stehlé et al., BMC Med. 9:87 (2011) |

No number is typed from a paper: every quantity is computed from the raw files by `a2lib.load`.
Pre-processing (frozen in `a2lib.load`): time → 20 s intervals; clock offsets LH10 +46800 s (t=0 is Monday
13:00), Thiers13 +3600 s (unix UTC→CET); self-pairs dropped; duplicate records merged; calendar days carrying
<1 % of all records dropped (weekends). Retained days are indexed 0..D-1.

Not used: the Malawi village and Hypertext-2009 files present in `data/raw` (kept as reserve, not part of any
decision); the co-presence files (`colocation.tar.gz`; exploratory only if used, labelled as such);
Copenhagen Networks Study and ATC trajectory data (large downloads on a disk-constrained shared machine) —
this is a scope limitation and is to be stated in the memo.

**Roles of datasets.** DEV = {InVS13, LyonSchool, LH10}; CONF = {InVS15, Thiers13, SFHH}. CONF datasets are not
touched by any model or statistic until `FIX_FREEZE.json` exists (section 8).

**Train/test split within every dataset** (`a2lib.split`): first ceil(D/2) retained days = train, last
floor(D/2) days = test. (InVS13/15: 5+5; Thiers13, LH10: 3+2; LyonSchool, SFHH: 1+1.)

## 3. Models (exact variants)

All models receive the same exogenous input on each day: who is present and when (each individual's presence
window = first to last own contact of that day; people with no contact that day are absent). Time step 20 s.

* **RW0 — lattice walk, homogeneous.** Rectangular lattice of M sites (Lx≈Ly), reflecting walls. Each step each
  present walker attempts a move to a uniformly chosen von Neumann neighbour, accepted with probability p
  (discrete-time version of the continuous-time walk). Start site uniform at arrival each day.
  Contact = same site at the sampling instant. Two parameters: M, p.
* **RWhet — lattice walk, heterogeneous two-zone.** As RW0 with a slow zone (fraction f_slow of columns, hop
  probability p_slow) and a fast zone (p_fast); interface acceptance = harmonic mean
  2 p p'/(p+p') (harmonic-mean rule for W(r→r')). Four parameters: M, p_fast, p_slow, f_slow.
* **WM — well-mixed null.** Every jointly present pair has the same contact probability c per interval; contact
  events start as a Poisson process, durations geometric with the training mean. Two parameters: c, mean duration.
* **BLOCK — constant within/between-group ratio null** (the analogue here of the "constant near/far risk ratio"
  competitor of the first validation): as WM with c_in for same-group pairs and c_out otherwise, groups = the
  recorded class / department / role. Three parameters. (SFHH has no groups: BLOCK ≡ WM there.)
* **REALtrain — noise ceiling.** The statistics of the real training days scored against the test days by the
  same rule; shows which statistics are stable enough to be predicted by anything.

## 4. Calibration (train days only; `a2lib.calibrate`)

* c = (contact pair-intervals) / (pair-intervals of joint presence); M = 1/c (so RW0's stationary co-location
  probability equals c).
* mean contact-event duration m → RW0's p from the same-site persistence probability q = (1-p)² + p²/4 = 1-1/m.
* RWhet: two-component geometric mixture fitted by EM to training durations → (q_slow, q_fast, weight) →
  p_slow, p_fast by the same inversion, f_slow from the time share of the slow component.
* BLOCK: c_in, c_out from within/between-group contact and joint-presence totals.
* Nothing is fitted on test days. For the epidemic comparison only, all models' contact probability is rescaled
  by one common factor so that the expected total contact time on the test days equals the real one
  ("matched mean contact rate", as planned for this analysis).

## 5. Predicted quantities and pass criteria (frozen in `a2lib.FAMILIES`, `metric_pass`, `epi_pass`)

Statistics on test days (model value = mean over 5 generated realisations):

| family | metrics | tolerance |
|---|---|---|
| S0 level | total contact pair-intervals | ratio in [0.8, 1.25] |
| S1 duration | mean event duration; share of contact time in events ≥ 15 intervals (5 min) | ratio in [0.8,1.25] (fractions: or abs. diff ≤ 0.03) |
| S2 inter-contact | burstiness B=(σ-μ)/(σ+μ) of an individual's gaps between contacts (abs. diff ≤ 0.10); recurrence = 1 − distinct pairs/events (ratio or ≤ 0.03) |
| S3 degree | mean and CV of distinct contacts per person-day (person-days present ≥ 30 min) | ratio |
| S4 heterogeneity | CV of individual total contact time | ratio |
| S5 persistence | fraction of a test day's contact pairs already in contact the previous day | ratio or ≤ 0.03 |
| S6 groups | within-group share of contact time (n/a for SFHH) | ratio or ≤ 0.03 |
| S7 density (descriptive, outside the ladder) | exponent α of contacts per 20-min bin vs number present; RW/WM predict mass action α≈2. Compared if the 90/10 occupancy range ≥ 1.5; "agrees" if |α_model − α_obs| ≤ 0.3 |

A family passes iff all its applicable metrics pass.

**E — epidemiological consequence.** SIR on the temporal network of the test days repeated periodically; one
uniformly random index case at a uniformly random time; per-contact-interval transmission probability β,
exponential infectious period of mean τ. Four scenarios (nominal R0, τ) ∈ {1.5, 3} × {1 d, 4 d}, β set from the
real mean contact rate (`epi_beta`) and used identically for real and model networks. Outputs: R_index (mean
number infected directly by the index case), P_major (final attack ≥ 10 %), mean attack rate. Real network:
4000 runs; each model: 10 generated networks × 400 runs. A scenario passes iff R_index ratio ∈ [0.9, 1/0.9] and
|ΔP_major| ≤ 0.05 and |Δattack| ≤ 0.05. **E passes iff ≥ 3 of 4 scenarios pass.**

## 6. Decision rule and stated expectations

Per model, over a set of n datasets (need = ceil(2n/3)):

* **A1 (the movement component is validated):** on ≥ need datasets the model passes E and ≥ 6 of the 7
  S-families (≥ 5 of 6 where S6 is n/a).
* **A2 (epidemiologically adequate and better than well-mixed):** not A1, but on ≥ need datasets the model
  passes E while WM fails E.
* **A0 (fail):** neither. If additionally WM or BLOCK does at least as well on E as the walk (number of
  datasets passing E), this is reported as "the lattice walk adds nothing over a non-spatial contact model".

Applied (i) to RW0 and RWhet over all six datasets — the test of the lattice-walk mechanism as formulated — and
(ii) to the minimal modification(s) over the three CONF datasets (section 8).
Per-family pass counts across datasets are reported for every model whatever the ladder says; a family is
called "reproduced" for a model if it passes on ≥ need datasets, "not reproduced" otherwise.

Expectations, written before any run (informed by the literature, see section 0): RW0/RWhet will pass S0 and
probably mean duration, and will FAIL S1's long-event share (geometric tail), S5 (no memory of who met whom),
S6 (no groups) and S2 recurrence; E is expected to fail in the direction "walk over-predicts outbreak
probability and size" because real contacts are concentrated on few repeated partners. I.e. the analyst expects
A0 for the lattice walk as formulated. Whatever happens is reported.

## 7. Power analysis (synthetic data only; `src/01_power.py`, `out/01_power_table.json`)

Synthetic "truths" with N=120, 4 groups, 4 days were generated from WM, BLOCK, RW0 and SYN (a structured
out-of-family truth: fixed log-normal pair affinities × group factor, Pareto durations), 4 replicates each, and
the frozen pipeline scored every candidate (2000/200 SIR runs; the real-data runs use 4000/400, i.e. less noise).

* Correct (or nesting) model on its own truth: all 7 S-families and E passed in 4/4 replicates for WM|WM,
  BLOCK|BLOCK, RW0|RW0, and for the nesting cases BLOCK-on-WM and RWhet-on-RW0 → the rule does not reject a true
  model through Monte-Carlo noise (0/20 false rejections; one-sided 95 % bound 14 %).
* Wrong, non-nested model reaching the per-dataset A1 criterion (E and ≥ 6 S-families): 0 of 40 model×replicate
  cases (max S-families passed by a wrong model: 4 of 7). One-sided 95 % upper bound on the per-dataset
  false-pass probability: 7 %; A1 additionally needs this on ≥ 2/3 of the datasets, so with independent datasets
  P(wrong model gets A1 | n=3) ≤ 3·0.07² ≈ 1.5 %, and far less for n=6.
  In particular the well-mixed null on a random-walk truth fails S2, S3, S4 (3/4), S5 and E (2/4); the
  constant-ratio BLOCK null on a random-walk truth fails S2, S3, S5.
* E ALONE is weakly discriminating: on a BLOCK truth, WM passed E 4/4 and RW0 3/4. Therefore E alone is never
  taken as support for the mechanism (hence A1 requires the S-families, and A2 requires WM to fail E).
* A2 can be reached by a wrong model: on the SYN truth RWhet passed E 4/4 while WM failed 4/4, although RWhet
  passed ≤ 3 S-families. A2 is therefore worded only as a comparative, epidemic-level statement, never as
  "the mechanism is validated".

## 8. Failure diagnosis and minimal modification (second freeze)

1. Stage 1: base models (WM, BLOCK, RW0, RWhet) on DEV only.
2. Stage 2 (exploratory, DEV only, every step logged in POSTHOC_LOG.md): for each failed family identify the
   smallest modification that keeps the framework (a Markov walk on the lattice with same-site contact)
   and makes that family pass on DEV. Candidate classes named now: (F1) home-anchored walk — each person has a
   fixed home site (desk/seat), groups' homes are spatially clustered, moves are biased towards home with
   probability b (a confining drift, still Markov); (F2) sticky sites / broader heterogeneity of hopping;
   (F3) a scheduled alternation between an anchored and a free phase. "Minimal" = fewest extra parameters.
   Their parameters must be fitted from training days only by a fixed, automatic procedure.
3. The chosen modification(s), the fitting procedure and code are frozen in `FIX_FREEZE.json` (SHA-256 +
   timestamp) BEFORE any model or section-5 statistic is computed on a CONF dataset.
4. Stage 3 (confirmatory, out-of-sample across settings): base models and frozen modification(s) on CONF.
   The ladder of section 6 is applied to the modification(s) with n=3, with BLOCK as the competitor that must
   be beaten: the modification is reported as "adding something beyond a constant-ratio block model" only if it
   passes E on at least as many CONF datasets as BLOCK and more S-families in total.

## 9. Reporting

All results, including every failed family and failed scenario, in `memo.md` with tables generated from
`out/*.json`. Seeds fixed (seed 20261001). Runs < 4 GB, < 20 min each, checkpointed per dataset.

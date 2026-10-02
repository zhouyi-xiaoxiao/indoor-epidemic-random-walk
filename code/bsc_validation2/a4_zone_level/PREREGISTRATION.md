# Pre-registration — Avenue 4: zone-level structure

Directory: `research/bsc_validation2/a4_zone_level/`. Written 2026-10-01, before any model was fitted to or
compared with outbreak outcomes. Frozen by `PREREG_FREEZE.json` (SHA-256 + timestamp). Every later choice goes to
`POSTHOC_LOG.md`.

## 0. Status of blindness (stated honestly)

* **T1 (Matsumoto):** before freezing, the outbreak file was only converted from JLD2 to CSV and summarised for
  *structure* (schools, grades, classes, respondent counts; `src/00_convert_jld2.py`). The only outcome numbers seen are
  the total (2,548 infected of 10,923, published in the source paper) and the minimum onset code (20). No class-, grade-
  or school-level outcome tabulation was made, no model was fitted. The power analysis (`src/02_power.py`) reads the
  structure columns only and uses synthetic epidemics.
* **Prior knowledge that cannot be unseen:** the source paper (Endo et al. 2021, PNAS) fitted its own within/between-class
  model to these data. I recall its qualitative conclusion (within-school reproduction number below one, about 0.8, with
  the class the dominant level, then grade, then school), but not its numerical shares. I have not opened the
  posterior-sample files in the repository (`output/*.csv`) and will not before T1 is run. Consequently the *direction* of
  Tier 1 below is foreseeable from the literature and carries little surprise value; the *quantitative* Tier 2 is not known to me.
* **T2 (Seoul):** data-aware (the floor counts are in the project's outbreak dataset and were re-verified today), model-blind
  (the zone model has not been evaluated on them). I did a back-of-envelope calculation while designing the check
  (ε = 2 %, κ = 1 gives about 2 expected off-floor cases versus 3 observed), so T2 is **not blind** and is registered only as
  a low-power consistency check, never as a test the model can "pass robustly".

## 1. What is being tested

The theory of the model [...] says that at the scale of zones the epidemic is governed by the
next-generation operator `K = βQ(γ−M)^{-1}`; with fast movement between zones attack rates equalise, with slow movement
they separate. In its zone-level (Lagrangian, time-share) reduction a person of group `i` spends a fraction `τ_iz` of the
time in zone `z`; with same-site, frequency-dependent transmission [...] the force
of infection on group `i` is

    λ_i(t) = β Σ_z q_z τ_iz · (Σ_j τ_jz I_j(t)) / (Σ_j τ_jz N_j).

**Identified component under test:** *transmission between groups is proportional to independently measured co-location
(contact) time, with frequency-dependent normalisation.* This is the same-site transmission kernel ("M1")
at zone resolution. The lattice random walk itself, the within-room proximity gradient and the airborne term ("M2") are
**not** tested here.

## 2. Test T1 (primary, well powered): class / grade / school structure of seasonal influenza, Matsumoto 2014/15

### 2.1 Data (outbreak)
Endo A, Uchida M, Hayashi N, Liu Y, Atkins KE, Kucharski AJ, Funk S. Within and between classroom transmission patterns
of seasonal influenza among primary school students in Matsumoto city, Japan. PNAS 2021;118(46):e2112605118,
doi:10.1073/pnas.2112605118. Open data (MIT licence): https://github.com/akira-endo/schooldynamics_FluMatsumoto14-15,
file `data/anonymizedstudents.jld2`, commit 658abfecc70440b0b849e179da457048082d4e72, accessed 2026-10-01,
SHA-256 55b454bfde0be0ab1e92830a495362ddf138b879d657fade944dd5195ff8496f; raw copy in `data/raw/matsumoto/`.
10,923 respondent students, 29 public primary schools, 6 grades, 455 classes; per student: infected yes/no, onset day
(day 1 = 1 Oct 2014), school, grade, class. Zones = classes; denominators = respondents per class.

### 2.2 Data (mixing, independent of any outbreak)
SocioPatterns Lyon primary-school contact data: Stehlé J et al. PLoS ONE 2011;6(8):e23176, doi:10.1371/journal.pone.0023176;
Gemmetto V, Barrat A, Cattuto C. BMC Infect Dis 2014;14:695, doi:10.1186/s12879-014-0695-9. Files
`primaryschool.csv.gz`, `primaryschool_metadata.txt` from http://www.sociopatterns.org/datasets/primary-school-temporal-network-data/
(CC BY-NC-SA), accessed 2026-10-01, SHA-256 5c93d9f5…47fc9 and 92844d12…f4329; raw copies in `data/raw/sociopatterns/`.
232 children (teachers dropped), 10 classes (2 per grade, 5 grades), 2 school days, 20-s proximity records.
From `src/01_lyon_mixing.py` (run before freeze; it touches no outbreak data): mean contact time per student per day is
7,851 s with own class, 1,076 s with the other class of the same grade, 1,376 s with other grades, i.e. **time shares**

    w_Lyon = (class 0.762, grade 0.104, school 0.134)

(day 1: 0.756/0.111/0.133; day 2: 0.768/0.098/0.134). Per-pair contact time 352 : 46.4 : 7.42 s/day.

### 2.3 Models (all discrete-time daily chain-binomial; identical likelihood code)
A susceptible student of class `c`, grade `g`, school `s` has onset on day `t` with probability `1 − exp(−λ_c(t))`,

    λ_c(t) = β · x_c(t) + α0 + α1 · P_city\s(t),

where, with the onset-to-onset kernel `k(τ)` (gamma, mean 1.7 d, sd 1 d, daily bins τ = 1…10; Vink et al. 2014,
Am J Epidemiol 180:865, the H3N2 serial interval, as used by Endo et al.) and `K_A(t) = Σ_{j∈A} k(t − t_j)`:

    P_c = K_c/(n_c − 1),  P_g = K_{g\c}/n_{g\c},  P_s = K_{s\g}/n_{s\g},  P_city\s = K_{city\s}/N_{city\s}
    (a term is zero when its stratum is empty; n = respondents).

Within-school terms are set to zero for onset days 28 Dec–8 Jan (documented winter break 27 Dec–7 Jan, shifted by one
day, as in Endo et al.); the city term is not. Weekends are ignored. Sample weights and covariates are not used.

| Model | `x_c(t)` | fitted parameters |
|---|---|---|
| **Z** — zone model, mixing from Lyon, frequency-dependent transfer | `0.762 P_c + 0.104 P_g + 0.134 P_s` | β, α0, α1 |
| **H** — homogeneous mixing within the school ("well-mixed") | `K_s/(N_s − 1)` | β, α0, α1 |
| **X** — each class independent | `P_c` | β, α0, α1 |
| **C** — generic constant own-class/rest-of-school per-pair risk ratio ρ | `(ρK_c + K_{s\c})/(ρ(n_c−1) + N_s − n_c)` | β, α0, α1, ρ |
| **F** — free shares (reference/ceiling) | `b_c P_c + b_g P_g + b_s P_s` | b_c, b_g, b_s, α0, α1 |

Nothing in Z is calibrated on the outbreak except the overall scale β and the two external-hazard parameters, exactly as
in H and X. The mixing structure of Z comes from Lyon alone. `ŵ = b/Σb` denotes the shares fitted by F.

### 2.4 Calibration, hold-out, metrics
* Two-fold cross-validation by school: calibrate on odd `schoolID`, predict even, and vice versa. Metric: held-out
  one-day-ahead log predictive likelihood, summed per school; `D_H = Σ_s [LL_Z − LL_H]`, `D_X` likewise; 95 % CI by
  bootstrap over the 29 schools (10,000 resamples). C is evaluated the same way (reported, not decisional).
* Full-data fits of Z and F: likelihood ratio `LR = 2(LL_F − LL_Z)` and the total-variation distance
  `TV = ½ Σ|ŵ − w_Lyon|`; bootstrap CI of ŵ over schools (300 resamples).
* Forward simulation (statement (iii) of the theory: mixing strength ↔ separation of attack rates): 500 whole-city
  stochastic simulations from the full-data fit of each of Z, H, X; statistic = pooled Pearson heterogeneity of class
  attack rates around their school's attack rate (≈1 if a school is well mixed), and the city attack rate.

### 2.5 Decision rule
* **Tier 1 (ranking, out of sample):** lower 95 % bootstrap bound of `D_H` > 0 **and** of `D_X` > 0.
* **Tier 2a (strict quantitative):** `LR ≤ 5.99` (χ²₂, 95 %).
* **Tier 2b (tolerance):** `TV ≤ 0.15`.
* **Tier 3 (forward dispersion):** observed heterogeneity statistic inside the central 95 % predictive interval of Z.
* **Verdict:** "pass" = Tier 1 ∧ Tier 2b ∧ Tier 3. Tier 2a is reported separately as the strict test; with ~2,500 cases it
  detects small deviations, and a French school is not a Japanese school, so failing 2a while passing 2b will be reported
  as "right structure within tolerance, not exact". Failing Tier 2b or Tier 3 is a failure of the quantitative
  prediction even if Tier 1 holds; Tier 1 alone is **not** to be reported as a validation (a generic model can pass it,
  see 2.6). The model's forward-simulated city attack rate is reported but is not a test of "levels": β and α are fitted.

### 2.6 Power analysis (done before freeze, synthetic data, real class structure; `src/02_power.py`, `results/power_*.jsonl`)
Synthetic epidemics on the actual 455-class structure; 40 replicates per truth per scenario; the identical pipeline
(`src/analysis.py`) is applied. Scenario A: β = 0.8, α0 = 2·10⁻⁵/d, α1 = 0.42 (city attack rate ≈ 0.2); scenario B: β = 0.5,
α1 = 0.70 (weaker school signal). Truths: Z (Lyon shares), H, X, C with ρ log-uniform on [1, 1000] ("generic constant
risk ratio"), F with shares uniform on the simplex ("generic random structure").

Probability of passing each tier (40 synthetic replicates per cell; binomial s.e. ≤ 0.08):

| scenario | truth | Tier 1 | Tier 2a | Tier 2b | Tier 3 | **verdict "pass"** |
|---|---|---|---|---|---|---|
| A | Z (Lyon zone model true) | 1.00 | 0.95 | 1.00 | 1.00 | **1.00** |
| A | H (well mixed) | 0.00 | 0.00 | 0.00 | 0.18 | **0.00** |
| A | X (independent classes) | 0.00 | 0.00 | 0.00 | 0.03 | **0.00** |
| A | C (generic constant risk ratio, ρ random) | 0.30 | 0.00 | 0.23 | 0.45 | **0.075** |
| A | F (generic random shares) | 0.63 | 0.00 | 0.03 | 0.33 | **0.025** |
| B | Z | 0.85 | 0.95 | 1.00 | 0.95 | **0.80** |
| B | H | 0.00 | 0.00 | 0.00 | 0.58 | **0.00** |
| B | X | 0.00 | 0.00 | 0.00 | 0.85 | **0.00** |
| B | C | 0.25 | 0.03 | 0.18 | 0.78 | **0.00** |
| B | F | 0.53 | 0.05 | 0.10 | 0.80 | **0.05** |

Pooled over scenarios, by the true ρ of the generic constant-risk-ratio model: ρ < 10: verdict pass 0/24; ρ 10–30: Tier 1
passes 12/12 but verdict 0/12; ρ 30–100: verdict 3/18; ρ ≥ 100: 0/26. Generic random shares farther than TV 0.15 from Lyon:
verdict pass 0/75. Under truth Z the 95th percentile of TV is 0.045 (A) / 0.071 (B) and of LR 4.7–4.8.
Tier 3 alone discriminates poorly when the school signal is weak (scenario B: the observed
statistic falls inside Z's interval for 58 % of H-truths and 85 % of X-truths), so Tier 3 is a necessary check, not evidence.

Reading: if the Lyon-calibrated zone model is true it passes with probability ≈ 1. If schools are well mixed (H) or
classes independent (X) the pass probability is ≈ 0. A generic constant-risk-ratio model passes Tier 1 in a broad band of ρ
but passes the full rule only when its ρ happens to mimic Lyon (ρ ≈ 30–100; Lyon's own-class vs rest-of-school per-pair
ratio is ≈ 30). A generic random structure passes only if it lies within TV ≤ 0.15 of Lyon, which is 11.4 % of the simplex.
These are false-pass rates *of the truth-generating alternatives*, under a correctly specified likelihood; real data are
over-dispersed (individual variation, recall error, class closures), which can only make passing harder for every model.

### 2.7 Pre-declared sensitivity analyses (non-decisional, full-data fits)
Lyon day 1 / day 2 shares; density-dependent transfer of Lyon per-pair times (per-pair weights 352 : 46.4 : 7.42 applied to
`K_c`, `K_{g\c}`, `K_{s\g}`); kernel mean 2.5 d, sd 1.5 d; no winter-break switch.

## 3. Test T2 (secondary, LOW power, not blind): Seoul call-centre building by floor

### 3.1 Data
Park SY et al. Coronavirus Disease Outbreak in Call Center, South Korea. Emerg Infect Dis 2020;26(8):1666–1670,
doi:10.3201/eid2608.201274, Table 1 (verified against the Europe PMC XML of PMC7392450, copy in `data/raw/park2020/`,
accessed 2026-10-01): 11th floor 94/216; floors 1–6 0/84; 7th 0/182; 8th 0/207; 9th 1/206; 10th 2/27; residents
(13th–19th) 0/201; visitors 0/20; total 97/1,143. Text: "Residents and employees in building X had frequent contact in the
lobby or elevators"; the first case by onset (22 Feb) worked on the 10th floor and "reportedly never went to 11th floor".

### 3.2 Model and independent information
Zones: each floor plus one shared zone (lifts + lobby). Every occupant spends a fraction ε of on-site time in the shared
zone. Frequency-dependent same-zone transmission with efficiency q on floors and κq in the shared zone. Cumulative hazard
of floor f: `Λ_f = c[(1−ε) z_f + κ ε z_tot]` (z = cumulative infected fraction). `c` is calibrated on the 11th floor alone
(`Λ_11 = −ln(1 − 94/216)`, with `z_tot = 94/1143`); prediction: first-generation spill-over to the 927 other occupants,
`E = 927 (1 − exp(−c κ ε · 94/1143))`. On-floor amplification on other floors (truncated by the closure on 9 March) is not
modelled, so the prediction is a first-order lower bound.
There are **no measured** movement data for this building. ε and κ are therefore priors, not calibrations:
shared-zone time 3–30 min per 9-h day (log-uniform; ε = 0.0056–0.056; motivated by lift-design criteria of ≤ 30 s
average waiting time and trips of 1–2 min, 4–8 trips per day — CIBSE Guide D as summarised at
designingbuildings.co.uk "Lifts for office buildings"), κ log-uniform on [1, 5] (crowded cars). This is the weak point of T2.

### 3.3 Decision rule and power
Statistic: number of cases among the 927 non-11th-floor occupants (observed 3; 2 if the 10th-floor first case is regarded
as a source rather than a recipient). "Consistent" if the observed count lies in the central 95 % prior-predictive
interval (Poisson mixture over the ε, κ prior). Competitors: well mixed (same attack rate everywhere: expected ≈ 79) and
independent floors (expected 0). **Power is low:** a generic predictor whose expected count is log-uniform between 0.01 and
400 has a 95 % Poisson interval containing 3 with probability ln(8.77/0.62)/ln(40000) ≈ 0.25, and the present model, whose
prior on κε spans a factor 50, will have a wide interval. T2 can therefore only be reported as "inconclusive / consistent
in order of magnitude" or "fail"; it cannot be reported as a robust pass.
The two wings of the 11th floor (79/137 vs 4/62 digitised seats): no independent measurement of between-wing mixing
exists, so **no test is registered**; the coupling that would reproduce the split may be reported as a descriptive,
fitted quantity only.

## 4. Other candidate datasets considered and not used (before freeze)
Diamond Princess (deck / crew vs passengers), USS Theodore Roosevelt, meat plant (Günther et al. 2020), prison
dormitories, nursing homes, school outbreaks by grade (e.g. Stein-Zamir et al. 2020): zone denominators exist in several,
but I found no *independently measured* between-zone mixing for any of them, so a zone model could only be fitted, not
tested. They are listed in the memo, not analysed.

## 5. Reporting
All tiers, all models, both folds, all sensitivity analyses and T2 are reported whatever the outcome. Seeds fixed
(20261001). Runs < 4 GB, < 20 min each.

# Reproducing the article

All commands are run from the repository root. Every script locates its inputs relative to its own position,
so the working directory does not matter for the Python commands; the shell lines below assume the root.
Figure and table numbers refer to the main text (`paper.pdf`) or, with the prefix S, to the Supplementary
Material (`supplement.pdf`).

## 0. Environment

```sh
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

* Python 3.14.6 with NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, mpmath 1.3.0, pandas 3.0.6 and Pillow 12.3.0
  produced the stored results.
* A C compiler available as `cc` is needed: the simulators (`code/bsc_sim/src/bsc_sim/gillespie.c`,
  `code/bsc_sim/verify/vsim.c`, `vwalk.c`) are compiled automatically with `cc -O2 -shared -fPIC` on first use.
* The Chinese figure variants (`*_zh.pdf`) need one of the fonts Songti SC, STSong, PingFang SC, Heiti SC,
  STHeiti or Arial Unicode MS; without one the figure scripts stop at the Chinese variant.
* LaTeX: TeX Live 2025. The main text and the Supplementary Material refer to each other through `xr-hyper`, so
  compile `latexmk -pdf main.tex`, then `latexmk -pdf supplement.tex`, then `latexmk -pdf main.tex` again.
* `code/bsc_<name>/data` etc. are symbolic links into `data/` (see `README.md`).

Conventions that matter when re-running:

* **Seeds.** Every stochastic script sets its seeds in its source, so a re-run reproduces the stored numbers
  exactly on the same platform (checked for `03b_outbreak_probability.py`, whose output is identical to the
  stored file). Fields that record run times, dates and speed-up factors change from run to run.
* **Checkpoints.** The simulation experiments save their JSON after each configuration and skip
  configurations that are already present. To force a full re-run, move the output file away first.
* **Generated tables.** The bodies of Tables 4, S8 and S9 are written to `data/s6_scenes_tables.tex` and pasted
  into `sections/m6_simulation.tex` and `sections/supp_sim.tex` between marker comments by
  `code/s6_scenes_paste_tables.py`; those of Tables S16, S17 and S18 are injected into
  `sections/supp_datasets.tex` by `code/appB_dataset_table.py`; the generated rows of Tables S27 to S29 must occur
  verbatim in `sections/supp_protocol1.tex`, which `code/appC_protocol_tables.py` checks. The bodies of Tables 8,
  S19, S31, S32, S36 and S37 are written to `data/v2_tab_*.tex` by `code/v2_numbers.py` and read by the LaTeX source
  directly. The rows of Tables 5, S10 to S14 are generated into `data/s7_interventions_tables.tex`, and
  `python code/s7_interventions_numbers.py --check` verifies that each occurs verbatim in the LaTeX source; the rows
  of Tables 2 and 3 are generated into `data/s2_model_table_rows.tex`, those of Tables S6 and S7 into
  `data/s5_pair_tables.tex`. All other tables were transcribed from the results files named in their `% src:`
  comments and checked against them.
* **Figures.** Each figure script writes an English and a Chinese variant (`_en`, `_zh`) as PDF and PNG.
  Article-level scripts (`code/*.py`) write directly into `figures/`. The research scripts write into
  `code/bsc_<name>/figures/`; copy the PDFs listed in Section 3 into `figures/`.

## 1. Quick check: article-level numbers and figures from the stored results

These commands read the stored results in `data/` and regenerate the numbers quoted in the text, the table bodies
and most of the figures.

```sh
python code/plan_theory_checks.py          # mean-field quantities on the four scenes -> data/plan_theory_checks.json, plan_theory_sweeps.npz (20 s to 4 min)
python code/s2_model_scene_table.py        # Tables 2, 3, S1; numbers of Section 2
python code/s3_r0_checks.py                # numbers of Section 3 and S2
python code/s4_geometry_checks.py          # Tables S3, S4; numbers of Sections 4 and S3
python code/s4_geometry_barriers_leak.py   # walls and barrier cells at D0 = 1, 10, 100 (Table S4, Fig. S2b); leak with heterogeneous q (about 1 min)
python code/s4_geometry_leak_kernels.py    # Remark S3.1 (leak formula and bounds with other contact kernels); office door against mobility
python code/plan_fig_theory.py             # Figs 2, 3, S1, S24 (also writes fig_meanfield_interventions, not used)
python code/s4_geometry_fig_interventions.py   # Fig. S2
python code/fig_finite_size.py             # Fig. 4
python code/s5_pair_checks.py              # Tables S6, S7; numbers of Sections 5 and S4
python code/s6_scenes_ode_pointseed.py     # last column of Table 4
python code/s6_scenes_tables.py            # Tables 4, S8, S9; numbers of Sections 6 and S5
python code/s6_scenes_paste_tables.py      # pastes them into sections/m6_simulation.tex and sections/supp_sim.tex
python code/s6_scenes_fig_spatial.py       # Fig. S6
python code/s7_interventions_numbers.py    # Tables 5, S10-S14; numbers of Section 6.4 and S5
python code/s7_interventions_fig_hetero.py # Fig. S7
python code/s8_outbreaks_checks.py         # model-free statistics of Section 7 (pooled near/far ratio with fixed and random effects)
python code/s8_default_level.py            # the bound on infections in a short event at the default parameters (Section 7.10)
python code/s8_outbreaks_numbers.py        # numbers of the first outbreak test (Tables S20 to S22)
python code/s8_outbreaks_fig_callcentre.py # Fig. S11
python code/s8_first_test_exact.py --force # first test with the round-off-safe propagator: Tables 7, S23 (A1), S27(b), S28; Figs S13, S14 (about 2 min; without --force the parts already stored in data/s8_first_test_exact.json are kept)
python code/v2_numbers.py                  # further analyses: every number of Sections 7 and S9, S10 -> data/v2_numbers.json; bodies of Tables 8, S19, S31, S32, S36, S37; re-checks the hashes of the four freeze records
python code/v2_figures.py                  # Figs S15 to S23
python code/v2_tidy_data.py                # data/outbreaks2/*.csv from the transcribed seat maps; copies of the tracer and school tables
python code/s9_discussion_numbers.py       # numbers of Section 8
python code/misc_checks.py                 # small facts (desk mapping, mobility gap, reference counts) -> data/misc_checks.json
python code/appA_numerics_checks.py derived    # derived numbers of Section S11
python code/appB_dataset_table.py          # Tables S16, S17, S18
python code/appB_dataset_fig_attack.py     # Fig. S10
python code/appC_protocol_tables.py        # Tables S27, S28, S29; re-checks the six frozen hashes
python code/plan_validation_checks.py      # data/validation_checks/* (call-centre probability by quadrature and its weight sensitivity, generic benchmark; --quadrature-only for the first part)
latexmk -pdf main.tex && latexmk -pdf supplement.tex && latexmk -pdf main.tex
```

Note on `appC_protocol_tables.py`: besides the SHA-256 hashes it records the modification times of the protocol
files. In a clone these are checkout times, so the block `freeze.file_times_local` of
`data/appC_protocol_numbers.json` cannot be reproduced; the stored values are the record of the author's working
tree (Section S8.1 says that this evidence is file metadata only). The same holds for the block `freeze.*.times` of
`data/v2_numbers.json` (Section S10.1): `code/v2_numbers.py` recomputes the SHA-256 of every file named in the freeze
records of the further analyses, and in a clone it records checkout times; hashes of source files that are not
redistributed cannot be recomputed.

## 2. Full re-run of the research code

### 2.1 Mean-field theory (`code/bsc_theory/`; about 20 minutes)

```sh
python code/bsc_theory/scripts/verify_theorems.py   # Table S40 -> data/bsc_theory/verify_theorems.json (15-40 s)
python code/bsc_theory/scripts/worked_examples.py   # Table S3, room-size numbers (2-7 min)
python code/bsc_theory/scripts/nonlinear_check.py   # Section S11.4 (nonlinear mean-field model; 1-4 min)
python code/bsc_theory/scripts/ibm_check.py         # annealed encounter formula against simulation (about 14 min)
python code/bsc_theory/scripts/make_figures.py 2    # writes data/bsc_theory/fig2_room_size_data.json, the input of Fig. S1
```

The re-checks are the scripts `code/bsc_theory/verify/v1*.py … v5*.py`; they are written to be run from inside
that directory (`cd code/bsc_theory/verify && python v1_theorems.py`), and their logs are stored next to them.
`python code/appA_numerics_checks.py adversarial` repeats the adversarial searches of Section S11.3 in 60-digit
arithmetic (about 1.5 minutes).

### 2.2 Exact simulation (`code/bsc_sim/`; 3 to 4 hours on one core, below 1.5 GB)

```sh
cd code/bsc_sim
python experiments/00_build_scenes.py          # scenes/*.json (Fig. 1, Tables 2, S1)
python experiments/01_validate.py              # validation tests A-G, Table S41 (about 8 min)
python experiments/02_r0_regimes.py            # Tables S6, S7; Figs 4, 5
python experiments/03_scenes.py                # Table 4; Figs S4, S5
python experiments/03b_outbreak_probability.py # Table S8 (about 2 min)
python experiments/04_spatial_pattern.py
python experiments/04b_spatial_first_generation.py   # Table S9, Fig. S6
python experiments/05_heterogeneity.py
python experiments/05_heterogeneity.py extra   # Table S10, Fig. S7
python experiments/06_obstacles.py             # Fig. S8
python experiments/07_occupancy.py             # Fig. S9
python experiments/08_interventions.py         # Table S12
python experiments/09_metro_time.py            # Table S13 (D0 = 1, 100)
python experiments/10_benchmark.py             # Table S15
python experiments/11_dwell_time.py            # Table 5
python experiments/20_figures.py scenes mobility validity_map scene_curves scene_attack obstacle_scan occupancy validation
cd ../..
cp code/bsc_sim/figures/{fig3_scene_maps,fig5_5_R_vs_mobility,fig4_meanfield_validity_map,fig5_3_scene_epidemic_curves,fig5_4_attack_rate_distributions,fig6_3_obstacle_scan,fig6_4_occupancy,figA_validation}_{en,zh}.pdf figures/
```

Fig. 4 is drawn by `code/fig_finite_size.py` from `data/bsc_sim/02_finite_size.json`.

Additional simulations made for the article (they use the same simulator):

```sh
python code/s7_interventions_barriers.py 1     # Table S11, D = 1: 40 layouts x 1,000 runs per arm (about 2 min)
python code/s7_interventions_barriers.py 10    # Table S11, D = 10 (about 8 min)
python code/s7_interventions_schedule.py       # Table 5 (exact pair-level duty-cycle values; office counts of 200,000 to 2,000,000 runs), Table S13 (D0 = 10)
python code/s5_pair_large_sample.py            # counted secondary cases of the office scene, 100,000 to 400,000 runs (about 25 min)
python code/s7_interventions_hetero_se.py      # Table S10 and Fig. S7d: standard errors; 2,000 epidemics per contrast at N = 1000 (about 40 min)
python code/s7_interventions_regime.py det     # Table S14 (about 3 min)
python code/s7_interventions_regime.py sim     # simulations at D0 = 100 quoted in Section 6.4 (about 12 min)
python code/appA_numerics_checks.py simulator  # large-sample final-size test of Section S11.6 (about 8 min)
```

The second simulator and its checks are in `code/bsc_sim/verify/` (`vsim.c`, `vlib.py`, `v01_det.py` …
`v14_schedule.py`); run them from inside that directory. Their stored outputs (`v*.json`, `v*.log`) are the
values quoted as "second simulator" or "re-check" in the article.

### 2.3 First outbreak dataset (`code/bsc_outbreaks/`)

```sh
cd code/bsc_outbreaks
python scripts/build_dataset.py       # scripts/records.py + mapping.py -> outbreaks.csv, outbreaks.json, sources.md, data/*.csv
python scripts/analyze_dataset.py     # model-free statistics, permutation test (seed 20261001) -> data/analysis_results.json
cd ../..
```

`scripts/records.py` is the single place where numbers transcribed from the sources are typed. Not reproducible
from this repository alone: `scripts/verify_dois.py`, `fetch_epmc.py`, `fetch_bioc.py` need network access;
`scripts/digitize_callcentre.py` needs the image of Fig. 2 of Park et al. (2020), which is not redistributed (its
output, `data/bsc_outbreaks/tables/park2020_callcentre_seats_digitized.csv`, is included). The re-computation of
the re-check is in `code/bsc_outbreaks/verify/scripts/` (`v02_stats.py`, `v03_callcentre.py`,
`v04_extra_and_figures.py`; the source texts they read are not redistributed), with outputs in `verify/results/`.

### 2.4 First outbreak test (`code/bsc_validation/`; order matters)

```sh
cd code/bsc_validation
python scripts/00_selftest.py          # synthetic inputs only
python scripts/01_calibrate.py         # calibration records only (trains, restaurant), all variants -> Table S20, Fig. S12
python scripts/02_power.py             # hold-out predictive distributions and power; reads no stratum outcome -> Table S27(a)
python scripts/02b_m1_exact.py         # same-site model in the exact simulator -> Table S21
python scripts/03_holdout_shape.py     # first script that reads data/heldout_outcomes.json (frozen propagator; the article prints the values of code/s8_first_test_exact.py)
python scripts/02c_mc_stability.py
python scripts/04_level_controls.py    # Tables S22, S29
python scripts/06_callcentre.py        # 40-draw Monte Carlo estimate of the call-centre test (the article reports the quadrature of code/plan_validation_checks.py)
python scripts/07_loo.py               # leave-one-event-out curves (frozen propagator)
python scripts/08_bruteforce_check.py  # brute-force check of exposures and conditional distributions
python scripts/09_sensitivity.py       # sensitivity analyses with the frozen propagator (the article prints the rows of code/s8_first_test_exact.py)
python scripts/10_decision.py          # decision rule with the frozen propagator and the 40-draw call-centre estimate
python scripts/20_figures.py           # validation figures
python scripts/30_tables.py            # data/tables.md
cd ../..
cp code/bsc_validation/figures/figV1_calibration_{en,zh}.pdf figures/
# Figs S13 and S14 are drawn by code/s8_first_test_exact.py with the round-off-safe propagator
```

The protocol documents `PREREGISTRATION.md`, `PREREG_ADDENDUM_A1_power.md`, `PREREG_FREEZE.json` and
`POSTHOC_LOG.md` are included as redacted copies (`code/REDACTIONS.md` lists every changed line);
`python code/appC_protocol_tables.py` checks each redacted copy against the hash listed in `code/REDACTIONS.json`
and compares the hashes of the originals listed there with the six SHA-256 hashes of the freeze record. The re-implementation of the re-check is in `code/bsc_validation/verify/` (`vlat.py`,
`vevents.py`, `v01_cal_T4.py` … `v13_c3.py`; run from inside that directory). `python code/s8_outbreaks_checks.py
--check` re-runs its sensitivity script (`v09_sens.py`, several minutes) and writes `data/s8_outbreaks_v09_sens.txt`.

### 2.5 Further analyses (`code/bsc_validation2/`, `code/checks/`)

Each of the four analyses has a `README.md` with its commands; the order inside an analysis matters, because the
protocol documents record what was computed before the outcomes were read.

```sh
python code/download_contact_data.py --schools       # third-party files that are not redistributed (network)

# second outbreak test (Section 7; Tables 8, S35, S36; Figs S15, S16, S22)
cd code/bsc_validation2/a1_more_outbreaks
python scripts/00_check_data.py && python scripts/00_selftest.py
sh scripts/run_tables.sh F1 F2 F3 F4 F5 F6 F7 F8 W1 W2 P1          # exposure tables (not stored: about 80 MB)
for e in F1 F2 F3 F4 F5 F6 F7 F8 W1 W2 P1; do python scripts/01b_tables_fix.py $e; done
python scripts/02_calibrate.py && python scripts/03_power.py      # calibration events only; power
python scripts/04_holdout.py                                      # held-out events, kernel as first frozen
python scripts/05_secondary.py; sh scripts/run_variants.sh; python scripts/06_variants.py
python scripts/07_descriptive.py; python scripts/09_numerics_check.py
python scripts/10_corrected.py                                    # round-off-safe kernel: the values of the article
cd ../../..

# contact records (Section 7; Tables S31, S37; Figs S17, S23)
cd code/bsc_validation2/a2_contact_mobility/src
python 02_run.py s1 InVS13,LyonSchool,LH10 WM,BLOCK,RW0,RWhet     # see the README for the other runs
cd ../../../..

# tracer measurements (Section 7; Tables S32, S38; Figs S18, S19)
sh code/bsc_validation2/a3_tracer_physics/run_all.sh

# zone level (Section 7; Table S33; Fig. S20)
cd code/bsc_validation2/a4_zone_level && python src/04_matsumoto.py && python src/05_seoul_floors.py && cd ../../..
```

The re-implementations of the re-checks are in `code/checks/` (see its `README.md`): for example
`cd code/checks/outbreaks2 && python t2_analysis.py` recomputes the second outbreak test with its own geometry code
and propagator (it rebuilds its exposure tables on the first run) and writes `t2_results.json`, and
`cd code/checks/zones && python v02_main.py` recomputes the zone-level test. The comparisons that the checks added
after the results were known (other near zones and strata, the exponential kernel, the static contact network, the
fits of the modified walks to their own synthetic data, guessed mixing shares) are produced by these scripts and by
no other.

### 2.6 Bibliography and numbers from the checks

Every entry of `refs.bib` was generated from, or checked against, the Crossref record of its DOI; the records are
stored in `data/plan_refs_crossref.json` (with the comparison in `data/plan_refs_report.json`),
`data/refs_search.json` (written by `code/refs_search.py`, network; it also stores the search queries behind the
statements of what is added and the manual checks of the patch-model literature) and `data/v2_refs_crossref.json`
(written by `code/v2_refs.py`, network). The files `data/checks/*.json` hold numbers recorded by the re-checks;
they were extracted from the records of the re-checks and are not regenerated by the scripts of this repository;
when the repository was assembled, each value was tested to be written in the item of `notes/VERIFICATION.md` it is
attributed to.

## 3. Where each figure comes from

| Figure | File in `figures/` | Generated by | From |
|---|---|---|---|
| 1 | `fig3_scene_maps` | `code/bsc_sim/experiments/20_figures.py scenes` | `code/bsc_sim/scenes/*.json` |
| 2 | `fig_r0_mobility` | `code/plan_fig_theory.py` | `data/plan_theory_checks.json`, `data/plan_theory_sweeps.npz` |
| 3 | `fig_hotspot_maps` | `code/plan_fig_theory.py` | `data/plan_theory_sweeps.npz` |
| 4 | `fig5_1_finite_size` | `code/fig_finite_size.py` | `data/bsc_sim/02_finite_size.json` |
| 5 | `fig5_5_R_vs_mobility` | `20_figures.py mobility` | `data/bsc_sim/02_mobility.json` |
| S1 | `fig_room_size_leak` | `code/plan_fig_theory.py` | `data/bsc_theory/fig2_room_size_data.json` |
| S2 | `fig_s4_interventions` | `code/s4_geometry_fig_interventions.py` | `data/plan_theory_checks.json`, `data/plan_corridor_sweep.json`, `data/s4_geometry_barriers_leak.json` |
| S3 | `fig4_meanfield_validity_map` | `20_figures.py validity_map` | closed form (Corollary 5.5) |
| S4 | `fig5_3_scene_epidemic_curves` | `20_figures.py scene_curves` | `data/bsc_sim/03_scenes.json`, `03_curves_*.npz` |
| S5 | `fig5_4_attack_rate_distributions` | `20_figures.py scene_attack` | same |
| S6 | `fig5_2_spatial_pattern` | `code/s6_scenes_fig_spatial.py` | `data/bsc_sim/04_spatial_office_D1.npz`, `04_spatial.json`, `data/s6_scenes_numbers.json` |
| S7 | `s7_interventions_hetero` | `code/s7_interventions_fig_hetero.py` | `data/bsc_sim/05_heterogeneity.json`, `05_heterogeneity_extra.json`, `data/s7_interventions_hetero_se.json` |
| S8 | `fig6_3_obstacle_scan` | `20_figures.py obstacle_scan` | `data/bsc_sim/06_obstacles.json` |
| S9 | `fig6_4_occupancy` | `20_figures.py occupancy` | `data/bsc_sim/07_occupancy.json` |
| S10 | `appB_dataset_attack` | `code/appB_dataset_fig_attack.py` | `data/appB_dataset_numbers.json` |
| S11 | `s8_outbreaks_callcentre` | `code/s8_outbreaks_fig_callcentre.py` | `data/bsc_outbreaks/tables/park2020_callcentre_seats_digitized.csv`, `callcentre_joincount_null.npy`; panel (c): `data/validation_checks/callcentre_quadrature.json`, `data/s8_outbreaks_numbers.json`, `data/bsc_validation2/a3_tracer_physics/04_evaluation.json`, `05_descriptive.json` |
| S12 | `figV1_calibration` | `code/bsc_validation/scripts/20_figures.py` | `data/bsc_validation/01_cal_T4_primary.json`, `01_cal_R1_primary.json` |
| S13 | `figV2_holdout_shape` | `code/s8_first_test_exact.py figures` | `data/s8_first_test_exact.json`, `data/bsc_validation/02_predictive.json`, `02b_m1_exact.json` |
| S14 | `figV8_mixing_heterogeneity` | `code/s8_first_test_exact.py figures` | `data/bsc_validation/07_loo.json`, classroom curve from `data/s8_first_test_exact.json` |
| S15 | `v2_rr_forest` | `code/v2_figures.py v2_rr_forest` | `data/bsc_validation2/a1_more_outbreaks/07_descriptive.json` |
| S16 | `v2_holdout_strata` | `code/v2_figures.py v2_holdout_strata` | `data/bsc_validation2/a1_more_outbreaks/10_corrected.json` |
| S17 | `v2_contacts_epidemic` | `code/v2_figures.py v2_contacts_epidemic` | `data/bsc_validation2/a2_contact_mobility/s*_<dataset>.json`, `code/checks/contacts/v4_static_*.json` |
| S18 | `v2_tracer_kernels` | `code/v2_figures.py v2_tracer_kernels` | `data/tracer/*.csv`, `data/bsc_validation2/a3_tracer_physics/01_kernels.json`, `05_descriptive.json` |
| S19 | `v2_tracer_predictions` | `code/v2_figures.py v2_tracer_predictions` | `data/bsc_validation2/a3_tracer_physics/04_evaluation.json` |
| S20 | `v2_zone` | `code/v2_figures.py v2_zone` | `data/bsc_validation2/a4_zone_level/matsumoto_primary.json`, `seoul_floors.json`, `data/schools/lyon_mixing.json` |
| S21 | `v2_summary_scores` | `code/v2_figures.py v2_summary_scores` | `data/v2_numbers.json`, `data/bsc_validation2/a1_more_outbreaks/10_corrected.json`, `a3_tracer_physics/04_evaluation.json`, `E2_ext_evaluation.json` |
| S22 | `v2_D_profiles` | `code/v2_figures.py v2_D_profiles` | `data/bsc_validation2/a1_more_outbreaks/10_corrected.json` |
| S23 | `v2_contacts_structure` | `code/v2_figures.py v2_contacts_structure` | `data/bsc_validation2/a2_contact_mobility/s*_<dataset>.json` |
| S24 | `fig_theorem_checks` | `code/plan_fig_theory.py` | `data/bsc_theory/T1_points.npy`, `T2_positions.npy`, `data/plan_theory_checks.json` |
| S25 | (algorithm box, typeset in `sections/supp_verification.tex`) | — | — |
| S26 | `figA_validation` | `20_figures.py validation` | `data/bsc_sim/01_validation.json`, `01_B_finalsize.npz`, `01_F_curves.npz` |

`20_figures.py` stands for `code/bsc_sim/experiments/20_figures.py`, which also loads `21_figures_more.py`.

## 4. Where each table comes from

| Table | Content | Results file(s) | Script that writes the entries |
|---|---|---|---|
| 1 | known results and what is added | as cited in the `% src:` comments of `sections/m1_intro.tex` | — |
| 2, 3, S1 | scenes; R0 by mobility and kernel; zones | `data/s2_model_scenes.json` (cross-checked against `data/plan_theory_checks.json`, `data/s3_r0_checks.json`) | `code/s2_model_scene_table.py` |
| 4 | scenes in simulation | `data/bsc_sim/03_scenes.json`, `data/s6_scenes_ode_pointseed.json` | `experiments/03_scenes.py`, `code/s6_scenes_ode_pointseed.py`, `code/s6_scenes_tables.py` |
| 5 | dwell time | `data/bsc_sim/11_dwell_time.json`, `data/s7_interventions_schedule.json` | `experiments/11_dwell_time.py`, `code/s7_interventions_schedule.py` |
| 6 | models | — (definitions) | — |
| 7 | hold-out events, first test | `data/s8_first_test_exact.json`, `data/bsc_validation/03_holdout_shape.json` | `code/s8_first_test_exact.py`; `scripts/02_power.py`, `02b_m1_exact.py`, `03_holdout_shape.py` |
| 8 | hold-out events, second test | `data/bsc_validation2/a1_more_outbreaks/10_corrected.json`, `07_descriptive.json` | `code/bsc_validation2/a1_more_outbreaks/scripts/10_corrected.py`, `code/v2_numbers.py` |
| 9 | components of the model after the tests | — (summary of Section 7) | — |
| S2 | notation | — (definitions) | — |
| S3 | room size and leaks | `data/bsc_theory/worked_examples.json`, `data/s4_geometry_checks.json` | `code/bsc_theory/scripts/worked_examples.py`, `code/s4_geometry_checks.py` |
| S4 | worked examples, office | `data/plan_theory_checks.json`, `data/s4_geometry_checks.json`, `data/s4_geometry_barriers_leak.json` | `code/plan_theory_checks.py`, `code/s4_geometry_checks.py`, `code/s4_geometry_barriers_leak.py` |
| S5 | targeted ventilation, office | `data/plan_theory_checks.json` | `code/plan_theory_checks.py` |
| S6, S7 | finite size; mobility | `data/bsc_sim/02_finite_size.json`, `02_mobility.json` | `experiments/02_r0_regimes.py`; rows in `data/s5_pair_tables.tex` by `code/s5_pair_checks.py` |
| S8 | outbreak probability | `data/bsc_sim/03b_outbreak_probability.json` | `experiments/03b_outbreak_probability.py`, `code/s6_scenes_tables.py` |
| S9 | infection maps | `data/bsc_sim/04_spatial.json`, `04b_first_generation.json`, `04_spatial_office_D*.npz` | `experiments/04_spatial_pattern.py`, `04b_spatial_first_generation.py`, `code/s6_scenes_tables.py` |
| S10 | two-zone room | `data/bsc_sim/05_heterogeneity.json`, `05_heterogeneity_extra.json`, `code/bsc_sim/verify/v09_het.json`, `data/s7_interventions_hetero_se.json` | `experiments/05_heterogeneity.py`, `code/s7_interventions_hetero_se.py`, `code/s7_interventions_numbers.py` |
| S11 | barriers | `data/s7_interventions_barriers.json`, `code/bsc_sim/verify/v08b_obstacles_many.json` | `code/s7_interventions_barriers.py`, `code/s7_interventions_numbers.py` |
| S12 | interventions | `data/bsc_sim/08_interventions.json` | `experiments/08_interventions.py`, `code/s7_interventions_numbers.py` |
| S13 | daily schedule | `data/bsc_sim/09_metro_time.json`, `data/s7_interventions_schedule.json` | `experiments/09_metro_time.py`, `code/s7_interventions_schedule.py` |
| S14 | effects against mobility | `data/s7_interventions_regime.json` | `code/s7_interventions_regime.py` |
| S15 | cost | `data/bsc_sim/10_benchmark.json`, `code/bsc_sim/verify/v12_benchmark_rerun.json`, `v13_bench.json` | `experiments/10_benchmark.py`, verification scripts `v12`, `v13` |
| S16, S17, S18 | first dataset; strata; controls | `data/bsc_outbreaks/outbreaks.csv`, `tables/near_far_risk_ratios.csv`, `tables/analysis_results.json` | `code/bsc_outbreaks/scripts/build_dataset.py`, `analyze_dataset.py`, `code/appB_dataset_table.py` |
| S19 | second dataset | `data/bsc_validation2/a1_more_outbreaks/07_descriptive.json`, `code/bsc_validation2/a1_more_outbreaks/PREREGISTRATION.md` | `code/v2_numbers.py` |
| S20 | calibration, first test | `data/bsc_validation/01_cal_T4_primary.json`, `01_cal_R1_primary.json`, `tables.md` | `scripts/01_calibrate.py`, `scripts/30_tables.py` |
| S21 | reach ceiling of the same-site model | `data/bsc_validation/02b_m1_exact.json` | `scripts/02b_m1_exact.py` |
| S22 | implied emission | `data/bsc_validation/04_level_controls.json`, `data/s8_outbreaks_v09_sens.txt` | `scripts/04_level_controls.py`, `code/s8_outbreaks_checks.py --check` |
| S23 | outcome of the first test | `data/bsc_validation/10_decision.json`, `data/validation_checks/*`, `data/s8_first_test_exact.json` | `scripts/10_decision.py`, `code/plan_validation_checks.py`, `code/s8_first_test_exact.py` |
| S24, S25, S26 | amendments; observation model; decision rule (first test) | `code/bsc_validation/PREREGISTRATION.md`, `PREREG_ADDENDUM_A1_power.md` | — (protocol text) |
| S27 | power; generic model | `data/bsc_validation/02_power.json`, `data/validation_checks/v05_generic.txt` | `scripts/02_power.py`, `code/appC_protocol_tables.py` |
| S28 | sensitivity analyses | `data/s8_first_test_exact.json` (sens), `data/bsc_validation/09_sensitivity.json` (frozen propagator) | `code/s8_first_test_exact.py`, `scripts/09_sensitivity.py`, `code/appC_protocol_tables.py` |
| S29 | level checks and controls | `data/bsc_validation/04_level_controls.json`, `data/validation_checks/v13_c3.txt` | `scripts/04_level_controls.py`, `code/appC_protocol_tables.py` |
| S30 | post-hoc log | `code/bsc_validation/POSTHOC_LOG.md` | — |
| S31 | contact records | `data/bsc_validation2/a2_contact_mobility/s*_<dataset>.json`, `code/checks/contacts/check_summary.json` | `code/bsc_validation2/a2_contact_mobility/src/02_run.py`, `code/checks/contacts/v4_static.py`, `code/v2_numbers.py` |
| S32 | tracer kernel on the first dataset | `data/bsc_validation2/a3_tracer_physics/04_evaluation.json` | `code/bsc_validation2/a3_tracer_physics/src/s04_evaluate.py`, `code/v2_numbers.py` |
| S33 | zone-level test | `data/bsc_validation2/a4_zone_level/matsumoto_primary.json`, `code/checks/zones/v_main.json`, `v_floor.json` | `code/bsc_validation2/a4_zone_level/src/04_matsumoto.py`, `code/checks/zones/v02_main.py`, `v03_floor_power.py` |
| S34 | freeze records of the further analyses | `code/bsc_validation2/*/PREREG_FREEZE.json` and the other freeze records; `data/v2_numbers.json` (`freeze`) | `code/v2_numbers.py` |
| S35 | rule and power, second outbreak test | `data/bsc_validation2/a1_more_outbreaks/03_power.json`, `10_corrected.json` | `scripts/03_power.py`, `10_corrected.py` of that analysis |
| S36 | mixing coefficient by event | `data/bsc_validation2/a1_more_outbreaks/10_corrected.json` | `code/v2_numbers.py` |
| S37 | contact datasets | `data/bsc_validation2/a2_contact_mobility/s1_*.json`, `s3a_*.json` | `code/v2_numbers.py` |
| S38 | rule and power, tracer analysis | `data/bsc_validation2/a3_tracer_physics/03_power.json` | `code/bsc_validation2/a3_tracer_physics/src/s03_power.py` |
| S39 | register of tests | — (summary; outcomes as in `notes/VERIFICATION.md`) | — |
| S40 | tests of the theorems | `data/bsc_theory/verify_theorems.json` | `code/bsc_theory/scripts/verify_theorems.py` |
| S41 | tests of the simulator | `data/bsc_sim/01_validation.json`, `data/appA_numerics_checks.json` | `experiments/01_validate.py`, `code/appA_numerics_checks.py` |

`experiments/…` stands for `code/bsc_sim/experiments/…` and `scripts/…` (Tables 7 and S20 to S29) for
`code/bsc_validation/scripts/…`.

## 5. What was re-run in the released layout

All commands of Section 1 were run in the working tree that the repository is assembled from, and the figure
scripts and table generators of Section 1 were re-run for the article. When the repository was
assembled, identifiers, comments and printed labels of the research code were reworded and a few outputs that
the article does not use were removed, from the code and from the stored results alike (keys of JSON files,
columns of tables, lines of logs), so that each stored file is what its script, as released, writes; the
protocol records are redacted copies (`code/REDACTIONS.md`). The
commands of Section 1, `code/bsc_theory/scripts/verify_theorems.py`, `code/bsc_sim/experiments/20_figures.py`,
`code/bsc_validation/scripts/20_figures.py`, `code/bsc_outbreaks/scripts/analyze_dataset.py` and
`code/bsc_sim/experiments/03b_outbreak_probability.py` (a simulation: identical output) were also run in a copy of
the repository layout, with outputs equal to the stored files up to time stamps. The other simulation experiments
and the verification scripts were run in the author's working tree, where the research code lived under
`research/bsc_<name>/` with the same internal structure, and were not repeated in the repository layout.

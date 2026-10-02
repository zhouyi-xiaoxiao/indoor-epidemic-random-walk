# Re-implementations of the analyses of the second set (re-checks)

Each analysis of `code/bsc_validation2/` was implemented a second time, with separately written code within the
same project and on the same machine, before its result was accepted. This is a re-check, not a third-party
verification. The scripts are kept as they were run; they are written to be run from inside their own
directory, and their stored outputs are next to them. What each check examined and found is recorded in
`notes/VERIFICATION.md`, Sections 4.6 to 4.9; the article reports the outcome after these checks.

| Directory | Analysis | Content |
|---|---|---|
| `outbreaks2/` | second outbreak test | `vlib.py` (own geometry, statistics and a third propagator, by uniformisation); `t0_kernel.py` (the three propagators compared); `t1_tables.py`, `t2_analysis.py` (the registered rule and its power, other near zones, an exponential kernel) -> `t2_results.json`; `t3_extra.py`; `t4_transfer.py` (train kernel on the flights); `t5_bins.py` (other strata); `t6_first_set_floor.py` (effect of the numerical defect on the restaurant calibration of the first test; its effect on the classroom is computed by `code/s8_first_test_exact.py`); `t7_m1.py` (same-site model) |
| `contacts/` | contact records | `vlib.py`, `v1.py` (own loader, statistics, walk simulator and SIR); `v2_rerun.py` (stored runs repeated); `v3_tables.py`; `v4_static.py` (static network of pair rates); `v5_selfpower.py`, `v5c_selfpower.py` (the modified walks fitted to their own synthetic data); `check_summary.json` |
| `tracer/` | tracer measurements | `v01_checks.py`, `v02_repro.py`, `v03_ext_checks.py` (cross-checks made with the analysis itself); `separate/` (own solver and statistics: `a01_kinahan.py … a07_ext.py`, `advlib.py`) |
| `zones/` | zone level | `vcore.py`, `v00_data.py`, `v01_lyon.py`, `v02_main.py` (tiers, guessed shares, other contact metrics), `v03_floor_power.py` (bound on the external hazard, building by floor, power), `v04_misc.py` -> `v_main.json`, `v_floor.json`, `v_power.json` |

Some scripts read source files that are not redistributed (downloaded articles, the raw contact records, the
tracer workbooks); `python code/download_contact_data.py --schools` fetches the contact records and the school
data into the directories of the analyses. The scripts of `contacts/` read the contact records from there. Those
of `zones/` expect their own copies, downloaded separately when the check was made, in a directory `raw/`
next to them: `raw/m.jld2` (the pupil-level file `data_anonymizedstudents.jld2`), `raw/primaryschool.csv.gz`
and `raw/meta.txt` (`primaryschool_metadata.txt`); `v00_data.py` then writes the table `v_students.csv` that
the other scripts read. The scripts of `tracer/separate/` read the workbooks of the cabin measurements
from the directory `data/raw/` of the tracer analysis (see `data/tracer/SOURCES.md` for their DOIs). Intermediate tables of the first directory (`cache/`) are
rebuilt on the first run.

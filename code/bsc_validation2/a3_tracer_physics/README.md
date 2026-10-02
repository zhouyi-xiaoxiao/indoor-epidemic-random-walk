# The airborne kernel fixed by tracer measurements

Section 7.8 of the article and Supplementary Sections S9.3 and S10.4. The mixing coefficient of the airborne term is estimated from
published tracer-gas and tracer-aerosol measurements (a coach, an aircraft cabin, a rail carriage, a
restaurant; a regression for rooms), with nothing fitted to outbreaks, and the kernel is then confronted with
the seven events of the first outbreak test.

| Path | Content |
|---|---|
| `PREREGISTRATION.md`, `PREREG_FREEZE.json` | protocol and freeze record (the record also lists the raw source files, which are not redistributed) |
| `PREREG_ADDENDUM_power.md`, `PREREG_FREEZE_2_before_unblinding.json` | power of the rule and the hashes of kernels and predictions, before any comparison with outbreak counts |
| `PREREG_EXT_flights.md`, `PREREG_EXT_FREEZE*.json` | exploratory extension to four further flights, written after the main results were known |
| `POSTHOC_LOG.md` | every choice made after the freeze |
| `src/p01_ou_seats.py`, `p02_kinahan_tidy.py`, `p03_woodward_digitise.py`, `k01_dump_kinahan.py` | extraction of the tracer tables from the source files (need the source files) |
| `src/a3lib.py`, `s01_kernels.py … s05_descriptive.py` | events, estimators and kernels; predictions; power; evaluation; descriptive tables |
| `src/e01_ext_predict.py`, `e02_ext_evaluate.py` | the extension |
| `data/derived/` (link to `data/tracer/`) | tracer tables, with `SOURCES.md` |
| `results/` (link to `data/bsc_validation2/a3_tracer_physics/`) | stored results |

## Running

```sh
sh code/bsc_validation2/a3_tracer_physics/run_all.sh
```

re-runs the analysis from the tracer tables (steps `s01` to `s05`, `e01`, `e02`; each less than 10 minutes and
1 GB), overwrites `results/` and then runs the two cross-checks of `code/checks/tracer/`. The figures of the
article are drawn by `code/v2_figures.py` from the stored results. The code of the first outbreak test
(`code/bsc_validation/src`) and the events of the second (`code/bsc_validation2/a1_more_outbreaks/src`) are
imported.

The cross-checks made with the analysis itself (`v01_checks.py`, `v02_repro.py`, `v03_ext_checks.py`) and the
re-implementation of the re-check with separately written code (`separate/`) are in `code/checks/tracer/`.

The protocol records of this directory (protocol, addenda, freeze records, the frozen code files and the post-hoc
log) are redacted copies: passages that do not concern the analyses of this article are omitted (marked `[...]`)
or reworded, identifiers are renamed as in the rest of the repository, and in the post-hoc log the words for the
organisation of the work are replaced as in the other notes ("re-check" for the second pass, "stretch of work"
for an interrupted or resumed working period, the names of the analyses for their labels). Every change is listed, line by line, in
`code/REDACTIONS.md`, with the SHA-256 hashes of the original and of the copy; the hashes of the freeze records
refer to the unredacted originals. The protocol texts call the four analyses of the second set avenues 1 to 4
(A1 further outbreaks, A2 contact records, A3 tracer measurements, A4 zone level) and the first and second sets
of tests rounds 1 and 2; the article does not use these labels.
Where they speak of a "verifier", they mean the re-check of the analysis (separately written code, within the same project
and on the same machine, `code/checks/`), not an independent party, and "the analyst" means the first analysis;
both passes were made by the author. "The memo" means the working notes of the analysis, which are not
distributed; their results are reported in the article.
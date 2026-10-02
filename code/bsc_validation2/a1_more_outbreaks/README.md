# Second outbreak test: eleven further records (ten outbreaks) with seat or position data

Sections 7.3 and 7.5 of the article and Supplementary Sections S9.1 and S10.2. Six outbreaks caused by pathogens other than SARS-CoV-2 calibrate
one mixing coefficient (and one constant near/far risk ratio) per setting class; five SARS-CoV-2 events are
held out.

| Path | Content |
|---|---|
| `PREREGISTRATION.md`, `PREREG_ADDENDUM_power.md`, `PREREG_FREEZE.json` | protocol, calibration result and power (written before any held-out stratum count was scored), freeze record with SHA-256 hashes; redacted copies, see `code/REDACTIONS.md` |
| `POSTHOC_LOG.md` | every choice made after the freeze; entry P5 is the numerical correction of the room kernel |
| `src/a1/digitized.py` | seat maps transcribed from the published figures (sources: `data/outbreaks2/SOURCES.md`) |
| `src/a1/events.py`, `lattice.py`, `engine.py`, `engine_fix.py`, `stats.py`, `analysis.py`, `pipeline.py` | geometry and strata of each event; lattice propagator (copied from `code/bsc_validation`); exposures; corrected propagator (image sum); conditional distributions; calibration and decision rule |
| `scripts/00_*.py … 10_corrected.py` | the analysis, in order |
| `out/` (link to `data/bsc_validation2/a1_more_outbreaks/`) | stored results: `04_holdout.json` (as first computed), `10_corrected.json` (corrected kernel; the values of the article), `05_secondary.json`, `06_variants.json`, `07_descriptive.json`, `09_numerics_check.json`, power and logs |

The protocol is called a pre-registration in its file name. It was written with every outbreak count known
and was not lodged with a registry; the article calls it pre-specified.

## Running

```sh
cd code/bsc_validation2/a1_more_outbreaks
python scripts/00_check_data.py        # transcribed seat maps against the totals printed in the sources
python scripts/00_selftest.py
sh scripts/run_tables.sh F1 F2 F3 F4 F5 F6 F7 F8 W1 W2 P1   # exposure tables on the grid of coefficients -> out/tab_*_primary.npz
for e in F1 F2 F3 F4 F5 F6 F7 F8 W1 W2 P1; do python scripts/01b_tables_fix.py $e; done   # corrected kernel -> out/tab_*_nfix.npz
python scripts/02_calibrate.py         # calibration events only
python scripts/03_power.py             # power of the rule; reads no held-out stratum count
python scripts/04_holdout.py           # first script that scores the held-out events
python scripts/05_secondary.py
sh scripts/run_variants.sh && python scripts/06_variants.py
python scripts/07_descriptive.py       # near/far risk ratios
python scripts/09_numerics_check.py
python scripts/10_corrected.py         # the whole rule again with the corrected kernel
```

The exposure tables (`out/tab_*.npz`, about 80 MB) are not stored in the repository; they are rebuilt by the
third and fourth commands (a few seconds to a few minutes per event). The figures of the article are drawn by
`code/v2_figures.py` from the stored results. Each script runs in less than 20 minutes and 4 GB.

The re-implementation of the re-check (separately written code, same project) is in `code/checks/outbreaks2/`.

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
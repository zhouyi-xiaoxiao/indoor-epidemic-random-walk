# Zone-level structure

Section 7.9 of the article and Supplementary Sections S9.4 and S10.5. Mixing shares between a pupil's own class, grade and school,
measured as contact durations in a primary school in Lyon, are used to predict the influenza epidemic of
2014/15 in the 29 primary schools of Matsumoto; a second, low-powered check uses the cases by floor of the
Seoul call-centre building.

| Path | Content |
|---|---|
| `PREREGISTRATION.md`, `PREREG_FREEZE.json` | protocol (with the power analysis on synthetic epidemics) and freeze record |
| `POSTHOC_LOG.md` | every choice made after the freeze |
| `src/00_convert_jld2.py`, `01_lyon_mixing.py` | conversion of the pupil-level data to a table; mixing shares from the contact records (both need the source files; their outputs are in `data/schools/`) |
| `src/zonecore.py`, `analysis.py` | chain-binomial likelihood of the five models, cross-validation, bootstrap, forward simulation |
| `src/02_power.py`, `03_power_summary.py` | power analysis |
| `src/04_matsumoto.py`, `05_seoul_floors.py`, `06_posthoc_descriptive.py` | the school test; the building by floor; descriptive additions |
| `data/derived/` (link to `data/schools/`) | pupil-level table (MIT licence of its authors), mixing shares, `SOURCES.md` |
| `results/` (link to `data/bsc_validation2/a4_zone_level/`), `logs/` | stored results and logs |

## Running

```sh
cd code/bsc_validation2/a4_zone_level
python src/04_matsumoto.py             # about 7 minutes; writes results/matsumoto_primary.json
python src/05_seoul_floors.py
python src/06_posthoc_descriptive.py
```

The scripts are run from this directory (their paths are relative to it). `python
code/download_contact_data.py --schools` fetches the two source files for `src/00_convert_jld2.py` (which
needs the package `h5py`) and `src/01_lyon_mixing.py`.

The re-implementation of the re-check (separately written code, same project) is in `code/checks/zones/`.

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
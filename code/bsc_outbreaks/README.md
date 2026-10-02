# First outbreak dataset

Supplementary Section S6.1 of the article. Twenty-seven records of documented indoor outbreaks and exposure
cohorts from 31 published sources, with counts, geometry, duration, ventilation, a quality grade and the mapping
to the lattice model.

| Path | Content |
|---|---|
| `outbreaks.csv`, `outbreaks.json` (links to `data/bsc_outbreaks/`) | the 27 records (51 fields) |
| `sources.md` | the 31 sources, with DOI, registry check and how each was read |
| `data/` (link to `data/bsc_outbreaks/tables/`) | auxiliary tables: distance and zone strata, risk ratios, the train attack-rate matrix, the restaurant tables, the digitised seat map and onset curve of the Seoul call centre, the record-to-model mapping, and `analysis_results.json` (model-free statistics) |
| `scripts/records.py`, `scripts/mapping.py` | the transcribed numbers (the only place where they are typed) and the mapping to the model |
| `scripts/build_dataset.py`, `scripts/analyze_dataset.py` | build the dataset files; model-free statistics and the permutation test |
| `PROTOCOL_PROPOSAL.md` | the validation protocol proposed with the dataset (redacted copy), on which the first outbreak test is based |
| `verify/` | the re-computation of the re-check (separately written code, same project, same machine) and its results |

## Running

```sh
cd code/bsc_outbreaks
python scripts/build_dataset.py       # scripts/records.py + mapping.py -> outbreaks.*, sources.md, data/*.csv
python scripts/analyze_dataset.py     # model-free statistics (seed 20261001) -> data/analysis_results.json
```

`run_all.sh` also redigitises the seat map of the call centre, which needs Fig. 2 of Park et al. (2020); that image
is third-party material and is not distributed, and `run_all.sh` writes its logs to `notes/`, which must be
created first. The scripts `verify_dois.py`, `fetch_epmc.py` and `fetch_bioc.py` need network access. See
`reproduce.md`, Section 2.3.

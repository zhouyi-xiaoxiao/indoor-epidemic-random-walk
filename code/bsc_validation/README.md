# First outbreak test

Sections 7.3 and 7.4 of the article and Supplementary Sections S7 and S8. Two spatial models, M1 (people walk on
the lattice, transmission only between people on the same site) and M2 (with a non-local airborne term: quanta
emitted by the index case walk on the lattice and are removed by ventilation), are compared with a well-mixed
Wells–Riley model M0. Calibration on the high-speed-train cohort (T4) and the Guangzhou restaurant (R1); held out:
the Zhejiang bus (T1), the Hunan coach (T2), flight VN54 (T5), a classroom (C1) and the Seoul call centre (O1), with
level checks (T3, C3), further level events (H1, O2, R2, C5) and negative controls (S1, C3, the premium economy of
T5, T7).

| Path | Content |
|---|---|
| `PREREGISTRATION.md`, `PREREG_ADDENDUM_A1_power.md`, `PREREG_FREEZE.json` | protocol, power addendum and freeze record with SHA-256 hashes (redacted copies, see below) |
| `POSTHOC_LOG.md` | every choice made after a model had been run against data |
| `src/valmod/` | event geometries, lattice kernels, statistics, the interface to the exact simulator (`exact.py`), predictions and the decision rule |
| `scripts/00_selftest.py … 30_tables.py` | the analysis, in order (`run_all.sh`) |
| `data/` (link to `data/bsc_validation/`) | stored results, including `tables.md` |
| `verify/` | the re-implementation of the re-check (separately written code, same project, same machine) |

The protocol was written with the outbreak counts known and was not lodged with a registry; the article calls it
pre-specified, although its file is named `PREREGISTRATION.md`. It is based on the protocol proposed with the
dataset, released as `../bsc_outbreaks/PROTOCOL_PROPOSAL.md`; Supplementary Table S24 lists every amendment.

## Running

```sh
cd code/bsc_validation
sh run_all.sh        # fixed seeds; writes only inside this directory; temporary files go to .tmp/
```

`reproduce.md`, Section 2.4, lists the scripts one by one and says which of their outputs the article uses (the
article prints the values obtained with the round-off-safe propagator of `code/s8_first_test_exact.py`).

## Protocol records

The protocol, the power addendum, the freeze record, the frozen code files and the post-hoc log are redacted
copies: passages that do not concern the analyses of this article are omitted (marked `[...]`) or reworded, and
identifiers are renamed as in the rest of the repository. Every change is listed, line by line, in
`code/REDACTIONS.md`, with the SHA-256 hashes of the original and of the copy; the hashes of `PREREG_FREEZE.json`
refer to the unredacted originals. Where these records speak of a "verifier", they
mean the re-check of this analysis (separately written code, within the same project and on the same machine,
`verify/`), not an independent party; "the analyst" means the first analysis; both passes were made by the author.
"The memo" means the working notes of the analysis, which are not distributed; their results are reported in the
article. The brute-force check that `POSTHOC_LOG.md` calls `08_independent_check.py` is
`scripts/08_bruteforce_check.py` here (results: `data/bsc_validation/08_bruteforce_check.json`).

#!/usr/bin/env python
"""download_contact_data.py -- fetch the third-party files that are not redistributed (needs network access).

1. The ten SocioPatterns files listed in data/contacts/datasets.csv are downloaded from sociopatterns.org into
   the raw-data directory of the contact-record analysis, and their SHA-256 is compared with the value recorded
   when the protocol of that analysis was frozen.  A file that is already present is only checked.
2. With --schools, the two raw files of the zone-level analysis are fetched as well: the primary-school
   contact file (the same file as in 1) and the pupil-level data of Endo et al. (2021) at the commit used.
   The tabular copy of the latter and the mixing shares derived from the former are already in data/schools/,
   so this step is needed only to re-derive them.

    python download_contact_data.py [--schools]
"""
from __future__ import annotations

import csv
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
LIST = HERE.parent / "data" / "contacts" / "datasets.csv"
RAW = ROOT / "code/bsc_validation2/a2_contact_mobility/data/raw"
ZONE_RAW = ROOT / "code/bsc_validation2/a4_zone_level/data/raw"
UA = {"User-Agent": "indoor-epidemic-random-walk/1.1 (mailto:zhouyixiaoxiao@gmail.com)"}
ENDO = ("https://raw.githubusercontent.com/akira-endo/schooldynamics_FluMatsumoto14-15/"
        "658abfecc70440b0b849e179da457048082d4e72/data/anonymizedstudents.jld2",
        "55b454bfde0be0ab1e92830a495362ddf138b879d657fade944dd5195ff8496f")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fetch(url: str, dst: Path, expected: str) -> bool:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists():
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r, open(dst, "wb") as f:
            shutil.copyfileobj(r, f)
    ok = sha256(dst) == expected
    print(("ok      " if ok else "MISMATCH") + f"  {dst.name}")
    return ok


if __name__ == "__main__":
    good = True
    with open(LIST, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            good &= fetch(row["url"], RAW / row["file"], row["sha256"])
    if "--schools" in sys.argv:
        for name in ("primaryschool.csv.gz", "primaryschool_metadata.txt"):
            src = RAW / name
            dst = ZONE_RAW / "sociopatterns" / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            if src.exists() and not dst.exists():
                shutil.copy2(src, dst)
        good &= fetch(ENDO[0], ZONE_RAW / "matsumoto" / "data_anonymizedstudents.jld2", ENDO[1])
    sys.exit(0 if good else 1)

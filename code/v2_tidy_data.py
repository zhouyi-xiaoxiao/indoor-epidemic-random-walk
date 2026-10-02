#!/usr/bin/env python
"""v2_tidy_data.py -- tidy tables of the data of the second set of tests.

Writes, from the transcriptions and derived files of the four analyses (nothing is typed here):

  data/outbreaks2/events.csv    the eleven outbreaks: setting, risk set, near/far counts, role, source
  data/outbreaks2/strata.csv    persons at risk and cases by distance stratum, per event
  data/outbreaks2/seats.csv     seat-map events: one row per person at risk and per source
                                            (lattice column and row, row distance to the nearest source, stratum, case)
  data/tracer/*.csv             tracer tables (coach and minibus seats, cabin sensor counts,
                                            digitised carriage profile), copied from the tracer analysis
  data/schools/matsumoto_students.csv, lyon_mixing.json, LICENSE_matsumoto.txt
  data/contacts/datasets.csv    the contact datasets: file, URL, SHA-256, licence (files not copied)

The sources files (reference and table or figure of origin of every number) are data/outbreaks2/SOURCES.md,
data/tracer/SOURCES.md, data/schools/SOURCES.md and data/contacts/SOURCES.md; they are written by hand and
kept next to the tables.

Run:  python code/v2_tidy_data.py
"""
from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
OUT = HERE.parent / "data"
A1 = ROOT / "code/bsc_validation2/a1_more_outbreaks"
A2 = ROOT / "code/bsc_validation2/a2_contact_mobility"
A3 = ROOT / "code/bsc_validation2/a3_tracer_physics"
A4 = ROOT / "code/bsc_validation2/a4_zone_level"
sys.path.insert(0, str(A1 / "src"))
from a1 import events as E  # noqa: E402

REF = {  # event: (first author and year, DOI, where the numbers come from)
    "F1": ("Olsen et al. 2003", "10.1056/NEJMoa031349", "seat map as redrawn in Hertzberg and Weiss 2016, Fig. 2 (10.1016/j.aogh.2016.06.003)"),
    "F2": ("Kenyon et al. 1996", "10.1056/NEJM199604113341501", "abstract: 4/13 within two rows, 2/55 in the rest of the section"),
    "F3": ("Baker et al. 2010", "10.1136/bmj.c2424", "seat-map figure of the rear section"),
    "F4": ("Young et al. 2014", "10.1111/irv.12181", "Fig. 1 (seat map)"),
    "F5": ("Toyokawa et al. 2022", "10.1111/irv.12913", "Fig. 2 (seat map)"),
    "F6": ("Speake et al. 2020", "10.3201/eid2612.203910", "Fig. 4 (seat map of the mid cabin)"),
    "F7": ("Swadi et al. 2021", "10.3201/eid2703.204714", "text and Fig. 3 (seats within rows 23-30)"),
    "F8": ("Hoehl et al. 2020", "10.1001/jamanetworkopen.2020.18044", "Fig. 1 (seat map)"),
    "W1": ("Wong et al. 2004", "10.3201/eid1002.030452", "text and table (3/3, 4/8, 0/8 by assigned bed); Fig. 4 (floor plan)"),
    "W2": ("Yu et al. 2005", "10.1086/428735", "attack rates by bay (13/20, 11/21, 6/33); Fig. 1"),
    "P1": ("Guenther et al. 2020", "10.15252/emmm.202013296", "Appendix Table S2 (cumulative counts by 1-m distance)"),
}
SETTING = {"F1": "Boeing 737-300", "F2": "wide-body aircraft", "F3": "Boeing 747-400, rear section", "F4": "Boeing 767",
           "F5": "Boeing 737-800", "F6": "Airbus A330-200, economy mid cabin", "F7": "Boeing 777-300ER",
           "F8": "Boeing 737-900", "W1": "hospital ward 8A, bedside session", "W2": "hospital ward 8A, inpatients",
           "P1": "beef-processing line, 32 m x 8.5 m hall"}
STRATA = {"air4": ["0-1 rows from the nearest source", "2 rows", "3-5 rows", "6 or more rows"],
          "F7": ["0-1 rows from the nearest source", "2 rows", "3 or more rows"],
          "F2": ["within 2 rows of the index case", "rest of the cabin section"],
          "W1": ["within 1 m of the index bed", "same cubicle", "elsewhere in the ward"],
          "W2": ["same bay as the index patient", "adjacent bay", "distant bays"],
          "P1": ["0-4 m from the index case", "4-8 m", "8-12 m", "more than 12 m"]}
DUR = dict(E.DUR) | {"W1": 40.0 / 60.0, "W2": None, "P1": 24.0}
des = json.loads((A1 / "out/07_descriptive.json").read_text())["events"]


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"{path.relative_to(OUT.parent)}: {len(rows)} rows")


# ---------------------------------------------------------------------------------------------- outbreaks
ev_rows, st_rows, seat_rows = [], [], []
for ev in E.EVENTS:
    meta, draws = E.build(ev)
    d = des[ev]
    role = "hold-out" if ev in sum(E.HOLD.values(), []) else "calibration"
    grade = "B" if ev == "W2" else "A"
    ev_rows.append([ev, E.LABEL[ev], E.PATH[ev], SETTING[ev], E.CLS[ev], DUR[ev] if DUR[ev] is not None else "",
                    sum(meta["sizes"]), meta["K"], d["near"][0], d["near"][1], d["far"][0], d["far"][1],
                    f"{d['RR']:.3f}", f"{d['lo']:.3f}", f"{d['hi']:.3f}", d["cc"], grade, role, REF[ev][0], REF[ev][1], REF[ev][2]])
    names = STRATA.get(ev, STRATA["air4"])
    for b, (n, k) in enumerate(zip(meta["sizes"], meta["obs"])):
        st_rows.append([ev, b, names[b], n, k, int(b in meta["near_bins"])])
    if "y" in meta:                                   # seat-map events
        lat = draws[0]["lat"]
        for s in draws[0]["src"]:
            seat_rows.append([ev, "source", int(s % lat.nx), int(s // lat.nx), 0, "", ""])
        for s, dist, b, y in zip(draws[0]["rec"], meta["dist"], draws[0]["bins"], meta["y"]):
            seat_rows.append([ev, "at risk", int(s % lat.nx), int(s // lat.nx), int(dist), int(b), int(bool(y))])
write(OUT / "outbreaks2/events.csv",
      ["event", "label", "pathogen", "setting", "class", "duration_h", "n_at_risk", "n_cases", "near_cases", "near_n",
       "far_cases", "far_n", "risk_ratio", "rr_lo95", "rr_hi95", "continuity_correction", "grade", "role", "source", "doi",
       "origin_of_numbers"], ev_rows)
write(OUT / "outbreaks2/strata.csv", ["event", "stratum", "definition", "n_at_risk", "n_cases", "in_near_zone"], st_rows)
write(OUT / "outbreaks2/seats.csv",
      ["event", "role", "lattice_column", "lattice_row", "rows_to_nearest_source", "stratum", "case"], seat_rows)
tot = {ev: (sum(r[3] for r in st_rows if r[0] == ev), sum(r[4] for r in st_rows if r[0] == ev)) for ev in E.EVENTS}
assert all(tot[r[0]] == (r[6], r[7]) for r in ev_rows)
for ev in ("F1", "F3", "F4", "F5", "F6", "F8"):
    n = sum(1 for r in seat_rows if r[0] == ev and r[1] == "at risk")
    k = sum(r[6] for r in seat_rows if r[0] == ev and r[1] == "at risk")
    assert (n, k) == tot[ev], (ev, n, k, tot[ev])

# ---------------------------------------------------------------------------------------------- tracer
def copy(src: Path, dst: Path):
    """Copy unless the source is missing (files that are not redistributed) or is the destination itself."""
    if src.exists() and src.resolve() != dst.resolve():
        shutil.copy2(src, dst)


(OUT / "tracer").mkdir(exist_ok=True)
for f in sorted((A3 / "data/derived").glob("*.csv")):
    copy(f, OUT / "tracer" / f.name)
print("tracer:", len(list((OUT / "tracer").glob("*.csv"))), "tables")

# ---------------------------------------------------------------------------------------------- schools
(OUT / "schools").mkdir(exist_ok=True)
copy(A4 / "data/derived/matsumoto_students.csv", OUT / "schools/matsumoto_students.csv")
copy(A4 / "data/derived/lyon_mixing.json", OUT / "schools/lyon_mixing.json")
copy(A4 / "data/raw/matsumoto/LICENSE", OUT / "schools/LICENSE_matsumoto.txt")
print("schools: matsumoto_students.csv, lyon_mixing.json, LICENSE_matsumoto.txt")

# ---------------------------------------------------------------------------------------------- contacts
USED = {"workplace_InVS_tij.dat.zip": "office building, 2013", "workplace_InVS_metadata.txt": "office building, 2013 (departments)",
        "workplace_InVS15_tij.dat.gz": "office building, 2015", "workplace_InVS15_metadata.txt": "office building, 2015 (departments)",
        "primaryschool.csv.gz": "primary school", "primaryschool_metadata.txt": "primary school (classes)",
        "hospital_lyon_contacts.dat.gz": "hospital ward", "HighSchool2013_proximity_net.csv.gz": "high school",
        "HighSchool2013_metadata.txt": "high school (classes)", "SFHH_tij.dat.gz": "scientific conference"}
rows = []
for line in (A2 / "data/raw_SHA256.txt").read_text().split("\n"):
    if not line.strip():
        continue
    h, name = line.split()
    if name in USED:
        rows.append([name, USED[name], "http://www.sociopatterns.org/assets/data/" + name, h,
                     "see http://www.sociopatterns.org/datasets/ (Creative Commons; not redistributed here)"])
write(OUT / "contacts/datasets.csv", ["file", "dataset", "url", "sha256", "licence"], rows)
assert len(rows) == len(USED)

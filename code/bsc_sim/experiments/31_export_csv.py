#!/usr/bin/env python
"""Flat CSV copies of the JSON result files (one row per configuration)."""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FILES = ["02_finite_size", "02_mobility", "03_scenes", "03b_outbreak_probability",
         "04b_first_generation", "05_heterogeneity", "05_heterogeneity_extra",
         "06_obstacles", "07_occupancy", "08_interventions", "09_metro_time",
         "11_dwell_time"]
for name in FILES:
    path = os.path.join(DATA, name + ".json")
    if not os.path.exists(path):
        continue
    d = json.load(open(path))
    rows = []
    for k, v in d.items():
        if isinstance(v, dict) and all(not isinstance(x, (dict, list)) for x in v.values()):
            rows.append(dict(key=k, **v))
    if not rows:
        continue
    cols = ["key"] + sorted({c for r in rows for c in r} - {"key"})
    with open(os.path.join(DATA, name + ".csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(name, len(rows), "rows")
# spatial correlations: one row per scene, mobility and generation
sp = json.load(open(os.path.join(DATA, "04_spatial.json")))
with open(os.path.join(DATA, "04_spatial.csv"), "w", newline="") as f:
    w = None
    for k, r in sp.items():
        for g in (1, 2, 3):
            row = dict(scene=r["scene"], D0=r["D0"], generation=g, **r[f"gen{g}"])
            if w is None:
                w = csv.DictWriter(f, fieldnames=list(row))
                w.writeheader()
            w.writerow(row)
# benchmark: one row per (D0, task)
bm = json.load(open(os.path.join(DATA, "10_benchmark.json")))
rows = []
for k, R in bm.items():
    for t, v in R.items():
        if isinstance(v, dict):
            rows.append(dict(config=k, task=t, **v))
cols = ["config", "task"] + sorted({c for r in rows for c in r} - {"config", "task"})
with open(os.path.join(DATA, "10_benchmark.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
print("done")

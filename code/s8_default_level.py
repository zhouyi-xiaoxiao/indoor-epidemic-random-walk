"""s8_default_level.py -- how many infections the default parameters allow in a short seated event.

Under the frequency-dependent pair hazard of Section 2.4 an infective at cell x infects at expected rate
beta*q_x when the other walkers are independent and uniformly distributed, which they are at every time because
the walks are stationary and do not depend on the epidemic.  Infections form chains of transmission between
distinct people.  Summing, over such chains, the expected number of ordered contact events along each chain (the
walkers move independently and remain at stationarity, and the pair hazard carries the factor 1/(N-1), so a chain
of k transmissions within time T through a given ordered choice of k distinct further people has expected weight
at most (beta*q_max*T/(N-1))^k / k!, and there are (N-1)_k such choices), the expected number of infections caused in an event of duration T by one index case is at most
sum_k (N-1)_k/(N-1)^k (beta*q_max*T)^k/k! <= exp(beta*q_max*T) - 1, where q_max is the largest infection
efficiency of the scenes.  (Counting every infection attempt instead, including attempts on a person's own
infector, would not give a bound: an infective is not independent of the position of its infector.)  Markov's inequality then bounds the probability of at least n infections by
that mean divided by n.  The script evaluates the bound at the default beta = 0.5/day for the two events with the
highest attack rates of short seated exposure in the first dataset (the Zhejiang bus and the Skagit choir) and
gives the beta that would make the bound equal to the observed count.

Inputs (read only):
    data/bsc_outbreaks/outbreaks.csv           (n_exposed, n_infected, n_infected_alt, exposure_duration_h)
    code/bsc_sim/scenes/{office,supermarket,classroom,metro}.json   (zone values of q)
Output: ../data/s8_default_level.json

    python s8_default_level.py
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
OB = ROOT / "data/bsc_outbreaks/outbreaks.csv"
SC = ROOT / "code/bsc_sim/scenes"
OUT = HERE.parent / "data" / "s8_default_level.json"

BETA = 0.5  # per day, default of Section 2.6

q_scene = {}
for name in ("office", "supermarket", "classroom", "metro"):
    z = json.loads((SC / f"{name}.json").read_text())["zones"]
    q_scene[name] = max(v["q"] for v in z.values())
q_max = max(q_scene.values())

rows = {r["id"]: r for r in csv.DictReader(OB.open(encoding="utf-8"))}
events = {}
for rid, label in (("T1_zhejiang_bus", "Zhejiang bus"), ("H1_skagit_choir", "Skagit choir")):
    r = rows[rid]
    T_day = float(r["exposure_duration_h"]) / 24.0
    mean_bound = math.expm1(BETA * q_max * T_day)
    counts = {"all cases": int(r["n_infected"])}
    if r["n_infected_alt"]:
        counts["confirmed cases only"] = int(r["n_infected_alt"])
    e = dict(id=rid, n_exposed=int(r["n_exposed"]), exposure_h=float(r["exposure_duration_h"]),
             expected_infections_bound=mean_bound, attack_rate_bound=mean_bound / int(r["n_exposed"]),
             cases={})
    for k, n in counts.items():
        # beta at which the bound equals the observed count: exp(beta*q_max*T) - 1 = n
        beta_req = math.log1p(n) / (q_max * T_day)
        e["cases"][k] = dict(n_infected=n, observed_attack_rate=n / int(r["n_exposed"]),
                             markov_bound_P_at_least_n=mean_bound / n, beta_needed_per_day=beta_req,
                             beta_needed_over_default=beta_req / BETA)
    events[label] = e

res = dict(beta_per_day=BETA, q_max_by_scene=q_scene, q_max=q_max,
           bound="E[infections in T] <= exp(beta*q_max*T) - 1; P(N >= n) <= E/n", events=events)
OUT.write_text(json.dumps(res, indent=1) + "\n")
print(json.dumps(res, indent=1))

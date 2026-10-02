"""s8_outbreaks_checks.py -- recompute the model-free statistics quoted in Section 7 and re-run one script
of the re-check whose output had not been stored.

Part A (a few seconds; always run).  From the outbreak dataset only:
  * near/far risk ratios of the four single-event shape tests, Cochran's Q, the fixed-effect pooled ratio and the
    DerSimonian-Laird random-effects pooled ratio (the estimator used for the second dataset);
  * the Seoul call centre: wing contrast (with the two extreme placements of the 10 unmapped cases) and the
    join-count statistic of the north wing (exact permutation moments, and the stored permutation sample of
    data/bsc_outbreaks/tables/callcentre_joincount_null.npy).
  Output: ../data/s8_outbreaks_checks.json

Part B (--check; several minutes).  Runs, unchanged and in a scratch copy inside article/data, the
re-check's script code/bsc_validation/verify/v09_sens.py, which gives
  * the emission rate (quanta/h) that model M2 needs for each hold-out event (5 %, 50 %, 95 % quantiles);
  * the posterior interval of D_air after tempering the train likelihood by its over-dispersion (deviance/d.f.);
  * the hold-out p-values under that tempered posterior and under ventilation priors scaled by 1/3 and 3.
  Output: ../data/s8_outbreaks_v09_sens.txt

Nothing is written to the research directories.

    python s8_outbreaks_checks.py [--check]
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                      # repository root
OB = ROOT / "code/bsc_outbreaks"
VER = ROOT / "code/bsc_validation/verify"
OUT = HERE.parent / "data"
OUT.mkdir(exist_ok=True)


def clopper_pearson(k: int, n: int, level: float = 0.95):
    a = 1 - level
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)


def risk_ratio(k1, n1, k0, n0):
    """Risk ratio with the log-normal (Katz) 95 % interval."""
    rr = (k1 / n1) / (k0 / n0)
    se = math.sqrt(1 / k1 - 1 / n1 + 1 / k0 - 1 / n0)
    return rr, rr * math.exp(-1.96 * se), rr * math.exp(1.96 * se), math.log(rr), se


res: dict = {}

# ---------------------------------------------------------------- A1. risk ratios and heterogeneity
rrf = pd.read_csv(OB / "data/near_far_risk_ratios.csv")
pick = {  # outbreak_id, contrast  ->  short name
    ("T5_vn54_flight", "within ~2 m vs beyond, same cabin"): "flight VN54",
    ("C1_marin_classroom", "rows nearest the index vs farther rows"): "Marin classroom",
    ("T1_zhejiang_bus", "within ~2 m (3 rows) vs beyond"): "Zhejiang bus",
    ("T2_hunan_coach", "index half vs other half of coach"): "Hunan coach",
}
rows, lg, se = [], [], []
for (oid, contrast), name in pick.items():
    r = rrf[(rrf.outbreak_id == oid) & (rrf.contrast == contrast)].iloc[0]
    rr, lo, hi, l, s = risk_ratio(int(r.near_k), int(r.near_n), int(r.far_k), int(r.far_n))
    rows.append(dict(event=name, near=f"{int(r.near_k)}/{int(r.near_n)}", far=f"{int(r.far_k)}/{int(r.far_n)}",
                     rr=rr, ci95=[lo, hi], fisher_p=float(stats.fisher_exact(
                         [[int(r.near_k), int(r.near_n - r.near_k)], [int(r.far_k), int(r.far_n - r.far_k)]])[1])))
    lg.append(l)
    se.append(s)
lg, se = np.array(lg), np.array(se)
w = 1 / se**2
pooled = float((w * lg).sum() / w.sum())
Q = float((w * (lg - pooled) ** 2).sum())
res["risk_ratios"] = rows
res["heterogeneity"] = dict(cochran_Q=Q, dof=len(lg) - 1, p=float(stats.chi2.sf(Q, len(lg) - 1)),
                            pooled_rr=math.exp(pooled),
                            pooled_ci95=[math.exp(pooled - 1.96 / math.sqrt(w.sum())),
                                         math.exp(pooled + 1.96 / math.sqrt(w.sum()))])
# the same four ratios pooled with the DerSimonian-Laird random-effects estimator, the estimator of the pooled
# ratio of the second dataset (code/bsc_validation2/a1_more_outbreaks/scripts/07_descriptive.py), so that
# the two pooled ratios can be compared on one estimator
dof = len(lg) - 1
tau2 = max(0.0, (Q - dof) / (w.sum() - (w ** 2).sum() / w.sum()))
wr = 1 / (1 / w + tau2)
mr = float((wr * lg).sum() / wr.sum())
sr = 1 / math.sqrt(wr.sum())
res["heterogeneity"].update(fixed_effect="inverse-variance weights", tau2_DL=tau2,
                            I2=max(0.0, (Q - dof) / Q) if Q > 0 else 0.0,
                            pooled_rr_random_DL=math.exp(mr),
                            pooled_ci95_random_DL=[math.exp(mr - 1.96 * sr), math.exp(mr + 1.96 * sr)])

# ---------------------------------------------------------------- A3. Seoul call centre
seats = pd.read_csv(OB / "data/park2020_callcentre_seats_digitized.csv")
nw = seats[seats.region == "north_wing"]
sw = seats[seats.region == "south_wing"]
kn, nn, ks, ns = int(nw.case.sum()), len(nw), int(sw.case.sum()), len(sw)
rr, lo, hi, _, _ = risk_ratio(kn, nn, ks, ns)
rr_s, lo_s, hi_s, _, _ = risk_ratio(kn, nn, ks + 10, ns)        # the 10 unmapped cases all on south-wing desks
rr_n, lo_n, hi_n, _, _ = risk_ratio(kn + 10, nn, ks, ns)        # ... all on north-wing desks
xy = nw[["x_pitch_units", "y_pitch_units"]].to_numpy()
d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
A = ((d > 0) & (d <= 1.25)).astype(float)                        # adjacency: centres within 1.25 desk pitches
c = nw.case.to_numpy().astype(float)
obs = float(c @ A @ c / 2)
n, k = nn, kn
m_pairs = A.sum() / 2
deg = A.sum(1)
p2 = k * (k - 1) / (n * (n - 1))
p3 = p2 * (k - 2) / (n - 2)
p4 = p3 * (k - 3) / (n - 3)
shared = (deg * (deg - 1)).sum() / 2                              # pairs of edges sharing one desk
disjoint = m_pairs * (m_pairs - 1) / 2 - shared
mean = m_pairs * p2
var = m_pairs * p2 + 2 * shared * p3 + 2 * disjoint * p4 - mean**2
null = np.load(OB / "data/callcentre_joincount_null.npy")
res["callcentre"] = dict(
    north=[kn, nn], south=[ks, ns], rr=rr, rr_ci95=[lo, hi],
    fisher_p=float(stats.fisher_exact([[kn, nn - kn], [ks, ns - ks]])[1]),
    rr_if_10_unmapped_cases_in_south=dict(rr=rr_s, ci95=[lo_s, hi_s]),
    rr_if_10_unmapped_cases_in_north=dict(rr=rr_n, ci95=[lo_n, hi_n]),
    case_seats_mapped=int(seats.case.sum()), cases_reported=94,
    joincount=dict(adjacent_pairs=int(m_pairs), observed=int(obs), exact_mean=float(mean), exact_sd=float(math.sqrt(var)),
                   z_exact=float((obs - mean) / math.sqrt(var)),
                   stored_permutations=int(null.size), perm_mean=float(null.mean()), perm_sd=float(null.std(ddof=1)),
                   z_perm=float((obs - null.mean()) / null.std(ddof=1)),
                   p_one_sided_clustering=float((np.sum(null >= obs) + 1) / (null.size + 1))))

(OUT / "s8_outbreaks_checks.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))

# ---------------------------------------------------------------- B. script v09_sens.py of the re-check
if "--check" in sys.argv:
    scratch = OUT / "_s8_scratch"
    if scratch.exists():
        shutil.rmtree(scratch)
    (scratch / "verify").mkdir(parents=True)
    for f in VER.iterdir():
        if f.suffix in (".py", ".npz", ".npy", ".json"):
            shutil.copy(f, scratch / "verify" / f.name)
    r = subprocess.run([sys.executable, "v09_sens.py"], cwd=scratch / "verify", capture_output=True, text=True)
    (OUT / "s8_outbreaks_v09_sens.txt").write_text(
        "# output of code/bsc_validation/verify/v09_sens.py (re-check), re-run by "
        "code/s8_outbreaks_checks.py --check\n" + r.stdout + ("\n[stderr]\n" + r.stderr if r.returncode else ""))
    print(r.stdout[-4000:])
    if r.returncode:
        print(r.stderr[-2000:])
    shutil.rmtree(scratch)

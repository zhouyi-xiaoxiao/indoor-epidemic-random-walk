"""plan_validation_checks.py -- re-run the re-check's cheap checks of the outbreak validation and
store their output inside the article, so that the statements of Section 7 have a source file.

It executes, unchanged and in a scratch copy, the re-check's scripts
    code/bsc_validation/verify/v05_generic.py      generic constant-risk-ratio model against the decision rule
    code/bsc_validation/verify/v10_c1_loo.py       classroom: implied emission over the room posterior; common-D test
    code/bsc_validation/verify/v11_breakdown.py    classroom/flight/bus/coach at the posterior mode
    code/bsc_validation/verify/v12_pooled_variants.py  pooled score without the flight / the classroom
    code/bsc_validation/verify/v13_c3.py           Utah classroom level check for mask factors 0.2-1.0
and recomputes the call-centre probability by quadrature from the re-check's stored grid
(verify/res_o1_grid.json) and the restaurant posterior of the first analysis (data/01_cal_R1_primary.json), with a
bound on the effect of the round-off of the re-check's propagator through the outbreak-size weights.

Nothing is written to the research directories.  Output: ../data/validation_checks/*.txt and
../data/validation_checks/callcentre_quadrature.json.

    python plan_validation_checks.py                    (2-6 minutes)
    python plan_validation_checks.py --quadrature-only  (part 1 only, a second)
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
VER = ROOT / "code/bsc_validation/verify"
VAL = ROOT / "code/bsc_validation"
OUT = HERE.parent / "data" / "validation_checks"
OUT.mkdir(parents=True, exist_ok=True)

# ---- 1. call centre: quadrature of P(z <= z_obs | D) over the restaurant posterior ---------------------
g = json.loads((VER / "res_o1_grid.json").read_text())
grid = np.array(sorted(float(k) for k in g))
Pz = np.array([g[str(k) if str(k) in g else repr(k)][0] for k in grid])
Kf = np.array([g[str(k) if str(k) in g else repr(k)][1] for k in grid])
cal = json.loads((VAL / "data/01_cal_R1_primary.json").read_text())
ll = np.array(cal["M2"]["profile_logL"])
logD = np.array(cal["logD"])
w = np.exp(ll - ll.max())
w /= w.sum()
Pg = np.interp(logD, grid, Pz)
Kg = np.interp(logD, grid, Kf)
res = {"description": "P(join-count z <= observed) for model M2, integrated over the restaurant-calibrated posterior of the "
                      "mixing coefficient (flat prior on log10 D between the grid limits)",
       "grid_log10D": grid.tolist(), "P_z_le_obs_given_D": Pz.tolist(), "fraction_outbreaks_retained": Kf.tolist(),
       "equal_weight": float((w * Pg).sum()), "pooled_weight": float((w * Kg * Pg).sum() / (w * Kg).sum()),
       "prior_upper_limit_sensitivity": {}}
for top in (2, 3, 4):
    m = logD <= top
    ww = w[m] / w[m].sum()
    res["prior_upper_limit_sensitivity"][f"1e{top}"] = dict(
        equal=float((ww * Pg[m]).sum()), pooled=float((ww * Kg[m] * Pg[m]).sum() / (ww * Kg[m]).sum()))
ext = np.concatenate([w, np.full(20, w[-1])])
ext /= ext.sum()
Pe = np.concatenate([Pg, np.full(20, Pg[-1])])
Ke = np.concatenate([Kg, np.full(20, Kg[-1])])
res["prior_upper_limit_sensitivity"]["1e6_flat_extension"] = dict(
    equal=float((ext * Pe).sum()), pooled=float((ext * Ke * Pe).sum() / (ext * Ke).sum()))
# The propagator of the re-check (a mode sum) suffers round-off below D = 0.1 m^2/h.  P(z <= z_obs | D) is 0 at
# every grid value below 5.6 m^2/h, so the defect can act only through the outbreak-size weights (pooled
# weighting).  Bound its effect: halve, or set to zero, the fraction of outbreaks retained at D <= 0.1 m^2/h.
P_zero_below = float(10 ** grid[np.nonzero(Pz)[0][0] - 1]) if Pz.any() else None
res["P_zero_at_all_grid_values_up_to_m2h"] = P_zero_below
res["weight_sensitivity_at_or_below_0.1"] = {}
for name, fac in (("retained_fraction_halved", 0.5), ("retained_fraction_zero", 0.0)):
    Kmod = np.where(logD <= -1.0 + 1e-9, Kg * fac, Kg)
    res["weight_sensitivity_at_or_below_0.1"][name] = float((w * Kmod * Pg).sum() / (w * Kmod).sum())
res["posterior_mass_at_or_below_0.1"] = float(w[logD <= -1.0 + 1e-9].sum())
(OUT / "callcentre_quadrature.json").write_text(json.dumps(res, indent=1))
print(json.dumps({k: v for k, v in res.items() if k not in ("grid_log10D", "P_z_le_obs_given_D", "fraction_outbreaks_retained")}, indent=1))

if "--quadrature-only" in sys.argv:
    sys.exit(0)

# ---- 2. scripts of the re-check, run in a scratch copy so that the research directory is not touched -------
scratch = OUT / "_scratch"
if scratch.exists():
    shutil.rmtree(scratch)
(scratch / "verify").mkdir(parents=True)
(scratch / "data").mkdir()
for f in VER.iterdir():
    if f.suffix in (".py", ".npz", ".npy", ".json"):
        shutil.copy(f, scratch / "verify" / f.name)
for f in ("02_predictive.json", "01_cal_R1_primary.json", "heldout_outcomes.json"):
    shutil.copy(VAL / "data" / f, scratch / "data" / f)
# the re-check's event builders read the outbreak dataset by relative path; keep the same relative layout
for script in ("v05_generic.py", "v12_pooled_variants.py", "v13_c3.py", "v11_breakdown.py", "v10_c1_loo.py"):
    r = subprocess.run([sys.executable, script], cwd=scratch / "verify", capture_output=True, text=True)
    (OUT / (script.replace(".py", "") + ".txt")).write_text(r.stdout + ("\n[stderr]\n" + r.stderr if r.returncode else ""))
    print("=" * 20, script, "exit", r.returncode)
    print(r.stdout[-3000:])
    if r.returncode:
        print(r.stderr[-2000:])
shutil.rmtree(scratch)

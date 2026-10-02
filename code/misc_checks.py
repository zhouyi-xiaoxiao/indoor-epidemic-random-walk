"""misc_checks.py -- small facts quoted in the article, recomputed from the project files (no
simulation; a few seconds).

  callcentre_lattice : the event builder of the call-centre test maps each desk to the nearest lattice site;
                       number of distinct sites and of sites that hold two desks
  mobility_gap       : range of the walking mobility over the four scenes, in cell^2/day and in orders of
                       magnitude above D0 = 1
  miller             : implied emission of the choir at the breathing rate used by Miller et al. (2021)
  jerusalem          : enrolled pupils by wing (grade denominators printed by the source) against the number tested
  flight_label       : the source of the flight record prints the near stratum as an underlined "<", i.e. "<= 2 seats away"
  coach              : the three denominators of the Hunan coach

    python misc_checks.py
Output: ../data/misc_checks.json
"""
import csv
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
DATA = HERE.parent / "data"
OUT = DATA / "misc_checks.json"
out = json.loads(OUT.read_text()) if OUT.exists() else {}

# ---- call centre: desks -> lattice sites (same arithmetic as code/bsc_validation/src/valmod/events.py, build_O1)
sys.path.insert(0, str(ROOT / "code" / "bsc_validation" / "src"))
import numpy as np                                    # noqa: E402
from valmod import events as ev                       # noqa: E402

xy, _ = ev.o1_desks()
x0, y0 = xy[:, 0].min() - 1.0, xy[:, 1].min() - 1.0
ix, iy = np.rint(xy[:, 0] - x0).astype(int), np.rint(xy[:, 1] - y0).astype(int)
occ = Counter(zip(ix.tolist(), iy.tolist()))
out["callcentre_lattice"] = dict(desks=int(len(xy)), distinct_sites=len(occ),
                                 sites_with_two_desks=sum(1 for v in occ.values() if v == 2),
                                 max_desks_per_site=max(occ.values()))

# ---- mobility gap
sc = json.loads((DATA / "s2_model_scenes.json").read_text())["scenes"]
lo = min(v["walking_5pct_D0_cell2_per_day"] for v in sc.values())
hi = max(v["walking_D0_cell2_per_day"] for v in sc.values())
out["mobility_gap"] = dict(walking_5pct_min=lo, walking_max=hi, log10_min=math.log10(lo), log10_max=math.log10(hi),
                           per_scene={k: [v["walking_5pct_D0_cell2_per_day"], v["walking_D0_cell2_per_day"]]
                                      for k, v in sc.items()})

# ---- choir emission at the breathing rate of Miller et al. (1.0 m^3/h instead of 0.5)
lev = json.loads((ROOT / "data/bsc_validation/04_level_controls.json").read_text())
E = lev["implied"]["H1"]["E_implied"]
out["miller"] = dict(E_implied_at_0p5_m3_per_h=E, E_at_1p0_m3_per_h=E * 0.5 / 1.0,
                     miller_reported="970 +- 390 quanta per hour (accepted manuscript, Table 1)")

# ---- Jerusalem school: grade denominators
rows = {r["id"]: r for r in csv.DictReader(open(ROOT / "data/bsc_outbreaks/outbreaks.csv", encoding="utf-8"))}
txt = rows["C2_jerusalem_highschool"]["attack_rate_reported"]
g = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in re.finditer(r"(\d+)th (\d+)/(\d+)", txt)}
out["jerusalem"] = dict(grades=g, junior_wing=[sum(g[k][0] for k in ("7", "8", "9")), sum(g[k][1] for k in ("7", "8", "9"))],
                        senior_wing=[sum(g[k][0] for k in ("10", "11", "12")), sum(g[k][1] for k in ("10", "11", "12"))],
                        tested=int(rows["C2_jerusalem_highschool"]["n_exposed"]))
assert out["jerusalem"]["junior_wing"] == [135, 581] and out["jerusalem"]["senior_wing"] == [18, 583]

# ---- flight: label of the near stratum in the stored source
xml = ROOT / "code/bsc_outbreaks/verify/src/khanh2020.xml"
if xml.exists():
    t = xml.read_text(encoding="utf-8", errors="replace")
    out["flight_label"] = dict(underlined_less_than_2_seats_away=t.count("<underline>&lt;</underline>2 seats away"),
                               greater_than_2_seats_away=t.count("&gt;2 seats away"),
                               reading="the source prints an underlined '<', i.e. '<= 2 seats away' against '> 2 seats away'")

# ---- Hunan coach: three denominators
t2 = rows["T2_hunan_coach"]
out["coach"] = dict(luo2020=int(t2["n_exposed"]), ou2022_dose_table="7/46" in t2["attack_rate_reported"] and 46,
                    ou2022_stratum_table=19 + 26, strata="3/19" in t2["spatial_pattern"] and "4/26" in t2["spatial_pattern"])

OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
print(json.dumps({k: v for k, v in out.items() if k != "mobility_gap"}, indent=1, ensure_ascii=False, default=str)[:2500])

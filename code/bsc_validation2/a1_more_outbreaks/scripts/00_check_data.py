"""Transcription checks against totals printed in the sources.  No model."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from a1 import events as E, digitized as dg
out = {}
for ev in E.EVENTS:
    m, d = E.build(ev, {"n_cfg": 2})
    out[ev] = dict(label=m["label"], sizes=m["sizes"], obs=m["obs"], K=m["K"], n=int(sum(m["sizes"])),
                   n_src=int(len(d[0]["src"])), near_bins=m["near_bins"])
    print(ev, out[ev])
chk = {}
m, _ = E.build("F1")
d = m["dist"]; y = m["y"]
# Olsen: 8/23 in the index row and the three rows in front, 10/88 elsewhere; Hertzberg: 9/29 within 2 rows, 9/75 beyond
lines = [dg.OLSEN[k] for k in "ABCDEF"]
front = [(i, r) for i in range(6) for r in range(10, 14) if lines[i][r] in "X/G"]
chk["F1_front3_plus_row"] = [sum(lines[i][r] == "G" for i, r in front), len(front)]
chk["F1_within2"] = [int(y[d <= 2].sum()), int((d <= 2).sum())]
chk["F1_beyond2"] = [int(y[d > 2].sum()), int((d > 2).sum())]
m, _ = E.build("F3"); chk["F3_within2"] = [int(m["y"][m["dist"] <= 2].sum()), int((m["dist"] <= 2).sum())]
chk["F3_beyond2"] = [int(m["y"][m["dist"] > 2].sum()), int((m["dist"] > 2).sum())]
m, _ = E.build("F4"); chk["F4_within2"] = [int(m["y"][m["dist"] <= 2].sum()), int((m["dist"] <= 2).sum())]
chk["F4_beyond2"] = [int(m["y"][m["dist"] > 2].sum()), int((m["dist"] > 2).sum())]
m, _ = E.build("F5"); chk["F5_within2"] = [int(m["y"][m["dist"] <= 2].sum()), int((m["dist"] <= 2).sum())]
chk["F5_beyond2"] = [int(m["y"][m["dist"] > 2].sum()), int((m["dist"] > 2).sum())]
m, _ = E.build("F6"); chk["F6_within2"] = [int(m["y"][m["dist"] <= 2].sum()), int((m["dist"] <= 2).sum())]
chk["F6_beyond2"] = [int(m["y"][m["dist"] > 2].sum()), int((m["dist"] > 2).sum())]
m, _ = E.build("F8"); chk["F8_within2"] = [int(m["y"][m["dist"] <= 2].sum()), int((m["dist"] <= 2).sum())]
chk["F8_beyond2"] = [int(m["y"][m["dist"] > 2].sum()), int((m["dist"] > 2).sum())]
for k, v in chk.items(): print(k, v)
json.dump(dict(events=out, checks=chk), open(os.path.join(os.path.dirname(__file__), "..", "out", "00_check_data.json"), "w"), indent=1)

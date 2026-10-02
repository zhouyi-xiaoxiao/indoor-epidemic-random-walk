#!/usr/bin/env python
"""s6_scenes_paste_tables.py -- paste the generated table bodies into sections/m6_simulation.tex and supp_sim.tex.

Reads data/s6_scenes_tables.tex (written by code/s6_scenes_tables.py) and replaces the text
between the marker lines "% BEGIN AUTO <label>" and "% END AUTO <label>" of sections/m6_simulation.tex (main
text) and sections/supp_sim.tex (Supplementary Material), so that no entry of these tables is typed by hand.

Run:  python code/s6_scenes_tables.py
      python code/s6_scenes_paste_tables.py
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "data" / "s6_scenes_tables.tex"
DSTS = [HERE.parent / "sections" / "m6_simulation.tex", HERE.parent / "sections" / "supp_sim.tex"]
LABELS = ["s6_scenes:tab:scenes", "s6_scenes:tab:pmajor", "s6_scenes:tab:hotspot:a", "s6_scenes:tab:hotspot:b"]

blocks, cur = [], None
for line in SRC.read_text().splitlines():
    if line.startswith("% ---- body of Table"):
        cur = []
        blocks.append(cur)
    elif cur is not None and line.strip():
        cur.append(line)
assert len(blocks) == len(LABELS), (len(blocks), len(LABELS))

texs = {d: d.read_text() for d in DSTS}
for label, body in zip(LABELS, blocks):
    if body and body[0].strip() == "\\addlinespace":      # no extra space directly below \midrule
        body = body[1:]
    pat = re.compile(r"(% BEGIN AUTO " + re.escape(label) + r"\n).*?(% END AUTO " + re.escape(label) + r"\n)",
                     re.S)
    hits = [d for d in DSTS if pat.search(texs[d])]
    assert len(hits) == 1, (label, hits)
    texs[hits[0]] = pat.sub(lambda m: m.group(1) + "\n".join(body) + "\n" + m.group(2), texs[hits[0]])
for d, tex in texs.items():
    d.write_text(tex)
print("pasted", len(LABELS), "table bodies into", [d.name for d in DSTS])

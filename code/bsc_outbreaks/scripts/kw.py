"""Keyword-in-context search over a saved full text. Usage: kw.py file kw1 kw2 ... (case-insensitive, prints +-width chars)"""
import sys, re, pathlib
p = pathlib.Path(__file__).resolve().parents[1] / "notes" / "fulltext" / sys.argv[1]
t = p.read_text()
W = 350
seen = []
for kw in sys.argv[2:]:
    for m in re.finditer(re.escape(kw), t, flags=re.I):
        a, b = max(0, m.start() - W), min(len(t), m.end() + W)
        if any(abs(a - s) < W for s in seen):
            continue
        seen.append(a)
        print(f"[{kw}] ..." + t[a:b].replace("\n", " ") + "...\n")

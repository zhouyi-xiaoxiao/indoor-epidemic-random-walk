"""Keyword-in-context helper: python kwic.py <file> <regex> [width] [max]"""
import re, sys
p = sys.argv[1]; pat = sys.argv[2]; w = int(sys.argv[3]) if len(sys.argv) > 3 else 180; mx = int(sys.argv[4]) if len(sys.argv) > 4 else 15
t = open(p, encoding="utf8", errors="replace").read()
t = re.sub(r"\s+", " ", t)
n = 0; last = -10**9
for m in re.finditer(pat, t, flags=re.I):
    if m.start() - last < w: continue
    last = m.start(); n += 1
    print(f"[{m.start()}] ...{t[max(0,m.start()-w):m.end()+w]}...")
    if n >= mx: break
print(f"-- {n} shown")

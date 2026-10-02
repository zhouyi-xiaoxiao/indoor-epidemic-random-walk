"""Predictive tables with the round-off-safe M2 kernel (POSTHOC_LOG P5).  M2 and M3 are recomputed;
M0, CRR and M1 are copied from the primary table.  Usage: 01b_tables_fix.py EVENT"""
import sys, os, time
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, engine as G, engine_fix as F
F.install()
ev = sys.argv[1]
fn = os.path.join(ROOT, "out", "tab_%s_nfix.npz" % ev)
if os.path.exists(fn):
    print("cached", fn); sys.exit()
t0 = time.time()
z = np.load(os.path.join(ROOT, "out", "tab_%s_primary.npz" % ev))
res = {k: z[k] for k in z.files}
meta, draws = E.build(ev)
assert (G.compositions(meta["sizes"], meta["K"]) == res["comps"]).all()
for m in ("M2", "M3"):
    res[m] = G.table(meta, draws, m, res["comps"])
    print(ev, "nfix", m, res[m].shape, "%.0fs" % (time.time() - t0), flush=True)
np.savez_compressed(fn, **res)

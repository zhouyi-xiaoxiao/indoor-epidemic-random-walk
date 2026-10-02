"""Predictive tables: for every event and model, the conditional pmf of the bin-count
vector at every grid value of the model parameter.  Reads the event total K only;
never reads the stratum counts.  Usage: 01_tables.py EVENT [variant-json] [tag]"""
import sys, os, json, time
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, engine as G
ev = sys.argv[1]
opt = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
tag = sys.argv[3] if len(sys.argv) > 3 else "primary"
models = opt.pop("models", ["M0", "CRR", "M2", "M1", "M3"])
fn = os.path.join(ROOT, "out", "tab_%s_%s.npz" % (ev, tag))
if os.path.exists(fn):
    print("cached", fn); sys.exit()
t0 = time.time()
meta, draws = E.build(ev, opt)
comps = G.compositions(meta["sizes"], meta["K"])
res = dict(comps=comps, sizes=np.array(meta["sizes"]), K=meta["K"])
for m in models:
    res[m] = G.table(meta, draws, m, comps)
    print(ev, tag, m, res[m].shape, "%.0fs" % (time.time() - t0), flush=True)
np.savez_compressed(fn, **res)

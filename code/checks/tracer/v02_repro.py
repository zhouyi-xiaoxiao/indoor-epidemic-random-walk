"""Reproducibility check: re-run s01..s05 with outputs redirected to verify/repro/, then compare every
number with results/ (the frozen files are not touched)."""
import os, sys, runpy, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "bsc_validation2", "a3_tracer_physics")
sys.path.insert(0, os.path.join(ROOT, "src"))
import a3lib as A
ORIG = A.RES
A.ROOT = os.path.join(HERE, "repro"); A.RES = os.path.join(A.ROOT, "results")
steps = sys.argv[1:] or ["s01_kernels", "s02_predict", "s03_power", "s04_evaluate", "s05_descriptive"]
for s in steps:
    t0 = time.time(); print("==", s, flush=True)
    sys.argv = [s]
    runpy.run_path(os.path.join(ROOT, "src", s + ".py"), run_name="__main__")
    print("== done", s, round(time.time() - t0, 1), "s", flush=True)

def walk(a, b, path, bad):
    if isinstance(a, dict):
        for k in a:
            if k == "stamp": continue
            if k not in b: bad.append((path + "/" + k, "missing")); continue
            walk(a[k], b[k], path + "/" + k, bad)
    elif isinstance(a, list):
        if len(a) != len(b): bad.append((path, "len")); return
        if a and all(isinstance(x, (int, float)) for x in a) and all(isinstance(x, (int, float)) for x in b):
            x, y = np.array(a, float), np.array(b, float)
            d = np.max(np.abs(x - y) / (1e-12 + np.maximum(np.abs(x), np.abs(y)))) if len(x) else 0
            if d > 1e-9: bad.append((path, float(d)))
        else:
            for i, (x, y) in enumerate(zip(a, b)): walk(x, y, path + "/%d" % i, bad)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if abs(a - b) > 1e-9 * max(1e-12, abs(a), abs(b)): bad.append((path, a, b))
    elif a != b: bad.append((path, a, b))
rep = {}
for f in ("01_kernels.json", "02_predictions.json", "03_power.json", "04_evaluation.json", "05_descriptive.json"):
    pa, pb = os.path.join(ORIG, f), os.path.join(A.RES, f)
    if not os.path.exists(pb): continue
    bad = []; walk(json.load(open(pa)), json.load(open(pb)), "", bad)
    rep[f] = dict(n_differences=len(bad), first=bad[:8]); print(f, len(bad), bad[:8])
json.dump(rep, open(os.path.join(A.ROOT, "repro_report.json"), "w"), indent=1, default=str)

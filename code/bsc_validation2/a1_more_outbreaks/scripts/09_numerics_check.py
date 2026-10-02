"""Diagnose the round-off floor of the frozen spectral kernel and validate the corrected kernel."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, engine as G, engine_fix as F
out = {}
for ev in E.EVENTS:
    meta, draws = E.build(ev, {"n_cfg": 7})
    rows = []
    for ld in G.LOGD[::4]:
        D = 10 ** ld; neg = 0; n = 0; relmax = 0.0; floor = 0; negM1 = 0
        for d in draws[:7:3]:
            Xs, xb = G.x_m2(d, D); Xi = F.x_images(d, D)[d["rec"]]
            X1, _ = G.x_m2(d, D, "M1")
            neg += int((Xs <= 0).sum()); negM1 += int((X1 <= 1e-12 * X1.max()).sum()); n += len(Xs)
            ok = Xs > 1e-6 * Xs.max()
            relmax = max(relmax, float(np.abs(Xi[ok] / Xs[ok] - 1).max()))
            floor += int((Xs < F.THRESH * Xs.max()).sum())
        rows.append(dict(D=float(D), frac_nonpos_spectral=neg / n, frac_replaced=floor / n, max_rel_diff_where_resolved=relmax, frac_unresolved_M1=negM1 / n))
    out[ev] = rows
    bad = [r["D"] for r in rows if r["frac_replaced"] > 0]
    print(ev, "max rel diff image vs spectral (resolved sites) %.2e" % max(r["max_rel_diff_where_resolved"] for r in rows),
          "| receivers replaced for D up to %s (max fraction %.2f)" % ("%.3g" % max(bad) if bad else "none", max(r["frac_replaced"] for r in rows)),
          "| M1 unresolved up to D=%.3g" % (max([r["D"] for r in rows if r["frac_unresolved_M1"] > 0] or [0])), flush=True)
json.dump(out, open(os.path.join(ROOT, "out", "09_numerics_check.json"), "w"), indent=1)

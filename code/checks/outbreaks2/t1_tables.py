import numpy as np, time, sys, json, os
import vlib as V
os.makedirs('cache',exist_ok=True)
evs = sys.argv[1:] or ["F1","F2","F3","F4","F5","F6","F7","F8","W1","W2","P1","W1x"]
for ev in evs:
    t=time.time()
    g=V.build(ev)
    pmf,comps=V.m2_table(g,cache='cache/m2_%s.npz'%ev)
    print(ev,g['sizes'],g['obs'],g['K'],pmf.shape,'%.0fs'%(time.time()-t),flush=True)

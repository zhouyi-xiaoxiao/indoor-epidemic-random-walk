import numpy as np, time, sys, json
import vlib as V
sys.path.insert(0,'../../bsc_validation2/a1_more_outbreaks/src')
from a1 import lattice as L, engine as G, engine_fix as F, events as E
out={}
# W2 geometry, compare three kernels at several D
g=V.build("W2")
meta,draws=E.build("W2")
d=draws[3]
for D in [0.02,0.5,2.5,16,100,794]:
    t=time.time()
    Xu=V.kern_table(g['nx'],g['ny'],1,1,g['srcsets'],100.0,[d['kappa']],[D])[0,0,0]
    Xe=V.kern_eigh(g['nx'],g['ny'],1,1,g['srcsets'],100.0,[d['kappa']],D)[0,0]
    Xs,_=G.x_m2(d,D)          # first-analysis spectral at receivers
    lat=d['lat']; 
    Xi=F.x_images(d,D)[d['rec']]
    r=d['rec']
    rel=lambda a,b: float(np.max(np.abs(a-b)/np.abs(b)))
    out[D]=dict(rel_img_vs_uni=rel(Xi,Xu[r]), rel_spec_vs_uni=rel(Xs,Xu[r]), rel_eigh_vs_uni=rel(Xe[r],Xu[r]),
                min_over_max=float(Xu[r].min()/Xu[r].max()), n_spec_nonpos=int((Xs<=0).sum()), t=time.time()-t)
    print(D,out[D],flush=True)
# F5 cabin
g=V.build("F5"); meta,draws=E.build("F5"); d=draws[3]
for D in [0.5,5,20,126,1000]:
    Xu=V.kern_table(g['nx'],g['ny'],.5,.8,g['srcsets'],2.0,[d['kappa']],[D])[0,0,0]
    Xs,_=G.x_m2(d,D); Xi=F.x_images(d,D)[d['rec']]; r=d['rec']
    same=np.array_equal(r,g['cfgs'][0][2])
    rel=lambda a,b: float(np.max(np.abs(a-b)/np.abs(b)))
    out['F5_%g'%D]=dict(rec_same=bool(same),rel_img_vs_uni=rel(Xi,Xu[r]), rel_spec_vs_uni=rel(Xs,Xu[r]),min_over_max=float(Xu[r].min()/Xu[r].max()))
    print('F5',D,out['F5_%g'%D],flush=True)
json.dump(out,open('t0_kernel.json','w'),indent=1)

import json, numpy as np
src=open('t2_analysis.py').read().split("# ------------------------------------------------------------------ 1. primary")[0]
exec(src)
res={}
for nm,edges in (("2 bins {0-2},{3+}",(2,)),("5 bins {0},{1},{2},{3-5},{6+}",(0,1,2,5)),("4 bins {0},{1-2},{3-5},{6+}",(0,2,5)),("3 bins {0-2},{3-5},{6+}",(2,5))):
    T2=dict(T); G2=dict(G)
    for e in ("F1","F3","F4","F5","F6","F8"):
        g=V.build(e,edges=edges)
        tag=''.join(map(str,edges))
        pmf,comps=V.m2_table(g,cache='cache/m2_%s_e%s.npz'%(e,tag))
        T2[e]=dict(comps=comps,M2=pmf,M0=V.m0_pmf(g,comps)[None,:],CRR=V.crr_table(g,comps),i=V.idx_of(comps,g['obs'])); G[e]=g
    r=run(CAL,HOLD,tabs=T2,label=nm)
    res[nm]=strip(r)
    print('    F5',G['F5']['sizes'],G['F5']['obs'],{m:(round(v['logp'],2),round(v['p_adeq'],4)) for m,v in r['events']['F5'].items() if m!='obs'})
json.dump(res,open('t5_bins.json','w'),indent=1,default=float)

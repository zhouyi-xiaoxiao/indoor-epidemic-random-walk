"""Does the frozen F3-O fitting procedure pass its own test when RWstickO IS the truth? (false-rejection rate of stage 3)"""
import sys,json,numpy as np
sys.path.insert(0,'../../bsc_validation2/a2_contact_mobility/src'); import a2lib as A, pipeline as PL, fix3
def synth_presence(rng,N,nd):
    pres=[]
    for d in range(nd):
        ids=np.flatnonzero(rng.random(N)<0.9); a=rng.integers(1440,1620,len(ids)); b=rng.integers(2880,3060,len(ids)); pres.append(dict(ids=ids,a=a,b=b))
    return pres
out=[]
for r in range(int(sys.argv[1])):
    rng=np.random.default_rng([r,555]); N=120; nd=4; groups=np.arange(N)%4
    pres=synth_presence(rng,N,nd)
    true=dict(M=1600.,p_fast=0.4,bias=0.2,sigma=2.0,rho=0.1)
    g=A.generate('RWstickO',rng,pres,N,groups,{}, list(range(nd)),extra=true)
    ds=dict(name=f'SYNO{r}',N=N,groups=groups,days=[g[d] for d in range(nd)])
    calfun=fix3.make_calfun(['RWstickO'],None,refine=True)
    res=PL.run_dataset(ds,['RWstickO'],G=5,G_epi=10,nrep_real=4000,nrep_model=400,seed=300+r,calfun=calfun,verbose=False)
    sc=res['models']['RWstickO']['score']; fit={k:v for k,v in res['extra']['RWstickO'].items() if k!='trace'}
    row=dict(rep=r,S=sum(bool(sc[f]) for f in A.FAMILIES),failed=[f for f in A.FAMILIES if sc[f] is False],E=res['models']['RWstickO']['epi_pass_count'],fit=fit,
             totratio=res['models']['RWstickO']['epi_total_ratio'],obs_c=res['cal']['c'],
             ratios={k:round(v['RWstickO']['R_index']/v['REAL']['R_index'],3) for k,v in res['epi'].items()},
             dA={k:round(v['RWstickO']['attack']-v['REAL']['attack'],3) for k,v in res['epi'].items()})
    print(json.dumps(row,default=float),flush=True); out.append(row)
    json.dump(out,open('v5_selfpower.json','w'),indent=1,default=float)

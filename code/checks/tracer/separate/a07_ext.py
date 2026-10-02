import sys,os,numpy as np,itertools
sys.path.insert(0,'../../../bsc_validation2/a1_more_outbreaks/src')
from a1 import events as E1, engine as G
from advlib import Lat
lg=lambda a:np.log(np.clip(a,1e-300,None))
kap=30+6*(np.arange(7)+.5)/7+E1.DEP; print('DEP',E1.DEP)
def run(Ds,opts={},label=''):
    tot={'M0':0,'CRR2':0,'CRR3':0}; s=[]
    for ev in ['F5','F6','F7','F8']:
        meta,draws=E1.build(ev,dict(opts.get(ev,{})))
        sizes,K=meta['sizes'],meta['K']; comps=G.compositions(sizes,K); b0=draws[0]['bins']; near=np.isin(b0,meta['near_bins'])
        cf=[draws[0]]*7 if len(draws)<=7 else draws
        acc=np.zeros(len(comps)); n=0
        for i,d in enumerate(cf):
            dd=dict(d); dd['kappa']=float(kap[i%7])
            for D in Ds:
                X,_=G.x_m2(dd,float(D),'M2'); acc+=G.cond_pmf(X,dd['bins'],sizes,K,comps); n+=1
        p=acc/n; i=int(np.flatnonzero((comps==np.array(meta['obs'])).all(1))[0])
        pv=p[p<=p[i]*(1+1e-12)].sum()
        cm={'M0':G.cond_pmf(np.ones(sum(sizes)),b0,sizes,K,comps),'CRR2':G.cond_pmf(np.where(near,2.,1.),b0,sizes,K,comps),'CRR3':G.cond_pmf(np.where(near,3.,1.),b0,sizes,K,comps)}
        for c in tot: tot[c]+=lg(p[i])-lg(cm[c][i])
        s.append('%s obs=%s p=%.3f'%(ev,meta['obs'],pv))
    print(label,'|',' ; '.join(s),'| delta',{c:round(v,2) for c,v in tot.items()})
run([707.9,112.2],label='P registered (708,112)')
run([281.8],label='single geomean D=282')
run([707.9,112.2,125.9],label='all 3 FWD releases (5L,5A,5G)')
run([125.9],label='5G only')
run([1259.],label='coach D 1259')
run([25.1],label='round-1 D 25')
run([707.9,112.2],{'F5':{'hh':True}},'P, Naha households removed (a1 option hh)')
run([707.9,112.2],{'F6':{'outcome':'wide'}},'P, Speake outcome=wide option')
run([707.9,112.2],{'F6':{'src':'alt1'}},'P, Speake src alt1')
run([707.9,112.2],{'F6':{'src':'alt2'}},'P, Speake src alt2')
# Naha check with the solver of the re-check
meta,draws=E1.build('F5'); d=draws[0]; lat=d['lat']; print('lat',lat.nx,lat.ny,lat.ax,lat.ay,'segs',d['segs'],'nsrc',len(d['src']),'K',meta['K'],'sizes',meta['sizes'])
my=Lat(lat.nx,lat.ny,lat.ax,lat.ay); X=sum(my.expo(112.2,33.93,d['segs'],int(s),d['rec']) for s in d['src'])
dd=dict(d); dd['kappa']=33.93; X2,_=G.x_m2(dd,112.2,'M2'); print('max rel diff exposure own vs a1:',np.abs(X/X2-1).max())

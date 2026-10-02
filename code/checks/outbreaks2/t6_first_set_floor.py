import sys, numpy as np, json
sys.path.insert(0,'../../bsc_validation/src')
import vlib as V
from valmod import events as E1, lattice as L1
rng=np.random.default_rng(1)
c=json.load(open('../../bsc_validation/data/01_cal_R1_primary.json'))
eps=np.array(c['M2']['eps']); logD=np.array(c['logD'])
res=[]
for trial in range(5):
    ev=E1.build_R1(rng); lat=ev['lat']
    src=(int(ev['src']%lat.nx), int(ev['src']//lat.nx))
    for D in [0.01,0.03,0.1,0.5,2,10,39]:
        e=eps[np.argmin(np.abs(logD-np.log10(D)))]
        worst=0; negs=0; dll=0
        for T in np.unique(ev['T']):
            sel=ev['T']==T
            Xs=L1.exposure(lat,'M2',D,ev['kappa'],[T],ev['src'],ev['sites'][sel])
            Xu=V.kern_table(lat.nx,lat.ny,1,1,[[src]],float(T),[ev['kappa']],[D])[0,0,0][ev['sites'][sel]]
            negs+=int((Xs<=0).sum())
            # impact on log-likelihood of non-cases at these tables: n * eps * |dX|
            dll+=float((ev['n'][sel]*e*np.abs(Xs-Xu)).sum())
            worst=max(worst,float(np.max(np.abs(Xs-Xu))/Xu.max()))
        res.append((trial,D,negs,worst,dll))
        print('layout',trial,'D',D,'n_nonpos',negs,'max abs err / max X %.1e'%worst,'bound on |dlogL| %.1e'%dll,'eps_grid %.3g'%e)
json.dump(res,open('t6_first_set_floor.json','w'))

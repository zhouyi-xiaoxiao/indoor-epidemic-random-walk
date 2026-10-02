import json,numpy as np,itertools
from advlib import p2
pred=json.load(open('../../../bsc_validation2/a3_tracer_physics/results/02_predictions.json'))['events']; EV=['T2','T1','T5','C1','R1']; TR=['P','M0','CRR2','CRR3']
pm={e:{m:np.array(pred[e]['pmf'][m]) for m in TR} for e in EV}
adv=json.load(open('pmfs.json'))
print('max |pmf diff| stored vs separate:',{e:{m:float(np.abs(pm[e][m]-np.array(adv[e][m])).max()) for m in TR} for e in EV})
tsp={e:np.array([p2(pm[e]['P'],k) for k in range(len(pm[e]['P']))]) for e in EV}
combos=list(itertools.product(*[range(len(pm[e]['P'])) for e in EV]))
W={t:np.array([np.prod([pm[e][t][k] for e,k in zip(EV,c)]) for c in combos]) for t in TR}
LG={t:np.array([sum(np.log(max(pm[e][t][k],1e-300)) for e,k in zip(EV,c)) for c in combos]) for t in TR}
nf=np.array([sum(tsp[e][k]<0.05 for e,k in zip(EV,c)) for c in combos])
def pval(j):  # exact one-sided p of delta under j for each combo
    d=LG['P']-LG[j]; o=np.argsort(-d); cw=np.cumsum(W[j][o]); p=np.empty(len(d)); 
    # ties: p = sum of w with d >= d_i
    ds=d[o]; 
    idx=np.searchsorted(-ds,-d-1e-12,side='right')-1; p=cw[idx]; return d,p
d0,p0=pval('M0'); d2,p2_=pval('CRR2'); d3,p3=pval('CRR3')
D2=(d0>0)&(p0<.05); D3=(d2>0)&(p2_<.05)&(d3>0)&(p3<.05)
# T4 simulation
rng=np.random.default_rng(11); t4=pred['T4']; piP=np.array(t4['P']['pi_draws']); pib=piP.mean(0)
def dev(k,p): e=k.sum(-1,keepdims=True)*p; return 2*np.where(k>0,k*np.log(np.where(k>0,k,1)/e),0).sum(-1)
null=dev(np.array([rng.multinomial(138,piP[i]) for i in rng.integers(len(piP),size=4000)]),pib); c95=np.quantile(null,.95)
pt={'P':piP,'M0':np.array(t4['M0']['pi'])[None],'CRR2':np.array(t4['CRR2']['pi'])[None],'CRR3':np.array(t4['CRR3']['pi'])[None]}
o1={'P':.95,'M0':.919,'CRR2':.919,'CRR3':.919}
for t in TR:
    k=np.array([rng.multinomial(138,pt[t][i]) for i in rng.integers(len(pt[t]),size=3000)])
    ok=dev(k,pib)<=c95; a=k @ np.log(piP).T; ll=np.log(np.exp(a-a.max(1,keepdims=True)).mean(1))+a.max(1); pos=(ll-k @ np.log(pt['M0'][0]))>0
    out=dict(S1=0,S2=0,S3=0,S4=0,S5=0)
    for a4 in(1,0):
      for ps in(1,0):
        p4=((ok==a4)&(pos==ps)).mean()
        for aO in(1,0):
            pO=o1[t] if aO else 1-o1[t]; ex=(1-a4)+(1-aO)
            for f in range(6):
                sel=nf==f; tot=f+ex; w=W[t][sel]
                if tot>=2: out['S5']+=p4*pO*w.sum()
                elif tot==1: out['S4']+=p4*pO*w.sum()
                else:
                    dd=D2[sel]&bool(ps); out['S1']+=p4*pO*w[dd&D3[sel]].sum(); out['S2']+=p4*pO*w[dd&~D3[sel]].sum(); out['S3']+=p4*pO*w[~dd].sum()
    print(t,{a:round(b,3) for a,b in out.items()},'T4 adequate',round(ok.mean(),3),'P(nofail5)',round(W[t][nf==0].sum(),3),'D2',round(W[t][D2].sum(),3),'D3',round(W[t][D3].sum(),3))
# Probability the actually observed qualitative pattern arises: "beats M0 (D2 stratum) but not both CRR"
for t in TR: print(t,'P(D2 & not D3)',round(W[t][D2&~D3].sum(),3),' P(D2&D3)',round(W[t][D2&D3].sum(),3))

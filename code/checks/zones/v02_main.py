import numpy as np, pandas as pd, json, time, vcore as v
t0=time.time(); rng=np.random.default_rng(777)
LY=np.array([0.7619836508613838,0.1044788607478434,0.13353748839077287])
cl,Y=v.load(); D=v.Data(cl,Y); S=D.S[:,1:]; Yd=Y[:,1:]
odd=cl.school.values%2==1; schools=np.unique(cl.school.values)
out={}
naive={'N_2/3_1/6_1/6':[2/3,1/6,1/6],'N_.8_.1_.1':[.8,.1,.1],'N_.7_.1_.2':[.7,.1,.2],'N_.5_.25_.25':[.5,.25,.25],'N_1/3each':[1/3]*3,'N_.9_.05_.05':[.9,.05,.05],'N_.6_.2_.2':[.6,.2,.2],
       'Lyon_distinctpairs':[.428,.1388,.4332],'Lyon_pairs>=60s':[.6777,.1362,.1861]}
specs={'Z':('Z',LY),'H':('H',None),'X':('X',None),'F':('F',None)}
for k,w in naive.items(): specs[k]=('Z',np.array(w))
ho={k:np.zeros(len(cl)) for k in list(specs)+['C']}; cv={}
for cal,name in ((odd,'cal_odd'),(~odd,'cal_even')):
    tst=~cal
    for k,(m,w) in specs.items():
        X=D.feats(m,w=w); th,l=v.fit(X[cal],S[cal],Yd[cal]); ho[k][tst]=v.ll_cls(th,X[tst],S[tst],Yd[tst]); cv[k+'_'+name]=th.tolist()
    th,l,rho=v.fitC(D,S,Yd,cal); ho['C'][tst]=v.ll_cls(th,D.feats('C',rho=rho)[tst],S[tst],Yd[tst]); cv['C_'+name]=th.tolist()+[rho]
ps={k:np.array([ho[k][cl.school.values==s].sum() for s in schools]) for k in ho}
out['heldout']={k:float(a.sum()) for k,a in ps.items()}; out['cv']=cv
idx=rng.integers(0,len(schools),(10000,len(schools)))
for k in ps:
    if k=='Z': continue
    d=ps['Z']-ps[k]; bs=d[idx].sum(1); out['D_Z_minus_'+k]=[float(d.sum()),float(np.quantile(bs,.025)),float(np.quantile(bs,.975)),int((d>0).sum())]
# stratified-by-fold bootstrap for robustness
full={}
for k,(m,w) in specs.items():
    th,l=v.fit(D.feats(m,w=w),S,Yd); full[k]=dict(theta=th.tolist(),ll=float(l))
th,l,rho=v.fitC(D,S,Yd,slice(None)); full['C']=dict(theta=th.tolist(),ll=float(l),rho=rho)
out['full']=full
b=np.array(full['F']['theta'][:3]); w=b/b.sum(); out['w_hat']=w.tolist(); out['TV']=float(0.5*np.abs(w-LY).sum()); out['LR']=2*(full['F']['ll']-full['Z']['ll'])
out['TV_naive']={k:float(0.5*np.abs(w-np.array(x)).sum()) for k,x in naive.items()}
out['LR_naive']={k:2*(full['F']['ll']-full[k]['ll']) for k in naive}
out['disp_obs']=float(v.disp(cl,Y)); out['AR']=float(Y.sum()/cl.n.sum())
json.dump(out,open('v_main.json','w'),indent=1); print('fits',time.time()-t0)
# bootstrap TV
rows={s:np.nonzero(cl.school.values==s)[0] for s in schools}; XF=D.feats('F'); W=[]
for i in range(300):
    ix=np.concatenate([rows[s] for s in rng.choice(schools,len(schools))]); th,_=v.fit(XF[ix],S[ix],Yd[ix]); W.append(th[:3]/th[:3].sum())
W=np.array(W); tvs=0.5*np.abs(W-LY).sum(1)
out['w_boot']=[[float(np.quantile(W[:,k],q)) for q in (.025,.975)] for k in range(3)]; out['TV_boot']=[float(np.quantile(tvs,q)) for q in (.025,.5,.975)]; out['P_TV_gt_0.15']=float((tvs>0.15).mean())
json.dump(out,open('v_main.json','w'),indent=1); print('boot',time.time()-t0)
# forward sims
for k in ('Z','H','X'):
    th=np.array(full[k]['theta']); ds=[];ars=[]
    for i in range(300):
        Ys=v.simulate(cl,k,th,rng,w=LY); ds.append(v.disp(cl,Ys)); ars.append(Ys.sum()/cl.n.sum())
    out['disp_'+k]=[float(np.quantile(ds,q)) for q in (.025,.5,.975)]; out['AR_'+k]=[float(np.quantile(ars,q)) for q in (.025,.5,.975)]
    json.dump(out,open('v_main.json','w'),indent=1); print(k,time.time()-t0)
th=np.array(full['C']['theta']); ds=[]
for i in range(300): ds.append(v.disp(cl,v.simulate(cl,'C',th,rng,rho=full['C']['rho'])))
out['disp_C']=[float(np.quantile(ds,q)) for q in (.025,.5,.975)]
json.dump(out,open('v_main.json','w'),indent=1); print('done',time.time()-t0)

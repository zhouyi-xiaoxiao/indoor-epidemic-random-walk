import json, glob, numpy as np, pandas as pd
rows=[json.loads(l) for f in sorted(glob.glob('results/power_*.jsonl')) for l in open(f)]
d=pd.DataFrame(rows)
g=d.groupby(['scenario','truth'])[['tier1','tier2a','tier2b','tier3','overall']].mean().round(3); g['n']=d.groupby(['scenario','truth']).size()
print(g.to_string())
out={'table':g.reset_index().to_dict('records')}
# false-pass detail for random generic truths
for t in ('C','F'):
    s=d[d.truth==t]
    if t=='F':
        near=s.TV_true<=0.15
        print('F truths: frac with true TV<=0.15:',near.mean().round(3),' overall pass | far:',s[~near].overall.mean().round(3),' tier2b pass | far:',s[~near].tier2b.mean().round(3),'| near: overall',s[near].overall.mean().round(3))
        out['F_far_overall']=float(s[~near].overall.mean()); out['F_far_n']=int((~near).sum())
    else:
        for lo,hi in ((1,3),(3,10),(10,30),(30,100),(100,1000)):
            q=s[(s.rho>=lo)&(s.rho<hi)]
            print(f'C truth rho in [{lo},{hi}): n={len(q)} tier1={q.tier1.mean():.2f} tier2b={q.tier2b.mean():.2f} tier3={q.tier3.mean():.2f} overall={q.overall.mean():.2f}')
zz=d[d.truth=='Z']; print('Z truth: LR quantiles',zz.LR_F_vs_Z.quantile([.5,.9,.95]).round(2).tolist(),' TV q95',zz.TV.quantile(.95).round(3), 'w_hat mean',np.mean(zz.w_hat.tolist(),0).round(3))
# analytic: prior volume of the simplex within TV<=0.15 of Lyon
rng=np.random.default_rng(0); w=rng.dirichlet([1,1,1],2_000_000); L=np.array([0.7619836508613838,0.1044788607478434,0.13353748839077287])
v=float((0.5*np.abs(w-L).sum(1)<=0.15).mean()); print('simplex volume fraction within TV<=0.15:',round(v,4)); out['simplex_fraction']=v
json.dump(out,open('results/power_summary.json','w'),indent=1)

import json, numpy as np
exec(open('t2_analysis.py').read().split("# ------------------------------------------------------------------ 1. primary")[0])
t4=json.load(open('../../bsc_validation/data/01_cal_T4_primary.json'))
ll=np.array(t4['M2']['profile_logL']); w=np.exp(ll-ll.max()); w=w.sum(axis=1) if w.ndim==2 else w; w/=w.sum()
print('T4 posterior', summ(w,V.LOGD))
air=[e for e in EV if e[0]=='F']
r=run({'air':CAL['air']},{'air':air},weights={('air','M2'):w},label='S-T4 transfer, 8 flights (CRR calibrated F1-F4)')
for e in r['events']: print('  ',e,{m:(round(v['logp'],2),round(v['p_adeq'],4)) for m,v in r['events'][e].items() if m!='obs'})
r=run({'air':CAL['air']},{'air':HOLD['air']},weights={('air','M2'):w},label='S-T4 transfer, covid flights')
json.dump({e:{m:v for m,v in r['events'][e].items()} for e in r['events']},open('t4_transfer.json','w'),default=float)

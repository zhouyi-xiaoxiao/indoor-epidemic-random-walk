import json, sys
import numpy as np
K=['S0_total','S1_mean_dur','S1_timefrac_ge15','S2_burst','S2_recurrence','S2_median_gap','S3_deg_mean','S3_deg_cv','S4_strength_cv','S5_persistence','S6_within_frac','S7_alpha','S7_n_range']
for p in sys.argv[1:]:
    r=json.load(open(p)); print('=====',p,'N',r['N'],'days',r['n_days'],'sec',round(r['seconds']))
    c=r['cal']; print('cal c=%.3g M=%.0f mean_dur=%.2f p=%.3f het=%s c_in=%.3g c_out=%.3g'%(c['c'],c['M'],c['mean_dur'],c['p'],{k:round(v,3) for k,v in c['het'].items()},c['c_in'],c['c_out']))
    if r.get('extra'): print('extra',r['extra'])
    print('%-18s'%'stat','%9s'%'OBS',*['%9s'%m[:9] for m in r['models']])
    for k in K:
        print('%-18s'%k,'%9.3g'%r['obs'][k],*['%9.3g'%r['models'][m]['stats'][k] for m in r['models']])
    for f in ['S0_level','S1_duration','S2_intercontact','S3_degree','S4_heterogeneity','S5_persistence','S6_groups','E_epidemic']:
        print('%-18s'%f,'%9s'%'',*['%9s'%{True:'pass',False:'FAIL',None:'n/a'}[r['models'][m]['score'].get(f)] for m in r['models']])
    for k,row in r.get('epi',{}).items():
        print(k,'beta=%.3g'%row['beta'])
        for m,e in row.items():
            if m=='beta': continue
            print('   %-8s R=%.2f Pmaj=%.3f AR=%.3f ARmaj=%.3f %s'%(m,e['R_index'],e['P_major'],e['attack'],e['attack_major'],e.get('pass','')))

import json, numpy as np, sys
d=json.load(open('../../bsc_validation2/a1_more_outbreaks/out/06_variants.json'))
for k,v in d.items():
    try:
        a=v['pooled']['M2']['air']; print('%-40s air dC %+.2f pB2 %.4f mirror %.4f'%(k,a['deltaC'],a['p_B2'],a['p_mirror']))
    except Exception as ex: print(k, 'ERR', list(v.keys())[:5] if isinstance(v,dict) else type(v))
s=json.load(open('../../bsc_validation2/a1_more_outbreaks/out/05_secondary.json'))
print(list(s.keys()))
a=s['CRR_2.4_primary_split']['pooled']['M2']; print('CRR2.4', {k:(round(v['deltaC'],2),round(v['p_B2'],4),round(v['p_mirror'],4)) for k,v in a.items() if k in('all','air','room')})

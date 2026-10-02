"""Descriptives needed to define the design (no model, no contact-structure statistic): N, days, counts."""
import a2lib as A, json
out={}
for n in A.SPEC:
    ds=A.load(n); tr,te=A.split(len(ds['days']))
    out[n]=dict(N=ds['N'],day_index=ds['day_index'],contacts_per_day=[int(len(d['s'])) for d in ds['days']],
                present_per_day=[int(len(A.presence(d)['ids'])) for d in ds['days']],
                first_last_interval=[[int(d['s'].min()),int(d['s'].max())] for d in ds['days']],
                n_groups=(len(ds['group_names']) if ds['group_names'] else 0),group_names=ds['group_names'],train=tr,test=te)
    print(n,out[n])
json.dump(out,open(A.ROOT+'/out/00_describe.json','w'),indent=1)

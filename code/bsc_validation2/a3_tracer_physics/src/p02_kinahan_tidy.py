"""P2: tidy table of the Kinahan et al. in-flight releases (integrated tracer counts per sensor seat).
Output data/derived/kinahan_tidy.csv: airframe, section, release, test, kind(B/C), mask, sensor, row, col, count
(blank cell -> NaN, zero kept as 0)."""
import csv, glob, os, re
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = []
for f in sorted(glob.glob(os.path.join(ROOT, 'data', 'derived', 'kinahan_7*_*.csv'))):
    b = os.path.basename(f)
    if 'tidy' in b: continue
    af = b.split('_')[1]
    sec = re.sub(r'^kinahan_7\d7_(767_)?', '', b).replace('.csv', '').replace('_-_Inflight', '').strip('_')
    rows = list(csv.reader(open(f)))
    hidx = [i for i, r in enumerate(rows) if r and r[0].startswith('Release')]
    h = rows[hidx[0]]
    n = next(k for k in range(7, len(h)) if h[k] == '')
    sensors = [re.sub(r'^S\d+\s+', '', c).strip() for c in h[7:n]]
    for r in rows[hidx[0] + 1: hidx[1]]:
        if len(r) < 8 or r[2].strip().upper() not in ('B', 'C'): continue
        for s, v in zip(sensors, r[7:n]):
            m = re.fullmatch(r'(\d+)([A-L])', s)
            out.append(dict(airframe=af, section=sec, release=r[0].strip(), test=int(r[1]), kind=r[2].strip().upper(),
                            mask=r[4].strip(), sensor=s, row=int(m.group(1)) if m else np.nan,
                            col=m.group(2) if m else '', count=float(v) if v != '' else np.nan))
df = pd.DataFrame(out)
df.to_csv(os.path.join(ROOT, 'data', 'derived', 'kinahan_tidy.csv'), index=False)
print(df.groupby(['airframe', 'section', 'release', 'kind', 'mask']).agg(tests=('test', 'nunique'), sensors=('sensor', 'nunique'),
      n_nan=('count', lambda x: int(x.isna().sum())), n_zero=('count', lambda x: int((x == 0).sum())),
      med=('count', 'median')).to_string())
g = df[(df.airframe == '777') & (df.section == 'FWD') & (df.kind == 'B') & (df['mask'] == 'No')]
for rel, d in g.groupby('release'):
    print('\n777 FWD release', rel, '(mean over replicates; rows = seat row, cols = A D G L)')
    print(d.groupby(['row', 'col'])['count'].mean().unstack().round(0).to_string())

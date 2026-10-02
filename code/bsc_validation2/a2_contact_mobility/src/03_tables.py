"""Markdown tables from out/*.json (no computation on data; pure formatting). usage: 03_tables.py > ../out/tables.md"""
import json, os, glob, sys
import numpy as np
import a2lib as A
OUT = A.ROOT + '/out'
FAM = ['S0_level', 'S1_duration', 'S2_intercontact', 'S3_degree', 'S4_heterogeneity', 'S5_persistence', 'S6_groups', 'E_epidemic']
KEYS = [('S0_total', 'total contact (pair-intervals)'), ('S1_mean_dur', 'mean event duration (x20 s)'), ('S1_timefrac_ge15', 'share of contact time in events >=5 min'),
        ('S2_burst', 'burstiness B of gaps'), ('S2_recurrence', 'recurrence 1-pairs/events'), ('S3_deg_mean', 'distinct contacts per person-day'),
        ('S3_deg_cv', 'CV of daily degree'), ('S4_strength_cv', 'CV of individual contact time'), ('S5_persistence', 'pairs repeated from previous day'),
        ('S6_within_frac', 'within-group share of contact time'), ('S7_alpha', 'density exponent alpha (descriptive)')]
SYM = {True: 'pass', False: 'FAIL', None: 'n/a'}


def merged(ds, tags):
    """Merge results of several tags for one dataset; later tags add models."""
    base = None
    for t in tags:
        p = f'{OUT}/{t}_{ds}.json'
        if not os.path.exists(p): continue
        r = json.load(open(p))
        if base is None: base = r; continue
        for m, v in r['models'].items():
            if m != 'REALtrain': base['models'][m] = v
        for k, row in r.get('epi', {}).items():
            for m, e in row.items():
                if m not in ('beta', 'REAL'): base['epi'][k][m] = e
            base['epi'][k]['REAL_' + t] = row['REAL']
        base.setdefault('extra', {}).update(r.get('extra', {}))
    return base


def g(x, f='%.3g'):
    return 'n/a' if x is None or (isinstance(x, float) and np.isnan(x)) else f % x


def tables(dss, tags, models, title):
    L = [f'### {title}', '']
    for ds in dss:
        r = merged(ds, tags)
        if r is None: continue
        ms = [m for m in models if m in r['models']]
        L += [f"**{ds}** — N={r['N']}, {r['n_days']} days (train {len(r['train'])}, test {len(r['test'])}), {r['n_groups']} groups; "
              f"calibrated c={r['cal']['c']:.3g} (M=1/c={r['cal']['M']:.0f} sites), mean event duration {r['cal']['mean_dur']:.2f} intervals, p={r['cal']['p']:.3f}", '']
        L += ['| statistic (test days) | observed | train days (ceiling) | ' + ' | '.join(ms) + ' |', '|---|---|---|' + '---|' * len(ms)]
        for k, lab in KEYS:
            L.append(f"| {lab} | {g(r['obs'].get(k))} | {g(r['obs_train'].get(k))} | " + ' | '.join(g(r['models'][m]['stats'].get(k)) for m in ms) + ' |')
        L += ['', '| family | train days (ceiling) | ' + ' | '.join(ms) + ' |', '|---|---|' + '---|' * len(ms)]
        for f in FAM:
            L.append(f"| {f} | {SYM[r['models']['REALtrain']['score'].get(f)]} | " + ' | '.join(SYM[r['models'][m]['score'].get(f)] for m in ms) + ' |')
        L.append('| S-families passed | ' + str(sum(bool(r['models']['REALtrain']['score'].get(f)) for f in FAM[:-1])) + ' | ' +
                 ' | '.join(str(sum(bool(r['models'][m]['score'].get(f)) for f in FAM[:-1])) for m in ms) + ' |')
        L += ['', '| SIR scenario | quantity | real network | ' + ' | '.join(ms) + ' |', '|---|---|---|' + '---|' * len(ms)]
        for k, row in r['epi'].items():
            for q, lab in [('R_index', 'R_index'), ('P_major', 'P(major)'), ('attack', 'mean attack')]:
                L.append(f"| {k} | {lab} | {row['REAL'][q]:.3f} | " + ' | '.join(f"{row[m][q]:.3f}" for m in ms) + ' |')
            L.append(f"| {k} | scenario passes | | " + ' | '.join(SYM[row[m]['pass']] for m in ms) + ' |')
        L.append('')
    return L


def summary(dss, tags, models, title):
    L = [f'### {title}', '', '| dataset | model | S-families passed / applicable | failed families | E scenarios passed (of 4) | mean R_index ratio model/real | mean dP(major) | mean d(attack) |', '|---|---|---|---|---|---|---|---|']
    for ds in dss:
        r = merged(ds, tags)
        if r is None: continue
        for m in models:
            if m not in r['models']: continue
            sc = r['models'][m]['score']
            app = [f for f in FAM[:-1] if sc.get(f) is not None]
            fail = [f.split('_')[0] for f in app if not sc[f]]
            rr = np.mean([row[m]['R_index'] / row['REAL']['R_index'] for row in r['epi'].values()])
            dp = np.mean([row[m]['P_major'] - row['REAL']['P_major'] for row in r['epi'].values()])
            da = np.mean([row[m]['attack'] - row['REAL']['attack'] for row in r['epi'].values()])
            L.append(f"| {ds} | {m} | {sum(sc[f] for f in app)}/{len(app)} | {', '.join(fail) or '-'} | {r['models'][m]['epi_pass_count']} | {rr:.2f} | {dp:+.3f} | {da:+.3f} |")
    L.append('')
    return L


if __name__ == '__main__':
    tags = ['s1', 's2a', 's2b', 's2d']
    M = ['WM', 'BLOCK', 'RW0', 'RWhet', 'RWhome', 'RWclus', 'RWstickZ', 'RWstickO']
    out = summary(A.DEV, tags, M, 'DEV summary') + tables(A.DEV, tags, M, 'DEV detail')
    ctags = ['s3a', 's3b', 's3c', 's3d']
    if glob.glob(f'{OUT}/s3a_*.json'):
        out += summary(A.CONF, ctags, M, 'CONF summary') + tables(A.CONF, ctags, M, 'CONF detail')
    print('\n'.join(out))

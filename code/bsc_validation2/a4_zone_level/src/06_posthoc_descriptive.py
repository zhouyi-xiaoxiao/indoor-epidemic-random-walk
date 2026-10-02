"""POST-HOC, descriptive (logged in POSTHOC_LOG.md, P1-P4). Not part of the registered decision rule."""
import sys, json, numpy as np
sys.path.insert(0, 'src'); import zonecore as z
R = json.load(open('results/matsumoto_primary.json')); P = R['primary']; out = {}
st, Y, info = z.load_outcomes()
# P1: recomputation, with separately written code, of headline quantities straight from the CSV (pandas, no zonecore features)
import pandas as pd
df = pd.read_csv('data/derived/matsumoto_students.csv')
out['check_counts'] = dict(n=int(len(df)), infected=int(df.isinfected.sum()), schools=int(df.schoolID.nunique()),
    classes=int(df.groupby(['schoolID', 'gradeID', 'classID']).ngroups))
g = df.groupby(['schoolID', 'gradeID', 'classID']).isinfected.agg(['sum', 'size']).reset_index()
num = 0; dof = 0
for s, gs in g.groupby('schoolID'):
    p = gs['sum'].sum()/gs['size'].sum()
    if len(gs) < 2 or p <= 0 or p >= 1: continue
    num += (((gs['sum']-gs['size']*p)**2)/(gs['size']*p*(1-p))).sum(); dof += len(gs)-1
out['check_dispersion'] = float(num/dof)
sch = df.groupby('schoolID').isinfected.agg(['sum', 'size']); ar = sch['sum']/sch['size']
out['school_attack_rate_range'] = [float(ar.min()), float(ar.median()), float(ar.max())]
out['school_size_range'] = [int(sch['size'].min()), float(sch['size'].median()), int(sch['size'].max())]
cl = g['sum']/g['size']; out['class_attack_rate_quantiles'] = [float(cl.quantile(q)) for q in (0, .25, .5, .75, 1)]
out['class_size_median'] = float(g['size'].median())
# P2: sign counts of the per-school held-out differences
for m in 'HXC':
    d = np.array(P[f'D_Z_minus_{m}_per_school']); out[f'schools_Z_ahead_of_{m}'] = [int((d > 0).sum()), len(d)]
    out[f'max_school_contrib_{m}'] = float(d.max()); out[f'sum_{m}'] = float(d.sum())
    # sign-flip-free check: leave-one-school-out range of the total
    out[f'leave_one_out_range_{m}'] = [float((d.sum()-d).min()), float((d.sum()-d).max())]
# P3: own-class / rest-of-school per-pair ratio implied by the Lyon shares in each Matsumoto school
#     (Z: per-pair weight to a classmate = w_c/(n_c-1); to a non-classmate of the school ~ (w_g+w_s)/(N_s-n_c))
w = z.LYON; n = st['n']; ns = st['n_school']
ok = (ns > n) & (n > 1); n, ns = n[ok], ns[ok]; rho_imp = (w[0]/(n-1))/((w[1]+w[2])/(ns-n))
out['rho_implied_by_Lyon_per_class_quantiles'] = [float(np.quantile(rho_imp, q)) for q in (.025, .25, .5, .75, .975)]
out['rho_implied_weighted_mean_log'] = float(np.exp(np.average(np.log(rho_imp), weights=n)))
out['rho_lyon_own_school'] = float(z.PERPAIR[0]/((1076.4655172413793+1375.8620689655172)/(232-23.2)))
out['rho_fitted_C'] = dict(full=P['full']['C']['rho'], cal_odd=P['cvfit']['C_cal_odd']['rho'], cal_even=P['cvfit']['C_cal_even']['rho'])
# P4: AIC/BIC on full data (k = number of fitted parameters)
k = dict(H=3, X=3, Z=3, C=4, F=5); nobs = info['n_inf']
out['AIC'] = {m: float(2*k[m]-2*P['full'][m]['ll']) for m in k}
# P5: Seoul 11th-floor wings, descriptive fitted coupling (NO test registered): north 79/137, south 4/62
LN, LS = -np.log(1-79/137), -np.log(1-4/62); zN, zS, zF = 79/137, 4/62, 83/199; r = LS/LN
d_first = r*zN/(zF+r*(zN-zF))                     # first-generation: L_S = c d zF ; L_N = c[(1-d) zN + d zF]
d_cum = (r*zN - zS)/((zF-zS) - r*(zF-zN))         # cumulative: L_S = c[(1-d) zS + d zF]
out['seoul_wings'] = dict(L_north=float(LN), L_south=float(LS), coupling_first_generation=float(d_first), coupling_cumulative=float(d_cum),
    note='cumulative form gives a negative coupling (boundary 0): not identifiable; first-generation form gives the share of time in common floor space that would reproduce the split. Fitted, not tested.')
json.dump(out, open('results/posthoc_descriptive.json', 'w'), indent=1); print(json.dumps(out, indent=1))

"""T2 (low-power consistency check): Seoul call-centre building, floor-level zone model with a shared lift/lobby zone.
All counts from Park et al. 2020 EID 26(8):1666, Table 1 (verified against Europe PMC XML, data/raw/park2020)."""
import numpy as np, json
from scipy import stats
N11, I11 = 216, 94
others = {'1st-6th': (84, 0), '7th': (182, 0), '8th': (207, 0), '9th': (206, 1), '10th': (27, 2), '13th-19th residents': (201, 0), 'visitors': (20, 0)}
No = sum(v[0] for v in others.values()); Io = sum(v[1] for v in others.values()); Ntot = N11+No
assert Ntot == 1143 and Io == 3 and No == 927
L11 = -np.log(1-I11/N11)
rng = np.random.default_rng(20261001); M = 200000
minutes = np.exp(rng.uniform(np.log(3), np.log(30), M)); eps = minutes/(9*60)
kappa = np.exp(rng.uniform(0, np.log(5), M))
def predict(eps, kappa):
    z11, ztot = I11/N11, I11/Ntot
    c = L11/((1-eps)*z11+kappa*eps*ztot)          # calibrated on the 11th floor only
    Lo = c*kappa*eps*ztot                         # first-generation spill-over hazard off the 11th floor
    return No*(1-np.exp(-Lo)), Lo
mu, Lo = predict(eps, kappa)
cnt = rng.poisson(mu)
out = dict(N_other=No, obs_other=Io, L11=float(L11), mu_quantiles={str(q): float(np.quantile(mu, q)) for q in (.025, .25, .5, .75, .975)},
           pred_count_interval95=[int(np.quantile(cnt, .025)), int(np.quantile(cnt, .975))], pred_count_median=float(np.median(cnt)),
           p_le_obs=float((cnt <= Io).mean()), p_ge_obs=float((cnt >= Io).mean()),
           central=dict(minutes=10, kappa=1, mu=float(predict(10/540, 1.0)[0])),
           well_mixed_expected=float(No*I11/N11*0+No*(I11+Io)/Ntot), independent_expected=0.0)
out['pass_obs3'] = bool(out['pred_count_interval95'][0] <= 3 <= out['pred_count_interval95'][1])
out['pass_obs2'] = bool(out['pred_count_interval95'][0] <= 2 <= out['pred_count_interval95'][1])
wm = No*(I11+Io)/Ntot; out['well_mixed_P_le_3'] = float(stats.binom.cdf(3, No, (I11+Io)/Ntot))
# per-floor expected counts at the prior median
for k, (n, i) in others.items():
    out.setdefault('per_floor', {})[k] = dict(n=n, obs=i, expected_median=float(np.median(n*(1-np.exp(-Lo)))), expected_95=[float(np.quantile(n*(1-np.exp(-Lo)), .025)), float(np.quantile(n*(1-np.exp(-Lo)), .975))])
# inverse (descriptive, NOT a test): shared-zone exposure kappa*eps that reproduces 3 (or 2) off-floor cases
from scipy.optimize import brentq
for k in (2, 3):
    f = lambda x: predict(x, 1.0)[0]-k
    out[f'implied_kappa_eps_for_{k}'] = float(brentq(f, 1e-6, 0.5))
json.dump(out, open('results/seoul_floors.json', 'w'), indent=1); print(json.dumps(out, indent=1))

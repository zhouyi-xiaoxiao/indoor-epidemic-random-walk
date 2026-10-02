"""Power analysis on SYNTHETIC epidemics using only the class structure of the Matsumoto respondents (no outcomes read).
usage: 02_power.py <scenario A|B> <nsim> ; checkpointed to results/power_<scenario>.jsonl"""
import sys, os, json, numpy as np
sys.path.insert(0, 'src'); import zonecore as z, analysis as an
from multiprocessing import Pool
SC = dict(A=dict(beta=0.8, a0=2e-5, a1=0.42), B=dict(beta=0.5, a0=2e-5, a1=0.70))
TRUTHS = ['Z', 'H', 'X', 'C', 'F']
def one(job):
    sc, truth, i = job; p = SC[sc]; st = z.load_structure()
    rng = np.random.default_rng([20261001, ord(sc), TRUTHS.index(truth), i])
    extra = {}
    if truth == 'C':
        rho = float(np.exp(rng.uniform(0, np.log(1000)))); extra['rho'] = rho
        Y = z.simulate(st, 'C', np.array([p['beta'], p['a0'], p['a1']]), rng, rho=rho)
    elif truth == 'F':
        w = rng.dirichlet([1, 1, 1]); extra['w'] = w.tolist(); extra['TV_true'] = float(z.tv(w, z.LYON))
        Y = z.simulate(st, 'F', np.r_[p['beta']*w, p['a0'], p['a1']], rng)
    else:
        Y = z.simulate(st, truth, np.array([p['beta'], p['a0'], p['a1']]), rng)
    r = an.run_tests(st, Y, rng, n_boot=1000, n_fwd=60)
    keep = {k: r[k] for k in ['tier1', 'tier2a', 'tier2b', 'tier3', 'overall', 'D_Z_minus_H', 'D_Z_minus_X', 'D_Z_minus_H_ci', 'D_Z_minus_X_ci', 'LR_F_vs_Z', 'TV', 'w_hat', 'disp_obs', 'disp_pred_Z', 'AR_obs']}
    return dict(scenario=sc, truth=truth, i=i, **extra, **keep)
if __name__ == '__main__':
    sc, nsim = sys.argv[1], int(sys.argv[2]); path = f'results/power_{sc}.jsonl'
    done = set()
    if os.path.exists(path):
        for l in open(path): d = json.loads(l); done.add((d['truth'], d['i']))
    jobs = [(sc, t, i) for t in TRUTHS for i in range(nsim) if (t, i) not in done]
    with Pool(4) as pool, open(path, 'a') as f:
        for r in pool.imap_unordered(one, jobs):
            f.write(json.dumps(r)+'\n'); f.flush()
    print('done', sc)

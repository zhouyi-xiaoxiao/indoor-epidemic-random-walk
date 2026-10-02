"""Shared analysis pipeline: identical code is applied to synthetic epidemics (power analysis) and to the real data."""
import numpy as np, zonecore as z

def school_folds(st):
    odd = (st['school'] % 2 == 1); return odd, ~odd

def run_tests(st, Y, rng, n_boot=2000, n_fwd=100, with_C=False, verbose=False):
    F = z.features(st, Y); odd, even = school_folds(st)
    models = ['H', 'X', 'Z'] + (['C'] if with_C else [])
    ho = {m: np.zeros(len(st['n'])) for m in models}; cvfit = {}
    for cal, tst, name in ((odd, even, 'cal_odd'), (even, odd, 'cal_even')):
        for m in models:
            r = z.fit_model(F, st, m, mask=cal); ho[m][tst] = z.heldout_ll(F, st, m, r, tst)
            cvfit[f'{m}_{name}'] = dict(theta=r['theta'].tolist(), rho=r['rho'], ll_cal=r['ll'])
    schools = np.unique(st['school'])
    per_school = {m: np.array([ho[m][st['school'] == s].sum() for s in schools]) for m in models}
    out = dict(heldout_ll={m: float(per_school[m].sum()) for m in models}, cvfit=cvfit)
    idx = rng.integers(0, len(schools), (n_boot, len(schools)))
    for m in models:
        if m == 'Z': continue
        d = per_school['Z'] - per_school[m]; bs = d[idx].sum(1)
        out[f'D_Z_minus_{m}'] = float(d.sum()); out[f'D_Z_minus_{m}_ci'] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
        out[f'D_Z_minus_{m}_per_school'] = d.tolist()
        # fold-specific
        o = (schools % 2 == 1); out[f'D_Z_minus_{m}_test_even'] = float(d[~o].sum()); out[f'D_Z_minus_{m}_test_odd'] = float(d[o].sum())
    out['tier1'] = bool(out['D_Z_minus_H_ci'][0] > 0 and out['D_Z_minus_X_ci'][0] > 0)
    # full-data fits
    full = {m: z.fit_model(F, st, m) for m in ['H', 'X', 'Z', 'F'] + (['C'] if with_C else [])}
    out['full'] = {m: dict(theta=r['theta'].tolist(), ll=float(r['ll']), rho=r['rho']) for m, r in full.items()}
    b = full['F']['theta'][:3]; w = b/b.sum()
    out['w_hat'] = w.tolist(); out['beta_F'] = float(b.sum())
    out['LR_F_vs_Z'] = float(2*(full['F']['ll']-full['Z']['ll'])); out['tier2a'] = bool(out['LR_F_vs_Z'] <= 5.991)
    out['TV'] = float(z.tv(w, z.LYON)); out['tier2b'] = bool(out['TV'] <= 0.15)
    # tier 3: forward simulation of the whole city from the full-data fits
    obs = z.dispersion_stat(st, Y); out['disp_obs'] = float(obs); out['AR_obs'] = float(Y.sum()/st['N'])
    for m in ['Z', 'H', 'X']:
        sims = [z.simulate(st, m, full[m]['theta'], rng) for _ in range(n_fwd)]
        d = np.array([z.dispersion_stat(st, s) for s in sims]); ar = np.array([s.sum()/st['N'] for s in sims])
        out[f'disp_pred_{m}'] = [float(np.quantile(d, q)) for q in (.025, .5, .975)]
        out[f'AR_pred_{m}'] = [float(np.quantile(ar, q)) for q in (.025, .5, .975)]
        out[f'disp_in_{m}'] = bool(out[f'disp_pred_{m}'][0] <= obs <= out[f'disp_pred_{m}'][2])
    out['tier3'] = out['disp_in_Z']
    out['overall'] = bool(out['tier1'] and out['tier2b'] and out['tier3'])
    return out

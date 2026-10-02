"""T1: real-data analysis, Matsumoto 2014/15 influenza. Run ONLY after PREREG_FREEZE.json exists."""
import sys, os, json, time, numpy as np
sys.path.insert(0, 'src'); import zonecore as z, analysis as an
assert os.path.exists('PREREG_FREEZE.json'), 'pre-registration not frozen'
t0 = time.time(); rng = np.random.default_rng(20261001)
st, Y, info = z.load_outcomes()
res = dict(info=info, n_students=float(st['N']), n_classes=int(len(st['n'])), onset_range=[int(np.nonzero(Y.sum(0))[0].min()), int(np.nonzero(Y.sum(0))[0].max())])
res['primary'] = an.run_tests(st, Y, rng, n_boot=10000, n_fwd=500, with_C=True)
json.dump(res, open('results/matsumoto_primary.json', 'w'), indent=1); print('primary done', time.time()-t0)
# school-bootstrap CI of the free-fit shares
F = z.features(st, Y); XF = z.design(F, 'F', st); S, Yd = F['S'][:, 1:], F['Y'][:, 1:]
schools = np.unique(st['school']); rows = {s: np.nonzero(st['school'] == s)[0] for s in schools}; W = []
for b in range(300):
    idx = np.concatenate([rows[s] for s in rng.choice(schools, len(schools))])
    th, _ = z.fit(XF[idx], S[idx], Yd[idx]); W.append(th[:3]/th[:3].sum())
W = np.array(W); res['w_hat_boot_ci'] = [[float(np.quantile(W[:, k], q)) for q in (.025, .975)] for k in range(3)]
res['TV_boot_ci'] = [float(np.quantile([z.tv(w, z.LYON) for w in W], q)) for q in (.025, .975)]
json.dump(res, open('results/matsumoto_primary.json', 'w'), indent=1)
# pre-declared sensitivity analyses (non-decisional)
sens = {}
def quick(tag, w=z.LYON, zmodel='Z'):
    F = z.features(st, Y); o = {}
    for m in ('H', 'X', zmodel, 'F'):
        r = z.fit_model(F, st, m, w=w); o[m] = dict(ll=float(r['ll']), theta=r['theta'].tolist())
    b = np.array(o['F']['theta'][:3]); o['w_hat'] = (b/b.sum()).tolist(); o['TV'] = float(z.tv(b/b.sum(), w)); o['LR_F_vs_Z'] = 2*(o['F']['ll']-o[zmodel]['ll'])
    o['dLL_Z_minus_H'] = o[zmodel]['ll']-o['H']['ll']; o['dLL_Z_minus_X'] = o[zmodel]['ll']-o['X']['ll']; sens[tag] = o
quick('base_full_data')
quick('lyon_day1', w=np.array([0.7557716936070056, 0.11143608736284655, 0.1327922190301478]))
quick('lyon_day2', w=np.array([0.7677978456305176, 0.09796711751842553, 0.13423503685105695]))
quick('density_dependent_transfer', zmodel='Zdd')
K0 = z.K.copy(); z.K = z.kernel(2.5, 1.5); quick('kernel_mean2.5_sd1.5'); z.K = K0
B0 = z.BREAK.copy(); z.BREAK = np.array([], int); quick('no_break_switch'); z.BREAK = B0
res['sensitivity'] = sens
json.dump(res, open('results/matsumoto_primary.json', 'w'), indent=1); print('all done', time.time()-t0)

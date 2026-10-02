"""Run the frozen pipeline on real datasets. usage: 02_run.py TAG ds1,ds2 model1,model2 [fixmodule]
Checkpoints one JSON per dataset in out/TAG_<ds>.json (skips if present)."""
import sys, os, importlib, time
import a2lib as A, pipeline as PL
tag, dss, models = sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(',')
calfun = importlib.import_module(sys.argv[4]).calfun if len(sys.argv) > 4 else None
for n in dss:
    path = f'{A.ROOT}/out/{tag}_{n}.json'
    if os.path.exists(path):
        print('skip', n); continue
    t0 = time.time(); print('==', n, flush=True)
    ds = A.load(n)
    res = PL.run_dataset(ds, models, G=5, G_epi=10, nrep_real=4000, nrep_model=400, seed=20261001, calfun=calfun)
    res['dataset'] = n
    PL.dump(res, path)
    print(n, f'{time.time()-t0:.0f}s', {m: res['models'][m]['score'] for m in res['models']}, flush=True)

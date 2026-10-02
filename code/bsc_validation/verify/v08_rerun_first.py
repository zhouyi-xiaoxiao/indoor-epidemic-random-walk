"""Re-run the first analysis's predictive code (current src) WITHOUT writing into ../data, compare with stored 02_predictive.json."""
import sys, os, json, importlib.util, numpy as np
os.environ['TMPDIR'] = os.path.abspath('.'); os.environ['MPLCONFIGDIR'] = os.path.abspath('.')
sys.path.insert(0, '../scripts')
import _common as C
saved = {}
C.save_json = lambda name, obj: saved.__setitem__(name, obj)
spec = importlib.util.spec_from_file_location('p02', '../scripts/02_power.py'); p02 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p02)
pred = p02.predictive_all()
old = json.load(open('../data/02_predictive.json'))
mx = 0
for ev in ('T1','T2','T5','C1'):
    for m in ('M0','M1','M2'):
        a = np.array(pred['events'][ev][m]['pmf']); b = np.array(old['events'][ev][m]['pmf']); mx = max(mx, np.abs(a-b).max())
print('max abs diff of predictive pmfs, current code vs stored file:', mx)
pw = p02.power(pred); oldpw = json.load(open('../data/02_power.json'))
print('power M2 truth M2', pw['M2']['M2']['P_A1'], oldpw['M2']['M2']['P_A1'], pw['M2']['M0']['P_A1'], oldpw['M2']['M0']['P_A1'])

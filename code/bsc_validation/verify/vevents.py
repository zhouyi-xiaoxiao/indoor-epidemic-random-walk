"""The re-check's own event builders, written from PREREGISTRATION.md section 4 (not from events.py)."""
import numpy as np
from vlat import Lat
_cache = {}
def L(nx, ny, ax, ay):
    key = (nx, ny, round(ax,6), round(ay,6))
    if key not in _cache: _cache[key] = Lat(nx, ny, ax, ay)
    return _cache[key]

def T1(rng):
    lat = L(6, 15, 2.5/6, 0.75)
    seats = [(c, r) for r in range(15) for c in (0,1,2,4,5)]
    idx = (1, 7)
    zoneA = [s for s in seats if 5 <= s[1] <= 9 and s != idx]       # rows 6-10: 24 seats
    zoneB = [s for s in seats if s[1] in (4, 10)]                    # rows 5, 11: 10 seats
    zoneC = [s for s in seats if s[1] < 4 or s[1] > 10]              # 40 seats
    a = [zoneA[i] for i in rng.choice(24, 23, replace=False)]
    c = [zoneC[i] for i in rng.choice(40, 34, replace=False)]
    occ = a + zoneB + c
    near = np.array([True]*33 + [False]*34)
    alt = np.array([True]*23 + [False]*44)
    rec = np.array([lat.site(x, y) for x, y in occ])
    kap = np.exp(rng.uniform(np.log(0.5), np.log(5))) + 0.93
    return dict(lat=lat, src=int(lat.site(*idx)), rec=rec, near=near, alt=alt, segs=[50/60, 50/60], kappa=kap, K=23)

def T2(rng, minutes=200):
    lat = L(5, 13, 0.5, 11.4/13)
    idx = (4, 11)
    drv = [(c, r) for r in range(13) for c in (3, 4) if (c, r) != idx and (c, r) != (3, 7)]
    oppF = [(c, r) for r in range(7) for c in (0, 1)]
    oppR = [(c, r) for r in range(7, 13) for c in (0, 1)]
    oF = [oppF[i] for i in rng.choice(14, 12, replace=False)]
    oR = [oppR[i] for i in rng.choice(12, 8, replace=False)]
    occ = drv + oF + oR + [(2, 12)]
    near = np.array([r >= 7 for (c, r) in occ])
    side = np.array([c >= 3 for (c, r) in occ]); valid = np.array([(c, r) != (2, 12) for (c, r) in occ])
    assert near.sum() == 19 and len(occ) == 45
    rec = np.array([lat.site(x, y) for x, y in occ])
    kap = 4.8*np.exp(0.2*rng.standard_normal()) + 0.93
    return dict(lat=lat, src=int(lat.site(*idx)), rec=rec, near=near, side=side, valid=valid, segs=[minutes/60], kappa=kap, K=7)

def T5(rng, index=None):
    lat = L(6, 46, 0.9, 1.1)
    seats = [(c, r+2) for r in range(7) for c in (0, 2, 3, 5)]
    i0 = rng.integers(28) if index is None else index
    others = [s for i, s in enumerate(seats) if i != i0]
    occ = [others[i] for i in rng.choice(27, 20, replace=False)]
    d = np.array([np.hypot((x-seats[i0][0])*0.9, (y-seats[i0][1])*1.1) for x, y in occ]) + 1e-6*rng.random(20)
    near = np.zeros(20, bool); near[np.argsort(d)[:12]] = True
    rec = np.array([lat.site(x, y) for x, y in occ])
    kap = np.exp(rng.uniform(np.log(10), np.log(30))) + 0.93
    return dict(lat=lat, src=int(lat.site(*seats[i0])), rec=rec, near=near, segs=[10.0], kappa=kap, K=12, dist=d)

def C1(rng, corner=False):
    lat = L(11, 13, 1.0, 1.0)
    front = [(x, y) for y in (2, 4) for x in (1, 3, 5, 7, 9)]
    back = [(x, y) for y in (6, 8, 10) for x in (1, 3, 5, 7, 9)]
    b = [back[i] for i in rng.choice(15, 14, replace=False)]
    occ = front + b
    near = np.array([True]*10 + [False]*14)
    rec = np.array([lat.site(x, y) for x, y in occ])
    kap = np.exp(rng.uniform(np.log(2), np.log(12))) + 0.93
    src = lat.site(0, 0) if corner else lat.site(5, 0)
    return dict(lat=lat, src=int(src), rec=rec, near=near, segs=[6.5, 6.5], kappa=kap, K=12)
B = dict(T1=T1, T2=T2, T5=T5, C1=C1)
OBS = dict(T1=(14, 33, 9, 34), T2=(3, 19, 4, 26), T5=(11, 12, 1, 8), C1=(8, 10, 4, 14))   # retyped by the re-check from sources

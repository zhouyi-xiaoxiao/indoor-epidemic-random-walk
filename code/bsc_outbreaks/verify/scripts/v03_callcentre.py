"""Seoul call centre (Park et al. 2020, Fig. 1 and Fig. 2): re-check.

The re-check read the published floor plan by eye from 3x enlarged crops (verify/src/crop_*.png) and typed
the seat map below, desk by desk. Nothing here uses the first analysis's image-segmentation code; their
digitised CSV is read only to compare.

Outputs: results/callcentre_check.json, results/callcentre_manual_seatmap.csv,
         results/callcentre_joincount_null_rookdiag.npy
Seed 771001.
"""
import json, math, csv
from pathlib import Path
import numpy as np
from scipy import stats

V = Path(__file__).resolve().parents[1]
ROOT = V.parent
RES = V / "results"
rng = np.random.default_rng(771001)
out = {}

# ---------------------------------------------------------------- manual seat map (north wing)
# left single column: 7 desks, top to bottom
LEFT = [0, 0, 0, 1, 1, 0, 0]
# islands 1..10 west -> east: head desk, then six rows of (left desk, right desk), top to bottom
ISL = [
    (1, [(1, 0), (0, 1), (1, 1), (1, 1), (0, 0), (1, 0)]),
    (1, [(0, 0), (0, 1), (1, 0), (0, 0), (0, 1), (1, 1)]),
    (1, [(1, 1), (0, 0), (0, 1), (1, 1), (0, 1), (0, 1)]),
    (0, [(0, 0), (1, 1), (1, 1), (0, 1), (1, 1), (1, 1)]),
    (1, [(1, 1), (0, 0), (0, 1), (0, 1), (0, 1), (0, 0)]),
    (1, [(0, 1), (0, 1), (0, 1), (1, 1), (1, 1), (1, 0)]),
    (1, [(1, 1), (1, 0), (1, 1), (0, 1), (1, 0), (1, 0)]),
    (0, [(1, 1), (0, 0), (1, 1), (0, 0), (1, 1), (0, 0)]),
    (1, [(1, 0), (0, 1), (1, 1), (1, 0), (0, 0), (1, 1)]),
    (1, [(0, 0), (1, 0), (1, 1), (0, 1), (1, 1), (1, 0)]),
]
# south wing: left single column of 6, islands of 2x5 + foot desk (11), the middle island 2x6 (12)
SOUTH_DESKS = 6 + 11 + 11 + 12 + 11 + 11
SOUTH_CASES = 2 + 0 + 1 + 1 + 0 + 0
EAST_OFFICE_CASES = 1

desks = []   # dict(id, island, col, row, case, x_px, y_px)
X0, DX, Y1, DY = 97.0, 82.3, 77.0, 32.5      # approximate pixel geometry of the 900 px figure
for r, c in enumerate(LEFT):
    desks.append(dict(island=0, col="S", row=r, case=c, x=23.0, y=47.0 + DY * r))
for i, (head, rows) in enumerate(ISL, start=1):
    X = X0 + DX * (i - 1)
    desks.append(dict(island=i, col="H", row=0, case=head, x=X, y=53.5))
    for r, (l, rr) in enumerate(rows, start=1):
        desks.append(dict(island=i, col="L", row=r, case=l, x=X - 9, y=Y1 + DY * (r - 1)))
        desks.append(dict(island=i, col="R", row=r, case=rr, x=X + 9, y=Y1 + DY * (r - 1)))
for j, d in enumerate(desks): d["id"] = j
n = len(desks); case = np.array([d["case"] for d in desks], bool); k = int(case.sum())
with open(RES / "callcentre_manual_seatmap.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "island", "col", "row", "case", "x", "y"]); w.writeheader(); w.writerows(desks)
out["manual_counts"] = dict(north_desks=n, north_cases=k, south_desks=SOUTH_DESKS, south_cases=SOUTH_CASES,
                            east_office_cases=EAST_OFFICE_CASES, total_case_seats=k + SOUTH_CASES + EAST_OFFICE_CASES,
                            cases_per_island=[int(sum(d["case"] for d in desks if d["island"] == i)) for i in range(11)],
                            desks_per_island=[sum(1 for d in desks if d["island"] == i) for i in range(11)])

# ---------------------------------------------------------------- compare with the first analysis's CSV
th = list(csv.DictReader(open(ROOT / "data" / "park2020_callcentre_seats_digitized.csv")))
thn = [t for t in th if t["region"] == "north_wing"]; ths = [t for t in th if t["region"] == "south_wing"]
cmp_ = dict(their_north=[sum(int(t["case"]) for t in thn), len(thn)], their_south=[sum(int(t["case"]) for t in ths), len(ths)],
            their_other=[sum(int(t["case"]) for t in th if t["region"] not in ("north_wing", "south_wing")),
                         sum(1 for t in th if t["region"] not in ("north_wing", "south_wing"))])
# nearest-desk matching between the grid of the re-check and the detections of the first analysis
check_xy = np.array([[d["x"], d["y"]] for d in desks]); their_xy = np.array([[float(t["x_px"]), float(t["y_px"])] for t in thn])
D = np.sqrt(((check_xy[:, None, :] - their_xy[None, :, :]) ** 2).sum(-1))
nn = D.argmin(1)
cmp_["max_match_distance_px"] = float(D.min(1).max())
cmp_["matching_is_bijection"] = bool(len(set(nn.tolist())) == n == len(thn))
their_case = np.array([int(thn[j]["case"]) for j in nn], bool)
cmp_["desks_with_different_case_flag"] = int((their_case != case).sum())
out["comparison_with_first_analysis"] = cmp_

# ---------------------------------------------------------------- wing contrast
def katz(k1, n1, k0, n0):
    rr = (k1 / n1) / (k0 / n0); se = math.sqrt(1 / k1 - 1 / n1 + 1 / k0 - 1 / n0)
    return rr, rr * math.exp(-1.959964 * se), rr * math.exp(1.959964 * se)


rr = katz(k, n, SOUTH_CASES, SOUTH_DESKS)
out["wing_contrast"] = dict(north=f"{k}/{n}", south=f"{SOUTH_CASES}/{SOUTH_DESKS}", rr=rr[0], lo=rr[1], hi=rr[2],
                            fisher_p=float(stats.fisher_exact([[k, n - k], [SOUTH_CASES, SOUTH_DESKS - SOUTH_CASES]])[1]),
                            north_ar=k / n, south_ar=SOUTH_CASES / SOUTH_DESKS)
# sensitivity: the 10 case seats missing from the figure. Worst case for the contrast: all 10 in the south wing
for lab, kn, ks in (("all10_south", k, SOUTH_CASES + 10), ("all10_north", k + 10, SOUTH_CASES), ("5_5", k + 5, SOUTH_CASES + 5)):
    r_ = katz(kn, n, ks, SOUTH_DESKS)
    out["wing_contrast"][f"sens_{lab}"] = dict(rr=r_[0], lo=r_[1], hi=r_[2],
        fisher_p=float(stats.fisher_exact([[kn, n - kn], [ks, SOUTH_DESKS - ks]])[1]))

# ---------------------------------------------------------------- adjacency (structural, three definitions)
idx = {(d["island"], d["col"], d["row"]): d["id"] for d in desks}
rook, diag, aisle = set(), set(), set()
def add(S, a, b):
    if a in idx and b in idx: S.add(tuple(sorted((idx[a], idx[b]))))
for r in range(6): add(rook, (0, "S", r), (0, "S", r + 1))
for i in range(1, 11):
    add(rook, (i, "H", 0), (i, "L", 1)); add(rook, (i, "H", 0), (i, "R", 1))
    for r in range(1, 7):
        add(rook, (i, "L", r), (i, "R", r))                       # facing across the desk
        add(rook, (i, "L", r), (i, "L", r + 1)); add(rook, (i, "R", r), (i, "R", r + 1))   # side by side
        add(diag, (i, "L", r), (i, "R", r + 1)); add(diag, (i, "R", r), (i, "L", r + 1))
        if i < 10: add(aisle, (i, "R", r), (i + 1, "L", r))       # back to back across the aisle
        if i == 1: add(aisle, (0, "S", r), (1, "L", r))
ADJ = {"rook": rook, "rook+diag": rook | diag, "rook+diag+aisle": rook | diag | aisle}
# geometric definition equivalent to the first analysis's (centres within 1.25 pitches)
xy = check_xy / 32.5
Dm = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
for rad in (1.0, 1.25, 1.5, 2.0):
    ADJ[f"radius{rad}"] = {(a, b) for a in range(n) for b in range(a + 1, n) if Dm[a, b] <= rad}

NPERM = 20000
perms = np.array([rng.permutation(case) for _ in range(NPERM)])
jc = {}
for name, S in ADJ.items():
    P = np.array(sorted(S)); a, b = P[:, 0], P[:, 1]
    obs = int((case[a] & case[b]).sum())
    null = (perms[:, a] & perms[:, b]).sum(1)
    E = len(P) * k * (k - 1) / (n * (n - 1))
    jc[name] = dict(pairs=int(len(P)), observed=obs, analytic_mean=E, null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                    z=float((obs - null.mean()) / null.std(ddof=1)), p_greater=float((np.sum(null >= obs) + 1) / (NPERM + 1)),
                    p_less=float((np.sum(null <= obs) + 1) / (NPERM + 1)))
    if name == "rook+diag": np.save(RES / "callcentre_joincount_null_rookdiag.npy", null.astype(np.int16))
out["joincount"] = jc

# island-level dispersion (cases per island, 13 desks each): observed variance vs permutation null
isl = np.array([d["island"] for d in desks])
def islvar(c): return np.var([c[isl == i].sum() for i in range(1, 11)], ddof=1)
v_obs = islvar(case); v_null = np.array([islvar(p) for p in perms[:5000]])
out["island_dispersion"] = dict(var_obs=float(v_obs), null_mean=float(v_null.mean()),
                                p_greater=float((np.sum(v_null >= v_obs) + 1) / 5001), p_less=float((np.sum(v_null <= v_obs) + 1) / 5001))
# gradients
xs = np.array([d["x"] for d in desks]); ys = np.array([d["y"] for d in desks])
r1 = stats.spearmanr(xs, case); r2 = stats.spearmanr(ys, case)
out["gradients"] = dict(west_east_rho=float(r1[0]), west_east_p=float(r1[1]), north_south_rho=float(r2[0]), north_south_p=float(r2[1]))

# ---------------------------------------------------------------- power of the join-count test
# Contagion alternative on the desk graph: desks are infected one at a time; the next case is drawn with
# probability proportional to (eps + number of already-infected graph neighbours). eps -> infinity is
# spatially random; eps -> 0 is pure neighbour-to-neighbour spread. Same k = 79 of n = 137.
P = np.array(sorted(ADJ["rook+diag"])); A = np.zeros((n, n), bool); A[P[:, 0], P[:, 1]] = True; A |= A.T
null_rd = (perms[:, P[:, 0]] & perms[:, P[:, 1]]).sum(1)
crit = np.quantile(null_rd, 0.95)
def simulate(eps):
    inf = np.zeros(n, bool); nb = np.zeros(n)
    for _ in range(k):
        w = np.where(inf, 0.0, eps + nb); j = rng.choice(n, p=w / w.sum())
        inf[j] = True; nb += A[j]
    return int((inf[P[:, 0]] & inf[P[:, 1]]).sum())
power = {}
for eps in (0.05, 0.2, 0.5, 1.0, 2.0, 5.0):
    sims = np.array([simulate(eps) for _ in range(400)])
    # share of infections that are 'local' at mid-epidemic ~ mean nb/(eps+mean nb); report simple summary
    power[str(eps)] = dict(mean_joins=float(sims.mean()), sd=float(sims.std(ddof=1)), power_at_alpha05=float((sims > crit).mean()),
                           frac_sims_at_or_below_observed=float((sims <= jc["rook+diag"]["observed"]).mean()))
out["power_contagion_alternative"] = dict(critical_value_95=float(crit), n_sims=400, by_eps=power,
    note="weight of an uninfected desk = eps + (# infected neighbours); eps is the background (non-local) weight")

# ---------------------------------------------------------------- onset curve (re-check's reading of Fig. 1)
ONSET11 = {  # day-of-year style key: (month, day): count on 11th floor
    (2, 25): 1, (2, 28): 3, (2, 29): 1, (3, 1): 1, (3, 2): 2, (3, 4): 3, (3, 5): 7, (3, 6): 16, (3, 7): 14, (3, 8): 12,
    (3, 9): 16, (3, 11): 2, (3, 13): 1, (3, 14): 1, (3, 15): 2, (3, 16): 2, (3, 19): 1, (3, 20): 1}
OTHER = {"floor10": [(2, 22), (2, 29), (3, 12)], "floor9": [(3, 2)]}
def dnum(m, d): return d - 25 if m == 2 else d + 4          # days since 25 Feb 2020 (leap year: Feb has 29 days)
series = np.zeros(25);
for (m, d), c in ONSET11.items(): series[dnum(m, d)] = c
th_epi = list(csv.DictReader(open(ROOT / "data" / "park2020_callcentre_epicurve_digitized.csv")))
their = np.zeros(25)
for t in th_epi:
    y, m, d = map(int, t["onset_date"].split("-"))
    if dnum(m, d) >= 0: their[dnum(m, d)] = int(t["floor11"])
out["onset"] = dict(total_floor11=int(series.sum()), peak_6_to_9_march=int(series[10:14].sum()),
                    identical_to_first_analysis=bool(np.array_equal(series, their)),
                    floor10_cases_in_figure=len(OTHER["floor10"]), floor10_cases_in_table=2)

def poisson_fit(t, y):
    b = np.array([math.log(max(y.mean(), 0.1)), 0.0])
    X = np.column_stack([np.ones_like(t), t])
    for _ in range(200):
        mu = np.exp(X @ b); g = X.T @ (y - mu); H = X.T @ (X * mu[:, None])
        step = np.linalg.solve(H, g); b = b + step
        if np.abs(step).max() < 1e-12: break
    cov = np.linalg.inv(X.T @ (X * np.exp(X @ b)[:, None]))
    mu = np.exp(X @ b)
    dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0) - (y - mu))
    return dict(r=float(b[1]), se=float(math.sqrt(cov[1, 1])), lo=float(b[1] - 1.959964 * math.sqrt(cov[1, 1])),
                hi=float(b[1] + 1.959964 * math.sqrt(cov[1, 1])), n_cases=int(y.sum()), deviance=float(dev), dof=len(y) - 2,
                dispersion=float(np.sum((y - mu) ** 2 / mu) / (len(y) - 2)))
fits = {}
for lab, end in (("25Feb-6Mar", 10), ("25Feb-7Mar", 11), ("25Feb-9Mar", 13)):
    t = np.arange(end + 1, dtype=float); fits[lab] = poisson_fit(t, series[: end + 1])
out["growth"] = fits
# quasi-Poisson widened CI
for lab in list(fits):
    f = fits[lab] if isinstance(fits[lab], dict) else None
    if f:
        s = f["se"] * math.sqrt(max(1.0, f["dispersion"]))
        f["quasi_lo"], f["quasi_hi"] = f["r"] - 1.959964 * s, f["r"] + 1.959964 * s
json.dump(out, open(RES / "callcentre_check.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))

"""Re-computation, with separately written code, of every derived number of the outbreak dataset.

The re-check's own implementation: nothing is imported from ../scripts. Counts are typed here from the
sources re-opened by the re-check (verify/src/*.txt), not copied from
records.py. Where the first analysis is read (outbreaks.json, data/*.csv) it is only to COMPARE.

Outputs (verify/results/):
  stats_check.json            all recomputed quantities + discrepancy list
  attack_rates_check.csv      Clopper-Pearson intervals, re-check vs first analysis
  risk_ratios_check.csv       RR / Katz CI / Fisher p, re-check vs first analysis
  train_matrix_reconstruction.csv   integer (k, n) per cell of Hu et al. Table 1 recovered from CIs
Deterministic: seed 771001 (different from the original 20261001).
"""
import json, math, csv
from pathlib import Path
import numpy as np
from scipy import stats, optimize

V = Path(__file__).resolve().parents[1]
ROOT = V.parent
RES = V / "results"; RES.mkdir(exist_ok=True)
SEED = 771001
out = {"seed": SEED}
disc = []          # discrepancies between the numbers of the re-check and those of the first analysis


def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)


def cp_upper_one_sided(k, n, a=0.05):
    return float(stats.beta.ppf(1 - a, k + 1, n - k))


def katz(k1, n1, k0, n0):
    a, b, c, d = k1, n1 - k1, k0, n0 - k0
    cc = 0.5 if min(a, b, c, d) == 0 else 0.0
    p1 = (a + cc) / (n1 + 2 * cc); p0 = (c + cc) / (n0 + 2 * cc)
    rr = p1 / p0
    se = math.sqrt(1 / (a + cc) - 1 / (n1 + 2 * cc) + 1 / (c + cc) - 1 / (n0 + 2 * cc))
    return rr, rr * math.exp(-1.959964 * se), rr * math.exp(1.959964 * se), cc > 0


# ------------------------------------------------------------------ 1. attack rates
# (id, k, n, T_hours or None, note) -- typed by the re-check from re-opened sources
REC = [
    ("O1_seoul_callcentre", 94, 216, None, "Park 2020 Table 1"),
    ("O2_zurich_office", 8, 12, 11.0, "Weissberg 2020 scenario 2 (67%); scenario 1 10/12"),
    ("O3_england_office", 22, 40, None, "Nicholls 2024 abstract"),
    ("O4_german_meatplant", 20, 78, None, "Guenther 2020 Methods: 20 out of 78 = 25.6%"),
    ("O5_tb_office", 27, 67, 160.0, "Nardell 1991 abstract; 160 h is an assumption (4 wk x 40 h)"),
    ("S1_liaocheng_supermarket", 11, 120, None, "Tian 2021 Results"),
    ("S2_tianjin_dept_store", 6, 194, None, "Li 2022 (6 employees; 194 quarantined staff); Chen 2020 gives 7"),
    ("S4_boston_grocery", 21, 104, None, "Lan 2021 abstract"),
    ("C1_marin_classroom", 12, 24, None, "Lam-Hine 2021"),
    ("C2_jerusalem_highschool", 153, 1161, None, "Stein-Zamir 2020 Table"),
    ("C3_utah_elementary", 5, 728, None, "Hershow 2021"),
    ("C5_cheonan_fitness", 57, 217, None, "Jang 2020: 26.3% of 217 (count derived)"),
    ("T1_zhejiang_bus", 23, 67, 100 / 60, "Shen 2020"),
    ("T2_hunan_coach", 7, 48, 2.5, "Luo 2020 (Ou 2022: 7/46, 200 min)"),
    ("T3_hunan_minibus", 2, 12, 1.0, "Luo 2020 (Ou 2022: 2/17)"),
    ("T4_china_hsr", 234, 72093, None, "Hu 2021"),
    ("T5_vn54_flight", 12, 20, 10.0, "Khanh 2020"),
    ("R1_guangzhou_restaurant", 5, 79, 82 / 60, "Li 2021 Table 2 (16+63 = 11+68 = 79 non-A patrons)"),
    ("R2_jeonju_restaurant", 2, 13, 21 / 60, "Kwon 2020"),
    ("H1_skagit_choir", 52, 60, 2.5, "Hamner 2020 (32 confirmed + 20 probable secondary)"),
    ("H2_sydney_church", 12, 508, 1.0, "Katelaris 2021"),
    ("X1_diamond_princess", 634, 3711, None, "Mizumoto 2020"),
    ("X2_singapore_dormitories", 17758, 295000, None, "Koh 2020"),
]
theirs = {r["id"]: r for r in json.load(open(ROOT / "outbreaks.json"))["records"]}
out["n_records_in_dataset"] = len(theirs)
out["grades"] = {g: sum(1 for r in theirs.values() if r["quality"] == g) for g in "ABCX"}
out["grade_A_ids"] = sorted(r["id"] for r in theirs.values() if r["quality"] == "A")
keys = set()
for r in theirs.values():
    for s in r["sources"]:
        keys.add(s["key"])
out["n_distinct_source_keys_in_records"] = len(keys)

rows = []
for rid, k, n, T, note in REC:
    lo, hi = cp(k, n)
    t = theirs[rid]
    row = dict(id=rid, k=k, n=n, ar=k / n, ci_lo=lo, ci_hi=hi, T_h=T,
               H=-math.log(1 - k / n), their_k=t.get("n_infected"), their_n=t.get("n_exposed"),
               their_ar=t.get("attack_rate"), their_lo=(t.get("attack_rate_ci95_exact") or [None, None])[0],
               their_hi=(t.get("attack_rate_ci95_exact") or [None, None])[1], their_T=t.get("exposure_duration_h"), note=note)
    row["hazard_per_h"] = row["H"] / T if T else None
    if (k, n) != (row["their_k"], row["their_n"]):
        disc.append(f"{rid}: counts differ re-check {k}/{n} first analysis {row['their_k']}/{row['their_n']}")
    if row["their_lo"] is not None and (abs(lo - row["their_lo"]) > 5e-5 or abs(hi - row["their_hi"]) > 5e-5):
        disc.append(f"{rid}: CP interval differs re-check {lo:.5f}-{hi:.5f} first analysis {row['their_lo']}-{row['their_hi']}")
    if T and row["their_T"] and abs(T - row["their_T"]) > 0.02:
        disc.append(f"{rid}: duration differs re-check {T} first analysis {row['their_T']}")
    rows.append(row)
with open(RES / "attack_rates_check.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
AR = {r["id"]: r for r in rows}

# ------------------------------------------------------------------ 2. near/far risk ratios
STRATA = [  # (id, label, k_near, n_near, k_far, n_far, source statement)
    ("O1", "north vs south wing (figure count)", 79, 137, 4, 62, "check's manual desk count of Park Fig. 2"),
    ("C1", "front 2 rows vs back 3 rows", 8, 10, 4, 14, "Lam-Hine: 80% (eight of 10) vs 28% (four of 14), p=0.036"),
    ("T1a", "rows 5-11 vs rows 1-4,12-15", 14, 33, 9, 34, "Shen: RR 1.6 (0.8-3.2)"),
    ("T1b", "rows 6-10 vs rows 1-5,11-15", 11, 23, 12, 44, "Shen: RR 1.8 (0.9-3.3)"),
    ("T2", "rear rows 8-13 vs front rows 1-7", 3, 19, 4, 26, "Ou 2022 Table 1"),
    ("T2side", "driver side vs opposite side", 6, 24, 1, 20, "Ou 2022 Table 1 (25% vs 5%, p=0.164)"),
    ("T5", "<2 seats vs >2 seats", 11, 12, 1, 8, "Khanh: RR 7.3 (1.2-46.2)"),
    ("R1a", "immediate vs remote tables", 5, 16, 0, 63, "Li 2021 Table 2"),
    ("R1b", "ABC zone vs non-ABC zone", 5, 11, 0, 68, "Li 2021 Table 2"),
    ("C2", "grades 7-9 vs 10-12 (persons)", 135, 581, 18, 583, "Stein-Zamir Table (persons)"),
    ("C2t", "grades 7-9 vs 10-12 (tested)", 135, 581, 18, 580, "Stein-Zamir Table (tested: 200+194+186)"),
    ("S1", "staff vs customers", 11, 120, 0, 8224, "Tian 2021"),
]
rr_rows = []
for sid, lab, k1, n1, k0, n0, src in STRATA:
    rr, lo, hi, cc = katz(k1, n1, k0, n0)
    p = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])[1]
    rr_rows.append(dict(id=sid, contrast=lab, near=f"{k1}/{n1}", far=f"{k0}/{n0}", rr=rr, lo=lo, hi=hi,
                        continuity=cc, fisher_p=float(p), source=src,
                        near_ci=cp(k1, n1), far_ci=cp(k0, n0)))
out["risk_ratios"] = rr_rows
their_rr = list(csv.DictReader(open(ROOT / "data" / "near_far_risk_ratios.csv")))
check_by_counts = {(r["near"], r["far"]): r for r in rr_rows}
for t in their_rr:
    key = (f"{t['near_k']}/{t['near_n']}", f"{t['far_k']}/{t['far_n']}")
    m = check_by_counts.get(key)
    if m is None:
        disc.append(f"RR row not reproduced (counts {key})"); continue
    for a, b in (("rr", "risk_ratio"), ("lo", "rr_ci_low"), ("hi", "rr_ci_high"), ("fisher_p", "fisher_p")):
        if abs(m[a] - float(t[b])) > 1e-3 * max(1, abs(m[a])):
            disc.append(f"RR {t['outbreak_id']} {key}: {a} re-check {m[a]:.5g} first analysis {float(t[b]):.5g}")
with open(RES / "risk_ratios_check.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rr_rows[0].keys())); w.writeheader(); w.writerows(rr_rows)

# ------------------------------------------------------------------ 3. single-exposure events with a duration
NC1_EVENTS = [("T1_zhejiang_bus", "metro"), ("T2_hunan_coach", "metro"), ("T3_hunan_minibus", "metro"),
              ("T5_vn54_flight", "metro"), ("R1_guangzhou_restaurant", "classroom"),
              ("R2_jeonju_restaurant", "classroom"), ("H1_skagit_choir", "classroom"),
              ("H2_sydney_church", "classroom"), ("O2_zurich_office", "office")]

# ------------------------------------------------------------------ 6. per-hour hazard heterogeneity
hz = {rid: AR[rid]["hazard_per_h"] for rid, _ in NC1_EVENTS}
v = np.array(list(hz.values()))
lg = np.log(v)
out["hazard_per_hour"] = dict(values=hz, min=float(v.min()), max=float(v.max()), ratio=float(v.max() / v.min()),
                              geometric_mean=float(np.exp(lg.mean())), gsd=float(np.exp(lg.std(ddof=1))),
                              gsd_ddof0=float(np.exp(lg.std(ddof=0))))
hz2 = dict(hz); hz2["H1_skagit_choir"] = -math.log(1 - 32 / 60) / 2.5
v2 = np.array(list(hz2.values()))
out["hazard_per_hour"]["confirmed_only_choir"] = dict(H1=hz2["H1_skagit_choir"], max=float(v2.max()),
    argmax=max(hz2, key=hz2.get), ratio=float(v2.max() / v2.min()), geometric_mean=float(np.exp(np.log(v2).mean())),
    gsd=float(np.exp(np.log(v2).std(ddof=1))))
out["hazard_per_hour"]["T2_over_T3"] = hz["T3_hunan_minibus"] / hz["T2_hunan_coach"]
# CI of each per-hour hazard (transforming the CP interval)
out["hazard_per_hour"]["ci"] = {rid: [-math.log(1 - AR[rid]["ci_lo"]) / AR[rid]["T_h"],
                                       -math.log(1 - AR[rid]["ci_hi"]) / AR[rid]["T_h"]] for rid, _ in NC1_EVENTS}

# ------------------------------------------------------------------ 7. negative controls
NEG = [("S1 customers", 0, 8224), ("T1 unexposed bus", 0, 60), ("R1 remote tables", 0, 63), ("R1 other zones", 0, 68),
       ("T5 premium economy", 0, 35), ("C3 Utah classroom contacts", 5, 728), ("O1 floors 7-9", 1, 595)]
out["negative_controls"] = [dict(name=a, k=k, n=n, upper95_one_sided=cp_upper_one_sided(k, n)) for a, k, n in NEG]

# ------------------------------------------------------------------ 8. Hu et al. train matrix
# values as returned twice by two separate page extractions (first analysis + re-check); 'printed' layout
printed = {  # row -> list under headers [same col, 1, 2, 3, 4, 5]
    0: [(3.53, 2.89, 4.31), (1.65, 1.18, 2.31), (0.38, 0.18, 0.78), (0.38, 0.19, 0.79), (0.29, 0.10, 0.85), None],
    1: [(0.21, 0.11, 0.38), (0.24, 0.14, 0.41), (0.14, 0.06, 0.32), (0.09, 0.03, 0.25), (0.03, 0.00, 0.16), (0.05, 0.00, 0.30)],
    2: [(0.25, 0.14, 0.45), (0.17, 0.09, 0.33), (0.23, 0.12, 0.46), (0.16, 0.07, 0.36), (0.09, 0.03, 0.27), (0.17, 0.06, 0.50)],
    3: [(0.05, 0.01, 0.18), (0.05, 0.01, 0.17), (0.13, 0.05, 0.33), (0.10, 0.03, 0.30), (0.10, 0.03, 0.30), (0.06, 0.00, 0.36)],
}
row_means = {0: (1.53, 1.30, 1.80), 1: (0.14, 0.10, 0.20), 2: (0.18, 0.13, 0.25), 3: (0.08, 0.05, 0.13)}
col_means = {0: (0.17, 0.12, 0.26), 1: (0.68, 0.56, 0.81), 2: (0.41, 0.31, 0.54), 3: (0.16, 0.10, 0.25),
             4: (0.12, 0.07, 0.20), 5: (0.13, 0.06, 0.25)}
overall = (0.32, 0.28, 0.36)


def wilson(k, n, z=1.959964):
    k = np.asarray(k, float); n = np.asarray(n, float)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def recover(p, lo, hi, kmax=320):
    """Best integer (k, n): rate rounds to the printed percentage and the Wilson 95% interval is closest to
    the printed one (sum of absolute errors in percentage points). Wilson reproduces the large cells to
    <= 0.005 points (results/train_ci_method_exploration.json); for k = 1 cells the printed lower limit 0.00
    is matched by Clopper-Pearson instead, which changes n by < 2%."""
    best = None
    for k in range(1, kmax):
        n_lo = int(math.floor(k / ((p + 0.005) / 100))); n_hi = int(math.ceil(k / (max(p - 0.005, 1e-4) / 100)))
        n = np.arange(max(n_lo, k + 1), n_hi + 1, dtype=float)
        if len(n) == 0 or len(n) > 40000: continue
        ok = np.round(100 * k / n, 2) == round(p, 2)
        if not ok.any(): continue
        n = n[ok]; l, h = wilson(k, n)
        err = np.abs(100 * l - lo) + np.abs(100 * h - hi); j = int(err.argmin())
        if best is None or err[j] < best[0]: best = (float(err[j]), k, int(n[j]))
    return best


cellfit = {}
rec_rows = []
for r in range(4):
    for c in range(6):
        cell = printed[r][c]
        if cell is None: continue
        e, k_, n_ = recover(*cell)
        cellfit[(r, c)] = (k_, n_)
        rec_rows.append(dict(row=r, printed_col=c, pct=cell[0], lo=cell[1], hi=cell[2], k=k_, n=n_, abs_err_pp=e))
agg = {}
for name, d in (("row_mean", row_means), ("col_mean", col_means)):
    for i, cell in d.items():
        e, k_, n_ = recover(*cell)
        agg[f"{name}_{i}"] = dict(pct=cell[0], k=k_, n=n_, abs_err_pp=e)
e, k_, n_ = recover(*overall)
agg["overall_table"] = dict(k=k_, n=n_, abs_err_pp=e, wilson_234_72093=[float(100 * x) for x in wilson(234, 72093)],
                            cp_234_72093=[100 * x for x in cp(234, 72093)])
with open(RES / "train_matrix_reconstruction.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rec_rows[0].keys())); w.writeheader(); w.writerows(rec_rows)


def colmeans(shift):
    """Pooled column rates when the five printed same-row cells are placed at columns shift..shift+4."""
    res = {}
    for c in range(6):
        K = N = 0.0
        for r in range(4):
            pc = c - shift if r == 0 else c
            if (r, pc) not in cellfit: continue
            K += cellfit[(r, pc)][0]; N += cellfit[(r, pc)][1]
        res[c] = 100 * K / N if N else None
    return res


A, B = colmeans(0), colmeans(1)
out["train_matrix"] = dict(
    cells=rec_rows, aggregates=agg,
    sum_k=sum(v[0] for v in cellfit.values()), sum_n=sum(v[1] for v in cellfit.values()),
    reported_total=[234, 72093],
    row_sums={r: [sum(cellfit[(r, c)][0] for c in range(6) if (r, c) in cellfit),
                  sum(cellfit[(r, c)][1] for c in range(6) if (r, c) in cellfit)] for r in range(4)},
    column_rates_printed_layout=A, column_rates_row0_shifted=B,
    printed_column_means={c: v[0] for c, v in col_means.items()},
    abs_err_printed_layout=sum(abs(A[c] - col_means[c][0]) for c in range(6) if A[c] is not None),
    abs_err_shifted=sum(abs(B[c] - col_means[c][0]) for c in range(6) if B[c] is not None),
    row_rates_from_cells={r: 100 * sum(cellfit[(r, c)][0] for c in range(6) if (r, c) in cellfit) /
                             sum(cellfit[(r, c)][1] for c in range(6) if (r, c) in cellfit) for r in range(4)},
    printed_row_means={r: v[0] for r, v in row_means.items()},
    adjacent_vs_far_ratio=dict(adjacent=3.53, rows1to3_min=0.03, rows1to3_max=0.25, rows2to3_mean=[0.18, 0.08]))

# ------------------------------------------------------------------ 9. Guangzhou restaurant table sums
tab = list(csv.DictReader(open(ROOT / "data" / "li2021_restaurant_tables.csv")))
out["restaurant"] = dict(
    patrons=sum(int(t["patrons"]) for t in tab), infected=sum(int(t["infected"]) for t in tab),
    nonA_patrons=sum(int(t["patrons"]) for t in tab if t["table"] != "TA"),
    nonA_infected=sum(int(t["infected"]) for t in tab if t["table"] != "TA"),
    immediate=[sum(int(t["infected"]) for t in tab if t["neighbour_class"] == "immediate"),
               sum(int(t["patrons"]) for t in tab if t["neighbour_class"] == "immediate")],
    abc_nonA=[sum(int(t["infected"]) for t in tab if t["zone"] == "ABC" and t["table"] != "TA"),
              sum(int(t["patrons"]) for t in tab if t["zone"] == "ABC" and t["table"] != "TA")],
    tracer_outside_zone=[float(t["norm_measured_tracer"]) for t in tab if t["zone"] != "ABC" and t["norm_measured_tracer"]],
    area_17x8_1=17 * 8.1, ach_to_Ls_per_patron=[x * 431 / 3.6 / 89 for x in (0.56, 0.77)])

out["discrepancies"] = disc
json.dump(out, open(RES / "stats_check.json", "w"), indent=1, default=float)
print("discrepancies:", len(disc))
for d in disc: print("  ", d)
print(json.dumps({k: out[k] for k in ("grades", "n_distinct_source_keys_in_records", "hazard_per_hour")}, indent=1, default=float))

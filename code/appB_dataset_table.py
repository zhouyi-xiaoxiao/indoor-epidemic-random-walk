"""appB_dataset_table.py -- tables and numbers of Supplementary Section S6 ("The outbreak records").

Reads the outbreak dataset and the re-check's results; nothing is typed
by hand except short English labels, and every number that occurs in such a label is checked against the
corresponding dataset fields (see ``_check_tokens``).

Inputs (relative to the repository root):
  data/bsc_outbreaks/outbreaks.csv                              27 records, 52 columns
  data/bsc_outbreaks/tables/near_far_risk_ratios.csv              near/far contrasts
  data/bsc_outbreaks/tables/analysis_results.json                 negative controls, per-hour hazards
  data/bsc_outbreaks/tables/park2020_callcentre_seats_digitized.csv   digitised seat map (Park et al., Fig. 2)
  code/bsc_outbreaks/verify/results/train_matrix_reconstruction.csv  re-check's recovered (k, n) per cell
  code/bsc_outbreaks/verify/results/callcentre_check.json       re-check's desk-by-desk reading
  code/bsc_outbreaks/verify/results/stats_check.json            re-check's recomputation (hazards)
  code/bsc_outbreaks/verify/results/extra_checks.json           Cochran Q, restaurant sensitivity

What is recomputed here (and asserted to agree with the stored values):
  * exact (Clopper-Pearson) 95 % intervals of all attack rates;
  * risk ratios with the Katz log-normal interval (0.5 added to every cell when one cell is zero) and
    Fisher's exact test;
  * one-sided exact 95 % upper bounds of the negative controls;
  * Wilson intervals of the re-check's recovered train cells against the percentages printed by Hu et al.;
  * wing counts of the digitised call-centre seat map and the wing contrast under the two extreme
    placements of the ten unmapped cases;
  * Cochran's Q for the four single-index near/far contrasts.

Outputs:
  ../data/appB_dataset_tables.tex     bodies of Tables appB_dataset:tab:records, :tab:strata, :tab:controls
  ../data/appB_dataset_numbers.json   every number quoted in the text of Supplementary Section S6
  ../sections/supp_datasets.tex       the three table bodies are injected between the marker comments
                                      "% >>> GENERATED <name>" and "% <<< GENERATED <name>" (if present)

    python appB_dataset_table.py            # a few seconds, < 200 MB
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                      # repository root
OB = ROOT / "code/bsc_outbreaks"
VER = OB / "verify/results"
OUT = HERE.parent / "data"
SEC = HERE.parent / "sections" / "supp_datasets.tex"
OUT.mkdir(exist_ok=True)

# ------------------------------------------------------------------------------------------------ statistics


def clopper_pearson(k: int, n: int, level: float = 0.95):
    a = 1 - level
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)


def upper_one_sided(k: int, n: int, level: float = 0.95):
    return float(stats.beta.ppf(level, k + 1, n - k))


def wilson(k: int, n: int, z: float = 1.959963984540054):
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return mid - half, mid + half


def risk_ratio(k1, n1, k0, n0):
    """Risk ratio, Katz 95 % interval; 0.5 is added to all four cells when a cell is zero."""
    cc = (k1 == 0) or (k0 == 0)
    a, b, c, d = k1, n1 - k1, k0, n0 - k0
    if cc:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    rr = (a / (a + b)) / (c / (c + d))
    se = math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
    p = float(stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])[1])
    return rr, rr * math.exp(-1.96 * se), rr * math.exp(1.96 * se), cc, p, math.log(rr), se


# ------------------------------------------------------------------------------------------------ formatting


def num(x: float) -> str:
    return f"{int(round(x)):,}".replace(",", "{,}")


def pct(x: float, small: bool) -> str:
    return f"{100 * x:.2f}" if small else f"{100 * x:.1f}"


def rrfmt(v: float) -> str:
    return num(v) if v >= 100 else f"{v:.1f}"


def pfmt(p: float) -> str:
    if p >= 0.995:
        return "1.0"
    if p >= 0.1:
        return f"{p:.2f}"
    if p >= 0.01:
        return f"{p:.2g}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


_NUMTOK = re.compile(r"\d+(?:\.\d+)?")


def _check_tokens(label: str, haystack: str, where: str):
    """Every number in a hand-written short label must occur in the dataset fields it abbreviates."""
    plain = re.sub(r"\$[\^_][^$]*\$", "", label)             # m$^2$, CO$_2$, h$^{-1}$ carry no data
    plain = plain.replace("{,}", ",").replace("--", "-")
    hay = haystack.replace(",", "")
    for tok in _NUMTOK.findall(plain.replace(",", "")):
        if not re.search(r"(?<![\d.])" + re.escape(tok) + r"(?!\d)", hay):
            # allow "11" for 11.0 and "0.9" inside e.g. "(0.9 m)"
            if not re.search(r"(?<![\d.])" + re.escape(tok) + r"(?:\.0)?(?![\d])", hay):
                raise AssertionError(f"{where}: number {tok!r} of label {label!r} not found in the dataset fields")


# ------------------------------------------------------------------------------------------------ labels
# Short English labels: (setting, exposure, ventilation).  Numbers are checked against the record's fields.
SHORT = {
    "O1": ("Call centre, 11th floor, open-plan desks (Seoul)", "working days over about 2 weeks", "not reported"),
    "O2": ("Team of 13 in an open-plan office (Zurich)", "11 h", "mechanical, about 1 air change per hour"),
    "O3": ("Public-facing office, three open-plan floors (England)", "weeks", r"CO$_2$ typically $\le$1000 ppm"),
    "O4": ("Beef-processing line, fixed work stations (Germany)", "3 shifts",
           "$<$1 air change per hour, recirculating cooling fans"),
    "O5": ("Office building, tuberculosis (Massachusetts)", "160 h (4 weeks)",
           "about 7 l/s outdoor air per occupant"),
    "S1": ("Supermarket, staff (Liaocheng)", "not reported", "not reported"),
    "S2": ("Department store, sales staff (Tianjin)", "not reported", "described as poorly ventilated"),
    "S3": ("Shopping mall, 8 floors (Wenzhou)", "not reported", "not reported"),
    "S4": ("Grocery store, staff, PCR screen (Massachusetts)", "not reported", "not reported"),
    "C1": ("Primary-school classroom, 24 pupils in 5 rows (Marin County)", "school days; hours not reported",
           "portable HEPA filter; doors and windows open"),
    "C2": ("High school, 35--38 pupils per classroom of 39--49 m$^2$ (Jerusalem)",
           "8--9 school days of 6.3--6.7 h", "split air conditioning, recirculating"),
    "C3": ("20 primary schools, masks worn, seats 0.9 m apart (Utah)", "school days; not quantified",
           "not quantified"),
    "C4": ("Primary school, 14 classrooms on one recirculating system, measles (New York State)",
           "not in abstract", "largely recirculated air"),
    "C5": ("Fitness dance classes in a 60 m$^2$ room (Cheonan)", "50 min per class", "not reported"),
    "T1": ("Bus, 15 rows of seats (Zhejiang)", "50 min each way", "air conditioning on recirculation"),
    "T2": ("Coach, 49 seats, windows closed (Hunan)", "2.5 h (200 min in the second report)",
           "1.72 l/s per person, measured"),
    "T3": ("Minibus, 18 seats (Hunan)", "60 min", "3.22 l/s per person, measured"),
    "T4": ("High-speed trains, pooled contact-tracing cohort (China)", "mean 2.1 h", "not reported"),
    "T5": ("Long-haul flight, business cabin (London--Hanoi)", "10 h", "aircraft cabin; not quantified"),
    "T6": ("Airliner held on the ground, influenza (Alaska)", "3 h", "none (system inoperative)"),
    "T7": ("National cluster surveillance; commuter rail as a venue (Japan)", "--", "--"),
    "R1": ("Restaurant, 18 tables, 89 patrons (Guangzhou)", "82 min (index table)",
           "0.56--0.77 air changes per hour, measured"),
    "R2": (r"Restaurant, 9.2 m $\times$ 10.5 m (Jeonju)", "overlaps of 5 and 21 min",
           "no ventilation system; air-conditioner jet"),
    "H1": ("Choir rehearsal in a church hall (Skagit County)", "2.5 h", "0.3--1.0 h$^{-1}$, estimated"),
    "H2": ("Church services, singer in a choir loft (Sydney)", "1 h per service",
           "minimal: fans off, doors and windows largely closed"),
    "X1": ("Cruise ship under quarantine (Yokohama)", "weeks", "--"),
    "X2": ("Migrant-worker dormitories (Singapore)", "weeks to months", "--"),
}
GROUPS = [("O", "Offices and workplaces"), ("S", "Supermarket and retail"), ("C", "Classrooms and class-like rooms"),
          ("T", "Vehicle cabins (stand-ins for the metro carriage)"), ("R", "Restaurants"),
          ("H", "Choir and church"), ("X", "Outside the scope of the model (recorded for reference)")]
# grade-A records that do not meet the stated grade-A criteria (re-check, section 4)
GRADE_A_NOT_MET = {"C1": "exposure duration in hours and room size not reported",
                   "T4": "2,334 index patients pooled; cell denominators not printed",
                   "H1": "no seat or zone positions published"}
# counts reported differently by a second source (string must occur in the record's attack_rate_reported)
ALT_DENOMINATOR = {"T2": "7/46", "T3": "2/17"}

# order matters: the three page-extraction sources also mention the abstract held in Europe PMC
ACCESS_CLASS = [("automated text extraction", "extraction"), ("accepted manuscript", "manuscript"), ("abstract", "abstract"),
                ("full text", "full")]
ACCESS_MARK = {"full": "", "abstract": r"$^{\mathrm{a}}$", "extraction": r"$^{\mathrm{e}}$", "manuscript": r"$^{\mathrm{m}}$"}


def access_class(text: str) -> str:
    for needle, cls in ACCESS_CLASS:
        if needle in text:
            return cls
    raise ValueError(text)


# ------------------------------------------------------------------------------------------------ records
df = pd.read_csv(OB / "outbreaks.csv")
df = df.drop(columns=[c for c in df.columns if c.endswith("_cited")])   # bookkeeping column of the working data, not distributed
assert len(df) == 27 and df.shape[1] == 51
df["sid"] = df["id"].str.split("_").str[0]
numbers: dict = {"n_records": int(len(df)), "n_columns": int(df.shape[1] - 1)}

source_access: dict[str, str] = {}
rec_rows: list[str] = []
rec_json: list[dict] = []
for prefix, title in GROUPS:
    rec_rows.append(r"\addlinespace[2pt]\multicolumn{8}{@{}l}{\emph{" + title + r"}}\\")
    for _, r in df[df.sid.str[0] == prefix].iterrows():
        sid = r.sid
        setting, expo, vent = SHORT[sid]
        hay = " | ".join(str(r[c]) for c in ("setting", "location", "exposure_type", "exposure_duration_h",
                                             "exposure_duration_note", "room_dimensions", "floor_area_m2",
                                             "ventilation", "occupancy_n", "geometry_note", "attack_rate_reported"))
        for lab in (setting, expo, vent):
            _check_tokens(lab, hay, sid)
        # counts and attack rate
        dag = r"$^{\dagger}$" if r.includes_index is True or str(r.includes_index) == "True" else ""
        if pd.notna(r.n_exposed) and pd.notna(r.n_infected):
            k, n = int(r.n_infected), int(r.n_exposed)
            lo, hi = clopper_pearson(k, n)
            assert abs(k / n - r.attack_rate) < 5e-7 and abs(lo - r.ci95_low) < 5e-7 and abs(hi - r.ci95_high) < 5e-7, sid
            small = k / n < 0.01
            counts = f"{num(k)}/{num(n)}{dag}"
            if pd.notna(r.n_infected_alt):
                if sid == "H1":                                   # 32 laboratory-confirmed + 20 probable cases
                    assert "32 confirmed + 20 probable" in r.attack_rate_reported
                    counts += f" ({int(r.n_infected_alt)}/{num(n)} confirmed)"
                else:
                    counts += f" (or {int(r.n_infected_alt)}/{num(n)})"
            if sid in ALT_DENOMINATOR:
                assert ALT_DENOMINATOR[sid] in r.attack_rate_reported, sid
                if sid == "T2":                                   # 45 passengers have a seat in the stratum table
                    assert "3/19" in r.spatial_pattern and "4/26" in r.spatial_pattern
                    counts += f" (or {ALT_DENOMINATOR[sid]}; 45 with seat positions)"
                else:
                    counts += f" (or {ALT_DENOMINATOR[sid]})"
            if sid == "C5":                                       # numerator derived from the printed attack rate
                assert "n_infected = round(0.263*217)" in str(r.derived_fields)
                counts += r"$^{\S}$"
            if sid == "R1":
                counts += r"$^{\ddagger}$"
            ar = f"{pct(k / n, small)} ({pct(lo, small)}--{pct(hi, small)})"
            alt = None
            if pd.notna(r.n_infected_alt):
                alt = [int(r.n_infected_alt), n]
            elif sid in ALT_DENOMINATOR:
                alt = [int(x) for x in ALT_DENOMINATOR[sid].split("/")]
            rj = dict(id=sid, k=k, n=n, attack_rate=k / n, ci95=[lo, hi], alt=alt)
        elif pd.notna(r.n_infected):
            counts, ar = f"{int(r.n_infected)} cases; no denominator", "--"
            rj = dict(id=sid, k=int(r.n_infected), n=None)
        elif sid == "T6":                                         # attack rate reported, counts not transcribed
            assert "72 per cent of the passengers became ill" in r.attack_rate_reported and int(r.occupancy_n) == 54
            counts, ar = "72\\,\\% of passengers; 54 aboard (counts not transcribed)", "--"
            rj = dict(id=sid, k=None, n=None, reported_attack_rate_pct=72, aboard=54)
        else:
            counts, ar = "no denominator", "--"
            rj = dict(id=sid, k=None, n=None)
        # grade
        grade = r.quality + (r"$^{*}$" if sid in GRADE_A_NOT_MET else "")
        if sid in GRADE_A_NOT_MET:
            assert r.quality == "A"
        # sources with access marks
        cites = []
        for part in str(r.source_access).split(" | "):
            key, acc = part.split(": ", 1)
            cls = access_class(acc)
            assert source_access.setdefault(key, cls) == cls, key
            cites.append((key, cls))
        src = []
        for cls in ("full", "manuscript", "extraction", "abstract"):
            keys = [k_ for k_, c_ in cites if c_ == cls]
            if keys:
                src.append(r"\citep{" + ",".join(keys) + "}" + ACCESS_MARK[cls])
        rj.update(grade=r.quality, grade_criteria_met=sid not in GRADE_A_NOT_MET, sources=[k_ for k_, _ in cites],
                  scene=r.scene, includes_index=str(r.includes_index) == "True",
                  exposure_duration_h=None if pd.isna(r.exposure_duration_h) else float(r.exposure_duration_h))
        rec_json.append(rj)
        rec_rows.append(" & ".join([sid, setting, counts, ar, expo, vent, grade, ", ".join(src)]) + r" \\")

assert len(source_access) == 31
acc_counts = {c: sorted(k for k, v in source_access.items() if v == c) for c in ("full", "manuscript", "extraction", "abstract")}
numbers["sources"] = {"n": len(source_access), **{c: {"n": len(v), "keys": v} for c, v in acc_counts.items()}}
numbers["grades"] = {g: int((df.quality == g).sum()) for g in "ABCX"}
numbers["grade_A"] = {"ids": sorted(df[df.quality == "A"].sid), "not_meeting_criteria": GRADE_A_NOT_MET,
                      "n_meeting_criteria": int((df.quality == "A").sum()) - len(GRADE_A_NOT_MET)}
numbers["records_by_group"] = {t: int((df.sid.str[0] == p).sum()) for p, t in GROUPS}
numbers["records_without_denominator"] = sorted(s for s in df[df.n_exposed.isna()].sid if s != "T6")
numbers["records_with_reported_rate_but_no_counts"] = ["T6"]
numbers["records_with_counts"] = int((df.n_exposed.notna() & df.n_infected.notna()).sum())
numbers["records"] = rec_json
numbers["other_pathogens"] = sorted(df[~df.pathogen.str.contains("SARS-CoV-2")].sid)
# abstract-only records must not be graded above C
for _, r in df.iterrows():
    accs = {access_class(p.split(": ", 1)[1]) for p in str(r.source_access).split(" | ")}
    if accs == {"abstract"}:
        assert r.quality == "C", r.sid
numbers["records_resting_on_abstracts_only"] = sorted(
    r.sid for _, r in df.iterrows()
    if {access_class(p.split(": ", 1)[1]) for p in str(r.source_access).split(" | ")} == {"abstract"})

# ------------------------------------------------------------------------------------------------ strata
CONTRAST = {  # (record, row index within the record) -> short label; numbers checked against near/far/contrast
    ("O1", 0): "north wing against south wing (mapped seats)",
    ("C1", 0): "front two rows against back three rows",
    ("T1", 0): "rows 5--11 (within about 2 m of the index) against rows 1--4 and 12--15",
    ("T1", 1): "rows 6--10 against rows 1--5 and 11--15 (alternative split)",
    ("T2", 0): "rear rows 8--13 (index zone) against front rows 1--7",
    ("T5", 0): "within 2 seats of the index against farther seats, business cabin",
    ("R1", 0): "neighbouring tables against remote tables",
    ("R1", 1): "same air-conditioning zone against other zones (alternative split)",
    ("C2", 0): "grades 7--9 against grades 10--12 (separate wings)",
    ("S1", 0): "staff against customers",
}
rr = pd.read_csv(OB / "data/near_far_risk_ratios.csv")
strata_rows, strata_json, seen = [], [], {}
for _, r in rr.iterrows():
    sid = r.outbreak_id.split("_")[0]
    j = seen.get(sid, 0)
    seen[sid] = j + 1
    lab = CONTRAST[(sid, j)]
    _check_tokens(lab, f"{r.contrast} | {r.near} | {r.far}", f"strata {sid}")
    k1, n1, k0, n0 = int(r.near_k), int(r.near_n), int(r.far_k), int(r.far_n)
    v, lo, hi, cc, p, lg, se = risk_ratio(k1, n1, k0, n0)
    assert abs(v / r.risk_ratio - 1) < 1e-9 and abs(lo / r.rr_ci_low - 1) < 1e-4 and abs(hi / r.rr_ci_high - 1) < 1e-4, sid
    assert cc == bool(r.continuity_corrected) and abs(p / r.fisher_p - 1) < 1e-6, sid
    star = r"$^{\mathrm{c}}$" if cc else ""
    mark = r"$^{\ddagger}$" if sid == "R1" else ""
    if sid == "C2":      # enrolled pupils by grade (581 and 583); 1,161 were tested, 580 of them in grades 10-12
        assert (n1, n0) == (581, 583)
        mark = r"$^{\mathrm{e}}$"
    strata_rows.append(" & ".join([
        sid, lab, f"{num(k1)}/{num(n1)}{mark if sid != 'C2' else ''} ({100 * k1 / n1:.1f})",
        f"{num(k0)}/{num(n0)}{mark if sid == 'C2' else ''} ({100 * k0 / n0:.1f})",
        f"{rrfmt(v)}{star} ({rrfmt(lo)}--{rrfmt(hi)})", pfmt(p)]) + r" \\")
    strata_json.append(dict(id=sid, contrast=lab, near=[k1, n1], far=[k0, n0], rr=v, ci95=[lo, hi],
                            continuity_corrected=cc, fisher_p=p, log_rr=lg, se=se))
numbers["strata"] = strata_json

# Cochran Q over the four single-index events with one near/far split each (classroom, bus split 1, coach, flight)
sel = [s for s in strata_json if (s["id"], s["contrast"]) in
       {("C1", CONTRAST[("C1", 0)]), ("T1", CONTRAST[("T1", 0)]), ("T2", CONTRAST[("T2", 0)]), ("T5", CONTRAST[("T5", 0)])}]
w = np.array([1 / s["se"] ** 2 for s in sel])
y = np.array([s["log_rr"] for s in sel])
mu = float((w * y).sum() / w.sum())
Q = float((w * (y - mu) ** 2).sum())
numbers["heterogeneity"] = dict(events=[s["id"] for s in sel], cochran_Q=Q, dof=3, p=float(stats.chi2.sf(Q, 3)),
                                pooled_rr=math.exp(mu), pooled_ci95=[math.exp(mu - 1.96 / math.sqrt(w.sum())),
                                                                     math.exp(mu + 1.96 / math.sqrt(w.sum()))])
ext = json.loads((VER / "extra_checks.json").read_text())
assert abs(Q - ext["heterogeneity"]["Q"]) < 1e-6
stc = json.loads((VER / "stats_check.json").read_text())
t2side = next(x for x in stc["risk_ratios"] if x["id"] == "T2side")
v, lo, hi, *_ = risk_ratio(6, 24, 1, 20)
assert abs(v - t2side["rr"]) < 1e-9
numbers["coach_side_contrast"] = dict(near=[6, 24], far=[1, 20], rr=v, ci95=[lo, hi])

# ------------------------------------------------------------------------------------------------ negative controls
CONTROL = {  # name in analysis_results.json -> (label, record, source key)
    "S1 customers": ("Supermarket customers", "S1", "tian2021"),
    "T1 unexposed bus": ("Other bus travelling to the same event", "T1", "shen2020"),
    "R1 remote tables": ("Restaurant, remote tables", "R1", "li2021"),
    "R1 other zones": ("Restaurant, other air-conditioning zones", "R1", "li2021"),
    "T5 premium economy": ("Flight, premium-economy cabin", "T5", "khanh2020"),
    "C3 Utah classroom contacts": ("Tested school contacts, masked classrooms", "C3", "hershow2021"),
    "O1 floors 7-9": ("Call-centre building, floors 7--9", "O1", "park2020"),
}
ana = json.loads((OB / "data/analysis_results.json").read_text())
ctrl_rows, ctrl_json = [], []
for c in ana["negative_controls"]:
    lab, sid, key = CONTROL[c["name"]]
    k, n = int(c["k"]), int(c["n"])
    ub = upper_one_sided(k, n)
    assert abs(ub - c["upper95_one_sided"]) < 1e-9, c["name"]
    assert key in source_access
    ubs = f"{100 * ub:.3f}" if ub < 0.001 else (f"{100 * ub:.2f}" if ub < 0.02 else f"{100 * ub:.1f}")
    ctrl_rows.append(" & ".join([lab, sid, f"{num(k)}/{num(n)}", ubs, r"\citep{" + key + "}" + ACCESS_MARK[source_access[key]]]) + r" \\")
    ctrl_json.append(dict(name=lab, record=sid, k=k, n=n, upper95_one_sided=ub, source=key))
numbers["negative_controls"] = ctrl_json
# the two restaurant rows are alternative splits of the same patrons; the second bus is in the T1 record's notes
t1 = df[df.sid == "T1"].iloc[0]
numbers["negative_controls_note"] = dict(t1_other_groups=str(t1.other_groups))

# ------------------------------------------------------------------------------------------------ call centre
seats = pd.read_csv(OB / "data/park2020_callcentre_seats_digitized.csv")
g = seats.groupby("region").case.agg(["sum", "count"])
kn, nn = int(g.loc["north_wing", "sum"]), int(g.loc["north_wing", "count"])
ks, ns = int(g.loc["south_wing", "sum"]), int(g.loc["south_wing", "count"])
mapped = int(seats.case.sum())
o1 = df[df.sid == "O1"].iloc[0]
reported = int(o1.n_infected)
unmapped = reported - mapped
cc_ver = json.loads((VER / "callcentre_check.json").read_text())
assert cc_ver["manual_counts"]["north_desks"] == nn and cc_ver["manual_counts"]["north_cases"] == kn
assert cc_ver["manual_counts"]["south_desks"] == ns and cc_ver["manual_counts"]["south_cases"] == ks
assert cc_ver["manual_counts"]["total_case_seats"] == mapped
v, lo, hi, *_ = risk_ratio(kn, nn, ks, ns)
vs, los, his, *_ = risk_ratio(kn, nn, ks + unmapped, ns)          # all unmapped cases in the south wing
vn, lon, hin, *_ = risk_ratio(kn + unmapped, nn, ks, ns)          # all in the north wing
assert abs(vs - cc_ver["wing_contrast"]["sens_all10_south"]["rr"]) < 1e-9
pitch_px = float((seats.x_px / seats.x_pitch_units).median())
numbers["callcentre"] = dict(
    north=[kn, nn], south=[ks, ns], east_offices_case_seats=int(g.loc["east_offices", "sum"]),
    desks_digitised_two_wings=nn + ns, case_seats_mapped=mapped, cases_reported=reported, cases_unmapped=unmapped,
    employees_floor=int(o1.n_exposed), desk_pitch_px=pitch_px,
    rr=v, rr_ci95=[lo, hi], rr_all_unmapped_south=dict(rr=vs, ci95=[los, his]),
    rr_all_unmapped_north=dict(rr=vn, ci95=[lon, hin]),
    check_match=dict(bijection=cc_ver["comparison_with_first_analysis"]["matching_is_bijection"],
                        desks_with_different_case_flag=cc_ver["comparison_with_first_analysis"]["desks_with_different_case_flag"],
                        max_position_difference_px=cc_ver["comparison_with_first_analysis"]["max_match_distance_px"]))

# ------------------------------------------------------------------------------------------------ train cohort
tm = pd.read_csv(VER / "train_matrix_reconstruction.csv")
worst = 0.0
for _, r in tm.iterrows():
    lo_, hi_ = wilson(int(r.k), int(r.n))
    err = max(abs(100 * r.k / r.n - r.pct), abs(100 * lo_ - r.lo), abs(100 * hi_ - r.hi))
    worst = max(worst, err)
assert worst < 0.0125, worst            # printed to two decimals: half a unit in the last place, plus slack
t4 = df[df.sid == "T4"].iloc[0]
numbers["train"] = dict(cells=int(len(tm)), recovered_cases=int(tm.k.sum()), recovered_contacts=int(tm.n.sum()),
                        reported_cases=int(t4.n_infected), reported_contacts=int(t4.n_exposed),
                        relative_difference_contacts=float(tm.n.sum() / t4.n_exposed - 1),
                        max_abs_error_pp_vs_printed=worst, cell_n_min=int(tm.n.min()), cell_n_max=int(tm.n.max()),
                        adjacent_seat=dict(k=int(tm.k.iloc[0]), n=int(tm.n.iloc[0]), pct=float(tm.pct.iloc[0])),
                        cells_with_at_most_3_cases=int((tm.k <= 3).sum()),
                        # total absolute error (percentage points) of the six column means recomputed from the
                        # recovered cells against the printed column means, for the two placements of the
                        # same-row cells (re-check, stats_check.json -> train_matrix)
                        column_mean_error_pp_as_printed=stc["train_matrix"]["abs_err_printed_layout"],
                        column_mean_error_pp_realigned=stc["train_matrix"]["abs_err_shifted"])

# ------------------------------------------------------------------------------------------------ hazards
hz = ana["hazard_per_hour_range"]
numbers["hazard_per_hour"] = dict(n_events=hz["n"], min=hz["min"], max=hz["max"], ratio=hz["ratio"],
                                  ratio_confirmed_choir_only=stc["hazard_per_hour"]["confirmed_only_choir"]["ratio"])

(OUT / "appB_dataset_numbers.json").write_text(json.dumps(numbers, indent=1))

# ------------------------------------------------------------------------------------------------ write / inject
blocks = {"records": rec_rows, "strata": strata_rows, "controls": ctrl_rows}
tex = ["% generated by code/appB_dataset_table.py -- do not edit"]
for name, rows in blocks.items():
    tex += [f"% ---- body of Table appB_dataset:tab:{name}"] + rows
(OUT / "appB_dataset_tables.tex").write_text("\n".join(tex) + "\n")

if SEC.exists():
    s = SEC.read_text()
    done = []
    for name, rows in blocks.items():
        a, b = f"% >>> GENERATED {name}", f"% <<< GENERATED {name}"
        if a in s and b in s:
            i, j = s.index(a) + len(a), s.index(b)
            s = s[:i] + " (code/appB_dataset_table.py; do not edit by hand)\n" + "\n".join(rows) + "\n" + s[j:]
            done.append(name)
    SEC.write_text(s)
    print("injected into sections/supp_datasets.tex:", done)

print(json.dumps({k: numbers[k] for k in ("sources", "grades", "grade_A", "records_by_group", "heterogeneity",
                                           "callcentre", "train", "hazard_per_hour")}, indent=1))

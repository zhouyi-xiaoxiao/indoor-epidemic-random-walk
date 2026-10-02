#!/usr/bin/env python
"""v2_numbers.py -- numbers and table bodies of the second set of tests (Sections 7.1, 7.5-7.10; Supplementary Sections S6, S9, S10).

The script reads the stored results of the four analyses of the second set (further outbreaks, contact
records, tracer measurements, zone level) and of their re-checks (separately written code), and the corrected
hold-out scores of the first test (s8_first_test_exact.json); the only quantities it computes itself are pooled
near/far ratios without one record and the probability of the pooled pattern of the tracer analysis.
It collects every number quoted in the article in

    data/v2_numbers.json          (value, with the file and key it comes from)

and writes the bodies of the generated tables to data/v2_tab_*.tex (read by \\input).

The freeze evidence of Supplementary Section S10 (file times and the hash check over all files, including third-party source
files that are not redistributed) can only be read in the author's working tree; it is recorded there with
--record-freeze in data/v2_freeze_evidence.json.  An ordinary run reads that record and re-checks the
hashes of the files that are present.

Run:  python code/v2_numbers.py [--record-freeze]
"""
from __future__ import annotations

import hashlib
import json
import os
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
OUT = HERE.parent / "data"

A1 = "code/bsc_validation2/a1_more_outbreaks/"
A2 = "code/bsc_validation2/a2_contact_mobility/"
A3 = "code/bsc_validation2/a3_tracer_physics/"
A4 = "code/bsc_validation2/a4_zone_level/"
# results, tracer and school tables, and re-checks (written out in full: the release build maps each
# of these directories to its place in the repository)
A1O = "data/bsc_validation2/a1_more_outbreaks/"
A1V = "code/checks/outbreaks2/"
A2O = "data/bsc_validation2/a2_contact_mobility/"
A2V = "code/checks/contacts/"
A2SUM = "code/checks/contacts/check_summary.json"
A3R = "data/bsc_validation2/a3_tracer_physics/"
A3D = "data/tracer/"
A4R = "data/bsc_validation2/a4_zone_level/"
A4D = "data/schools/"
A4V = "code/checks/zones/"
R1D = "data/bsc_validation/"


def load(rel):
    return json.loads((ROOT / rel).read_text())


def frozen_sha256(path) -> str:
    """SHA-256 of a frozen protocol file.  In the public repository the protocol records are redacted copies
    (code/REDACTIONS.md): such a copy is first checked against the hash listed for it in code/REDACTIONS.json,
    and the hash of the unredacted original listed there is returned, so that the comparison with the freeze
    record is a comparison of the originals (as recorded) and the copy is checked against the list."""
    path = Path(path)
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    listing = ROOT / "code" / "REDACTIONS.json"
    if listing.exists():
        red = json.loads(listing.read_text())
        try:
            key = Path(os.path.normpath(path.absolute())).relative_to(os.path.normpath(ROOT.absolute())).as_posix()
        except ValueError:
            return h
        if key in red:
            assert h == red[key]["sha256_copy"], f"redacted copy differs from code/REDACTIONS.json: {key}"
            return red[key]["sha256_original"]
    return h


def sha(rel):
    return frozen_sha256(ROOT / rel)


N = {}          # the collected numbers


def put(key, value, src):
    N[key] = {"value": value, "src": src}
    return value


# ------------------------------------------------------------------------------------------- formatting
def sgn(x, nd=2):
    return f"${'+' if x >= 0 else '-'}{abs(x):.{nd}f}$"


def pval(p):
    """p-value for a LaTeX table."""
    if p is None:
        return "--"
    if p >= 0.995:
        return "1.00"
    if p >= 0.0995:
        return f"{p:.2f}"
    if p >= 0.0095:
        return f"{p:.3f}"
    if p >= 0.00095:
        return f"{p:.4f}"
    e = int(math.floor(math.log10(p)))
    m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = 1.0, e + 1
    return f"${m:.1f}\\times10^{{{e}}}$"


def vec(v, nd=1):
    return ", ".join(f"{x:.{nd}f}" if not float(x).is_integer() or nd else f"{int(x)}" for x in v)


def ivec(v):
    return ", ".join(str(int(x)) for x in v)


# =================================================================================== 1. further outbreaks
cor = load(A1O + "10_corrected.json")
reg = load(A1O + "04_holdout.json")
sec = load(A1O + "05_secondary.json")
des = load(A1O + "07_descriptive.json")
var = load(A1O + "06_variants.json")
num = load(A1O + "09_numerics_check.json")
chk = load(A1V + "t2_results.json")
chk_bins = load(A1V + "t5_bins.json")
chk_m1 = load(A1V + "t7_m1.json")
chk_tr = load(A1V + "t4_transfer.json")
P = cor["primary_corrected"]

LABEL = {
    "F1": "Hong Kong--Beijing, CA112", "F2": "Chicago--Honolulu", "F3": "Los Angeles--Auckland",
    "F4": "Cancun--Birmingham", "F5": "flight to Naha", "F6": "Sydney--Perth", "F7": "Dubai--Auckland, EK448",
    "F8": "Tel Aviv--Frankfurt", "W1": "ward 8A, students", "W2": "ward 8A, inpatients",
    "P1": "meat-processing line"}
CITE = {"F1": "olsen2003transmission", "F2": "kenyon1996transmission", "F3": "baker2010transmission",
        "F4": "young2014international", "F5": "toyokawa2022transmission", "F6": "speake2020flight",
        "F7": "swadi2021genomic", "F8": "hoehl2020assessment", "W1": "wong2004cluster", "W2": "yu2005temporal",
        "P1": "guenther2020"}
ORDER = ["F1", "F2", "F3", "F4", "W1", "W2", "F5", "F6", "F7", "F8", "P1"]
HOLD = ["F5", "F6", "F7", "F8", "P1"]

# --- descriptive near/far risk ratios
src = A1O + "07_descriptive.json"
for e in ORDER:
    r = des["events"][e]
    put(f"out2.rr.{e}", dict(near=r["near"], far=r["far"], RR=r["RR"], lo=r["lo"], hi=r["hi"],
                             pathogen=r["pathogen"], role=r["role"], cc=r["cc"],
                             pred_M2=r["pred_RR"]["M2"], pred_CRR=r["pred_RR"]["CRR"]), src + f" (events.{e})")
for k, v in des["pooled"].items():
    put(f"out2.rr_pooled.{k}", dict(n=v["n"], RR_random=v.get("RR_random"), lo=v.get("lo_r"), hi=v.get("hi_r"),
                                    RR_fixed=v["RR_fixed"], Q=v.get("Q"), df=v.get("df"), p_Q=v.get("p_Q"),
                                    I2=v.get("I2")), src + f" (pooled.{k})")
put("out2.rr_pooled.check", chk["RR"], A1V + "t2_results.json (RR)")


def _dl(evs):
    """DerSimonian-Laird pool of the stored log risk ratios (same estimator as 07_descriptive.py)."""
    y = np.array([des["events"][e]["logRR"] for e in evs])
    w = 1 / np.array([des["events"][e]["se"] for e in evs]) ** 2
    m = (w * y).sum() / w.sum()
    Q = float((w * (y - m) ** 2).sum())
    df = len(evs) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (1 / w + tau2)
    mr = (wr * y).sum() / wr.sum()
    sr = 1 / np.sqrt(wr.sum())
    return dict(n=len(evs), RR_random=float(np.exp(mr)), lo=float(np.exp(mr - 1.96 * sr)),
                hi=float(np.exp(mr + 1.96 * sr)), Q=Q, df=df, p_Q=float(stats.chi2.sf(Q, df)))


_all = list(des["events"])
_chk = _dl(_all)
assert abs(_chk["RR_random"] - des["pooled"]["all 11 events"]["RR_random"]) < 1e-9, "pooling not reproduced"
put("out2.rr_pooled.without P1", _dl([e for e in _all if e != "P1"]),
    A1O + "07_descriptive.json (events.*.logRR, se), pooled here without the processing line (P1) with the "
    "DerSimonian-Laird estimator of 07_descriptive.py, which reproduces its pooled value for all 11 records")
put("out2.rr_pooled.without W2", _dl([e for e in _all if e != "W2"]),
    A1O + "07_descriptive.json (events.*.logRR, se), pooled here without the second ward record (W2)")

# --- calibration
src = A1O + "10_corrected.json (primary_corrected.post)"
for k, v in P["post"].items():
    put(f"out2.cal.{k}", v, src)
put("out2.cal.registered_room_M2", reg["post"]["room|M2"] if "post" in reg else sec["CRR_2.4_primary_split"]["post"]["room|M2"],
    A1O + "05_secondary.json (CRR_2.4_primary_split.post.room|M2: the posterior as first computed)")
put("out2.cal.maxloglik_check", chk["calib_maxloglik"], A1V + "t2_results.json (calib_maxloglik)")

# --- hold-out events
src = A1O + "10_corrected.json (primary_corrected.events)"
for e in HOLD:
    r = P["events"][e]
    put(f"out2.hold.{e}", dict(obs=r["obs"], sizes=r["sizes"],
                               **{m: dict(logp=r[m]["logp"], p=r[m]["p_adeq"], expected=r[m]["expected"])
                                  for m in ("M0", "CRR", "M2", "M1")},
                               d0=r["M2"]["logp"] - r["M0"]["logp"], dC=r["M2"]["logp"] - r["CRR"]["logp"]),
        src + f".{e}")
for grp in ("all", "air", "room"):
    put(f"out2.pooled.M2.{grp}", P["pooled"]["M2"][grp], A1O + f"10_corrected.json (primary_corrected.pooled.M2.{grp})")
    put(f"out2.pooled.M1.{grp}", P["pooled"]["M1"][grp], A1O + f"10_corrected.json (primary_corrected.pooled.M1.{grp})")
    put(f"out2.pooled_check.M2.{grp}", chk["primary"][grp] if grp in chk["primary"] else chk["primary"]["all"],
        A1V + f"t2_results.json (primary.{grp})")
put("out2.pooled.M2.decision", P["pooled"]["M2"]["decision"], A1O + "10_corrected.json (primary_corrected.pooled.M2.decision)")
put("out2.pooled.M1.decision", P["pooled"]["M1"]["decision"], A1O + "10_corrected.json (primary_corrected.pooled.M1.decision)")
put("out2.pooled.M3.all", P["pooled"]["M3"]["all"], A1O + "10_corrected.json (primary_corrected.pooled.M3.all)")
put("out2.pooled.M3.decision", P["pooled"]["M3"]["decision"], A1O + "10_corrected.json (primary_corrected.pooled.M3.decision)")
put("out2.pooled_registered.M2.all", reg["pooled"]["M2"]["all"], A1O + "04_holdout.json (pooled.M2.all; kernel as first frozen)")
put("out2.pooled_registered.P1", {m: reg["events"]["P1"][m] for m in ("M0", "CRR", "M2")},
    A1O + "04_holdout.json (events.P1)")
plant_share = P["pooled"]["M2"]["room"]["deltaC"] / P["pooled"]["M2"]["all"]["deltaC"]
put("out2.pooled.plant_share_of_deltaC", plant_share, "ratio of the two entries above (room / all)")
put("out2.check.M1", {k: chk_m1[k] for k in ("all", "air", "room", "decision")}, A1V + "t7_m1.json")
put("out2.check.M1_events", {e: chk_m1["events"][e]["M1"] for e in HOLD}, A1V + "t7_m1.json (events)")

# --- power
put("out2.power.registered", {t: {k: v for k, v in reg_p.items() if k.startswith("P_")}
                              for t, reg_p in load(A1O + "03_power.json")["rules"]["M2"]["truth"].items()},
    A1O + "03_power.json (rules.M2.truth)")
put("out2.power.registered_air", load(A1O + "03_power.json")["rules"]["M2"]["air_only"],
    A1O + "03_power.json (rules.M2.air_only)")
put("out2.power.corrected", cor["power_corrected"]["M2"], A1O + "10_corrected.json (power_corrected.M2)")
put("out2.power.check", chk["power"], A1V + "t2_results.json (power)")

# --- sensitivity to the definition of the competitor and of the strata (re-check)
for k, name in (("CRR_near<=1row / P1<=4m", "near1_4m"), ("CRR_near<=5rows / P1<=12m", "near5_12m"),
                ("CRR_fixed2.4", "rho2.4"), ("vs_trueCRR", "registered_form"), ("S_hh", "households"),
                ("generic_EXP_target", "exp_vs_crr"), ("M2_vs_EXP", "m2_vs_exp"), ("swap", "swap"),
                ("S_noW2", "noW2")):
    x = chk[k]
    put(f"out2.check.{name}", {g: x[g] for g in ("all", "air", "room", "decision") if g in x}
        | ({"post": x["post"]} if "post" in x else {}), A1V + f"t2_results.json ({k})")
for k, x in chk_bins.items():
    put(f"out2.check.bins.{k}", {g: x[g] for g in ("all", "air", "decision")}, A1V + f"t5_bins.json ({k})")
put("out2.check.heterogeneity", chk["heterogeneity"], A1V + "t2_results.json (heterogeneity)")
put("out2.check.train_kernel_events", {e: chk_tr[e]["M2"] for e in chk_tr}, A1V + "t4_transfer.json")
put("out2.check.loo_air", {k: chk["LOO_air"][k] for k in ("delta0", "deltaC")}, A1V + "t2_results.json (LOO_air)")

# --- secondary analyses of the protocol
put("out2.sec.rho2.4", sec["CRR_2.4_primary_split"]["pooled"]["M2"], A1O + "05_secondary.json (CRR_2.4_primary_split.pooled.M2)")
t4 = cor["S-T4_corrected"]
put("out2.sec.train_kernel.flights", sec["S-T4_newflights_only"]["pooled"]["M2"]["all"],
    A1O + "05_secondary.json (S-T4_newflights_only.pooled.M2.all)")
put("out2.sec.train_kernel.padeq", {e: sec["S-T4_newflights_only"]["events"][e]["M2"]["p_adeq"]
                                    for e in ("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8")},
    A1O + "05_secondary.json (S-T4_newflights_only.events.*.M2.p_adeq)")
put("out2.sec.train_kernel.posterior", sec["S-T4"]["round1_posterior"]["air|M2"], A1O + "05_secondary.json (S-T4.round1_posterior)")
put("out2.sec.train_kernel.rooms", {e: t4["events"][e]["M2"]["p_adeq"] for e in ("W1", "W2", "P1")},
    A1O + "10_corrected.json (S-T4_corrected.events.*.M2.p_adeq)")
sw = cor["S-swap_corrected"]["pooled"]["M2"]
put("out2.sec.swap", {g: sw[g] for g in ("all", "air", "room", "decision")}, A1O + "10_corrected.json (S-swap_corrected.pooled.M2)")
put("out2.sec.swap_W2", cor["S-swap_corrected"]["events"]["W2"]["M2"]["p_adeq"], A1O + "10_corrected.json (S-swap_corrected.events.W2)")
put("out2.sec.loo_air", cor["S-LOO_pooled_air_corrected"]["M2"], A1O + "10_corrected.json (S-LOO_pooled_air_corrected.M2)")
put("out2.sec.loo_room", cor["S-LOO_pooled_room_corrected"]["M2"], A1O + "10_corrected.json (S-LOO_pooled_room_corrected.M2)")
put("out2.sec.seat", sec["S-seat"]["pooled"], A1O + "05_secondary.json (S-seat.pooled)")
vv = {}
for k, x in var.items():
    if k.startswith("primary"):
        continue
    try:
        m = x["pooled"]["M2"]
        vv[k] = dict(d0=m["all"]["delta0"], dC=m["all"]["deltaC"], dC_air=m["air"]["deltaC"],
                     pB2_air=m["air"].get("p_B2"), pmirror_air=m["air"].get("p_mirror"),
                     outcome=m["decision"]["outcome"], fails=m["decision"]["B3_failures"])
    except (KeyError, TypeError):
        pass
put("out2.variants", vv, A1O + "06_variants.json (<variant>.pooled.M2)")
put("out2.variants.summary", dict(n=len(vv), d0_min=min(v["d0"] for v in vv.values()), d0_max=max(v["d0"] for v in vv.values()),
                                  dC_min=min(v["dC"] for v in vv.values()), dC_max=max(v["dC"] for v in vv.values()),
                                  dC_air_min=min(v["dC_air"] for v in vv.values()), dC_air_max=max(v["dC_air"] for v in vv.values()),
                                  n_V5=sum(v["outcome"] == "V5" for v in vv.values()),
                                  n_V3=sum(v["outcome"] == "V3" for v in vv.values()),
                                  air_sig_M2_better=sorted(k for k, v in vv.items() if v["pB2_air"] < 0.05),
                                  air_sig_CRR_better=sorted(k for k, v in vv.items() if v["pmirror_air"] < 0.05)),
    A1O + "06_variants.json (summary over the variants)")
# the variants were computed with the room kernel as first frozen (before the numerical correction): their reference
ref = var["primary (reference)"]["pooled"]["M2"]
put("out2.variants.reference", dict(d0=ref["all"]["delta0"], dC=ref["all"]["deltaC"], dC_air=ref["air"]["deltaC"]),
    A1O + "06_variants.json (primary (reference).pooled.M2: the primary analysis with the kernel as first frozen)")
put("out2.sec.train_kernel.pooled", {g: t4["pooled"]["M2"][g] for g in ("all", "air", "room", "decision")},
    A1O + "10_corrected.json (S-T4_corrected.pooled.M2)")
put("out2.sec.loo_both", cor["S-LOO_pooled_both_corrected"]["M2"] if "M2" in cor["S-LOO_pooled_both_corrected"]
    else cor["S-LOO_pooled_both_corrected"], A1O + "10_corrected.json (S-LOO_pooled_both_corrected)")

# --- mixing coefficient preferred by each event
pd_ = cor["preferred_D_corrected"]
for cls in ("air", "room"):
    x = pd_[cls]
    put(f"out2.D.{cls}", dict(D_common=x["D_common"], D_common_ci=x["D_common_ci"], LR=x["LR"], df=x["df"], p=x["p"],
                              per_event={e: {k: v[k] for k in ("D_hat", "D_lo", "D_hi", "gain_over_wellmixed")}
                                         for e, v in x["per_event"].items()},
                              by_pathogen={k: {kk: v[kk] for kk in ("D_hat", "D_lo", "D_hi", "gain_over_wellmixed")}
                                           for k, v in x["by_pathogen"].items()}),
        A1O + f"10_corrected.json (preferred_D_corrected.{cls})")
put("out2.numerics", {k: num[k] for k in list(num)[:6]} if isinstance(num, dict) else num, A1O + "09_numerics_check.json")

# --- table: hold-out events (Section 7.5)
rows = []
for e in HOLD:
    r = N[f"out2.hold.{e}"]["value"]
    near, far = des["events"][e]["near"], des["events"][e]["far"]
    nf = f"{near[0]}/{near[1]} vs {far[0]}/{far[1]}"
    rows.append(" & ".join([
        LABEL[e][0].upper() + LABEL[e][1:] + f" \\citep{{{CITE[e]}}}", nf, ivec(r["obs"]) + f" ({ivec(r['sizes'])})",
        vec(r["CRR"]["expected"]), vec(r["M2"]["expected"]),
        f"{pval(r['M0']['p'])} / {pval(r['CRR']['p'])} / {pval(r['M2']['p'])}",
        sgn(r["d0"]), sgn(r["dC"])]) + " \\\\")
pm = P["pooled"]["M2"]
rows.append("\\midrule")
rows.append(f"Four flights & & & & & & {sgn(pm['air']['delta0'])} & {sgn(pm['air']['deltaC'])} \\\\")
rows.append(f"All five events & & & & & & {sgn(pm['all']['delta0'])} & {sgn(pm['all']['deltaC'])} \\\\")
(OUT / "v2_tab_holdout2.tex").write_text(
    "% generated by code/v2_numbers.py from " + A1O + "10_corrected.json (primary_corrected) and out/07_descriptive.json\n"
    + "\n".join(rows) + "\n")

# --- table: the eleven events (Supplementary Section S6)
META = {  # setting and duration; risk set: code/bsc_validation2/a1_more_outbreaks/PREREGISTRATION.md section 2.1
    "F1": ("SARS-CoV-1", "B737-300, 3\\,h", "full seat map (redrawn in \\citep{hertzberg2016risk})", "A"),
    "F2": ("\\emph{M.\\ tuberculosis}", "wide-body, 8.75\\,h", "two strata of the abstract; no seat map", "A"),
    "F3": ("influenza A(H1N1)\\allowbreak pdm09", "B747-400, rear section, 13\\,h", "seat map; 9 sources", "A"),
    "F4": ("influenza A(H1N1)\\allowbreak pdm09", "B767, 9.5\\,h", "seat map; 6 sources", "A"),
    "W1": ("SARS-CoV-1", "hospital ward, 40\\,min", "three strata by bed; floor plan", "A"),
    "W2": ("SARS-CoV-1", "same ward, days", "three strata by bay; floor plan", "B"),
    "F5": ("SARS-CoV-2", "B737-800, 2\\,h", "seat map; index in 23C", "A"),
    "F6": ("SARS-CoV-2", "A330-200, mid cabin, 5\\,h", "seat map; 4 sources", "A"),
    "F7": ("SARS-CoV-2", "B777-300ER, 18\\,h", "seats within rows 23--30; 2 sources", "A"),
    "F8": ("SARS-CoV-2", "B737-900, 4.7\\,h", "seat map; 7 sources", "A"),
    "P1": ("SARS-CoV-2", "$32\\times8.5$\\,m hall, 3 shifts", "counts by 1-m distance band", "A"),
}
rows = []
for e in ORDER:
    r = des["events"][e]
    path, setting, pos, grade = META[e]
    cc = "$^{a}$" if r["cc"] else ""
    rows.append(" & ".join([
        e, LABEL[e][0].upper() + LABEL[e][1:] + f" \\citep{{{CITE[e]}}}", path, setting, pos,
        f"{r['near'][0]}/{r['near'][1]} vs {r['far'][0]}/{r['far'][1]}",
        f"{r['RR']:.2f}{cc} ({r['lo']:.2f}--{r['hi']:.1f})" if r["hi"] >= 10 else f"{r['RR']:.2f}{cc} ({r['lo']:.2f}--{r['hi']:.2f})",
        grade, "calibration" if r["role"].startswith("calibration") else "hold-out"]) + " \\\\")
    if e == "W2":
        rows.append("\\addlinespace")
(OUT / "v2_tab_events2.tex").write_text(
    "% generated by code/v2_numbers.py from " + A1O + "07_descriptive.json (events); setting and position data: "
    + A1 + "PREREGISTRATION.md section 2.1\n" + "\n".join(rows) + "\n")

# --- table: mixing coefficient by event (Section 7.6)
rows = []
for cls in ("air", "room"):
    for e, v in pd_[cls]["per_event"].items():
        lo, hi = v["D_lo"], v["D_hi"]
        dh = v["D_hat"]
        if dh <= 0.011:
            d = "lower edge"
        elif dh >= 9999:
            d = "upper edge"
        else:
            d = f"{dh:,.0f}".replace(",", "{,}")
        if lo <= 0.011:
            ci = f"$\\le{hi:,.0f}$".replace(",", "{,}")
        elif hi >= 9999:
            ci = f"$\\ge{lo:,.0f}$".replace(",", "{,}")
        else:
            ci = f"{lo:,.0f}--{hi:,.0f}".replace(",", "{,}")
        rows.append(f"{LABEL[e][0].upper() + LABEL[e][1:]} & {des['events'][e]['pathogen']} & {d} & {ci} & "
                    f"{v['gain_over_wellmixed']:.2f} \\\\")
    if cls == "air":
        rows.append("\\addlinespace")
(OUT / "v2_tab_Dprofiles.tex").write_text(
    "% generated by code/v2_numbers.py from " + A1O + "10_corrected.json (preferred_D_corrected.*.per_event)\n"
    + "\n".join(rows).replace("M. tuberculosis", "\\emph{M.\\ tuberculosis}") + "\n")

# ====================================================================================== 2. contact records
DS = ["InVS13", "LyonSchool", "LH10", "InVS15", "Thiers13", "SFHH"]
DSNAME = {"InVS13": "office, 2013", "LyonSchool": "primary school", "LH10": "hospital ward",
          "InVS15": "office, 2015", "Thiers13": "high school", "SFHH": "conference"}
STAGE = {"InVS13": "s1", "LyonSchool": "s1", "LH10": "s1", "InVS15": "s3a", "Thiers13": "s3a", "SFHH": "s3a"}
FAM = ["S0_level", "S1_duration", "S2_intercontact", "S3_degree", "S4_heterogeneity", "S5_persistence", "S6_groups"]
MODS = {"WM": {}, "BLOCK": {}, "RW0": {}, "RWhet": {}, "RWclus": {}, "RWstickZ": {}, "RWstickO": {}}
EXTRA = {"RWclus": {"dev": "s2b", "conf": "s3b"}, "RWstickZ": {"dev": "s2d", "conf": "s3c"},
         "RWstickO": {"dev": "s2d", "conf": "s3d"}}
con = {}
for ds in DS:
    base = load(A2O + f"{STAGE[ds]}_{ds}.json")
    con[ds] = dict(N=base["N"], n_days=base["n_days"], n_train=len(base["train"]), n_test=len(base["test"]),
                   n_groups=base["n_groups"], M=base["cal"]["M"], p=base["cal"]["p"], mean_dur=base["cal"]["mean_dur"],
                   obs={k: base["obs"][k] for k in base["obs"]}, models={})
    real = {sc: v["REAL"] for sc, v in base["epi"].items()}
    con[ds]["real"] = real

    def harvest(d, model):
        m = d["models"][model]
        fams = [f for f in FAM if f in m["score"] and m["score"][f] is not None]
        if ds == "SFHH":
            fams = [f for f in fams if f != "S6_groups"]
        ratios, dA, dP, passed = [], [], [], 0
        for sc, v in d["epi"].items():
            ratios.append(v[model]["R_index"] / v["REAL"]["R_index"])
            dA.append(v[model]["attack"] - v["REAL"]["attack"])
            dP.append(v[model]["P_major"] - v["REAL"]["P_major"])
            passed += bool(v[model]["pass"])
        return dict(S_pass=sum(bool(m["score"][f]) for f in fams), S_n=len(fams),
                    S_failed=[f.split("_")[0] for f in fams if not m["score"][f]],
                    E_scen=passed, ratios=ratios, dA=dA, dP=dP, stats=m["stats"],
                    contact_time_ratio=m.get("epi_total_ratio"))

    for model in ("WM", "BLOCK", "RW0", "RWhet"):
        con[ds]["models"][model] = harvest(base, model)
    con[ds]["models"]["REALtrain"] = dict(
        S_pass=sum(bool(base["models"]["REALtrain"]["score"].get(f)) for f in FAM
                   if not (ds == "SFHH" and f == "S6_groups")))
    for model, tags in EXTRA.items():
        tag = tags["dev"] if STAGE[ds] == "s1" else tags["conf"]
        d = load(A2O + f"{tag}_{ds}.json")
        con[ds]["models"][model] = harvest(d, model)
for ds in DS:
    put(f"con.dataset.{ds}", {k: con[ds][k] for k in ("N", "n_days", "n_train", "n_test", "n_groups", "M", "p", "mean_dur")},
        A2O + f"{STAGE[ds]}_{ds}.json (N, n_days, train, test, n_groups, cal)")
    put(f"con.obs.{ds}", con[ds]["obs"], A2O + f"{STAGE[ds]}_{ds}.json (obs)")
summ = {}
for model in MODS:
    cells_ratio = [r for ds in DS for r in con[ds]["models"][model]["ratios"]]
    cells_dA = [r for ds in DS for r in con[ds]["models"][model]["dA"]]
    cells_dP = [r for ds in DS for r in con[ds]["models"][model]["dP"]]
    per_ds = [sum(con[ds]["models"][model]["ratios"]) / 4 for ds in DS]
    summ[model] = dict(
        S=[con[ds]["models"][model]["S_pass"] for ds in DS],
        E_scen=[con[ds]["models"][model]["E_scen"] for ds in DS],
        E_datasets=sum(con[ds]["models"][model]["E_scen"] >= 3 for ds in DS),
        scen_total=sum(con[ds]["models"][model]["E_scen"] for ds in DS),
        scen_conf=sum(con[ds]["models"][model]["E_scen"] for ds in DS[3:]),
        ratio_per_dataset=per_ds, ratio_min=min(per_ds), ratio_max=max(per_ds),
        ratio_cell_min=min(cells_ratio), ratio_cell_max=max(cells_ratio),
        n_cells_over=sum(r > 1 for r in cells_ratio), n_cells=len(cells_ratio),
        mean_abs_log_ratio=sum(abs(math.log(r)) for r in cells_ratio) / len(cells_ratio),
        mean_abs_dA=sum(abs(x) for x in cells_dA) / len(cells_dA),
        mean_abs_dP=sum(abs(x) for x in cells_dP) / len(cells_dP),
        dA_per_dataset=[sum(con[ds]["models"][model]["dA"]) / 4 for ds in DS],
        contact_time_ratio=[con[ds]["models"][model]["contact_time_ratio"] for ds in DS])
    put(f"con.model.{model}", summ[model], A2O + "{s1,s2b,s2d,s3a,s3b,s3c,s3d}_<dataset>.json (models.<model>.score, epi)")
put("con.realtrain_S", [con[ds]["models"]["REALtrain"]["S_pass"] for ds in DS], A2O + "*.json (models.REALtrain.score)")
put("con.ladder", load(A2O + "ladder.json"), A2O + "ladder.json")
# statistics quoted in the text (observed against the homogeneous walk)
stat = {}
for ds in DS:
    o, m = con[ds]["obs"], con[ds]["models"]["RW0"]["stats"]
    stat[ds] = {k: (o[k], m[k]) for k in ("S3_deg_mean", "S4_strength_cv", "S6_within_frac", "S2_burst", "S5_persistence",
                                           "S1_timefrac_ge15", "S7_alpha", "S7_n_range") if k in o and k in m}
put("con.stats_obs_vs_RW0", stat, A2O + "{s1,s3a}_<dataset>.json (obs, models.RW0.stats)")
# re-check
vs = load(A2SUM)
put("con.check.static", vs["static_null_E_scenarios_passed"], A2SUM + " (static_null_E_scenarios_passed)")
put("con.check.static_total", dict(scenarios=sum(vs["static_null_E_scenarios_passed"].values()),
                                   datasets=sum(v >= 3 for v in vs["static_null_E_scenarios_passed"].values())),
    "sums of the entry above")
put("con.check.train_replay", vs["train_days_replay_E_scenarios_passed"], A2SUM + " (train_days_replay_E_scenarios_passed)")
put("con.check.train_replay_datasets", sum(v >= 3 for v in vs["train_days_replay_E_scenarios_passed"].values()), "count of the entry above")
put("con.check.self_power", vs["self_power_true_model_refit"], A2SUM + " (self_power_true_model_refit)")
pw = load(A2O + "01_power_table.json")
true_cells = [k for k in pw if k.split("|")[0] == k.split("|")[1] or k in ("BLOCK|WM", "RWhet|RW0")]
put("con.power.cells", {k: dict(nS=v["nS"], E=v["epi_pass_count"]) for k, v in pw.items()}, A2O + "01_power_table.json")

# table: contact records (Section 7.7)
NAME = {"WM": ("Well mixed", 2), "BLOCK": ("Constant within/between-group ratio", 3),
        "RW0": ("Lattice walk, homogeneous", 2), "RWhet": ("Lattice walk, two zones", 4),
        "RWclus": ("Walk anchored to clustered home sites$^{b}$", 4),
        "RWstickZ": ("\\quad with a slow seating zone$^{b}$", 5), "RWstickO": ("\\quad with a slow own seat$^{b}$", 5)}
rows = []
for model in ("WM", "BLOCK", "RW0", "RWhet"):
    s = summ[model]
    rows.append(f"{NAME[model][0]} & {NAME[model][1]} & {', '.join(map(str, s['S']))} & {s['scen_total']} & "
                f"{s['ratio_min']:.2f}--{s['ratio_max']:.2f} & {s['mean_abs_dA']:.3f} \\\\")
rows.append("\\addlinespace")
st = N["con.check.static_total"]["value"]
rows.append(f"Static network of training-day pair rates$^{{a}}$ & -- & -- & {st['scenarios']} & -- & -- \\\\")
for model in ("RWclus", "RWstickZ", "RWstickO"):
    s = summ[model]
    rows.append(f"{NAME[model][0]} & {NAME[model][1]} & {', '.join(map(str, s['S']))} & {s['scen_total']} & "
                f"{s['ratio_min']:.2f}--{s['ratio_max']:.2f} & {s['mean_abs_dA']:.3f} \\\\")
(OUT / "v2_tab_contacts.tex").write_text(
    "% generated by code/v2_numbers.py from " + A2O + "*.json and " + A2SUM + "\n" + "\n".join(rows) + "\n")

# table: the six datasets (Supplementary Section S9)
DCITE = {"InVS13": "genois2015data", "InVS15": "genois2018can", "LyonSchool": "stehle2011high",
         "Thiers13": "mastrandrea2015contact", "LH10": "vanhems2013estimating", "SFHH": "genois2018can"}
rows = []
for ds in DS:
    c = con[ds]
    o = c["obs"]
    rw = c["models"]["RW0"]["stats"]
    alpha = f"{o['S7_alpha']:.2f}" if o.get("S7_n_range", 0) >= 1.5 else "n.a."
    rows.append(f"{DSNAME[ds][0].upper() + DSNAME[ds][1:]} \\citep{{{DCITE[ds]}}} & "
                f"{'dev.' if STAGE[ds] == 's1' else 'conf.'} & {c['N']} & {c['n_train']}+{c['n_test']} & "
                f"{c['M']:,.0f} & {o['S3_deg_mean']:.1f} / {rw['S3_deg_mean']:.1f} & "
                f"{o['S6_within_frac']:.2f} / {rw['S6_within_frac']:.2f} & {o['S2_burst']:.2f} / {rw['S2_burst']:.2f} & "
                f"{alpha} / {rw['S7_alpha']:.2f} \\\\".replace(",", "{,}", 0))
rows = [r.replace("nan / nan", "-- / --") for r in rows]
(OUT / "v2_tab_datasets.tex").write_text(
    "% generated by code/v2_numbers.py from " + A2O + "{s1,s3a}_<dataset>.json (N, train, test, cal.M, obs, models.RW0.stats)\n"
    + "\n".join(rows) + "\n")

# ================================================================================== 3. tracer measurements
ker = load(A3R + "01_kernels.json")
ev3 = load(A3R + "04_evaluation.json")
de3 = load(A3R + "05_descriptive.json")
pw3 = load(A3R + "03_power.json")
ex3 = load(A3R + "E2_ext_evaluation.json")
put("tr.D.coach", dict(D=ker["coach"]["D_hat"], lo=ker["coach"]["boot"]["q05"], hi=ker["coach"]["boot"]["q95"],
                       n=ker["coach"]["n_measured"], D_cfd=ker["coach"]["D_hat_cfd_S4A"]), A3R + "01_kernels.json (coach)")
put("tr.D.cabin", dict(D_5L=ker["cabin"]["zeros_kept"]["fits"]["5L"]["D_hat"], D_5A=ker["cabin"]["zeros_kept"]["fits"]["5A"]["D_hat"],
                       D_geomean=ker["cabin"]["zeros_kept"]["D_geomean"],
                       D_5L_zeros_missing=ker["cabin"]["zeros_missing"]["fits"]["5L"]["D_hat"],
                       D_5A_zeros_missing=ker["cabin"]["zeros_missing"]["fits"]["5A"]["D_hat"]),
    A3R + "01_kernels.json (cabin)")
put("tr.D.train", dict(middle=ker["train"]["B_middle"]["D_hat"], middle_lo=ker["train"]["B_middle"]["boot"]["q05"],
                       middle_hi=ker["train"]["B_middle"]["boot"]["q95"], end=ker["train"]["A_end"]["D_hat"],
                       end_lo=ker["train"]["A_end"]["boot"]["q05"], end_hi=ker["train"]["A_end"]["boot"]["q95"],
                       len_middle=ker["train"]["B_middle"]["mixing_length_m"], len_end=ker["train"]["A_end"]["mixing_length_m"]),
    A3R + "01_kernels.json (train)")
put("tr.D.table", de3["D_table"], A3R + "05_descriptive.json (D_table)")
put("tr.restaurant_tracer", de3["R1_tracer"], A3R + "05_descriptive.json (R1_tracer)")
put("tr.dose_ratio", de3["dose_ratio_needed"], A3R + "05_descriptive.json (dose_ratio_needed)")
put("tr.emission", de3["implied_emission_quanta_per_h"], A3R + "05_descriptive.json (implied_emission_quanta_per_h)")
EV3 = ["T2", "T1", "C1", "T5", "R1"]
for e in EV3:
    x = ev3["events"][e]
    put(f"tr.event.{e}", dict(n_near=x["n_near"], n_far=x["n_far"], K=x["K"], k_obs=x["k_obs"], rr_obs=x["rr_obs"],
                              models={m: {k: v[k] for k in ("p_two_sided", "logscore") if k in v} | (
                                  {"rr_expected": v["rr_expected"]} if "rr_expected" in v else {}) for m, v in x["models"].items()}),
        A3R + f"04_evaluation.json (events.{e})")
put("tr.restaurant_byK", {k: {m: v["models"][m]["p_two_sided"] for m in v["models"]} for k, v in ev3["events"]["R1_byK"].items()},
    A3R + "04_evaluation.json (events.R1_byK)")
for k in ev3:
    if k not in ("events", "stamp", "observed_near"):
        put(f"tr.eval.{k}", ev3[k], A3R + f"04_evaluation.json ({k})")
put("tr.power", dict(truths=pw3["truths"], c95=pw3["c95"]), A3R + "03_power.json (truths, c95)")
put("tr.trains", {m: {k: v[k] for k in ("loglik", "deviance", "p_boot", "delta_vs_M0") if k in v}
                  for m, v in ev3["events"]["T4"]["models"].items()}, A3R + "04_evaluation.json (events.T4.models)")
put("tr.trains_rows", {r: {k: v[k] for k in ("k", "n", "P", "M0", "CRR3", "M2cal_mode")} for r, v in ev3["events"]["T4"]["row_profile"].items()},
    A3R + "04_evaluation.json (events.T4.row_profile: observed k and expected cases by rows apart)")
ver3 = load(A3R + "05_verify.json")
put("tr.pooled_without", {k: ver3[k] for k in ("pooled_P_vs_M0", "pooled_P_vs_M0_without_T5", "pooled_P_vs_M0_without_C1",
                                               "pooled_P_vs_M0_without_R1")}, A3R + "05_verify.json")
put("tr.callcentre", {m: {k: v[k] for k in ("P_le_obs_pooled", "P_le_obs_equal_weight", "int95", "median")}
                      for m, v in ev3["events"]["O1"]["models"].items()} | {"z_obs": ev3["events"]["O1"]["z_obs"]},
    A3R + "04_evaluation.json (events.O1)")
put("tr.ext.events", {e: dict(near=v["near"], far=v["far"], obs=v["obs"], sizes=v["sizes"],
                              models={m: dict(p=w["p_two_sided"], logp=w["logp"], rr=w["rr_expected"]) for m, w in v["models"].items()})
                      for e, v in ex3["events"].items()}, A3R + "E2_ext_evaluation.json (events)")
put("tr.ext.pooled", {m: {c: {k: ex3["pooled"][m][c][k] for k in ("delta", "p_under_competitor", "per_event")}
                          for c in ("M0", "CRR2", "CRR3")} | {"decision": ex3["pooled"][m]["decision"]}
                      for m in ("P", "P_econ")}, A3R + "E2_ext_evaluation.json (pooled)")
for k in ex3:
    if k not in ("events", "pooled", "stamp"):
        put(f"tr.ext.{k}", ex3[k], A3R + f"E2_ext_evaluation.json ({k})")

# probability of the pooled pattern that was observed (better than well mixed, D2, and not better than both
# constant ratios, not D3) under each truth: exact enumeration over the joint outcomes of the five stratum events
# (the computation of the check script verify/adv/a06_power.py, vectorised)
_pr = load(A3R + "02_predictions.json")["events"]
_EV, _TR = ["T2", "T1", "T5", "C1", "R1"], ["P", "M0", "CRR2", "CRR3"]
_pm = {e: {m: np.array(_pr[e]["pmf"][m], float) for m in _TR} for e in _EV}
_g = np.meshgrid(*[np.arange(len(_pm[e]["P"])) for e in _EV], indexing="ij")
_W = {m: np.prod([_pm[e][m][_g[a]] for a, e in enumerate(_EV)], axis=0).ravel() for m in _TR}
_LG = {m: sum(np.log(np.maximum(_pm[e][m][_g[a]], 1e-300)) for a, e in enumerate(_EV)).ravel() for m in _TR}


def _pv(j):
    d = _LG["P"] - _LG[j]
    o = np.argsort(-d)
    cw = np.cumsum(_W[j][o])
    idx = np.searchsorted(-d[o], -d - 1e-12, side="right") - 1
    return d, cw[idx]


_d0, _p0 = _pv("M0"); _d2, _p2 = _pv("CRR2"); _d3, _p3 = _pv("CRR3")
_D2 = (_d0 > 0) & (_p0 < 0.05)
_D3 = (_d2 > 0) & (_p2 < 0.05) & (_d3 > 0) & (_p3 < 0.05)
_pat = {m: float(_W[m][_D2 & ~_D3].sum()) for m in _TR}
put("tr.power.pattern", dict(P_D2_not_D3=_pat, ratio_CRR2_over_P=_pat["CRR2"] / _pat["P"],
                             ratio_CRR3_over_P=_pat["CRR3"] / _pat["P"]),
    A3R + "02_predictions.json (events.*.pmf), enumerated as in code/checks/tracer/separate/a06_power.py "
    "(added by the re-check after the results were known)")

# table: tracer kernel on the first set of events (Section 7.8); rr_pred of the kernel from 05_descriptive / 04_evaluation
def rr_pred(e):
    return ev3["events"][e]["models"]["P"].get("rr_expected")


TN = {"T2": "Hunan coach", "T1": "Zhejiang bus", "C1": "Marin classroom", "T5": "Flight VN54",
      "R1": "Restaurant ($K=2$)"}
KSRC = {"T2": "tracer-validated CFD, same coach", "T1": "$\\Dair$ from the coach tracer", "C1": "room relation",
        "T5": "cabin tracer, two releases", "R1": "tracer gas, same room"}
rows = []
for e in EV3:
    x = ev3["events"][e]
    m = x["models"]
    rp = rr_pred(e)
    obs = f"{x['k_obs']}/{x['n_near']} vs {x['K'] - x['k_obs']}/{x['n_far']}"
    rows.append(" & ".join([
        TN[e], obs, (f"{rp:.2f}" if rp is not None else "--") + f" / {x['rr_obs']:.2f}" + ("$^{a}$" if e == "R1" else ""),
        pval(m["P"]["p_two_sided"]), pval(m["M0"]["p_two_sided"]),
        sgn(m["P"]["logscore"] - m["M0"]["logscore"]), sgn(m["P"]["logscore"] - m["CRR2"]["logscore"]),
        sgn(m["P"]["logscore"] - m["CRR3"]["logscore"])]) + " \\\\")
(OUT / "v2_tab_tracer.tex").write_text(
    "% generated by code/v2_numbers.py from " + A3R + "04_evaluation.json (events) and results/05_descriptive.json\n"
    + "\n".join(rows) + "\n")
tot = {c: sum(ev3["events"][e]["models"]["P"]["logscore"] - ev3["events"][e]["models"][c]["logscore"] for e in EV3)
       for c in ("M0", "CRR2", "CRR3")}
put("tr.pooled_sum", tot, A3R + "04_evaluation.json (sum over T2, T1, T5, C1, R1 of logscore differences)")
put("tr.crr_vs_M0", {c: sum(ev3["events"][e]["models"][c]["logscore"] - ev3["events"][e]["models"]["M0"]["logscore"] for e in EV3)
                     for c in ("CRR2", "CRR3")}, A3R + "04_evaluation.json (same sum for the constant ratios)")

# ======================================================================================== 4. zone level
mp = load(A4R + "matsumoto_primary.json")
sf = load(A4R + "seoul_floors.json")
ph = load(A4R + "posthoc_descriptive.json")
ps = load(A4R + "power_summary.json")
vm = load(A4V + "v_main.json")
vf = load(A4V + "v_floor.json")
pr = mp["primary"]
put("zone.counts", dict(n=mp["n_students"], infected=mp["info"]["n_inf"], classes=mp["n_classes"], schools=ph["check_counts"]["schools"]),
    A4R + "matsumoto_primary.json (n_students, info.n_inf, n_classes); results/posthoc_descriptive.json (check_counts)")
put("zone.heldout_ll", pr["heldout_ll"], A4R + "matsumoto_primary.json (primary.heldout_ll)")
for k in ("D_Z_minus_H", "D_Z_minus_X", "D_Z_minus_C"):
    put(f"zone.{k}", dict(value=pr[k], ci=pr[k + "_ci"]), A4R + f"matsumoto_primary.json (primary.{k}, {k}_ci)")
rest = {k: v for k, v in pr.items() if not k.startswith("D_Z") and k not in ("heldout_ll", "cvfit")
        and not (isinstance(v, list) and len(v) > 12)}
put("zone.primary_other", rest, A4R + "matsumoto_primary.json (primary.*)")
for k in mp:
    if k not in ("primary", "info", "n_students", "n_classes", "onset_range"):
        put(f"zone.{k}", mp[k], A4R + f"matsumoto_primary.json ({k})")
put("zone.posthoc", {k: ph[k] for k in ("schools_Z_ahead_of_H", "schools_Z_ahead_of_X", "schools_Z_ahead_of_C",
                                         "rho_fitted_C", "rho_lyon_own_school", "rho_implied_by_Lyon_per_class_quantiles",
                                         "AIC", "check_dispersion", "leave_one_out_range_C", "max_school_contrib_C")},
    A4R + "posthoc_descriptive.json")
put("zone.power", ps["table"], A4R + "power_summary.json (table)")
put("zone.seoul", {k: sf[k] for k in ("N_other", "obs_other", "pred_count_median", "pred_count_interval95",
                                       "well_mixed_expected", "well_mixed_P_le_3", "mu_quantiles", "central")},
    A4R + "seoul_floors.json")
put("zone.check.heldout", vm["heldout"], A4V + "v_main.json (heldout)")
for k in vm:
    if k not in ("heldout", "cv") and not (isinstance(vm[k], list) and len(vm[k]) > 40):
        put(f"zone.check.{k}", vm[k], A4V + f"v_main.json ({k})")
put("zone.check.hazard_bound", {b: {k: v[k] for k in ("D_H", "D_X", "D_N702", "D_N811", "D_N2/3", "TV", "LR") if k in v}
                                for b, v in vf.items() if isinstance(v, dict)}, A4V + "v_floor.json")

# ===================================================================================== 5. freeze records
import time  # noqa: E402


def mtime(rel):
    """Modification time of a file of the working tree (local time); None in a copy without the original times."""
    p = ROOT / rel
    return time.strftime("%H:%M:%S", time.localtime(p.stat().st_mtime)) if p.exists() else None


def check(base, hashes):
    """Hashes of a freeze record against the files; source files that are not redistributed may be absent."""
    present = {k: h for k, h in hashes.items() if (ROOT / base / k).exists()}
    return dict(n_files=len(hashes), n_present=len(present), all_match=all(sha(base + k) == h for k, h in present.items()))


# The file times and the full hash check exist only in the author's working tree (a copy or a clone does not keep
# modification times, and some hashed files are third-party sources that are not redistributed).  They are therefore
# recorded once, with --record-freeze, in data/v2_freeze_evidence.json; an ordinary run reads that record and
# re-checks the hashes of the files that are present.
EVID = OUT / "v2_freeze_evidence.json"


def build_freeze():
    frz = {}
    f1 = load(A1 + "PREREG_FREEZE.json")
    frz["outbreaks"] = dict(frozen_at=f1["frozen_at"], addendum_at=f1["addendum"]["written_at"],
                            protocol=check(A1, f1["files"]), addendum=check(A1, f1["addendum"]["files"]),
                            times=dict(protocol=mtime(A1 + "PREREGISTRATION.md"), calibration=mtime(A1O + "02_calibration.json"),
                                       power=mtime(A1O + "03_power.json"), holdout=mtime(A1O + "04_holdout.json"),
                                       corrected=mtime(A1O + "10_corrected.json")))
    f2 = load(A2 + "PREREG_FREEZE.json")
    f2b = load(A2 + "FIX_FREEZE.json")
    frz["contacts"] = dict(frozen_at_utc=f2["frozen_at_utc"], second_freeze_utc=f2b["frozen_at_utc"],
                           protocol=check(A2, f2["sha256"]),
                           second=check(A2, {k: h for k, h in f2b["sha256"].items() if k != "POSTHOC_LOG.md"}),
                           note="the second record also holds the hash of the post-hoc log up to its item 12; the log was "
                                "continued afterwards, so that hash is not compared here",
                           times=dict(protocol=mtime(A2 + "PREREGISTRATION.md"), first_result=mtime(A2O + "s1_InVS13.json"),
                                      first_confirmatory=mtime(A2O + "s3a_InVS15.json")))
    f3 = load(A3 + "PREREG_FREEZE.json")
    f3b = load(A3 + "PREREG_FREEZE_2_before_unblinding.json")
    f3c = load(A3 + "PREREG_EXT_FREEZE.json")
    f3d = load(A3 + "PREREG_EXT_FREEZE_2_before_unblinding.json")
    frz["tracer"] = dict(frozen_at=f3["frozen_at"], before_unblinding_at=f3b["written_at"],
                         protocol=check(A3, {**f3["files"], "PREREGISTRATION.md": f3["prereg_sha256"]}),
                         before_unblinding=check(A3, f3b["files"]),
                         extension_at=f3c["frozen_at"], extension=check(A3, {**f3c["files"], "PREREG_EXT_flights.md": f3c["prereg_ext_sha256"]}),
                         extension_before_unblinding_at=f3d["written_at"], extension_before_unblinding=check(A3, f3d["files"]),
                         times=dict(protocol=mtime(A3 + "PREREGISTRATION.md"), kernels=mtime(A3R + "01_kernels.json"),
                                    evaluation=mtime(A3R + "04_evaluation.json"),
                                    extension_protocol=mtime(A3 + "PREREG_EXT_flights.md"),
                                    extension_evaluation=mtime(A3R + "E2_ext_evaluation.json")))
    f4 = load(A4 + "PREREG_FREEZE.json")
    frz["zones"] = dict(frozen_at=f4["frozen_at_local"], protocol=check(A4, f4["sha256"]),
                        times=dict(protocol=mtime(A4 + "PREREGISTRATION.md"), schools=mtime(A4R + "matsumoto_primary.json"),
                                   floors=mtime(A4R + "seoul_floors.json")))
    return frz


if "--record-freeze" in sys.argv:
    frz = build_freeze()
    assert all(v["n_present"] == v["n_files"] and v["all_match"]
               for a in frz.values() for v in a.values() if isinstance(v, dict) and "n_files" in v), "freeze check"
    EVID.write_text(json.dumps(frz, indent=1, ensure_ascii=False) + "\n")
    print("freeze evidence recorded ->", EVID)
else:
    now = build_freeze()
    bad = [(a, k) for a, rec in now.items() for k, v in rec.items()
           if isinstance(v, dict) and "n_files" in v and not v["all_match"]]
    assert not bad, f"hash mismatch against the freeze records: {bad}"
    n_here = sum(v["n_present"] for rec in now.values() for v in rec.values() if isinstance(v, dict) and "n_files" in v)
    n_all = sum(v["n_files"] for rec in now.values() for v in rec.values() if isinstance(v, dict) and "n_files" in v)
    red = (ROOT / "code" / "REDACTIONS.json").exists()
    print(f"freeze records: {n_here} of {n_all} hashed files present, all match"
          + (" (redacted copies: through the hashes of the originals listed in code/REDACTIONS.json, against which "
             "each copy was checked; the originals themselves are not in this repository)" if red else
             " (hashes recomputed on the original files)"))
    frz = json.loads(EVID.read_text())
put("freeze", frz, "data/v2_freeze_evidence.json <- PREREG_FREEZE.json, FIX_FREEZE.json, "
    "PREREG_FREEZE_2_before_unblinding.json and PREREG_EXT_FREEZE*.json of the four analyses; SHA-256 of every hashed "
    "file recomputed on the unredacted original in the author's working tree, and modification times of the files there "
    "(British Summer Time) read, when that record was written (v2_numbers.py --record-freeze); the protocol files of the "
    "public repository are redacted copies whose own hashes differ (code/REDACTIONS.md)")

# ============================================================== 6. first set: numbers used by the summary figure
h1 = load(R1D + "03_holdout_shape.json")
# the hold-out scores of the first test with the round-off-safe propagator (only the classroom changes; the
# other events are identical to the stored ones): data/s8_first_test_exact.json
fx = load("data/s8_first_test_exact.json")["primary"]["events"]
for e in h1["events"]:
    assert abs(fx[e]["M0"]["log_score"] - h1["events"][e]["M0"]["log_score"]) < 1e-12
    h1["events"][e]["M2"]["log_score"] = fx[e]["M2"]["log_score"]
first = {}
for e, a3e in (("T1", "T1"), ("T2", "T2"), ("T5", "T5"), ("C1", "C1")):
    m2 = h1["events"][e]["M2"]["log_score"]
    m0 = h1["events"][e]["M0"]["log_score"]
    c2 = ev3["events"][a3e]["models"]["CRR2"]["logscore"]
    c3 = ev3["events"][a3e]["models"]["CRR3"]["logscore"]
    assert abs(m0 - ev3["events"][a3e]["models"]["M0"]["logscore"]) < 1e-9          # same test in both analyses
    first[e] = dict(d0=m2 - m0, dC2=m2 - c2, dC3=m2 - c3)
put("first.per_event", first, "data/s8_first_test_exact.json (primary.events.*.M2, M0: log_score; equal to "
    + R1D + "03_holdout_shape.json except for the classroom); constant ratios 2 and 3 on the same "
    "strata: " + A3R + "04_evaluation.json (events.*.models.CRR2, CRR3: logscore)")
put("first.pooled", dict(d0=sum(v["d0"] for v in first.values()), dC2=sum(v["dC2"] for v in first.values()),
                         dC3=sum(v["dC3"] for v in first.values())), "sums of the entry above")

(OUT / "v2_numbers.json").write_text(json.dumps(N, indent=1, ensure_ascii=False, default=float) + "\n")
print(f"{len(N)} entries -> {OUT / 'v2_numbers.json'}")
for f in sorted(OUT.glob("v2_tab_*.tex")):
    print("table body:", f.name)

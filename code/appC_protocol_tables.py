"""appC_protocol_tables.py -- numbers and table rows of Supplementary Section S8 (protocol of the first outbreak test, sensitivity analyses,
post-hoc log).

Nothing is simulated here.  The script
  1. recomputes the SHA-256 of the six files frozen in code/bsc_validation/PREREG_FREEZE.json and compares
     them with the recorded hashes, and reads the file-system times that are the only evidence for the order
     "protocol frozen -> first model run -> power addendum -> hold-out comparison";
  2. reads the stored results of the first outbreak test (data/bsc_validation/*.json), of its
     re-run with the round-off-safe propagator (data/s8_first_test_exact.json, written by
     s8_first_test_exact.py: the M2 row of Table (b) of tab:power and every row of tab:sensitivity come from
     it, because the first analysis was affected by a numerical defect in the classroom) and of the
     re-check (data/validation_checks/*) and writes
        ../data/appC_protocol_numbers.json   every number quoted in sections/supp_protocol1.tex
        ../data/appC_protocol_tables.tex     the LaTeX rows of Tables appC_protocol:tab:power, :tab:sensitivity
                                             and :tab:levels (pasted into the section)
No table entry of the appendix is typed by hand.  The script fails (exit status 1) if a generated row is not
found verbatim in sections/supp_protocol1.tex.  Block (a) of tab:power is the power computed before unblinding
(first analysis); its values recomputed with the corrected propagator are stored as power_corrected.

    python appC_protocol_tables.py            (about 1 s)
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                       # repository root
VAL = ROOT / "code/bsc_validation"
VER = HERE.parent / "data" / "validation_checks"
OUT = HERE.parent / "data"


def load(name):
    return json.loads((VAL / "data" / name).read_text())


EXACT = json.loads((OUT / "s8_first_test_exact.json").read_text())   # first test with the corrected propagator


def stamp(ts):
    return dt.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")


def fmt_p(p):
    """p-value for a table cell: two decimals above 0.1, two significant digits below, a x 10^b below 1e-3."""
    if p >= 0.995:
        return "1.00"
    if p >= 0.0995:
        return f"{p:.2f}"
    if p >= 0.00995:
        return f"{p:.3f}"
    if p >= 0.000995:
        return f"{p:.4f}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def fmt_p3(p):
    """p-value printed by the re-check with three decimals (no further digits are available)."""
    if p < 0.0005:
        return "$<0.001$"
    if p >= 0.995:
        return "1.00"
    return f"{p:.2f}" if p >= 0.0995 else f"{p:.3f}"


def fmt_rr(x):
    return f"{x:.1f}" if x < 9.95 else f"{x:.0f}"


def fmt_delta(x):
    return f"${'+' if x >= 0 else '-'}{abs(x):.2f}$"


def fmt_D(x):
    if x >= 1000:
        return f"{x:,.0f}".replace(",", "{,}")
    if x >= 10:
        return f"{x:.0f}"
    if x >= 1:
        return f"{x:.1f}"
    return f"{x:.2g}"


def frozen_sha256(path) -> tuple[str, str | None]:
    """SHA-256 of a frozen protocol file, and the hash of its unredacted original if the file is a redacted copy.
    In the author's working tree the files are the originals: the first value is their hash and the second is
    None.  In the public repository the protocol records are redacted copies (code/REDACTIONS.md): the copy is
    checked against the hash listed for it in code/REDACTIONS.json, and the hash of the original listed there is
    returned as the second value; the comparison with the freeze record then rests on that list."""
    path = Path(path)
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    listing = ROOT / "code" / "REDACTIONS.json"
    if listing.exists():
        red = json.loads(listing.read_text())
        try:
            key = Path(os.path.normpath(path.absolute())).relative_to(os.path.normpath(ROOT.absolute())).as_posix()
        except ValueError:
            return h, None
        if key in red:
            assert h == red[key]["sha256_copy"], f"redacted copy differs from code/REDACTIONS.json: {key}"
            return h, red[key]["sha256_original"]
    return h, None


num = {}
tex = []

# ---------------------------------------------------------------------------------------------------------
# 1. freeze record: hashes and file-system times
# ---------------------------------------------------------------------------------------------------------
frz = json.loads((VAL / "PREREG_FREEZE.json").read_text())
hashes = {}
for rel, h in frz["sha256"].items():
    own, listed = frozen_sha256(VAL / rel)
    if listed is None:                        # the original itself (author's working tree)
        hashes[rel] = dict(recorded=h, sha256_recomputed_on_original=own, match=bool(own == h))
    else:                                     # redacted copy (public repository)
        hashes[rel] = dict(recorded=h, sha256_copy_recomputed=own, sha256_original_listed=listed,
                           match=bool(listed == h))
basis = ("hashes recomputed on the unredacted original files in the author's working tree"
         if all("sha256_recomputed_on_original" in v for v in hashes.values()) else
         "redacted copies: each copy recomputed and checked against code/REDACTIONS.json; the comparison with the "
         "freeze record uses the hashes of the originals listed there, which cannot be recomputed from the repository")
times = {}
for rel in ["PREREGISTRATION.md", "PREREG_FREEZE.json", "PREREG_ADDENDUM_A1_power.md", "POSTHOC_LOG.md",
            "src/valmod/events.py", "src/valmod/lattice.py", "src/valmod/stats.py", "src/valmod/exact.py",
            "src/valmod/predict.py", "src/valmod/decision.py", "src/valmod/level.py",
            "data/heldout_outcomes.json", "logs/01_primary.log", "data/01_cal_T4_primary.json",
            "data/02_predictive.json", "data/02_power.json", "data/03_holdout_shape.json",
            "scripts/02_power.py", "scripts/03_holdout_shape.py"]:
    st = os.stat(VAL / rel)
    times[rel] = dict(created=stamp(getattr(st, "st_birthtime", st.st_mtime)), modified=stamp(st.st_mtime))
first_log_line = (VAL / "logs/01_primary.log").read_text().splitlines()[0]
num["freeze"] = dict(
    frozen_at=frz["frozen_at"], verified_unchanged_at=frz["verified_unchanged_at"],
    hashes=hashes, hash_basis=basis, all_hashes_match=all(v["match"] for v in hashes.values()), n_hashed_files=len(hashes),
    file_times_local=times, first_calibration_log_line=first_log_line,
    note="file-system times as read when this script was run; they are metadata, not a version-control record")

# ---------------------------------------------------------------------------------------------------------
# 2. power of the decision rule (computed before unblinding) and the generic constant-risk-ratio model
# ---------------------------------------------------------------------------------------------------------
pw = load("02_power.json")
num["power"] = {j: {t: {k: pw[j][t][k] for k in ("P_A1", "P_A2_nofail", "P_A2_le1", "P_W1_like")}
                    for t in ("M0", "M1", "M2")} for j in ("M1", "M2")}
num["power"]["per_event_reject_M0_if_M2_true"] = {e: pw["M2_per_event_power"][e]["P_reject_M0_if_model_true"]
                                                  for e in ("T1", "T2", "T5", "C1")}
num["power"]["level_dispersion"] = pw["level_dispersion"]
num["power"]["stamp"] = pw["stamp"]
num["power_corrected"] = {j: {t: {k: EXACT["power"][j][t][k] for k in ("P_A1", "P_A2_nofail", "P_A2_le1", "P_W1_like")}
                              for t in ("M0", "M1", "M2")} for j in ("M1", "M2")}


def fmt_prob(p):
    if p < 5e-5:
        return "$<10^{-4}$"
    if p < 0.01:
        return f"{p:.4f}"
    return f"{p:.3f}"


rows = ["% --- Table appC_protocol:tab:power, block (a): data/bsc_validation/02_power.json"]
labels = [("P_A1", "A1 passes"), ("P_A2_nofail", "no A2 failure"), ("P_A2_le1", "at most one A2 failure"),
          ("P_W1_like", "A1 and no A2 failure")]
for j, name in (("M1", "$\\Mone$"), ("M2", "$\\Mtwo$")):
    for i, (k, lab) in enumerate(labels):
        rows.append(f"    {name if i == 0 else ''} & {lab} & " +
                    " & ".join(fmt_prob(pw[j][t][k]) for t in ("M0", "M1", "M2")) + " \\\\")
    if j == "M1":
        rows.append("    \\addlinespace")
tex += rows

gen = []
txt = (VER / "v05_generic.txt").read_text()
m = re.search(r"pooled \(np\.float64\(([-\d.e]+)\), np\.float64\(([-\d.e]+)\)\)", txt)
num["generic"] = dict(M2_check_reimplementation=dict(delta=float(m.group(1)), p=float(m.group(2))), models=[])
for line in txt.splitlines():
    mm = re.match(r"constant RR=([\d.]+): pooled delta ([+-][\d.]+), p under M0 ([\d.]+); A2 p-values \[(.*?)\] "
                  r"failures \[(.*?)\]", line)
    if mm:
        ps = [float(x) for x in mm.group(4).split()]
        fails = [s.strip(" '") for s in mm.group(5).split(",") if s.strip()]
        num["generic"]["models"].append(dict(rr=float(mm.group(1)), delta=float(mm.group(2)), p=float(mm.group(3)),
                                             p_T1=ps[0], p_T2=ps[1], p_T5=ps[2], p_C1=ps[3], failures=fails))
sh = load("03_holdout_shape.json")
dec = load("10_decision.json")
EVN = {"T1": "bus", "T2": "coach", "T5": "flight", "C1": "classroom"}
rows = ["% --- Table appC_protocol:tab:power, block (b): data/validation_checks/v05_generic.txt;",
        "%     M2 row: data/s8_first_test_exact.json (primary.pooled.M2.all, primary.events.*.M2.p_two_sided)"]
m2 = dec["models"]["M2"]                                   # first analysis (call centre, A4): not affected
pm2 = EXACT["primary"]["pooled"]["M2"]
fails_m2 = [EVN[e] for e in EVN if EXACT["primary"]["events"][e]["M2"]["p_two_sided"] < 0.05]
rows.append("    $\\Mtwo$ (calibrated) & " + fmt_delta(pm2["all"]["delta"]) + " & " + fmt_p(pm2["all"]["p_under_M0"]) + " & " +
            " & ".join(fmt_p(EXACT["primary"]["events"][e]["M2"]["p_two_sided"]) for e in EVN) + " & " +
            (", ".join(fails_m2) or "none") + " \\\\")
num["M2_first_test_corrected"] = dict(delta=pm2["all"]["delta"], p=pm2["all"]["p_under_M0"],
                                      p_events={e: EXACT["primary"]["events"][e]["M2"]["p_two_sided"] for e in EVN},
                                      first_analysis=dict(delta=m2["A1"]["delta"], p=m2["A1"]["p"]))
rows.append("    \\addlinespace")
for g in num["generic"]["models"]:
    rr = f"{g['rr']:g}"
    pp = f"{g['p']:.5f}"                       # printed with five decimals in the source
    m_, e_ = f"{g['p']:.0e}".split("e")
    rows.append(f"    constant ratio {rr} & " + fmt_delta(g["delta"]) + f" & ${m_}\\times10^{{{int(e_)}}}$ & " +
                " & ".join(fmt_p3(g[k]) for k in ("p_T1", "p_T2", "p_T5", "p_C1")) + " & " +
                (", ".join(EVN[f] for f in g["failures"]) or "none") + " \\\\")
tex += rows
passing = [g for g in num["generic"]["models"] if g["delta"] > 0 and g["p"] < 0.05 and not g["failures"]]
num["generic"]["ratios_passing_A1_with_no_A2_failure"] = [g["rr"] for g in passing]
num["generic"]["delta_range_of_those"] = [min(g["delta"] for g in passing), max(g["delta"] for g in passing)]

# ---------------------------------------------------------------------------------------------------------
# 3. sensitivity analyses (12 variants, primary included)
# ---------------------------------------------------------------------------------------------------------
sens_first = load("09_sensitivity.json")["variants"]       # first analysis (classroom affected by the defect)
sens = EXACT["sens"]                                         # corrected: these rows are printed
assert set(sens) == set(sens_first), "variant sets differ"
cal = {v: json.loads((VAL / "data" / f"01_cal_T4_{v}.json").read_text())
       for v in ("primary", "S-adj", "S-row", "S-logit", "S-pitch", "S-aniso")}
calr = {v: json.loads((VAL / "data" / f"01_cal_R1_{v}.json").read_text())
        for v in ("primary", "S-R1-upper", "S-R1-lower")}
num["calibration_variants"] = {
    **{v: dict(n_cells=c.get("n_cells"), M2_D_hat=c["M2"]["D_hat"], M1_D_hat=c["M1"]["D_hat"],
               M2_theta_hat=c["M2"].get("theta_hat"), M1_theta_hat=c["M1"].get("theta_hat")) for v, c in cal.items()},
    **{"R1_" + v: dict(M2_D_hat=c["M2"]["D_hat"], M1_D_hat=c["M1"]["D_hat"]) for v, c in calr.items()}}
nu = sens["S-nu"]["events"]["T1"]
num["S_nu_factor"] = dict(M1=nu["M1_nu"], M2=nu["M2_nu"])
DESC = {
    "primary (closed form for M1)": ("primary", "protocol as frozen; $\\Dair=%s$ (vehicles), $%s$ (rooms)"
                                     % (fmt_D(cal["primary"]["M2"]["D_hat"]), fmt_D(calr["primary"]["M2"]["D_hat"]))),
    "S-adj": ("S-adj", "adjacent-seat cell put back (%d cells); $\\Dair=%s$"
              % (cal["S-adj"]["n_cells"], fmt_D(cal["S-adj"]["M2"]["D_hat"]))),
    "S-row": ("S-row", "whole index row excluded (%d cells); $\\Dair=%s$"
              % (cal["S-row"]["n_cells"], fmt_D(cal["S-row"]["M2"]["D_hat"]))),
    "S-logit": ("S-logit", "logit-normal likelihood on printed intervals (%d cells); $\\Dair=%s$"
                % (cal["S-logit"]["n_cells"], fmt_D(cal["S-logit"]["M2"]["D_hat"]))),
    "S-pitch": ("S-pitch", "train row pitch 0.4\\,m (source) instead of 1.0\\,m; $\\Dair=%s$"
                % fmt_D(cal["S-pitch"]["M2"]["D_hat"])),
    "S-aniso": ("S-aniso", "fitted seat-back factor on hops between rows ($\\theta=%.2f$); $\\Dair=%s$"
                % (cal["S-aniso"]["M2"]["theta_hat"], fmt_D(cal["S-aniso"]["M2"]["D_hat"]))),
    "S-R1-upper": ("S-R1-upper", "restaurant counts at their upper bound (3, 2); $\\Dair=%s$"
                   % fmt_D(calr["S-R1-upper"]["M2"]["D_hat"])),
    "S-R1-lower": ("S-R1-lower", "restaurant counts at their lower bound (1, 1); $\\Dair=%s$"
                   % fmt_D(calr["S-R1-lower"]["M2"]["D_hat"])),
    "S-5K": ("S-5K", "flight: index case in seat 5K (recalled, not verified)"),
    "S-Luo": ("S-Luo", "coach: 150\\,min and 7/48 \\citep{luo2020}"),
    "S-corner": ("S-corner", "classroom: teacher's desk in a front corner"),
    "S-nu": ("S-nu", "train adjacent-seat excess (factor %.1f) applied to adjacent seats of hold-out events"
             % nu["M2_nu"]),
}
rows = ["% --- Table appC_protocol:tab:sensitivity: data/s8_first_test_exact.json (sens.*.events.*.M2, sens.*.M2,",
        "%     sens.*.M1; first analysis: data/bsc_validation/09_sensitivity.json); calibrated coefficients:",
        "%     data/bsc_validation/01_cal_*.json (M2.D_hat)"]
summ = dict(n_variants=len(sens), M2_A1_pass=[], M2_A1_fail=[], coach_fail=[], coach_pass=[], bus_fail=[],
            flight_fail=[], classroom_fail=[], M1_A1_pass=[], M2_p_under_M0=[], M2_delta={}, M1_delta={})
for tag, v in sens.items():
    short, desc = DESC[tag]
    cells = []
    for e in EVN:                       # bold: two-sided predictive p < 0.05, i.e. an A2 failure
        c = f"{fmt_rr(v['events'][e]['M2']['rr_pred'])}; {fmt_p(v['events'][e]['M2']['p'])}"
        assert (e in v["M2"]["failed"]) == (v["events"][e]["M2"]["p"] < 0.05)
        cells.append(f"{{\\bfseries\\boldmath {c}}}" if e in v["M2"]["failed"] else c)
    rows.append(f"    {short} & {desc} & " + " & ".join(cells) + " & " + fmt_delta(v["M2"]["delta"]) + " \\\\")
    (summ["M2_A1_pass"] if v["M2"]["A1"] else summ["M2_A1_fail"]).append(short)
    (summ["coach_fail"] if "T2" in v["M2"]["failed"] else summ["coach_pass"]).append(short)
    for e, key in (("T1", "bus_fail"), ("T5", "flight_fail"), ("C1", "classroom_fail")):
        if e in v["M2"]["failed"]:
            summ[key].append(short)
    if v["M1"]["A1"]:
        summ["M1_A1_pass"].append(short)
    summ["M2_p_under_M0"].append(v["M2"]["p_under_M0"])
    summ["M2_delta"][short] = v["M2"]["delta"]
    summ["M1_delta"][short] = v["M1"]["delta"]
summ["M2_p_under_M0_range"] = [min(summ["M2_p_under_M0"]), max(summ["M2_p_under_M0"])]
summ["n_M2_A1_pass"] = len(summ["M2_A1_pass"])
summ["n_coach_fail"] = len(summ["coach_fail"])
summ["coach_p"] = {DESC[t][0]: v["events"]["T2"]["M2"]["p"] for t, v in sens.items()}
summ["classroom_p"] = {DESC[t][0]: v["events"]["C1"]["M2"]["p"] for t, v in sens.items()}
summ["bus_p"] = {DESC[t][0]: v["events"]["T1"]["M2"]["p"] for t, v in sens.items()}
num["sensitivity"] = summ
num["sensitivity_first_analysis"] = dict(
    M2_delta={DESC[t][0]: v["M2"]["delta"] for t, v in sens_first.items()},
    M2_failed={DESC[t][0]: v["M2"]["failed"] for t, v in sens_first.items()},
    classroom_p={DESC[t][0]: v["events"]["C1"]["M2"]["p"] for t, v in sens_first.items()},
    note="data/bsc_validation/09_sensitivity.json, computed with the mode-sum propagator")
tex += rows

# ---------------------------------------------------------------------------------------------------------
# 4. Monte-Carlo stability (five other seeds), leave-one-event-out, common mixing coefficient
# ---------------------------------------------------------------------------------------------------------
mc = load("02c_mc_stability.json")
num["mc_stability"] = dict(
    seeds=mc["seeds"],
    M2={e: [min(r["p"][e] for r in mc["M2"]), max(r["p"][e] for r in mc["M2"])] for e in EVN},
    M2_delta=[min(r["delta"] for r in mc["M2"]), max(r["delta"] for r in mc["M2"])],
    M2_p_pooled=[min(r["p_pooled"] for r in mc["M2"]), max(r["p_pooled"] for r in mc["M2"])],
    M1cf_delta=[min(r["delta"] for r in mc["M1cf"]), max(r["delta"] for r in mc["M1cf"])],
    primary=dict(delta=sh["pooled"]["M2"]["delta"] if "pooled" in sh else m2["A1"]["delta"], p=m2["A1"]["p"]))

from math import exp  # noqa: E402

from scipy.stats import chi2  # noqa: E402

loo = load("07_loo.json")["classes"]
num["loo"] = {}
NAME = {"T4": "trains", "T1": "bus", "T2": "coach", "T5": "flight", "R1": "restaurant", "C1": "classroom"}
# leave-one-event-out numbers are quoted in the text (stored in the JSON), not tabulated
for cls, evs in (("vehicle", ("T4", "T1", "T2", "T5")), ("room", ("R1", "C1"))):
    c = loo[cls]["M2"]
    lr = 2 * (c["sum_own_max"] - c["joint_max"])
    df = len(evs) - 1
    num["loo"][cls] = dict(events={e: c["loo"][e] for e in evs}, joint_D_mode=c["joint_D_mode"],
                           sum_own_max=c["sum_own_max"], joint_max=c["joint_max"], LR=lr, df=df,
                           p=float(chi2.sf(lr, df)))
    c1 = loo[cls]["M1"]
    lr1 = 2 * (c1["sum_own_max"] - c1["joint_max"])
    num["loo"][cls + "_M1"] = dict(LR=lr1, df=df, p=float(chi2.sf(lr1, df)), joint_D_mode=c1["joint_D_mode"])
import numpy as np  # noqa: E402

_l = load("07_loo.json")
_logD = np.array(_l["logD"])
_c1 = np.array(_l["classes"]["room"]["M2"]["curves"]["C1"])
_loc = [i for i in range(1, len(_c1) - 1) if _c1[i] > _c1[i - 1] and _c1[i] > _c1[i + 1]]
num["loo"]["classroom_curve_local_maxima"] = [dict(D=float(10 ** _logD[i]), logL=float(_c1[i])) for i in _loc]
_t2 = np.array(_l["classes"]["vehicle"]["M2"]["curves"]["T2"])
num["loo"]["coach_curve"] = dict(logL_at_grid_top=float(_t2[-1]), max=float(_t2.max()),
                                 logL_drop_from_1e3_4_to_top=float(_t2.max() - _t2[np.argmin(abs(_logD - 3.4))]),
                                 grid_top_D=float(10 ** _logD[-1]))
# the re-check's recomputation of the common-D test (vehicles)
t10 = (VER / "v10_c1_loo.txt").read_text()
mm = re.search(r"sum of maxima ([-\d.]+) common-D max ([-\d.]+) at D ([\d.]+) LR ([\d.]+) p\(3 df\) ([\d.]+)", t10)
num["loo"]["vehicle_check"] = dict(sum_own_max=float(mm.group(1)), joint_max=float(mm.group(2)),
                                      joint_D=float(mm.group(3)), LR=float(mm.group(4)), p=float(mm.group(5)))
num["loo"]["vehicle_check"]["D_alone"] = {
    k: float(v) for k, v in re.findall(r"^(T\d) D preferred alone ([\d.]+)", t10, flags=re.M)}
mm = re.search(r"C1 with the train-calibrated D_air.*?p=([\d.]+)", (VER / "v11_breakdown.txt").read_text())
num["classroom_with_train_D_posthoc_p"] = float(mm.group(1))

# ---------------------------------------------------------------------------------------------------------
# 5. level checks and negative controls
# ---------------------------------------------------------------------------------------------------------
lv = load("04_level_controls.json")
ver_c3 = {}
for line in (VER / "v13_c3.txt").read_text().splitlines():
    mm = re.match(r"mask ([\d.]+): mean ([\d.]+), P\(X>=5\)=([\d.]+)", line)
    if mm:
        ver_c3[float(mm.group(1))] = dict(mean=float(mm.group(2)), P_ge5=float(mm.group(3)))
num["levels"] = dict(
    tied={m: dict(mean=lv["tied_T2_T3"][m]["mean"], pi90=lv["tied_T2_T3"][m]["pi90"],
                  inside=lv["tied_T2_T3"][m]["inside"]) for m in ("M0", "M1", "M2")},
    C3={k: dict(mean=v["mean"], pi90=v["pi90"], inside=v["inside"], P_exceed_bound=v["P_exceed_bound"],
                ar_mean_pct=v["ar_mean_pct"]) for k, v in lv["C3"].items() if isinstance(v, dict)},
    C3_bound_pct=lv["C3"]["bound_pct"], C3_check=ver_c3,
    S1=lv["S1"], T7=lv["T7"], T5_PE=lv["T5_PE"], dispersion=lv["dispersion"], original_slope=lv["original_slope"],
    choir=dict(E_probable_included=lv["implied"]["H1"]["E_implied"], E_confirmed_only=lv["implied"]["H1"]["E_implied_alt"],
               E_coach=lv["implied"]["T2"]["E_implied"]),
    A4_P_exceed=dec["models"]["M2"]["A4"]["P_exceed"])
num["levels"]["choir"]["spread_probable"] = lv["implied"]["H1"]["E_implied"] / lv["implied"]["T2"]["E_implied"]
others = [lv["implied"][k]["E_implied"] for k in ("O2", "R2", "C5", "T2", "T3")] + [lv["implied"]["H1"]["E_implied_alt"]]
num["levels"]["choir"]["spread_confirmed_only"] = max(others) / min(others)


def pi(v):
    return f"{v[0]:.0f}--{v[1]:.0f}"


rows = ["% --- Table appC_protocol:tab:levels: data/bsc_validation/04_level_controls.json (tied_T2_T3, C3, S1, T5_PE, T7);",
        "%     last column of the classroom rows: data/validation_checks/v13_c3.txt"]
t = lv["tied_T2_T3"]
for i, (m, nm) in enumerate((("M0", "$\\Mzero$"), ("M1", "$\\Mone$"), ("M2", "$\\Mtwo$"))):
    rows.append(f"    {'A3: coach $\\to$ minibus, 17 exposed, 2 observed' if i == 0 else ''} & {nm} & "
                f"{t[m]['mean']:.2f} ({pi(t[m]['pi90'])}) & {'inside' if t[m]['inside'] else 'outside'} & \\\\")
rows.append("    \\addlinespace")
c3 = lv["C3"]
for i, (key, nm, f) in enumerate((("M2_masks_0.35", "$\\Mtwo$, factor 0.35", 0.35), ("M2_masks_0.2", "$\\Mtwo$, 0.2", 0.2),
                                  ("M2_masks_0.6", "$\\Mtwo$, 0.6", 0.6), ("M2_no_masks", "$\\Mtwo$, no masks", 1.0),
                                  ("M1_masks_0.35", "$\\Mone$, factor 0.35", None))):
    last = f"$\\Prob(X\\ge5)={ver_c3[f]['P_ge5']:.3f}$" if f is not None else ""
    rows.append(f"    {'A3: trains $\\to$ Utah classrooms, 728 contacts, 5 observed' if i == 0 else ''} & {nm} & "
                f"{c3[key]['mean']:.2f} ({pi(c3[key]['pi90'])}) & {'inside' if c3[key]['inside'] else 'outside'} & {last} \\\\")
rows.append("    \\addlinespace")
s1, pe, t7 = lv["S1"], lv["T5_PE"], lv["T7"]
rows.append(f"    A4: supermarket customers, 0/8{{,}}224 (bound {s1['bound_pct']:.4f}\\,\\%) & $\\Mzero$/$\\Mtwo$ & "
            f"{s1['M2_single_index']['median_pct']:.4f}\\,\\% (95th percentile {s1['M2_single_index']['q95_pct']:.4f}\\,\\%) & pass & "
            f"$\\Prob(\\text{{exceed}})={s1['M2_single_index']['P_exceed']:.4f}$ \\\\")
rows.append(f"    A4: masked Utah classrooms, 5/728 (bound {c3['bound_pct']:.2f}\\,\\%) & $\\Mzero$/$\\Mtwo$, 0.35 & "
            f"{c3['M2_masks_0.35']['ar_mean_pct']:.2f}\\,\\% & pass & $\\Prob(\\text{{exceed}})={c3['M2_masks_0.35']['P_exceed_bound']:.4f}$ \\\\")
rows.append(f"    A4: premium economy of the flight, 0/35 (bound {pe['bound_pct']:.1f}\\,\\%) & $\\Mtwo$ & "
            f"{pe['M2']['ar_median_pct']:.2f}\\,\\% (95th percentile {pe['M2']['ar_q95_pct']:.1f}\\,\\%) & pass & "
            f"$\\Prob(\\text{{exceed}})={pe['M2']['P_exceed']:.2f}$ \\\\")
rows.append(f"    A4: one 20-minute metro ride, 309 other riders & $\\Mzero$/$\\Mtwo$ & {t7['M2']['mean_cases']:.3f} expected cases & "
            f"pass & {int(round(t7['M2']['P_ge5'] * 20000))} of 20{{,}}000 rides with $\\ge5$ cases \\\\")
tex += rows

# ---------------------------------------------------------------------------------------------------------
# 6. call centre: first analysis (40 draws) against quadrature; software verification
# ---------------------------------------------------------------------------------------------------------
q = json.loads((VER / "callcentre_quadrature.json").read_text())
num["callcentre"] = dict(
    first_analysis_40_draws=dict(pooled=m2["O1"]["P_z_le_obs_pooled"], equal=m2["O1"]["P_z_le_obs_equal_weight"]),
    quadrature=dict(equal=q["equal_weight"], pooled=q["pooled_weight"]),
    prior_upper_limit=q["prior_upper_limit_sensitivity"], n_grid=len(q["grid_log10D"]),
    M1_first_analysis=dict(pooled=dec["models"]["M1"]["O1"]["P_z_le_obs_pooled"],
                           equal=dec["models"]["M1"]["O1"]["P_z_le_obs_equal_weight"]))
st = load("00_selftest.json")
ic = load("08_bruteforce_check.json")["events"]
num["software"] = dict(
    kernel_vs_quadrature=dict(M1=st["A_max_abs_err_M1"], M2=st["A_max_abs_err_M2"]),
    M2_Dinf_minus_M0=st["B_M2_Dinf_minus_M0"], mean_M2_minus_M0=st["B_mean_M2_minus_M0"],
    cond_pmf_vs_enumeration=st["C_cond_pmf_err"], hypergeom=st["C_hypergeom_err"],
    joincount_mean=st["D_mean_theory_sim"], joincount_var=st["D_var_theory_sim"], M1_pair=st["E_pair"],
    brute_force={e: dict(rel_err_exposure=v["rel_err_exposure"], max_abs_diff_pmf=v["max_abs_diff_pmf"],
                         mc_se_max=v["mc_se_max"], n_accepted=v["n_accepted"]) for e, v in ic.items()},
    brute_force_max_rel_err=max(v["rel_err_exposure"] for v in ic.values()),
    brute_force_pmf_diff_range=[min(v["max_abs_diff_pmf"] for v in ic.values()),
                                max(v["max_abs_diff_pmf"] for v in ic.values())])
num["decision_as_first_run"] = {k: dict(outcome_as_run=v["outcome_as_run"], outcome_conservative=v["outcome_conservative"],
                                        A2_failures_as_run=v["A2_failures_as_run"]) for k, v in dec["models"].items()}

num["_generated"] = dict(script="code/appC_protocol_tables.py",
                         at=dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
(OUT / "appC_protocol_numbers.json").write_text(json.dumps(num, indent=1))
(OUT / "appC_protocol_tables.tex").write_text(
    "% generated by code/appC_protocol_tables.py -- do not edit; rows are pasted into sections/supp_protocol1.tex\n"
    + "\n".join(tex) + "\n")

sec = HERE.parent / "sections" / "supp_protocol1.tex"
if sec.exists():
    body = sec.read_text()
    missing = [r for r in tex if not r.lstrip().startswith("%") and r.strip() != "\\addlinespace" and r.strip() not in body]
    print("generated rows not found verbatim in sections/supp_protocol1.tex:", len(missing))
    for r in missing:
        print("   MISSING:", r)
    if missing:
        sys.exit(1)
print("hashes match:", num["freeze"]["all_hashes_match"])
print("sensitivity: A1 passes in", summ["n_M2_A1_pass"], "of", summ["n_variants"], "; coach fails in", summ["n_coach_fail"])
print("vehicles common-D: LR %.2f p %.4f (first analysis); LR %.2f p %.4f (re-check)" % (
    num["loo"]["vehicle"]["LR"], num["loo"]["vehicle"]["p"], num["loo"]["vehicle_check"]["LR"],
    num["loo"]["vehicle_check"]["p"]))
print("\n".join(tex))

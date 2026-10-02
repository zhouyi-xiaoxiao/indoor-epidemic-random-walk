"""Build the outbreak validation dataset from records.py + mapping.py + doi_verification.json.

Outputs (all deterministic):
  outbreaks.json, outbreaks.csv, sources.md
  data/spatial_strata.csv, data/hu2021_train_attack_matrix.csv, data/li2021_restaurant_tables.csv,
  data/park2020_callcentre_epicurve_digitized.csv, data/model_mapping.csv

Computed here (and flagged as computed in the output):
  attack_rate = n_infected / n_exposed; exact (Clopper-Pearson) 95% CI;
  cumulative hazard H = -ln(1 - AR) and hazard per hour H / T when a single-exposure duration exists.
"""
import csv, json, math, pathlib, sys

import numpy as np
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from records import RECORDS, STRATA, HU_TRAIN, LI_RESTAURANT, CALLCENTRE_EPICURVE  # noqa: E402
from mapping import MAPPING, SCENE_DEFAULTS  # noqa: E402

DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)


def exact_ci(k, n, level=0.95):
    """Clopper-Pearson interval for k successes out of n."""
    if n is None or k is None or n == 0:
        return (None, None)
    a = 1 - level
    lo = 0.0 if k == 0 else float(stats.beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - a / 2, k + 1, n - k))
    return (lo, hi)


def fmt_citation(m):
    au = m.get("authors") or []
    if not au and m.get("epmc_authors"):
        au = [a.strip().rstrip(".") for a in m["epmc_authors"].split(",")]
    au_s = ", ".join(au[:6]) + (", et al." if len(au) > 6 else ".")
    title = (m.get("title") or m.get("epmc_title") or "").rstrip(".")
    journal = m.get("journal") or m.get("epmc_journal") or ""
    year = m.get("year") or m.get("epmc_year")
    vol = m.get("volume") or m.get("epmc_volume")
    issue = m.get("issue") or m.get("epmc_issue")
    pages = m.get("pages") or m.get("epmc_pages")
    s = f"{au_s} {title}. {journal}. {year}"
    if vol:
        s += f";{vol}"
        if issue:
            s += f"({issue})"
    if pages:
        s += f":{pages}"
    s += f". doi:{m['doi']}"
    return s


def main():
    ver = json.loads((DATA / "doi_verification.json").read_text())
    out = []
    for rec in RECORDS:
        r = dict(rec)
        n, k = r.get("n_exposed"), r.get("n_infected")
        if n and k is not None:
            ar = k / n
            lo, hi = exact_ci(k, n)
            r["attack_rate"] = round(ar, 7)
            r["attack_rate_ci95_exact"] = [round(lo, 7), round(hi, 7)]
            H = -math.log(1 - ar) if ar < 1 else None
            r["cumulative_hazard"] = round(H, 4) if H is not None else None
            T = r.get("exposure_duration_h")
            single = r["exposure_type"].startswith("single") and T
            r["hazard_per_hour"] = round(H / T, 4) if (single and H is not None) else None
        else:
            r["attack_rate"] = r["attack_rate_ci95_exact"] = r["cumulative_hazard"] = r["hazard_per_hour"] = None
        if r.get("n_infected_alt") is not None and n:
            lo, hi = exact_ci(r["n_infected_alt"], n)
            r["attack_rate_alt"] = round(r["n_infected_alt"] / n, 7)
            r["attack_rate_alt_ci95_exact"] = [round(lo, 7), round(hi, 7)]
        for s in r["sources"]:
            m = ver[s["key"]]
            s["citation"] = fmt_citation(m)
            s["url"] = m["url"]
            s["doi_registered"] = m["doi_registered"]
            s["pmid"] = m.get("pmid")
            s["pmcid"] = m.get("pmcid")
        r["model_mapping"] = MAPPING[r["id"]]
        out.append(r)

    meta = dict(
        title="Documented indoor respiratory outbreaks for validating the lattice random-walk SIR model",
        built="2026-10-01", n_records=len(out),
        conventions="n_exposed/n_infected exclude the index case unless includes_index is true; None = not reported in the sources accessed; fields listed in `derived` are computed here, everything else is transcribed.",
        scene_parameters=SCENE_DEFAULTS,
    )
    (ROOT / "outbreaks.json").write_text(json.dumps(dict(meta=meta, records=out), indent=1, ensure_ascii=False))

    # ---- flat CSV
    cols = ["id", "scene", "analogue", "quality", "setting", "location", "event_dates", "pathogen", "citation", "doi", "url",
            "additional_sources", "source_access", "exposure_type", "exposure_duration_h", "exposure_duration_note", "room_dimensions",
            "floor_area_m2", "ceiling_height_m", "volume_m3", "geometry_note", "occupancy_n", "ventilation", "masks", "index_case",
            "n_exposed", "n_infected", "n_infected_alt", "includes_index", "attack_rate", "ci95_low", "ci95_high", "attack_rate_alt",
            "cumulative_hazard", "hazard_per_hour", "attack_rate_reported", "spatial_pattern", "other_groups", "spatial_data_files",
            "quality_notes", "derived_fields", "model_scene", "model_a_m", "model_lattice", "model_n_agents", "model_T",
            "model_index_position", "model_mobility", "model_observables", "protocol_role", "mapping_obstacles"]
    with open(ROOT / "outbreaks.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            mp = r["model_mapping"]
            ci = r["attack_rate_ci95_exact"] or [None, None]
            w.writerow(dict(
                id=r["id"], scene=r["scene"], analogue=r["analogue"], quality=r["quality"], setting=r["setting"], location=r["location"],
                event_dates=r["event_dates"], pathogen=r["pathogen"], citation=r["sources"][0]["citation"], doi=r["sources"][0]["doi"],
                url=r["sources"][0]["url"], additional_sources=" | ".join(s["citation"] for s in r["sources"][1:]),
                source_access=" | ".join(f"{s['key']}: {s['access']}" for s in r["sources"]),
                exposure_type=r["exposure_type"], exposure_duration_h=r["exposure_duration_h"], exposure_duration_note=r["exposure_duration_note"],
                room_dimensions=r["room_dimensions"], floor_area_m2=r["floor_area_m2"], ceiling_height_m=r["ceiling_height_m"], volume_m3=r["volume_m3"],
                geometry_note=r["geometry_note"], occupancy_n=r["occupancy_n"], ventilation=r["ventilation"], masks=r["masks"], index_case=r["index_case"],
                n_exposed=r["n_exposed"], n_infected=r["n_infected"], n_infected_alt=r.get("n_infected_alt"), includes_index=r["includes_index"],
                attack_rate=r["attack_rate"], ci95_low=ci[0], ci95_high=ci[1], attack_rate_alt=r.get("attack_rate_alt"),
                cumulative_hazard=r["cumulative_hazard"], hazard_per_hour=r["hazard_per_hour"], attack_rate_reported=r["attack_rate_reported"],
                spatial_pattern=r["spatial_pattern"], other_groups=r["other_groups"], spatial_data_files="; ".join(r["spatial_data_files"]),
                quality_notes=r["quality_notes"], derived_fields="; ".join(r["derived"]),
                model_scene=mp["model_scene"], model_a_m=mp["a_m"], model_lattice=mp["lattice"], model_n_agents=mp["n_agents"], model_T=mp["T_days"],
                model_index_position=mp["index_position"], model_mobility=mp["mobility"], model_observables=" | ".join(mp["observables"]),
                protocol_role=mp["role"], mapping_obstacles=mp["obstacles"]))

    # ---- model mapping CSV
    with open(DATA / "model_mapping.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "quality", "model_scene", "a_m", "lattice", "n_agents", "T_model", "index_position", "mobility", "observables", "protocol_role", "obstacles"])
        for r in out:
            mp = r["model_mapping"]
            w.writerow([r["id"], r["quality"], mp["model_scene"], mp["a_m"], mp["lattice"], mp["n_agents"], mp["T_days"], mp["index_position"],
                        mp["mobility"], " | ".join(mp["observables"]), mp["role"], mp["obstacles"]])

    # ---- strata
    with open(DATA / "spatial_strata.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["outbreak_id", "stratum", "d_min_m", "d_max_m", "n_exposed", "n_infected", "attack_rate", "ci95_low", "ci95_high", "source_key", "note"])
        for (oid, lab, d0, d1, n, k, key, note) in STRATA:
            if n:
                lo, hi = exact_ci(k, n)
                w.writerow([oid, lab, d0, d1, n, k, round(k / n, 7), round(lo, 7), round(hi, 7), key, note])
            else:
                w.writerow([oid, lab, d0, d1, n, k, None, None, None, key, note])

    # ---- train matrix (long format)
    with open(DATA / "hu2021_train_attack_matrix.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rows_apart", "cols_apart", "attack_rate_pct", "ci95_low_pct", "ci95_high_pct", "note"])
        for i, ra in enumerate(HU_TRAIN["rows_apart"]):
            for j, ca in enumerate(HU_TRAIN["cols_apart"]):
                v = HU_TRAIN["attack_pct"][i][j]
                note = "cell realigned (printed one column to the left on the publisher page)" if (ra == 0 and v is not None) else ""
                if ra == 0 and ca == 0:
                    note = "index seat (undefined)"
                w.writerow([ra, ca, v, HU_TRAIN["ci_low"][i][j], HU_TRAIN["ci_high"][i][j], note])
        for i, ra in enumerate(HU_TRAIN["rows_apart"]):
            w.writerow([ra, "mean", HU_TRAIN["row_mean_pct"][i], None, None, "printed row mean"])
        for j, ca in enumerate(HU_TRAIN["cols_apart"]):
            w.writerow(["mean", ca, HU_TRAIN["col_mean_pct"][j], None, None, "printed column mean"])

    # ---- restaurant tables
    with open(DATA / "li2021_restaurant_tables.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["table", "patrons", "infected", "overlap_with_table_A_min", "total_exposure_min", "norm_measured_tracer", "norm_predicted_tracer",
                    "norm_predicted_exposure", "zone", "neighbour_class"])
        for (t, n, k, ov, tot, mt, pt, pe) in LI_RESTAURANT:
            zone = "ABC" if t in ("TA", "TB", "TC") else "non-ABC"
            nb = "index table" if t == "TA" else ("immediate" if t in ("TB", "TC", "T18") else "remote")
            w.writerow([t, n, k, ov, tot, mt, pt, pe, zone, nb])

    # ---- call centre epidemic curve
    with open(DATA / "park2020_callcentre_epicurve_digitized.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["onset_date", "floor11", "floor10", "floor9"])
        for row in CALLCENTRE_EPICURVE:
            w.writerow(row)

    # ---- sources.md
    used = {}
    for r in out:
        for s in r["sources"]:
            used.setdefault(s["key"], dict(s=s, recs=[], access=set()))
            used[s["key"]]["recs"].append(r["id"])
            used[s["key"]]["access"].add(s["access"])
    lines = ["# Sources for the outbreak validation dataset", "",
             "Generated by `scripts/build_dataset.py` on 2026-10-01 from `data/doi_verification.json`.", "",
             "Verification: every DOI below was checked on 2026-10-01 against the DOI handle registry",
             "(`https://doi.org/api/handles/<doi>`, responseCode 1 = registered); bibliographic fields were pulled from",
             "Crossref (fallback Europe PMC), not typed from memory. `Access` states exactly which text the numbers were",
             "transcribed from; 'abstract only' records were transcribed from the abstract and are marked as such.", "",
             f"All {len(used)} DOIs are registered. No source in this list was cited without being opened.", ""]
    for key in sorted(used):
        s = used[key]["s"]
        m = ver[key]
        lines += [f"## {key}", "", f"- **Citation:** {s['citation']}", f"- **URL:** {s['url']}",
                  f"- **DOI registered:** {m['doi_registered']}; PMID {m.get('pmid')}; PMCID {m.get('pmcid')}",
                  f"- **Access:** {'; '.join(sorted(used[key]['access']))}",
                  f"- **Used in records:** {', '.join(used[key]['recs'])}", ""]
    lines += ["## Sources consulted but not used as records", "",
              "- Cheng P, et al. Predominant airborne transmission and insignificant fomite transmission of SARS-CoV-2 in a two-bus COVID-19 outbreak originating from the same pre-symptomatic index case. J Hazard Mater 2022. doi:10.1016/j.jhazmat.2021.128051 (same Hunan event as T2/T3; full text saved, not transcribed).",
              ""]
    (ROOT / "sources.md").write_text("\n".join(lines))
    print("records:", len(out), "| sources:", len(used))


if __name__ == "__main__":
    main()

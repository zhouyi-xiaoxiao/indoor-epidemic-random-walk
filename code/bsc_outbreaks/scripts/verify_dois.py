"""Verify that every DOI used in records.py exists, and pull bibliographic metadata.

For each DOI:
  1. https://doi.org/api/handles/<doi>   -> responseCode 1 means the DOI is registered
     (and gives the resolver target URL);
  2. https://api.crossref.org/works/<doi> -> title, authors, journal, year, volume, issue, pages;
  3. Europe PMC search (fallback for metadata, also gives PMID/PMCID).
Output: data/doi_verification.json (one entry per source key).  Network access required;
results are cached in the JSON so downstream scripts are deterministic.
"""
import json, pathlib, sys, time, urllib.parse, urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from records import RECORDS  # noqa: E402

UA = {"User-Agent": "outbreak-dataset-verification/1.0 (mailto:research@example.org)"}


def get_json(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def main():
    keys = {}
    for rec in RECORDS:
        for s in rec["sources"]:
            keys[s["key"]] = s["doi"]
    out = {}
    for key, doi in sorted(keys.items()):
        entry = dict(doi=doi, url=f"https://doi.org/{doi}", checked="2026-10-01")
        try:
            h = get_json("https://doi.org/api/handles/" + urllib.parse.quote(doi))
            entry["doi_registered"] = h.get("responseCode") == 1
            vals = [v for v in h.get("values", []) if v.get("type") == "URL"]
            entry["resolves_to"] = vals[0]["data"]["value"] if vals else None
        except Exception as e:  # noqa: BLE001
            entry["doi_registered"] = False
            entry["handle_error"] = repr(e)
        try:
            m = get_json("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
            entry["title"] = (m.get("title") or [None])[0]
            entry["journal"] = (m.get("container-title") or [None])[0]
            entry["volume"], entry["issue"], entry["pages"] = m.get("volume"), m.get("issue"), m.get("page")
            dp = (m.get("published-print") or m.get("published") or m.get("issued") or {}).get("date-parts", [[None]])
            entry["year"] = dp[0][0]
            entry["authors"] = [f"{a.get('family', '')} {''.join(w[0] for w in a.get('given', '').replace('-', ' ').split())}".strip()
                                for a in m.get("author", [])]
            entry["crossref_ok"] = True
        except Exception as e:  # noqa: BLE001
            entry["crossref_ok"] = False
            entry["crossref_error"] = repr(e)
        try:
            q = urllib.parse.quote(f'DOI:"{doi}"')
            r = get_json(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=lite")
            res = r["resultList"]["result"]
            if res:
                entry["pmid"], entry["pmcid"] = res[0].get("pmid"), res[0].get("pmcid")
                entry["epmc_title"] = res[0].get("title")
                entry["epmc_authors"] = res[0].get("authorString")
                entry["epmc_journal"] = res[0].get("journalTitle")
                entry["epmc_year"] = res[0].get("pubYear")
                entry["epmc_volume"], entry["epmc_issue"], entry["epmc_pages"] = res[0].get("journalVolume"), res[0].get("issue"), res[0].get("pageInfo")
        except Exception as e:  # noqa: BLE001
            entry["epmc_error"] = repr(e)
        out[key] = entry
        print(key, doi, "registered=", entry.get("doi_registered"), "crossref=", entry.get("crossref_ok"), "pmid=", entry.get("pmid"))
        time.sleep(0.4)
    (ROOT / "data" / "doi_verification.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    bad = [k for k, v in out.items() if not v.get("doi_registered")]
    print("unregistered:", bad)


if __name__ == "__main__":
    main()

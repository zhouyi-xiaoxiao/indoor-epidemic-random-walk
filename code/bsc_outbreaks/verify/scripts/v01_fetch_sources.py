"""Second fetch of every cited source (re-check's own code).

For each DOI: (1) doi.org handle registry, (2) Crossref metadata, (3) Europe PMC search by DOI
(-> PMID/PMCID, abstract, open-access flag), (4) Europe PMC full-text XML when available.
Outputs: verify/results/source_check.json and verify/src/<key>.txt (plain text of the JATS XML),
verify/src/<key>.abstract.txt. Network step; not deterministic in the sense that remote servers
may change, but all downstream scripts read the saved copies.
"""
import json, re, sys, time, urllib.request, urllib.parse, html, os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"; SRC.mkdir(exist_ok=True)
RES = ROOT / "results"; RES.mkdir(exist_ok=True)

DOIS = {
 "cai2020": "10.3201/eid2606.200412",
 "chen2020tianjin": "10.1016/j.ajic.2020.06.006",
 "furuse2020": "10.3201/eid2609.202272",
 "guenther2020": "10.15252/emmm.202013296",
 "hamner2020": "10.15585/mmwr.mm6919e6",
 "hershow2021": "10.15585/mmwr.mm7012e3",
 "hu2021": "10.1093/cid/ciaa1057",
 "jang2020": "10.3201/eid2608.200633",
 "jiang2021deptstore": "10.1007/s11783-021-1386-6",
 "katelaris2021": "10.3201/eid2706.210465",
 "khanh2020": "10.3201/eid2611.203299",
 "koh2020": "10.1136/oemed-2020-106626",
 "kwon2020": "10.3346/jkms.2020.35.e415",
 "lamhine2021": "10.15585/mmwr.mm7035e2",
 "lan2021": "10.1136/oemed-2020-106774",
 "li2021": "10.1016/j.buildenv.2021.107788",
 "li2022twoclusters": "10.3390/ijerph19084876",
 "lu2020": "10.3201/eid2607.200764",
 "luo2020": "10.1093/ofid/ofaa430",
 "miller2021": "10.1111/ina.12751",
 "mizumoto2020": "10.2807/1560-7917.ES.2020.25.10.2000180",
 "moser1979": "10.1093/oxfordjournals.aje.a112781",
 "nardell1991": "10.1164/ajrccm/144.2.302",
 "nicholls2024": "10.1093/occmed/kqad100",
 "ou2022": "10.1016/j.buildenv.2021.108414",
 "park2020": "10.3201/eid2608.201274",
 "riley1978": "10.1093/oxfordjournals.aje.a112560",
 "shen2020": "10.1001/jamainternmed.2020.5225",
 "steinzamir2020": "10.2807/1560-7917.ES.2020.25.29.2001352",
 "tian2021": "10.1186/s12889-021-10713-z",
 "weissberg2020": "10.1186/s13756-020-00861-z",
}
UA = {"User-Agent": "outbreak-dataset-check/1.0 (mailto:zhouyixiaoxiao@gmail.com)"}

def get(url, timeout=40, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except Exception as e:  # noqa
            last = e; time.sleep(1.5 * (i + 1))
    return None, repr(last).encode()

def jats_to_text(xml: str) -> str:
    xml = re.sub(r"<xref[^>]*>.*?</xref>", "", xml, flags=re.S)
    xml = re.sub(r"</(p|title|sec|tr|caption|label|abstract|td|th)>", "\n", xml)
    txt = re.sub(r"<[^>]+>", " ", xml)
    txt = html.unescape(txt)
    txt = re.sub(r"[ \t]+", " ", txt)
    txt = re.sub(r"\n\s*\n+", "\n", txt)
    return txt

out = {}
only = sys.argv[1:] or list(DOIS)
prev = {}
if (RES / "source_check.json").exists():
    prev = json.load(open(RES / "source_check.json"))
for key in DOIS:
    if key not in only:
        if key in prev: out[key] = prev[key]
        continue
    doi = DOIS[key]; rec = {"doi": doi}
    st, body = get("https://doi.org/api/handles/" + urllib.parse.quote(doi, safe="/"))
    try:
        j = json.loads(body); rec["handle_responseCode"] = j.get("responseCode")
        urls = [v["data"]["value"] for v in j.get("values", []) if v.get("type") == "URL"]
        rec["handle_url"] = urls[0] if urls else None
    except Exception:
        rec["handle_responseCode"] = None; rec["handle_error"] = body[:200].decode("utf8", "replace")
    st, body = get("https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/"))
    try:
        m = json.loads(body)["message"]
        rec["crossref"] = dict(
            title=(m.get("title") or [None])[0],
            first_author=(m.get("author") or [{}])[0].get("family"),
            n_authors=len(m.get("author") or []),
            journal=(m.get("container-title") or [None])[0],
            year=(m.get("issued", {}).get("date-parts") or [[None]])[0][0],
            volume=m.get("volume"), issue=m.get("issue"), page=m.get("page") or m.get("article-number"),
        )
    except Exception:
        rec["crossref"] = None; rec["crossref_error"] = body[:200].decode("utf8", "replace")
    q = urllib.parse.quote(f'DOI:"{doi}"')
    st, body = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&resultType=core&format=json")
    try:
        res = json.loads(body)["resultList"]["result"]
        r0 = res[0] if res else {}
        rec["epmc"] = dict(pmid=r0.get("pmid"), pmcid=r0.get("pmcid"), title=r0.get("title"),
                           isOpenAccess=r0.get("isOpenAccess"), inEPMC=r0.get("inEPMC"),
                           journal=(r0.get("journalInfo") or {}).get("journal", {}).get("title"),
                           year=r0.get("pubYear"), authorString=(r0.get("authorString") or "")[:120])
        if r0.get("abstractText"):
            (SRC / f"{key}.abstract.txt").write_text(html.unescape(re.sub(r"<[^>]+>", " ", r0["abstractText"])))
            rec["abstract_saved"] = True
    except Exception:
        rec["epmc"] = None
    pmcid = (rec.get("epmc") or {}).get("pmcid")
    rec["fulltext_chars"] = 0
    if pmcid:
        st, body = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML", timeout=60)
        if st == 200 and len(body) > 2000:
            (SRC / f"{key}.xml").write_bytes(body)
            txt = jats_to_text(body.decode("utf8", "replace"))
            (SRC / f"{key}.txt").write_text(txt)
            rec["fulltext_chars"] = len(txt)
        else:
            rec["fulltext_status"] = st
    out[key] = rec
    print(key, rec.get("handle_responseCode"), (rec.get("crossref") or {}).get("title", "")[:60] if rec.get("crossref") else None,
          pmcid, rec["fulltext_chars"], flush=True)
    json.dump(out, open(RES / "source_check.json", "w"), indent=1, ensure_ascii=False)

"""v2_refs.py -- references of the second outbreak test and of the further analyses (needs network access).

Fetches the Crossref record of every reference of Sections 7.5-7.9 and Supplementary Sections S6, S9 and S10 in refs.bib,
stores it in ../data/v2_refs_crossref.json and prints BibTeX entries built from the records (authors, title,
journal, year, volume, issue, pages or article number, DOI), so that no bibliographic field is typed by hand.

    python v2_refs.py            # fetch and write ../data/v2_refs_crossref.json, print the entries
"""
import json
import os
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "v2_refs_crossref.json")
UA = {"User-Agent": "indoor-epidemic-random-walk/1.0 (mailto:zhouyixiaoxiao@gmail.com)"}

DOIS = {
    # eleven further outbreaks
    "olsen2003transmission": ("10.1056/NEJMoa031349", "flight CA112, SARS-CoV-1 (F1)"),
    "kenyon1996transmission": ("10.1056/NEJM199604113341501", "Chicago-Honolulu flight, tuberculosis (F2)"),
    "baker2010transmission": ("10.1136/bmj.c2424", "Los Angeles-Auckland flight, influenza (F3)"),
    "young2014international": ("10.1111/irv.12181", "Cancun-Birmingham flight, influenza (F4)"),
    "toyokawa2022transmission": ("10.1111/irv.12913", "flight to Naha, SARS-CoV-2 (F5)"),
    "speake2020flight": ("10.3201/eid2612.203910", "Sydney-Perth flight, SARS-CoV-2 (F6)"),
    "swadi2021genomic": ("10.3201/eid2703.204714", "flight EK448, SARS-CoV-2 (F7)"),
    "hoehl2020assessment": ("10.1001/jamanetworkopen.2020.18044", "Tel Aviv-Frankfurt flight, SARS-CoV-2 (F8)"),
    "wong2004cluster": ("10.3201/eid1002.030452", "ward 8A, medical students (W1)"),
    "yu2005temporal": ("10.1086/428735", "ward 8A, inpatients (W2)"),
    "hertzberg2016risk": ("10.1016/j.aogh.2016.06.003", "two-row rule on aircraft; redrawn seat map of CA112"),
    # contact records
    "genois2018can": ("10.1140/epjds/s13688-018-0140-1", "office 2015 and conference contact data"),
    "mastrandrea2015contact": ("10.1371/journal.pone.0136497", "high-school contact data"),
    "vanhems2013estimating": ("10.1371/journal.pone.0073970", "hospital-ward contact data"),
    "stehle2011simulation": ("10.1186/1741-7015-9-87", "SEIR on a recorded conference contact network"),
    "cattuto2010dynamics": ("10.1371/journal.pone.0011596", "RFID proximity sensing; broad contact-duration distributions"),
    "starnini2013modeling": ("10.1103/PhysRevLett.110.168701", "random-walk model of face-to-face interactions with attractiveness (anchored walkers)"),
    # tracer measurements
    "kinahan2021aerosol": ("10.1371/journal.pone.0246916", "aerosol tracer tests in wide-body cabins"),
    "woodward2022evaluation": ("10.1111/ina.13121", "salt-aerosol tracer in a rail carriage"),
    # zone level
    "endo2021within": ("10.1073/pnas.2112605118", "influenza in 29 primary schools, Matsumoto"),
    "gemmetto2014mitigation": ("10.1186/s12879-014-0695-9", "class and school mixing from the Lyon contact data"),
    "vink2014serial": ("10.1093/aje/kwu209", "serial interval of influenza"),
}


# Fields that Crossref does not settle and that were taken from the publisher's page (1 October 2026).
MANUAL = {
    "hertzberg2016risk": "Annals of Global Health 82(5), September-October 2016, pp. 819-823; Crossref lists the "
                         "online date (2017) and the first page only. Authors as printed: V. S. Hertzberg, H. Weiss.",
    "woodward2022evaluation": "article number e13121 (in the DOI; Crossref has no page field).",
    "baker2010transmission": "article number c2424; the Crossref issue field ('may21 1') is not used.",
}


def fetch(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["message"]


def bibtex(key, m, citedfor):
    au = []
    for a in m.get("author", []):
        if "family" in a:
            au.append(f"{a['family']}, {a.get('given', '')}".strip().rstrip(","))
        elif "name" in a:
            au.append("{" + a["name"] + "}")
    if len(au) > 10:
        au = au[:10] + ["others"]
    year = (m.get("published-print") or m.get("published") or m.get("issued"))["date-parts"][0][0]
    pages = m.get("page") or m.get("article-number") or ""
    fields = [("author", " and ".join(au)), ("title", m["title"][0]),
              ("journal", (m.get("container-title") or [""])[0]), ("year", str(year)),
              ("volume", m.get("volume", "")), ("number", m.get("issue", "")), ("pages", pages.replace("-", "--")),
              ("doi", m["DOI"]), ("citedfor", citedfor)]
    body = ",\n".join(f"  {k:<9} = {{{v}}}" for k, v in fields if v)
    return f"@article{{{key},\n{body}\n}}"


if __name__ == "__main__":
    rec, out = {}, []
    for key, (doi, citedfor) in DOIS.items():
        m = fetch(doi)
        rec[key] = {k: m.get(k) for k in ("DOI", "title", "author", "container-title", "volume", "issue", "page",
                                          "article-number", "published", "published-print", "issued", "type")}
        rec[key]["citedfor"] = citedfor
        out.append(bibtex(key, m, citedfor))
        time.sleep(0.3)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump({"fetched": time.strftime("%Y-%m-%d"), "source": "api.crossref.org/works/<doi>", "records": rec,
                   "manual_checks": MANUAL},
                  f, indent=1, ensure_ascii=False)
    print("\n\n".join(out))

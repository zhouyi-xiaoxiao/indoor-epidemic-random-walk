"""refs_search.py -- Crossref records of selected references, and the targeted literature
search behind the statements of what is added in Sections 1, 3 and 5 (needs network access).

1. Crossref records of eight cited references and of further entries whose bibliographic data (title,
   authors, journal, volume, issue, pages, year) were checked against Crossref -> "crossref".
2. Targeted search: the Crossref and arXiv queries listed in QUERIES, with the number of records returned and
   the first 20 titles of each -> "searches".  This is a targeted search, not a systematic review.
3. "manual_checks": what was read by hand on 1 October 2026 for the points that Crossref does not settle
   (page ranges missing from Crossref, theorem numbers of the published version of Gao & Dong, abstracts).
   These entries are a record typed by the author; the script only stores them.

    python refs_search.py
Output: ../data/refs_search.json
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "refs_search.json")
UA = {"User-Agent": "indoor-epidemic-random-walk/1.0 (mailto:zhouyixiaoxiao@gmail.com)"}

DOIS = {
    # references whose Crossref records are stored
    "kesten2005spread": "10.1214/009117905000000413",
    "kesten2006phase": "10.1215/ijm/1258059486",
    "kesten2008shape": "10.4007/annals.2008.167.701",
    "maia2007diffusive": "10.1088/0953-8984/19/6/065143",
    "tien2015disease": "10.1007/s00285-014-0791-x",
    "bichara2015sis": "10.1007/s11538-015-0113-5",
    "peng2022practical": "10.1021/acs.est.1c06531",
    "carvalho2025generalized": "10.1088/1742-5468/ad9f4d",     # found by the targeted search below
    # bibliographic data checked against Crossref
    "cohen1981convexity": "10.1090/s0002-9939-1981-0601750-2",
    "shen2020": "10.1001/jamainternmed.2020.5225",
    "mizumoto2020": "10.2807/1560-7917.ES.2020.25.10.2000180",
    "steinzamir2020": "10.2807/1560-7917.ES.2020.25.29.2001352",
    # attribution re-checked
    "gao2020fast": "10.1090/proc/14868",
    "allen2007asymptotic": "10.1137/060672522",
    "gao2011sis": "10.1016/j.mbs.2011.05.001",
}

QUERIES = [
    ("crossref", "spread of a rumor or infection in a moving population random walks"),
    ("crossref", "diffusive epidemic process"),
    ("crossref", "infection independent random walkers lattice phase transition recovery"),
    ("crossref", "epidemic random walkers lattice basic reproduction number"),
    ("crossref", "two random walkers transmission of infection confinement"),
    ("crossref", "SIS patch model basic reproduction number monotone dispersal rate"),
    ("crossref", "basic reproduction number fast dispersal limit patch model expansion"),
    ("crossref", "residence time patch epidemic model virtual dispersal"),
    ("crossref", "indoor airborne transmission outbreaks compilation risk indicators"),
    ("arxiv", 'all:"diffusive epidemic process"'),
    ("arxiv", 'ti:infection AND ti:"random walks"'),
    ("arxiv", 'ti:epidemic AND abs:"random walkers" AND abs:lattice'),
]

MANUAL_CHECKS = {
    "date": "2026-10-01",
    "kesten2005spread": "Pages 2402-2462: journal reference of arXiv:math/0312496 (arXiv API). Crossref has no page field.",
    "kesten2006phase": "Pages 547-634: Project Euclid page of the DOI. Model as stated in the abstract of "
                       "arXiv:math/0410371: independent continuous-time random walks on Z^d, a healthy particle "
                       "becomes infected when it meets an infected one, infected particles recuperate at a rate "
                       "lambda; there is a critical recuperation rate below which the infection survives with "
                       "positive probability and above which it dies out.",
    "maia2007diffusive": "Abstract (OpenAlex): continuous absorbing-state phase transition of the one-dimensional "
                         "diffusive epidemic process studied by mean-field theory and Monte Carlo simulation; two "
                         "species hop on a lattice, with recovery and infection on contact; particle number conserved.",
    "carvalho2025generalized": "Abstract (arXiv:2407.08175): walkers on a two-dimensional lattice, infection "
                               "between individuals on the same node, recovery with permanent immunity; a phase "
                               "transition in the dynamic-percolation universality class controlled by the "
                               "population size. Found by the arXiv query on the diffusive epidemic process; "
                               "journal data from Crossref.",
    "tien2015disease": "Cited for the 1/d expansion of R0 at fast dispersal on the authority of the Discussion of "
                       "Gao & Dong (2020, published version): Tien et al. derived the limit of R0 as the dispersal "
                       "rate tends to infinity and a difference of order 1/d_I through a Laurent series. The paper "
                       "of Tien et al. itself was not read (title, authors, journal, pages from Crossref).",
    "bichara2015sis": "Abstract (arXiv:1503.08881): multi-group epidemic framework in which the risk of infection "
                      "is a function of residence times in patches; R0 computed as a function of the "
                      "residence-times matrix.",
    "peng2022practical": "Abstract (OpenAlex): risk parameters combining emission rate, breathing rate, masking, "
                         "ventilation and removal rates, occupancy and duration, applied to documented COVID-19 "
                         "outbreaks and to measles, influenza and tuberculosis outbreaks.",
    "gao2020fast": "Published version (Proc. Amer. Math. Soc. 148 (2020) 1709-1722, PDF on ams.org) read: "
                   "Theorem 2.1 (symmetric connectivity), Theorem 2.3 (R0 strictly decreasing and strictly convex "
                   "in the dispersal rate, irreducible connectivity), Theorem 3.3 (min patch number < R0(infinity) "
                   "< R0(d) < max patch number; its proof credits the min/max statement to Gao & Ruan 2011), "
                   "Corollary 3.5 (spectral bound s(F - V), 'a generalization of Lemma 3.4 in Allen et al.'), "
                   "Introduction (Allen et al.: global stability of the disease-free equilibrium for R0 < 1, "
                   "unique endemic equilibrium for R0 > 1, conjecture that R0 decreases with the dispersal rate; "
                   "Gao 2019 proved it for symmetric connectivity), assumption (B2) (Proposition 2.2 of Gao & "
                   "Ruan: R0 independent of movement when the patch numbers coincide). The numbering is the same "
                   "as in arXiv:1907.12229v1.",
    "allen2007asymptotic": "Not read (no access). Statements attributed to it in the article follow the account "
                           "in Gao & Dong (2020).",
    "cohen1981convexity": "Pages 657-658 on the publisher's page (ams.org); Crossref gives 657-657.",
    "shen2020": "Pages 1665-1671 (PubMed 32870239); Crossref gives the first page only.",
    "mizumoto2020": "Article number 2000180 (PubMed 32183930).",
    "steinzamir2020": "Article number 2001352 (PubMed 32720636).",
    "moser1979": "PubMed 463858 (abstract): a jet airliner with 54 persons aboard, delayed three hours; within "
                 "72 hours 72 per cent of the passengers became ill.",
    "miller2021": "Accepted manuscript, Table 1: probability of infection sampled uniformly on 53-87 %, volumetric "
                  "breathing rate uniform on 0.65-1.38 m3/h, loss rates 0.3-1.0 (ventilation), 0.3-1.5 "
                  "(deposition), 0-0.63 (inactivation) per hour; result 970 +- 390 quanta per hour; sensitivity "
                  "runs at 1.0 m3/h.",
}


def get(url, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:                      # noqa: BLE001
            err = str(e)
            time.sleep(2 + 2 * k)
    return json.dumps({"error": err})


def crossref_record(doi):
    m = json.loads(get("https://api.crossref.org/works/" + urllib.parse.quote(doi))).get("message", {})
    return dict(doi=doi, title=(m.get("title") or [""])[0],
                authors=[f"{a.get('family', '')}, {a.get('given', '')}" for a in m.get("author", [])],
                journal=(m.get("container-title") or [""])[0], volume=m.get("volume"), issue=m.get("issue"),
                pages=m.get("page"), article_number=m.get("article-number"),
                year=(m.get("issued", {}).get("date-parts") or [[None]])[0][0], type=m.get("type"))


def crossref_search(q, rows=20):
    d = json.loads(get("https://api.crossref.org/works?rows=%d&query.bibliographic=%s"
                       % (rows, urllib.parse.quote(q)))).get("message", {})
    items = [dict(title=(i.get("title") or [""])[0], doi=i.get("DOI"),
                  year=(i.get("issued", {}).get("date-parts") or [[None]])[0][0],
                  journal=(i.get("container-title") or [""])[0]) for i in d.get("items", [])]
    return dict(engine="crossref", query=q, total_results=d.get("total-results"), first=items)


def arxiv_search(q, rows=20):
    t = get("https://export.arxiv.org/api/query?max_results=%d&search_query=%s" % (rows, urllib.parse.quote(q)))
    tot = re.search(r"<opensearch:totalResults[^>]*>(\d+)<", t)
    items = [dict(id=re.search(r"<id>(.*?)</id>", e).group(1),
                  title=" ".join(re.search(r"<title>(.*?)</title>", e, re.S).group(1).split()))
             for e in re.findall(r"<entry>(.*?)</entry>", t, re.S)]
    return dict(engine="arxiv", query=q, total_results=int(tot.group(1)) if tot else None, first=items)


def main():
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    out["retrieved"] = time.strftime("%Y-%m-%d")
    out.setdefault("crossref", {})
    for key, doi in DOIS.items():
        if key not in out["crossref"] or not out["crossref"][key].get("title"):
            out["crossref"][key] = crossref_record(doi)
            json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
            time.sleep(1)
    out.setdefault("searches", [])
    done = {(s["engine"], s["query"]) for s in out["searches"] if s.get("total_results") is not None}
    for engine, q in QUERIES:
        if (engine, q) in done:
            continue
        out["searches"] = [s for s in out["searches"] if (s["engine"], s["query"]) != (engine, q)]
        out["searches"].append(crossref_search(q) if engine == "crossref" else arxiv_search(q))
        json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
        time.sleep(3)
    out["manual_checks"] = MANUAL_CHECKS
    json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
    print(json.dumps({k: v["title"] for k, v in out["crossref"].items()}, indent=1, ensure_ascii=False))
    for s in out["searches"]:
        print(s["engine"], "|", s["query"], "|", s["total_results"], "|", len(s["first"]))


if __name__ == "__main__":
    main()

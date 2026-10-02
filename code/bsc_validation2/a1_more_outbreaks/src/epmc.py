"""Europe PMC helper: search, fetch abstract + OA full-text XML into data/raw/."""
import sys, json, urllib.request, urllib.parse, os, re, time
RAW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "raw")
def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research; outbreak data)"})
            return urllib.request.urlopen(req, timeout=40).read()
        except Exception as e:
            err = e; time.sleep(2)
    raise err
def search(q, n=5):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&pageSize=%d&query=%s" % (n, urllib.parse.quote(q))
    return json.loads(get(u))["resultList"]["result"]
def fetch(key, q):
    res = search(q, 3)
    if not res: print(key, "NO RESULT"); return
    r = res[0]
    meta = {k: r.get(k) for k in ("id","source","pmid","pmcid","doi","title","authorString","journalTitle","pubYear","journalVolume","issue","pageInfo","isOpenAccess","abstractText")}
    meta["journalTitle"] = (r.get("journalInfo") or {}).get("journal", {}).get("title")
    meta["volume"] = (r.get("journalInfo") or {}).get("volume"); meta["issue"]=(r.get("journalInfo") or {}).get("issue")
    meta["query"] = q; meta["accessed"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    json.dump(meta, open(os.path.join(RAW, key + "_meta.json"), "w"), indent=1, ensure_ascii=False)
    got = "no-fulltext"
    if r.get("pmcid"):
        try:
            x = get("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML" % r["pmcid"])
            open(os.path.join(RAW, key + "_fulltext.xml"), "wb").write(x); got = "xml %d bytes" % len(x)
        except Exception as e:
            got = "fulltext fail %s" % e
    print(key, "|", r.get("pmid"), r.get("pmcid"), r.get("doi"), "|", (r.get("title") or "")[:90], "|", r.get("pubYear"), "|", got)
if __name__ == "__main__":
    for line in open(sys.argv[1]):
        line=line.strip()
        if not line or line.startswith("#"): continue
        key, q = line.split("\t", 1)
        if os.path.exists(os.path.join(RAW, key + "_meta.json")) and "--force" not in sys.argv: print(key, "cached"); continue
        try: fetch(key, q)
        except Exception as e: print(key, "ERROR", e)

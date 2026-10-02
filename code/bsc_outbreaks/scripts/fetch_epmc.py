"""Fetch Europe PMC metadata/full text for a DOI or PMCID and save plain text.

Usage: python fetch_epmc.py <DOI-or-PMCID> [outname]
Writes notes/fulltext/<outname>.txt (plain text of the JATS XML body, tables included)
and prints metadata. Used only to *read* sources; nothing here is modelling.
"""
import json, re, sys, urllib.request, urllib.parse, pathlib, xml.etree.ElementTree as ET

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
OUT = pathlib.Path(__file__).resolve().parents[1] / "notes" / "fulltext"

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "research-script/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()

def meta(q):
    if q.upper().startswith("PMC"):
        query = f"PMCID:{q}"
    else:
        query = f'DOI:"{q}"'
    url = f"{BASE}/search?query={urllib.parse.quote(query)}&format=json&resultType=core"
    d = json.loads(get(url))
    return d["resultList"]["result"]

def xml_to_text(xml_bytes):
    root = ET.fromstring(xml_bytes)
    parts = []
    def walk(el, depth=0):
        tag = el.tag.split('}')[-1]
        if tag in ("title",):
            parts.append("\n## " + "".join(el.itertext()).strip() + "\n")
            return
        if tag in ("p", "caption"):
            parts.append("".join(el.itertext()).strip() + "\n")
            return
        if tag == "tr":
            parts.append(" | ".join("".join(c.itertext()).strip() for c in el) + "\n")
            return
        for c in el:
            walk(c, depth + 1)
    walk(root)
    return "\n".join(parts)

if __name__ == "__main__":
    q = sys.argv[1]
    name = sys.argv[2] if len(sys.argv) > 2 else re.sub(r"[^A-Za-z0-9]+", "_", q)
    res = meta(q)
    if not res:
        print("NO RESULT for", q); sys.exit(1)
    r = res[0]
    info = {k: r.get(k) for k in ("title", "doi", "pmid", "pmcid", "journalTitle", "pubYear", "authorString", "isOpenAccess", "inEPMC")}
    info["journalTitle"] = r.get("journalInfo", {}).get("journal", {}).get("title")
    print(json.dumps(info, ensure_ascii=False, indent=1))
    OUT.mkdir(parents=True, exist_ok=True)
    txt = "# " + json.dumps(info, ensure_ascii=False) + "\n\nABSTRACT: " + (r.get("abstractText") or "") + "\n\n"
    pmcid = r.get("pmcid")
    if pmcid:
        try:
            x = get(f"{BASE}/{pmcid}/fullTextXML")
            txt += xml_to_text(x)
            print("fulltext chars:", len(txt))
        except Exception as e:
            print("fulltext failed:", e)
    (OUT / f"{name}.txt").write_text(txt)
    print("saved", OUT / f"{name}.txt")

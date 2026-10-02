"""Fallback full-text fetch via NCBI BioC PMC API (open-access subset). Usage: fetch_bioc.py PMCID outname"""
import json, sys, urllib.request, pathlib
OUT = pathlib.Path(__file__).resolve().parents[1] / "notes" / "fulltext"
pmcid, name = sys.argv[1], sys.argv[2]
url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmcid}/unicode"
req = urllib.request.Request(url, headers={"User-Agent": "research-script/1.0"})
d = json.loads(urllib.request.urlopen(req, timeout=60).read())
docs = d if isinstance(d, list) else [d]
lines = []
for coll in docs:
    for doc in coll["documents"]:
        for p in doc["passages"]:
            t = p.get("infons", {}).get("type", "")
            txt = p.get("text", "")
            if t.startswith("title"):
                lines.append("\n## " + txt)
            else:
                lines.append(txt)
path = OUT / f"{name}.txt"
prev = path.read_text() if path.exists() else ""
path.write_text(prev + "\n\n# BIOC FULLTEXT\n" + "\n".join(lines))
print(name, "chars", len("\n".join(lines)))

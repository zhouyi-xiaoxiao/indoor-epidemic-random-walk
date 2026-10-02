import sys, re, os, xml.etree.ElementTree as ET
RAW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "raw")
def txt(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()
def conv(key):
    p = os.path.join(RAW, key + "_fulltext.xml")
    root = ET.parse(p).getroot()
    out = []
    for el in root.iter():
        tag = el.tag.split('}')[-1]
        if tag in ("article-title",) and not out: out.append("# " + txt(el))
        if tag == "abstract": out.append("## ABSTRACT\n" + txt(el))
        if tag == "title": out.append("\n### " + txt(el))
        if tag == "p": out.append(txt(el))
        if tag == "table-wrap":
            out.append("\n[TABLE] " + " ".join(txt(c) for c in el if c.tag.split('}')[-1] in ("label","caption")))
            for tr in el.iter():
                if tr.tag.split('}')[-1] == "tr":
                    out.append(" | ".join(txt(c) for c in tr))
        if tag == "fig":
            g = [x.attrib.get("{http://www.w3.org/1999/xlink}href") for x in el.iter() if x.tag.split('}')[-1]=="graphic"]
            out.append("\n[FIG] " + txt(el) + " GRAPHIC=" + str(g))
        if tag == "supplementary-material":
            out.append("[SUPP] " + txt(el) + " " + str([ (x.attrib) for x in el.iter() if x.tag.split('}')[-1]=="media"]))
    # dedupe consecutive
    res = []
    for o in out:
        if not res or res[-1] != o: res.append(o)
    open(os.path.join(RAW, key + "_fulltext.txt"), "w").write("\n".join(res))
    return len(res)
if __name__ == "__main__":
    keys = sys.argv[1:] or [f[:-13] for f in os.listdir(RAW) if f.endswith("_fulltext.xml")]
    for k in keys: print(k, conv(k))

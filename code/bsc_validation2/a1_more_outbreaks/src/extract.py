import sys,re,os
RAW=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","raw")
OUT=os.path.join(os.path.dirname(RAW),"..","notes")
pat=re.compile(sys.argv[1],re.I)
for k in sys.argv[2:]:
    lines=open(os.path.join(RAW,k+"_fulltext.txt")).read().split("\n")
    seen=set();out=[]
    for i,l in enumerate(lines):
        if i<9: continue
        if (pat.search(l) or l.startswith("[TABLE]") or " | " in l) and l not in seen:
            seen.add(l); out.append("%d: %s"%(i,l[:2200]))
    open(os.path.join(OUT,"extract_%s.txt"%k),"w").write("\n".join(out))
    print(k,len(out))

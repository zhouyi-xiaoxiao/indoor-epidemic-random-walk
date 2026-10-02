import sys,re,html,os
RAW=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","raw")
for k in sys.argv[1:]:
    s=open(os.path.join(RAW,k+"_pmc.html"),errors="ignore").read()
    m=re.search(r"<main.*?</main>",s,re.S); s=m.group(0) if m else s
    s=re.sub(r"<(script|style).*?</\1>","",s,flags=re.S)
    s=re.sub(r"</(td|th)>"," | ",s); s=re.sub(r"</(p|tr|h\d|div|li|figcaption)>","\n",s)
    s=re.sub(r'<img[^>]*src="([^"]+)"[^>]*>',r" [IMG \1] ",s)
    s=html.unescape(re.sub(r"<[^>]+>","",s)); s=re.sub(r"[ \t]+"," ",s); s=re.sub(r"\n\s*\n+","\n",s)
    open(os.path.join(RAW,k+"_fulltext.txt"),"w").write(s); print(k,len(s))

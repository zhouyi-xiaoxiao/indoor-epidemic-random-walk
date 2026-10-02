import re,html,sys
for p in sys.argv[1:]:
    x=open(p).read()
    x=re.sub(r'<xref[^>]*>.*?</xref>','',x,flags=re.S)
    x=re.sub(r'</(p|title|tr|caption|sec|label)>','\n',x); x=re.sub(r'</t[dh]>',' | ',x); x=re.sub(r'<[^>]+>','',x); x=html.unescape(x)
    open(p.replace('.xml','.txt'),'w').write(x); print(p,len(x))

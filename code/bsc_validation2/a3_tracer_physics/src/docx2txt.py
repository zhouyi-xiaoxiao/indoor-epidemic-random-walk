import zipfile,re,sys,html
def docx_text(p):
    z=zipfile.ZipFile(p); x=z.read('word/document.xml').decode('utf8')
    x=re.sub(r'</w:tc>','\t',x); x=re.sub(r'</w:tr>','\n',x); x=re.sub(r'</w:p>','\n',x)
    x=re.sub(r'<[^>]+>','',x); return html.unescape(x), [n for n in z.namelist() if 'media' in n]
if __name__=='__main__':
    for p in sys.argv[1:]:
        t,m=docx_text(p); open(p+'.txt','w').write(t); print('=====',p,len(t),m)

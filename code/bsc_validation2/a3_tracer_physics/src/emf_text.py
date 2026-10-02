"""Minimal EMF text extractor with world-transform tracking (positions in device units)."""
import struct,sys
def mul(a,b):  # apply a then b ; xform = (m11,m12,m21,m22,dx,dy): x'=x*m11+y*m21+dx ; y'=x*m12+y*m22+dy
    a11,a12,a21,a22,adx,ady=a; b11,b12,b21,b22,bdx,bdy=b
    return (a11*b11+a12*b21, a11*b12+a12*b22, a21*b11+a22*b21, a21*b12+a22*b22,
            adx*b11+ady*b21+bdx, adx*b12+ady*b22+bdy)
def texts(p):
    b=open(p,'rb').read(); out=[]; off=0; W=(1,0,0,1,0,0); stack=[]
    while off+8<=len(b):
        t,sz=struct.unpack_from('<II',b,off)
        if sz<8: break
        if t==33: stack.append(W)
        elif t==34:
            if stack: W=stack.pop()
        elif t==35: W=struct.unpack_from('<6f',b,off+8)
        elif t==36:
            X=struct.unpack_from('<6f',b,off+8); mode=struct.unpack_from('<I',b,off+32)[0]
            if mode==1: W=(1,0,0,1,0,0)
            elif mode==2: W=mul(X,W)
            elif mode==3: W=mul(W,X)
            elif mode==4: W=X
        elif t in (84,83):
            x,y,n,offs=struct.unpack_from('<iiII',b,off+36)
            s=b[off+offs:off+offs+2*n].decode('utf-16le','replace') if t==84 else b[off+offs:off+offs+n].decode('latin1')
            X=x*W[0]+y*W[2]+W[4]; Y=x*W[1]+y*W[3]+W[5]
            out.append((round(X,1),round(Y,1),round(W[0],2),round(W[1],2),s))
        off+=sz
    return out
if __name__=='__main__':
    for p in sys.argv[1:]:
        print('=====',p)
        for r in texts(p): print(r)

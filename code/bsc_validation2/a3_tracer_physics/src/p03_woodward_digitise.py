"""P3: digitise Woodward et al. 2022 Fig. 9 (steady-state PM2.5 along a rail carriage).
Markers by colour segmentation (red: y<0 side, blue: y>=0 side), axes calibrated on the plot frame
(x/L from -0.5 to 0.5, ordinate 0 to 3.5).  Also the solid 1D-model curve.
Output data/derived/woodward2022_fig9_points.csv, woodward2022_fig9_curve.csv"""
import os
import numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage as ndi
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
im = np.asarray(Image.open(os.path.join(ROOT, 'data', 'raw', 'woodward2022_supp', 'INA-32-0-g001.jpg')).convert('RGB')).astype(int)
H, W, _ = im.shape
R, G, B = im[..., 0], im[..., 1], im[..., 2]
grey = (R + G + B) < 450          # frame lines are grey
dark = ((R + G + B) < 400) & (abs(R - B) < 40)   # curve and text
print('image', W, H)
# frames: rows / columns with long dark runs
rowsum = grey.sum(axis=1); colsum = grey.sum(axis=0)
hr = np.flatnonzero(rowsum > 0.6 * W); vc = np.flatnonzero(colsum > 0.3 * H)
def groups(a):
    g = np.split(a, np.flatnonzero(np.diff(a) > 2) + 1); return [float(x.mean()) for x in g if len(x)]
hl, vl = groups(hr), groups(vc)
print('horizontal frame lines', hl); print('vertical frame lines', vl)
assert len(hl) == 4 and len(vl) == 2
x0, x1 = vl
panels = {'A_end': (hl[0], hl[1], -0.392), 'B_middle': (hl[2], hl[3], -0.135)}
red = (R > 170) & (G < 110) & (B < 110)
blue = (B > 170) & (R < 110) & (G < 110)
pts, cur = [], []
for name, (ytop, ybot, xsrc) in panels.items():
    fx = lambda px: -0.5 + (px - x0) / (x1 - x0)
    fy = lambda py: 3.5 * (ybot - py) / (ybot - ytop)
    for colour, mask in (('red', red), ('blue', blue)):
        m = mask.copy(); m[:int(ytop) + 2] = False; m[int(ybot) - 1:] = False
        m[:, :int(x0) + 2] = False
        # legend region: x/L > 0.2 and ordinate > 1.0 contains the legend symbols -> excluded
        lab, n = ndi.label(ndi.binary_dilation(m, iterations=1))
        for i in range(1, n + 1):
            yy, xx = np.nonzero(lab == i)
            if len(yy) < 6: continue
            cx, cy = fx(xx.mean()), fy(yy.mean())
            if cx > 0.2 and cy > 1.0: continue
            h = fy(yy.min()) - fy(yy.max())
            pts.append(dict(panel=name, colour=colour, x_over_L=round(cx, 4), c_norm=round(cy, 3),
                            c_lo=round(fy(yy.max()), 3), c_hi=round(fy(yy.min()), 3), npix=len(yy),
                            near_source=abs(cx - xsrc) < 0.02))
    # 1D-model curve: dark pixels strictly inside the frame, away from the dotted source line and legend text
    for px in range(int(x0) + 3, int(x1) - 2, 3):
        col = dark[int(ytop) + 3:int(ybot) - 3, px]
        ys = np.flatnonzero(col)
        if len(ys) == 0: continue
        g = groups(ys)
        # the curve is the lowest short dark run in the column (text/legend sits above 1.0 on the right)
        cand = [v for v in g]
        yv = [fy(v + int(ytop) + 3) for v in cand]
        xv = fx(px)
        ok = [v for v in yv if (xv < 0.18 or v < 1.0)]
        if ok and abs(xv - xsrc) > 0.012:
            cur.append(dict(panel=name, x_over_L=round(xv, 4), c_norm=round(min(ok), 3), n_runs=len(g)))
P = pd.DataFrame(pts).sort_values(['panel', 'x_over_L']); C = pd.DataFrame(cur)
P.to_csv(os.path.join(ROOT, 'data', 'derived', 'woodward2022_fig9_points.csv'), index=False)
C.to_csv(os.path.join(ROOT, 'data', 'derived', 'woodward2022_fig9_curve.csv'), index=False)
pd.set_option('display.width', 200); pd.set_option('display.max_rows', 300)
print(P.to_string()); print(C.groupby('panel').apply(lambda d: d.iloc[::12][['x_over_L', 'c_norm']].to_string()).to_string())

"""P1: seat-level tracer values of Ou et al. 2022 Fig. S3 / S4, read from the EMF text records.
Output: data/derived/ou2022_B1_seats.csv, ou2022_B2_seats.csv
 columns: seat, row, col, x_fig, y_fig, s3_measured, s3_cfd (test conditions), s4_cfd (infection conditions)"""
import os, re, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from emf_text import texts
MED = os.path.join(ROOT, 'data', 'raw', 'ou2022_supp', 'media')

def parse(fn):
    rec = texts(os.path.join(MED, fn))
    lab = [(x, y, s.strip()) for x, y, a, b, s in rec if (a, b) == (0.0, 1.0)]
    val = [(x, y, s.strip()) for x, y, a, b, s in rec if (a, b) == (1.0, 0.0)]
    # join split seat labels (same x, next piece 74 units further along)
    L = []
    i = 0
    while i < len(lab):
        x, y, s = lab[i]
        if i + 1 < len(lab) and abs(lab[i + 1][0] - x) < 2 and 60 < lab[i + 1][1] - y < 90 and re.fullmatch(r'\d+', s):
            s = s + lab[i + 1][2]; i += 1
        L.append((x, y, s)); i += 1
    # join split numbers (same y, consecutive records)
    V = []
    i = 0
    while i < len(val):
        x, y, s = val[i]
        j = i + 1
        while j < len(val) and abs(val[j][1] - y) < 2 and 0 < val[j][0] - x < 260 and re.fullmatch(r'[\d.]+', val[j][2]) and re.fullmatch(r'[\d.]+', s) and not re.fullmatch(r'\d\.\d\d', s):
            s += val[j][2]; j += 1
        V.append((x, y, s)); i = j
    return L, V

def assign(L, V, paired):
    seats = [(x, y, s) for x, y, s in L if re.fullmatch(r'\d+[A-E]|[CD]', s)]
    nums = [(x, y, float(re.search(r'(\d\.\d\d)', s).group(1)), s) for x, y, s in V if re.search(r'\d\.\d\d', s)]
    out = {}
    used = set()
    if paired:   # measured above predicted: same x, dy ~ 160-185
        pairs = []
        for i, (x, y, v, s) in enumerate(nums):
            for j, (x2, y2, v2, s2) in enumerate(nums):
                if j != i and abs(x2 - x) < 2 and 150 < y2 - y < 195 and i not in used and j not in used:
                    pairs.append((x, 0.5 * (y + y2), v, v2)); used |= {i, j}
        singles = [(x, y, None, v) for i, (x, y, v, s) in enumerate(nums) if i not in used]
        items = pairs + singles
    else:
        items = [(x, y, None, v) for x, y, v, s in nums]
    return seats, items

def match(seats, items, off):
    rows = []
    sx = np.array([[x, y] for x, y, s in seats])
    taken = {}
    for x, y, m, p in items:
        d = np.hypot(sx[:, 0] + off[0] - x, sx[:, 1] + off[1] - y)
        k = int(np.argmin(d))
        rows.append(dict(seat=seats[k][2], x_fig=seats[k][0], y_fig=seats[k][1], measured=m, cfd=p, dist=float(d[k])))
    return pd.DataFrame(rows)

def run(bus, f3, f4, off3, off4):
    L, V = parse(f3); seats, items = assign(L, V, True); a = match(seats, items, off3)
    L4, V4 = parse(f4); seats4, items4 = assign(L4, V4, False); b = match(seats4, items4, off4)
    print(bus, 'S3 items', len(a), 'max dist', a.dist.max().round(1), '| S4 items', len(b), 'max dist', b.dist.max().round(1))
    for nm, df in (('S3', a), ('S4', b)):
        dup = df[df.seat.duplicated(keep=False)]
        if len(dup): print('  DUPLICATE seat assignment in', nm, '\n', dup)
    allseats = pd.DataFrame([dict(seat=s, x_fig=x, y_fig=y) for x, y, s in seats4])
    df = allseats.merge(a[['seat', 'measured', 'cfd']].rename(columns={'measured': 's3_measured', 'cfd': 's3_cfd'}), on='seat', how='left')
    df = df.merge(b[['seat', 'cfd']].rename(columns={'cfd': 's4_cfd'}), on='seat', how='left')
    df['row'] = df.seat.str.extract(r'(\d+)')[0].astype(float)
    df['col'] = df.seat.str.extract(r'([A-E])$')[0]
    df = df.sort_values(['row', 'col']).reset_index(drop=True)
    df.to_csv(os.path.join(ROOT, 'data', 'derived', f'ou2022_{bus}_seats.csv'), index=False)
    return df, a, b

if __name__ == '__main__':
    pd.set_option('display.width', 200); pd.set_option('display.max_rows', 200)
    d1, a1, b1 = run('B1', 'image7.emf', 'image9.emf', (180, 160), (180, 195))
    print(d1.to_string())
    # B2 (minibus): the automatic pairing fails where a CFD-only value of an A seat sits directly above the
    # measured/CFD pair of the B seat.  The 30+17 text records were assigned by hand from
    # data/raw/ou2022_supp/emf_text_dump.txt (label position + vertical offset: measured ~ +50..+97,
    # CFD ~ +230..+280) and checked against the rendered figure (conv/mmc1.pdf, page S13-S14).
    B2 = {  # seat: (s3_measured, s3_cfd, s4_cfd)
        '1A': (None, 0.56, 0.61), '1B': (0.36, 0.77, 0.78), '2A': (None, 0.60, 0.56), '2B': (0.46, 0.72, 0.65),
        '3A': (None, 0.69, 0.59), '3B': (0.69, 0.74, 0.68), '4A': (None, 0.96, 0.67), '4A_mother': (None, None, 1.10),
        '4B': (1.00, 1.00, 1.00), '5A': (0.61, 0.65, 0.68), '5B': (0.73, 0.63, None), '5C': (0.84, 0.75, 0.62),
        '5D': (0.74, 0.86, 0.57), '5E': (None, 0.65, 0.52), '4D': (0.56, 0.59, None), '3C': (0.66, 0.72, 0.63),
        '2C': (1.00, 0.87, 0.60), 'D': (0.26, 0.86, 0.88), 'C': (None, 0.77, 0.91)}
    d2 = pd.DataFrame([dict(seat=k, s3_measured=v[0], s3_cfd=v[1], s4_cfd=v[2]) for k, v in B2.items()])
    d2.to_csv(os.path.join(ROOT, 'data', 'derived', 'ou2022_B2_seats.csv'), index=False)
    print(d2.to_string())

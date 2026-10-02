"""Digitize the seat-level infection map of the Seoul call-centre outbreak.

Source figure: Park SY et al., Emerg Infect Dis 2020;26(8):1666-1670, Figure 2
("Floor plan of the 11th floor of building X ... Blue indicates the seating places of
persons with confirmed cases."), https://wwwnc.cdc.gov/eid/images/20-1274-F2.jpg
(downloaded to notes/figs_src/20-1274-F2.jpg).

Method (deterministic, no randomness):
  * pixel classes: blue fill (case seat) and white fill (non-case seat);
  * a desk is a connected component with bounding box ~31x17 px (vertical desk) or
    ~17x36 px (end-of-column horizontal desk) and fill fraction > 0.8;
  * one chair-shaped component (~13x9 px) accompanies every blue desk, which is used as a
    consistency check (#blue chairs == #blue desks);
  * the east-wing blue desk (17x31 px, a different furniture type) is added explicitly by
    its colour (it is the only blue component of desk size outside the two wings);
  * hatched wall strips on the west side produce desk-sized white components and are
    excluded by position (x < 40 and 300 < y < 640).
Outputs: data/park2020_callcentre_seats_digitized.csv, figures/callcentre_digitization_check.{png,pdf}

Caveats: the figure has no scale bar; physical coordinates are given only
in units of the desk pitch. Seat occupancy on the exposure days is not reported; we assume
one employee per drawn desk. The figure marks 84 case seats whereas the text reports 94 cases
on the floor, so 10 cases have no seat in the figure.
"""
import pathlib
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
IMG = ROOT / "notes" / "figs_src" / "20-1274-F2.jpg"


def components(mask):
    lab, n = ndi.label(mask)
    out = []
    for i, o in enumerate(ndi.find_objects(lab)):
        s = int((lab[o] == i + 1).sum())
        h = o[0].stop - o[0].start
        w = o[1].stop - o[1].start
        out.append(dict(y0=o[0].start, y1=o[0].stop, x0=o[1].start, x1=o[1].stop,
                        h=h, w=w, area=s, fill=s / (h * w)))
    return out


def is_desk(c):
    return (((28 <= c["h"] <= 33 and 15 <= c["w"] <= 20) or
             (15 <= c["h"] <= 20 and 30 <= c["w"] <= 39)) and c["fill"] > 0.8)


def is_chair(c):
    return (11 <= c["h"] <= 14 and 7 <= c["w"] <= 10) or (7 <= c["h"] <= 10 and 11 <= c["w"] <= 14)


def region(xc, yc):
    if yc < 300:
        return "north_wing"
    if yc > 640 and xc < 700:
        return "south_wing"
    if xc > 700:
        return "east_offices"
    return "other"


def main():
    im = np.asarray(Image.open(IMG).convert("RGB")).astype(int)
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    blue = (B > 200) & (R < 190) & (G > 150) & (G < 225) & (B - R > 40)
    white = (R > 215) & (G > 215) & (B > 215)
    cb, cw = components(blue), components(white)
    rows = []
    for colour, comps in (("blue", cb), ("white", cw)):
        for c in comps:
            if not is_desk(c):
                continue
            xc, yc = (c["x0"] + c["x1"]) / 2, (c["y0"] + c["y1"]) / 2
            reg = region(xc, yc)
            if colour == "white" and xc < 40 and 300 < yc < 640:
                continue  # hatched wall strip, not a desk
            if colour == "white" and reg not in ("north_wing", "south_wing"):
                continue  # other furniture types are not digitized for non-cases
            rows.append(dict(x_px=xc, y_px=yc, w_px=c["w"], h_px=c["h"], case=int(colour == "blue"),
                             region=reg))
    df = pd.DataFrame(rows).sort_values(["region", "y_px", "x_px"]).reset_index(drop=True)
    n_blue_chairs = sum(is_chair(c) for c in cb)
    # desk pitch (vertical desks): median spacing of y centres in the north wing
    nw = df[(df.region == "north_wing") & (df.h_px > 25)]
    ys = np.sort(nw.y_px.unique())
    pitch = float(np.median(np.diff(ys)[np.diff(ys) > 5]))
    df["x_pitch_units"] = df.x_px / pitch
    df["y_pitch_units"] = df.y_px / pitch
    out = ROOT / "data" / "park2020_callcentre_seats_digitized.csv"
    df.to_csv(out, index=False)
    summ = df.groupby("region").case.agg(["sum", "count"])
    summ["attack_rate_by_seat"] = summ["sum"] / summ["count"]
    print(summ)
    print("total blue desks", int(df.case.sum()), "blue chairs", n_blue_chairs, "desk pitch px", pitch)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(im.astype(np.uint8))
    for _, r in df.iterrows():
        ax.plot(r.x_px, r.y_px, "x" if r.case else "o", ms=4, mfc="none",
                color="red" if r.case else "black", mew=0.8)
    ax.set_title("Digitization check: Park et al. 2020 Fig. 2\nred x = case seat (blue), black o = non-case desk",
                 fontsize=9)
    ax.set_axis_off()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"callcentre_digitization_check.{ext}", dpi=200)
    return df


if __name__ == "__main__":
    main()

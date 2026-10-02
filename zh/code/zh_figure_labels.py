"""zh_figure_labels.py -- redraw one Chinese-label figure with the terminology of the Chinese edition.

The Chinese edition uses the *_zh.pdf variants written by the figure scripts of the English article
(../figures).  One of the five used in the main text carries a label that differs from the term used in the Chinese text:

    fig_hotspot_maps_zh      panel (b): 感染效率 q      -> 传染效率 q          (infection efficiency q)

This script imports the drawing function of the English figure (../code/plan_fig_theory.py), replaces only this
Chinese string, and writes the PDF into ../figures of the Chinese edition.  Nothing is recomputed: the function
draws stored results, so data and numbers are those of the English figure.  The English article's figure
directory is not written to.

    python code/zh_figure_labels.py        (run from the zh directory or anywhere)
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ART_CODE = HERE.parents[1] / "code"
OUT = HERE.parent / "figures"
sys.path.insert(0, str(ART_CODE))

import figstyle  # noqa: E402  (sets the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

REPLACE = {
    r"感染效率 $q$": r"传染效率 $q$",
}


def patched_tr(orig):
    def tr(lang, en, zh):
        for a, b in REPLACE.items():
            zh = zh.replace(a, b)
        return orig(lang, en, zh)
    return tr


def main():
    import plan_fig_theory as pft
    pft.tr = patched_tr(pft.tr)
    OUT.mkdir(exist_ok=True)
    for name, draw in (("fig_hotspot_maps", pft.fig_hotspot_maps),):
        figstyle.setup("zh")
        fig = draw("zh")
        path = OUT / f"{name}_zh.pdf"
        fig.savefig(str(path), bbox_inches="tight")
        plt.close(fig)
        print(path.name)


if __name__ == "__main__":
    main()

"""Bilingual figure helpers: every figure is produced twice, with English
labels (suffix _en) and Simplified-Chinese labels (suffix _zh), as vector PDF
and PNG."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

FIG_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures"))

_ZH_CANDIDATES = ["Songti SC", "STSong", "PingFang SC", "Heiti SC", "STHeiti",
                  "Arial Unicode MS"]

# colour-blind-safe palette (Okabe-Ito)
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
         purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", black="#000000",
         grey="#7F7F7F")


def _zh_font():
    names = {f.name for f in font_manager.fontManager.ttflist}
    for c in _ZH_CANDIDATES:
        if c in names:
            return c
    return None


def setup(lang: str):
    plt.rcdefaults()
    rc = {
        "figure.dpi": 100, "savefig.dpi": 200, "font.size": 9,
        "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.unicode_minus": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "legend.frameon": False, "lines.linewidth": 1.4,
        "mathtext.fontset": "dejavusans",
    }
    if lang == "zh":
        f = _zh_font()
        if f is None:
            raise RuntimeError("no CJK font found")
        rc["font.family"] = [f, "DejaVu Sans"]
        rc["font.sans-serif"] = [f, "DejaVu Sans"]
    else:
        rc["font.family"] = ["DejaVu Sans"]
    plt.rcParams.update(rc)


def tr(lang: str, en: str, zh: str) -> str:
    return en if lang == "en" else zh


def save(fig, name: str, lang: str):
    os.makedirs(FIG_DIR, exist_ok=True)
    base = os.path.join(FIG_DIR, f"{name}_{lang}")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    return base


def bilingual(draw, name: str):
    """Call ``draw(lang)`` (which must return a Figure) for 'en' and 'zh'."""
    out = []
    for lang in ("en", "zh"):
        setup(lang)
        fig = draw(lang)
        out.append(save(fig, name, lang))
    return out

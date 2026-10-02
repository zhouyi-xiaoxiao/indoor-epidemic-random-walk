"""figstyle.py -- shared bilingual figure helper for every figure generated inside article/code.

Same visual style as the research directories (code/bsc_sim/src/bsc_sim/plotting.py and
code/bsc_validation/src/valmod/plotting.py): DejaVu Sans, Okabe-Ito colours, vector PDF + PNG.

Every figure script must call ``bilingual(draw, name)``; ``draw(lang)`` returns a Figure and uses
``tr(lang, english, chinese)`` for every string, so that each figure exists as <name>_en.pdf and
<name>_zh.pdf (the Chinese edition of the article uses the _zh files).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
_ZH_CANDIDATES = ["Songti SC", "STSong", "PingFang SC", "Heiti SC", "STHeiti", "Arial Unicode MS"]

# colour-blind-safe palette (Okabe-Ito)
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00", purple="#CC79A7",
         sky="#56B4E9", yellow="#F0E442", black="#000000", grey="#7F7F7F")


def _zh_font():
    names = {f.name for f in font_manager.fontManager.ttflist}
    for c in _ZH_CANDIDATES:
        if c in names:
            return c
    return None


def setup(lang: str):
    plt.rcdefaults()
    rc = {
        "figure.dpi": 100, "savefig.dpi": 200, "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
        "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "axes.unicode_minus": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "legend.frameon": False, "lines.linewidth": 1.4, "mathtext.fontset": "dejavusans",
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
    FIG_DIR.mkdir(exist_ok=True)
    base = FIG_DIR / f"{name}_{lang}"
    fig.savefig(str(base) + ".pdf", bbox_inches="tight")
    fig.savefig(str(base) + ".png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    return str(base)


def bilingual(draw, name: str):
    """Call ``draw(lang)`` (which must return a Figure) for 'en' and 'zh' and save both."""
    out = []
    for lang in ("en", "zh"):
        setup(lang)
        out.append(save(draw(lang), name, lang))
    return out

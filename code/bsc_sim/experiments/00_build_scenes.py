#!/usr/bin/env python
"""Build the four scenes as data (JSON) from ASCII floor plans.

Each scene is a rectangular lattice. One character = one lattice cell of side
``a`` metres. Zone letters:

  office      W workstations, C corridor, M meeting room, K kitchen, B barrier
  supermarket S shelf lanes, A main aisles, F fresh produce, Q checkout queue,
              E entrance/exit, B shelf racks / counters (barrier)
  classroom   d desks (window block), D desks (inner block), P platform,
              I aisles, O obstacle
  metro       T seats, N standing area, G door vestibule, R grab rails

Zone values of the relative diffusion coefficient D/D0 and of the infection
efficiency q, the floor plans and the numbers of people are illustrative
assumptions; they are not measurements and are not derived from a standard.

Run:  python experiments/00_build_scenes.py
Writes scenes/<name>.json and prints zone counts.
"""
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "scenes")

# ----------------------------------------------------------------------------
# Office: 30 m x 20 m, a = 1.5 m, 20 x 13 = 260 cells, N = 100.
# Area fractions: W 60 %, C 15 %, M 10 %, K 5 %, B 10 %
#   -> 156 / 39 / 26 / 13 / 26 cells (reproduced exactly).
OFFICE = [
    "BWWWWWWWWCWWWWBMMMMM",
    "WWWWWWWWWCWWWWBMMMMM",
    "WWWWBWWWWCWWBWBMMMMM",
    "WWWWWWWWWCWWWWBMMMMM",
    "WWWWWWWBWCWWWWBMMMMM",
    "WWWWWWWWWCWWWWBBMBBB",
    "CCCCCCCCCCCCCCCCCCCC",
    "WWWWWWWWWCWWWWCCWWWW",
    "WWWWWWWBWCWWWWCWWWBW",
    "WWWWWWWWWCWWWWCBBBBB",
    "WWWWBWWWWCWWBWCBKKKK",
    "WWWWWWWWWCWWWWCKKKKK",
    "BWWWWWWWWCWWWWCBKKKK",
]
OFFICE_ZONES = {
    "W": {"name_en": "workstations", "name_zh": "工位区", "D": 0.3, "q": 1.4},
    "C": {"name_en": "corridor", "name_zh": "走廊", "D": 1.5, "q": 0.8},
    "M": {"name_en": "meeting room", "name_zh": "会议室", "D": 0.5, "q": 1.5},
    "K": {"name_en": "kitchen", "name_zh": "茶水间", "D": 1.0, "q": 1.2},
    "B": {"name_en": "barrier", "name_zh": "障碍物", "D": 0.0, "q": 0.0,
          "barrier": True},
}

# ----------------------------------------------------------------------------
# Supermarket: 50 m x 30 m, a = 2 m, 25 x 15 = 375 cells, N = 100-200 (150).
# Zones S, A, F, Q, E; the floor plan is a design choice.
# Shelf racks and checkout counters are modelled as barriers (B).
MARKET = [
    "FFFFFFFFFFFFAAAAAAAAAAAAA",
    "FFFFFFFFFFFFAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "ASBSBSBSBSBSASBSBSBSBSBSA",
    "AAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAQBQBQBQBQBQBQBQAAAAAA",
    "AAAAQBQBQBQBQBQBQBQAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAA",
    "EEEEEEBBBBBBBBBBBBBEEEEEE",
]
MARKET_ZONES = {
    "S": {"name_en": "shelf lanes", "name_zh": "货架区", "D": 0.4, "q": 1.0},
    "A": {"name_en": "main aisles", "name_zh": "通道", "D": 2.0, "q": 0.6},
    "F": {"name_en": "fresh produce", "name_zh": "生鲜区", "D": 0.3, "q": 1.8},
    "Q": {"name_en": "checkout queue", "name_zh": "收银台", "D": 0.1, "q": 2.0},
    "E": {"name_en": "entrance/exit", "name_zh": "出入口", "D": 1.5, "q": 0.8},
    "B": {"name_en": "shelf racks / counters", "name_zh": "货架/柜台",
          "D": 0.0, "q": 0.0, "barrier": True},
}

# ----------------------------------------------------------------------------
# Classroom: 15 m x 10 m, a = 1 m, 15 x 10 = 150 cells, 80 seats, N = 80.
# x runs front (platform) -> back, y runs from the window wall (row 0).
# 'd' = desk in the window-side block (q = 1.2), 'D' = desk in the inner
# block (q = 2.0): 40 + 40 = 80 desk cells.
CLASSROOM = [
    "PPIddddddddddII",
    "PPIddddddddddII",
    "PPIddddddddddII",
    "PPIddddddddddII",
    "POIIIIIIIIIIIII",
    "POIDDDDDDDDDDII",
    "PPIDDDDDDDDDDII",
    "PPIDDDDDDDDDDII",
    "PPIDDDDDDDDDDII",
    "PPIIIIIIIIIIIII",
]
CLASSROOM_ZONES = {
    "d": {"name_en": "desks (window side)", "name_zh": "课桌区（靠窗）",
          "D": 0.05, "q": 1.2},
    "D": {"name_en": "desks (inner)", "name_zh": "课桌区（中央）",
          "D": 0.05, "q": 2.0},
    "P": {"name_en": "platform", "name_zh": "讲台", "D": 0.8, "q": 1.5},
    "I": {"name_en": "aisles", "name_zh": "过道", "D": 0.2, "q": 0.8},
    "O": {"name_en": "obstacle", "name_zh": "障碍物", "D": 0.0, "q": 0.0,
          "barrier": True},
}

# ----------------------------------------------------------------------------
# Metro car (Chinese type-A): 22.5 m x 3 m, a = 0.5 m, 45 x 6 = 270 cells,
# crush load 310 (4.6 persons/m^2).  Five door pairs, 1.5 m wide, 4.5 m apart.
# D depends on the crowd density rho (persons/m^2):
#   T 0.02/rho, N 0.1/rho^2, G 0.5/rho, R 0.05/rho   (times D0).
def _metro_rows():
    nx, ny = 45, 6
    doors = [4, 13, 22, 31, 40]
    door_cols = set()
    for c in doors:
        door_cols.update([c - 1, c, c + 1])
    poles = set(doors) | {8, 17, 26, 35}
    rows = []
    for y in range(ny):
        row = ""
        for x in range(nx):
            if x in door_cols and y in (0, 1, 4, 5):
                ch = "G"
            elif y in (0, 5):
                ch = "T"
            elif x in poles and y in (2, 3):
                ch = "R"
            else:
                ch = "N"
            row += ch
        rows.append(row)
    return rows


METRO = _metro_rows()
METRO_ZONES = {
    "T": {"name_en": "seats", "name_zh": "座位", "D_coef": 0.02, "D_exp": 1,
          "q": 2.5},
    "N": {"name_en": "standing area", "name_zh": "站立区", "D_coef": 0.1,
          "D_exp": 2, "q": 2.8},
    "G": {"name_en": "doors", "name_zh": "车门", "D_coef": 0.5, "D_exp": 1,
          "q": 2.0},
    "R": {"name_en": "grab rails", "name_zh": "扶手", "D_coef": 0.05,
          "D_exp": 1, "q": 2.5},
}

SCENES = {
    "office": {
        "title_en": "Open-plan office", "title_zh": "开放式办公室",
        "a_m": 1.5, "size_m": [30.0, 20.0], "N_default": 100,
        "map": OFFICE, "zones": OFFICE_ZONES,
        "dwell": "8 h",
        "source": "illustrative assumption; area fractions "
                  "W/C/M/K/B = 60/15/10/5/10 % reproduced exactly",
    },
    "supermarket": {
        "title_en": "Supermarket", "title_zh": "中型超市",
        "a_m": 2.0, "size_m": [50.0, 30.0], "N_default": 150,
        "map": MARKET, "zones": MARKET_ZONES,
        "dwell": "27 min",
        "source": "illustrative assumption; floor plan and "
                  "barrier cells (shelf racks, counters) are a design choice",
    },
    "classroom": {
        "title_en": "Classroom", "title_zh": "教室",
        "a_m": 1.0, "size_m": [15.0, 10.0], "N_default": 80,
        "map": CLASSROOM, "zones": CLASSROOM_ZONES,
        "dwell": "50-75 min",
        "source": "illustrative assumption; 80 desk cells, "
                  "window block q=1.2, inner block q=2.0",
    },
    "metro": {
        "title_en": "Metro carriage", "title_zh": "地铁车厢",
        "a_m": 0.5, "size_m": [22.5, 3.0], "N_default": 310,
        "map": METRO, "zones": METRO_ZONES,
        "dwell": "20 min",
        "density_dependent_D": True,
        "source": "illustrative assumption; type-A car, 5 door "
                  "pairs; D evaluated at rho = N / 67.5 m^2",
    },
}


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, sc in SCENES.items():
        rows = sc["map"]
        ny, nx = len(rows), len(rows[0])
        assert all(len(r) == nx for r in rows), name
        sc["nx"], sc["ny"] = nx, ny
        cnt = Counter("".join(rows))
        sc["zone_counts"] = dict(cnt)
        with open(os.path.join(OUT, name + ".json"), "w", encoding="utf-8") as f:
            json.dump(sc, f, ensure_ascii=False, indent=1)
        tot = nx * ny
        print(f"{name}: {nx} x {ny} = {tot} cells, a = {sc['a_m']} m")
        for z, c in sorted(cnt.items()):
            print(f"   {z}: {c:4d}  ({100 * c / tot:5.1f} %)")


if __name__ == "__main__":
    main()

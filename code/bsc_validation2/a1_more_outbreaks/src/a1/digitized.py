"""Seat maps transcribed by eye from the published figures (see sources.md).

Every map is a list of strings, one per lateral seat line, one character per
aircraft row.  No model code reads the case symbols: `events.py` turns the
maps into geometry (occupied seats, sources, strata) and `outcomes()` returns
the case indicators, which only calibration/evaluation scripts call.

Common symbols
  '-' no seat        'e' empty seat
  'o' susceptible, not a case
  'C' susceptible, case (primary outcome)
  'P' susceptible, case only under the wider (sensitivity) outcome
  'I' infectious source (primary source set)
  'j' other person infectious/ill on board: not susceptible, not a source in
      the primary analysis (source in a sensitivity analysis)
  'x' on board but not in the risk set (immune, unknown status, not
      interviewed, travel-group member exposed before the flight)
"""

# ---------------------------------------------------------------------------
# Olsen 2003, flight CA112 Hong Kong-Beijing 15 Mar 2003, Boeing 737-300.
# Figure as redrawn in Hertzberg & Weiss 2016 Fig. 2 (gr2_lrg.jpg).
# Rows 1..22.  Lines: A, B, C | D, E, F (window, middle, aisle | aisle, middle, window).
# 'o' = no illness (X in the figure) or not interviewed (/ in the figure).
# symbols: X no illness, / not interviewed, G case, I index, . empty, - no seat
OLSEN_ROWS = list(range(1, 23))
OLSEN = {
    # row:        1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
    "A": list("--" ".""/X///GGXX" "GXXXXXGX/."),
    "B": list(".X" ".""//////GXX" "GIXX//X/X."),
    "C": list("--" "X""X//X/GX/G" "XX/X/XX/X/"),
    "D": list("--" ".""XX/GXXXGG" "G/X//XX/X."),
    "E": list(".." ".""/X/G/////" "/G/GGXX//."),
    "F": list("--" ".""/X/X///XG" "XXXGXX////"),
}

# ---------------------------------------------------------------------------
# Toyokawa 2022, domestic flight to Naha 23 Mar 2020, Boeing 737-800 (Fig. 2).
# Rows 1,2,3,5..31 (no row 4).  Seats A B C | F G H.  Index 23C.
# C = confirmed case, P = probable case (grey square), o = other passenger.
TOYO_ROWS = [1, 2, 3] + list(range(5, 32))
TOYO = {
    #  row   A B C F G H
    1: "oee---", 2: "ooeeee", 3: "oeeooo", 5: "ooooeo", 6: "ooooeo", 7: "Ceoooo",
    8: "oooooo", 9: "ooPooo", 10: "oooeoo", 11: "ooeooo", 12: "eeePoo", 13: "oPoPoo",
    14: "PCoCCC", 15: "oeoooe", 16: "ooeooo", 17: "ooooeo",
    18: "oooooo", 19: "oooooo", 20: "oPooeo", 21: "ooCCoo", 22: "CeCoeo", 23: "PeIooo",
    24: "ooeeoC", 25: "ooeCCo", 26: "oooCoo", 27: "ooeoeo", 28: "ooeeeo", 29: "CooPeo",
    30: "oeoeoo", 31: "ooeooo",
}
TOYO_HOUSEHOLDS = [["14B", "14F", "14G", "14H"], ["25F", "25G"]]

# ---------------------------------------------------------------------------
# Speake 2020, Sydney-Perth 19 Mar 2020, Airbus A330-200, economy mid cabin (Fig. 4).
# 17 seat rows (columns of the figure, front to rear); 8 seat lines, 2-4-2.
# I = infectious Ruby Princess passenger (A2 lineage; the paper's reference
#     set for its "within 2 rows" statement): three symptomatic, culture-positive
#     ('S' in SPEAKE_SYMPT) and one presymptomatic.
# j = other primary cases (B.1 lineage infectious, or non-infectious primary).
# C = flight-associated secondary case, P = possibly flight-associated secondary case.
SPEAKE = [
    #12345678901234567
    "eoooPooCCPeoCooej",   # line 1 (window)
    "ooooojoooPeooeoej",   # line 2
    "-oejeCo-ooooooooo",   # line 3
    "-eeoooIooeooooooo",   # line 4
    "-jejooIooooeooooo",   # line 5
    "-oojooIoCoooeoeeo",   # line 6
    "oooIooCoooooeoooe",   # line 7
    "oCejeoooCoeoooooo",   # line 8 (window)
]
SPEAKE_SYMPT = [(3, 6), (4, 6), (5, 6)]          # (line index, column index), zero-based
SPEAKE_B1_INFECTIOUS = [(1, 5), (4, 1)]

# ---------------------------------------------------------------------------
# Baker 2010, Los Angeles-Auckland 25 Apr 2009, Boeing 747-400 rear section (Fig.).
# 15 seat rows; 10 seat lines, 3-4-3.
# I = laboratory-confirmed symptomatic case during flight (9)
# j = suspected symptomatic case during flight (3)
# x = immune (2) or unknown status (5)
# C = laboratory-confirmed post-flight case (3), P = suspected post-flight case (1)
BAKER = [
    #123456789012345
    "---Ioooooooxooo",   # T1 window
    "oxjIooooooooojo",   # T2
    "oxIxIoCoooooooo",   # T3 aisle
    "oooxoIjoooo----",   # M1
    "oooPooIoooo----",   # M2
    "oCooooooooo----",   # M3
    "oooCooooooo----",   # M4
    "ooooIoxooxooooo",   # B1 aisle
    "ooooIeoooooooeo",   # B2
    "---Iooooooooooo",   # B3 window
]

# ---------------------------------------------------------------------------
# Young 2014, Cancun-Birmingham 2009, Boeing 767 (Fig. 1).  Rows 1..47; seats a b | c d e f | g h.
# I = infectious in flight (6), C = infected in flight (9 with known seat),
# o = passenger with data, x = passenger without data ("X") or immune.
YOUNG = {
    1: "--Coe---", 2: "ooeeo-oo", 3: "ooooe-oo", 4: "ooooo-oo", 5: "Ioooo-xx",
    6: "oooooooo", 7: "oooooooo", 8: "ooooIIoo", 9: "ooooooxo", 10: "oooooooo",
    14: "--oooo--", 15: "--xxxe--", 16: "--oooe--", 17: "eeooooxx", 18: "oooooooo",
    19: "oooxIooo", 20: "ooooooxo", 21: "ooooooxo", 22: "oooooooo", 23: "ooooooxo",
    24: "oooooooo", 25: "ooxooooo", 26: "ooxxxxxx", 27: "xxooxoxx", 28: "oooooooo",
    32: "xxxxxxxx", 33: "--xxxx--", 34: "xeoeeeee", 35: "oCxooxIo", 36: "oCooooox",
    37: "xxoooooo", 38: "ooCIxooC", 39: "ooxooooo", 40: "Cooooooo", 41: "oCoooooo",
    42: "ooxxoxxx", 43: "ooCoooox", 44: "oooooooo", 46: "ooxxoxxe", 47: "--oooC--",
}

# ---------------------------------------------------------------------------
# Hoehl 2020, Tel Aviv-Frankfurt 9 Mar 2020, Boeing 737-900 (Fig. 1).
# Row labels as printed; seats K J H | C B A.
# I = index case (7 members of the tourist group), x = tourist-group member tested
# negative at the airport or person not interviewed, o = interviewed passenger,
# C = likely transmission (tested positive).
HOEHL_ROWS = [1, 2, 3, 4] + list(range(21, 29)) + list(range(34, 53))
HOEHL = {
    1: "e-ee-e", 2: "o-ee-e", 3: "e-ee-e", 4: "e-ee-e",
    21: "oex---", 22: "oeeeoe", 23: "oooeeo", 24: "ooeeee", 25: "xxeeee", 26: "eeeeee",
    27: "eeoooo", 28: "oeoeeo", 34: "eeeoee", 35: "eeeeee", 36: "oeoooo", 37: "ooxxeo",
    38: "ooooxo", 39: "oooxxx", 40: "oooooo", 41: "oooooo", 42: "ooooeo", 43: "oexeeo",
    44: "oeeoeo", 45: "xIIxoo", 46: "xxxxoC", 47: "oIxCox", 48: "IIIxex", 49: "xIxoex",
    50: "ooooex", 51: "ooeeeo", 52: "eeeeee",
}

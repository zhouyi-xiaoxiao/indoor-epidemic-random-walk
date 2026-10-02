# Sources of the second set of outbreak records

Tables: `events.csv` (one row per outbreak), `strata.csv` (persons at risk and cases by distance stratum),
`seats.csv` (seat-map events: one row per person at risk and per source). They are written by
`code/v2_tidy_data.py` from the transcriptions in `code/bsc_validation2/a1_more_outbreaks/src/a1/digitized.py` and
`events.py`. All sources were retrieved on 1 October 2026 through Europe PMC. Seat maps were transcribed by eye
from the published figures; copies of the source texts and figures are not redistributed.

| Event | Reference | DOI | Where the numbers come from | Check against printed totals |
|---|---|---|---|---|
| F1 | Olsen SJ et al., N Engl J Med 2003;349:2416-22 | 10.1056/NEJMoa031349 | seat map as redrawn in Hertzberg and Weiss 2016, Fig. 2 (the article itself is not open access); abstract: 8/23 in the three rows in front against 10/88 elsewhere | 112 on board, 8/23, and 9/29 within two rows: exact |
| F2 | Kenyon TA et al., N Engl J Med 1996;334:933-8 | 10.1056/NEJM199604113341501 | abstract: 4/13 within two rows against 2/55 in the rest of the section; no seat map (seats placed at random inside the two strata) | typed from the abstract |
| F3 | Baker MG et al., BMJ 2010;340:c2424 | 10.1136/bmj.c2424 | seat-map figure of the rear section | 128 seats, 126 occupied, 103 non-cases: exact; 67 susceptibles within two rows (Hertzberg and Weiss: 67; Baker et al.: 57) |
| F4 | Young N et al., Influenza Other Respir Viruses 2014;8:66-73 | 10.1111/irv.12181 | Fig. 1 (seat map) | 216 passengers at risk with a seat against 220 in the text; 5/126 and 4/90 within and beyond two rows against the printed 5/131 and 4/89 |
| F5 | Toyokawa T et al., Influenza Other Respir Viruses 2022;16:63-71 | 10.1111/irv.12913 | Fig. 2 (seat map) | 141 persons on the map against 142 passengers in the text; 14 confirmed cases: exact |
| F6 | Speake H et al., Emerg Infect Dis 2020;26:2872-80 | 10.3201/eid2612.203910 | Fig. 4 (seat map of the mid cabin) | 28 window seats occupied with 7 cases, 8 of 11 secondary cases within two rows: exact; 110 mid-cabin passengers against 111-112 implied by the text |
| F7 | Swadi T et al., Emerg Infect Dis 2021;27:687-93 | 10.3201/eid2703.204714 | text and Fig. 3 (seats of everyone within rows 23-30; the others are known only to sit elsewhere and are placed at random) | 2 sources, 4 in-flight infections, 84 others: as in the text |
| F8 | Hoehl S et al., JAMA Netw Open 2020;3:e2018044 | 10.1001/jamanetworkopen.2020.18044 | Fig. 1 (seat map) | 71 interviewed, 17 group members negative, 7 not interviewed, 7 index cases: exact |
| W1 | Wong TW et al., Emerg Infect Dis 2004;10:269-76 | 10.3201/eid1002.030452 | text and table: 3/3 within 1 m of the index bed, 4/8 in the same cubicle, 0/8 elsewhere; Fig. 4 (floor plan) | typed from the text |
| W2 | Yu IT et al., Clin Infect Dis 2005;40:1237-43 | 10.1086/428735 | attack rates by bay: 13/20, 11/21, 6/33; Fig. 1 | typed from the text |
| P1 | Guenther T et al., EMBO Mol Med 2020;12:e13296 | 10.15252/emmm.202013296 | Appendix Table S2 (cumulative count and positives by 1-m distance): 5/9 within 4 m, 12/17 at 4-8 m, 1/22 at 8-12 m, 2/30 beyond | re-read from the appendix text |

Cross-check source and published constant near/far ratio 2.4: Hertzberg VS, Weiss H, Ann Glob Health
2016;82:819-23, doi:10.1016/j.aogh.2016.06.003.

Limits of the transcription. (i) The re-check re-read every seat symbol of F1, F5, F6, F7 and F8
against the published figures and re-typed the counts of P1; the maps of F3 and F4 were not re-read seat by
seat. (ii) About 22 passengers of F5 who were lost to follow-up are counted as non-cases. (iii) Most non-cases
of F8 were never tested. (iv) For W1 the published figure gives each student's bed; the analysis placed the
students at random beds within the three strata.

Events found and not used (reason fixed before any model was run): Murphy et al. 2020 (no identified source on
board); Foxwell et al. 2011 (survey response 42 %, no denominators by seat); Desenclos et al. 2004 (two cases);
Chen et al. 2020 (exposure before the flight); Choi et al. 2020, Bae et al. 2020, Eldin et al. 2020 (one or two
cases); Cho et al. 2016 (zone level, index case moved between zones); Moser et al. 1979 and Riley et al. 1978
(no positions published); Khanh et al. 2020 and Katelaris et al. 2021 (used in the first test).

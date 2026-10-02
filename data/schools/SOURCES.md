# Sources of the school data

| File | Reference | What it holds | Origin and licence |
|---|---|---|---|
| `matsumoto_students.csv` | Endo A, Uchida M, Hayashi N, Liu Y, Atkins KE, Kucharski AJ, Funk S. Within and between classroom transmission patterns of seasonal influenza among primary school students in Matsumoto city, Japan. PNAS 2021;118(46):e2112605118. doi:10.1073/pnas.2112605118 | one row per pupil (10,923 pupils, 29 schools, 455 classes): infected or not, day of onset (day 1 = 1 October 2014), school, grade, class, sex, class size, number of classes in the grade | tabular copy of `data/anonymizedstudents.jld2` of https://github.com/akira-endo/schooldynamics_FluMatsumoto14-15 (commit 658abfecc70440b0b849e179da457048082d4e72), written by `code/bsc_validation2/a4_zone_level/src/00_convert_jld2.py`; MIT licence, copy in `LICENSE_matsumoto.txt` |
| `lyon_mixing.json` | Stehle J, Voirin N, Barrat A, et al. High-resolution measurements of face-to-face contact patterns in a primary school. PLoS ONE 2011;6(8):e23176. doi:10.1371/journal.pone.0023176; Gemmetto V, Barrat A, Cattuto C. BMC Infect Dis 2014;14:695. doi:10.1186/s12879-014-0695-9 | contact time per pupil and day with the own class, the other class of the grade and the other grades, and the three shares (0.762, 0.104, 0.134) | computed by `code/bsc_validation2/a4_zone_level/src/01_lyon_mixing.py` from the SocioPatterns primary-school file, which is not redistributed (see `data/contacts/SOURCES.md`) |

Seoul call-centre building (cases by floor): Park SY et al., Emerg Infect Dis 2020;26(8):1666-70,
doi:10.3201/eid2608.201274, Table 1 (11th floor 94/216; floors 1-6 0/84; 7th 0/182; 8th 0/207; 9th 1/206;
10th 2/27; residents 0/201; other 0/20). The numbers are typed into
`code/bsc_validation2/a4_zone_level/src/05_seoul_floors.py`.

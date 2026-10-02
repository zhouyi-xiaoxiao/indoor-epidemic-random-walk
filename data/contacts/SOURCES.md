# Sources of the contact records

The contact records are the SocioPatterns face-to-face proximity datasets (RFID sensors, 20 s resolution),
distributed at http://www.sociopatterns.org/datasets/ under Creative Commons licences, some of which exclude
commercial use. They are **not redistributed** in this repository. `datasets.csv` lists the ten files used,
with their URL and the SHA-256 recorded when the protocol was frozen (the complete list of the files then
present is `code/bsc_validation2/a2_contact_mobility/data/raw_SHA256.txt`). `code/download_contact_data.py` downloads them
and checks the hashes.

| Dataset | Reference |
|---|---|
| office building, 2013 | Genois M, Vestergaard CL, Fournet J, Panisson A, Bonmarin I, Barrat A. Network Science 2015;3:326-47. doi:10.1017/nws.2015.10 |
| office building, 2015; scientific conference | Genois M, Barrat A. EPJ Data Science 2018;7:11. doi:10.1140/epjds/s13688-018-0140-1 |
| primary school | Stehle J, Voirin N, Barrat A, et al. PLoS ONE 2011;6(8):e23176. doi:10.1371/journal.pone.0023176 |
| high school | Mastrandrea R, Fournet J, Barrat A. PLoS ONE 2015;10(9):e0136497. doi:10.1371/journal.pone.0136497 |
| hospital ward | Vanhems P, Barrat A, Cattuto C, et al. PLoS ONE 2013;8(9):e73970. doi:10.1371/journal.pone.0073970 |

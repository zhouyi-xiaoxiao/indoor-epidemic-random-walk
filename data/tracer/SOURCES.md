# Sources of the tracer tables

All sources were accessed on 1 October 2026. The tables are numbers read from the sources; copies of the
source files are not redistributed, except as noted.

| Table | Reference | What it holds | Origin |
|---|---|---|---|
| `ou2022_B1_seats.csv`, `ou2022_B2_seats.csv` | Ou C, Hu S, Luo K, et al. Insufficient ventilation led to a probable long-range airborne transmission of SARS-CoV-2 on two buses. Build Environ 2022;207:108414. doi:10.1016/j.buildenv.2021.108414 | coach (B1) and minibus (B2): ethane tracer measured at 8 and 14 seats, tracer-validated flow simulation at all seats | Supplementary Fig. S3 and S4 (values read from the text records of the figures and checked on the rendered page); air change rates from Table S2 |
| `kinahan_tidy.csv`, `kinahan_777_*.csv`, `kinahan_767_*.csv` | Kinahan SM, Silcott DB, Silcott BE, et al. Aerosol tracer testing in Boeing 767 and 777 aircraft to simulate exposure potential of infectious aerosol such as SARS-CoV-2. PLoS ONE 2021;16(12):e0246916. doi:10.1371/journal.pone.0246916 | time-integrated counts of 1 um tracer particles per breathing-zone sensor and release, in-flight tests | figshare workbooks doi:10.6084/m9.figshare.13537349.v1 (777) and doi:10.6084/m9.figshare.13537319.v1 (767), licence CC BY 4.0; the tables are extracts of those workbooks |
| `woodward2022_fig9_points.csv`, `woodward2022_fig9_curve.csv` | Woodward H, de Kreij RJB, Kruger ES, et al. An evaluation of the risk of airborne transmission of COVID-19 on an inter-city train carriage. Indoor Air 2022;32(10):e13121. doi:10.1111/ina.13121 | steady concentration of nebulised salt aerosol along the saloon, release in the middle and at the end | digitised from Fig. 9 by colour segmentation (`code/bsc_validation2/a3_tracer_physics/src/p03_woodward_digitise.py`); saloon dimensions and 11-15 air changes per hour from the text |

Used without a table of their own: Li Y, Qian H, Hang J, et al., Build Environ 2021;196:107788,
doi:10.1016/j.buildenv.2021.107788, Tables 1 and 3 (tracer gas at the restaurant tables; transcribed in
`data/bsc_outbreaks/tables/li2021_restaurant_tables.csv`); Cheng KC, Acevedo-Bolton V, Jiang RT, et al.,
Environ Sci Technol 2011;45:4016-22, doi:10.1021/es103080p (K/L^2 = 0.52 ACH + 0.31 per hour, p. 4020);
Khanh NC et al., Emerg Infect Dis 2020;26:2617-24, doi:10.3201/eid2611.203299, Fig. 1 (seat map of the
business class of flight VN54, typed into `code/bsc_validation2/a3_tracer_physics/src/a3lib.py`).

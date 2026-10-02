"""Mapping of outbreak records onto the lattice model (scene, lattice, occupancy, duration).

Default scene definitions (illustrative values, as in the scene files of code/bsc_sim/scenes/):

  scene        room (m)    a (m)  lattice   people   dwell        D/D0 range         q range
  office       30 x 20     1.5    20 x 13   100      8 h          0.3 .. 1.5         0.8 .. 1.5
  supermarket  50 x 30     2.0    25 x 15   100-200  27 +/- 12 min 0.1 .. 2.0        0.6 .. 2.0
  classroom    15 x 10     1.0    15 x 10   80       50-75 min    0.05 .. 0.8        0.8 .. 2.0
  metro        22.5 x 3    0.5    45 x 6    310      20 min       0.02/rho .. 0.5/rho 2.0 .. 2.8

with D0 = 1 lattice^2/day, beta = 0.5/day, gamma = 0.14/day.

Lattice sizes below are nx = round(L/a), ny = round(W/a) with the default a for the mapped
scene.  `T_days` is the exposure duration expressed in the model's time unit (days).
`role` is the record's role in the validation protocol (code/bsc_outbreaks/PROTOCOL_PROPOSAL.md).
"""

SCENE_DEFAULTS = dict(
    beta_per_day=0.5, gamma_per_day=0.14, D0_lattice2_per_day=1.0,
    scenes=dict(
        office=dict(a_m=1.5, lattice=(20, 13), n=100, dwell_h=8.0, D_seated=0.3, D_max=1.5, q_max=1.5),
        supermarket=dict(a_m=2.0, lattice=(25, 15), n=150, dwell_h=0.45, D_seated=0.1, D_max=2.0, q_max=2.0),
        classroom=dict(a_m=1.0, lattice=(15, 10), n=80, dwell_h=1.0, D_seated=0.05, D_max=0.8, q_max=2.0),
        metro=dict(a_m=0.5, lattice=(45, 6), n=310, dwell_h=0.33, D_seated=0.02, D_max=0.5, q_max=2.8),
    ),
)

# id -> mapping
MAPPING = {
    "O1_seoul_callcentre": dict(
        model_scene="office", a_m=1.5, lattice="28 x 28 desk-pitch units (floor plan 900 x 897 px, desk pitch 32.5 px); physical pitch not reported (plausibly 1.0-1.6 m), central service core = obstacle block",
        n_agents=216, T_days="11-12 working days x ~8 h (21 Feb - 8 Mar 2020); multi-generation",
        index_position="unknown (seed anywhere in north wing; marginalise)", mobility="seated at desk (W zone), corridors, shared lobby/elevators",
        observables=["floor attack rate 94/216", "wing contrast 79/137 vs 4/62 (seat-resolved map)", "onset curve (86 dated 11th-floor cases)", "within-wing clustering statistic (join counts)", "other floors as negative control (1/595 on floors 7-9)"],
        role="held-out multi-generation test (day-scale SIR is appropriate here)",
        obstacles="no room scale, no ventilation, no hours; SIR lacks a latent period although onsets are delayed ~5 d; 10 of 94 cases have no seat",
    ),
    "O2_zurich_office": dict(
        model_scene="office", a_m=1.5, lattice="9 x 9 (175 m2 team area ~ 13.2 m square) + conference room 4 x 4 (30 m2)",
        n_agents=13, T_days="11 h = 0.46 d (single index, first generation only)", index_position="known (Fig. 1 of source, not digitized)",
        mobility="seated + frequent shared-desk work + 1-h meeting", observables=["first-generation attack rate 8/12 to 10/12"],
        role="held-out level test (overall attack rate)", obstacles="behavioural close contact dominates; n = 12"),
    "O3_england_office": dict(
        model_scene="office", a_m=1.5, lattice="not mappable (no geometry)", n_agents=40, T_days="weeks; multi-generation",
        index_position="unknown", mobility="open-plan office", observables=["cumulative attack rate 22/40"],
        role="qualitative only", obstacles="abstract-level data"),
    "O4_german_meatplant": dict(
        model_scene="office (analogue: fixed work stations)", a_m=1.5, lattice="21 x 6 (32 m x 8.5 m line area)",
        n_agents=78, T_days="3 shifts (~8 h each) = ~1 d of contact; first generation", index_position="known to authors (proximal half), not public",
        mobility="fixed stations", observables=["attack rate among fixed-station workers 20/78", "radial profile of cumulative attack rate vs distance (source Fig. 3B, to be digitized)"],
        role="held-out shape test (distance kernel) once Appendix Table S2 is transcribed", obstacles="index position withheld; cold recirculated air unlike an office"),
    "O5_tb_office": dict(
        model_scene="office", a_m=1.5, lattice="not mappable (no geometry)", n_agents=68, T_days="4 weeks (~160 working hours)",
        index_position="unknown", mobility="office work", observables=["attack rate 27/67"], role="qualitative only (different pathogen)",
        obstacles="abstract-level data"),
    "S1_liaocheng_supermarket": dict(
        model_scene="supermarket", a_m=2.0, lattice="scene default 25 x 15 (store size not reported)", n_agents="120 staff + through-flow of customers (scene: 100-200 present)",
        T_days="customers: one visit (scene dwell 27 min = 0.019 d); staff: daily shifts over >= 1 week", index_position="unknown (staff)",
        mobility="customers walking (Robin in/out boundary); staff semi-fixed", observables=["customer attack rate 0/8224 (upper bound)", "staff attack rate 11/120", "staff/customer contrast"],
        role="negative control for short-dwell risk + level test for staff", obstacles="no geometry, no dwell times, primary case unidentified"),
    "S2_tianjin_dept_store": dict(
        model_scene="supermarket (analogue)", a_m=2.0, lattice="29 x 29 if the 3,400 m2 floor is taken as square (shape not reported)", n_agents="~200 staff + customers",
        T_days="staff: ~7 days of shifts (20-26 Jan 2020)", index_position="unknown", mobility="staff at counters, customers walking",
        observables=["staff attack rate 6-7/194", "ratio of customer to staff cases 21-23 : 6-7 (no customer denominator)"],
        role="qualitative / order-of-magnitude level check", obstacles="inconsistent counts across sources; no customer denominator"),
    "S3_wenzhou_mall": dict(
        model_scene="supermarket (analogue)", a_m=2.0, lattice="not mappable", n_agents=None, T_days="multi-day", index_position="unknown", mobility="mixed",
        observables=["existence of customer infections without direct contact"], role="qualitative only", obstacles="no denominators"),
    "S4_boston_grocery": dict(
        model_scene="supermarket", a_m=2.0, lattice="not mappable", n_agents=104, T_days="cumulative occupational", index_position="n/a", mobility="staff",
        observables=["staff point prevalence 21/104", "customer-facing vs not (OR 5.1)"], role="qualitative only", obstacles="prevalence, not attack rate"),
    "C1_marin_classroom": dict(
        model_scene="classroom", a_m=1.0, lattice="~10 x 10 (5 rows of desks at 1.83 m spacing => ~9 m x 9 m; room size not reported)",
        n_agents=25, T_days="2-5 school days x ~6.5 h (17-21 May 2021); effectively first generation", index_position="known: teacher's desk at front",
        mobility="students seated (D zone), teacher in podium zone (P)", observables=["class attack rate 12/24 (12/22 tested)", "front-two-rows vs back-three-rows 8/10 vs 4/14"],
        role="calibration candidate for the classroom kernel OR held-out shape test", obstacles="hours not reported; masks on students; Delta variant; HEPA + open windows"),
    "C2_jerusalem_highschool": dict(
        model_scene="classroom", a_m=1.0, lattice="7 x 6 to 7 x 7 per classroom (39-49 m2); 30+ classrooms coupled through teachers/schoolyard",
        n_agents="35-38 students + teacher per class; 1,164 students, 152 staff", T_days="8-9 school days x 6.5 h; 1-2 generations",
        index_position="unknown", mobility="seated; between-class coupling via teachers and breaks",
        observables=["school attack rate 153/1161", "grade profile (20.3, 17.3, 32.6, 4.5, 3.1, 1.6 %)", "max cases per class 20, 14, 13, 13"],
        role="held-out multi-generation test of the multi-room (scene-coupling) extension", obstacles="per-class denominators and positions not reported"),
    "C3_utah_elementary": dict(
        model_scene="classroom", a_m=1.0, lattice="scene default scaled to ~25 pupils (seat spacing 0.9 m)", n_agents="~20 contacts per index (1041/51)",
        T_days="1-2 school days before isolation of the index", index_position="unknown", mobility="seated, masked",
        observables=["per-index secondary attack rate 5/728 = 0.7%"], role="negative control (masked classrooms): model with mask factor must not exceed ~2%",
        obstacles="exposure hours and ventilation not quantified"),
    "C4_measles_school": dict(
        model_scene="classroom (multi-room)", a_m=1.0, lattice="not mappable from abstract", n_agents=None, T_days="school days of one infectious period",
        index_position="known classroom", mobility="seated", observables=["28 first-generation cases spread over 14 classrooms"], role="qualitative only (shows need for a non-local/ventilation coupling term)",
        obstacles="abstract-level data; different pathogen"),
    "C5_cheonan_fitness": dict(
        model_scene="classroom (analogue)", a_m=1.0, lattice="8 x 8 (60 m2)", n_agents="5-22 per class", T_days="50 min per class (0.035 d), 2 classes per week",
        index_position="instructor at front", mobility="vigorous movement in place", observables=["attack rate 57/217", "no cases in classes with < 5 participants or low-intensity classes"],
        role="held-out level test (density and activity dependence)", obstacles="per-class denominators in appendix (not transcribed)"),
    "T1_zhejiang_bus": dict(
        model_scene="metro (analogue: bus)", a_m=0.5, lattice="23 x 5 (15 rows x 0.75 m = 11.25 m; width ~2.5 m assumed, not reported)",
        n_agents=68, T_days="100 min = 0.069 d", index_position="known: row 8, middle seat of 3-seat side", mobility="seated (T zone); no movement",
        observables=["bus attack rate 23/67", "near/far rows 14/33 vs 9/34 (and 11/23 vs 12/44)", "unexposed bus 0/60", "window-seat sparing (qualitative)"],
        role="held-out shape + level test", obstacles="bus width/height not printed; seat-level map not digitized"),
    "T2_hunan_coach": dict(
        model_scene="metro (analogue: coach)", a_m=0.5, lattice="23 x 5 (11.3 m x 2.5 m)", n_agents=49, T_days="150-200 min = 0.10-0.14 d",
        index_position="known: second rear row, window", mobility="seated", observables=["attack rate 7/48 (7/46)", "rear vs front rows 3/19 vs 4/26", "driver side vs opposite side 6/24 vs 1/20"],
        role="held-out shape + level test (measured ventilation 1.72 L/s per person)", obstacles="two reports disagree on denominator and duration"),
    "T3_hunan_minibus": dict(
        model_scene="metro (analogue: minibus)", a_m=0.5, lattice="11 x 5 (5.5 m x 2.5 m)", n_agents="13-18", T_days="60 min = 0.042 d",
        index_position="known", mobility="seated", observables=["attack rate 2/12 (2/17)", "infectee distances 1.5 m and 4.5 m"], role="held-out level test (same index as T2: tests duration + ventilation scaling with infectiousness fixed)",
        obstacles="n = 12-17"),
    "T4_china_hsr": dict(
        model_scene="metro (analogue: train coach)", a_m=0.5, lattice="6 columns (A B C aisle D F) x rows of one coach (rows per coach not reported); 1 site per seat (0.5 m wide); row spacing per source 0.4 m",
        n_agents="full coach (seat count not reported in source)", T_days="0.13-13.8 h; mean 2.1 h = 0.088 d", index_position="origin of the (d_row, d_col) offset matrix",
        mobility="seated; occasional aisle walks", observables=["23-cell attack-rate matrix by (rows apart, columns apart)", "per-hour slopes by seat class", "seat re-use 1/1342"],
        role="CALIBRATION of the distance kernel and of the duration dependence", obstacles="cell denominators not transcribed; 0.4 m row spacing implausible as seat pitch"),
    "T5_vn54_flight": dict(
        model_scene="metro (analogue: aircraft cabin)", a_m=0.5, lattice="business cabin 28 seats (layout in source Fig. 1, not digitized)", n_agents=21, T_days="10 h = 0.42 d",
        index_position="known seat", mobility="seated; aisle/toilet trips", observables=["near/far 11/12 vs 1/8", "adjacent cabins 0/35 and 2/145"], role="held-out shape test (long duration)",
        obstacles="cabin geometry not transcribed; cabin air exchange far higher than metro"),
    "T6_influenza_airliner": dict(
        model_scene="metro (analogue)", a_m=0.5, lattice="not mappable from abstract", n_agents=54, T_days="3 h = 0.125 d", index_position="unknown", mobility="seated/moving",
        observables=["72% clinical attack rate with zero ventilation"], role="qualitative only", obstacles="abstract-level data; different pathogen"),
    "T7_japan_no_train_clusters": dict(
        model_scene="metro", a_m=0.5, lattice="scene default 45 x 6", n_agents="up to 310", T_days="20 min = 0.014 d (scene dwell)", index_position="n/a", mobility="standing/seated",
        observables=["no rail cluster of >= 5 cases detected nationally in 11 weeks"], role="weak negative control for the metro scene",
        obstacles="non-detection is not absence"),
    "R1_guangzhou_restaurant": dict(
        model_scene="classroom-type seated room (no restaurant scene)", a_m=1.0, lattice="17 x 8 (17 m x 8.1 m)", n_agents=89,
        T_days="per-table exposure 48-89 min (index table 82 min = 0.057 d)", index_position="known: seat A1 at table A (window side)", mobility="seated at tables; staff walking",
        observables=["per-table attack rates (18 tables) vs per-table exposure time", "zone contrast 5/11 vs 0/68", "neighbour contrast 5/16 vs 0/63", "CFD/tracer-normalised exposure per table as an independent covariate"],
        role="CALIBRATION candidate (best single event) or primary held-out test", obstacles="airflow zoning (not distance) organises the risk; table A confounded by household exposure"),
    "R2_jeonju_restaurant": dict(
        model_scene="classroom-type seated room", a_m=1.0, lattice="9 x 11 (9.2 m x 10.5 m)", n_agents="14 (13 exposed + infector)", T_days="5 min and 21 min overlaps (0.0035 d, 0.015 d)",
        index_position="known (table near door 2)", mobility="seated", observables=["2/13 infected at 4.8 m and 6.5 m within 5-21 min"], role="stress test: advective transport (not representable by isotropic diffusion)",
        obstacles="requires a drift (directed-airflow) term"),
    "H1_skagit_choir": dict(
        model_scene="classroom-type seated room", a_m=1.0, lattice="15 x 12 (180 m2; shape not reported)", n_agents=61, T_days="2.5 h = 0.104 d", index_position="known to investigators, withheld",
        mobility="seated in rows; one re-seating for sectionals; 15-min break", observables=["attack rate 32/60 to 52/60", "absence of a spatial pattern"], role="held-out level test + well-mixed limit check",
        obstacles="seat map not public; 20 probable cases"),
    "H2_sydney_church": dict(
        model_scene="classroom-type seated room", a_m=1.0, lattice="not mappable (dimensions not reported)", n_agents="120-215 per service", T_days="1 h per service = 0.042 d",
        index_position="known (choir loft)", mobility="seated", observables=["per-service attack rates 5/215, 7/120, 0/173", "all cases within one 70-degree sector up to 15 m"], role="qualitative shape test (anisotropy)",
        obstacles="no sector denominators"),
    "X1_diamond_princess": dict(
        model_scene="none", a_m=None, lattice="not mappable", n_agents=3711, T_days="weeks", index_position="n/a", mobility="n/a", observables=["cumulative positives 634/3711"],
        role="outside model scope; recorded for reference", obstacles="ship-wide, quarantine intervention"),
    "X2_singapore_dormitories": dict(
        model_scene="none", a_m=None, lattice="not mappable", n_agents=295000, T_days="months", index_position="n/a", mobility="n/a", observables=["17,758/295,000 by 6 May 2020; one dormitory 19.4%"],
        role="outside model scope; recorded for reference", obstacles="residential, months-long"),
}

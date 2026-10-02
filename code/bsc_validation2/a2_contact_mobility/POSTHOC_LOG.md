# POSTHOC LOG — the contact-record analysis
Every choice made after `PREREG_FREEZE.json`.

1. (erratum, before any real-data run) PREREGISTRATION §7 says "0 of 40" wrong non-nested model×replicate cases; the
   table actually contains 11 such model×truth cells × 4 replicates = 44 (0 of 44 reach the A1 criterion; bound 6.6 %).
2. Stage 1 (base models, DEV) ran as frozen: `out/s1_{InVS13,LyonSchool,LH10}.json`. All four models fail E on all three.
3. Stage 2 begins (exploratory, DEV only). F1 = home-anchored walk `src/fix1.py`: homes fixed across days, groups'
   homes in contiguous strips (`a2lib.make_homes`, written before the freeze), with probability b the attempted move
   is the step towards home. Fitting rule chosen BEFORE seeing any F1 output on test days: (M, p) iterated to match
   training total contact time and mean duration; b from a fixed grid to match the training mean daily degree.
   Runs tagged `s2a_*`.
4. F1 result on DEV (`out/s2a_*.json`): epidemic error roughly halved, E still fails on all three (LH10 passes 2/4
   scenarios); S3 mean degree now matched (fitted), S5 passes on LH10 and LyonSchool and fails on InVS13; still failing: S2 burstiness,
   S4 strength heterogeneity, S6 within-group share (strips are too thin relative to the localisation length), S1 tail.
5. F2 = F1 with spatially CLUSTERED homes (`src/fix2.py`, model `RWclus`): group centres on a regular grid, each home =
   centre + Gaussian(σ) rounded/clipped. One more parameter σ, fitted with b on a fixed 2-D grid to the training mean
   degree and training within-group share (objective |log deg ratio| + |log within ratio|). Datasets without groups keep
   uniform homes (σ irrelevant). Rationale: desk-/class-mates share or neighbour home sites → long, repeated contacts.
6. F2 result on DEV (`out/s2b_*.json`): S3 mean degree, S6 within-group share and S2 recurrence now pass (fitted or
   induced); E error reduced again (mean R_index ratio 1.24 / 1.15 / 1.13 for InVS13 / LyonSchool / LH10; base walk
   1.53 / 1.25 / 1.36) but E still fails on all three (0, 1, 1 scenarios of 4). Still failing: S1 long-contact share
   (model ≈ 0, data 0.03–0.14), S2 burstiness (0.14–0.27 vs 0.36–0.61), S4 strength heterogeneity.
7. (resumed stretch of work, 2026-10-01) Interim memo v0 written from stages 1–2b. `src/03_tables.py` added (formatting only).
8. F3 = F2 + slow home sites (`src/fix3.py`), the pre-registered candidate class "sticky sites / broader heterogeneity
   of hopping". Two variants, both tried on DEV, decided BEFORE seeing either output:
   * `RWstickZ` (zone): every site that is somebody's home is a slow site with hop probability ρ·p; all other sites
     have p; moves between sites of different speed use the harmonic-mean rule. This is exactly the model's
     heterogeneous D(r) with the seating area as the slow zone.
   * `RWstickO` (own): a walker standing on its OWN home site leaves with probability ρ·p; elsewhere p (person-specific,
     not a D(r) field — outside the model's formulation, recorded as such).
   One more parameter ρ ∈ {1, 0.3, 0.1, 0.03, 0.01}. Automatic fit on training days: stage A = the F2 fit of (b, σ)
   with ρ = 1; stage B = scan ρ at that (b, σ); stage C = rescan the 3×3 grid neighbourhood of (b, σ) at the chosen ρ.
   (M, p) are re-iterated inside every fit to match training total contact time and mean event duration. Objective =
   |log ratio| of mean daily degree + within-group share (if groups) + share of contact time in events ≥ 5 min
   (floored at 0.002). Rationale: people sit at desks/seats for long stretches → long contacts with seat neighbours,
   bursty gaps, unequal contact time. Selection rule between Z and O fixed now: prefer Z (inside the model's framework)
   unless O passes strictly more DEV family×dataset cells in total (S-families + E).
9. F3 result on DEV (`out/s2c_*.json`, log `out/s2c_dev.log`). Cells passed (S-families + E, three datasets):
   RWstickZ 11, RWstickO 15 → by the rule of item 8 the person-specific variant O is the primary modification.
   RWstickO: InVS13 6/7 S-families and E 4/4 scenarios; LH10 3/7 and E 3/4; LyonSchool 4/7 and E 0/4 — on the school
   the model now UNDER-predicts outbreaks (mean attack 0.05 vs 0.11 at R0=1.5), i.e. the fix over-corrects there; the
   training fit is poor (fit error 0.63: degree 31 vs 47) because the ρ grid is coarse (factor 3) and the outcome is
   steep in ρ. RWstickZ (the model's own D(r) mechanism) fails E on all three (2, 0, 1 scenarios).
10. Fitting refinement, decided after item 9 and before any CONF run: stage D = at the chosen (b, σ) additionally try
   ρ·√3 and ρ/√3 (geometric midpoints of the grid) and keep the best training error (`fix3.make_calfun(..., refine=True)`,
   tag `s2d`). No change to the model. This is the last exploratory step; whatever s2d gives is frozen.
11. Stage D result on DEV (`out/s2d_*.json`): the refined fit changed only InVS13/RWstickO (ρ 0.1 → 0.173, training
   error 0.46 → 0.23) and made its TEST performance worse: E 4/4 → 0/4 scenarios (R_index ratio 1.19), S4 fails again.
   LH10 (E 3/4) and LyonSchool (E 0/4, under-prediction) unchanged. So with the better-fitting parameters RWstickO passes
   E on 1 of 3 DEV datasets, with the coarse grid on 2 of 3; the epidemic outcome is steep in ρ and ρ is only weakly
   pinned by the training targets. As committed in item 10, the refined procedure is what is frozen. Reported in full.
12. SECOND FREEZE (`FIX_FREEZE.json`), before any model or statistic touches a CONF dataset. Carried to CONF:
   base models WM, BLOCK, RW0, RWhet and three modifications, each with its automatic train-only fit:
   F2 `RWclus` (fix2.calfun), F3-Z `RWstickZ` and F3-O `RWstickO` (fix3 with stage D; modules fix3rZ / fix3rO, which
   reuse the stored F2 fit of the same dataset). PRIMARY modification for the ladder of PREREG §6/§8: RWstickO;
   RWstickZ is reported as "the modification that stays inside the model's D(r) formulation"; RWclus as the
   lower-parameter alternative. Run plan (only to respect the 20-min limit; no analytic content): tags s3a = base
   models, s3b = RWclus, s3c = RWstickZ, s3d = RWstickO; same seed 20261001, same real-network SIR in every tag.
   Expectation written now: A0 for every model on CONF; RWstickO passes E on at most 1 of 3 and fewer than 6 S-families
   everywhere; S1 (long contacts), S2 (burstiness), S4 fail for all walk variants; SFHH (no groups, conference)
   is where the home-anchored variants should be least appropriate.
13. Stage 3 (confirmatory, CONF = InVS15, Thiers13, SFHH) run exactly as frozen: `out/s3{a,b,c,d}_*.json`, logs alongside.
   Ladder (`out/ladder.json`): A0 for every model, over all six datasets and over CONF. E passed on 0 of 3 CONF datasets
   by every model including the primary modification RWstickO (2, 0, 0 of 4 scenarios). The real-network SIR results are
   identical in the four tags (checked).
14. Degenerate fit, reported not repaired: Thiers13 / RWstickO — the automatic fit ended at p = 0.8 (the cap) and
   M = 26 912 sites; on the test days it produces 2.2× the observed contact time (S0 FAIL), within-group share 1.00,
   and in the epidemic comparison its contact time is 2.34× the real one, so that comparison is NOT at matched contact
   rate. Its Thiers13 epidemic numbers are therefore not interpretable as support or refutation of the mechanism.
15. PREREG §8 item 4 ("adds something beyond BLOCK" iff E on at least as many CONF datasets as BLOCK and more S-families
   in total): RWstickO 0 vs BLOCK 0 datasets on E, 10 vs 6 S-family passes. The letter of the criterion is met only
   because both fail E everywhere; the memo reports this as "criterion met vacuously", not as a success.
16. Not pre-registered, descriptive only (memo §5): mean |log R_index ratio|, mean |Δattack|, mean |ΔP_major| over the
   24 dataset×scenario cells per model; figures `figs/fig{1,2,3}_*_{en,zh}.{pdf,png}` (`src/04_figs.py`; Fig. 2 uses one
   regenerated realisation per model, seed fixed, cached in `out/fig2_dist.npz`).
17. Scope not covered (as declared in PREREG §2): Copenhagen Networks Study, ATC / retail trajectory data, zone-occupancy
   statistics (SocioPatterns files carry no positions), Malawi and Hypertext files (downloaded, unused).

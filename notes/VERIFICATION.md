# Verification

How the results of the article (`paper.pdf`) and of its Supplementary Material (`supplement.pdf`) were checked,
and what the checks found. Section 1 is a summary. Section 4 goes through the parts of the work one by one and
lists every item that was examined, with its outcome; the `% src:` comments of the LaTeX files cite these items by
the heading of the part and the number or title of the item. Numbers that the article takes from these items are
also collected, under descriptive keys and with the item that reports each of them, in `data/checks/<part>.json`.

## 1. How the results were checked

Each part of the work was re-checked with separately written code, within the same project and on the same
machine. This is not an external or third-party verification, and the article has not been through external
peer review. The parts are: theory, simulations, references, the outbreak dataset, the first outbreak test, and
the four further analyses (second outbreak test, contact records, tracer measurements, zone level). In each check

* proofs were derived again and theorems tested on random instances, with a search for counter-examples;
* computations were implemented again: a second exact simulator, and second implementations of the primary
  analyses, with their own kernels and random numbers; re-runs of the stored analyses, comparisons of kernels
  and comparisons added by the checks partly reuse modules of the first analysis (data loaders of the contact
  and tracer analyses, event geometries, the transcribed seat maps, the simulator for timing benchmarks and schedules);
* numbers taken from the literature were typed again from the sources, and seat maps were read again from
  the published figures;
* the freeze records of the protocols were compared with the hashes and modification times of the files;
* generic competitors and analyses that the first analysis had not run were added where they could change
  the reading of a result; these additions are post hoc and are marked as such in the article.

Each item examined received one of four outcomes: *confirmed*, *partially confirmed*, *refuted* or
*unverifiable*. Where a check corrected a result of the first analysis, the article reports the corrected result.

## 2. Code of the checks

The code of the checks is in `code/bsc_theory/verify/`, `code/bsc_sim/verify/`, `code/bsc_outbreaks/verify/`
and `code/bsc_validation/verify/` and in `code/checks/` (the four further analyses, with a `README.md`). The
stored outputs of the checks are next to their code.

## 3. Numbers that come from the checks

Numbers that the article takes from a check rather than from a results file of the analyses are collected in
`data/checks/*.json`, one file per part, each value with the item of Section 4 that reports it. When the
repository is assembled, every such value is tested to be written in that item.

## 4. The checks, part by part

In the tables below, "the first analysis" is the analysis under examination and "the check" is its
re-examination with separately written code. Numbers in brackets after a quantity of the check are those of the
first analysis where the two are compared. File names such as `verify/v1_theorems.py` are relative to the
directory of the check named in Section 1. The items are records of the checks, reworded impersonally; where an
item differs from the article, the article prevails, as the notes in Sections 4.5, 4.6 and 4.8 explain. Items that
concern only material that is not part of the article are omitted; the items and further points shown are
numbered consecutively in each section. Each section ends with a summary written for this repository.

Terms. "Table Vn" is a table of the first outbreak test in `data/bsc_validation/tables.md`. The protocols are
called pre-registrations in their files, but they were not registered with any registry: they were frozen by
SHA-256 hashes of files on the author's machine, written with the outbreak data known (data-aware), and their
timing is self-attested. The protocol files of this repository are redacted copies (`code/REDACTIONS.md`); where
an item below says that the frozen hashes match, they were recomputed on the unredacted originals, which the
author can provide to the editor and referees, and the released copies hash differently.

### 4.1 Theory

*How it was checked.* Every proof was derived again symbolically, line by line, without code. Each theorem was then tested with separate code
on several hundred to several thousand random instances (symmetric and non-reversible movement, sparse and
dense generators), and counter-examples to the bounds and to the monotonicity in mobility were searched for
with an optimiser; where double precision was inconclusive the instances were evaluated again in 60-digit
arithmetic. All worked examples were recomputed.

*What was found.* No theorem was refuted. Apparent violations found by the optimiser were traced to the
numerics of the search and vanish in exact arithmetic. Two items were partially confirmed: the agreement of
the annealed encounter-factor formula with simulation is within about 5 %, not exact; and the comparison
of targeted ventilation with "highest efficiency first" depends on a random tie-break.

*Which office layout.* The worked examples of the check (items 2, 5, 6, 7, 11, 12, 13, 14, 16
and 17 and further point 2) were computed on layout O1, one of the two office plans of the first analysis with the cell counts
of the office scene but a different arrangement of the zones (R0 = 5.109 at D0 = 1; the second plan, O2, gives
5.041). The article uses the office scene of `code/bsc_sim/scenes/office.json` (R0 = 5.107),
whose values were computed again by `code/plan_theory_checks.py`, `code/s3_r0_checks.py` and
`code/s4_geometry_checks.py`; that is where the numbers of the article come from (for example R0 exceeds the
well-mixed value by 0.06 % at D0 = 100, the meeting room has an exit time of 114 days and 98.1 % of the elasticity,
the spectral gap is 0.0615, the total-variation distance is 0.996 and targeting gains 0.35 to 0.48). The per-visit
numbers of further point (1) are those of the uniform-leak formula βq̄/(γ + 1/T); the per-visit numbers of the
article are counted in simulations of a single stay.

*Items: 17* (15 confirmed, 2 partially confirmed).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | R0 = ρ(βQ(γI−M)^-1), Laplace-propagator interpretation, sign(R0−1) = sign s(M+βQ−γI), unique root of s(M+βQ/R−γI) = 0, comparison bound I(t) ≤ e^{Jt}I(0) (Theorem 3.1(a)–(d) as stated) | confirmed | Proof re-derived line by line (M-matrix facts P1–P5 used correctly; claim h(c)<0 ⇔ cρ(K)<1 holds in both directions; comparison principle valid because S/N≤1 and J is Metzler). Own implementation (verify/v1_theorems.py, 600 random instances, 4 generator families incl. non-reversible and dense): max \|s(M+βQ/R0−γ)\| = 2.1e-14, 0 sign failures, quadrature of ∫e^{-γt}e^{Mt}dt vs (γ−M)^-1 rel. err 1.2e-14. Nonlinear ODE re-implemented (verify/v3_dynamics.py): final attack 3.86e-6, 9.88e-6, 5.09%, 19.83%, 72.62% for R0 = 0.80…2.00 — identical to the reported values. Part (d) is a comparison bound: for R0 < 1 the infectives decay exponentially from t = 0 (in a closed SIR model I → 0 for every R0, so (d) is not an extinction statement); the check did not redo the 4e6-lifetime Monte-Carlo sojourn test (the identity it tests was checked directly). |
| 2 | Rayleigh formula for reversible movement; R0 = λmax[βΛ_nm/√((γ+μ_n)(γ+μ_m))]; the largest diagonal modal term is only a lower bound; Galerkin truncation monotone (Supplementary Theorem S2.4, Corollary S2.5) | confirmed | Symmetrisation by Π^{1/2}, Dirichlet-form identity and Cauchy interlacing checked by hand; 'exact iff q constant' is correct (a diagonal Q commuting with an irreducible A must be constant, and φ_n, n≥1, cannot be the Perron vector). Numerics on 300 symmetric/departure instances: Rayleigh vs NGM rel. err 1.6e-14; π ∝ 1/D for the departure rule asserted; diagonal modal max never above R0 (max rel. excess −8e-11); 0 Galerkin monotonicity failures. Office layouts O1 / O2 (not the office scene of the article) reproduced exactly: diagonal max 4.9848 / 4.7319, 3 modes 5.0887 / 5.0160, 10 modes 5.1064 / 5.0364, exact 5.1090 / 5.0407. |
| 3 | Β⟨q⟩_π/γ ≤ R0 ≤ βq_max/γ and the λ1 analogue, equality iff q constant, for every irreducible generator (Theorem 3.2); nonlinear demonstration (well-mixed 0.900, NGM 1.0757, attack 12.4%) | confirmed | Upper bound (1ᵀV^-1 = 1ᵀ/γ + entrywise comparison) and lower bound (Cohen convexity, tangent at v=0 with gradient π) are correct; the strictness argument (convex + analytic ⇒ g≡0 ⇒ v≤0 ⇒ v=0) is valid. 3000 random instances: max relative violation 9e-16 (lower), none (upper). An adversarial Nelder–Mead search over non-reversible rates and q initially reported −6e-2 and −3e-8; both were traced to its own numerics (stationary vector taken from an ill-conditioned eigen-solve; rates spanning 1e-12…1e18 so the float-built generator no longer had zero column sums). With a robust solver and rates clipped to e^±12 the minimum of R0/lower−1 over 120 optimised instances is −7e-13, i.e. zero. Nonlinear ODE re-run: β=0.10568, R0=1.0757, attack 0.1244, 97.25% of infections in the workspace, homogenised room attack 1.0e-5 — all identical. |
| 4 | R0→β⟨q⟩_π/γ (D→∞), →βq_max/γ (D→0), λ1(D) convex non-increasing, R0(D) non-increasing and strictly decreasing for non-constant q, any irreducible generator (Theorem 3.3) | confirmed | Perspective-function argument and the strictness reduction to Theorem 3.2 re-derived; both are sound. 400 instances × 49 log-spaced D: no increase (max relative step −1.8e-12), 0 strict-decrease failures, λ1 convexity/monotonicity violation none. Adversarial search for R0(D2)>R0(D1) with non-reversible generators: float64 'violations' up to 5e-4 appeared only at extreme rates/D (cond(V)~1e12) and vanished when the same instances were re-evaluated in 60-digit arithmetic with an exactly conservative generator (max increase 8e-59 and 0.0 over 80 optimised instances; verify/v1d_adv3_exact.log). Prior art confirmed by a search: Gao & Dong 2020 (PAMS 148:1709–1722, 'strictly decreasing and strictly convex… solves a conjecture of Allen et al.'), Chen–Shi–Shuai–Wu 2020 (JMB 80:2327–2361), Altenberg 2012 (PNAS 109:3705–3710). |
| 5 | Asymptotic expansions R0 ≈ β⟨q⟩/γ + βC/(⟨q⟩D) and (βq_max/γ)(1 − Dμ_Z/γ) (Proposition 3.4) | confirmed | Second-order and degenerate first-order perturbation derivations re-done by hand (group inverse, constants drop out, δ = βC/(⟨q⟩D)). Own check with an explicitly constructed group inverse on 300 instances: error ratio per decade of D median 99.9 (large D) and 99.5 (small D); C<0 never; symmetric closed form agrees to 2e-14. Office O1/O2 coefficients reproduced: 4.6429 + 0.4418/D0, 5.3571(1 − 0.0646 D0); 4.6429 + 0.2455/D0, 5.3571(1 − 0.0981 D0). Proposition 3.4 is stated with O/o remainders and proved by perturbation theory; no explicit remainder constants are given. |
| 6 | Zone lower bound R0 ≥ β·min_Z q/(γ+μ_Z) and block closed forms (Proposition 3.5) | confirmed | Schur-complement proof is correct (S ≤ V_ZZ entrywise, both M-matrices, so (V^-1)_ZZ ≥ V_ZZ^-1). 1800 random zones over all generator types: no violation (closest −1.8e-3 relative). Dirichlet block eigenvalue 2w[2−cos(π/(ℓx+1))−cos(π/(ℓy+1))] and the one-wall variant 2w[1−cos(π/(2ℓ+1))] verified to 1e-14 for 7×5, 10×2, 3×12. Office O1 meeting-room bound 5.0322 (μ_Z = 0.00904/day) reproduced. |
| 7 | Edge-wise monotonicity for symmetric rates; no monotonicity under the departure-site rule; corridor-D numbers (Theorem 3.6) | confirmed | One-line Rayleigh proof is valid. 1200 single-edge increases: max change +1.4e-14. Departure rule, doubling D at one cell: R0 up 191 / down 109 of 300 (the first analysis reports 201/98 for a ×3 change). Corridor-D table reproduced digit for digit: O1 symmetric 5.1090/5.2188/5.3252, O1 departure 5.1602/5.1235/5.0892, two-zone aisles departure 5.1427/4.1796/2.9401. |
| 8 | No finite-size correction with reflecting walls and uniform q (Theorem 4.1) | confirmed | -- |
| 9 | Leaky room: R0_room = βq/(γ+μ_κ), 0<μ_κ≤⟨κ⟩_π, monotone, absorbing limit, ν-weighted bounds, corridor with absorbing ends, room ratios and per-visit numbers (Theorem 4.2) | confirmed | All four parts of the proof check out (including the limit-point argument for k→∞). 400 instances: uniform-q identity 4.9e-14; μ ≤ ⟨κ⟩_π, monotonicity and both heterogeneous bounds never violated; absorbing limit reached. Corridor: μ = 2w[1−cos(π/(L+1))] exact for L=3…100 (note μ/(π²D/L²) = 0.53 at L=3, 0.86 at L=13, 0.98 at L=100: the continuum formula is reached only in the continuum limit). Room ratios reproduced: one door 0.7767, 0.9887, 0.9955, 0.9990; open side 0.5857, 0.9539, 0.9795, 0.9947. Per-visit values 0.1911, 0.0112, 0.0311, 0.0173 are arithmetically right; they treat every departure as permanent (see the further points). |
| 10 | One-door square room: μ_abs ≈ πD/(L² ln(L/c)), c ≈ 0.82–0.91; series formula for finite door rate | confirmed | Own sparse eigen-solves for L = 8…128 give L²μ_abs = 1.3808, 1.0762, 0.8755, 0.7363, 0.6347, c = 0.822, 0.864, 0.885, 0.898, 0.907 and series errors −1.51…−0.28% — identical. c is still drifting with L, so this is a fit, not an asymptotic constant, and the article states it as a numerical fit. Schuss–Singer–Holcman 2007 (PNAS 104:16098–16103) exists. |
| 11 | Hotspot maps u, w, z, elasticity e_r = z_r w_r/(z·w), Σe=1, reversible case e ∝ q x²/π (Theorem 4.4) | confirmed | Derivations checked (∂K/∂q_r, w ∝ Qx, z = Π^-1 x under detailed balance). 240 instances: \|Σe−1\| ≤ 4e-16; finite-difference elasticity error 3e-9; w ∝ Qx to 3e-15; reversible formula to 1.6e-14; e^{Jt}I0 → u to 1.3e-14. Office O1 elasticity shares W 2.21%, C 0.54%, M 97.23%, K 0.02% reproduced. |
| 12 | Transient caveat: at D0=1 the Perron map is not reached before saturation from a far seed | confirmed | Re-done with a separately written method (spectral solution of the linearised system): D0=1 gap 0.0561, far-desk seed fitted rate 0.5485, TV to u 0.981 at prevalence 1e-4 and 0.978 at 1e-2; meeting-room seed rate 0.6024, TV 0.002; D0=100 gap 0.580, rate 0.5129, TV ≈ 0. Matches the nonlinear-ODE values reported. |
| 13 | Recomputed two-zone and full-office R0 (5.090/4.371/4.360/4.381; bounds [4.258, 5.357] and [4.643, 5.357]; O1 5.109, O2 5.041; λ1 ≈ 0.60/day) | confirmed | verify/v2_examples.py (own generator and eigen code; only the zone arrays taken from scripts/layouts.py, cell counts re-checked: W156 C39 M26 K13 B26; two-zone 180/80, scatter 182/78) reproduces every number: 5.0895, 4.3709, 4.3598, 4.3805 (symmetric); 5.2320, 5.1427, 5.1420, 5.1484 (departure); O1 5.1090 / λ1 0.6021, O2 5.0407 / 0.5994; departure 5.1602, 5.1582; ⟨q⟩ = 1.3000 (area-weighted) and 1.3787 (π ∝ 1/D). Layouts are the first analysis's own constructions, as disclosed. |
| 14 | Partitions change R0 by +0.05…+0.27%, obstacles +0.18%; masks+ventilation 0.77 | confirmed | Separate random draws of the check on layout O1 (fewer samples): walls on 10/30/50% of desk–desk edges give 5.1109 / 5.1157 / 5.1226 (+0.04/+0.13/+0.27%), maximal walls (125 edges) 5.1220 (+0.25%), 39 desks → obstacles 5.1184 ± 0.0027 (+0.18%); never below baseline, as Theorem 3.6 requires. q_W→0.5q_W gives 5.0905 at D0=1 and 3.0312 (from 4.6473) at D0=100. 2.56 and 0.77 are exact consequences of R0 being linear in β and q. |
| 15 | Individual-based rule βq·n_I/n_tot is sub-critical; encounter factor c = 1 − [1−(1−π)^N]/(Nπ); simulated offspring numbers | partially confirmed | E[k/(k+1)] formula is correct; c = 0.1848, 0.1685, 0.1142 recomputed. Own simulation (verify/v3_dynamics.py, v4_ibm_big.py), local rule, N=100: D0=1 0.754 ± 0.022 (the first analysis's value 0.761), D0=0.1 0.545 ± 0.017 (the first analysis's value 0.564), mean-occupancy rule 4.44 ± 0.12 (the first analysis's value 4.50). At D0=10 the check's first run gave 0.769 ± 0.023 and a 9000-replicate run 0.825 ± 0.013 (integrated hazard 0.831 ± 0.009) against the formula 0.843: consistent within about 5 %; the room relaxation time at D0=10 is ≈4 days, so D0=10 is not yet the annealed limit. N=400: the check's value 2.57 ± 0.09, the first analysis's value 2.33 ± 0.06, formula 2.45. Qualitative conclusion (sub-critical by a factor ≈5–6 at 100 people) is solid; the quantitative match of the annealed formula is only to ≈±5%. |
| 16 | At realistic mobility the closed room is well mixed (+10.0% at D0=1, +0.095% at D0=100, +0.004% at 2.4e3) | confirmed | Reproduced for O1: 5.1090 (+10.04%), 4.6473 (+0.095%), 4.64304 (+0.0039%), 4.64287 at 4.87e4. O2 gives +0.053% at D0=100. Meeting-room exit time 1/μ_M = 110.6 days reproduced (layout-specific: one door cell). |
| 17 | Targeted ventilation: elasticity-guided placement beats random by 0.26–0.39 at D0=1; one-shot ranking worse than random from 20% | partially confirmed | Layout O1: adaptive and one-shot values reproduced exactly (4.7499, 4.5745, 4.2688, 3.5624; one-shot 4.7107 at 20%, 4.5195 at 50%). The 'highest q first' rule depends on a random tie-break among equal-q desk cells: one tie-break gives 4.484 and 3.583, another 4.503 and 3.629, and at 50% the adaptive value at the nominal area (3.593) is not better than the q-ranked 3.583 of the first tie-break (the article gives means over 30 random tie-breaks). Elasticity-guided placement is better than random; it is better than treating the highest-q cells first at 10–20% of the floor but not at 50%. |

Further points:

* (1) Per-visit reproduction numbers (office 0.19, supermarket 0.011, classroom 0.031, metro 0.017) treat every departure as permanent.
* (2) All layout-specific numbers of the check (meeting room as hotspot with 97% of the elasticity, 110-day exit time, +0.05…+0.27% for partitions) were computed on layout O1 with a one-cell door; layout O2 gives different secondary numbers (gap 0.0185 vs 0.0561, +0.053% vs +0.095% at D0=100). The article's values come from the office scene (see the note above).
* (3) Numerical hygiene for anyone reusing the scripts: with rates or D spanning more than ~10 orders of magnitude, float64 generators lose exact conservation and produce spurious violations of Theorems 3.2–3.3 of order 1e-5…1e-2 (the check hit this in its own adversarial tests); the scripts of the article stay inside safe ranges.

Summary:

Every proof was derived again and no theorem is refuted. Theorems 3.1–3.3, 3.6, 4.1, 4.2 and 4.4 and
Proposition 3.5 are proved in full; they rest on standard M-matrix facts and on Cohen's convexity theorem, which
are cited rather than proved again. Proposition 3.4 is stated with O/o remainders and proved by perturbation
theory, without explicit remainder constants; the one-door formula is a numerical fit. The optimiser-driven
searches for counter-examples first showed apparent violations (−6e-2 for the lower bound, +5e-4 for the
monotonicity in D) that came from the numerics of the search; with a robust stationary solver and 60-digit
arithmetic on an exactly conservative generator the worst cases are −7e-13 and 8e-59. The worked examples
reproduce to the printed digits. Two items are partially confirmed: the annealed encounter-factor formula agrees
with simulation to about 5 %, and the comparison of targeted ventilation with "highest q first" depends on a
random tie-break. Theorems 3.1–3.3 are known results for patch models and are credited to that literature in
the article.

### 4.2 Simulations

*How it was checked.* A second exact simulator was written separately (`code/bsc_sim/verify/vsim.c`),
together with separate solvers for the next-generation matrix and for the pair-level equations. Stored
results were re-run from their seeds, and the key statistics were simulated again with the second simulator.

*What was found.* The deterministic values and every statistic that was simulated again agree within Monte
Carlo error, and the stored results are reproducible bit for bit. One statement was refuted: the lower bound
of the reproduction number does not hold for contact kernels other than the same-cell kernel (the article
states it for the same-cell kernel and gives the counter-example in Proposition 3.7). Three items are partially
confirmed: corridors are not the hotspots of infection (the workstation cells that border them are); barriers
delay the peak, by about 0.8 day; and the speed-up of the spectral method over Monte Carlo depends on the
implementation. The effect of heterogeneity on the counted number of secondary cases at high
occupancy holds only for an index case placed with the Perron vector.

*Items: 16* (12 confirmed, 3 partially confirmed, 1 refuted).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | Simulator is correct (agrees with pure-Python reference, exact well-mixed final-size law, NGM via exposure identity, lattice ODE at large occupancy, schedule tests) | confirmed | The check wrote a separate brute-force exact simulator (verify/vsim.c: different event selection, RNG, no incremental bookkeeping) and separate theory code (verify/vlib.py). The check's values vs exact final-size law N=40: mean 27.808 exact vs 27.79-27.87 over 5x200k runs. Exposure identity, office uniform index: pooled 4.638 +- 0.009 vs 4.6429 (400k runs); Perron index 5.133 +- 0.025 vs 5.107. N=1000, D0=10 mean curve peak 0.431 vs ODE 0.454 (first analysis: 0.428/0.454). The check's simulator and the first analysis's C core agree in ~25 configurations. The first analysis's stored numbers are bitwise reproducible from its seeds (03_scenes office/metro, 01 tests B and E; test F re-run with a new seed gives max diff 0.007). |
| 2 | Pair-level theory: exact pair-lattice solve, closed form R1 = (beta q/gamma)/(1+(beta q/rho) g0), elliptic-integral g0 on the infinite lattice; rises with D and room size | confirmed | Derivation re-done (Feynman-Kac for the killed pair walk; translation invariance gives p0 = h g0/(1+h g0); sum over y gives the closed form). Closed form equals the check's exact pair solve on a torus to 8 digits (L=6: 2.44337619 both). g0 = (2/pi)K(k)/(gamma+8D) with modulus k = 8D/(gamma+8D) matches the k-integral: 0.241939 (D=1), R1_inf = 2.710 / 2.303 / 3.421. The check's exact solve vs the check's counts (30k-200k runs): L=4 1.8821 vs 1.8825 +- 0.0037; L=10 2.4652 vs 2.468 +- 0.009; L=20 2.596 vs 2.612 +- 0.012; D=10 L=10 3.141 vs 3.159 +- 0.011; D=0.51 L=20 2.199 vs 2.204 +- 0.010; L=40 closed 2.651 vs 2.629 +- 0.017. Small discrepancies: the reflecting-wall closed form is off by 0.7-0.8 % at D=0.51, L=4-5, and two of 24 points lie outside their CI against the exact solve. |
| 3 | Mean-field R0 overstates discrete secondary cases ~2x at the default parameters; converge with mobility; opposite mobility dependence; gen-2/gen-1 = 1.63 | confirmed | Office, the check's code: R0 = 5.348 / 5.107 / 4.676 / 4.646 at D0 = 0.03 / 1 / 10 / 100. rho(K1) from a separately written forward formulation (occupation density of the surviving pair, 234 CG solves): 0.7095 / 2.3212 / 3.8025 / 4.3432 (first analysis 0.710 / 2.321 / 3.803 / 4.343). Counted with K1-Perron index: 0.708 +- 0.005 / 2.335 +- 0.012 / 3.813 +- 0.020 / 4.317 +- 0.052. Gen-2/gen-1: 0.356 / 1.62 / 3.07 / 3.49. The check's uniform-index count at D0=10 is 3.787 +- 0.020 (theory 3.793), supporting the first analysis's explanation that its 3.697 was a chance-low batch. |
| 4 | Scene results at D0 = 1 with 1 m kernel; D0 = 10 near mean field | confirmed | The check's simulator (3000-4000 runs): office R0 5.107, pair R1 2.266, counted 2.235 +- 0.039, P(major) 0.498 (0.482-0.513), attack\|major 44.8 %, peak 12.9 % (0.5-day grid) on day 23.5; supermarket 4.394 / 1.990 / 1.973 +- 0.011, 0.429, 42.0 %, 10.0 %, day 30.9; classroom 5.866 / 2.581 / 2.576 +- 0.013, 0.603, 39.3 %, 14.1 %, day 16.1; metro 9.241 / 5.584 / 5.559 +- 0.023, 0.864 (0.851-0.876), 94.7 %, 25.1 %, day 22.8. Office D0=10: 0.734 (0.714-0.752), 96.3 %, 39.0 %, day 16.9; branching 0.783; ODE 0.989/0.451/11.5 all recomputed. Metro P(major) and peak day and office D0=10 P(major) differ from the first analysis by 1.6-2.3 sigma, i.e. sampling noise. |
| 5 | Contact kernel matters: strict same-cell rule gives no metro outbreak at D0 = 1; classroom 0.276 vs 0.614 | confirmed | The check's simulator: metro same cell 0 of 3000 runs reach 10 % (counted 1.535 +- 0.030, pair 1.514, gen-2/gen-1 0.61); classroom same cell P(major) 0.270 (0.257-0.284) vs 0.603 with the 1 m kernel. Same-cell R0: metro 9.836, classroom 6.182. Note that rho(K1) = 2.24 in the metro same-cell case while no outbreak occurs: R1 is not a threshold quantity. |
| 6 | Hotspots: mean-field map fails in the office at D0 = 1, pair-level map works; where in the office infections concentrate | partially confirmed | Correlations reproduce with the check's simulator and the check's K1 (40,000 outbreaks, gmax=2): mean field 0.105 / -0.234 / -0.366, pair level 0.878 / 0.824 / 0.726; first generation alone r = 0.963 vs 0.972 ceiling; zone shares corridor 15.4 % observed, 15.3 % pair, 10.3 % mean field. The corridor is 16.7 % of the floor, so it is below average per cell (relative rate 0.96 in generation 1, 1.09 in gen 2, 1.15 in gen 3). The hottest cells are workstation cells bordering a corridor (1.15 / 1.26 / 1.36 vs 0.94 / 0.88 / 0.87 for interior workstations); the top 15 cells are all workstation cells (one meeting-room cell in gen 1, one corridor cell in gen 3). |
| 7 | Barriers with identical starts: mean-field R0 unchanged, outcomes fall modestly; timing of the peak | partially confirmed | Confirmed: R0 = 4.2857 at every density; pair R1 falls (uniform index 2.92 -> 2.83 -> 2.73); with its own 40 nested layouts x 1000 runs, identical starts: mean attack 0.612 -> 0.564 -> 0.376 (differences -0.047 +- 0.004 and -0.235 +- 0.011, layout-level SE), peak of major outbreaks 27.0 % -> 24.6 % -> 17.0 %, P(major) 0.702 -> 0.683 -> 0.627. The peak is delayed: peak day moves 21.46 -> 22.28, i.e. +0.82 +- 0.10 d at 10 % barriers and +0.71 +- 0.30 d at 30 %. The first analysis's own paired runs give +0.67 +- 0.31 d at 10 %; an across-layout CI of +-1.1 d cannot resolve it. The delay is about 0.8 day. A paired CI that treats runs as independent (+-0.017) ignores between-layout variance; the check's first 8 layouts gave -0.025 +- 0.007. |
| 8 | Heterogeneity in q raises mean-field R0; simulation confirms rise at N = 1000; at 100 people realised R and attack rate fall | confirmed | Two-zone case recomputed: R0/(beta/gamma) = 1.498 / 1.485 / 1.392 / 1.240 / 1.204 at D0 = 0.01 / 0.1 / 1 / 10 / 100. Lower bound proven for the same-cell kernel by a Rayleigh quotient with the constant vector. N=1000, D0=1, the check's K1 and the check's simulator: rho(K1) 3.967 -> 4.335 -> 4.817, counted (Perron index) 3.92 +- 0.08 -> 4.33 +- 0.08 -> 4.78 +- 0.09. N=100: rho(K1) 2.747 -> 2.115, counted 2.774 +- 0.021 -> 2.113 +- 0.017, mean attack 0.341 -> 0.239, P(major) 0.527 -> 0.476. Controls: uniform D 2.777 -> 2.822; hot corridor c=-2: R0 9.141, rho(K1) 4.972, counted 4.973. Caveat: at N=1000 the rise appears only for a Perron-weighted index; for a uniformly placed index the count falls (3.96 -> 3.86 -> 3.82) and outbreak probability and attack rate fall. The fall of -24 % in rho(K1) at 100 people is -15 % (uniform index) to -28 % (gen ratio) with other measures. |
| 9 | Under the frequency-dependent law mean-field R0 is independent of occupancy; capacity limits act through discreteness unless density dependence is assumed | confirmed | R0 = 5.107 for N = 50, 100, 400 by construction (rho_bar normalisation). The check's pair solver and simulator, office D0=1: FD N=50 R1 1.521, counted 1.556 +- 0.028, P(major) 0.383, peak day 13.4; N=100 2.266 / 2.235 +- 0.038; N=400 3.650 / 3.54 +- 0.06, P 0.729; DD N=50 R0 2.528, R1 1.122, counted 1.109 +- 0.022; DD N=400 R0 20.58, R1 9.133, counted 8.96 +- 0.13. Elasticities are arithmetic on these values (ln(2.72/1.95)/ln 2 = 0.48). |
| 10 | Interventions: halving q halves R0; combination gives R < 1 only under density dependence at D0 = 1 | confirmed | The check's code: ventilation R0 2.554, R1 1.509, counted 1.491 +- 0.029, P(major) 0.303; masks 1.532 / 1.049 / 0.153; capacity 50 % peak day 13.4 vs 23.5 baseline. Combination with three of its own partition layouts (39 workstation cells): R0 2.55 (FD) / 1.52 (DD), R1 1.10 / 0.816, counted 1.11-1.13 / 0.79-0.84. Classroom combination R1 1.20 / 0.85; office D0=10 1.78-1.80 / 1.16; classroom D0=10 1.69-1.70 / 1.09-1.10. R1 is not a threshold quantity, so R1 < 1 is not a criterion of control (see the further points). |
| 11 | A daily schedule behaves like its time average; realistic dwell times shrink the closed-room reproduction numbers | confirmed | Test in the check's simulator (office, half-day hop x1.8, m=1.3 / half-day hop x0.2, m=0.7 vs constants, 6000 runs each): P(major) 0.491 vs 0.498, attack 0.448 vs 0.446, peak 12.83 % vs 12.80 %, day 23.7 vs 23.3. The first analysis's metro arms re-run with a new seed: varying vs mean P 0.844 vs 0.841, peak 22.8 % vs 23.0 %, day 25.1 vs 25.0. Duty cycle: R0 1.754, pair R1 0.889 recomputed; counted with explicit schedule 0.928 +- 0.003 pooled over both simulators (the first analysis's 0.918 is reproducible from seed 1100 but is a low batch). P(major) with schedule 0.065. Per visit: office 0.1750 +- 0.0007 (mean field 0.2117), supermarket 0.0085, lecture 0.0273, metro ride 0.0175. D0=10 duty-cycle count (1.408) not re-simulated. |
| 12 | Fair benchmark: speed-up of the spectral R0 over Monte Carlo at equal accuracy; mean-field route biased for whole-epidemic outputs | partially confirmed | Re-running the first analysis's script reproduces its factors within ~20 % (9.2x and ~113x at 1 %; pair solve 22x and 97x; biases +0.283 / +0.030 and +0.546 / +0.026). The method is sound (same quantity, cost scaled as accuracy^-2), but the Monte Carlo side is vectorised numpy. With the check's compiled C walker estimator the factors are 3.1x (D0=1) and 86x (D0=10) at 1 %, 307x and 8,590x at 0.1 % (spectral 4.1 ms and 0.95 ms). At 1 % accuracy the speed-up is therefore a few-fold to about 10^2, and the factors depend on the implementation by 2-3x. |
| 13 | Lower bound β⟨q⟩/γ ≤ R0 of Theorem 3.2 for contact kernels other than the same-cell kernel (the article states the bound for the same-cell kernel and gives the counter-example below in Proposition 3.7(iv)) | refuted | The upper bound holds (column sums). The lower bound is proven for the same-cell kernel (Theorem 3.2). Counter-example with the kernel construction of the article: 3 cells in a row, radius-1 kernel, q = (1, 0, 1), D -> 0 gives R0 = 0.500 beta/gamma < <q> beta/gamma = 0.667; q = (1, 0.2, 1) gives 0.620 < 0.733. In the four scenes the bound happens to hold (classroom 5.87-6.54, metro 9.24-9.26 above 4.71 and 9.06). |
| 14 | Scene summary values (area-mean q 1.30 / 0.91 / 1.32 / 2.54; metro mean D/D0 0.028) and scene data | confirmed | Recomputed from scenes/*.json with its own parser: mean q 1.300, 0.914, 1.318, 2.536; metro mean D 0.0282 (zones 0.0044 / 0.0047 / 0.109 / 0.011); office zone counts 156/39/26/13/26 reproduce the 60/15/10/5/10 % fractions. Floor plans of the other three scenes are design choices, as the first analysis states. |
| 15 | Results: 22 figures each in _en and _zh, PDF + PNG, from real runs | confirmed | 88 files present (22 x 2 languages x 2 formats). Inspected fig5_3 (zh), fig6_3 (en), fig5_2 (zh), fig4_meanfield_validity_map (en), fig6_1 (zh): CJK text renders, curves and values match the JSON data. |
| 16 | Spatial branching and lattice ODE predictions quoted in the scene tables | confirmed | Its own fixed-point and LSODA implementations give identical values: branching P(major) 0.781 / 0.667 / 0.767 / 0.888 at D0=1 and 0.783 / 0.685 / 0.777 / 0.888 at D0=10; ODE attack/peak/day 0.989/0.451/11.5, 0.936/0.289/17.0, 0.984/0.426/11.0, 1.000/0.641/7.0. Galton-Watson with the simulated offspring law: office 0.612 (mean 2.267, variance 6.09, P0 0.269), metro 0.881. |

Further points:

* (1) The lower bound beta<q>/gamma <= R0 is false in general once a contact kernel P != I is used (3-cell counter-example: 0.50 vs 0.667).
* (2) R1 and rho(K1) are not threshold quantities. Metro same-cell at D0=1: rho(K1) = 2.24, R1 = 1.51, yet gen-2/gen-1 = 0.61 and 0 of 3000 outbreaks; classroom same-cell: rho(K1) = 2.27 with gen ratio 0.90. No quantity computed here predicts the actual threshold of the discrete system.
* (3) Corridors are not the hotspots. The corridor receives 15.4 % of first-generation infections on 16.7 % of the floor; the hottest cells are workstation cells bordering corridors (relative rate 1.15-1.36 against 0.87-0.94 for interior workstations).
* (4) Benchmark factors depend on the Monte Carlo implementation: with a compiled walker estimator the equal-accuracy speed-up is 3x (D0=1) and 86x (D0=10) at 1 %.
* (5) The evidence that heterogeneity raises the counted number of secondary cases at N=1000 rests on a Perron-weighted index (and the gen ratio).
* (6) Physical regime: every spatial effect of the scenes lives at D0 = 1 cell^2/day (a person diffuses about 1.5 m per day) and a closed room occupied 24 h a day. With an 8-hour day the office R is already 0.9-1.4.
* (7) Barrier design details: a paired CI (+-0.017) that treats runs as independent ignores layout-to-layout variance; starting everyone on cells walkable in the densest arm is not the stationary law of the sparser arms (the check measured attack 0.600 vs 0.589 at zero barriers, within noise, so the effect is small).
* (8) Metro time-dependence keeps N = 310 while the 'density' signal varies, so only mobility and the multiplier change; it is not a variable-occupancy experiment.
* (9) Small discrepancies: closed-form accuracy for reflecting walls is 0.7-0.8 % at D = 0.51 in 4x4 and 5x5 rooms; two of 24 finite-size points lie outside their CI against the exact solve; the duty-cycle count 0.918 of seed 1100 is a low batch (0.928 +- 0.003 pooled).

Summary:

All recomputed deterministic values and every statistic that was simulated again with the second simulator agree
within Monte Carlo error, and the stored results that were re-run are reproducible bit for bit from their seeds.
The lower bound of Theorem 3.2 needs the same-cell kernel (counter-example 0.50 against 0.667 for a kernel over
neighbouring cells). The pair-level numbers R1 and ρ(K1) are not thresholds: the metro with the same-cell kernel
has ρ(K1) = 2.24 and no outbreak in 3000 runs. Corridors are not the hotspots; the workstation cells bordering
them are. Barriers delay the peak by about 0.8 day. The speed-up of the spectral method depends on the Monte
Carlo implementation (3× and 86× at 1 % accuracy with a compiled estimator). The rise of the counted number of
secondary cases with heterogeneity at N = 1000 holds for a Perron-weighted index only. The spatial effects of the
scenes are large only at low mobility (D0 = 1 cell²/day) in a room occupied around the clock. Not re-simulated
with the second simulator: the scene rows at D0 = 100, and test F with N = 5000 (re-run with the same code).

### 4.3 References

*How it was checked.* Every entry of `refs.bib` was generated from, or compared with, the Crossref (or DataCite)
record of its DOI; the records are stored in `data/plan_refs_crossref.json`, `data/refs_search.json` and
`data/v2_refs_crossref.json`. Books without a DOI were checked against publisher or library catalogue records.
The theorems of the patch-model literature quoted in the article were checked against the published versions of
Gao and Dong (2020), Gao and Ruan (2011) and Tien et al. (2015); for Allen et al. (2007) and Gao (2019) only the
abstracts were read, and the article relies on the account of Gao and Dong (2020) for them.

*What was found.* The references exist as given, and the theorems are as quoted.

### 4.4 Outbreak dataset

*How it was checked.* The sources were fetched again, the counts typed again, and every derived number
(attack rates and intervals, risk ratios, the join-count test of the call centre) was recomputed with separate code. The seat map of the call centre was read again desk by desk from
enlarged crops of the published figure.

*What was found.* The numbers and the sources stand. Three of the eight records graded A do not meet the
stated criteria; the "three regimes" of proximity gradient are a reading of point estimates; and the protocol
first proposed for a test would have let a weak model pass. In addition, the five secondary cases of the
Guangzhou restaurant are an upper bound. That no outbreak in a metro
carriage with a known number of exposed riders exists could not be verified.

*Items: 4* (3 partially confirmed, 1 unverifiable).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | A validation dataset of 27 outbreak records from 31 sources was assembled (all DOIs registered); 8 records are grade A (C1, T1, T2, T4, T5, R1, R2, H1). | partially confirmed | Counts and sources hold: 27 records, 31 distinct source keys, grades 8/9/8/2, 52 CSV columns (51 in this repository, without a bookkeeping column). Its own registry query found all 31 DOIs registered with Crossref titles matching. The check re-fetched 21 full texts from Europe PMC plus the Miller accepted-manuscript PDF and found every quoted count in them. The three sources not served as raw text (Shen, Hu, Hamner) were re-extracted separately with the same values. Two internal checks support them: risk ratios recomputed from Shen's zone counts reproduce the published 1.6 (0.8-3.2) and 1.8 (0.9-3.3), and cell counts recovered from Hu's printed intervals sum to 230/72,692 against the reported 234/72,093. Only 5 of the 8 grade-A records meet the stated grade criteria: C1 has no exposure duration in hours, T4 pools 2,334 index cases, and H1 has no published seat positions. |
| 2 | Documented outbreaks fall into near-field, near-well-mixed and airflow-zoned patterns (RRs: flight 7.3, classroom 2.8, bus 1.6, coach 1.0; restaurant 5/11 vs 0/68; train adjacent 3.53% vs 0.05-0.25%), and per-hour risk varies 34-fold between single-index events (0.024 to 0.81 per hour, GM 0.13, GSD 3.1). | partially confirmed | Every number reproduces. The three-regime classification is not statistically established. The four single-event risk ratios are mutually compatible (Cochran Q = 3.77, 3 d.f., p = 0.29; pooled 2.0, 1.2-3.3). The bus and coach intervals reach 3.2 and 4.1, and Ou et al. report a side-to-side contrast on the coach (6/24 vs 1/20). The choir's 'well mixed' label is the authors' description with no seat data. The train gradient is subject to companion confounding that Hu et al. acknowledge. The restaurant zone effect depends on counting all five cases as restaurant-acquired. The 34-fold spread rests on 20 probable choir cases; with confirmed cases only it is 20-fold and GSD 2.6. |
| 3 | No documented metro-carriage outbreak with an exposed denominator was found; supermarket evidence is one record with denominators and no geometry (Liaocheng: staff 11/120, customers 0/8,224). | unverifiable | The Liaocheng numbers are confirmed in the Tian et al. full text, and Furuse et al. report 61 clusters with one transport-related (an airplane). The check's two additional searches also returned no metro-carriage outbreak with a denominator. Absence cannot be proven by search, as the result itself says. |
| 4 | Mapping and validation protocol written; no model-versus-data comparison was run. | partially confirmed | The mapping table and protocol exist as described and no model run is present. (1) The distance kernel is calibrated on the train matrix without addressing companion confounding of the adjacent-seat cell. (2) The shape criterion is passed by any constant near/far ratio between 1.16 and 3.19 for all four single events. (3) The level criterion's predictive distribution for infectiousness is undefined from the calibration set (one event-specific value); with the observed GSD of 3.1 the 90% band spans a factor 6.4 either way, and all nine observed hazards fall inside it, so a constant-hazard model with no spatial structure would score full coverage. (4) Restaurant calibration rests on 2 to 5 cases. Only condition 4 of the decision rule (beat the well-mixed null by two standard errors) has real power. |

Further points:

* (1) The Guangzhou restaurant counts are an upper bound. Lu et al. 2020 say at least one member of each of families B and C was infected at the lunch and the others may be within-family (onsets 27 Jan to 5 Feb for a 24 Jan lunch). The dataset uses 5/79, 5/16 and 5/11 throughout.
* (2) Three of the eight grade-A records do not meet the stated grade-A criteria: C1 Marin classroom (no exposure duration in hours, no room size), T4 trains (2,334 pooled index cases, no printed cell denominators), H1 Skagit choir (no seat or zone positions). Five records satisfy the definition.
* (3) The train cohort is the designated calibration set for the distance kernel, but Hu et al. state that family members and friends may inflate the on-train risk, most of all in the adjacent-seat cell; the first outbreak test therefore excludes that cell from the calibration.
* (4) The protocol's decision rule has low power. Any constant near/far ratio in 1.16-3.19 passes all four single-event shape tests. The level test's infectiousness distribution is undefined from the calibration set, and with the observed spread (GSD 3.1, computed from the held-out events themselves) a spatially blind constant-hazard model would achieve full coverage.
* (5) The 'three regimes' of proximity gradient are a reading of point estimates, not a statistical finding between single events (Cochran Q p = 0.29). The bus and coach are compatible with gradients as steep as the classroom's.
* (6) The share of employees in the north wing of the call centre assumes one employee per drawn desk in the two wings; about 66 further seats in the east rooms and offices were not digitised, so the share of the floor's employees seated in the north wing lies between 52 and 69%. No conclusion depends on it.
* (7) The wing contrast of the call centre applies to the 84 mapped case seats; with the 10 unmapped cases all in the south wing it would be 2.6-fold (CI 1.6-4.1). The onset growth-rate interval that ignores over-dispersion (2.07) is too narrow; with it the interval is 0.18-0.64.
* (8) The 34-fold per-hour hazard spread depends on 20 probable (not laboratory-confirmed) choir cases. With confirmed cases only it is 20-fold and the maximum becomes a 2/13 event.
* (9) Minor: Luo 2020 and Ou 2022 disagree on where the coach cases sat (farthest about 4.5 m vs 9.46 m). The Jerusalem senior-wing denominator is 583 persons but 580 tested. One duration is rounded (1.37 h for 82 min).
* (10) Three sources (Shen zone counts, Hu Table 1, Hamner schedule and chair layout) could again be read only through a page-extraction tool. Europe PMC, NCBI efetch and BioC do not serve their full text. They pass internal consistency checks but were not compared with the PDFs.

Summary:

The dataset is numerically sound and its sources exist: the sources were fetched again, the counts typed again
and every derived number recomputed with separate code, and all agree; all 31 DOIs of the first dataset resolve
with matching titles. No model is fitted to the dataset in this part (the tests are in Sections 4.5 to 4.9).
Three of the eight grade-A records do not meet the grade criteria; the "three regimes" of proximity gradient are
a reading of point estimates; the restaurant counts are an upper bound; the train calibration set is subject to
companion confounding, which its source acknowledges; and the protocol first proposed with the dataset
(`code/bsc_outbreaks/PROTOCOL_PROPOSAL.md`) has discriminating power only in its comparison with the well-mixed
model.

### 4.5 First outbreak test

*How it was checked.* The calibration, the hold-out predictions, the power table, the ceiling of the
same-site model and the call-centre test were implemented again with separate code
(`code/bsc_validation/verify/`), the counts were typed again from the sources, and the freeze record was
compared with the hashes and times of the files.

*What was found.* Every headline number was reproduced. One reading was refuted: the protocol estimates the
call-centre probability from 40 posterior draws (0.031 and 0.017), and that estimate is Monte Carlo noise;
integrating the posterior by quadrature gives 0.007 and 0.013, below the limit 0.025, so M2 fails on the call
centre and its outcome is W4. The quadrature was added by the check, after the results. The classroom agreement
is not a prediction, the pooled test rests on one flight, and a constant near/far risk ratio of 2 to 3 without
any calibration passes the same rule with a higher score.

*Note on the items below.* The propagator of the first test, as frozen, returns round-off noise at very small
mixing coefficients, and the check's code shared it, so the numbers in the items below are those of the frozen
propagator. The article reports the values of a round-off-safe propagator: the classroom has p = 0.068 (frozen
propagator: 0.11); less than 0.1 % of its predictive probability comes from coefficients below 0.1 m²/h, a
regime that is numerical, not a property of the lattice; its median implied emission is 8.1 × 10⁷ quanta per
hour (frozen: 1.4 × 10⁷); and the pooled score is +2.83 (p = 2.5 × 10⁻⁴; frozen: +3.87). The pooled score without
the flight is −3.19 (frozen: −2.16); with the adjacent-seat cell included in the calibration it is −4.86, and with
the whole index row excluded +4.91 (frozen, as computed by the check: −4.0 and +6.0); and the constant ratios of 2
to 3 score +7.1 to +7.7 against +2.83 for M2 (frozen: +3.8). The outcome W4 is the same with both propagators.
The call-centre probabilities were computed with a propagator that shares the defect. The probability of a join
count as small as observed is 0 at every mixing coefficient below 5.6 m²/h, so the defect can act only through
the outbreak-size weights; with zero weight at coefficients up to 0.1 m²/h the pooled value would be 0.015, still
below the limit 0.025. For the same-site model the expected number of cases stays below the observed one in three
of four held-out events, but the observed totals are reached with probability 0.0002 (bus), 0.025 (flight) and
0.081 (classroom); under rule P4 of the post-hoc log, adopted after the results, only the bus total counts as
unattainable. The classroom draws of the same-site model take their mobility from the restaurant posterior, not
from the trains.

*Items: 14* (7 confirmed, 6 partially confirmed, 1 refuted).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | Pre-registration written before results and followed (timestamps, hashes, post-hoc log) | partially confirmed | All file metadata is consistent with the claim, but it is the only evidence (no version control). PREREGISTRATION.md was created and last modified at 05:04:36; PREREG_FREEZE.json was created at 05:12:03; the first calibration log starts at 05:13:30; 02_predictive.json (05:15:43) and the power addendum (05:16:17) precede 03_holdout_shape.json (05:19:04). All six frozen hashes match (recomputed on the unredacted originals; the released protocol files are redacted copies whose hashes differ, see `code/REDACTIONS.md`). The freeze file itself was rewritten at 05:29, so 'frozen_at' is self-reported. The hold-out code (predict.py, decision.py, 02/03 scripts) was written after calibration results were seen and is not hashed; predict.py was edited at 05:21:59, after unblinding. The check re-ran the current code without writing to data/ and it reproduces the stored predictive pmfs with zero difference, so that edit did not change the primary predictions. The protocol was followed with small logged deviations (exact enumeration instead of Monte Carlo for A1; 300 draws instead of 400). One pre-registered sensitivity analysis, the choir with 32 confirmed cases, was computed (1,034 quanta/h). The pre-registration is data-aware: strata, ventilation priors and the exclusion of the adjacent-seat cell were chosen knowing the outcomes. |
| 2 | Decision rule and its power, computed before unblinding | partially confirmed | The power table reproduces from its own predictive distributions (M2 judged, M2 true: P(A1) 0.998, P(A1 and no A2 failure) 0.877; M0 true: P(A1) 0.0019). The rule is strengthened only against the well-mixed null. A calibration-free model with one constant near/far risk ratio of 2 to 3 in every event passes A1 with a larger pooled score than M2 (+7.1 to +7.7 against +3.8, p = 4e-5) and has no A2 failure on the four single events. Passing A1 therefore shows that a proximity gradient exists in the data, not that the lattice kernel or its calibration adds anything. |
| 3 | Outbreak counts retyped from sources for the events used | confirmed | Zhejiang bus re-fetched from the JAMA Internal Medicine page and PMC: 14/33 (rows 5-11) vs 9/34, RR 1.6 (0.8-3.2); 11/23 vs 12/44; index in row 8, middle seat of the 3-seat side; 2 x 50 min. Train table re-fetched from the CID page: all 23 cell percentages and intervals match, as do 2,334 index cases, 72,093 contacts, 2.1 h (SD 1.8), 0.4 m row spacing and the companion caveat. The reconstructed counts reproduce the printed percentages to within 0.005 points. From the stored full texts: Hunan coach Table 1 (3/19 rear, 4/26 front; 6/24 vs 1/20 by side; 4.8 and 9.6 ACH; 60.42 and 21.69 m3), flight VN54 Table 2 (11/12 vs 1/8, RR 7.3; 28 business seats, 21 occupied; 0/35 premium economy), Marin classroom (8/10 vs 4/14 of 24 students, 22 tested), restaurant table sizes and overlap times. heldout_outcomes.json matches. Details not modelled: two coach passengers boarded 40 minutes late and four left early; VN54 carried four companion couples seated together in business class, three of them case-pairs. |
| 4 | Calibration fits on trains and restaurant | confirmed | Re-implemented with a different propagator method (numerical eigendecomposition of the generator, no import from valmod). Trains, 22 cells: M0 -99.704; M1 -73.116 at D_p 0.501 (90% 0.31-0.79); M2 -72.965 at D_air 25.1, kappa 13.12 (90% 9.7-43); deviances 126.1, 72.9, 72.6, p = 5.4e-17, 6.0e-8, 6.7e-8; fitted cell rates identical to Table V2. S-adj and S-row modes reproduce (10 and 100 m2/h). Restaurant, with its own random layouts: M0 -6.079; M1 -2.50 at 0.25; M2 -2.51 at 0.50; 90% interval 0.02-41 (first analysis: -2.47, -2.48; Monte Carlo layout difference). The calibrated model is rejected by its own calibration data. |
| 5 | M2 beats the well-mixed model out-of-sample and reproduces the flight and classroom gradients | partially confirmed | The numbers reproduce with its own event builders and seeds: pooled Delta +3.76 to +3.79, p = 1.5e-4 to 1.7e-4; flight p = 0.61 (RR 3.6); classroom p = 0.11-0.14 (RR 5.0-5.2); bus p = 0.38 (RR 2.4). Two qualifications apply. (1) The pooled pass rests on the flight alone: without it Delta is -2.16 and A1 fails; without the classroom it is +2.41, p = 0.0024. (2) The classroom agreement is produced by the width of the restaurant posterior, not by the calibrated kernel. At the posterior mode (0.5 m2/h) the observed 8/10 vs 4/14 is rejected at p = 0.0009, because the model puts all 10 front desks among the cases. The 63% of posterior mass between 0.1 and 5.6 m2/h contributes 5% of the predictive probability; the rest comes from D < 0.1 (a regime where the frozen propagator returns round-off noise; see the note above) and from the diffuse tail D > 5.6. In 80% of the posterior the fitted classroom emission exceeds 1e4 quanta/h (median 1.4e7, 95th percentile 1.6e16), which is not physical. The flight result is robust: p = 0.60 at the train posterior mode, and 0.20-1.0 when the ventilation prior is scaled by 1/3 or 3. |
| 6 | M2 fails on the Hunan coach | confirmed | Reproduced: model RR 10.8-11.1, predictive p = 0.006, M0 p = 1.0. The failure survives checks that were not in the first analysis's list: widening the train posterior for its over-dispersion (deviance/d.f. 3.6) gives p = 0.026; cutting the coach ventilation to a third gives p = 0.040. S-row (whole index row dropped) gives RR 2.6, p = 0.31, as reported. At the train posterior mode alone the coach is rejected at p = 0.001. |
| 7 | M1 (same-site transmission) fails | confirmed | The check wrote its own event-driven simulator of the per-cell rule in the limit of infinite transmission rate (infection at first co-location, seats reset per segment). Bus: mean total 9.1 at D_p = 0.5 (6.2 at 0.31, 12.2 at 0.79) against 23 observed; no replicate reached 23. Business cabin: 5.4 (5.0-7.1) against 12; P(total >= 12) = 0-0.08. Classroom: 1.1 at the room mode against 12. These match Table V8. Closed-form M1 reproduces too (bus RR 20, coach p = 1e-4, pooled -12.4 against the first analysis's -11.1). The classroom ceiling depends on the weakly identified room posterior: at D_p = 14 m2/h the ceiling is 23 of 24. The verdict is specific to events in which everyone is seated. |
| 8 | Seoul call centre, M2: probability of a join count as small as observed, estimated by the protocol from 40 posterior draws (0.031 pooled, 0.017 equal weight) | refuted | The observed statistic reproduces (91 case-case pairs vs 94.6 +/- 5.6, z = -0.65), as do M0 (P = 0.28) and M2's median z of about 11. The 40-draw estimate does not: it is Monte Carlo noise, three of the 40 draws having D > 100 m2/h, where the posterior expects 1.2. Integrating P(z <= observed \| D) over the same posterior on a 25-point grid gives 0.007 with equal weight and 0.013 with pooled weight. Both are below 0.025, so the observation is outside the central 95% interval under either convention. The check's separate 40-draw run gave 0.0008 and 0.0014, with 1 of 40 draws compatible. The value also depends on the arbitrary upper limit of the flat prior on log D (10^3: 0.004/0.007; 10^4: 0.007/0.013; 10^6: 0.013/0.024), because the restaurant likelihood is flat towards the well-mixed limit. The call centre is a clear A2 failure for M2. Rounding 137 desks onto 113 lattice sites is immaterial (checked on 2x and 3x finer lattices). |
| 9 | One mixing coefficient per setting class is rejected; restaurant fit conflicts with tracer data | confirmed | Recomputed for vehicles: preferred D_air alone is 0.06 (flight), 25 (trains), 63 (bus) and 6,300 m2/h (coach, grid edge); sum of maxima -76.76, common-D maximum -81.54 at 40 m2/h; likelihood ratio 9.56 on 3 d.f., p = 0.023. This analysis is post hoc, as logged. The measured tracer ratios in Table V20 match the stored Li et al. table; the check did not recompute the fitted-model column. |
| 10 | Levels not predicted; two level checks pass; no negative control fails | partially confirmed | The numbers reproduce: choir 2,732 and coach 36 quanta/h by hand; train cohort mean 1.34 (0.80-2.31); coach to minibus 1.18 expected cases; Utah classrooms at mask factor 0.35 mean 2.1, 90% interval 0-5, P(X >= 5) = 0.09. The passes are weak. The Utah check fails at mask factor 0.2 (P(X >= 5) = 0.016), inside the stated range 0.2-0.6, and the factor was chosen knowing the observed 5/728. A negative control fails only if more than half the predictive mass exceeds the 95% bound, which almost nothing can fail: premium economy passes although the 95th percentile of the prediction (18.6%) is more than twice the bound (8.2%). Table V10 does not list the emission M2 needs for the events where it is not rejected: flight median 1,600 quanta/h (435-11,800), classroom median 1.4e7. The 36-2,700 spread uses the 20 probable choir cases. |
| 11 | Decision outcome (M1: W4; M2: W2 with the protocol's 40-draw estimate of the call-centre probability, W4 with that probability computed by quadrature) | partially confirmed | M1 W4 is confirmed. For M2, the outcome W2 rests on the Monte Carlo noise of the 40-draw estimate: with the call-centre probability computed by quadrature there are two A2 failures (coach, call centre) under either weighting, so the outcome is W4 with no judgement call. The classroom gradient is rejected at the calibrated mode, the call centre is a failure, and a common mixing coefficient for vehicles is rejected at p = 0.023. |
| 12 | Code checked separately from the outbreak data | confirmed | The check's separately written implementation reproduces the train calibration to three decimals in log-likelihood, the hold-out p-values to Monte Carlo accuracy, the power table and the M1 ceilings. The first analysis's current code regenerates the stored predictive pmfs exactly. Frozen hashes are unchanged (recomputed on the unredacted originals; the released protocol files are redacted copies whose hashes differ, see `code/REDACTIONS.md`), and no file in bsc_outbreaks, bsc_sim or the first analysis's data/ was modified by the check's run. |
| 13 | Leakage between calibration and hold-out | partially confirmed | No numerical leakage found. 01_calibrate.py and 02_power.py do not read heldout_outcomes.json; events.py contains only stratum sizes and totals; the per-event infectiousness is profiled from the event total, and the shape test conditions on that total. Leakage through knowledge of the outcomes cannot be excluded: the split, the strata, the hold-out ventilation priors and the choice to exclude only the adjacent-seat cell were fixed knowing the outcomes. That last choice decides the result: with the cell included the bus fails and pooled Delta is -4.0; with the whole index row excluded nothing fails and pooled Delta is +6.0 (both reproduced). |
| 14 | Overfitting (parameters vs data points) | confirmed | Vehicle class: 2 fitted parameters on 22 cells (138 cases among 70,086 contacts after exclusion); no parameter is fitted to hold-out strata. Room class: 2 parameters on 2-5 informative cases, giving a posterior that spans more than three decades. That is under-identification, not overfitting, and it is what makes the classroom result uninterpretable. The S-aniso variant (3 parameters) is not used for the primary result. |

Further points:

* (1) The protocol's 40-draw estimate of the call-centre probability is Monte Carlo noise. Quadrature over the same posterior gives P = 0.007 (equal weight) and 0.013 (pooled), a clear A2 failure for M2 under both conventions, so the outcome is W4 without any judgement call. The value also depends on the arbitrary upper limit of the prior on D.
* (2) The classroom agreement comes from the width of a posterior calibrated on 2-5 restaurant cases. At the posterior mode the classroom is rejected (p = 0.0009), and 80% of the posterior needs more than 1e4 quanta/h (median 1.4e7 with the frozen propagator).
* (3) A constant near/far risk ratio of 2-3 with no calibration scores +7.1 to +7.7 on the pooled test (M2: +3.8 with the frozen propagator, +2.83 with the round-off-safe one) and has no failure on the four single events. The rule therefore does not show that the lattice kernel or its calibration adds predictive value.
* (4) The positive result rests on one event, not two. Without the flight the pooled Delta is negative and A1 fails. The flight has 20 exposed passengers, including four companion couples seated side by side (three are case-pairs). The companion-confounded cell is excluded in the train data, but no equivalent treatment exists for the flight; within-couple dependence makes the exact p-values (M0 0.0008, pooled 1.4e-4) too small.
* (5) Everyone is seated in every event used. A pass for M2 supports standard aerosol-diffusion physics, not the heterogeneity framework.
* (6) The emission M2 needs for the hold-out events where it is not rejected is about 1,600 quanta/h on the flight (range 435-11,800) and physically impossible in the classroom (median 1.4e7 with the frozen propagator); attack-rate levels are not predicted.
* (7) The pre-specified sensitivity analysis S-H1-confirmed (choir 32/60, 1,034 quanta/h) reduces the 76-fold spread of fitted emissions to about 53-fold.
* (8) The negative-control criterion (fail only if more than half the predictive mass exceeds the bound) has almost no power, and the Utah level check becomes a failure at mask factor 0.2, inside the stated range; passes of the controls and of A3 carry little weight.
* (9) Exploratory, post hoc, not usable as a claim: with the train-calibrated coefficient (about 25 m2/h, a physically plausible eddy diffusivity) the classroom is predicted well (p = 0.59, 8.7 of 10 front desks expected).
* (10) Not modelled: the Zhejiang bus passengers also attended a 150-minute event at which other attendees had a 4.1% attack rate (about 3 of the 23 bus cases could be background); two coach passengers boarded 40 minutes late and four left early.
* (11) The posterior for D_air ignores the over-dispersion of the train fit (deviance 72.6 on 20 d.f.). The check verified that tempering the likelihood by 3.6 does not change any verdict (coach p = 0.026, pooled p = 9e-5); the untempered interval 9.7-43 m2/h is narrower than the tempered one (8-96), as the Supplementary Material states.

Summary:

The calibration, the hold-out predictions, the power table, the computations for the same-site model and the
call-centre test reproduce with separately written code, and the counts reproduce from the sources. The protocol
was frozen at 05:12, before the first model run at 05:13, on the evidence of file times and hashes only, and it
was written with the counts known. The same-site model M1 fails (outcome W4). M2 fails on the coach (p = 0.006)
and on the call centre (0.007 and 0.013 by quadrature), so its outcome is W4. Its one robust positive result is
the flight (p = 0.61; well-mixed rejected at 0.0008), which carries the pooled test; the classroom agreement
depends on a coefficient that the restaurant calibration leaves undetermined over three decades. A generic
constant near/far ratio passes the same rule with a higher score.

### 4.6 Second outbreak test

*How it was checked.* The geometry of the eleven events, the statistics and the power analysis were
implemented again (`code/checks/outbreaks2/`), with a third lattice propagator built by uniformisation that
shares no code with the two propagators of the analysis. Every seat symbol of the four held-out flights and
of one calibration flight was read again from the published figures, and the other counts were typed again.

*What was found.* Every number that decides the outcome was reproduced (pooled scores +14.30 and −9.20 in the
check, +14.35 and −9.15 in the analysis), and so was the numerical correction of the room kernel. Three
qualifications apply: nine-tenths of the pooled deficit against the constant ratio comes from the
meat-processing line, whose near zone was selected by the source's authors on the same data; on aircraft the sign
of the comparison with the constant ratio depends on the sensitivity variant; and the constant ratio and the
well-mixed model also fail the adequacy test on some events. A kernel that decays exponentially with
distance, added by the check, performs as the lattice kernel does.

*Note on the items below.* The pooled gain over the well-mixed model is +14.35. The closed form of the same-site model, scored by the same rule, gains +15.45 over the well-mixed model
and has the same formal outcome V5, with pooled −8.06 against the constant ratio and one more adequacy failure.
The protocol had fixed in advance that this closed form overstates the reach of the same-site model (the first
test showed it off by up to 9.9 log-score units on one event) and that its pass would not rehabilitate that
model, so the gain is not a result for the same-site model. For M2 the pre-specified pooled comparison with the
constant ratio is significant in favour of the constant ratio (−9.15, mirror p below 1.5 × 10⁻⁵), mostly through
the processing line. On the four flights the comparison with the constant ratio does not have a determined sign: over the
15 listed sensitivity variants it runs from −3.42 to +1.21; three variants that widen the outcome definitions
make the constant ratio significantly better (mirror p = 0.013, 0.026, 0.005) and one makes the kernel
significantly better. The variants were computed with the room kernel as first frozen.

*Items: 16* (12 confirmed, 3 partially confirmed, 1 unverifiable).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | Pre-registration frozen before results and followed | partially confirmed | SHA-256 of PREREGISTRATION.md, digitized.py, lattice.py, stats.py and of the addendum files all match PREREG_FREEZE.json (recomputed on the unredacted originals; the released protocol files are redacted copies whose hashes differ, see `code/REDACTIONS.md`). File times are in the right order: prereg 09:41:14, first kernel table 09:49:53, calibration 09:51:33, power 09:52:12, addendum 09:52:41, hold-out scored 09:53:08. Split, strata, decision rule and outcome V5 were applied as written. Weaknesses: (1) the freeze is self-attested (no external timestamp; PREREG_FREEZE.json itself was rewritten at 09:52:41); (2) the first analysis had seen all outcomes, as disclosed; (3) events.py, engine.py, analysis.py and pipeline.py were written after the freeze (pipeline.py 4 s before the hold-out output); (4) CRR is coded as p = 1-exp(-eps*rho), not the registered p = rho*q; with the exact form pooled DeltaC is -7.19 (coded form -9.20), same outcome; |
| 2 | Data values re-typed from sources | confirmed | The check re-read the Naha map (all 30 rows), the Sydney-Perth mid-cabin map, the EK448 figure and seat table, the Tel Aviv-Frankfurt map and the CA112 map (Hertzberg redraw) against the stored figures; every seat symbol and every stratum count matches ([3,4,1,6], [4,4,2,1], [1,3,0], [2,0,0,0], [4,5,7,2]). Günther Table S2 gives 9/5, 26/17, 48/18, 78/20, so the bins 5/9, 12/17, 1/22, 2/30 are right. Kenyon 4/13 vs 2/55 and Olsen 8/23 vs 10/88 were re-fetched from PubMed; Yu 65/52/18% (30 of 74) and Wong 3/3, 4/8, 0/8 from Europe PMC and the stored full text. Baker and Young maps were not re-read seat by seat (calibration only). Caveats: about 22 Naha passengers lost to follow-up are counted as non-cases; most Tel Aviv-Frankfurt non-cases were never tested; Wong Fig. 4 prints the exact bed of each student, which the first analysis randomised instead (using exact beds changes pooled Delta0 by +0.13, nothing else). |
| 3 | B1: M2 vs well-mixed, pooled hold-out | confirmed | Re-implementation with separately written code (own geometry code, own statistics, a third kernel by uniformization): Delta0 = +14.30 (first analysis +14.35), p < 2.5e-6; aircraft alone +10.22. Per-event log scores agree to 0.01 on flights and 0.08 on the meat plant. As the first analysis says, this is weak evidence for the kernel: a CRR world passes B1 with probability 0.73, and a physics-free kernel exp(-distance/2.5 m) calibrated the same way scores Delta0 = +18.0. |
| 4 | B2: M2 vs constant risk ratio, pooled hold-out | partially confirmed | Numbers reproduce: DeltaC = -9.20, p(M2 better) 0.91, mirror p 2.5e-6. But -8.4 of this comes from the meat plant, where CRR's near zone (8 m) is the radius the source authors selected as the most significant on this very outcome. With a 4 m zone M2 beats CRR on the plant by about +2.2; with 12 m it is about -0.6. With near = within 1 row and 4 m, pooled DeltaC is +10.0; with within 5 rows and 12 m, +1.54 (p 0.033). 'CRR predicts better' holds only for the registered zone definitions, and on the plant that definition is not an out-of-sample competitor. |
| 5 | B2 on aircraft only (four held-out SARS-CoV-2 flights) | partially confirmed | As registered it reproduces: DeltaC = -0.74, p 0.157 for M2 better, 0.067 for CRR better, power 0.82; leave-one-out over eight flights gives DeltaC +0.02. Two sensitivity variants are significantly positive: S-hh has aircraft DeltaC +1.21 with p_B2 = 0.044 (out/06_variants.json), and CRR with the published rho = 2.4 has aircraft DeltaC +1.18 with p_B2 = 0.028 (out/05_secondary.json; pooled -4.44). Additional checks: two strata {0-2 rows},{3+} give +2.44 (p 0.003); CRR zone within 1 row +7.8; within 5 rows +2.2 (p 0.016). In the specifications tried by the check M2 ranges from tie to modestly better than CRR on aircraft; over the 15 listed sensitivity variants the comparison runs from -3.42 to +1.21 (see the note above). |
| 6 | B3: adequacy of M2 on each held-out event | confirmed | Naha p = 0.047, Sydney-Perth 0.30, EK448 0.21, Tel Aviv-Frankfurt 1.0, meat plant below 1e-5; outcome V5. Naha is fragile: p = 0.25 with households collapsed (then outcome V3), 0.018-0.037 under other binnings, 0.026 in leave-one-out. In addition, CRR itself fails adequacy on EK448 (p = 0.0195) and the well-mixed model fails on three events. |
| 7 | Meat-processing line | confirmed | M2 expects 2.9, 4.8, 5.7, 6.6 against 5, 12, 1, 2; log p -15.32 (first analysis -15.24) vs CRR -6.86 and M0 -19.40. At its own best D (40 m2/h) M2 reaches -9.46. The failure is real, but mostly a failure of the ward-to-plant calibration transfer (the room posterior is nearly well mixed because of W2); the CRR comparison is biased by the 8 m zone as noted above. |
| 8 | Common mixing coefficient within a setting class | confirmed | Rooms: LR 17.7 on 2 d.f., p = 1.5e-4 (first analysis 17.4). Aircraft: common D 158 m2/h, LR 15.0 on 7 d.f., p = 0.036. Per-event preferred D and intervals match the table of the first analysis. The four SARS-CoV-2 flights alone are also borderline heterogeneous (LR 7.9 on 3 d.f., p = 0.049). |
| 9 | S-T4: frozen transfer of the first-test train posterior to eight flights | confirmed | First-test posterior D 20 (9.7-43) reproduced; Delta0 +7.56, DeltaC -8.29; adequacy failures CA112 0.020, Cancun-Birmingham 0.004, Naha 0.0002, Sydney-Perth 0.015; EK448 0.72 and Tel Aviv-Frankfurt 1.0 pass. Caveat: the CRR comparator here was calibrated on F1-F4, which are among the predicted flights, so DeltaC over eight flights is unfair to M2; the adequacy failures do not depend on that. |
| 10 | S-swap: calibrate on SARS-CoV-2 events, predict pre-2020 events | confirmed | Aircraft Delta0 +6.46 (p 1e-4), DeltaC +3.62 (p 0.0037 under CRR); rooms DeltaC about -3.8 with W2 failing (p < 1e-5); pooled DeltaC -0.22, outcome V3. 'pass_weak' is a fair label: registered as secondary, and the win comes from CRR's rho learned on 2020 flights over-predicting CA112 and the gradient-free influenza flight. |
| 11 | Descriptive near/far risk ratios | confirmed | Recomputed from its own 2x2 tables: random-effects 3.98 (2.32-6.82), fixed 3.16, Q 18.9 on 10 d.f. (p 0.041); aircraft 3.75 (2.01-7.00); SARS-CoV-2 6.66 (3.38-13.1); SARS-CoV-1 2.34 (1.53-3.59); influenza 1.15 (0.35-3.72); tuberculosis 8.46 (1.73-41.3). All match. |
| 12 | M1cf (same-site model, closed form) on the hold-out | confirmed | Recomputed with a positive (uniformization) kernel: adequacy failures Naha 0.008, EK448 0.013, plant 0.004; aircraft DeltaC -4.62, mirror p 0.0014; outcome V5. The M1cf room numbers do not depend on the kernel: the positive kernel gives the same room posterior (mode 1.58, 0.54-23 m2/h) and the same plant p-value. |
| 13 | Numerical defect in the frozen spectral kernel and the post-hoc fix | confirmed | The check's uniformization kernel (non-negative series, independent of both kernels of the first analysis) agrees with the image-sum kernel to 1e-14 at small D, while the spectral kernel is wrong by many orders of magnitude and returns non-positive exposures at up to 32 of 74 ward receivers. The fix is legitimate and the corrected room posterior (mode 794, 305-4530 m2/h) reproduces; aircraft numbers do not change. The image kernel is itself inaccurate at large D (8% at D = 1000 in a cabin) but is only substituted where the spectral value is below 1e-8 of the maximum, i.e. small D, so no harm. First test: on the first-test restaurant calibration the floor changes the log-likelihood by less than 1e-7 at every D in its posterior range, so that calibration is not affected; first-test hold-out predictions conditioned on the event total were not re-run. |
| 14 | Power analysis of the decision rule | confirmed | 400,000 independent draws, corrected kernel: P(B1 and B2) = 0.927 if M2 is true, 0.048 under CRR, 0.010 if well mixed; P(V1) = 0.76, 0.037, 0.007. Aircraft only: 0.822, 0.050, 0.003. Critical values agree. |
| 15 | Does a generic model do as well as M2? | confirmed | Yes. A static kernel exp(-distance/ell) with one length per class (aircraft posterior ell 2.5 m, 1.6-5.7), no lattice and no airborne removal, behaves like M2: as target against CRR it gives aircraft DeltaC -1.06 (M2: -0.74) and fails the same two events; M2 against it is +0.31 on aircraft (not significant) and -3.7 pooled (the generic kernel is better on the plant). |
| 16 | Cherry-picking of events | unverifiable | Exclusion reasons were fixed in the pre-registration and are plausible (no identified source, one or two cases, pre-flight exposure, no positions). The check found no excluded event that would obviously favour or hurt M2, but the check did not run a literature search of its own. One inconsistency: AF171 was dropped as 'secondary tabulation only' while CA112 and the tuberculosis flight are themselves taken from a redraw and an abstract; with 2 cases it would not change anything. |

Further points:

* (1) The pooled 'CRR predicts better than M2 (p < 5e-6)' is 91% the meat plant, where CRR's 8 m near zone was chosen by the source authors as the most significant radius on the same outcome.
* (2) Two sensitivity variants make M2 significantly better than CRR on aircraft: households collapsed (DeltaC +1.21, p 0.044) and published rho = 2.4 (DeltaC +1.18, p 0.028); three that widen the outcome definitions make CRR significantly better (see the note above).
* (3) The near-zone boundary of CRR and the strata are forking paths with large leverage on aircraft: DeltaC runs from -0.74 (registered) to +2.4 (two strata, p 0.003) and +7.8 (zone within 1 row). The registered two-row zone is the conventional contact-tracing rule, itself derived from earlier outbreaks including the Kenyon flight used for calibration.
* (4) A physics-free exponential distance kernel does as well as M2, so nothing in this analysis tests the lattice-walk or airborne-removal structure as such.
* (5) CRR also fails adequacy on one held-out event (EK448, p = 0.0195); B3 is applied only to M2.
* (6) CRR is implemented in complementary log-log form, not the registered p = rho*q. Exact form: pooled DeltaC -7.19 (coded form -9.20); conclusion unchanged.
* (7) Naha: about 22 passengers lost to follow-up are counted as non-cases (mostly in far rows, which flatters the gradient). The paper's own 'within two rows' group is 14 persons, not the 23 in rows 21-25, so its near definition differs from the one used here.
* (8) Wong 2004 Fig. 4 gives each student's exact bed; the analysis randomises over bed lists instead. Effect negligible.

Summary:

A re-implementation with separately written code and a third kernel reproduces every headline number and the
registered outcome V5 (pooled Δ0 = +14.30 against +14.35, ΔC = −9.20 against −9.15; four held-out flights
ΔC = −0.74; power P(B1 and B2) = 0.927 if M2 is true, 0.048 under the constant ratio, 0.010 if well mixed). The
seat maps of the hold-out flights reproduce, and the post-hoc numerical fix of the room kernel is correct. The
pooled advantage of the constant ratio comes mostly from the meat plant, whose near zone was chosen on the same
outcome; on aircraft M2 and the constant ratio are not distinguished; and M2 is not distinguished from a generic
exponential distance decay, which the check added after the results. The freeze is self-attested and
data-aware. The Baker and Young calibration maps were not read again seat by seat, the hold-out predictions of
the first test were not re-run, and no literature search for omitted outbreaks was made.

### 4.7 Contact records

*How it was checked.* The ten data files were downloaded again and their hashes compared; the loader, the
contact statistics, the walk simulator and the epidemic simulation were written again
(`code/checks/contacts/`); three stages of the stored analysis were re-run from the frozen code.

*What was found.* The failure of the lattice walk on all six datasets was reproduced. The confirmatory test
of the walks anchored to home sites is uninformative, because the frozen fitting procedure does not recover
its own synthetic truth; the tolerance of the epidemic criterion is at the level of the noise; and a static
network of pair contact rates, added by the check, does better than any walk.

*Items: 10* (6 confirmed, 4 partially confirmed).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | Pre-registration frozen before results, and followed | confirmed | All PREREG_FREEZE hashes match (recomputed on the unredacted originals; the released protocol files are redacted copies whose hashes differ, see `code/REDACTIONS.md`). The freeze is stamped 08:23:27Z and the first real-data model output (s1_InVS13.json) is 33 s later, consistent with its 19 s run time. FIX_FREEZE (09:34:24Z) precedes every s3 output, and the post-hoc log up to item 12 matches its frozen hash. Re-running the frozen code for s1/InVS13, s2d/InVS13 and s3d/InVS15 reproduces the stored JSON bit for bit. DEV/CONF roles and tolerances were set in the frozen file. Residual weakness: no version control, so integrity rests on file times and hashes; pre-freeze descriptive counts on the confirmatory datasets are disclosed. |
| 2 | Data values and sources | confirmed | The check re-downloaded the ten files used from sociopatterns.org/assets/data; all SHA-256 match the raw copies. No numbers are typed from papers. Its own loader (verify/vlib.py, which does not import a2lib) gives identical N, day counts, calibration (c, mean duration, p) and identical observed statistics on InVS13, Thiers13 and LH10. No sign of cherry-picked datasets: the unused Malawi and Hypertext files were declared unused before the freeze. |
| 3 | Test 1: lattice walk (RW0, RWhet) vs contact structure and SIR outcomes, six datasets — reported fail (A0) | confirmed | The check's separately written walk simulator and SIR simulation reproduce the result within Monte-Carlo noise. Examples (check vs first analysis): InVS13 partners per day 5.8 vs 6.0 (observed 4.6); Thiers13 33.8 vs 33.8 (observed 14.6); Thiers13 index-case secondary cases at R0=1.5, 1 day: model 1.40 vs 1.39, real 0.86 vs 0.83. Recomputed from the stored JSON: S-families 1,2,2,1,1,1 (RW0) and 2 each (RWhet); E 0/4 on all six; secondary cases too high in 24/24 cells (ratios 1.15–1.86); WM and BLOCK also 0/6. Contact time in the epidemic comparison is matched (ratio 0.98–1.01), so the over-prediction is not a rate artefact. Two qualifications: it holds for one identification (site = contact range, M = 1/c, synchronous 20 s steps), and the size of the gap depends on the high per-interval transmission probability that the nominal-R0 design forces on sparse datasets (0.06–0.12 on InVS13). |
| 4 | Test 2: modifications on held-out settings (F3-O primary, F3-Z, F2) — reported fail (A0) | partially confirmed | The numbers reproduce (RWstickO E scenarios 2, 0, 0; S-families 4/7, 4/7, 2/6; RWstickZ and RWclus 0/12), and A0 is the correct reading of the rule. But the test cannot confirm a correct model. When the check generates synthetic data from RWstickO itself and applies the frozen fit, it passes E in 1 of 4 replicates and reaches the per-dataset A1 criterion in 0 of 4; three fits are degenerate (M of 9,000–160,000 sites, contact time 2.9–4.9× the target). RWclus on its own truth passes E in 0 of 3. Cause: the home-cluster width is in lattice sites, so total contact barely depends on M and the M iteration diverges; the common rescaling then no longer matches contact rate. The Thiers13 degeneracy is this defect. So the outcome is a failure of the fitting procedure, not evidence against home anchoring. |
| 5 | Test 3: descriptive reduction of epidemic error by home anchoring | partially confirmed | All table values recompute exactly (attack-rate error 0.167 → 0.041 / 0.038 / 0.050; log secondary-case error 0.361 → 0.191 / 0.166 / 0.108; scenarios 2, 3, 5 of 24). 'More than half' holds for attack rate; for secondary cases under F2 the reduction is 47%. In addition, RWclus and RWstickZ still over-predict secondary cases in 24/24 cells (ratios 1.06–1.36). Contact time for the anchored variants is not exactly matched (0.92–1.11, plus 2.34 on Thiers13). A non-spatial static-network model does better on the pre-registered rule (see next item). Post hoc and label-dependent. |
| 6 | Does a generic model do as well? (well-mixed, constant ratio, simple static network) | confirmed | WM and BLOCK match the lattice walk on E (0/6 each), as reported. The check added a static-network null with no space and no walk: each pair's contact rate taken from the training days, geometric durations, one common rescaling. It passes 12 of 24 scenarios and E on 3 of 6 datasets (LyonSchool 4/4, InVS15 3/4, SFHH 3/4; InVS13 0, LH10 0, Thiers13 2), against 0/6 for every walk variant. On the confirmatory set it passes E on 2 of 3 where WM fails, which is A2 under the frozen ladder. Caveats: it was added after seeing the results, its contact time is 2–9% low, and it is a crude estimator (2.3× contact time unscaled on InVS13). The conclusion that the walk adds nothing over a non-spatial model is therefore conservative. |
| 7 | Test 4: density scaling (S7) | partially confirmed | Exponents recompute as stated (observed 0.45, 0.95, 0.76, 0.81, 2.83; RW0 1.95–2.20; 0 of 5 within 0.3; LyonSchool correctly excluded for range 1.27). The estimator drops bins with zero contacts and occupancy is inferred from each person's own first and last contact, so the absolute exponents are biased on sparse datasets and confounded with schedule. Model and data go through the same estimator, so the qualitative contrast (mass action vs sub-linear in four settings) stands as a description only. |
| 8 | Test 5: burstiness and day-to-day repetition not reproduced by any walk variant | confirmed | S2 passes on at most 1 of 6 datasets for every model, while the real training days pass S2 on 6/6. Persistence S5 is too low for the free walk (1/6) and too high for anchored variants. The check's separately computed statistics agree (InVS13 walk burstiness 0.24 vs observed 0.36; Thiers13 0.16 vs 0.51). The check's static-network null also fails S2 on all but InVS13, so this deficit is shared by memoryless generators in general. |
| 9 | Test 6: power of the decision rule — reported pass_robust | partially confirmed | Recomputed from the stored 01_power.json: true or nesting model passes everything in 20/20; wrong non-nested model reaches A1 in 0/44; a wrong model reaches the A2-type outcome in 4 cases, as disclosed. The check did not re-run the synthetic generation. The label 'pass_robust' reads stronger than the evidence: (1) The 44 cases are 11 cells × 4 replicates and are not independent, so the 6.6% bound is optimistic. (2) The power analysis covers only the base models; the modified models fail on their own truth. (3) On real data the E tolerance sits at the noise ceiling: the real training days, run as a predictor of the test days at equal nominal R0, pass E on only 4 of 6 datasets (InVS13 2/4, Thiers13 1/4), and two identical 4,000-run SIRs on the real network differ by up to 7% in secondary cases against a 10% tolerance. A1 on at least 4 of 6 had no slack for any predictor. The lattice walk, whose errors are 25–72%, fails regardless; but it makes the borderline failures of the modified walks uninformative. |
| 10 | Leakage, forking paths, overfitting | confirmed | No leakage in favour of the reported failures. Presence windows come from the test days for every model alike, and the epidemic rescaling uses the test total for every model alike. Families and tolerances are fixed in frozen code; the unused metrics are in no decision. The stage-D refinement was committed before its result and kept although it hurt (InVS13 4/4 → 0/4), which the check reproduced. |

Further points:

* (1) The confirmatory test of the modifications has almost no power. Synthetic data generated from RWstickO or RWclus and refitted by the frozen procedure fail E in 6 of 7 replicates and reach A1 in 0 of 7, with degenerate fits in 4 of 7. The stage-3 'fail' says the fit procedure is broken (M is not identifiable when cluster width is in lattice sites); it does not say home anchoring is wrong.
* (2) The E tolerance is at the noise ceiling on real data. Real training days, used as a predictor of the test days, pass E on only 4 of 6 datasets (InVS13 2/4, Thiers13 1/4 scenarios). Repeating the 4,000-run real-network SIR with another seed moves index-case secondary cases by up to 7%.
* (3) A static-network null with no space and no walk (pair rates from training days) passes E on 3 of 6 datasets and 12 of 24 scenarios, and would be A2 on the confirmatory set. Every walk variant passes E on 0 of 6; WM and BLOCK fail like the walk.
* (4) RWclus and RWstickZ still over-predict index-case secondary cases in 24 of 24 cells (ratios 1.06–1.36). Home anchoring reduces the error by more than half for the attack rate; for secondary cases under F2 the reduction is 47%.
* (5) The epidemic comparison for the anchored variants is not at matched contact rate (generated/real contact time 0.92–1.11, besides the disclosed 2.34 on Thiers13). The common rescaling M → M/scale only works when contact is proportional to 1/M.
* (6) The verdict on the lattice walk holds for one identification: site = contact range, M = 1/c (120–1,400 sites), synchronous 20 s steps, hop probability 0.2–0.4. Walkers do not mix across the lattice within a day; other site sizes or the continuous-time walk were not tried.
* (7) The size of the over-prediction depends on the nominal-R0 design. On sparse datasets the per-interval transmission probability is large (0.06–0.12 on InVS13), so saturation on repeated contacts dominates; the gap is smallest where it is low (LyonSchool, ratio 1.15–1.3).
* (8) The power analysis has 4 replicates per cell and 11 dependent cells, covers only the base models, and was not extended to the modified models before the second freeze.
* (9) The S7 exponent fit drops 20-minute bins with zero contacts, which biases the observed exponent on sparse datasets.

Summary:

The failure of the lattice walk on all six datasets reproduces with separately written code, and every table
value recomputes from the stored results; the well-mixed and block models fail the epidemic criterion as the walk
does. The confirmatory test of the home-anchored walks is uninformative: with these models as synthetic truth
the frozen fit fails its own epidemic test in 6 of 7 replicates, with degenerate fits in 4 of 7, so the outcome
is a failure of the fitting procedure, not evidence against home anchoring. The tolerance of the epidemic
criterion is at the level of the noise (the real training days pass it on 4 of 6 datasets). A static network of
pair contact rates, added after the results, passes the epidemic test on 3 of 6 datasets against 0 of 6 for
every walk variant; it is exploratory.

### 4.8 Tracer measurements

*How it was checked.* The tracer values were typed again from the sources, the mixing coefficients were
fitted again with a separate solver, and the predictions, the pooled scores and the power table were computed
again (`code/checks/tracer/`); a complete re-run of the stored analysis reproduced every result file.

*What was found.* The computations stand. Four results need a narrow reading: the non-rejections on
the coach and the call centre are not informative, because a well-mixed model fits as well; flight VN54 and
the restaurant are not reproduced; the cabin coefficient depends on how zero sensor readings are treated; and
the pooled pattern is more probable under a constant risk ratio than under the kernel. The extension to four
flights is exploratory and fragile.

*Note on the items below.* In the restaurant the tracer was measured at seven tables
and taken from the source's flow simulation at the other nine; the remote mean 0.55 (0.32 to 0.86) mixes the
two, and its lowest value is simulated. Comparisons below of the tracer values with the 0.1 to 25 m²/h
fitted to outbreaks hold for the coefficients fitted to the first set only (25 m²/h
for vehicles; the room value, 0.5 m²/h, is compared with a regression on two houses, since rooms have no
tracer measurement). The aircraft coefficient of the second outbreak test, 126 m²/h, lies inside the range
measured in cabins, and its room coefficient, 794 m²/h, above the regression values; neither makes an
outbreak-fitted coefficient a measurement of air mixing. The tracer values in vehicles are of order 100 to
1,000 m²/h (the coach, 1,259 m²/h, lies above that range).

*Items: 14* (7 confirmed, 7 partially confirmed).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | Pre-registration freeze, timestamps and adherence (main protocol) | partially confirmed | All 33 frozen SHA-256 hashes match (recomputed on the unredacted originals; the released protocol files are redacted copies whose hashes differ, see `code/REDACTIONS.md`). File times are consistent with the stated order: protocol 09:18:00, freeze 09:18:16, kernels 09:28, predictions 09:36, power 09:43, second freeze 09:44:50, evaluation script 09:48:44, evaluation 09:49:08. s02/s03 do not read stratum outcomes (T4 counts and O1 case labels are loaded but unused). Limits: the freeze is self-attested (local mtimes, no external timestamp); the protocol is openly data-aware on both sides (outbreak counts, first-test failures and tracer tables had been seen). Deviations are logged (no 7A release exists, so 2 kernels instead of 3; T4 uses 20 D-quantiles x 7 nodes). The rule was followed as written. |
| 2 | Raw data re-typed from sources | confirmed | Ou 2022: The check rendered supplementary pages S13/S14 at 500 dpi and re-typed all 8 measured Fig. S3A points and all 45 Fig. S4A seat values plus the 7 case seats; all match the derived CSV (3/19 rear, 4/26 front). Kinahan: both workbooks' MD5 match figshare's API (queried again); the check re-read the 777 FWD sheet with openpyxl and reproduces every 5L and 5A-mirrored kernel value. Li 2021 Tables 1 and 3 re-fetched from Europe PMC: all 16 tables match (16 near / 63 remote patrons). Khanh 2020 text and Fig. 1: 20 susceptibles, 11/12 vs 1/8, seat map matches. Woodward XML re-fetched: identical hash. Train table sums to 138. Not re-done: second digitisation of Woodward Fig. 9; Cheng PDF (regression coefficients taken from the stored text only). |
| 3 | Re-implementation of the calibration (D fits) with separately written code | confirmed | With its own finite-difference Laplacian and eigen-solver (no shared code) the check gets exactly the same grid optima: coach 1,259 m2/h, cabin 5L 708 and 5A 112, carriage 891 (middle) and 355 (end). Exposure fields agree with the authors' mode-sum to 1e-11. |
| 4 | Hunan coach (pass_weak) | confirmed | Reproduced: kernel p = 1.00, log score -1.19 vs well-mixed -1.14, predicted RR 1.29 vs 1.03 observed; first-test D = 25 gives p = 0.0025. The pass is real but carries no evidence for the kernel: well-mixed fits marginally better. The 8 tracer points are nearly flat (0.70-1.03), so D is only bounded below. |
| 5 | Seoul call centre (pass_weak) | partially confirmed | Its own percolation simulation reproduces z_obs = -0.64 and P(z <= obs) of about 0.09-0.10 for the kernel and 0.25 for well-mixed. But the kernel predicts positive clustering (median z +0.6 to +0.9) where the observation is negative; well-mixed describes it better. It is a non-rejection, not support. D comes from a regression on two residential rooms. |
| 6 | Zhejiang bus (inconclusive) | confirmed | Reproduced with own geometry draws: kernel p = 0.22, well-mixed 0.20, predicted RR 1.03 vs 1.60. Non-discriminating; constant ratio 2 scores 0.63 log units better. |
| 7 | Marin classroom (pass_weak) | confirmed | Reproduced: kernel p = 0.42, RR 1.76 vs 2.8; well-mixed rejected at 0.036. This is the one event where the tracer-based kernel beats well-mixed without an artefact. Constant ratios 2 and 3 score better (-0.34, -0.59), and D is transferred from Cheng's residential regression, with the teacher position assumed. |
| 8 | Flight VN54 (inconclusive / formal pass p = 0.14) | partially confirmed | Numbers reproduce exactly (mixture 0.136; 5L alone 0.0007; 5A-mirrored 0.23; zeros-as-missing 0.004; diffusion form 0.045). The pass is an artefact of two things: averaging the predictive pmfs of two kernels that disagree, and exact-zero sensor readings (3L, 4L, 1D in the 5A release, while neighbouring 1L, 2L, 6L read hundreds). If the two normalised kernels are averaged before prediction, p = 0.044. VN54 is not reproduced. |
| 9 | Guangzhou restaurant (pass_weak) | partially confirmed | Numbers reproduce (p = 0.095, 0.027, 0.007, 0.002 for K = 2..5). The source's own count is 5 of 16 vs 0 of 63 (Li Table 2), at which the kernel is rejected at 0.002. It passes only under the lenient registered rule and only at K = 2. The restaurant is adequate only at the lower bound K = 2. |
| 10 | High-speed trains (fail) | confirmed | Own computation: deviance 112.7 on 21 d.f., 14.6 expected same-row cases vs 50 observed. Every model fails, including constant ratios (best at ratio 5, deviance 54) and a row-saturated model (45.2 on 18 d.f., p = 0.0004), so the table has structure no row-distance model captures. Counts are reconstructed from published percentages. |
| 11 | Pooled decision rule and power table | confirmed | Own enumeration: +7.86 vs well-mixed (p = 1.4e-6), -0.81 vs ratio 2, -0.79 vs ratio 3; without the flight +3.40 (p = 0.0012). Power table reproduced (S4: 0.20 / 0.40 / 0.45 / 0.20). Constant ratios 2 and 3 beat well-mixed by more (+8.67, +8.65) than the kernel does. Added by the check: the observed pattern 'beats well-mixed but not both constant ratios' has probability 0.23 if the kernel is true, 0.84 under ratio 2 and 1.00 under ratio 3, so the pooled evidence favours a generic ratio. Formal outcome S4; outcome S5 if the flight is also counted as a failure (trains plus flight). |
| 12 | Tracer D vs outbreak-calibrated D (descriptive) | partially confirmed | The direction holds: tracer D is of order 100-1,000 m2/h against 0.1-25 fitted to outbreaks. Three qualifications. (1) The cabin estimator drops all zero readings (rows 9-10 and others); linear or square-root least squares keeping zeros gives 56-450 m2/h, so 708 is partly an estimator artefact, and the diffusion form fits the cabin field poorly (residual SD about 1.1 log units for 5L). (2) Cabin decay lengths are 1.8-4.5 m, so the decay lengths range from about 2 to 14 m. (3) At the low cabin value the tracer D is 4.5 times the 25 m2/h fitted to vehicles. The unused third release (5G) gives 126. |
| 13 | Extension EXT-F: four further flights (pass_weak, outcome E2) | partially confirmed | All numbers reproduce (p = 0.054, 0.26, 0.25, 1.00; +10.36, +2.83, -0.16), and its own solver matches the exposure to 1e-11. Concerns: (a) protocol 10:24:56, script 10:25:33, results 10:25:35, evaluation 10:26:08 - written within 72 seconds, so the freeze is nominal and the test is fully outcome-aware, as labelled. (b) 'No rejection' is knife-edge: Naha p = 0.054 with the 708/112 mixture, 0.044 with the single geometric-mean D of 282, 0.035 if the third release (126) is included, 0.003 with 126 alone - each gives E4. (c) the second outbreak test's seat maps are unverified and carry options of the first analysis that move results (Naha households removed: p = 0.18 and +1.5 vs ratio 3; Speake alternative outcome: p = 0.052). Robust part: any tracer-range D beats well-mixed by 8-10 log units and is about equal to a constant ratio of 3. |
| 14 | Leakage, forking paths, cherry-picking | partially confirmed | No outcome leakage into kernels: D and the direct kernels use tracer data only, and the check re-derived them without outbreak files. Forking paths that favour a pass: mixture-of-pmfs over disagreeing kernels (VN54, extension), the R1 'fails only if rejected at every K' rule, positive-only log fit for the cabin. Event selection is inherited from the first test with exclusions fixed in the protocol; the extension then reverses two of those exclusions post hoc (disclosed). No evidence of tuning to fit; sensitivity variants that reject are reported in full. |

Further points:

* (1) The pooled result is itself evidence for a generic model: the pattern actually obtained (beats well-mixed, does not beat both constant ratios) has probability 0.23 under the physics kernel, 0.84 under a constant ratio of 2 and about 1.00 under a ratio of 3.
* (2) Constant risk ratios of 2 and 3 beat the well-mixed model by more (+8.67, +8.65 log units) than the physics kernel does (+7.87) on the same five events, so the gain over well-mixed (p = 1e-6) is not specific to the kernel.
* (3) Cabin D is estimator-dependent: the registered log least-squares fit drops every zero sensor reading (9 of 39 in release 5L, 14 of 39 in 5A). Fits that keep zeros give 56-450 m2/h. The Kinahan field is not diffusive (strong forward bias, hard cut-off aft of row 8), so 'tracer-measured D' for the cabin is a loose summary, and the extension rests on it.
* (4) Sensor 7A reads about 21,000 counts in the 5A release (six times the adjacent 6A) and is included in the D fit; excluding the four extra sensors moves the 5A fit from 112 to 178 and the 5L fit from 708 to 891.
* (5) Adequacy on all four flights of the extension fails under three equally defensible physics choices (single geometric-mean D 282: Naha p = 0.044; all three forward releases: 0.035; centre release alone: 0.003). The main protocol itself defined D_cabin as the geometric mean.
* (6) The extension protocol, its prediction script and its results carry timestamps within 72 seconds of each other, so the code was drafted together with or before the protocol; the second 'before unblinding' freeze adds nothing because the near/far counts were already known.
* (7) For the call centre the kernel predicts the wrong sign of clustering (median z about +0.7 vs observed -0.64); it is no longer rejected, nothing more.
* (8) For the restaurant, the source's own table gives 5 of 16 vs 0 of 63, at which the tracer kernel is rejected at p = 0.002; the K = 2 reading is carried over from the first test.
* (9) VN54 was flown on a different aircraft type from the 777-200 tested, and Naha is a narrow-body 3-3 cabin while all tracer data are wide-body; the transfer is untested in both directions.
* (10) Freezes are self-attested local files with no external timestamp or version control; hashes prove the files are unchanged since, not when they were written.

Summary:

Every source value that was typed again agrees, and the re-implementation reproduces all kernels, p-values,
pooled scores and the power table to 2–3 digits. The registered outcome is S4 as computed; it is S5 if the flight,
whose formal pass depends on zero sensor readings and on averaging two contradictory kernels, is counted as a
failure. The tracer-calibrated kernel is not shown to be a spatial predictor: its non-rejections on the coach
and the call centre are at events where a well-mixed model does as well or better; its one gain over well-mixed
(the classroom) is matched by a constant ratio; and the observed pooled pattern is 3–4 times more probable under
a constant ratio of 2–3 than under the kernel. The extension to four flights is outcome-aware, uses the seat maps
of the second outbreak test, and has one rejection (Naha) under the geometric-mean cabin coefficient that the
main protocol specified. The robust part is that a cabin coefficient of order 100 m²/h beats well-mixed and ties
with a constant ratio of 3.

### 4.9 Zone level

*How it was checked.* The school data were downloaded again and compared; the likelihoods, the
cross-validation, the bootstrap and the forward simulations were implemented again (`code/checks/zones/`),
and the power of the rule was examined for the generic model at the ratio that the data favour.

*What was found.* The full-data likelihoods were reproduced to six digits, and the outcomes of the three
tiers stand. Two points bound what the pass means: a zone model with guessed shares predicts held-out schools
at least as well as the measured shares, and a generic model with its fitted ratio earns the verdict in about
38 % of synthetic replicates.

*Items: 14* (7 confirmed, 5 partially confirmed, 1 refuted, 1 unverifiable).

| No. | Item | Outcome | Notes |
|---|---|---|---|
| 1 | Pre-registration frozen before results and followed | partially confirmed | SHA-256 of PREREGISTRATION.md and all 13 frozen files equal to PREREG_FREEZE.json (recomputed on the unredacted originals; the released protocol files are redacted copies whose hashes differ, see `code/REDACTIONS.md`). File times are consistent with the claim: prereg last modified 09:48:58, freeze 09:48:59, matsumoto_primary.json created 09:50:13, seoul_floors.json 10:21:08; scripts 04/05, zonecore and analysis unchanged since freeze. The registered pipeline was followed (same tiers, thresholds, folds, seeds). Limits: file-system evidence only, no version control; the CSV with outcomes existed 44 minutes before the freeze, so blindness to class-level outcomes rests on the author's statement; the 0.15 tolerance is a judgement roughly 2 to 3 times the sampling noise under a true Z (95th percentile TV 0.045 to 0.071). |
| 2 | Data values re-typed from original sources | confirmed | Re-downloaded all three sources: Matsumoto JLD2 from GitHub commit 658abfe (hash identical), SocioPatterns primaryschool.csv.gz and metadata (hashes identical), Europe PMC XML of PMC7392450 (hash identical). Parsed the JLD2 separately using the compound-type member offsets: 10,923 students, 2,548 infected, 29 schools, 455 classes, onset days 20 to 151; school, grade, class, infection and onset columns identical to the first analysis's CSV. Park Table 1 re-read: 84/0, 182/0, 207/0, 206/1, 27/2, 216/94, 201/0, 20/0, total 1,143/97 all correct (the 20 are labelled 'Other' in the table). Lyon shares recomputed: 0.762/0.1045/0.1335, day 1 0.756/0.111/0.133, day 2 0.768/0.098/0.134. The 11th-floor wing counts (79/137, 4/62) are digitised from a figure and were not checked. |
| 3 | T1 Tier 1: out-of-sample ranking against well-mixed and independent-class models | confirmed | Re-implementation with separately written code reproduces Z-H = +798.4 (95% 563 to 1,074) and Z-X = +211.1 (109 to 345) with the same lower bound on the baseline hazard; Z ahead in 27 and 23-24 of 29 schools. One fragility found: in the odd-school calibration fold the baseline hazard sits on its lower bound (1e-9), and held-out scores depend on that bound; varying it from 1e-9 to 1e-5 moves Z-X only from 211 to 203 and Z-H not at all, and letting it go to zero gives Z-H = +745 (531 to 999). The tier holds in every variant. The check's power replicates: well-mixed truth passes 0/10, independent-class truth 0/10. As the first analysis says, the direction was foreseeable and this tier alone is weak evidence. |
| 4 | T1 Tier 2a: strict equality of fitted and Lyon shares | confirmed | Fail reproduced exactly: full-data log likelihoods identical to the first analysis's to 6 digits (Z -14,795.14, F -14,773.78), LR = 42.71 against 5.99; fitted shares 0.659/0.104/0.236. Under a true Z the check's synthetic replicates give LR 0.7 to 4.2, so the observed 42.7 rejects Lyon shares as the exact structure. |
| 5 | T1 Tier 2b: fitted shares within tolerance 0.15 of Lyon shares | partially confirmed | TV = 0.1028 reproduced; the check's school bootstrap gives 0.052 to 0.169, with 11% of resamples above 0.15. Sensitivity variants reproduced: no winter-break switch 0.162 (outside), kernel mean 2.5 d 0.094; two variants the check added: shorter kernel (mean 1.2 d) 0.123, break widened by 3 days 0.090. The pass is mechanical and not secure, as the first analysis states. The tier is also weakly discriminating: round-number guesses with no measurement pass it too and sit closer to the fitted shares than Lyon does: (0.7, 0.1, 0.2) TV 0.04, (2/3, 1/6, 1/6) TV 0.07, (0.6, 0.2, 0.2) TV 0.10, (0.8, 0.1, 0.1) TV 0.14. |
| 6 | T1 Tier 3: forward-simulated between-class heterogeneity | confirmed | Observed statistic 3.383 reproduced. The check's 300 simulations per model: Z 2.72 to 3.65 (contains it), well-mixed 0.87 to 1.14 (excluded), independent classes 3.03 to 4.09 (contains it). The first analysis's reading that this tier excludes only the well-mixed school is correct. Added by the check (post hoc): the generic constant-ratio model with its fitted ratio gives 2.32 to 3.11 and so does not contain the observed value; this is the one place where Z does better than that competitor, but it is a post-hoc observation of the check and not a registered test. |
| 7 | T1 registered verdict (weak pass) | partially confirmed | The mechanical verdict 'pass' and the first analysis's downgrade to weak are both reproduced and appropriate. Two points make it weak. (1) The false-pass probability of the power table of 0 to 0.075 for a constant-ratio truth averages over ratios log-uniform on 1 to 1,000. At the ratio the real data actually favour (37), the check's replicates pass the full registered rule 6/16 times (about 38%; all pass Tier 1 and Tier 3, 6 pass Tier 2b); at ratio 10 and 100 they pass 0/8 and 0/8. The power table of the protocol shows 3/18 for ratios 30 to 100. (2) Lyon-measured shares are not shown to be better than an unmeasured guess: see the naive-shares check. |
| 8 | T1 comparison with generic constant own-class/rest ratio model | confirmed | Held-out Z-C = +18.7 (95% -8.9 to +53.8), Z ahead in 16/29 schools; full-data log likelihood favours C by 13.6 (C -14,781.55, ratio 36.6). Reproduced. 'Inconclusive / competitor not beaten' is the right label. |
| 9 | Naive-shares competitor (added by the check): do the measured Lyon shares earn the result? | refuted | The reading that the separately measured Lyon mixing earns the result does not hold. Under the same two-fold cross-validation, a zone model with guessed shares (0.7, 0.1, 0.2) beats the Lyon-share model out of sample by 19 log units (95% 3 to 38); (2/3, 1/6, 1/6) and (0.6, 0.2, 0.2) are ahead by 9 and about 25, not significantly; (0.8, 0.1, 0.1) is behind by 19 (8 to 32), (0.9, 0.05, 0.05) by 82, equal thirds by 170. So Lyon shares are in the right region but are not distinguishable from, and are beaten by one of, several 'mostly own class' guesses. |
| 10 | T1 sensitivity: frequency- versus density-dependent transfer | unverifiable | Not re-implemented; the density-dependent variant contains an unexplained normalising constant (22.2) that does not affect the likelihood after refitting the rate. Stored numbers are internally consistent with the reproduced base fit. Descriptive only. |
| 11 | T2: Seoul call-centre building by floor | confirmed | Recomputed with separately written code: expected off-floor cases median 3.97 (95% 0.86 to 18.1), predictive count median 4, interval 0 to 19; well-mixed expectation 78.7, P(3 or fewer) = 2.2e-31; generic log-uniform guess passes with probability 0.250. 'Inconclusive' is the correct label. The outcome depends on where the assumed prior is placed: shared-zone time 10 to 60 min/day gives 1 to 39 (still contains 3), 30 to 120 min gives 7 to 81 (does not), and the prior was set after a back-of-envelope look at the counts. It rules out a well-mixed building, which needs no model. |
| 12 | Seoul 11th-floor wings and early-growth / Perron statements (not run) | confirmed | The reason given (no outbreak-independent mixing measurement) is sound. The check did not search further for such datasets, so that none exists is not established by the check. |
| 13 | Power analysis | partially confirmed | Replicates of the check in scenario A (its own simulator, 60 forward runs each): true Z passes 8/10 (two Tier 3 misses; the first analysis reports 1.00 from 40), well-mixed 0/10, independent classes 0/10, consistent with the first analysis's table. The table's constant-ratio row is correct for a ratio drawn log-uniformly from 1 to 1,000 but understates the risk for the relevant competitor: about 38% at ratio 37 (6/16). A truth equal to the fitted shares (0.66/0.10/0.24) passes 8/10 while failing the strict test 10/10, which is what was observed. |
| 14 | Leakage, forking paths, cherry-picking | partially confirmed | No leakage found: Lyon shares are computed without outbreak data, cross-validation is by school, and the same code path is used for all models. Model form, kernel, break window and denominators (respondents, 82% of enrolled pupils) were fixed before the freeze. Only one outbreak dataset was analysed for T1, chosen with knowledge of the source paper's qualitative conclusion (disclosed). The external-hazard design was changed before the freeze using synthetic data only (logged). |

Further points:

* (1) A zone model with guessed, unmeasured shares (0.7, 0.1, 0.2) predicts held-out schools better than the Lyon-measured shares by 19 log units (95% 3 to 38), and (2/3, 1/6, 1/6) and (0.6, 0.2, 0.2) do at least as well; the result therefore does not show that independently measured contact time adds anything beyond 'most transmission is within the class'.
* (2) At the constant ratio the data actually favour (about 37), a generic constant-ratio world passes the full registered rule in about 38% of synthetic replicates (6/16); the 0 to 7.5% of the power table averages over ratios from 1 to 1,000.
* (3) Held-out scores depend on an arbitrary lower bound (1e-9) on the baseline hazard, which is active in the odd-school calibration fold; with the bound removed Z minus well-mixed is +745 instead of +798. Conclusions are unchanged, but the reported magnitudes are not bound-free.
* (4) The forward-simulated heterogeneity of the constant-ratio model (2.32 to 3.11) excludes the observed 3.38, while the zone model's interval contains it: a small point in favour of the zone model (post hoc).
* (5) The Lyon mixing shares depend strongly on the contact metric: by duration 0.76/0.10/0.13, by distinct contact pairs 0.43/0.14/0.43; with the latter the tolerance test would fail (TV 0.23). Duration was pre-declared.
* (6) Denominators are survey respondents (82% of enrolled pupils); infections among non-respondents are unobserved and their effect on the fitted shares was not examined.
* (7) In 36 single-class grades (453 pupils) the grade term is set to zero rather than redistributed, so the zone model's shares do not sum to one there; pre-declared, effect not assessed.
* (8) The Seoul 'consistent' outcome flips to 'inconsistent' if the assumed shared-zone time is 30 to 120 minutes a day instead of 3 to 30; the prior was chosen after a back-of-envelope calculation on the observed counts.

Summary:

The data values match their re-downloaded sources, and the re-implementation reproduces the full-data
likelihoods to six digits and the outcomes of all tiers (Tier 1 pass, strict test fail, tolerance and dispersion
weak pass, generic competitor not beaten, Seoul inconclusive). Two findings of the check, both added after the
results, limit the weak pass: a zone model with guessed shares (0.7, 0.1, 0.2) beats the Lyon-share model out of
sample by 19 log units (95 % 3 to 38), so measured mixing is not shown to matter beyond "mostly within the
class"; and a constant-ratio world at the ratio the data favour (about 37) passes the full registered rule in 6
of 16 synthetic replicates. The analysis supports class-dominated three-level mixing against a well-mixed school
or isolated classes; it does not involve the lattice walk, the within-room gradient, the airborne term or
attack-rate levels.

### 4.10 Overall statement of the tests against data

The statement that the evidence of both sets of tests supports, after the checks:

> Tested against documented outbreaks, tracer measurements and recorded human contacts, the model does not
> show general agreement, and what it gets right is not specific to it. The model with people
> performing a random walk on the lattice and transmission only between people on the same site is not
> supported: with its mobility calibrated (train cohort, restaurant) its expected attack rate falls short of the
> observed one in three of four held-out seated outbreaks of the first test, however large the transmission
> rate (the observed totals have probability 0.0002, 0.025 and 0.081), its closed form fails on three of five
> held-out events of the second test, and the contacts generated by its walk are spread over too many, too
> equal partners, so that simulated outbreaks are over-predicted, no better than by a well-mixed model. In the
> second outbreak test the framework with the non-local airborne term predicts distance gradients better than a
> well-mixed model (log-score gain +14.35 on five held-out SARS-CoV-2 events), but not better than a constant
> near/far risk ratio: pooled over those five events the pre-specified comparison favours the constant ratio, mostly
> through one event whose near zone its source chose on the same data; on four held-out flights the two are
> statistically indistinguishable (a test with 82 % power), and in a comparison added afterwards so is an
> exponential decay without physics. It fails on several events and admits no single mixing coefficient.
> Tracer-measured values in vehicles, of order 100 to 1,000 m²/h, are far above the 25 m²/h fitted to the
> vehicles of the first dataset (the room value of that dataset, 0.5 m²/h, is a few hundred times below a
> regression estimate for rooms); the aircraft value of the second test, 126 m²/h, lies inside the range
> measured in cabins, which does not make it a measurement of air mixing. The data therefore support a smooth
> one-parameter distance gradient (pooled near/far risk ratio with random effects 3.98 over the eleven records
> of the second dataset, 2.05 over the four held-out events of the first test with distance strata) and not the
> specific shape, calibration or physics of the lattice kernel. At the scale of zones, mixing shares measured in a different school predict
> an influenza epidemic in 29 schools far better than a well-mixed school or independent classes and within a
> pre-set tolerance (a weak pass), but exact agreement is rejected, a generic model with one parameter does as
> well, and rough shares of 0.7, 0.1 and 0.2 (compared after the results were known) do better; this supports a class-dominated zone structure and not the value of measured contact time.
> Attack-rate levels were never predicted, and heterogeneous mobility was not tested.

Reservations that apply to all of it:

* Every protocol was written with the data known (it is data-aware and model-blind), and the freezes are
  files on the author's machine with matching hashes, without an external time stamp or version-control
  history.
* Several findings of the checks were themselves added after the results were known and are exploratory: the
  static network of pair contact rates, the exponential kernel, the guessed mixing shares, the alternative
  near zones and strata.
* One infectiousness parameter was fitted per outbreak, so no test bears on attack-rate levels.
* People are seated or at fixed stations in every outbreak that publishes positions, and the contact records
  carry no positions: spatially heterogeneous mobility, the defining feature of the model, remains untested.

Wordings that the evidence does not support on any reading: "excellent", "good" or "close" agreement;
"validated in N settings"; any claim that the lattice kernel outperforms a generic distance rule; any claim
that a calibrated coefficient transfers between settings or that the fitted coefficient is a physical
air-mixing coefficient; "tracer calibration repairs the coach and call-centre failures" without saying that a
well-mixed model does as well; "reproduces" for flight VN54 or the restaurant with the tracer kernel; any
claim that attack rates are predicted; any validation of the illustrative parameter values of the scenes; any claim
that the random-walk contact mechanism reproduces real contact patterns; "home anchoring was tested and
failed" or "home anchoring works"; "measured contact time predicts transmission between classes" as a
quantitative law; any statement that heterogeneous mobility has been tested.

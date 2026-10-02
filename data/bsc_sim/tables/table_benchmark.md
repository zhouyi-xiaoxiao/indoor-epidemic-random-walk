
**office, D0=1** — single-core CPU seconds

| quantity | deterministic method | value | CPU (s) | Monte Carlo estimator | MC value (rel. SE of the run) | MC CPU for the target | speed-up at equal accuracy | deterministic − MC |
|---|---|---|---|---|---|---|---|---|
| mean-field R0 = ρ(K) | sparse LU + power iteration | 5.1073 | 4.63e-03 | particle power method (83 burn-in generations) | 5.108 (0.12 %) | 3.88e-02 (1 %) | 8× (1 %); 837× (0.1 %) | -0.001 |
| cases of a case born at one cell, (1ᵀK)_x | one sparse solve | 5.3345 | 2.71e-04 | single-walker path integral | 5.357 (0.22 %) | 9.41e-03 (1 %) | 35× | -0.022 |
| discrete-people R1 (uniform index) | CG on the pair lattice | 2.2661 | 6.84e-02 | individual-based runs, index alone | 2.265 (0.54 %) | 1.24e+00 (1 %, 11626 runs) | 18× | +0.001 |
| P(major outbreak) | mean-field branching fixed point | 0.781 | 9.41e-03 | full epidemics | 0.498 (SE 0.006) | 1.29e+00 (SE 0.01, 2500 runs) | 137× | +0.283 |
| attack rate of a major outbreak | mean-field lattice ODE | 0.989 | 4.04e-01 | full epidemics | 0.443 (0.80 %) | 2.61e+00 (1 %, 5063 runs) | 6.4× | +0.546 |

full epidemic: 5949 events, C core 5.15e-04 s per run; pure-Python reference 0.30 ms per event, i.e. 1.8 s per run (3480× slower).

**office, D0=10** — single-core CPU seconds

| quantity | deterministic method | value | CPU (s) | Monte Carlo estimator | MC value (rel. SE of the run) | MC CPU for the target | speed-up at equal accuracy | deterministic − MC |
|---|---|---|---|---|---|---|---|---|
| mean-field R0 = ρ(K) | sparse LU + power iteration | 4.6762 | 1.40e-03 | particle power method (27 burn-in generations) | 4.667 (0.12 %) | 1.86e-01 (1 %) | 133× (1 %); 13317× (0.1 %) | +0.009 |
| cases of a case born at one cell, (1ᵀK)_x | one sparse solve | 5.0720 | 2.85e-04 | single-walker path integral | 5.082 (0.22 %) | 9.43e-02 (1 %) | 331× | -0.010 |
| discrete-people R1 (uniform index) | CG on the pair lattice | 3.7934 | 1.37e-01 | individual-based runs, index alone | 3.794 (0.53 %) | 1.16e+01 (1 %, 11346 runs) | 85× | -0.000 |
| P(major outbreak) | mean-field branching fixed point | 0.783 | 9.73e-03 | full epidemics | 0.754 (SE 0.005) | 1.09e+01 (SE 0.01, 1856 runs) | 1125× | +0.030 |
| attack rate of a major outbreak | mean-field lattice ODE | 0.990 | 1.54e+00 | full epidemics | 0.964 (0.06 %) | 1.92e-01 (1 %, 33 runs) | 0.12× | +0.026 |

full epidemic: 72773 events, C core 5.90e-03 s per run; pure-Python reference 0.12 ms per event, i.e. 8.9 s per run (1516× slower).

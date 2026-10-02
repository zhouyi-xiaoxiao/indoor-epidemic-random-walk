# Addendum to the pre-registration: calibration result and power of the decision rule

Written 2026-10-01 09:52 BST, after calibration (`scripts/02_calibrate.py`) and before any hold-out stratum count was scored.
`scripts/03_power.py` reads the calibration posteriors and the hold-out predictive tables, which depend on hold-out geometry, stratum sizes and event totals only. 200,000 Monte Carlo replicates per truth, fixed seeds.

## Calibrated parameters (posterior on the grid, flat prior)

| Class | Model | Posterior mode | 5% - 50% - 95% | max log-likelihood (calibration events) |
|---|---|---|---|---|
| aircraft (F1-F4) | M0 | - | - | -19.90 |
| aircraft (F1-F4) | CRR | rho = 3.16 | 1.61 - 2.90 - 5.32 | -15.04 |
| aircraft (F1-F4) | M2 | D = 126 m2/h | 51.6 - 128 - 438 | -13.35 |
| aircraft (F1-F4) | M1 | D = 1.26 m2/h | 0.647 - 1.29 - 4.2 | -13.52 |
| aircraft (F1-F4) | M3 | D = 158 m2/h, w = 0.00 | D: 16.3 - 85.9 - 349; w: 0.00 - 0.06 - 0.47 | -13.39 |
| rooms (W1, W2) | M0 | - | - | -16.23 |
| rooms (W1, W2) | CRR | rho = 3.55 | 1.83 - 3.21 - 5.53 | -9.88 |
| rooms (W1, W2) | M2 | D = 0.02 m2/h | 0.0108 - 0.0188 - 2.3e+03 | -10.67 |
| rooms (W1, W2) | M1 | D = 1.58 m2/h | 0.542 - 2.32 - 22.5 | -9.72 |
| rooms (W1, W2) | M3 | D = 25.1 m2/h, w = 0.50 | D: 0.245 - 5.4 - 667; w: 0.01 - 0.19 - 0.68 | -11.08 |

The room-class posterior of M2 is bimodal and reaches the lower edge of the D grid: the students (W1) ask for a steep kernel and the inpatients (W2) for a flat one; no model fits W2 (in-sample adequacy p between 0.001 and 0.025). This is recorded before the hold-out room event (P1) is scored.

## Probability of each outcome of the decision rule, by which model is true

### Rule applied to M2

Critical values: B1 requires Delta0 > max(0, -3.07); B2 requires DeltaC > max(0, -0.80); mirror requires DeltaC < min(0, 0.47).

| Truth | P(B1) | P(B2) | P(B1 and B2) | P(mirror) | P(B3, no failure) | P(V1) | P(V2) | P(V3) | P(V4) | P(V5) |
|---|---|---|---|---|---|---|---|---|---|---|
| well mixed (M0) | 0.007 | 0.096 | 0.005 | 0.904 | 0.162 | 0.003 | 0.000 | 0.003 | 0.575 | 0.418 |
| constant risk ratio (CRR) | 0.652 | 0.029 | 0.028 | 0.971 | 0.222 | 0.020 | 0.000 | 0.452 | 0.204 | 0.316 |
| M2 itself | 0.989 | 0.961 | 0.956 | 0.039 | 0.809 | 0.784 | 0.000 | 0.031 | 0.010 | 0.015 |

Aircraft hold-out events only (F5-F8): P(B1 and B2) = 0.821 if M2 is true, 0.050 if CRR is true, 0.004 if M0 is true.

### Rule applied to M3

Critical values: B1 requires Delta0 > max(0, -2.05); B2 requires DeltaC > max(0, 0.40); mirror requires DeltaC < min(0, -0.82).

| Truth | P(B1) | P(B2) | P(B1 and B2) | P(mirror) | P(B3, no failure) | P(V1) | P(V2) | P(V3) | P(V4) | P(V5) |
|---|---|---|---|---|---|---|---|---|---|---|
| well mixed (M0) | 0.013 | 0.187 | 0.007 | 0.558 | 0.343 | 0.005 | 0.003 | 0.003 | 0.754 | 0.234 |
| constant risk ratio (CRR) | 0.813 | 0.050 | 0.048 | 0.878 | 0.474 | 0.038 | 0.064 | 0.611 | 0.154 | 0.124 |
| M3 itself | 0.980 | 0.892 | 0.885 | 0.050 | 0.800 | 0.718 | 0.051 | 0.041 | 0.019 | 0.018 |

Aircraft hold-out events only (F5-F8): P(B1 and B2) = 0.671 if M3 is true, 0.050 if CRR is true, 0.006 if M0 is true.

### Rule applied to M1

Critical values: B1 requires Delta0 > max(0, -4.08); B2 requires DeltaC > max(0, -1.03); mirror requires DeltaC < min(0, 0.76).

| Truth | P(B1) | P(B2) | P(B1 and B2) | P(mirror) | P(B3, no failure) | P(V1) | P(V2) | P(V3) | P(V4) | P(V5) |
|---|---|---|---|---|---|---|---|---|---|---|
| well mixed (M0) | 0.003 | 0.043 | 0.001 | 0.957 | 0.042 | 0.001 | 0.000 | 0.002 | 0.392 | 0.605 |
| constant risk ratio (CRR) | 0.653 | 0.025 | 0.025 | 0.975 | 0.177 | 0.017 | 0.000 | 0.448 | 0.191 | 0.337 |
| M1 itself | 0.994 | 0.966 | 0.964 | 0.034 | 0.808 | 0.790 | 0.000 | 0.028 | 0.005 | 0.016 |

Aircraft hold-out events only (F5-F8): P(B1 and B2) = 0.930 if M1 is true, 0.046 if CRR is true, 0.004 if M0 is true.

## Reading

- If outbreaks follow a constant near/far risk ratio, M2 passes B1 with probability 0.65: beating the well-mixed model alone is not evidence for the kernel (the round-1 criticism). It passes B1 and B2 together with probability 0.028, and reaches outcome V1 with probability 0.020.
- If outbreaks are well mixed, M2 reaches V1 with probability 0.003.
- If M2 is true, B1 and B2 are both met with probability 0.956 (the power condition of V1, at least 0.8, is satisfied) and V1 is reached with probability 0.78; the gap is the chance that one of five events fails its own adequacy test at the 0.05 level.
- The test therefore discriminates M2 from both generic competitors.

## Predicted bin counts for the hold-out events (no outcome read)

| Event | bin sizes | K | M0 | CRR | M2 | M3 |
|---|---|---|---|---|---|---|
| F5 | [12, 11, 30, 87] | 14 | 1.2, 1.1, 3.0, 8.7 | 2.7, 2.4, 2.3, 6.6 | 3.8, 2.7, 4.5, 3.1 | 3.6, 2.4, 4.1, 3.8 |
| F6 | [32, 14, 22, 30] | 11 | 3.6, 1.6, 2.5, 3.4 | 5.5, 2.4, 1.3, 1.8 | 5.9, 1.9, 2.0, 1.2 | 5.7, 1.8, 2.0, 1.5 |
| F7 | [4, 9, 71] | 4 | 0.2, 0.4, 3.4 | 0.4, 1.0, 2.6 | 0.8, 1.5, 1.7 | 0.7, 1.3, 1.9 |
| F8 | [15, 5, 17, 34] | 2 | 0.4, 0.1, 0.5, 1.0 | 0.8, 0.3, 0.3, 0.6 | 1.1, 0.3, 0.4, 0.2 | 1.0, 0.3, 0.4, 0.3 |
| P1 | [9, 17, 22, 30] | 20 | 2.3, 4.4, 5.6, 7.7 | 4.0, 7.6, 3.5, 4.8 | 6.7, 8.3, 2.3, 2.7 | 6.3, 5.2, 4.0, 4.5 |

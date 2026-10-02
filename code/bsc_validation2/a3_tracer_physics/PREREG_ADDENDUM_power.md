# Addendum: power of the decision rule (computed before unblinding)

Written 2026-10-01 09:43:16 BST by `src/s03_power.py` from `results/02_predictions.json` (kernels from physics data,
strata sizes and event totals). No stratum outcome was read. Exact enumeration over the joint distribution
of the five stratum events (57024 outcomes); T4 by 2,000 simulated tables per truth; O1 from the simulated
join-count distributions.

## Probability of each outcome of the decision rule

| Truth | S1 | S2 | S3 | S4 (one failure) | S5 (two or more) | no failure among 5 stratum events | D2 (stratum part) | D3 | T4 adequate | T4 score above M0 | O1 adequate |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 0.213 | 0.117 | 0.453 | 0.196 | 0.021 | 0.865 | 0.612 | 0.705 | 0.953 | 0.674 | 0.950 |
| M0 | 0.009 | 0.004 | 0.488 | 0.401 | 0.099 | 0.581 | 0.050 | 0.905 | 0.937 | 0.303 | 0.919 |
| CRR2 | 0.017 | 0.204 | 0.112 | 0.455 | 0.211 | 0.464 | 0.875 | 0.050 | 0.782 | 0.802 | 0.919 |
| CRR3 | 0.000 | 0.022 | 0.002 | 0.204 | 0.773 | 0.100 | 0.998 | 0.001 | 0.254 | 0.965 | 0.919 |

Reading: row P is the probability that the physics kernel earns each outcome if it is true; rows M0, CRR2,
CRR3 are the probabilities that it earns them when the truth is well mixed or a constant near/far risk ratio
of 2 or 3 (for these truths O1 is generated under M0 and T4 with near = same row).

Critical values of the pooled log-score difference (95th percentile under the competitor): M0: +1.17, CRR2: +1.16, CRR3: -0.96.

## Which events can discriminate

Probability that the model in the column is rejected by the two-sided 5% test when the model in the row is true.

| Event | Truth | reject P | reject M0 | reject CRR2 | reject CRR3 |
|---|---|---|---|---|---|
| T2 | P | 0.044 | 0.044 | 0.060 | 0.232 |
| T2 | M0 | 0.031 | 0.031 | 0.112 | 0.359 |
| T2 | CRR2 | 0.135 | 0.135 | 0.031 | 0.080 |
| T2 | CRR3 | 0.295 | 0.295 | 0.064 | 0.023 |
| T1 | P | 0.044 | 0.044 | 0.503 | 0.846 |
| T1 | M0 | 0.039 | 0.039 | 0.535 | 0.868 |
| T1 | CRR2 | 0.450 | 0.450 | 0.033 | 0.176 |
| T1 | CRR3 | 0.836 | 0.836 | 0.091 | 0.043 |
| T5 | P | 0.010 | 0.285 | 0.070 | 0.252 |
| T5 | M0 | 0.054 | 0.019 | 0.260 | 0.612 |
| T5 | CRR2 | 0.004 | 0.356 | 0.010 | 0.064 |
| T5 | CRR3 | 0.032 | 0.708 | 0.033 | 0.039 |
| C1 | P | 0.044 | 0.262 | 0.044 | 0.426 |
| C1 | M0 | 0.340 | 0.036 | 0.340 | 0.893 |
| C1 | CRR2 | 0.022 | 0.388 | 0.022 | 0.275 |
| C1 | CRR3 | 0.104 | 0.803 | 0.104 | 0.039 |
| R1 | P | 0.000 | 0.095 | 0.000 | 0.000 |
| R1 | M0 | 0.000 | 0.039 | 0.000 | 0.000 |
| R1 | CRR2 | 0.000 | 0.111 | 0.000 | 0.000 |
| R1 | CRR3 | 0.000 | 0.186 | 0.000 | 0.000 |

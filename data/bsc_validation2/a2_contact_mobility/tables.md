### DEV summary

| dataset | model | S-families passed / applicable | failed families | E scenarios passed (of 4) | mean R_index ratio model/real | mean dP(major) | mean d(attack) |
|---|---|---|---|---|---|---|---|
| InVS13 | WM | 1/7 | S1, S2, S3, S4, S5, S6 | 0 | 1.69 | +0.208 | +0.239 |
| InVS13 | BLOCK | 2/7 | S1, S2, S3, S4, S5 | 0 | 1.53 | +0.180 | +0.194 |
| InVS13 | RW0 | 1/7 | S1, S2, S3, S4, S5, S6 | 0 | 1.53 | +0.178 | +0.192 |
| InVS13 | RWhet | 2/7 | S1, S3, S4, S5, S6 | 0 | 1.44 | +0.156 | +0.172 |
| InVS13 | RWhome | 2/7 | S1, S2, S4, S5, S6 | 0 | 1.31 | +0.122 | +0.100 |
| InVS13 | RWclus | 4/7 | S1, S4, S5 | 0 | 1.24 | +0.100 | +0.070 |
| InVS13 | RWstickZ | 5/7 | S1, S4 | 2 | 1.11 | +0.046 | +0.018 |
| InVS13 | RWstickO | 5/7 | S1, S4 | 0 | 1.19 | +0.072 | +0.043 |
| LyonSchool | WM | 2/7 | S2, S3, S4, S5, S6 | 0 | 1.25 | +0.074 | +0.125 |
| LyonSchool | BLOCK | 4/7 | S2, S3, S4 | 0 | 1.20 | +0.067 | +0.103 |
| LyonSchool | RW0 | 2/7 | S2, S3, S4, S5, S6 | 0 | 1.25 | +0.078 | +0.125 |
| LyonSchool | RWhet | 2/7 | S2, S3, S4, S5, S6 | 0 | 1.25 | +0.076 | +0.121 |
| LyonSchool | RWhome | 3/7 | S2, S3, S4, S6 | 0 | 1.19 | +0.068 | +0.097 |
| LyonSchool | RWclus | 3/7 | S2, S3, S4, S5 | 1 | 1.15 | +0.045 | +0.056 |
| LyonSchool | RWstickZ | 3/7 | S1, S2, S3, S4 | 1 | 1.16 | +0.051 | +0.072 |
| LyonSchool | RWstickO | 4/7 | S2, S3, S5 | 0 | 0.97 | -0.061 | -0.087 |
| LH10 | WM | 2/7 | S1, S2, S3, S4, S6 | 0 | 1.37 | +0.089 | +0.119 |
| LH10 | BLOCK | 3/7 | S1, S2, S3, S4 | 0 | 1.32 | +0.076 | +0.098 |
| LH10 | RW0 | 2/7 | S1, S2, S3, S4, S6 | 0 | 1.36 | +0.089 | +0.114 |
| LH10 | RWhet | 2/7 | S1, S2, S3, S4, S6 | 0 | 1.32 | +0.084 | +0.110 |
| LH10 | RWhome | 3/7 | S1, S2, S3, S4 | 2 | 1.12 | +0.030 | +0.007 |
| LH10 | RWclus | 3/7 | S1, S2, S3, S4 | 1 | 1.13 | +0.034 | +0.020 |
| LH10 | RWstickZ | 3/7 | S1, S2, S3, S4 | 0 | 1.19 | +0.049 | +0.048 |
| LH10 | RWstickO | 3/7 | S1, S2, S3, S4 | 3 | 1.01 | -0.009 | -0.038 |

### DEV detail

**InVS13** — N=92, 10 days (train 5, test 5), 5 groups; calibrated c=0.000785 (M=1/c=1273 sites), mean event duration 2.04 intervals, p=0.303

| statistic (test days) | observed | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|---|
| total contact (pair-intervals) | 5.16e+03 | 4.66e+03 | 4.54e+03 | 4.3e+03 | 4.57e+03 | 4.65e+03 | 4.31e+03 | 4.3e+03 | 4.31e+03 | 4.79e+03 |
| mean event duration (x20 s) | 2.24 | 2.04 | 2.03 | 2.01 | 2.08 | 2.04 | 2.01 | 2.01 | 1.96 | 2.05 |
| share of contact time in events >=5 min | 0.133 | 0.0819 | 0.00135 | 0.00237 | 0.00227 | 0.0703 | 0.000732 | 0 | 0.0756 | 0.0209 |
| burstiness B of gaps | 0.359 | 0.365 | 0.00155 | 0.025 | 0.228 | 0.263 | 0.227 | 0.265 | 0.259 | 0.275 |
| recurrence 1-pairs/events | 0.691 | 0.672 | 0.159 | 0.317 | 0.592 | 0.584 | 0.654 | 0.67 | 0.687 | 0.686 |
| distinct contacts per person-day | 4.58 | 4.98 | 12.5 | 9.72 | 5.96 | 6.3 | 4.92 | 4.67 | 4.54 | 4.86 |
| CV of daily degree | 0.568 | 0.536 | 0.488 | 0.494 | 0.58 | 0.701 | 0.575 | 0.667 | 0.683 | 0.665 |
| CV of individual contact time | 1.05 | 0.899 | 0.544 | 0.559 | 0.577 | 0.601 | 0.629 | 0.728 | 0.797 | 0.784 |
| pairs repeated from previous day | 0.417 | 0.342 | 0.188 | 0.315 | 0.0811 | 0.0848 | 0.287 | 0.317 | 0.368 | 0.343 |
| within-group share of contact time | 0.817 | 0.824 | 0.247 | 0.814 | 0.239 | 0.242 | 0.589 | 0.805 | 0.877 | 0.834 |
| density exponent alpha (descriptive) | 0.452 | 0.554 | 2.08 | 1.92 | 1.97 | 2.22 | 1.88 | 1.94 | 1.81 | 1.95 |

| family | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| S0_level | pass | pass | pass | pass | pass | pass | pass | pass | pass |
| S1_duration | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S2_intercontact | pass | FAIL | FAIL | FAIL | pass | FAIL | pass | pass | pass |
| S3_degree | pass | FAIL | FAIL | FAIL | FAIL | pass | pass | pass | pass |
| S4_heterogeneity | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S5_persistence | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass | pass |
| S6_groups | pass | FAIL | pass | FAIL | FAIL | FAIL | pass | pass | pass |
| E_epidemic | n/a | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S-families passed | 6 | 1 | 2 | 1 | 2 | 2 | 4 | 5 | 5 |

| SIR scenario | quantity | real network | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|---|
| R0=1.5,tau=1.0d | R_index | 0.710 | 1.248 | 1.145 | 1.137 | 1.043 | 0.971 | 0.908 | 0.810 | 0.907 |
| R0=1.5,tau=1.0d | P(major) | 0.117 | 0.329 | 0.315 | 0.296 | 0.263 | 0.243 | 0.225 | 0.179 | 0.201 |
| R0=1.5,tau=1.0d | mean attack | 0.042 | 0.202 | 0.160 | 0.151 | 0.131 | 0.085 | 0.082 | 0.060 | 0.073 |
| R0=1.5,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=1.0d | R_index | 1.049 | 1.875 | 1.651 | 1.609 | 1.476 | 1.361 | 1.301 | 1.135 | 1.177 |
| R0=3.0,tau=1.0d | P(major) | 0.271 | 0.486 | 0.454 | 0.443 | 0.401 | 0.393 | 0.369 | 0.318 | 0.328 |
| R0=3.0,tau=1.0d | mean attack | 0.099 | 0.410 | 0.363 | 0.336 | 0.298 | 0.225 | 0.188 | 0.128 | 0.149 |
| R0=3.0,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass | FAIL |
| R0=1.5,tau=4.0d | R_index | 0.850 | 1.341 | 1.195 | 1.273 | 1.200 | 1.117 | 1.058 | 0.960 | 1.017 |
| R0=1.5,tau=4.0d | P(major) | 0.169 | 0.365 | 0.313 | 0.354 | 0.337 | 0.288 | 0.268 | 0.212 | 0.238 |
| R0=1.5,tau=4.0d | mean attack | 0.055 | 0.208 | 0.151 | 0.191 | 0.175 | 0.106 | 0.099 | 0.068 | 0.084 |
| R0=1.5,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=4.0d | R_index | 1.298 | 2.139 | 1.989 | 1.940 | 1.919 | 1.635 | 1.573 | 1.437 | 1.504 |
| R0=3.0,tau=4.0d | P(major) | 0.386 | 0.596 | 0.581 | 0.563 | 0.566 | 0.509 | 0.480 | 0.420 | 0.462 |
| R0=3.0,tau=4.0d | mean attack | 0.171 | 0.502 | 0.466 | 0.455 | 0.452 | 0.350 | 0.278 | 0.184 | 0.233 |
| R0=3.0,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass | FAIL |

**LyonSchool** — N=242, 2 days (train 1, test 1), 11 groups; calibrated c=0.00155 (M=1/c=644 sites), mean event duration 1.62 intervals, p=0.417

| statistic (test days) | observed | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|---|
| total contact (pair-intervals) | 6.52e+04 | 6.06e+04 | 5.38e+04 | 5.62e+04 | 5.34e+04 | 5.29e+04 | 5.64e+04 | 5.63e+04 | 5.59e+04 | 5.94e+04 |
| mean event duration (x20 s) | 1.62 | 1.62 | 1.62 | 1.63 | 1.67 | 1.68 | 1.63 | 1.63 | 2.18 | 1.88 |
| share of contact time in events >=5 min | 0.029 | 0.0307 | 0 | 0 | 0.000184 | 0.0197 | 0.000117 | 0 | 0.0698 | 0.0311 |
| burstiness B of gaps | 0.58 | 0.595 | -0.0471 | -0.0434 | 0.0827 | 0.163 | 0.0922 | 0.142 | 0.153 | 0.199 |
| recurrence 1-pairs/events | 0.862 | 0.842 | 0.43 | 0.706 | 0.665 | 0.643 | 0.808 | 0.853 | 0.778 | 0.884 |
| distinct contacts per person-day | 46.6 | 50 | 159 | 85 | 90.4 | 94.7 | 55.7 | 42.8 | 47.8 | 30.7 |
| CV of daily degree | 0.426 | 0.378 | 0.138 | 0.161 | 0.211 | 0.299 | 0.285 | 0.28 | 0.298 | 0.285 |
| CV of individual contact time | 0.577 | 0.502 | 0.214 | 0.273 | 0.239 | 0.244 | 0.288 | 0.35 | 0.342 | 0.537 |
| pairs repeated from previous day | 0.583 | n/a | 0.738 | 0.503 | 0.426 | 0.45 | 0.709 | 0.738 | 0.677 | 0.752 |
| within-group share of contact time | 0.729 | 0.722 | 0.0969 | 0.736 | 0.0972 | 0.0956 | 0.327 | 0.705 | 0.597 | 0.822 |
| density exponent alpha (descriptive) | 1.22 | 0.743 | 2.19 | 1.6 | 2.12 | 2.13 | 1.88 | 1.56 | 1.62 | 1.56 |

| family | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| S0_level | pass | pass | pass | pass | pass | pass | pass | pass | pass |
| S1_duration | pass | pass | pass | pass | pass | pass | pass | FAIL | pass |
| S2_intercontact | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S3_degree | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S4_heterogeneity | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| S5_persistence | n/a | FAIL | pass | FAIL | FAIL | pass | FAIL | pass | FAIL |
| S6_groups | pass | FAIL | pass | FAIL | FAIL | FAIL | pass | pass | pass |
| E_epidemic | n/a | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S-families passed | 6 | 2 | 4 | 2 | 2 | 3 | 3 | 3 | 4 |

| SIR scenario | quantity | real network | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|---|
| R0=1.5,tau=1.0d | R_index | 1.248 | 1.435 | 1.424 | 1.437 | 1.472 | 1.437 | 1.376 | 1.381 | 1.147 |
| R0=1.5,tau=1.0d | P(major) | 0.271 | 0.329 | 0.318 | 0.323 | 0.340 | 0.324 | 0.289 | 0.306 | 0.172 |
| R0=1.5,tau=1.0d | mean attack | 0.106 | 0.192 | 0.162 | 0.184 | 0.192 | 0.156 | 0.115 | 0.136 | 0.046 |
| R0=1.5,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | pass | pass | FAIL |
| R0=3.0,tau=1.0d | R_index | 1.948 | 2.554 | 2.367 | 2.574 | 2.497 | 2.351 | 2.272 | 2.291 | 1.869 |
| R0=3.0,tau=1.0d | P(major) | 0.533 | 0.603 | 0.592 | 0.619 | 0.599 | 0.609 | 0.588 | 0.577 | 0.498 |
| R0=3.0,tau=1.0d | mean attack | 0.403 | 0.552 | 0.529 | 0.564 | 0.546 | 0.539 | 0.494 | 0.497 | 0.286 |
| R0=3.0,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=1.5,tau=4.0d | R_index | 1.204 | 1.454 | 1.438 | 1.450 | 1.449 | 1.432 | 1.389 | 1.368 | 1.238 |
| R0=1.5,tau=4.0d | P(major) | 0.255 | 0.340 | 0.328 | 0.330 | 0.329 | 0.318 | 0.299 | 0.303 | 0.181 |
| R0=1.5,tau=4.0d | mean attack | 0.098 | 0.199 | 0.173 | 0.186 | 0.185 | 0.155 | 0.120 | 0.137 | 0.050 |
| R0=1.5,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=4.0d | R_index | 1.946 | 2.580 | 2.441 | 2.603 | 2.584 | 2.334 | 2.313 | 2.368 | 1.888 |
| R0=3.0,tau=4.0d | P(major) | 0.563 | 0.646 | 0.651 | 0.659 | 0.655 | 0.643 | 0.625 | 0.637 | 0.527 |
| R0=3.0,tau=4.0d | mean attack | 0.431 | 0.596 | 0.586 | 0.605 | 0.600 | 0.577 | 0.533 | 0.555 | 0.306 |
| R0=3.0,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |

**LH10** — N=75, 5 days (train 3, test 2), 4 groups; calibrated c=0.00834 (M=1/c=120 sites), mean event duration 2.26 intervals, p=0.266

| statistic (test days) | observed | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|---|
| total contact (pair-intervals) | 1.28e+04 | 1.96e+04 | 1.39e+04 | 1.33e+04 | 1.37e+04 | 1.38e+04 | 1.32e+04 | 1.36e+04 | 1.35e+04 | 1.34e+04 |
| mean event duration (x20 s) | 2.4 | 2.26 | 2.27 | 2.26 | 2.44 | 2.57 | 2.23 | 2.27 | 3.01 | 2.19 |
| share of contact time in events >=5 min | 0.14 | 0.082 | 0.00243 | 0.00314 | 0.00853 | 0.106 | 0.008 | 0.0093 | 0.171 | 0.069 |
| burstiness B of gaps | 0.612 | 0.545 | 0.119 | 0.101 | 0.25 | 0.328 | 0.26 | 0.239 | 0.247 | 0.246 |
| recurrence 1-pairs/events | 0.86 | 0.873 | 0.734 | 0.737 | 0.773 | 0.78 | 0.866 | 0.872 | 0.822 | 0.886 |
| distinct contacts per person-day | 16.1 | 16.3 | 35.3 | 33.4 | 27.6 | 25.8 | 17.2 | 16.6 | 17.3 | 15.2 |
| CV of daily degree | 0.566 | 0.571 | 0.249 | 0.269 | 0.344 | 0.377 | 0.401 | 0.41 | 0.43 | 0.421 |
| CV of individual contact time | 1.13 | 1.07 | 0.561 | 0.65 | 0.574 | 0.576 | 0.731 | 0.726 | 0.725 | 0.888 |
| pairs repeated from previous day | 0.65 | 0.49 | 0.81 | 0.785 | 0.687 | 0.646 | 0.674 | 0.663 | 0.647 | 0.692 |
| within-group share of contact time | 0.574 | 0.586 | 0.326 | 0.55 | 0.329 | 0.325 | 0.647 | 0.609 | 0.579 | 0.552 |
| density exponent alpha (descriptive) | 2.83 | 2.36 | 2.06 | 1.76 | 2.2 | 2.32 | 1.54 | 1.62 | 1.73 | 1.6 |

| family | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| S0_level | FAIL | pass | pass | pass | pass | pass | pass | pass | pass |
| S1_duration | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S2_intercontact | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S3_degree | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S4_heterogeneity | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S5_persistence | FAIL | pass | pass | pass | pass | pass | pass | pass | pass |
| S6_groups | pass | FAIL | pass | FAIL | FAIL | pass | pass | pass | pass |
| E_epidemic | n/a | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| S-families passed | 4 | 2 | 3 | 2 | 2 | 3 | 3 | 3 | 3 |

| SIR scenario | quantity | real network | WM | BLOCK | RW0 | RWhet | RWhome | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|---|
| R0=1.5,tau=1.0d | R_index | 0.957 | 1.250 | 1.224 | 1.266 | 1.205 | 1.058 | 1.089 | 1.181 | 0.946 |
| R0=1.5,tau=1.0d | P(major) | 0.297 | 0.366 | 0.367 | 0.371 | 0.362 | 0.304 | 0.325 | 0.344 | 0.275 |
| R0=1.5,tau=1.0d | mean attack | 0.143 | 0.214 | 0.208 | 0.213 | 0.209 | 0.120 | 0.145 | 0.167 | 0.106 |
| R0=1.5,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | pass | FAIL | FAIL | pass |
| R0=3.0,tau=1.0d | R_index | 1.331 | 1.950 | 1.795 | 1.908 | 1.825 | 1.525 | 1.601 | 1.609 | 1.409 |
| R0=3.0,tau=1.0d | P(major) | 0.428 | 0.526 | 0.507 | 0.528 | 0.521 | 0.479 | 0.489 | 0.487 | 0.447 |
| R0=3.0,tau=1.0d | mean attack | 0.291 | 0.443 | 0.417 | 0.440 | 0.429 | 0.327 | 0.336 | 0.361 | 0.261 |
| R0=3.0,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| R0=1.5,tau=4.0d | R_index | 0.987 | 1.283 | 1.307 | 1.311 | 1.310 | 1.093 | 1.105 | 1.139 | 0.990 |
| R0=1.5,tau=4.0d | P(major) | 0.302 | 0.383 | 0.378 | 0.384 | 0.386 | 0.318 | 0.328 | 0.346 | 0.279 |
| R0=1.5,tau=4.0d | mean attack | 0.142 | 0.218 | 0.207 | 0.215 | 0.217 | 0.122 | 0.146 | 0.169 | 0.106 |
| R0=1.5,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | pass | FAIL | FAIL | pass |
| R0=3.0,tau=4.0d | R_index | 1.522 | 2.149 | 2.011 | 2.077 | 2.003 | 1.702 | 1.628 | 1.776 | 1.505 |
| R0=3.0,tau=4.0d | P(major) | 0.505 | 0.616 | 0.585 | 0.606 | 0.602 | 0.552 | 0.529 | 0.552 | 0.498 |
| R0=3.0,tau=4.0d | mean attack | 0.336 | 0.513 | 0.472 | 0.502 | 0.494 | 0.371 | 0.362 | 0.407 | 0.284 |
| R0=3.0,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | pass | FAIL | FAIL |

### CONF summary

| dataset | model | S-families passed / applicable | failed families | E scenarios passed (of 4) | mean R_index ratio model/real | mean dP(major) | mean d(attack) |
|---|---|---|---|---|---|---|---|
| InVS15 | WM | 1/7 | S1, S2, S3, S4, S5, S6 | 0 | 1.58 | +0.186 | +0.208 |
| InVS15 | BLOCK | 2/7 | S1, S2, S3, S4, S5 | 0 | 1.43 | +0.158 | +0.164 |
| InVS15 | RW0 | 1/7 | S1, S2, S3, S4, S5, S6 | 0 | 1.46 | +0.162 | +0.177 |
| InVS15 | RWhet | 2/7 | S2, S3, S4, S5, S6 | 0 | 1.39 | +0.143 | +0.157 |
| InVS15 | RWclus | 4/7 | S1, S2, S5 | 0 | 1.16 | +0.059 | +0.006 |
| InVS15 | RWstickZ | 5/7 | S2, S5 | 0 | 1.12 | +0.044 | -0.004 |
| InVS15 | RWstickO | 4/7 | S1, S2, S5 | 2 | 0.97 | -0.034 | -0.059 |
| Thiers13 | WM | 1/7 | S1, S2, S3, S4, S5, S6 | 0 | 1.79 | +0.322 | +0.319 |
| Thiers13 | BLOCK | 2/7 | S1, S2, S3, S4, S5 | 0 | 1.63 | +0.263 | +0.234 |
| Thiers13 | RW0 | 1/7 | S1, S2, S3, S4, S5, S6 | 0 | 1.72 | +0.303 | +0.292 |
| Thiers13 | RWhet | 2/7 | S2, S3, S4, S5, S6 | 0 | 1.64 | +0.289 | +0.276 |
| Thiers13 | RWclus | 2/7 | S1, S2, S3, S4, S5 | 0 | 1.31 | +0.021 | -0.000 |
| Thiers13 | RWstickZ | 2/7 | S1, S2, S3, S4, S5 | 0 | 1.30 | +0.038 | +0.006 |
| Thiers13 | RWstickO | 4/7 | S0, S2, S5 | 0 | 1.26 | -0.076 | -0.014 |
| SFHH | WM | 2/6 | S1, S2, S3, S4 | 0 | 1.50 | +0.120 | +0.157 |
| SFHH | BLOCK | 2/6 | S1, S2, S3, S4 | 0 | 1.48 | +0.119 | +0.157 |
| SFHH | RW0 | 1/6 | S1, S2, S3, S4, S5 | 0 | 1.35 | +0.078 | +0.103 |
| SFHH | RWhet | 2/6 | S1, S2, S4, S5 | 0 | 1.24 | +0.016 | +0.027 |
| SFHH | RWclus | 1/6 | S1, S2, S3, S4, S5 | 0 | 1.28 | +0.053 | +0.070 |
| SFHH | RWstickZ | 2/6 | S2, S3, S4, S5 | 0 | 1.21 | +0.019 | +0.037 |
| SFHH | RWstickO | 2/6 | S2, S3, S4, S5 | 0 | 1.13 | -0.008 | +0.019 |

### CONF detail

**InVS15** — N=217, 10 days (train 5, test 5), 12 groups; calibrated c=0.00072 (M=1/c=1388 sites), mean event duration 2.27 intervals, p=0.264

| statistic (test days) | observed | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| total contact (pair-intervals) | 3.47e+04 | 4.35e+04 | 3.49e+04 | 3.47e+04 | 3.44e+04 | 3.49e+04 | 3.74e+04 | 3.78e+04 | 3.97e+04 |
| mean event duration (x20 s) | 2.39 | 2.27 | 2.27 | 2.26 | 2.31 | 2.28 | 2.26 | 2.71 | 2.51 |
| share of contact time in events >=5 min | 0.163 | 0.168 | 0.00258 | 0.00283 | 0.00336 | 0.173 | 0.00473 | 0.177 | 0.0807 |
| burstiness B of gaps | 0.447 | 0.453 | 0.0252 | 0.0847 | 0.205 | 0.319 | 0.292 | 0.269 | 0.312 |
| recurrence 1-pairs/events | 0.775 | 0.783 | 0.149 | 0.398 | 0.607 | 0.603 | 0.772 | 0.739 | 0.814 |
| distinct contacts per person-day | 8.6 | 10.1 | 34.7 | 24.4 | 15.5 | 16.2 | 9.95 | 9.62 | 7.8 |
| CV of daily degree | 0.65 | 0.73 | 0.428 | 0.494 | 0.467 | 0.703 | 0.675 | 0.631 | 0.716 |
| CV of individual contact time | 0.821 | 0.832 | 0.417 | 0.602 | 0.444 | 0.477 | 0.779 | 0.784 | 0.95 |
| pairs repeated from previous day | 0.337 | 0.399 | 0.22 | 0.43 | 0.101 | 0.102 | 0.548 | 0.502 | 0.584 |
| within-group share of contact time | 0.797 | 0.735 | 0.136 | 0.734 | 0.135 | 0.139 | 0.65 | 0.653 | 0.672 |
| density exponent alpha (descriptive) | 0.947 | 0.471 | 2.05 | 2.02 | 2.06 | 2.05 | 2.01 | 2.04 | 2.02 |

| family | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|
| S0_level | FAIL | pass | pass | pass | pass | pass | pass | pass |
| S1_duration | pass | FAIL | FAIL | FAIL | pass | FAIL | pass | FAIL |
| S2_intercontact | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S3_degree | pass | FAIL | FAIL | FAIL | FAIL | pass | pass | pass |
| S4_heterogeneity | pass | FAIL | FAIL | FAIL | FAIL | pass | pass | pass |
| S5_persistence | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S6_groups | pass | FAIL | pass | FAIL | FAIL | pass | pass | pass |
| E_epidemic | n/a | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S-families passed | 6 | 1 | 2 | 1 | 2 | 4 | 5 | 4 |

| SIR scenario | quantity | real network | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| R0=1.5,tau=1.0d | R_index | 0.855 | 1.392 | 1.242 | 1.288 | 1.153 | 1.043 | 0.986 | 0.907 |
| R0=1.5,tau=1.0d | P(major) | 0.116 | 0.333 | 0.293 | 0.300 | 0.257 | 0.200 | 0.177 | 0.126 |
| R0=1.5,tau=1.0d | mean attack | 0.038 | 0.209 | 0.163 | 0.166 | 0.136 | 0.064 | 0.057 | 0.038 |
| R0=1.5,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| R0=3.0,tau=1.0d | R_index | 1.261 | 2.171 | 1.975 | 1.978 | 1.868 | 1.449 | 1.415 | 1.185 |
| R0=3.0,tau=1.0d | P(major) | 0.306 | 0.489 | 0.480 | 0.478 | 0.445 | 0.365 | 0.355 | 0.268 |
| R0=3.0,tau=1.0d | mean attack | 0.167 | 0.433 | 0.403 | 0.406 | 0.370 | 0.180 | 0.171 | 0.099 |
| R0=3.0,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=1.5,tau=4.0d | R_index | 1.020 | 1.441 | 1.361 | 1.363 | 1.337 | 1.152 | 1.175 | 1.021 |
| R0=1.5,tau=4.0d | P(major) | 0.165 | 0.351 | 0.326 | 0.329 | 0.321 | 0.226 | 0.236 | 0.154 |
| R0=1.5,tau=4.0d | mean attack | 0.054 | 0.208 | 0.170 | 0.186 | 0.175 | 0.073 | 0.075 | 0.045 |
| R0=1.5,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| R0=3.0,tau=4.0d | R_index | 1.588 | 2.450 | 2.177 | 2.260 | 2.283 | 1.822 | 1.683 | 1.422 |
| R0=3.0,tau=4.0d | P(major) | 0.456 | 0.616 | 0.575 | 0.586 | 0.592 | 0.489 | 0.452 | 0.360 |
| R0=3.0,tau=4.0d | mean attack | 0.297 | 0.540 | 0.474 | 0.506 | 0.504 | 0.261 | 0.235 | 0.136 |
| R0=3.0,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |

**Thiers13** — N=327, 5 days (train 3, test 2), 9 groups; calibrated c=0.000918 (M=1/c=1089 sites), mean event duration 2.79 intervals, p=0.206

| statistic (test days) | observed | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| total contact (pair-intervals) | 7.22e+04 | 1.16e+05 | 7.14e+04 | 7.25e+04 | 7.11e+04 | 7.18e+04 | 7.39e+04 | 7.26e+04 | 1.61e+05 |
| mean event duration (x20 s) | 2.79 | 2.79 | 2.78 | 2.79 | 2.86 | 2.84 | 2.82 | 2.98 | 2.8 |
| share of contact time in events >=5 min | 0.243 | 0.231 | 0.0111 | 0.0116 | 0.0171 | 0.251 | 0.0182 | 0.147 | 0.251 |
| burstiness B of gaps | 0.507 | 0.513 | -0.000351 | 0.0145 | 0.153 | 0.33 | 0.267 | 0.218 | 0.385 |
| recurrence 1-pairs/events | 0.836 | 0.833 | 0.155 | 0.618 | 0.607 | 0.619 | 0.848 | 0.82 | 0.923 |
| distinct contacts per person-day | 14.6 | 15.3 | 74.7 | 34.2 | 33.8 | 33.4 | 13.7 | 15.1 | 15.2 |
| CV of daily degree | 0.518 | 0.488 | 0.261 | 0.199 | 0.308 | 0.635 | 0.382 | 0.331 | 0.421 |
| CV of individual contact time | 0.926 | 0.776 | 0.321 | 0.362 | 0.344 | 0.387 | 0.635 | 0.618 | 1.09 |
| pairs repeated from previous day | 0.461 | 0.471 | 0.257 | 0.714 | 0.116 | 0.115 | 0.692 | 0.688 | 0.78 |
| within-group share of contact time | 0.943 | 0.93 | 0.116 | 0.931 | 0.119 | 0.117 | 0.937 | 0.934 | 1 |
| density exponent alpha (descriptive) | 0.759 | 0.719 | 2.06 | 1.86 | 2.12 | 2.09 | 1.8 | 1.86 | 1.92 |

| family | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|
| S0_level | FAIL | pass | pass | pass | pass | pass | pass | FAIL |
| S1_duration | pass | FAIL | FAIL | FAIL | pass | FAIL | FAIL | pass |
| S2_intercontact | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S3_degree | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| S4_heterogeneity | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | pass |
| S5_persistence | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S6_groups | pass | FAIL | pass | FAIL | FAIL | pass | pass | pass |
| E_epidemic | n/a | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S-families passed | 6 | 1 | 2 | 1 | 2 | 2 | 2 | 4 |

| SIR scenario | quantity | real network | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| R0=1.5,tau=1.0d | R_index | 0.828 | 1.404 | 1.371 | 1.391 | 1.333 | 1.086 | 1.117 | 1.113 |
| R0=1.5,tau=1.0d | P(major) | 0.023 | 0.327 | 0.232 | 0.308 | 0.278 | 0.034 | 0.041 | 0.015 |
| R0=1.5,tau=1.0d | mean attack | 0.016 | 0.194 | 0.076 | 0.163 | 0.142 | 0.023 | 0.025 | 0.025 |
| R0=1.5,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=1.0d | R_index | 1.239 | 2.463 | 2.161 | 2.304 | 2.233 | 1.639 | 1.634 | 1.453 |
| R0=3.0,tau=1.0d | P(major) | 0.234 | 0.568 | 0.558 | 0.553 | 0.549 | 0.274 | 0.301 | 0.115 |
| R0=3.0,tau=1.0d | mean attack | 0.071 | 0.512 | 0.468 | 0.488 | 0.476 | 0.069 | 0.078 | 0.043 |
| R0=3.0,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=1.5,tau=4.0d | R_index | 0.891 | 1.413 | 1.340 | 1.389 | 1.338 | 1.177 | 1.136 | 1.193 |
| R0=1.5,tau=4.0d | P(major) | 0.022 | 0.336 | 0.226 | 0.308 | 0.301 | 0.044 | 0.042 | 0.019 |
| R0=1.5,tau=4.0d | mean attack | 0.017 | 0.197 | 0.075 | 0.166 | 0.158 | 0.025 | 0.025 | 0.026 |
| R0=1.5,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=4.0d | R_index | 1.380 | 2.625 | 2.237 | 2.433 | 2.299 | 1.756 | 1.725 | 1.636 |
| R0=3.0,tau=4.0d | P(major) | 0.298 | 0.635 | 0.613 | 0.621 | 0.603 | 0.311 | 0.346 | 0.126 |
| R0=3.0,tau=4.0d | mean attack | 0.094 | 0.571 | 0.513 | 0.550 | 0.526 | 0.080 | 0.093 | 0.049 |
| R0=3.0,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |

**SFHH** — N=403, 2 days (train 1, test 1), 0 groups; calibrated c=0.00082 (M=1/c=1219 sites), mean event duration 2.82 intervals, p=0.203

| statistic (test days) | observed | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| total contact (pair-intervals) | 2.45e+04 | 4.58e+04 | 2.7e+04 | 2.71e+04 | 2.7e+04 | 2.67e+04 | 2.76e+04 | 2.65e+04 | 2.75e+04 |
| mean event duration (x20 s) | 2.49 | 2.82 | 2.83 | 2.82 | 2.89 | 3.11 | 2.79 | 2.84 | 2.71 |
| share of contact time in events >=5 min | 0.221 | 0.252 | 0.01 | 0.0108 | 0.0198 | 0.282 | 0.0125 | 0.246 | 0.253 |
| burstiness B of gaps | 0.499 | 0.522 | 0.0203 | 0.0201 | 0.167 | 0.306 | 0.165 | 0.243 | 0.235 |
| recurrence 1-pairs/events | 0.536 | 0.641 | 0.0918 | 0.092 | 0.574 | 0.575 | 0.63 | 0.667 | 0.62 |
| distinct contacts per person-day | 26.5 | 29.9 | 51 | 51.3 | 23.4 | 21.6 | 21.5 | 18.3 | 22.7 |
| CV of daily degree | 0.666 | 0.741 | 0.322 | 0.323 | 0.36 | 0.71 | 0.379 | 0.396 | 0.488 |
| CV of individual contact time | 1.13 | 1.24 | 0.446 | 0.443 | 0.496 | 0.561 | 0.499 | 0.556 | 0.887 |
| pairs repeated from previous day | 0.18 | n/a | 0.183 | 0.182 | 0.0812 | 0.0768 | 0.368 | 0.35 | 0.323 |
| within-group share of contact time | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| density exponent alpha (descriptive) | 0.805 | 1.25 | 1.99 | 2.1 | 1.95 | 2.08 | 2.08 | 1.94 | 2.19 |

| family | train days (ceiling) | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|
| S0_level | FAIL | pass | pass | pass | pass | pass | pass | pass |
| S1_duration | pass | FAIL | FAIL | FAIL | FAIL | FAIL | pass | pass |
| S2_intercontact | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S3_degree | pass | FAIL | FAIL | FAIL | pass | FAIL | FAIL | FAIL |
| S4_heterogeneity | pass | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S5_persistence | n/a | pass | pass | FAIL | FAIL | FAIL | FAIL | FAIL |
| S6_groups | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| E_epidemic | n/a | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| S-families passed | 4 | 2 | 2 | 1 | 2 | 1 | 2 | 2 |

| SIR scenario | quantity | real network | WM | BLOCK | RW0 | RWhet | RWclus | RWstickZ | RWstickO |
|---|---|---|---|---|---|---|---|---|---|
| R0=1.5,tau=1.0d | R_index | 0.952 | 1.406 | 1.404 | 1.310 | 1.198 | 1.292 | 1.194 | 1.072 |
| R0=1.5,tau=1.0d | P(major) | 0.215 | 0.326 | 0.327 | 0.256 | 0.205 | 0.236 | 0.184 | 0.156 |
| R0=1.5,tau=1.0d | mean attack | 0.076 | 0.186 | 0.186 | 0.110 | 0.073 | 0.078 | 0.055 | 0.048 |
| R0=1.5,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=1.0d | R_index | 1.505 | 2.383 | 2.405 | 2.070 | 1.852 | 1.917 | 1.819 | 1.741 |
| R0=3.0,tau=1.0d | P(major) | 0.428 | 0.552 | 0.555 | 0.538 | 0.465 | 0.514 | 0.496 | 0.476 |
| R0=3.0,tau=1.0d | mean attack | 0.269 | 0.463 | 0.465 | 0.432 | 0.313 | 0.394 | 0.354 | 0.333 |
| R0=3.0,tau=1.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=1.5,tau=4.0d | R_index | 0.999 | 1.402 | 1.371 | 1.329 | 1.278 | 1.258 | 1.234 | 1.143 |
| R0=1.5,tau=4.0d | P(major) | 0.229 | 0.342 | 0.333 | 0.281 | 0.224 | 0.236 | 0.196 | 0.165 |
| R0=1.5,tau=4.0d | mean attack | 0.081 | 0.195 | 0.192 | 0.122 | 0.081 | 0.082 | 0.060 | 0.052 |
| R0=1.5,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| R0=3.0,tau=4.0d | R_index | 1.623 | 2.467 | 2.408 | 2.119 | 1.920 | 2.006 | 1.850 | 1.787 |
| R0=3.0,tau=4.0d | P(major) | 0.480 | 0.612 | 0.610 | 0.589 | 0.521 | 0.578 | 0.552 | 0.522 |
| R0=3.0,tau=4.0d | mean attack | 0.304 | 0.515 | 0.514 | 0.478 | 0.369 | 0.454 | 0.411 | 0.373 |
| R0=3.0,tau=4.0d | scenario passes | | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |


#!/bin/sh
# Reproduce the validation.  Order matters; seeds are fixed (master seed 20261001).
set -e
cd "$(dirname "$0")"
PY=python
export TMPDIR="$PWD/.tmp" MPLCONFIGDIR="$PWD/.mplcache"
mkdir -p .tmp .mplcache logs data figures
$PY scripts/00_selftest.py                 # synthetic inputs only
$PY scripts/01_calibrate.py                # calibration records only (T4, R1), all variants
$PY scripts/02_power.py                    # hold-out predictive distributions + power (no stratum outcomes read)
$PY scripts/02b_m1_exact.py                # pre-registered M1 check against the exact simulator (../bsc_sim)
$PY scripts/03_holdout_shape.py            # first script that reads data/heldout_outcomes.json
$PY scripts/02c_mc_stability.py
$PY scripts/04_level_controls.py
$PY scripts/06_callcentre.py
$PY scripts/07_loo.py
$PY scripts/08_bruteforce_check.py
$PY scripts/09_sensitivity.py
$PY scripts/10_decision.py
$PY scripts/20_figures.py
$PY scripts/30_tables.py

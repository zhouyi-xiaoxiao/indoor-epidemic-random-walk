#!/bin/sh
# Re-runs the whole analysis from the tracer tables. Each step < 10 min, < 1 GB. Overwrites results/ (time stamps change).
cd "$(dirname "$0")/src" || exit 1
PY=python
for s in s01_kernels s02_predict s03_power s04_evaluate s05_descriptive e01_ext_predict e02_ext_evaluate; do $PY $s.py || exit 1; done
(cd ../../../checks/tracer && $PY v01_checks.py && $PY v03_ext_checks.py) || exit 1

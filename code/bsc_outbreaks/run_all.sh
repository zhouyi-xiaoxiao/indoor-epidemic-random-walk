#!/bin/sh
# Rebuild everything except the network steps (DOI check, source downloads). Deterministic.
set -e
cd "$(dirname "$0")"
PY=python
$PY scripts/digitize_callcentre.py
$PY scripts/build_dataset.py
$PY scripts/analyze_dataset.py > notes/analyze_stdout.txt

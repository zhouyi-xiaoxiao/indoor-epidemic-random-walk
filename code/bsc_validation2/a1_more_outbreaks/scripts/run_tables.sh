#!/bin/bash
# usage: run_tables.sh EVENT...   (primary tables)
cd "$(dirname "$0")/.."
PY=python
for ev in "$@"; do $PY scripts/01_tables.py $ev >> out/01_tables.log 2>&1; done
echo done "$@" >> out/01_tables.log

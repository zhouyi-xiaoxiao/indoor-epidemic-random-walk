#!/bin/bash
cd "$(dirname "$0")/.."
PY=python
run() { $PY scripts/01_tables.py "$1" "$2" "$3" >> out/01_tables_variants.log 2>&1; }
for e in F5 F6 F3 F7; do run $e '{"outcome":"wide"}' wide; done
for e in F6 F3 F7; do run $e '{"src":"alt1"}' src1; done
run F6 '{"src":"alt2"}' src2
run F5 '{"hh":true}' hh
run F1 '{"interviewed":true}' int
for e in F1 F2 F3 F4 F5 F6 F7 F8; do run $e '{"ay":0.75}' ay075; run $e '{"ay":0.86}' ay086; run $e '{"cabin_ach":[5,15]}' ach5; done
echo ALLDONE >> out/01_tables_variants.log

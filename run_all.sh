#!/usr/bin/env bash
# Reproduce everything in output/ from a clean checkout.
# Tests run first: if the measurement logic is broken there is no point fitting
# 100 models to find out.
set -euo pipefail
cd "$(dirname "$0")"

R="${RSCRIPT:-Rscript}"
mkdir -p output

echo "### tests"
"$R" tests.R

echo
echo "### 00 download and verify vintage"
"$R" 00_download.R

echo
echo "### 01 audit"
"$R" 01_audit.R | tee output/01_audit.log

echo
echo "### 02 invariance sequence, partial invariance, consequence"
"$R" 02_invariance.R | tee output/02_invariance.log

echo
echo "### 03 negative controls (takes about a minute, fits ~130 models)"
"$R" 03_negative_control.R | tee output/03_negative_control.log

echo
echo "### 04 figures and tables"
"$R" 04_report.R

echo
echo "### 05 graded response model and DIF by sex"
"$R" 05_grm_dif.R | tee output/05_grm_dif.log

echo
echo "### 06 DIF error rate (the slow one: ~700 multigroup GRM fits, parallel)"
echo "###    set N_SPLITS / N_REPS / N_CORES to shrink it for a smoke run"
"$R" 06_dif_error_rate.R | tee output/06_dif_error_rate.log

echo
echo "### 07 IRT vs CFA comparison"
"$R" 07_compare.R | tee output/07_compare.log

echo
echo "### 09 local dependence, and whether the DIF flags survive it"
"$R" 09_local_dependence.R | tee output/09_local_dependence.log

echo
echo "### 08 IRT figures and tables"
"$R" 08_dif_figures.R

echo
echo "done. output/:"
ls -la output/

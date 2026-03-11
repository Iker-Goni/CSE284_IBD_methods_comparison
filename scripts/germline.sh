#!/bin/bash
# Script 2: GERMLINE IBD detection
set -e
set -o pipefail

PHASE_PREFIX="lwk_phased"
GERMLINE_OUT="lwk_ibd"

echo "=== Step 1: Run GERMLINE ==="
../germline/germline \
  -input ${PHASE_PREFIX}.ped ${PHASE_PREFIX}.map \
  -output ${GERMLINE_OUT} \
  -bits 64 \   # adjust bits if you have memory issues
  -min_m 3

echo "=== Step 2: Calculate average segment length (MB) excluding self-matches ==="
awk '$1 != $3 {sum += ($7-$6)/1e6; n++} END {print "Average segment length (MB):", sum/n}' ${GERMLINE_OUT}.match

echo "GERMLINE IBD detection completed. Matches in ${GERMLINE_OUT}.match"

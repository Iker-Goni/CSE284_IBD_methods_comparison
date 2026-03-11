#!/bin/bash
# Script 2: GERMLINE IBD detection
set -e
set -o pipefail

PHASE_PREFIX="data/lwk_phased"
GERMLINE_OUT="data/lwk_ibd"

# adjust bits if you have memory issues
echo "=== Step 1: Run GERMLINE ==="
tools/germline/germline \
  -input ${PHASE_PREFIX}.ped ${PHASE_PREFIX}.map \
  -output ${GERMLINE_OUT} \
  -bits 8 \
  -min_m 3

echo "=== Step 2: Calculate average segment length (MB) excluding self-matches ==="
awk '$1 != $3 {sum += ($7-$6)/1e6; n++} END {print "Average segment length (MB):", sum/n}' ${GERMLINE_OUT}.match

echo "GERMLINE IBD detection completed. Matches in ${GERMLINE_OUT}.match"

#!/bin/bash
# Script 2: GERMLINE IBD detection per chromosome
set -e
set -o pipefail

PHASE_PREFIX="data/lwk_phased"
GERMLINE_OUT_PREFIX="data/lwk_ibd"
BITS=8
MIN_M=3

# Create output directory if it doesn't exist
mkdir -p data

echo "=== Running GERMLINE on each chromosome ==="
for chr in {1..22}; do
    echo "Processing chromosome $chr..."
    
    # Input PED/MAP for this chromosome
    PED_FILE="${PHASE_PREFIX}_chr${chr}.ped"
    MAP_FILE="${PHASE_PREFIX}_chr${chr}.map"
    
    # Output prefix for GERMLINE
    OUT_FILE="${GERMLINE_OUT_PREFIX}_chr${chr}"
    
    # Run GERMLINE
    tools/germline/germline \
      -input $PED_FILE $MAP_FILE \
      -output $OUT_FILE \
      -bits $BITS \
      -min_m $MIN_M
done

echo "=== Combine all chromosome match files ==="
COMBINED_MATCH="${GERMLINE_OUT_PREFIX}_allchr.match"
cat ${GERMLINE_OUT_PREFIX}_chr*.match > $COMBINED_MATCH
echo "Combined match file created: $COMBINED_MATCH"

echo "=== Calculate average segment length (MB) excluding self-matches ==="
awk '$1 != $3 {sum += ($7-$6)/1e6; n++} END {print "Average segment length (MB):", sum/n}' $COMBINED_MATCH

echo "GERMLINE IBD detection completed."

#!/bin/bash

# Exit if any command fails
set -e

# Input prefix (no extension)
INPUT_PREFIX="../data/ps2_ibd.lwk"

# Output prefix
OUTPUT_PREFIX="../results/plink_ibd"

echo "Running PLINK IBD calculation..."
echo "Input:  $INPUT_PREFIX"
echo "Output: $OUTPUT_PREFIX"

# Run PLINK genome-wide IBD estimation
plink \
  --bfile $INPUT_PREFIX \
  --genome full \
  --out $OUTPUT_PREFIX

echo "PLINK IBD computation complete."
echo "Output file: ${OUTPUT_PREFIX}.genome"

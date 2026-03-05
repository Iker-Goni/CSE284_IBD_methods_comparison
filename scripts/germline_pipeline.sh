#!/bin/bash

# Step 1: Convert VCF to PED/MAP (GERMLINE-friendly)
plink \
  --vcf data/ps2_ibd_phased.vcf.gz \
  --biallelic-only strict \
  --geno 0 \
  --snps-only just-acgt \
  --keep-allele-order \
  --recode ped \
  --out data/germline_input

# Step 2: Run GERMLINE on converted data
./tools/germline/germline \
  -input data/germline_input.ped data/germline_input.map \
  -output results/germline_out \
  -min_m 3 \
  -bits 16

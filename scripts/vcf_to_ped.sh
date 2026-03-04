#!/bin/bash

for chr in {1..22}; do
    ../../plink_2/plink2 \
        --vcf ../data/ps2_ibd_phased.vcf.gz \
        --chr ${chr} \
        --recode ped \
        --out ../data/germline/peds/ps2_ibd_chr${chr}_germline \
        --keep-allele-order \
        --snps-only just-acgt \
        --max-alleles 2 \
        --min-alleles 2 \
        --geno 0 \
        --mind 0 \
        --allow-extra-chr

done
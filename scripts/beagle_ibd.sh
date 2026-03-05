#!/bin/bash

for chr in {1..22}; do
  java -jar tools/beagle/beagle.27Jan18.7e1.jar \
    gt=data/ps2_ibd_phased.vcf.gz \
    chrom=$chr \
    map=data/maps/chr${chr}.map \
    out=data/beagle/ps2_ibd_chr${chr} \
    ibd=true
done
#!/bin/bash

# This is code for if the vcf -> ped conversion actually works, which it isn't rn
for chr in {1..22}; do
    ~/Project/germline/germline \
        -input ../data/germline/peds/ps2_ibd_chr${chr}_germline.ped \
               ../data/germline/peds/ps2_ibd_chr${chr}_germline.map \
        -output ../data/germline/peds/ps2_ibd_chr${chr} \
        -map ../data/maps/chr${chr}.map
done
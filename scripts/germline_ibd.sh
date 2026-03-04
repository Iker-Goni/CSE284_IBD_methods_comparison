#!/bin/bash


for chr in {21..22}; do
    ~/Project/germline/germline \
        -input ../data/germline/peds/ps2_ibd_chr${chr}_germline.ped \
               ../data/germline/peds/ps2_ibd_chr${chr}_germline.map \
        -output ../data/germline/peds/ps2_ibd_chr${chr} \
        -map ../data/maps/chr${chr}.map
done
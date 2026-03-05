# CSE284_IBD_methods_comparison

## Overview
The goal of this project is to compare the IBD analysis of 3 tools: plink 1.9, germline, and Beagle. Plink utilizes a Method of Moments estimator to compute the probability that a pair of individuals shares 0, 1, or 2 alleles IBD across the genome. Germline uses a hashing based technique to run in linear time with respect to the number of individuals. Beagle uses a refined IBD algorithm, which utilizes a HMM and computes a log odds score (LOD score) of how likely a segment is to be in IBD. 

## Data
The dataset that we will be analyzing is the 1000 Genomes Phase 3 release, consisting of 2504 individuals from 26 different populations. Our analysis will focus on the subset of individuals from the LWK population. This dataset uses the GRCh37 reference genome and includes VCF files for every chromosome containing genotype information regarding variants for every individual. Individuals were sequenced with whole-genome sequencing with a mean depth of 7.4x and targeted exome sequencing with a mean depth of 65.7x. 

## Repository structure

## Dependencies

## Instructions to reproduce results


## Results so far

## Remaining work to complete
There is currently a large discrepancy between the total number of segments, mean segment length, median segment length, and max segment length computed by germline as opposed to Beagle. We believe it is an issue with the input .ped file and will work to correct it. Some other future steps include:
- Comparing runtime of each method
- Estimating pi_hat from germline and Beagle results and comparing to plink

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
import glob
import os
import seaborn as sns
from collections import defaultdict


BEAGLE_DIR = "../data.beagle"
OUTPUT_DIR = "../data/beagle_analysis"
MIN_LOD = 3.0
os.makedirs(OUTPUT_DIR, exist_ok=True)


colnames = ['sample1', 'hap1', 'sample2', 'hap2', 'chrom', 'start', 'end', 'lod']
all_segments = []

file_pattern = os.path.join(BEAGLE_DIR, "ps2_ibd_chr*.ibd.gz")
file_list = sorted(glob.glob(file_pattern))

if not file_list:
    print(f"Error: No files found in {BEAGLE_DIR}")
    exit(1)

print(f"Found IBD files for {len(file_list)} chromosomes")

for file_path in file_list:
    chromosome = file_path.split('chr')[-1].split('.')[0]
    
# Unfinished script
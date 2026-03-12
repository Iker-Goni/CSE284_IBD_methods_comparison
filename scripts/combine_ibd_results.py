#!/usr/bin/env python3
import csv
from collections import defaultdict

# -------- FILE PATHS --------
plink_file = "results/plink_ibd.genome"
germline_file = "data/lwk_ibd_allchr.match"
beagle_file = "results/beagle/pairwise_pi_hat.csv"
output_file = "results/combined_ibd_table.csv"

# -------- DATA STRUCTURES --------
plink_ibd = {}
germline_ibd = defaultdict(float)
beagle_ibd = {}

# -------- READ PLINK --------
with open(plink_file) as f:
    next(f)  # skip header
    for line in f:
        parts = line.split()
        id1 = parts[1]
        id2 = parts[3]
        pihat = float(parts[9])
        pair = tuple(sorted([id1, id2]))
        plink_ibd[pair] = pihat

# ---------------- GERMLINE ----------------
with open(germline_file) as f:
    for line in f:
        parts = line.split()

        id1 = parts[0]
        id2 = parts[2]

        start = int(parts[5])
        end = int(parts[6])

        length_mb = (end - start) / 1e6

        pair = tuple(sorted([id1, id2]))
        germline_ibd[pair] += length_mb

GENOME_MB = 3200.0
for pair in germline_ibd:
    germline_ibd[pair] /= GENOME_MB

# ---------------- BEAGLE ----------------
with open(beagle_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        id1 = row["sample_a"].split("_")[0]
        id2 = row["sample_b"].split("_")[0]
        pihat = float(row["pi_hat"])
        pair = tuple(sorted([id1, id2]))
        beagle_ibd[pair] = pihat

# -------- COMBINE PAIRS --------
all_pairs = set(plink_ibd) | set(germline_ibd) | set(beagle_ibd)

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["ID1", "ID2", "PLINK_PIHAT", "GERMLINE_IBD", "BEAGLE_PIHAT"])

    for id1, id2 in sorted(all_pairs):
        writer.writerow([
            id1,
            id2,
            plink_ibd.get((id1, id2), 0),
            germline_ibd.get((id1, id2), 0),
            beagle_ibd.get((id1, id2), 0)
        ])

print("Combined table written to:", output_file)

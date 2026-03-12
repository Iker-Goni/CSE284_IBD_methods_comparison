#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---------------- LOAD DATA ----------------
df = pd.read_csv("results/combined_ibd_table.csv")

# ---------------- APPLY THRESHOLD ----------------
threshold = 0.01
# Only keep rows where all tools are above threshold
df_filtered = df[(df["PLINK_PIHAT"] >= threshold) & 
                 (df["GERMLINE_IBD"] >= threshold) & 
                 (df["BEAGLE_PIHAT"] >= threshold)]

print(f"Number of pairs above threshold {threshold}: {len(df_filtered)}")

# ---------------- SCATTER PLOTS ----------------
plt.figure(figsize=(18, 5))  # wider figure for 3 plots

# Scatter: PLINK vs GERMLINE
plt.subplot(1, 3, 1)
plt.scatter(df_filtered["PLINK_PIHAT"], df_filtered["GERMLINE_IBD"], alpha=0.5, label="Pairs")
plt.plot([0, max(df_filtered["PLINK_PIHAT"].max(), df_filtered["GERMLINE_IBD"].max())],
         [0, max(df_filtered["PLINK_PIHAT"].max(), df_filtered["GERMLINE_IBD"].max())],
         color='red', linestyle='--', label="y = x")
plt.xlabel("PLINK IBD")
plt.ylabel("GERMLINE IBD")
plt.title("PLINK vs GERMLINE (thresholded)")
plt.legend()

# Scatter: PLINK vs BEAGLE
plt.subplot(1, 3, 2)
plt.scatter(df_filtered["PLINK_PIHAT"], df_filtered["BEAGLE_PIHAT"], alpha=0.5, label="Pairs")
plt.plot([0, max(df_filtered["PLINK_PIHAT"].max(), df_filtered["BEAGLE_PIHAT"].max())],
         [0, max(df_filtered["PLINK_PIHAT"].max(), df_filtered["BEAGLE_PIHAT"].max())],
         color='red', linestyle='--', label="y = x")
plt.xlabel("PLINK IBD")
plt.ylabel("BEAGLE IBD")
plt.title("PLINK vs BEAGLE (thresholded)")
plt.legend()

# Scatter: BEAGLE vs GERMLINE
plt.subplot(1, 3, 3)
plt.scatter(df_filtered["BEAGLE_PIHAT"], df_filtered["GERMLINE_IBD"], alpha=0.5, label="Pairs")
plt.plot([0, max(df_filtered["BEAGLE_PIHAT"].max(), df_filtered["GERMLINE_IBD"].max())],
         [0, max(df_filtered["BEAGLE_PIHAT"].max(), df_filtered["GERMLINE_IBD"].max())],
         color='red', linestyle='--', label="y = x")
plt.xlabel("BEAGLE IBD")
plt.ylabel("GERMLINE IBD")
plt.title("BEAGLE vs GERMLINE (thresholded)")
plt.legend()

plt.tight_layout()
plt.savefig("results/ibd_scatter_plots_3way_thresholded.png", dpi=300)
plt.show()

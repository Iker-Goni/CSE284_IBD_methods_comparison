import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../results/germline_out.match", sep=r"\s+", header=None)

df.columns = [
    "id1","hap1","id2","hap2","chr","start","end",
    "snp1","snp2","nsnp","length_mb","unit",
    "err1","err2","err3"
]

segment_lengths = df["length_mb"]

plt.figure()
plt.hist(segment_lengths, bins=50)
plt.xlabel("IBD Segment Length (MB)")
plt.ylabel("Count")
plt.title("GERMLINE IBD Segment Length Distribution")
plt.savefig("germline_segment_hist.png", dpi=300)

pair_ibd = df.groupby(["id1","id2"])["length_mb"].sum()

plt.figure()
plt.hist(pair_ibd, bins=40)
plt.xlabel("Total Shared IBD Length (MB)")
plt.ylabel("Number of Pairs")
plt.title("GERMLINE Total Shared IBD per Pair")
plt.savefig("germline_pair_hist.png", dpi=300)

plt.show()

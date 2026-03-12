import pandas as pd

df = pd.read_csv("results/combined_ibd_table.csv")

# compute pairwise differences
df["plink_vs_germline"] = abs(df["PLINK_PIHAT"] - df["GERMLINE_IBD"])
df["plink_vs_beagle"] = abs(df["PLINK_PIHAT"] - df["BEAGLE_PIHAT"])
df["germline_vs_beagle"] = abs(df["GERMLINE_IBD"] - df["BEAGLE_PIHAT"])

# compute averages
avg_pg = df["plink_vs_germline"].mean()
avg_pb = df["plink_vs_beagle"].mean()
avg_gb = df["germline_vs_beagle"].mean()

print("\nAverage IBD differences between tools:\n")
print("PLINK vs GERMLINE :", avg_pg)
print("PLINK vs BEAGLE   :", avg_pb)
print("GERMLINE vs BEAGLE:", avg_gb)

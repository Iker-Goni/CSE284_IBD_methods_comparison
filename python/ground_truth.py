import pandas as pd

# Loading the PLINK genome output
df = pd.read_csv("../results/plink_ibd.genome", delim_whitespace=True)

threshold = 0.125

df["relationship"] = df["PI_HAT"].apply(
    lambda x: "related" if x >= threshold else "unrelated"
)

# Now only keep the relevant columns
truth = df[["IID1", "IID2", "PI_HAT", "relationship"]]

truth.to_csv("../results/ground_truth.csv", sep="\t", index=False)

print("Ground truth (binary) generated.\n")

print(truth["relationship"].value_counts())

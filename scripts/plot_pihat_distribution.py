import pandas as pd
import matplotlib.pyplot as plt

# load the PLINK IBD results
df = pd.read_csv("results/plink_ibd.genome", delim_whitespace=True)

threshold = 0.125

# Plot the histogram
plt.figure()
plt.hist(df["PI_HAT"], bins=50)
plt.axvline(threshold)
plt.xlabel("PI_HAT")
plt.ylabel("Number of Pairs")
plt.title("Distribution of Pairwise PI_HAT Values")

plt.savefig("results/pihat_distribution.png", dpi=300)

print("Plot saved as png")

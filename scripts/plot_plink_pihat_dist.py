import pandas as pd
import matplotlib.pyplot as plt

# load the PLINK IBD results
df = pd.read_csv("results/plink_ibd.genome", delim_whitespace=True)

# Plot the histograms
fig, axes = plt.subplots(2, 1, figsize=(14, 12))
axes[0].hist(df["PI_HAT"], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('PI_HAT', fontsize=12)
axes[0].set_ylabel('Number of Pairs', fontsize=12)
axes[0].set_title("Plink PI_HAT Values", fontsize=14)
axes[0].grid(True, alpha=0.3)
axes[0].axvline(df['PI_HAT'].median(), color='red', linestyle='--',
                label=f"Median: {df['PI_HAT'].median():.4f}")
axes[0].axvline(df['PI_HAT'].median(), color='green', linestyle='--',
                label=f"Median: {df['PI_HAT'].median():.4f}")
axes[0].legend()


axes[1].hist(df["PI_HAT"], bins=50, color='coral', edgecolor='black', alpha=0.7)
axes[1].set_xlabel('PI_HAT', fontsize=12)
axes[1].set_ylabel('Number of Pairs', fontsize=12)
axes[1].set_title("Plink PI_HAT Values (Log scale)", fontsize=14)
axes[1].set_yscale('log')
axes[1].grid(True, alpha=0.3)

plt.savefig("results/plink_pihat_distribution.png", dpi=150)

print("Plot saved as png")

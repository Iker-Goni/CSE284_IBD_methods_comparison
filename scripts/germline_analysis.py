import pandas as pd
import matplotlib.pyplot as plt
import gzip
import glob
import os
import numpy as np

# Configure paths and gather files
GERMLINE_FILE = "data/lwk_ibd_allchr.match"
OUTPUT_DIR = "results/germline"
os.makedirs(OUTPUT_DIR, exist_ok=True)

colnames = [
    "fid1","iid1","fid2","iid2","chr","start","end",
    "snp1","snp2","nsnp","error","unit",
    "err1","err2","err3"
           ]

df = pd.read_csv(GERMLINE_FILE, sep=r"\s+", names=colnames)
df['length_mb'] = (df["end"] - df["start"]) / 1000000
print(f"\nTotal segments read: {len(df):,}")
print(f"Mean segment length: {df['length_mb'].mean():.3f} Mb")
print(f"Median segment length: {df['length_mb'].median():.3f} Mb")
print(f"Max segment length: {df['length_mb'].max():.3f} Mb")


# Ensures that pairs are not double counted
def merge_intervals(df):
    df_sorted = df.sort_values('start')
    
    # Iteratively merge overlapping or adjacent intervals
    merged = []
    curr_start, curr_end = df_sorted.iloc[0]['start'], df_sorted.iloc[0]['end']
    
    for _, row in df_sorted.iloc[1:].iterrows():
        if row['start'] <= curr_end:
            # Overlap found: extend the current end if this segment goes further
            curr_end = max(curr_end, row['end'])
        else:
            # No overlap: save finished segment and start new one
            merged.append(curr_end - curr_start)
            curr_start, curr_end = row['start'], row['end']
            
    merged.append(curr_end - curr_start)
    
    # Return total physical length in Mb
    return sum(merged) / 1000000

print("Calculating IBD sharing")

# Group by pair and chromosome to merge within genomic regions
df = df[df['iid1'] != df['iid2']]
df['pair_key'] = list(map(tuple, np.sort(df[['iid1','iid2']].values, axis=1)))
pair_chrom_groups = df.groupby(['pair_key', 'chr'])
merged_data = pair_chrom_groups.apply(merge_intervals).reset_index(name='merged_mb')

# Sum the merged chromosome totals for each pair
pair_ibd = merged_data.groupby('pair_key')['merged_mb'].sum().reset_index()
pair_ibd.columns = ['pair', 'total_ibd_mb']

# Splitting keys back into sample columns
pair_ibd['sample_a'] = pair_ibd['pair'].apply(lambda x: x[0])
pair_ibd['sample_b'] = pair_ibd['pair'].apply(lambda x: x[1])
pair_ibd = pair_ibd.sort_values('total_ibd_mb', ascending=False).reset_index(drop=True)

print(f"Total unique pairs with IBD sharing: {len(pair_ibd):,}")
print(f"\nTop 10 pairs by IBD sharing:")
print(pair_ibd.head(10))

# Estimate pi_hat value from IBD results
genome_length = 3200 # in Mb
pair_ibd['pi_hat'] = pair_ibd['total_ibd_mb'] / genome_length
print(f"\nPi_hat statistics:")
print(f"  Mean pi_hat: {pair_ibd['pi_hat'].mean():.6f}")
print(f"  Median pi_hat: {pair_ibd['pi_hat'].median():.6f}")
print(f"  Max pi_hat: {pair_ibd['pi_hat'].max():.6f}")


print("\nCreating distribution plots...")
fig, axes = plt.subplots(2, 1, figsize=(14, 12))

# Histogram with cumulative IBD per pair
axes[0].hist(pair_ibd['total_ibd_mb'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Cumulative IBD Sharing (Mb)', fontsize=12)
axes[0].set_ylabel('Number of Pairs', fontsize=12)
axes[0].set_title('Germline Total IBD Sharing per Pair', fontsize=14)
axes[0].grid(True, alpha=0.3)
axes[0].axvline(pair_ibd['total_ibd_mb'].median(), color='red', linestyle='--', 
                   label=f"Median: {pair_ibd['total_ibd_mb'].median():.2f} Mb")
axes[0].axvline(pair_ibd['total_ibd_mb'].mean(), color='green', linestyle='--', 
                   label=f"Mean: {pair_ibd['total_ibd_mb'].mean():.2f} Mb")
axes[0].legend()

# Log scale histogram
axes[1].hist(pair_ibd['total_ibd_mb'], bins=50, color='coral', edgecolor='black', alpha=0.7)
axes[1].set_xlabel('Cumulative IBD Sharing (Mb)', fontsize=12)
axes[1].set_ylabel('Number of Pairs (log scale)', fontsize=12)
axes[1].set_title('Germline Total IBD Sharing per Pair (Log Scale)', fontsize=14)
axes[1].set_yscale('log')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'cumulative_ibd_distribution.png'), dpi=150)
print(f"Distribution plots saved to {OUTPUT_DIR}/")


print("\nCalculating per-chromosome IBD contributions...")

chrom_ibd = df.groupby('chr')['length_mb'].sum().reset_index()
chrom_ibd.columns = ['chromosome', 'total_ibd_mb']
chrom_ibd['percentage'] = (chrom_ibd['total_ibd_mb'] / chrom_ibd['total_ibd_mb'].sum() * 100).round(2)

print("\nTotal IBD by chromosome:")
print(chrom_ibd)

# Create chromosome contribution plot
fig, ax = plt.subplots(figsize=(12, 6))
bars = ax.bar(chrom_ibd['chromosome'].astype(str), chrom_ibd['total_ibd_mb'], 
              color='darkorange', alpha=0.7, edgecolor='black')
ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('Total IBD (Mb)', fontsize=12)
ax.set_title('Germline Total IBD Contribution by Chromosome', fontsize=14)
ax.grid(True, alpha=0.3, axis='y')

# Add percentage labels
for bar, pct in zip(bars, chrom_ibd['percentage']):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 5,
            f'{pct}%', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'ibd_by_chromosome.png'), dpi=150)

# Create pi_hat distribution plots (normal and log)
fig, axes = plt.subplots(2, 1, figsize=(14, 12))
axes[0].hist(pair_ibd['pi_hat'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Pi_hat', fontsize=12)
axes[0].set_ylabel('Number of Pairs', fontsize=12)
axes[0].set_title('Germline Pi_hat Values', fontsize=14)
axes[0].grid(True, alpha=0.3)
axes[0].axvline(pair_ibd['pi_hat'].median(), color='red', linestyle='--', 
                label=f"Median: {pair_ibd['pi_hat'].median():.4f}")
axes[0].axvline(pair_ibd['pi_hat'].mean(), color='green', linestyle='--', 
                label=f"Mean: {pair_ibd['pi_hat'].mean():.4f}")
axes[0].legend()

axes[1].hist(pair_ibd['pi_hat'], bins=50, color='coral', edgecolor='black', alpha=0.7)
axes[1].set_xlabel('Pi_hat', fontsize=12)
axes[1].set_ylabel('Number of Pairs', fontsize=12)
axes[1].set_title('Germline Pi_hat Distribution (Log Scale)', fontsize=14)
axes[1].set_yscale('log')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'pi_hat_distribution.png'), dpi=150)
print(f"Pi_hat distribution plots saved to {OUTPUT_DIR}/pi_hat_distribution.png")


print("\nExporting results...")

# Save full pair IBD data
pair_ibd.to_csv(os.path.join(OUTPUT_DIR, 'pairwise_ibd_sharing.csv'), index=False)
print(f"  Pairwise IBD saved to {OUTPUT_DIR}/pairwise_ibd_sharing.csv")

# Save per-chromosome summary
chrom_ibd.to_csv(os.path.join(OUTPUT_DIR, 'chromosome_ibd_summary.csv'), index=False)

# Create summary statistics file
summary_file = os.path.join(OUTPUT_DIR, 'ibd_summary_stats.txt')
with open(summary_file, 'w') as f:
    f.write("Germline IBD Analysis Summary\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
    
    f.write("SEGMENT STATISTICS\n")
    f.write("-" * 30 + "\n")
    f.write(f"Total segments detected: {len(df):,}\n")
    f.write(f"Mean segment length: {df['length_mb'].mean():.3f} Mb\n")
    f.write(f"Median segment length: {df['length_mb'].median():.3f} Mb\n")
    f.write(f"Total IBD across all segments: {df['length_mb'].sum():.2f} Mb\n\n")
    
    f.write("PAIRWISE STATISTICS\n")
    f.write("-" * 30 + "\n")
    f.write(f"Total pairs with IBD: {len(pair_ibd):,}\n")
    f.write(f"Mean IBD per pair: {pair_ibd['total_ibd_mb'].mean():.3f} Mb\n")
    f.write(f"Median IBD per pair: {pair_ibd['total_ibd_mb'].median():.3f} Mb\n")
    f.write(f"Max IBD per pair: {pair_ibd['total_ibd_mb'].max():.3f} Mb\n\n")
    
    f.write("IBD PERCENTILES (per pair)\n")
    f.write("-" * 30 + "\n")
    for p in [10, 25, 50, 75, 90, 95, 99]:
        percentile_val = pair_ibd['total_ibd_mb'].quantile(p/100)
        f.write(f"  {p}th percentile: {percentile_val:.3f} Mb\n")
        
    f.write("\n\nPI_HAT STATISTICS\n")
    f.write("-" * 30 + "\n")
    f.write(f"Genome length used for pi_hat: {genome_length:.0f} Mb\n")
    f.write(f"Mean pi_hat: {pair_ibd['pi_hat'].mean():.6f}\n")
    f.write(f"Median pi_hat: {pair_ibd['pi_hat'].median():.6f}\n")
    f.write(f"Max pi_hat: {pair_ibd['pi_hat'].max():.6f}\n")
    f.write(f"Min pi_hat (>0): {pair_ibd[pair_ibd['pi_hat'] > 0]['pi_hat'].min():.6f}\n\n")
    
    f.write("PI_HAT PERCENTILES\n")
    f.write("-" * 30 + "\n")
    for p in [10, 25, 50, 75, 90, 95, 99, 99.5, 99.9]:
        percentile_val = pair_ibd['pi_hat'].quantile(p/100)
        f.write(f"  {p}th percentile: {percentile_val:.6f}\n")
        
    pair_ibd[['sample_a', 'sample_b', 'total_ibd_mb', 'pi_hat']].to_csv(
    os.path.join(OUTPUT_DIR, 'pairwise_pi_hat.csv'), index=False)

print(f"  Summary statistics saved to {summary_file}")









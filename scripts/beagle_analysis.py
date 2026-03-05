import pandas as pd
import matplotlib.pyplot as plt
import gzip
import glob
import os

# Configure paths and gather files
BEAGLE_DIR = "data/beagle"
OUTPUT_DIR = "results/beagle"
MIN_LOD = 3.0 # Standard value
os.makedirs(OUTPUT_DIR, exist_ok=True)

colnames = ['sample1', 'hap1', 'sample2', 'hap2', 'chrom', 'start', 'end', 'lod']
all_segments = []

# Gather .ibd files and combine into a single df
file_pattern = os.path.join(BEAGLE_DIR, "ps2_ibd_chr*.ibd.gz")
file_list = sorted(glob.glob(file_pattern))

if not file_list:
    print(f"Error: No files found in {BEAGLE_DIR}")
    exit(1)

print(f"Found IBD files for {len(file_list)} chromosomes")

for file_path in file_list:
    chromosome = file_path.split('chr')[-1].split('.')[0]
    with gzip.open(file_path, 'rt') as f:
        df_chr = pd.read_csv(f, sep='\s+', names=colnames, comment='#')
        all_segments.append(df_chr)
        
df_all = pd.concat(all_segments, ignore_index=True)
print(f"\nTotal segments read: {len(df_all):,}")

# Filter for segments with a LOD above MIN_LOD
print(f"\nFiltering segments (LOD >= {MIN_LOD})...")
df_filtered = df_all[df_all['lod'] >= MIN_LOD].copy()
df_filtered['length_mb'] = (df_filtered['end'] - df_filtered['start']) / 1000000

print(f"Mean segment length: {df_filtered['length_mb'].mean():.3f} Mb")
print(f"Median segment length: {df_filtered['length_mb'].median():.3f} Mb")
print(f"Max segment length: {df_filtered['length_mb'].max():.3f} Mb")


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
df_filtered['pair_key'] = df_filtered.apply(lambda row: tuple(sorted([row['sample1'], row['sample2']])), axis=1)
pair_chrom_groups = df_filtered.groupby(['pair_key', 'chrom'])
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



print("\nCreating distribution plots...")
fig, axes = plt.subplots(2, 1, figsize=(14, 12))

# Histogram with cumulative IBD per pair
axes[0].hist(pair_ibd['total_ibd_mb'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Cumulative IBD Sharing (Mb)', fontsize=12)
axes[0].set_ylabel('Number of Pairs', fontsize=12)
axes[0].set_title('Distribution of Total IBD Sharing per Pair', fontsize=14)
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
axes[1].set_title('Distribution of Total IBD Sharing per Pair (Log Scale)', fontsize=14)
axes[1].set_yscale('log')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'cumulative_ibd_distribution.png'), dpi=150)
print(f"Distribution plots saved to {OUTPUT_DIR}/")


print("\nCalculating per-chromosome IBD contributions...")

chrom_ibd = df_filtered.groupby('chrom')['length_mb'].sum().reset_index()
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
ax.set_title('Total IBD Contribution by Chromosome', fontsize=14)
ax.grid(True, alpha=0.3, axis='y')

# Add percentage labels
for bar, pct in zip(bars, chrom_ibd['percentage']):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 5,
            f'{pct}%', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'ibd_by_chromosome.png'), dpi=150)


print("\nExporting results...")

# Save full pair IBD data
pair_ibd.to_csv(os.path.join(OUTPUT_DIR, 'pairwise_ibd_sharing.csv'), index=False)
print(f"  Pairwise IBD saved to {OUTPUT_DIR}/pairwise_ibd_sharing.csv")

# Save per-chromosome summary
chrom_ibd.to_csv(os.path.join(OUTPUT_DIR, 'chromosome_ibd_summary.csv'), index=False)

# Create summary statistics file
summary_file = os.path.join(OUTPUT_DIR, 'ibd_summary_stats.txt')
with open(summary_file, 'w') as f:
    f.write("BEAGLE 4.1 IBD Analysis Summary\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"Filter settings: LOD ≥ {MIN_LOD}\n\n")
    
    f.write("SEGMENT STATISTICS\n")
    f.write("-" * 30 + "\n")
    f.write(f"Total segments detected: {len(df_filtered):,}\n")
    f.write(f"Mean segment length: {df_filtered['length_mb'].mean():.3f} Mb\n")
    f.write(f"Median segment length: {df_filtered['length_mb'].median():.3f} Mb\n")
    f.write(f"Total IBD across all segments: {df_filtered['length_mb'].sum():.2f} Mb\n\n")
    
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

print(f"  Summary statistics saved to {summary_file}")









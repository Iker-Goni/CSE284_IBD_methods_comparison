#!/bin/bash
# Script 1: Preprocessing and phasing
set -e
set -o pipefail

# Define prefixes
RAW_PREFIX="data/ps2_ibd.lwk"
QC_PREFIX="data/lwk_qc"
PRUNE_PREFIX="data/lwk_prune"
LD_PREFIX="data/lwk_ld"
PHASE_PREFIX="data/lwk_phased"

echo "=== Step 1: QC filtering ==="
plink \
  --bfile ${RAW_PREFIX} \
  --geno 0.05 \
  --maf 0.05 \
  --make-bed \
  --out ${QC_PREFIX}

echo "=== Step 2: LD pruning ==="
plink \
  --bfile ${QC_PREFIX} \
  --indep-pairwise 50 5 0.2 \
  --out ${PRUNE_PREFIX}

echo "=== Step 3: Extract pruned SNPs ==="
plink \
  --bfile ${QC_PREFIX} \
  --extract ${PRUNE_PREFIX}.prune.in \
  --make-bed \
  --out ${LD_PREFIX}

echo "=== Step 4: Convert to VCF for phasing ==="
plink \
  --bfile ${LD_PREFIX} \
  --recode vcf \
  --out ${LD_PREFIX}

echo "=== Step 5: Phase with Beagle ==="
java -jar tools/beagle/beagle.27Feb25.75f.jar \
  gt=${LD_PREFIX}.vcf \
  out=${PHASE_PREFIX}

echo "=== Step 6: Filter multi-allelic SNPs ==="
bcftools view -m2 -M2 -v snps ${PHASE_PREFIX}.vcf.gz -Oz -o ${PHASE_PREFIX}_clean.vcf.gz
tabix -p vcf ${PHASE_PREFIX}_clean.vcf.gz

echo "=== Step 7: Convert phased VCF back to PED/MAP ==="
plink \
  --vcf ${PHASE_PREFIX}_clean.vcf.gz \
  --recode ped \
  --out ${PHASE_PREFIX}

echo "Preprocessing and phasing completed. Output: ${PHASE_PREFIX}.ped and ${PHASE_PREFIX}.map"

# --- Subset phased PED/MAP per chromosome ---
for chr in {1..22}; do
    echo "Subsetting chromosome $chr..."
    plink \
      --file ${PHASE_PREFIX} \
      --chr $chr \
      --recode ped \
      --out ${PHASE_PREFIX}_chr${chr}
done

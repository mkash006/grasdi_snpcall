#!/bin/bash
#SBATCH --job-name=GenomicsDB_Import
#SBATCH --partition=epyc
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --array=0-23
#SBATCH --output=logs/gatk_genomicsdb_%A_%a.out
#SBATCH --error=logs/gatk_genomicsdb_%A_%a.err

# Load necessary modules
module load gatk/4.3.0.0

# Configuration
GVCF_DIR="results/vcf"
SAMPLE_MAP="sample_map.txt"
DB_DIR="results/genomicsdb_per_chr"
SCRATCH_DIR="${SCRATCH:-/scratch/$USER/genomicsdb_tmp}"
THREADS=8

CHROM_LIST=(
  "LG1_scaffold_3" "LG2_1_scaffold_12" "LG2_2_scaffold_24"
  "LG3_scaffold_2" "LG4_scaffold_10" "LG5_scaffold_8"
  "LG6_scaffold_11" "LG7_scaffold_7" "LG8_scaffold_19"
  "LG9_scaffold_4" "LG10_scaffold_9" "LG11_scaffold_16"
  "LG12_scaffold_17" "LG13_scaffold_5" "LG14_scaffold_15"
  "LG15_scaffold_14" "LG16_scaffold_5" "LG17_scaffold_13"
  "LG18_scaffold_20" "LG19_scaffold_18" "LG20_scaffold_21"
  "LG21_scaffold_22" "LG22_scaffold_1" "LG23_scaffold_23"
)

CHROM=${CHROM_LIST[$SLURM_ARRAY_TASK_ID]}

set -euo pipefail

mkdir -p "$DB_DIR"
TMP_PATH="${SCRATCH_DIR}/${CHROM}_tmp"
mkdir -p "$TMP_PATH"

# Generate sample map if needed
if [ ! -s "$SAMPLE_MAP" ]; then
  echo "Creating sample map: $SAMPLE_MAP"
  > "$SAMPLE_MAP"
  for gvcf in "$GVCF_DIR"/*.g.vcf.gz; do
    [ -e "$gvcf" ] || continue
    idx="${gvcf}.tbi"
    if [[ ! -f "$idx" ]]; then
      echo "Index missing for: $gvcf — indexing now"
      gatk IndexFeatureFile -I "$gvcf"
    fi
    sample_name=$(basename "$gvcf" .g.vcf.gz)
    echo -e "${sample_name}\t${gvcf}" >> "$SAMPLE_MAP"
  done
fi

# Clean Windows-style line endings
sed -i 's/\r$//' "$SAMPLE_MAP"

echo "==== Starting $CHROM with $THREADS threads ===="

gatk --java-options "-Xmx100G" GenomicsDBImport \
  --genomicsdb-workspace-path "${DB_DIR}/${CHROM}" \
  --sample-name-map "$SAMPLE_MAP" \
  --reader-threads "$THREADS" \
  --intervals "$CHROM" \
  --tmp-dir "$TMP_PATH"

echo "==== Finished $CHROM ===="

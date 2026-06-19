#!/bin/bash -l
#SBATCH --job-name=snpcall_init
#SBATCH --output=logs/snakemake_init_%j.out
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load snakemake/7.18
module load bwa
module load samtools
module load gatk
export PYTHONNOUSERSITE=1
mkdir -p logs

snakemake \
    --snakefile Snakefile \
    --jobs 100 \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb} --time=24:00:00 --output=logs/slurm-%j.out"

echo "Snakemake workflow submitted."

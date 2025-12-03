# GRAS-Di SNP Calling Pipeline

This repository contains a GATK-based variant calling pipeline designed for GRAS-Di (Genotyping by Random Amplicon Sequencing–Direct) sequencing data using a combination of `Snakemake` and `bash`. The pipeline is tailored for use on an HPC cluster with SLURM and follows **GATK Best Practices** for variant calling.

**Note:** All steps assume input files have already been **trimmed to remove adapter contamination**. Raw FASTQ files must be pre-processed using tools such as **Trimmomatic** or **fastp** before running this pipeline.

---

##  Pipeline Overview

###  Snakefile: Per-sample GVCF Calling
The [`Snakefile`](Snakefile) performs **per-sample variant calling** and generates **GVCF files** via GATK HaplotypeCaller. Key steps include:

1. Reference Indexing
2. Alignment to Reference (BWA-MEM)
3. BAM Sorting & Marking Duplicates
4. Read Group Assignment
5. GATK HaplotypeCaller in GVCF mode

All outputs are organized in `results/`, and logs are saved in the `logs/` directory.

---

### `.sh` Scripts for initiating snakemake and joint genotyping using GATK gvcf mode

- [Snakemake_init.sh](Snakemake_init.sh)
- [make_genomicsdb.sh](make_genomicsdb.sh)
- [genotype_gvcf.sh](genotype_gvcf.sh)
- [merge_vcf_allchr.sh](merge_vcf_allchr.sh)

Making Genomicsdb and GenotypeGVCFs are separated intentionally from Snakemake because:
- Selecting scaffolds or chromosomes of interest can be manually performed using [[make_genomicsdb.sh]], [[genotype_gvcf.sh]] and [[merge_vcf_allchr.sh]]
- This workflow allows to accumulation of g.vcf files from separate sequencing runs and then combines the  haplotypecaller step for separate runs
---

##  Requirements

- Snakemake
- BWA
- SAMtools
- GATK 
- Python 3
- SLURM (for HPC execution)

Optional tools:
- Trimmomatic or fastp (for trimming before this pipeline)
- bcftools, plink (for downstream filtering)

---

##  Usage Instructions

1. Make sure all FASTQ files in `data/` are **adapter-trimmed**.
2. The reference genome should be placed in `reference/` before initiating Snakemake
3. Run the Snakemake pipeline by submitting Snakemake_init.sh to an HPC using Slurm job scheduler or via the following terminal command if running locally. 

```bash
snakemake --jobs 100 --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb} --time=24:00:00" --use-conda

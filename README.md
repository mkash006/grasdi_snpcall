# GRAS-Di SNP Calling Pipeline

GATK Best Practices variant calling for GRAS-Di sequencing data, built with
Snakemake + bash for SLURM on UCR HPCC. Input reads are assumed already
adapter-trimmed (Trimmomatic/fastp). **The reference genome is assumed
pre-indexed** — this workflow does not include reference indexing.

## Layout

```
.
├── Snakefile                 # per-sample GVCF generation
├── config/config.yaml        # paths, read naming, resources
├── Snakemake_init.sh         # SLURM submission wrapper
├── make_genomicsdb.sh        # joint: GenomicsDBImport
├── genotype_gvcf.sh          # joint: GenotypeGVCFs
├── merge_vcf_allchr.sh       # joint: concat per-chr VCFs
├── data/                     # trimmed paired FASTQ
├── reference/                # pre-indexed genome (.fai .dict bwa indices)
├── results/                  # outputs (bam, metrics, gvcf, vcf)
└── logs/
```

## Per-sample workflow (Snakefile)

1. **bwa_map** — BWA-MEM alignment with inline read groups, coordinate-sorted
2. **adding read groups** — GATK add read groups
3. **haplotype_caller** — HaplotypeCaller in `-ERC GVCF` mode

Output: `results/gvcf/<sample>.g.vcf.gz`

## Reference prerequisites

Place the indexed genome in `reference/`:
```
P_prolifica.ref.fasta       P_prolifica.ref.fasta.amb
P_prolifica.ref.fasta.fai   P_prolifica.ref.fasta.ann
P_prolifica.ref.dict        P_prolifica.ref.fasta.bwt
                            P_prolifica.ref.fasta.pac
                            P_prolifica.ref.fasta.sa
```

## Configuration

Edit `config/config.yaml` — sample names are auto-detected from
`data/*<r1_suffix>`. Set `r1_suffix` / `r2_suffix` to match your trimmed
FASTQ naming before running.

## Running

Dry run:
```
snakemake -n -p
```

Submit on slurm HPCC:
```
sbatch Snakemake_init.sh
```

## Joint genotyping (run after all per-sample GVCFs exist)

Separated from Snakemake so GVCFs from multiple sequencing runs can be
accumulated, and so scaffolds/chromosomes of interest can be selected
manually via the interval list.

```
# 1. list chromosomes/scaffolds of interest in reference/intervals.list
sbatch make_genomicsdb.sh
sbatch genotype_gvcf.sh
sbatch merge_vcf_allchr.sh   # only if genotyping was split by interval
```

## Requirements

`module load`: snakemake/7.18, bwa, samtools, gatk, bcftools

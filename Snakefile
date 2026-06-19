import os
from glob import glob

# -----------------------------
# Samples and reference
# -----------------------------
SAMPLES = sorted(set(
    os.path.basename(f).replace('_R1_paired.fq.gz', '')
    for f in glob("data/*_R1_paired.fq.gz")
))

REF = "reference/P_prolifica.ref.fasta"

# -----------------------------
# Resources
# threads / mem_mb sized for a 720 Mbp reference and light GRAS-Di samples.
# Kept lean so SLURM can pack many of the 1400 jobs in parallel.
# -----------------------------
ALIGN_THREADS, ALIGN_MEM = 12, 24000
RG_THREADS,    RG_MEM    = 2,  8000
INDEX_THREADS            = 2
HC_THREADS,    HC_MEM    = 8,  32000

SORT_MEM_PER_THREAD = "1500M"   # caps samtools sort so it cannot balloon

# -----------------------------
# Final targets
# -----------------------------
rule all:
    input:
        expand("results/vcf/{sample}_snp_calls.g.vcf.gz", sample=SAMPLES),
        expand("results/bam/{sample}_rg.bam", sample=SAMPLES),
        expand("results/bam/{sample}_rg.bam.bai", sample=SAMPLES)

# -----------------------------
# Align + sort
# Reference is already indexed (.fai, .dict, bwa indices in reference/)
# -----------------------------
rule align_and_sort:
    input:
        r1 = "data/{sample}_R1_paired.fq.gz",
        r2 = "data/{sample}_R2_paired.fq.gz",
        ref = REF
    output:
        bam = "results/bam/{sample}_aligned.bam"
    threads: ALIGN_THREADS
    resources: mem_mb=ALIGN_MEM
    log:
        "logs/align/{sample}.log"
    shell:
        """
        set -euo pipefail
        module load bwa/0.7.17
        module load samtools/1.18

        echo "[`date`] Aligning {wildcards.sample}" > {log}
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2>> {log} | \
        samtools sort -@ {threads} -m %s -o {output.bam} - 2>> {log}
        echo "[`date`] Finished alignment {wildcards.sample}" >> {log}
        """ % SORT_MEM_PER_THREAD

# -----------------------------
# Add read group
# -----------------------------
rule add_readgroup:
    input:
        bam = "results/bam/{sample}_aligned.bam"
    output:
        bam = "results/bam/{sample}_rg.bam"
    threads: RG_THREADS
    resources: mem_mb=RG_MEM
    log:
        "logs/readgroup/{sample}.log"
    shell:
        """
        set -euo pipefail
        module load gatk/4.3.0.0

        ID="RG1"
        LB="Lib1"
        PL="ILLUMINA"
        PU="Unit1"
        SM="{wildcards.sample}"

        gatk AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -ID $ID -LB $LB -PL $PL -PU $PU -SM $SM \
            >> {log} 2>&1
        """

# -----------------------------
# Index read-group BAM
# -----------------------------
rule index_bam:
    input:
        "results/bam/{sample}_rg.bam"
    output:
        "results/bam/{sample}_rg.bam.bai"
    threads: INDEX_THREADS
    log:
        "logs/index_bam/{sample}.log"
    shell:
        """
        set -euo pipefail
        module load samtools/1.18

        samtools index {input} >> {log} 2>&1
        """

# -----------------------------
# HaplotypeCaller
# Heap is derived from mem_mb so the SLURM request and -Xmx never drift apart.
# -----------------------------
rule haplotype_caller:
    input:
        bam = "results/bam/{sample}_rg.bam",
        bai = "results/bam/{sample}_rg.bam.bai",
        ref = REF
    output:
        gvcf = "results/vcf/{sample}_snp_calls.g.vcf.gz",
        bam = "results/bam_snpcall/{sample}_hpcall.bam"
    threads: HC_THREADS
    resources: mem_mb=HC_MEM
    params:
        heap = lambda wc, resources: int(resources.mem_mb * 0.85)
    log:
        "logs/haplotypecaller/{sample}.log"
    shell:
        """
        set -euo pipefail
        module load gatk/4.3.0.0

        gatk --java-options "-Xmx{params.heap}m -XX:ParallelGCThreads=2" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            --sample-name {wildcards.sample} \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads {threads} \
            --bam-output {output.bam} \
            -O {output.gvcf} >> {log} 2>&1
        """

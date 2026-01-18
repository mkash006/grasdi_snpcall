import os
from glob import glob

# -----------------------------
# Samples and reference
# -----------------------------
SAMPLES = sorted(set(
    os.path.basename(f).replace('_R1_paired.fq.gz', '')
    for f in glob("data/*_R1_paired.fq.gz")
))

REF = "reference/p_reticulata.ref.fasta"

ALIGN_THREADS = 20
RG_THREADS = 4
HC_THREADS = 20

# -----------------------------
# Final targets
# -----------------------------
rule all:
    input:
        expand("results/vcf/{sample}_snp_calls.g.vcf.gz", sample=SAMPLES),
        expand("results/bam/{sample}_rg.bam", sample=SAMPLES),
        expand("results/bam/{sample}_rg.bam.bai", sample=SAMPLES)

# -----------------------------
# Index reference
# -----------------------------
rule index_reference:
    input:
        REF
    output:
        REF + ".bwt",
        REF + ".fai",
        REF.replace(".fasta", ".dict")
    threads: 4
    resources: mem_mb=32000
    log:
        "logs/index_reference.log"
    shell:
        """
        set -euo pipefail
        module load bwa/0.7.17
        module load samtools/1.18
        module load gatk/4.3.0.0

        echo "[`date`] Indexing reference" > {log}
        bwa index {input} >> {log} 2>&1
        samtools faidx {input} >> {log} 2>&1
        gatk CreateSequenceDictionary -R {input} >> {log} 2>&1
        """

# -----------------------------
# Align + sort
# -----------------------------
rule align_and_sort:
    input:
        index = rules.index_reference.output,
        r1 = "data/{sample}_R1_paired.fq.gz",
        r2 = "data/{sample}_R2_paired.fq.gz",
        ref = REF
    output:
        bam = "results/bam/{sample}_aligned.bam"
    threads: ALIGN_THREADS
    resources: mem_mb=160000
    log:
        "logs/align/{sample}.log"
    shell:
        """
        set -euo pipefail
        module load bwa/0.7.17
        module load samtools/1.18

        echo "[`date`] Aligning {wildcards.sample}" > {log}
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2>> {log} | \
        samtools view -Sb - 2>> {log} | \
        samtools sort -@ {threads} -o {output.bam} 2>> {log}
        echo "[`date`] Finished alignment {wildcards.sample}" >> {log}
        """

# -----------------------------
# Add read group
# -----------------------------
rule add_readgroup:
    input:
        bam = "results/bam/{sample}_aligned.bam"
    output:
        bam = "results/bam/{sample}_rg.bam"
    threads: RG_THREADS
    resources: mem_mb=160000
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
# NEW: Index read-group BAM
# -----------------------------
rule index_bam:
    input:
        "results/bam/{sample}_rg.bam"
    output:
        "results/bam/{sample}_rg.bam.bai"
    threads: 2
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
    resources: mem_mb=220000
    log:
        "logs/haplotypecaller/{sample}.log"
    shell:
        """
        set -euo pipefail
        module load gatk/4.3.0.0

        gatk --java-options "-Xmx120G -XX:ParallelGCThreads={threads}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            --sample-name {wildcards.sample} \
            --emit-ref-confidence GVCF \
            --bam-output {output.bam} \
            -O {output.gvcf} >> {log} 2>&1
        """


# snakemake -s piccard_markduplicates.smk --configfile config.yml --profile slurm --nodes 50 --use-conda
ssample = config["samples"].keys()

rule all:
    input: 
        expand("/scratch/tbrittoborges/mage_dedup/{sample}.bam", sample=sample),
        expand("/scratch/tbrittoborges/mage_dedup/{sample}.bam.bai", sample=sample)

rule mark_duplicates:
    input:
        "mappings/{sample}.bam"
    output:
        bam="/scratch/tbrittoborges/mage_dedup/{sample}.bam",
        metrics="/scratch/tbrittoborges/mage_dedup/{sample}.metrics.txt"
    log:
        "/scratch/tbrittoborges/mage_dedup/{sample}.log"
    resources:
        mem_mb = 32000
    params:
        "REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate TMP_DIR=/scratch/tbrittoborges/tmp/"
    wrapper:
        "0.66.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "/scratch/tbrittoborges/mage_dedup/{sample}.bam"
    output:
        "/scratch/tbrittoborges/mage_dedup/{sample}.bam.bai"
    log:
        "/scratch/tbrittoborges/mage_dedup/index_{sample}.log"
    params:
        "-@ 10" 
    wrapper:
        "0.66.0/bio/samtools/index"
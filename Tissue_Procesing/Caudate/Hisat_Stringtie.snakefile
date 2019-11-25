configfile: "/scratch/mjpete11/GTEx/Configs/Sex_Sorted_Quantification.config.json"

import os

# TOOLS
# Variables should point to full paths to tools.
# Currently assumes all are available in one's PATH
# (e.g., installed with Bioconda)
HISAT2_PATH = "hisat2"
SAMTOOLS_PATH = "samtools"
STRINGTIE_PATH = "stringtie"

# Directory variables
FASTQ_DIRECTORY = "/mnt/storage/public/dbgap-8834/brain_caudate/"

# Samples
XX_SAMPLES = config["Caudate_RNA_Trimmed"]["Female"]
XY_SAMPLES = config["Caudate_RNA_Trimmed"]["Male"]
SAMPLES = XX_SAMPLES + XY_SAMPLES

rule all:
    input:
        expand("stringtie_results/{sample}/{sample}_GRCh38.fully_covered_transcripts.secondpass.gtf", sample=SAMPLES)
        
rule hisat2_align_reads:
    input:
        paired_1 = "/scratch/mjpete11/GTEx/Caudate/Quality_Control/trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_1.fastq.gz",
        paired_2 = "/scratch/mjpete11/GTEx/Caudate/Quality_Control/trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_2.fastq.gz",
        xx = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.8.ht2",
        xy = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/GRCh38.p12.genome.XY/GRCh38.p12.genome.XY.8.ht2"
    output:
        "bams/{sample}_GRCh38.sorted.bam"
    params:
        hisat2 = HISAT2_PATH,
        samtools = SAMTOOLS_PATH,
        threads = 4,
        base_name_xx = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly",
        base_name_xy = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/GRCh38.p12.genome.XY/GRCh38.p12.genome.XY"
    run:
        if wildcards.sample in XY_SAMPLES:
            shell(
                "{params.hisat2} --dta -p {params.threads} -x {params.base_name_xy} "
                "-1 {input.paired_1} -2 {input.paired_2} | "
                "{params.samtools} view -b - | "
                "{params.samtools} sort -O bam -o {output} -")
        else:
            shell(
                "{params.hisat2} --dta -p {params.threads} -x {params.base_name_xx} "
                "-1 {input.paired_1} -2 {input.paired_2} | "
                "{params.samtools} view -b - | "
                "{params.samtools} sort -O bam -o {output} -")

rule stringtie_first_pass:
    input:
        bam = "bams/{sample}_GRCh38.sorted.bam",
        gff = "/scratch/mjpete11/GTEx/gencode.v29.annotation.gtf"
    output:
        "stringtie_results/{sample}/{sample}_GRCh38.assembled_transcripts.firstpass.gtf"
    threads: 4
    params:
        stringtie = STRINGTIE_PATH,
        threads = 4
    shell:
        "{params.stringtie} {input.bam} -o {output} -p {params.threads} "
        "-G {input.gff}"

rule create_stringtie_merged_list:
    input: 
         expand("stringtie_results/{sample}/{sample}_GRCh38.assembled_transcripts.firstpass.gtf", sample = SAMPLES)
    output: 
        "stringtie_results/GRCh38_gtflist.txt"
    run:
        shell("echo -n > {output}")
        for i in input:
            shell("echo {} >> {{output}}".format(i))

rule stringtie_merge:
    input: 
        stringtie_list = "stringtie_results/GRCh38_gtflist.txt",
        gff = "/scratch/mjpete11/GTEx/gencode.v29.annotation.gtf"
    output:
        "stringtie_results/GRCh38.merged.gtf"
    threads: 4
    params:
        stringtie = STRINGTIE_PATH, 
        threads = 4
    shell:
        "{params.stringtie} --merge {input.stringtie_list} -o {output} "
        "-p {params.threads} -G {input.gff}"
       
       
rule stringtie_second_pass:
    input: 
        bam = "bams/{sample}_GRCh38.sorted.bam",
        gff = "stringtie_results/GRCh38.merged.gtf"
    output:
        assembled_transcripts = "stringtie_results/{sample}/{sample}_GRCh38.assembled_transcripts.secondpass.gtf",
        gene_abundances = "stringtie_results/{sample}/{sample}_GRCh38.gene_abundances.secondpass.txt",
        fully_covered_transcripts = "stringtie_results/{sample}/{sample}_GRCh38.fully_covered_transcripts.secondpass.gtf"
    threads: 4
    params:
        stringtie = STRINGTIE_PATH,
        threads = 4
    shell:
        "{params.stringtie} {input.bam} -p {params.threads} "
        "-G {input.gff} -B -e "
        "-o {output.assembled_transcripts} "
        "-A {output.gene_abundances} "
        "-C {output.fully_covered_transcripts}"



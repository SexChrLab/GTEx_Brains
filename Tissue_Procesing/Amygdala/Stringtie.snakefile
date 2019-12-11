configfile: "/scratch/mjpete11/GTEx/Configs/Sex_Sorted_Quantification.config.json"

import os

# TOOLS
# Variables should point to full paths to tools.
# Currently assumes all are available in one's PATH
# (e.g., installed with Bioconda)
STRINGTIE_PATH = "stringtie"

# Directory variables
BAM_DIRECTORY = "/mnt/storage/SAYRES/Mollie/GTEx_Brains/GTEx_Brains_092719/GTEx/Amygdala/Hisat_Stringtie/bams/"

# Samples
XX_SAMPLES = config["Amygdala_RNA_Trimmed"]["Female"]
XY_SAMPLES = config["Amygdala_RNA_Trimmed"]["Male"]
SAMPLES = XX_SAMPLES + XY_SAMPLES

rule all:
    input:
        expand("Annotated_Stringtie/{sample}/{sample}_GRCh38.fully_covered_transcripts.gtf", sample=SAMPLES)

rule stringtie:
    input:
        bam = "/mnt/storage/SAYRES/Mollie/GTEx_Brains/GTEx_Brains_092719/GTEx/Amygdala/Hisat_Stringtie/bams/{sample}_GRCh38.sorted.bam",
        gff = "/scratch/mjpete11/GTEx/gencode.v29.annotation.gtf"
    output:
        assembled_transcripts = "Annotated_Stringtie/{sample}/{sample}_GRCh38.assembled_transcripts.gtf",
        gene_abundances = "Annotated_Stringtie/{sample}/{sample}_GRCh38.gene_abundances.txt",
        fully_covered_transcripts = "Annotated_Stringtie/{sample}/{sample}_GRCh38.fully_covered_transcripts.gtf"
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

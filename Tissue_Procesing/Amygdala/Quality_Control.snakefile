configfile: "/scratch/mjpete11/GTEx/GTEx_Brain_v1.8.config"

from os.path import join

# TOOLS
# Variables should point to full paths to tools.
# Currently assumes all are available in one's PATH
# (e.g., installed with Bioconda)
FASTQC_PATH = "fastqc"
MULTIQC_PATH = "multiqc"
TRIMMOMATIC_PATH = "trimmomatic"

# Directories
FASTQ_DIRECTORY = "/data/storage/public/dbgap-8834/brain_amygdala" 

# Samples
XX_SAMPLES = config["Brain_Amygdala_RNA_Female"]
XY_SAMPLES = config["Brain_Amygdala_RNA_Male"]
SAMPLES = XX_SAMPLES + XY_SAMPLES


rule all:
        input: 
            "multiqc_trimmed_results/multiqc_report.html",
            expand("fastqc_results/{sample_female}_fixed_1_fastqc.html", sample_female = XX_SAMPLES), #Specify output of first rule that uses defined wilcards (e.g. "sample" is a regex, "sample_female" is constrained to a pre-defined list). 
            expand("fastqc_results/{sample_female}_fixed_2_fastqc.html", sample_female = XX_SAMPLES), #Have to do this so snakemake doesn't look for {male} samples in {female} directories.
            expand("fastqc_results/{sample_male}_fixed_1_fastqc.html", sample_male = XY_SAMPLES),
            expand("fastqc_results/{sample_male}_fixed_2_fastqc.html", sample_male = XY_SAMPLES)


rule fastqc_analysis:
        input:
            fq1 = os.path.join(FASTQ_DIRECTORY, "{sample}/{sample}_fixed_1.fastq"),
            fq2 = os.path.join(FASTQ_DIRECTORY, "{sample}/{sample}_fixed_2.fastq")
        output:
            ofq1 = "fastqc_results/{sample}_fixed_1_fastqc.html", #output files need to match input pattern; else use mv to overwrite name.
            ofq2 = "fastqc_results/{sample}_fixed_2_fastqc.html"
        params:
            FASTQC = FASTQC_PATH
        shell:
            "{params.FASTQC} -o fastqc_results {input.fq1} {input.fq2}"


rule multiqc:
        input:
            expand("fastqc_results/{sample}_{num}_fastqc.html", sample=SAMPLES, num=[1, 2])
        output:
            "multiqc_results/multiqc_report.html"
        params:
            MULTIQC = MULTIQC_PATH
        shell:
            "{params.MULTIQC} --interactive fastqc_results -o multiqc_results"


rule trimmomatic:
        input:
            fq1 = os.path.join(FASTQ_DIRECTORY, "{sample}/{sample}_fixed_1.fastq"),
            fq2 = os.path.join(FASTQ_DIRECTORY, "{sample}/{sample}_fixed_2.fastq"),
            ADAPTER_FASTA =  "/scratch/mjpete11/GTEx/adapter_sequences.fa"
        output:
            paired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_1.fastq.gz",
            paired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_2.fastq.gz",
            unpaired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_unpaired_1.fastq.gz",
            unpaired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_unpaired_2.fastq.gz",
            logfile = "logfiles/{sample}_trimmomatic.log"
        params:
            threads = 4,
            seed_mismatches = 2,
            palindrome_clip_threshold = 30,
            simple_clip_threshold = 10,
            leading = 3,
            trailing = 3,
            winsize = 4,
            winqual = 30,
            minlen = 50
        shell:
            "trimmomatic PE -threads {params.threads} -trimlog {output.logfile} "
            "{input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} "
            "{output.paired_2} {output.unpaired_2} "
            "ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} "
            "LEADING:{params.leading} TRAILING:{params.trailing} "
            "SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}"


rule fastqc_analysis_trimmed:
        input:
            expand("trimmed_fastqs/{{sample}}_trimmomatic_trimmed_paired_{num}.fastq.gz", num=[1, 2])
        output:
            ofq1 = "fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_1_fastqc.html",
            ofq2 = "fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_2_fastqc.html"
        params:
            FASTQC = FASTQC_PATH
        shell:
            "{params.FASTQC} -o fastqc_trimmed_results {input}"


rule multiqc_trimmed_paired:
        input:
            expand("fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_{num}_fastqc.html", sample=SAMPLES, num=[1, 2])
        output:
            "multiqc_trimmed_results/multiqc_report.html"
        params:
            MULTIQC = MULTIQC_PATH
        shell:
            "{params.MULTIQC} --interactive fastqc_trimmed_results -o multiqc_trimmed_results"


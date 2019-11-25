from os.path import join

# Config contains list of files that passed quality control.
configfile: "Quantification.config.json"

# Directories
INPUT_DIR = "/scratch/mjpete11/GTEx/Hippocampus/trimmed_fastqs/"
SALMON_DIR = "/scratch/mjpete11/GTEx/Hippocampus/salmon_quants/"
HISAT_DIR = "/scratch/mjpete11/GTEx/Hippocampus/hisat_quants/"

# Tools
SALMON = "salmon"

# Index
SALMON_INDEX = config["Salmon_hg38_transcriptome_index_path"]

# Samples
SAMPLES = config["Hippocampus_RNA_Trimmed"]

rule all:
	input:
	    expand(SALMON_DIR + "{sample}", sample=SAMPLES)

rule salmon_quant_paired:
	input:
		fq1 = os.path.join(INPUT_DIR, "{sample}_trimmomatic_trimmed_paired_1.fastq.gz"),
		fq2 = os.path.join(INPUT_DIR, "{sample}_trimmomatic_trimmed_paired_2.fastq.gz")
	output:
		OUTPUT = os.path.join(SALMON_DIR, "{sample}"),
	params:
		SALMON = SALMON,
		SALMON_INDEX = SALMON_INDEX,
		LIBTYPE = "A", # LIBTYPE A for automatic detection of library type
		threads = 8
	run:
		shell("{params.SALMON} quant -i {params.SALMON_INDEX} -l {params.LIBTYPE} -1 {input.fq1} -2 {input.fq2} -o {output.OUTPUT}")
		


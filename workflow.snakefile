# snakefile to execute all scripts

import os

# Constants
BASE = "/scratch/mjpete11/human_monkey_brain/"

rule all:
	input: 
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/Isc_Dim4.pdf")

# Rules
# Rule 1: Generate metadata to use to subset counts
rule metadata:
	input: 
		os.path.join(BASE, "data/metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
		os.path.join(BASE, "data/metadata/GTEx_Analysis_v8_Anndata/metadata/otations_SubjectPhenotypesDS.txt")
	output:
		os.path.join(BASE, "data/metadata/older_metadata.csv"),
		os.path.join(BASE, "data/metadata/metadata.csv")
	script:
		"/scratch/mjpete11/human_monkey_brain/data/metadata/metadata.R"

# Rule 2: Gene expression plots to determine filtering threshold(s)
rule plot_genes:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
	output:
		os.path.join(BASE, "data_exploration/plots/sex_filtered_hist.pdf"),
		os.path.join(BASE, "data_exploration/plots/region_TPM1_hist.pdf"),
		os.path.join(BASE, "data_exploration/plots/region_TPM5_hist.pdf"),
		os.path.join(BASE, "data_exploration/plots/region_TPM10_hist.pdf"),
		os.path.join(BASE, "data_exploration/plots/filter_by_region.csv")
	script:
		"/scratch/mjpete11/human_monkey_brain/data_exploration/gene_expression_plots.R"

# Rule 3: multidimensional scaling by region
rule mds:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"),
		os.path.join(BASE, "data/gene_annotation/gencodeGenes_Xchr.txt"),
		os.path.join(BASE, "data/gene_annotation/gencodeGenes_Ychr.txt")
	output:
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/Sex_Dim12.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/RIN_Dim1.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/RIN_Dim2.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/RIN_Dim3.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/RIN_Dim4.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/Isc_Dim1.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/Isc_Dim2.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/Isc_Dim3.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/plots/Isc_Dim4.pdf")
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/multidimensional_scaling/multidimensional_scaling.R"

# Rule 4: uniform manifold approximation and projection
rule umap:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"),
		os.path.join(BASE, "data/gene_annotation/gencodeGenes_Xchr.txt"),
		os.path.join(BASE, "data/gene_annotation/gencodeGenes_Ychr.txt")
	output:
		os.path.join(BASE, "dimension_redcution/umap_plots/umap_NoSexChr.txt")
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/umap/umap.R"

# Rule 5: limma
#rule limma:
#	input:
#		os.path.join(BASE, "data/metadata/metadata.csv"),
#		os.path.join(BASE, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
#	output:
## Tables
#	script:
#		"/scratch/mjpete11/human_monkey_brain/limma/limma.R"

# Rule 6: multivariate adaptive shrinkage
# rule mashr:
# 	input: 

	
		


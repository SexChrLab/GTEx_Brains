# snakefile to execute all scripts

import os

# Constants
BASE = "/scratch/mjpete11/human_monkey_brain/"

rule all:
	input: 
		os.path.join(BASE, "dimension_reduction/umap/no_sex/no_sex_projection.csv") 

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

# Rule 4: Write list of genes determined via various filtering thresholds

# Rule 5: Write filtered/normalized count expression matrices 

# Rule 6: multidimensional scaling with sex chr
rule mds_with_sex:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv"),
	output:
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/Sex_Dim12.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/RIN_Dim1.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/RIN_Dim2.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/RIN_Dim3.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/RIN_Dim4.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/Isc_Dim1.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/Isc_Dim2.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/Isc_Dim3.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/sex_chr/Isc_Dim4.pdf")
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/multidimensional_scaling/multidimensional_scaling.R"

# Rule 7: multidimensional scaling without sex chr
rule mds_no_sex:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_no_sex_chr.csv"),
	output:
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/Sex_Dim12.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/RIN_Dim1.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/RIN_Dim2.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/RIN_Dim3.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/RIN_Dim4.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/Isc_Dim1.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/Isc_Dim2.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/Isc_Dim3.pdf"),
		os.path.join(BASE, "dimension_reduction/multidimensional_scaling/no_sex/Isc_Dim4.pdf")
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/multidimensional_scaling/multidimensional_scaling.R"

# Rule 5: uniform manifold approximation and projection; write projections
rule umap_with_sex:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv"),
	output:
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/sex_projection.csv"),
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/umap/umap.R"

rule umap_no_sex:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_no_sex_chr.csv"),
	output:
		os.path.join(BASE, "dimension_reduction/umap/no_sex/no_sex_projection.csv"),
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/umap/umap.R"

# Rule 6: umap plots
rule umap_plots:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/sex_projection.csv"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/no_sex_projection.csv")
	output:
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/umap_tissue.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_tissue.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/umap_sex.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_sex.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/umap_age.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_age.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/umap_isc.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_isc.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/umap_rin.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_rin.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/umap_indiv.pdf"),
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_indiv.pdf")
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/umap/umap_plots.R"


# Rule 6: limma
# This produces 11 tables
rule limma:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
	output:
	script:
		"/scratch/mjpete11/human_monkey_brain/limma/limma.R"

#  Rule 7: multivariate adaptive shrinkage
# rule mashr:
# 	input: 

	
		


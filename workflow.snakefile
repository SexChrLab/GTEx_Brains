# Snakefile to execute all scripts

import os

# Constants
BASE = "/scratch/mjpete11/human_monkey_brain/"

# Output
rule all:
	input: 
		os.path.join(BASE, "dimension_reduction/umap/no_sex/umap_indiv.pdf")

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

# Rule 3: Write list of genes determined via various filtering thresholds
 rule write_genes:
 	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
	output:
		os.path.join(BASE, "data/expression_matrices/output/filtered_by_sex.csv")
		os.path.join(BASE, "data/expression_matrices/output/region_tally.csv")
		os.path.join(BASE, "data/expression_matrices/output/filtered_by_region.csv")
		os.path.join(BASE, "data/expression_matrices/output/union_region_filtered.csv")
	script:
		"/scratch/mjpete11/human_monkey_brain/data/expression_matrices/tpm_filtering.R"
		
# Rule 4: Write filtered/normalized count expression matrices 
rule process_counts:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv")
		os.path.join(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
		os.path.join(BASE, "data/gene_annotation/gencodeGenes_Xchr.txt")
		os.path.join(BASE, "data/gene_annotation/gencodeGenes_Ychr.txt")
		os.path.join(BASE, "data/expression_matrices/output/union_region_filtered.csv")
	output:
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv")
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_no_sex_chr.csv")
	script:
		"scratch/mjpete11/human_monkey_brain/data/expression_matrices/expression_matrices.R""

# Rule 5: multidimensional scaling with sex chr
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

# Rule 6: multidimensional scaling without sex chr
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

# Rule 7: umap; write projections from counts with sex chr
rule umap_with_sex:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv"),
	output:
		os.path.join(BASE, "dimension_reduction/umap/sex_chr/sex_projection.csv"),
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/umap/umap.R"

# Rule 8: umap; write projections from counts without sex chr
rule umap_no_sex:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "data/expression_matrices/output/processed_counts_no_sex_chr.csv"),
	output:
		os.path.join(BASE, "dimension_reduction/umap/no_sex/no_sex_projection.csv"),
	script:
		"/scratch/mjpete11/human_monkey_brain/dimension_reduction/umap/umap.R"

# Rule 9: umap plots; with and without sex chr; color each point be feature
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

# Rule 10: limma
# This produces 11 tables
rule limma:
	input:
		os.path.join(BASE, "data/metadata/metadata.csv"),
		os.path.join(BASE, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
	output:
	script:
		"/scratch/mjpete11/human_monkey_brain/limma/limma.R"

#  Rule 11: multivariate adaptive shrinkage
# rule mashr:
# 	input: 

	
		


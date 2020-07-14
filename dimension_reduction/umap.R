#!/usr/bin/env Rscript

# Uniform manifold approximation and projection for dimension reduction

# Load libraries
library(umap)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(stringr)
library(edgeR)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")

# Read in files
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_counts <- data.frame(fread(COUNTS))

#_______________________________________________________________________________
# Sample pre-processing
#_______________________________________________________________________________
# Drop gene name and ID from gene_counts df
gene_counts <- gene_counts[, 3:ncol(gene_counts)]

# Replace . to - in colnames
colnames(gene_counts) <- str_replace_all(colnames(gene_counts),
										 pattern = "\\.", "-")

# Number of samples in count df
ncol(gene_counts) # 17,382

# Number of samples in meta (brain only)
nrow(meta) # 2,570

# Drop samples in metadata that do not have count data
select_samples <- colnames(gene_counts)[colnames(gene_counts) %in% meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples, ]

# Number of samples in meta
nrow(meta) # 2,146

# Set rownames of metadata object equal to sample names
rownames(meta) <- meta$Sample_ID

# Subset gene_counts to only samples present in metadata
gene_counts <- gene_counts[, select_samples]

# Number of samples
ncol(gene_counts) # 2,146

# Check that the count and meta data have the same samples in the same order
identical(colnames(gene_counts), rownames(meta)) # TRUE

#_______________________________________________________________________________
# Filter genes by expression in each sex and voom normalize 
#_______________________________________________________________________________
# Make factor indicating sex of samples for filtering
#sex <- factor(meta$Sex)

# Make design matrix to use with voom
design <- model.matrix(~meta$Sex)
rownames(design) <- colnames(gene_counts)

# How many genes are there before filtering?
nrow(gene_counts) # 56,200

# Remove genes with cpm < 1 in each sex
keep <- filterByExpr(gene_counts, design = design,
					 min.count = 1, min.prop = 0.5)
gene_counts <- gene_counts[keep, ]

# How many genes are left after filtering?
nrow(gene_counts) # 35,610

# limma-voom normalization
gene_counts <- voom(gene_counts, design=design)

# Convert gene_counts (EList object) to matrix for use w/ umap
gene_counts <- data.matrix(gene_counts)

#_______________________________________________________________________________
# umap 
#_______________________________________________________________________________
# gene_counts is an expression matrix, with genes as rows and samples (libraries) as columns. 
# For this particular matrix, we used the partial residuals of a model after removing technical covariates. 
# If this doesn't apply to your data, use your voom-normalized expression matrix
# n_neighbors and min_dist are flexible parameters: 
# see https://umap-learn.readthedocs.io/en/latest/parameters.html for explanation
start_t <- Sys.time()
a <- umap(gene_counts, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time:", end_t-start_t))

# Write umap object to file
write.table(a$layout, file="umap.txt", sep=" ")

# extract the projections (a$layout) and join it with metadata (particularly the 'Region' column)
# full.covariates is a vector of model covariate names (in case you want to plot by those as well)
#a.umap = data.frame(as.data.frame(a$layout),meta[,c(full.covariates,'Region')])
#
## Plot UMAP projections and color by Region.
## region.colors is a custom vector of pretty colors.
## batch.variable is a length=1 character vector that contains the name of our batch variable
#p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape=batch.variable)) +
#	geom_point(size=1.25) +
#	scale_color_manual(values=region.colors) +
#	theme_classic(base_size=18) +
#	xlab('UMAP 1') +
#	ylab('UMAP 2') +
#	coord_fixed() +
#	theme(axis.ticks=element_blank(),axis.text=element_blank())
#ggsave(p,file='figures/data_visualization_dimension_reduction_umap_batched.pdf',useDingbats=FALSE)
#
## Redo the plot, this time without highlighting batch
#p = ggplot(a.umap,aes_string('V1','V2',color='Region')) +
#	geom_point(size=1.25) +
#	scale_color_manual(values=region.colors) +
#	theme_classic(base_size=18) +
#	xlab('UMAP 1') +
#	ylab('UMAP 2') +
#	coord_fixed() +
#	theme(axis.ticks=element_blank(),axis.text=element_blank())
#ggsave(p,file='figures/data_visualization_dimension_reduction_umap.pdf',useDingbats=FALSE)
#
## Redo the plot, this time showing color AND shape to better differentiate points since many colors look similar
#p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape='Region')) +
#	geom_point(size=1.25) +
#	scale_color_manual(values=region.colors) +
#	scale_shape_manual(values=region.shapes) +
#	theme_classic(base_size=18) +
#	xlab('UMAP 1') +
#	ylab('UMAP 2') +
#	coord_fixed() +
#	guides(shape = guide_legend(override.aes = list(size = 3))) +
#	theme(axis.ticks=element_blank(),axis.text=element_blank())
#ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes.pdf',useDingbats=FALSE)

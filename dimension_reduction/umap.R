#!/usr/bin/env Rscript

# Uniform manifold approximation and projection for dimension reduction

# Load libraries
library(umap)
library(ggplot2)
library(data.table)
library(stringr)
library(edgeR)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
CHRX <- file.path(BASE, "data/gene_annotation/gencodeGenes_Xchr.txt")
CHRY <- file.path(BASE, "data/gene_annotation/gencodeGenes_Ychr.txt")

# Output
PROJECTION <- file.path(BASE, "dimension_redcution/umap_plots/umap_NoSexChr.txt")

# Read in files
xchr <- read.table(CHRX, sep = "")
ychr <- read.table(CHRY, sep = "")
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_counts <- data.frame(fread(COUNTS))

#_______________________________________________________________________________
# Drop X and Y-linked genes to check if clustering is driven by sex-linked loci
#_______________________________________________________________________________
# Original number of genes
nrow(gene_counts) # 56,200

# Remove X and Y-linked genes
rows_to_drop <- intersect(xchr$V6, gene_counts$Name)
rows_to_drop <- c(rows_to_drop, intersect(ychr$V6, gene_counts$Name))

gene_counts <- gene_counts[-which(gene_counts[, 1] %in% rows_to_drop), ]

# Are the expected number of genes remaining?; Yes
nrow(gene_counts) # 53,325
nrow(xchr) + nrow(ychr) # 2,875

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
nrow(gene_counts) # 34,341

# limma-voom normalization
gene_counts <- voom(gene_counts, design = design)

# Convert gene_counts (EList object) to matrix for use w/ umap
gene_counts <- data.matrix(gene_counts)

#_______________________________________________________________________________
# umap
#_______________________________________________________________________________
# gene_counts is an expression matrix: genes as rows and samples as columns.
start_t <- Sys.time()
a <- umap(gene_counts, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time:", end_t - start_t))

# Write umap object to file
write.table(a$layout, file = PROJECTION, sep = " ")

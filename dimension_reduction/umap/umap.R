#!/usr/bin/env Rscript

# Script to plot the umap projections

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)

# Input
METADATA <- snakemake@input[[1]]
COUNTS <- snakemake@input[[2]]

# Output
PROJ <- snakemake@output[[1]]

# Read in files
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
counts <- read.table(COUNTS, header = TRUE, sep = ",", stringsAsFactors = FALSE)

#_______________________________________________________________________________
# umap
#_______________________________________________________________________________
# Replaces '.' to '-' in the sample IDs for the projections
colnames(counts) <- str_replace_all(colnames(counts), pattern = "\\.", "-")

# Convert gene_counts to matrix and transpose for use w/ umap
res_matrix <- as.matrix(t(counts))

# Apply umap to an expression matrix: genes as rows and samples as columns.
# processed count matrix with sex chr
start_t <- Sys.time()
res_umap <- umap(res_matrix, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time to do umap with sex chr counts:", end_t - start_t))

# Make df of projections from umap object
res_proj <- data.frame(res_umap$layout)

#_______________________________________________________________________________
# Sanity check 
#_______________________________________________________________________________
# Number of samples in meta
nrow(meta) # 2,146

# Number of samples in projection
nrow(res_proj) == nrow(meta) # TRUE

# Check that the samples in the projections are in the same order as in meta
identical(row.names(res_proj), meta$Sample_ID) # TRUE

# Write umap objects to file
write.csv(res_proj, file = PROJ, row.names = TRUE)


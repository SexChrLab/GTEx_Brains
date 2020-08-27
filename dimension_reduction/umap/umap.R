#!/usr/bin/env Rscript

# Script to plot the umap projections

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
SEX_COUNTS <- file.path(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv")
NO_SEX_COUNTS <- file.path(BASE, "data/expression_matrices/output/processed_counts_no_sex_chr.csv")

# Output
SEX_PROJ <- file.path(BASE, "dimension_reduction/umap/sex_chr/sex_projection.csv")
NO_SEX_PROJ <- file.path(BASE, "dimension_reduction/umap/no_sex/no_sex_projection.csv")

# Read in files
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
sex_counts <- read.table(SEX_COUNTS, header = TRUE, sep = ",", stringsAsFactors = FALSE)
no_sex_counts <- read.table(NO_SEX_COUNTS, header =TRUE, sep = ",", stringsAsFactors = FALSE)

#_______________________________________________________________________________
# umap
#_______________________________________________________________________________
# Replaces '.' to '-' in the sample IDs for the projections
colnames(sex_counts) <- str_replace_all(colnames(sex_counts), pattern = "\\.", "-")
colnames(no_sex_counts) <- str_replace_all(colnames(no_sex_counts), pattern = "\\.", "-")

# Convert gene_counts to matrix and transpose for use w/ umap
sex_matrix <- as.matrix(t(sex_counts))
no_sex_matrix <- as.matrix(t(no_sex_counts))

# Apply umap to an expression matrix: genes as rows and samples as columns.
# processed count matrix with sex chr
start_t <- Sys.time()
sex_umap <- umap(sex_matrix, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time to do umap with sex chr counts:", end_t - start_t))

# processed count matrix without sex chr
start_t <- Sys.time()
no_sex_umap <- umap(no_sex_matrix, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time to do umap without sex chr counts:", end_t - start_t))

# Make df of projections from umap object
sex_proj <- data.frame(sex_umap$layout)
no_sex_proj <- data.frame(no_sex_umap$layout)

#_______________________________________________________________________________
# Sanity check 
#_______________________________________________________________________________
# Number of samples in meta
nrow(meta) # 2,146

# Number of samples in projection
nrow(sex_proj) == nrow(meta) # TRUE
nrow(no_sex_proj) == nrow(meta) # TRUE

# Check that the samples in the projections are in the same order as in meta
identical(row.names(sex_proj), meta$Sample_ID) # TRUE
identical(row.names(no_sex_proj), meta$Sample_ID) # TRUE

# Write umap objects to file
write.csv(sex_proj, file = SEX_PROJ, row.names = TRUE)
write.csv(no_sex_proj, file = NO_SEX_PROJ, row.names = TRUE)


#!/usr/bin/env Rscript

# Purpose: Generate filtered and normalized  count matrices
# Write two expression matrices, with and without sex chr linked genes

# Load libraries
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
REGION <- file.path(BASE, "data/expression_matrices/output/union_region_filtered.csv")
# SEX <- file.path(BASE, "data/expression_matrices/output/filtered_by_sex.csv")

# Output
TABLE1 <- file.path(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv")
TABLE2 <- file.path(BASE, "data/expression_matrices/output/processed_counts_no_sex_chr.csv")

# Read in data
xchr <- read.table(CHRX, sep = "")
ychr <- read.table(CHRY, sep = "")
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_counts <- data.frame(fread(COUNTS))
region <- read.csv(REGION)
# sex <- read.csv(SEX) # In case we want to filter by sex instead

#_______________________________________________________________________________
# Drop X and Y-linked genes to check if clustering is driven by sex-linked loci
#_______________________________________________________________________________
# Original number of genes
nrow(gene_counts) # 56,200

# Remove X and Y-linked genes
rows_to_drop <- intersect(xchr$V6, gene_counts$Name)
rows_to_drop <- c(rows_to_drop, intersect(ychr$V6, gene_counts$Name))

counts_no_sex <- gene_counts[-which(gene_counts[, 1] %in% rows_to_drop), ]

# Are the expected number of genes remaining?; Yes
nrow(counts_no_sex) # 53,325
nrow(xchr) + nrow(ychr) # 2,875

#_______________________________________________________________________________
# Sample pre-processing
#_______________________________________________________________________________
# Drop gene name and ID from gene_counts df
#gene_counts <- gene_counts[, 3:ncol(gene_counts)]
#counts_no_sex <- counts_no_sex[, 3:ncol(counts_no_sex)]

# Replace . to - in colnames
colnames(gene_counts) <- str_replace_all(colnames(gene_counts),
                                         pattern = "\\.", "-")
colnames(counts_no_sex) <- str_replace_all(colnames(counts_no_sex),
                                         pattern = "\\.", "-")

# Number of samples in count df
ncol(gene_counts) # 17,382
ncol(counts_no_sex) # 17,382

# Number of samples in meta (brain only)
nrow(meta) # 2,146

# Define samples to keep from metadata
select_samples <- c("Name", "Description", colnames(gene_counts)[colnames(gene_counts) %in% meta$Sample_ID])

# Set rownames of metadata object equal to sample names
rownames(meta) <- meta$Sample_ID

# Subset gene_counts to only samples present in metadata
gene_counts <- gene_counts[, select_samples]
counts_no_sex <- counts_no_sex[, select_samples]

# Number of samples
ncol(gene_counts) - 2 # 2,146
ncol(counts_no_sex) - 2 # 2,146

# Check that the count and meta data have the same samples in the same order
identical(colnames(gene_counts)[-c(1,2)], rownames(meta)) # TRUE
identical(colnames(counts_no_sex)[-c(1,2)], rownames(meta)) # TRUE

#_______________________________________________________________________________
# Filter genes by expression in each sex and voom normalize
#_______________________________________________________________________________
# Make design matrix to use with voom
sex_design <- model.matrix(~meta$Sex)
no_sex_design <- model.matrix(~meta$Sex)
rownames(sex_design) <- colnames(gene_counts)[-c(1,2)]
rownames(no_sex_design) <- colnames(counts_no_sex)[-c(1,2)]

# How many genes are there before filtering?
nrow(gene_counts) # 56,200
nrow(counts_no_sex) # 53,325

# Originally used cpm to filter; keeping in case I need it...
# Remove genes with cpm < 1 in each sex
# keep <- filterByExpr(gene_counts, design = design,
#                     min.count = 1, min.prop = 0.5)
#gene_counts <- gene_counts[keep, ]

# Use predefined lists of genes to keep to filter
gene_counts <- gene_counts[which(gene_counts$Name %in% region$gene_id), ]
counts_no_sex <- counts_no_sex[which(counts_no_sex$Name %in% region$gene_id), ]

# How many genes are left after filtering?
nrow(gene_counts) # 13,955
nrow(counts_no_sex) # 13,468

# limma-voom normalization
gene_counts <- voom(gene_counts[3:ncol(gene_counts)], design = sex_design)
counts_no_sex <- voom(counts_no_sex[3:ncol(counts_no_sex)], design = no_sex_design)

# Write to file
write.csv(gene_counts$E, TABLE1, row.names = FALSE)
write.csv(counts_no_sex$E, TABLE2, row.names = FALSE)


#!/usr/bin/env Rscript

# Purpose:
# 1) Filter samples by sex with median(TPM >= 1)
# 2) Filter samples by regions with median(TPM >= 1, 5, and 10)
# 3) Write df with list of genes to drop depending on filtering method;
#    Will use this to drop genes if analyzing count data

# Load libraries
library(data.table)
library(stringr)
library(zoo)
library(purrr)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
TPM <- file.path(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")

# Output
DF1 <- file.path(BASE, "data/expression_matrices/output/filtered_by_sex.csv")
DF2 <- file.path(BASE, "data/expression_matrices/output/region_tally.csv")
DF3 <- file.path(BASE, "data/expression_matrices/output/filtered_by_region.csv")
DF4 <- file.path(BASE, "data/expression_matrices/output/union_region_filtered.csv")

# Read in files
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
unfiltered_mat <- data.frame(fread(TPM))

#_______________________________________________________________________________
# Sample pre-processing
#_______________________________________________________________________________
# Store gene name and ID columns (will be dropped and needs to be appended)
genes <- unfiltered_mat[,1:2]

# Replace . to - in colnames
colnames(unfiltered_mat) <- str_replace_all(colnames(unfiltered_mat), pattern = "\\.", "-")

# Number of samples in count df
ncol(unfiltered_mat) # 17,384

# Number of samples in meta (brain only)
nrow(meta) # 2,570

# Drop samples in metadata that do not have count data
select_samples <- colnames(unfiltered_mat)[colnames(unfiltered_mat) %in% meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples, ]

# Number of samples in meta
nrow(meta) # 2,146

# Set rownames of metadata object equal to sample names
rownames(meta) <- meta$Sample_ID

# Subset unfiltered_mat to only samples present in metadata
unfiltered_mat <- unfiltered_mat[, select_samples]

# Number of samples
ncol(unfiltered_mat) # 2,146

# Check that the count and meta data have the same samples in the same order
identical(colnames(unfiltered_mat), rownames(meta)) # TRUE

# Append gene name and gene ID columns back to the expression matrix
unfiltered_mat <- cbind(genes, unfiltered_mat)

#_______________________________________________________________________________
# Filter genes by expression in each sex 
#_______________________________________________________________________________
# How many genes  are there before filtering?
nrow(unfiltered_mat) # 56,200

# Get list of female sample IDs 
fem_samples <- meta$Sample_ID[which(meta$Sex == "Female")]
male_samples <- meta$Sample_ID[which(meta$Sex == "Male")]

# Split expression matrix into two by sex
fem_counts <- unfiltered_mat[, c(1:2, which(colnames(unfiltered_mat) %in% fem_samples))]
male_counts <- unfiltered_mat[, c(1:2, which(colnames(unfiltered_mat) %in% male_samples))]

# Median filter: Function to filter rows that have a median >= condition
median_filter <- function(DF, TPM) {
    DF <- DF[which(apply(DF[, -c(1:2)], 1, median) >= TPM), ]
    return(DF)
}

# Remove genes with TPM < 5
fem_counts <- median_filter(DF = fem_counts, TPM = 1)
male_counts <- median_filter(DF = male_counts, TPM = 1)

# How many genes are left after filtering in each sex?
nrow(fem_counts) # 15,611
nrow(male_counts) # 15,569

#_______________________________________________________________________________
# Make df listing genes to drop depending on filter 
#_______________________________________________________________________________
# Make df of gene name and ID to drop in each sex 
fem_genes_id <- fem_counts$Name
fem_genes_name <- fem_counts$Description
male_genes_id <- male_counts$Name
male_genes_name <- male_counts$Description

# Set all vectors to the max vector length;
# This will allow cbind to fill empty cells with NA 
n <- max(length(fem_genes_id), length(male_genes_id))
length(fem_genes_id) <- n
length(fem_genes_name) <- n
length(male_genes_id) <- n
length(male_genes_name) <- n
sex_filtered_genes <- cbind(fem_genes_id, fem_genes_name,
                            male_genes_id, male_genes_name)
# Write filtered_by_sex.csv
write.csv(sex_filtered_genes, DF1, row.names=FALSE)

#_______________________________________________________________________________
# In case I want the TPM gene expression matrix of union of both sexes  
#_______________________________________________________________________________
# Make a gene expression matrix that is the union of m/f genes remaining
# Use zoo to combine dfs with different num rows and different cols
# The larger df has to go first
sex_filtered_mat <- data.frame(fem_counts, cbind(zoo(, 1:nrow(fem_counts)),
                                                 as.zoo(male_counts)))

# Replace . to - in colnames
colnames(sex_filtered_mat) <- str_replace_all(colnames(sex_filtered_mat), 
                                              pattern = "\\.", "-")

# Check for the expected number og rows/cols
nrow(sex_filtered_mat) == nrow(fem_counts) # TRUE 
ncol(sex_filtered_mat) == ncol(fem_counts) + ncol(male_counts) # TRUE

# Write to file
# write.csv(sex_filtered_mat, 
#            "/scratch/mjpete11/human_monkey_brain/data/expression_matrices/output/filtered_by_sex.csv")

#_______________________________________________________________________________
# Keep genes with median TPM >= 1, 5, 10 in each region 
#_______________________________________________________________________________
# For each tissue, make a df of samples from that tissue and store in a list
types <- unique(meta$Tissue)

tissue_lst <- list()
for (i in types) {
    tissue_lst[[i]] <- meta[meta$Tissue == i,]
}

# Sort expression mat into list of dfs by tissue type
sort_counts <- function(x) {
	lst <- which(colnames(unfiltered_mat) %in% x$Sample_ID)
	res <- unfiltered_mat[, c(1:2, lst)]
	return(res)
}
counts_by_tissue <- lapply(tissue_lst, sort_counts)

# Apply filtering threshold
TPM_5 <- Map(median_filter, TPM = 5, counts_by_tissue)

# Name each df in list by  which tissue type it is
names(TPM_5) <- names(tissue_lst)

# How many genes left per regions with TPM >= 5?
res <- data.frame(regions = names(TPM_5),
				  num_genes_remaining = as.numeric(lapply(TPM_5, function(x) nrow(x))))

write.csv(res, DF2, row.names = FALSE)
# Code to write results to file, in case I want it at some point
# sapply(names(TPM_5), function(x) write.table(TPM_5[[x]], 
#       file=paste0("/scratch/mjpete11/human_monkey_brain/data/expression_matrices/output", x, ".csv"))) 

#_______________________________________________________________________________
# Write df where cols list which genes to keep in each region with TPM >= 5 
#_______________________________________________________________________________
result <- do.call(rbind, lapply(1:11, function(i){
				  data.frame(TPM_5[[i]][1], TPM_5[[i]][2], names(TPM_5)[i])}))

colnames(result) <- c("gene_id", "gene_name", "region")
write.csv(result, DF3, row.names=FALSE)

#_______________________________________________________________________________
# Write df of the union of genes with TPM >= 5 across regions 
#_______________________________________________________________________________
lst <- lapply(TPM_5, function(x) x$Name)
res <- Reduce(union, lst)
res <- data.frame(gene_id = res)
write.csv(res, DF4, row.names=FALSE)

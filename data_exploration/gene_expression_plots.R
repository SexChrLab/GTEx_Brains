#!/usr/bin/env Rscript

# Purpose:
# 1) Filter gene expression matrix by various methods 
# 2) Plot histograms of gene TPMs to determine filtering parameters

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
SEX_PLT <- file.path(BASE, "data_exploration/plots/sex_filtered_hist.pdf")
REGION_PLT_1 <- file.path(BASE, "data_exploration/plots/region_TPM1_hist.pdf")
REGION_PLT_5 <- file.path(BASE, "data_exploration/plots/region_TPM5_hist.pdf")
REGION_PLT_10 <- file.path(BASE, "data_exploration/plots/region_TPM10_hist.pdf")

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
fem_counts <- unfiltered_mat[, c(which(colnames(unfiltered_mat) 
											%in% fem_samples))]
male_counts <- unfiltered_mat[, c(which(colnames(unfiltered_mat) 
											 %in% male_samples))]

# Function to filter rows that do not meet condition
filter_func <- function(DF, TPM) {
	DF <- DF[apply(DF, 1, function(x) all(x > TPM)),]
	return(DF)
}

# Remove genes with TPM < 1
fem_counts <- filter_func(DF = fem_counts, TPM = 1)
male_counts <- filter_func(DF = male_counts, TPM = 1)

# How many genes are left after filtering?
nrow(fem_counts) # 7,175
nrow(male_counts) # 6,136

# Make a gene expression matrix that is the union of m/f genes remaining
# Use zoo to combine dfs with different num rows and different cols
sex_filtered_mat <- data.frame(fem_counts, cbind(zoo(1:nrow(fem_counts)), 
							   as.zoo(male_counts)))

#_______________________________________________________________________________
# Plot histogram of the median log(TPM) across all regions 
#_______________________________________________________________________________
# Plot
pdf(SEX_PLT)
hist(log10(apply(sex_filtered_mat, 1, median)), 100,
	main = "Histogram of gene expression filtered by sex", 
   	xlab = "Mean log10(TPM)",
	ylab = "Frequency")
dev.off()

#_______________________________________________________________________________
# Filter genes by expression (TPM > 1, 5, 10) in each region 
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
	res <- unfiltered_mat[, c(lst)]
	return(res)
}
counts_by_tissue <- lapply(tissue_lst, sort_counts)

# Apply different filtering thresholds
TPM_1 <- Map(filter_func, TPM = 1, counts_by_tissue)
TPM_5 <- Map(filter_func, TPM = 5, counts_by_tissue)
TPM_10 <- Map(filter_func, TPM = 10, counts_by_tissue)

#_______________________________________________________________________________
# Plot histogram of the median log(TPM) in each region separately 
#_______________________________________________________________________________
# Function to plot histograms on one page
plot_func <- function(mat, tissue, thresh) {
	hist(log10(apply(mat, 1, median)), 100, main = tissue)
    mtext(paste("Histogram of gene expression with TPM >", thresh, "filtered by region"), 
		  side = 3, outer = TRUE,  cex=1.2, line=3)
	mtext("Median log10(TPM)", side = 1, outer = TRUE, cex = 0.8, line = 1)
	mtext("Frequency", side = 2, outer = TRUE, cex = 0.8, line = 2)
}

# Plot filter by regions with TPM > 1
pdf(REGION_PLT_1)
par(mfrow = c(4, 3), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
pmap(list(mat = TPM_1, tissue = names(tissue_lst)), plot_func, thresh = 1) 
dev.off()

# Plot filter by regions with TPM > 5
pdf(REGION_PLT_5)
par(mfrow = c(4, 3), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
pmap(list(mat = TPM_5, tissue = names(tissue_lst)), plot_func, thresh = 5) 
dev.off()

# Plot filter by regions with TPM > 10
pdf(REGION_PLT_10)
par(mfrow = c(4, 3), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
pmap(list(mat = TPM_10, tissue = names(tissue_lst)), plot_func, thresh = 10) 
dev.off()

#!/usr/bin/env Rscript

# Purpose:
# 1) Filter samples with meidan(TPM > 1) in each sex separately
# 2) Filter samples by regions with median(TPM > 1, 5, and 10)
# 3) Plot frequency of log10(median(TPM)) 

# Load libraries
library(data.table)
library(stringr)
library(zoo)
library(purrr)

# Snakemake constants
# Input
METADATA <- snakemake@input[[1]]
TPM <- snakemake@input[[2]]

# Output
SEX_PLT <- snakemake@output[[1]]
REGION_PLT_1 <- snakemake@output[[2]]
REGION_PLT_5 <- snakemake@output[[3]]
REGION_PLT_10 <- snakemake@output[[4]]
TABLE <- snakemake@output[[5]]

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
fem_counts <- unfiltered_mat[, c(which(colnames(unfiltered_mat) %in% fem_samples))]
male_counts <- unfiltered_mat[, c(which(colnames(unfiltered_mat) %in% male_samples))]

# Median filter: Function to filter rows that do not have a median >= condition
median_filter <- function(DF, TPM) {
    DF <- DF[which(apply(DF[, -c(1:2)], 1, median) > TPM), ]
    return(DF)
}

# Remove genes with TPM < 1
fem_counts <- median_filter(DF = fem_counts, TPM = 1)
male_counts <- median_filter(DF = male_counts, TPM = 1)

# How many genes are left after filtering in each sex?
nrow(fem_counts) # 15,604
nrow(male_counts) # 15,568

# Make a gene expression matrix that is the union of m/f genes remaining
# Use zoo to combine dfs with different num rows and different cols
sex_filtered_mat <- data.frame(fem_counts, cbind(zoo(, 1:nrow(fem_counts)),
                                                 as.zoo(male_counts)))

# Check for the expected number og rows/cols
nrow(sex_filtered_mat) == nrow(fem_counts) # TRUE 
ncol(sex_filtered_mat) == ncol(fem_counts) + ncol(male_counts) # TRUE

#_______________________________________________________________________________
# Plot histogram of the median log(TPM) across all regions 
#_______________________________________________________________________________
# Plot
pdf(SEX_PLT)
hist(log10(apply(sex_filtered_mat, 1, median)), 100,
	main = "Histogram of gene expression with median TPM >= 1 filtered by sex", 
   	xlab = "Median log10(TPM)",
	ylab = "Frequency")
dev.off()

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
	res <- unfiltered_mat[, c(lst)]
	return(res)
}
counts_by_tissue <- lapply(tissue_lst, sort_counts)

# Apply different filtering thresholds
TPM_1 <- Map(median_filter, TPM = 1, counts_by_tissue)
TPM_5 <- Map(median_filter, TPM = 5, counts_by_tissue)
TPM_10 <- Map(median_filter, TPM = 10, counts_by_tissue)

# Name each df in list by  which tissue type it is
names(TPM_1) <- names(tissue_lst)
names(TPM_5) <- names(tissue_lst)
names(TPM_10) <- names(tissue_lst)

# Number of genes left after filtering by region with various thresholds
res1 <- lapply(TPM_1, nrow)
res2 <- lapply(TPM_5, nrow)
res3 <- lapply(TPM_10, nrow)

# Write summary to file
res <- do.call(rbind, Map(data.frame, thresh_1=res1, thresh_5=res2, thresh_10=res3))
write.csv(res, file = TABLE)

#_______________________________________________________________________________
# Plot histogram of the median log(TPM) in each region separately 
#_______________________________________________________________________________
# Function to plot histograms on one page
plot_func <- function(mat, tissue, thresh) {
	hist(log10(apply(mat, 1, median)), 100, main = tissue)
    mtext(paste("Histograms of gene expression with median TPM >=", thresh, "filtered by region"), 
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

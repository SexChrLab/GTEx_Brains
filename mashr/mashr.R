#!/usr/bin/env Rscript

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/expression_matrices/output/filtered_counts_with_sex_chr.csv")

# load libraries
library(mashr)
library(readr)
library(data.table)

# Read in files
meta <- read.csv(METADATA, header = TRUE, stringsAsFactors = FALSE)
counts <- read.csv(COUNTS, header = TRUE, stringsAsFactors = FALSE)

# _____________________________________________________________________________
# Convert test results (csv files) into effect and standard error matrices
# _____________________________________________________________________________  
# Path to limma results; female gene expression regressed onto male
p1 <- "/scratch/mjpete11/human_monkey_brain/limma/tables"

# Read in csv files in directory as list of data frames
limma_lst <- Map(fread, header=TRUE, list.files(path = p1, pattern = "*.csv", full.names=T))

# Assign names to dfs in list
names(limma_lst) <- levels(as.factor(meta$Tissue))

# Add gene name and id column
#limma_lst <- lapply(limma_lst, function(x) cbind(counts[, 1], x))

# Combine the first column on every df in list into one matrix 
# Name columns in each df, then choose which to drop 
name_drop <- function(x, d){
    x <- data.frame(x)
    colnames(x) <- c("gene_ID", "beta", "st_dev")
    drp <- c(d)
    x <- x[,!colnames(x) %in% drp]
    return(x)
}    
effect_lst <- Map(name_drop, x=limma_lst, d="st_dev")
st_dev_lst <- Map(name_drop, x=limma_lst, d="beta")

# Combine the betas into one matrix
effect_df <- Reduce(function(...) merge(..., by = "gene_ID", all = TRUE), effect_lst)
colnames(effect_df) <- c("gene_ID", levels(as.factor(meta$Tissue)))
effect_mat <- as.matrix(effect_df[, 2:ncol(effect_df)])

# The gene order will be different in the output
# Check a random gene and make sure the corresponding values are correct
#lapply(limma <- lst, function(x) x[which(x$V1=="ENSG00000279457.4"),2])

# Cobine the sqrt(SEs) into one matrix
st_dev_df <- Reduce(function(...) merge(..., by = "gene_ID", all = TRUE), st_dev_lst)
colnames(st_dev_df) <- c("gene_ID", levels(as.factor(meta$Tissue)))
st_dev_mat <- as.matrix(st_dev_df[, 2:ncol(st_dev_df)])

# Add rownames to the effect (Bhat) and st dev (Shat) matrices so they are included in the output
rownames(effect_mat) <- effect_df$gene_ID
rownames(st_dev_mat) <- effect_df$gene_ID

#________________________________________________________________________________________________________  
# sanity check 
#________________________________________________________________________________________________________  
sanity <- function(e, s){
    # Do effect and st_dev_matrices have same nrow and rownames in same order?
    print(nrow(s) == nrow(e)) # TRUE
    print(identical(rownames(s), rownames(e))) # TRUE

    # Check for missing values
    # TRUE= All cols contain finite entries
    print(all(sapply(e, function(x) all(is.finite(x))))) 
    print(all(sapply(s, function(x) all(is.finite(x))))) 

    # FALSE = No cols contain missing values
    print(any(sapply(e, function(x) any(is.nan(x))))) 
    print(any(sapply(s, function(x) any(is.nan(x))))) 
   
    }
# Expect: 4T, 2F
sanity(e=effect_mat, s=st_dev_mat)

# ________________________________________________________________________________________________________  
# Mashr analysis; NSM lab approach
# ________________________________________________________________________________________________________  
# Create the mashr data object
mash_limma <- mash_set_data(Bhat=effect_mat, Shat=st_dev_mat)

# Apply multivariate adaptive shrinkage method in initial mode ("naive" run)
# Compute canonical covariance matrix
U_c <- cov_canonical(mash_limma)
m_c <- mash(mash_limma, U_c) 

# Save initial results and pairwise sharing matrix
# New vignette uses get fsr() (false sign rate); investigate later
strong_limma <- get_significant_results(m_c, thresh=0.05, sig_fn=ashr::get_lfdr)
shared_limma <- get_pairwise_sharing(m_c, factor=0.05)

# Print IDs of strong genes to standard output
print(strong_limma)

# Get a randomized set of results with the same length as the significant results
#random_limma <- sample(1:nrow(mash_limma$Bhat), length(strong_limma))

# Increase size of random subset to 50% of input data because we do not expect
# most genes to be significant
random_limma <- sample(1:nrow(mash_limma$Bhat), round(nrow(effect_mat)/2))

# Estimate null correlation on random subset of the data
# The resulting correlation matrix (vhat_limma) will be used to construct the "random" and "strong"mashr datasets
temp_limma <- mash_set_data(mash_limma$Bhat[random_limma,], mash_limma$Shat[random_limma,])
temp_U_c_limma <- cov_canonical(temp_limma)
vhat_limma <- estimate_null_correlation(temp_limma, temp_U_c_limma)

# This step was failing...works with new values for missing values
estimate_null_correlation(mash_limma, U_c)

# Create mashr datasets on the sig results ("strong" set) and random results ("random" set)
mash_limma_random <- mash_set_data(mash_limma$Bhat[random_limma,],mash_limma$Shat[random_limma,],V=vhat_limma)
mash_limma_strong <- mash_set_data(mash_limma$Bhat[strong_limma,],mash_limma$Shat[strong_limma,],V=vhat_limma)

# Perform PCA and extreme deconvolution
U_pca_limma <- cov_pca(mash_limma_strong, 5)
U_ed_limma <- cov_ed(mash_limma_strong, U_pca_limma)

# Compute canonical covariance matrix on random subset
U_c_random <- cov_canonical(mash_limma_random)

# Perform mashr on random subset for the purpose of estimating g
# In contrast to the earlier "naive" model, empirical covariance matrices from the extreme deconvolution and canonical methods are used
m_limma <- mash(mash_limma_random, Ulist=c(U_ed_limma, U_c_random), outputlevel=1)

# Run mashr on the strong subset, this time substituting the empirical value of g
#m2_limma <- mash(mash_limma_strong, g=get_fitted_g(m_limma), fixg=TRUE)

# Use covar matrix learned using strong effects, but run on the entire mash data set.
m2_limma <- mash(mash_limma, g=get_fitted_g(m_limma), fixg=TRUE)

# Calculate pairwise sharing
limma_share <- get_pairwise_sharing(m2_limma, factor=0.5)

# Posterior local false sign rate
lfsr_limma <- get_lfsr(m2_limma)

# Set a flase sign rate cutoff
fsr_cutoff <- 0.2

# Get number of significant genes in each region
sig_genes <- apply(lfsr_limma, 2, function(x) sum(x < fsr_cutoff))
sig_genes <- data.frame(region = names(sig_genes), num_sig_genes = sig_genes) 

# Get posterior values
# Posterior betas
beta_limma <- as.data.frame(get_pm(m2_limma))

# Write to file
saveRDS(m2_limma, file = "/scratch/mjpete11/human_monkey_brain/mashr/output/mashr_results.rds")
write.csv(limma_share, file = "/scratch/mjpete11/human_monkey_brain/mashr/output/limma_share.csv")
write.csv(sig_genes, file = "/scratch/mjpete11/human_monkey_brain/mashr/output/sig_genes.csv", row.names = FALSE)
write.csv(beta_limma, file = "/scratch/mjpete11/human_monkey_brain/mashr/output/beta_limma.csv")
write.csv(lfsr_limma, file = "/scratch/mjpete11/human_monkey_brain/mashr/output/lfsr_limma.csv")

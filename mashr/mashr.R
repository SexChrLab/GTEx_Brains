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
#names(limma_lst) <- c("Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex",
#                      "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", 
#                      "Putamen", "Spinal_Cord", "Substantia_Nigra")
names(limma_lst) <- levels(as.factor(meta$Tissues))
names(limma_lst)

# Add gene name and id column
limma_lst <- lapply(limma_lst, function(x) cbind(counts[, 1], x))

# Combine the first column on every df in list into one matrix 
# Name columns in each df, then choose which to drop 
name_drop <- function(x, d){
    x <- data.frame(x)
    colnames(x) <- c("V1", "beta", "sqrt_SE")
    drp <- c(d)
    x <- x[,!colnames(x) %in% drp]
    return(x)
}    
effect_lst <- Map(name_drop, x=limma_lst, d="sqrt_SE")
sqrt_s2_lst <- Map(name_drop, x=limma_lst, d="beta")

# Iteratively join dfs in list by gene_ID col
# supress warning about duplicate cols when merging; keep all rows, including those not present in both
# extremely hacky
lst_to_mat <- function(l, m){
    # supress warning about duplicate cols when merging; keep all rows, including those not present in both
    x <- suppressWarnings(Reduce(function(df1, df2) merge(df1, df2, by="V1", all=TRUE), l))
    rownames(x) <- x[["V1"]] # join by gene_ID col, then drop
    x[["V1"]] <- NULL
#    colnames(x)[1:ncol(x)] <- levels(as.factor(meta$Tissue)) # explicitly name conditions
    x[is.na(x)] <- m # Replace NA (row present in a df but not in all others) with value of choice
    x <- as.matrix(x) # mashr only accepts matrices
    print(head(x)) # visually inspect df
    print(tail(x))
    return(x)
}
effect_mat <- lst_to_mat(l=effect_lst, m=0)
sqrt_s2_mat <- lst_to_mat(l=sqrt_s2_lst, m=100)

# temporary hack
effect_mat <- apply(effect_mat, 2, as.numeric)
head(effect_mat)

# ________________________________________________________________________________________________________  
# sanity check 
# ________________________________________________________________________________________________________  
sanity <- function(e, s){
    # Do effect and sqrt_s2_matrices have same nrow and rownames in same order?
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
sanity(e=effect_mat, s=sqrt_s2_mat)

# ________________________________________________________________________________________________________  
# Mashr analysis; NSM lab approach
# ________________________________________________________________________________________________________  
# Create the mashr data object
mash.Limma <- mash_set_data(Bhat=effect_mat, Shat=sqrt_s2_mat)

# Apply multivariate adaptive shrinkage method in initial mode ("naive" run)
# Compute canonical covariance matrix
U.c <- cov_canonical(mash.Limma)
m.c <- mash(mash.Limma, U.c) # error: chol(): decompisition failed

# Save initial results and pairwise sharing matrix
strong.Limma <- get_significant_results(m.c, thresh=0.05, sig_fn=ashr::get_lfdr)
shared.Limma <- get_pairwise_sharing(m.c, factor=0.05)

# Get a randomized set of results with the same length as the significant results
random.Limma <- sample(1:nrow(mash.Limma$Bhat), length(strong.Limma))

# Estimate null correlation on random subset of the data
# The resulting correlation matrix (Vhat.Limma) will be used to construct the "random" and "strong"mashr datasets
temp.Limma <- mash_set_data(mash.Limma$Bhat[random.Limma,], mash.Limma$Shat[random.Limma,])
temp.U.c.Limma <- cov_canonical(temp.Limma)
Vhat.Limma <- estimate_null_correlation(temp.Limma, temp.U.c.Limma)

# This step was failing...works with new values for missing values
estimate_null_correlation(mash.Limma, U.c)

# Create mashr datasets on the sig results ("strong" set) and random results ("random" set)
mash.Limma.random <- mash_set_data(mash.Limma$Bhat[random.Limma,],mash.Limma$Shat[random.Limma,],V=Vhat.Limma)
mash.Limma.strong <- mash_set_data(mash.Limma$Bhat[strong.Limma,],mash.Limma$Shat[strong.Limma,],V=Vhat.Limma)

# Perform PCA and extreme deconvolution
U.pca.Limma <- cov_pca(mash.Limma.strong, 5)
U.ed.Limma <- cov_ed(mash.Limma.strong, U.pca.Limma)

# Compute canonical covariance matrix on random subset
U.c.random <- cov_canonical(mash.Limma.random)

# Perform mashr on random subset for the purpose of estimating g
# In contrast to the earlier "naive" model, empirical covariance matrices from the extreme deconvolution and canonical methods are used
m.Limma <- mash(mash.Limma.random, Ulist=c(U.ed.Limma, U.c.random), outputlevel=1)

# Now run mashr on the strong subset, this time substituting the empirical value of g
m2.Limma <- mash(mash.Limma.strong, g=get_fitted_g(m.Limma), fixg=TRUE)

# Calculate pairwise sharing
Limma.share <- get_pairwise_sharing(m2.Limma, factor=0.5)

# Set a flase sign rate cutoff
fsr.cutoff <- 0.2

# Get number of significant genes in each region
apply(get_lfsr(m2.Limma), 2, function(x) sum(x < fsr.cutoff))

# Get posterior values
# Posterior betas
beta.Limma <- get_pm(m2.Limma)

# Posterior local false sign rate
lfsr.Limma <- get_lfsr(m2.Limma)

#!/usr/bin/env Rscript

# DGX with limma-voom

#Load packages                                                                 
library(limma)
library(edgeR) 
library(readr)
library(stringr)
library(rjson)
library(dplyr)

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/Metadata.csv"
COUNTS = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/Gene_Salmon_CountMatrix.tsv"

#---------------------------------------------------------------------------------------------------------------------
# Data Preprocessing and normalization 
#---------------------------------------------------------------------------------------------------------------------
# Read Metadata CSV.                                                            
Samples <- read.csv(METADATA, header = TRUE, stringsAsFactors=FALSE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample   

# Set tissues column as factor to add levels attribute
Samples$Tissue <- as.factor(Samples$Tissue)

# Read in count matrix
cts <- read.table(COUNTS, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Identify sample in metadata that does not have count data and drop from metadata
setdiff(Samples$Sample, colnames(cts)) # "GTEX-13N2G-0011-R2a-SM-5MR4Q"
Samples <- Samples[!grepl("GTEX-13N2G-0011-R2a-SM-5MR4Q", Samples$Sample),]

# Number of samples with count data after subsetting with metadata
ncol(cts) == nrow(Samples) # 1,340 samples total

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
      Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Function to split count matrix into list of dfs
Tissues <- names(Meta)

Split_Cols <- function(w, z){
  names <- which(colnames(w) %in% subset(Samples[['Sample']], Samples[['Tissue']] == z))
  cts_lst <- cts[, names]
  return(cts_lst)
}

# Split count matrix into list of dfs
cts_lst <- list()
for (i in Tissues){
      cts_lst[[i]] <- Split_Cols(w=cts, z=i)
}

#---------------------------------------------------------------------------------------------------------------------
# Sanity check; make sure samples are in same order in count dfs and metadata dfs
#---------------------------------------------------------------------------------------------------------------------
Check_Order <- function(x, z){
    cts_lst <- identical(colnames(x), rownames(z))
    return(cts_lst)
}
all(Map(Check_Order, x=cts_lst, z=Meta)) # TRUE

# Sort cols in cts in same order as rows in metadata anyways, just to be safe
# If the order of samples in metadata and count data don't match, 
# samples will be labelled incorrectly in the design matrix
Sort_Cols <- function(x, z){
    x <- x[, match(rownames(z), colnames(x))]
    return(x)
}
cts_lst <- Map(Sort_Cols, x=cts_lst, z=Meta)
all(Map(Check_Order, x=cts_lst, z=Meta)) # TRUE

#---------------------------------------------------------------------------------------------------------------------
# Make list of DGEList objects with corresponding design matrix 
#---------------------------------------------------------------------------------------------------------------------
# Create design matrix: Step 1
# Create sex factor
Factor_Func <- function(x){
  cts_lst <- factor(x[['Sex']])
  return(cts_lst)
}
Groups <- lapply(Meta, Factor_Func)

# Combine count dfs and factors to make list of DGEList objects; 
# required object class for limma
DGE_Func <- function(d, l){
    cts_lst <- DGEList(d, group=l)
    return(cts_lst)
}
DGE_Lst <- Map(DGE_Func, d=cts_lst, l=Groups)

# Create design matrix: Step 2; Part 1
Model_Func <- function(fc, d){
    x <- model.matrix(~0 + fc, data = d[['samples']])
    colnames(x) <- gsub("fc", "", colnames(x)) # adds var name to col name; drop it
    return(x)
}
Design <- Map(Model_Func, fc=Groups, d=DGE_Lst)

# Create design matrix: Step 2; Part 2
Set_Levels <- function(x, z){
    colnames(x) <- levels(z[['samples']][['group']])
    return(x)
}
Design <- Map(Set_Levels, x=Design, z=DGE_Lst)

#---------------------------------------------------------------------------------------------------------------------
# Filter genes by cpm expression level 
#---------------------------------------------------------------------------------------------------------------------
# Keep only genes with cpm > 1 in at least half the samples for each tissue type
# How many cols have cpm > 1 in each row; return TRUE if that value >= half the number of samples 
Keep <- lapply(DGE_Lst, function(x){rowSums(cpm(x[['counts']]) > 1) >= ceiling(ncol(x[['counts']])/2)})

# Do all genes pass filter?
all(lapply(Keep, function(x) all(x) == TRUE)) # FALSE

# Range of genes left across conditions
sapply(Keep, function(x) sum(x)) # 16,486-18,273 

# How many genes are left is I use a more stringent filter: cpm > 1 in all samples 
Keep.test <- lapply(DGE_Lst, function(x){rowSums(cpm(x[['counts']]) > 1) >= ceiling(ncol(x[['counts']]))})

# Range of genes left across conditions after applying more stringent filter
sapply(Keep.test, function(x) sum(x)) # 10-12,539 

# Apply filter (cpm >1 in at least half the samples) to list of DGEList object
Filter_Func <- function(x, k){
      x <- x[k, , keep.lib.sizes=FALSE]
}
DGE_Lst <- Map(Filter_Func, x=DGE_Lst, k=Keep)

# TMM Normalization
DGE_Lst <- lapply(DGE_Lst, calcNormFactors)

#---------------------------------------------------------------------------------------------------------------------
# Voom
#---------------------------------------------------------------------------------------------------------------------
# Apply voom (rather than plain limma) since library sizes are highly variable
Voom_Func <- function(x, d){
    v <- voom(x, d, plot=TRUE)
    fit <- lmFit(v, d)
    fit <- eBayes(fit)
    return(fit)
}
Voom_Res <- Map(Voom_Func, x=DGE_Lst, d=Design)

# Combine the betas (second col in coeff matrix) with the standard error (var^2)
Make_Tables <- function(x){
    cts_lst <- data.frame(Beta=x[['coefficients']][,2], SE=x[['s2.post']])
    return(cts_lst)
}
Table_Lst <- lapply(Voom_Res, Make_Tables)

# Write tables 
sapply(names(Table_Lst), function(x) write.table(Table_Lst[[x]], file=paste(x, "csv", sep="."))) 

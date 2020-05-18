# DGX with limma-voom

#Load packages                                                                 
library(limma)
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(rjson)
library(dplyr)

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/Metadata.csv"
COUNTS = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/Gene_Salmon_CountMatrix.tsv"

#---------------------------------------------------------------------------------------------------------------------
# Data Preprocessing, experimental design, and normalization 
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

# Remove sample counts that do not have metadata
# cts <- cts[intersect(colnames(cts), as.character(Samples$Sample))]

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
  res <- cts[, names]
  return(res)
}

# Split count matrix into list of dfs
res <- list()
for (i in Tissues){
      res[[i]] <- Split_Cols(w=cts, z=i)
}

# Create design matrix: Step 1
# Create sex and tissue factor
Factor_Func <- function(x){
  res <- factor(paste(x$Tissue, x$Sex, sep="."))
  return(res)
}
Groups <- lapply(Meta, Factor_Func)

# Add column with new factor
Col_Bind <- function(df, fc){
      cbind(df, fc)
}
Meta <- Map(Col_Bind, Meta, Groups)

# Check that the count samples are in the same order as the metadata 
Check_Order <- function(x, z){
    res <- identical(colnames(x), rownames(z))
    return(res)
}
all(Map(Check_Order, x=res, z=Meta)) # TRUE

# Sort cols in cts in same order as rows in Meta
Sort_Cols <- function(x, z){
  x <- x[, match(rownames(z), colnames(x))]
  return(x)
}
cts <- Map(Sort_Cols, x=res, z=Meta)

# Create list of DGEList objects
DGE_Func <- function(df, lst){
      res <- DGEList(df, group=lst)
  return(res)
}
DGE_Lst <- Map(DGE_Func, cts, Groups)

# Create design matrix: Step 2; Part 1
Model_Func <- function(fc, df){
      res <- model.matrix(~0 + fc, data = df$samples)
  return(res)
}
Design <- Map(Model_Func, Groups, DGE_Lst)

# Remove "fc" from colnames
Rename_Cols <- function(x){
      colnames(x) <- gsub("fc", "", colnames(x))
  return(x)
}
Design <- lapply(Design, Rename_Cols)

# Create design matrix: Step 2; Part 2
Set_Levels <- function(x, z){
      colnames(x) <- levels(z$samples$group)
  return(x)
}
Design <- Map(Set_Levels, x=Design, z=DGE_Lst)

# Keep only genes expressed in at least half the samples for each tissue type
# How many cols have cpm > 1 in each row; return TRUE if that value >= half the number of samples 
Keep <- lapply(DGE_Lst, function(x){rowSums(cpm(x[['counts']]) > 1) >= ceiling(ncol(x[['counts']])/2)})

# Every single gene passed the filter...
# Are all values in each list TRUE?
all(lapply(Keep, function(x) all(x) == TRUE)) # TRUE

# Do all genes pass the filter if they mucst be expressed in all samples?
Keep.test <- lapply(DGE_Lst, function(x){rowSums(cpm(x[['counts']]) > 1) >= ceiling(ncol(x[['counts']]))})

# Are all values in each list TRUE?
all(lapply(Keep.test, function(x) all(x) == TRUE)) # FALSE

# Apply expression threshold filter
Filter_Func <- function(x, k){
      x <- x[k, , keep.lib.sizes=FALSE]
}
DGE_Lst <- Map(Filter_Func, x=DGE_Lst, k=Keep)

# TMM Normalization
DGE_Lst <- lapply(DGE_Lst, calcNormFactors)

#---------------------------------------------------------------------------------------------------------------------
# Voom
#---------------------------------------------------------------------------------------------------------------------
# Chose voom since library sizes are > 3 fold diff between largest and smallest
v <- voom(dge, design, plot=TRUE)

# Option to apply between-array quantile normalization for very noisy data
#v <- voom(counts, design, plot=TRUE, normalize="quantile")

# Fit model
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef=ncol(design))

# Option to give more weight to the log-fold changes in the rankings
#fit <- treat(fit, lfc=log2(1.2))
#topTreat(fit, coef=ncol(design))

#---------------------------------------------------------------------------------------------------------------------
# Write results to disk 
#---------------------------------------------------------------------------------------------------------------------




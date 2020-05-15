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
Samples = read.csv(METADATA, header = TRUE, stringsAsFactors=FALSE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample   

# Set tissues column as factor to add levels attribute
Samples$Tissue <- as.factor(Samples$Tissue)

# Drop samples < 55 
Samples <- subset(Samples, Samples[,4] > 55)

# Total number of samples after filtering by age
nrow(Samples) # 932

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
      Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Read in count matrix
cts <- read.table(COUNTS, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Remove Samples not present in metadata
cts <- cts[intersect(colnames(cts), as.character(Samples$Sample))]

# Function to split count matrix into list of dfs
Tissues <- names(Meta)

Split_Cols <- function(w, z){
  names <- which(colnames(w) %in% subset(Samples[['Sample']], Samples[['Tissue']] == z))
  res <- cts[, names]
  return(res)
}
#Split_Cols <- function(w, z){
#  names <- which(colnames(w) %in% Samples[['Sample']] & Samples[['Tissue']] == z)
#  res <- cts[, names]
#  return(res)
#}

# Split count matrix into list of dfs
res <- list()
for (i in Tissues){
      res[[i]] <- Split_Cols(w=cts, z=i)
}

### There is a sample in Meta[[13]] that is not present in res[[13]]!!!
# Check if columns were subset correctly
check <- list()
for (i in 1:13){
      check[[i]] <- colnames(res[[i]]) %in% subset(Samples$Sample, Samples$Tissue==Tissues[[i]])
}
all(as.logical(lapply(check, all)))

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
Keep <- lapply(y, function(x){
                     rowSums(cpm(x[['counts']])>1) >= ncol(x[['counts']]) 
})

Filter_Func <- function(x, k){
      x <- x[k, , keep.lib.sizes=FALSE]
}
DGE_Lst <- Map(Filter_Func, DGE_Lst, Keep)

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




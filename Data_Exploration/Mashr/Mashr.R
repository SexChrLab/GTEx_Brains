#!/usr/bin/env Rscript

setwd("/scratch/mjpete11/GTEx/Data_Exploration/Mashr/")

library(mashr)
library(readr)
library(data.table)

# ________________________________________________________________________________________________________  
# Function to convert test results (csv files) into effect matrices
# ________________________________________________________________________________________________________  
# Paths to csv files with significant differentially expressed genes
# Add transcripts later
#BASE  <- "/scratch/mjpete11/GTEx/Differential_Expression/"
#EXACT  <- "Exact_Test/"
#FTEST <- "F_Test/"
#RATIO <- "Ratio_Test/"
#SAL  <- "Salmon/"
#HIS  <- "Hisat/"
#AGE <- "Age_Matched/Gene/"
#MATCH <- "Matched/Gene/"
#
## Function to generate paths
#Paths <- function(test, aligner, sample_set){
#    res <- list()
#    for (i in test){
#        for (j in aligner){
#            for (k in sample_set){
#                path <- paste0(BASE, i, j, k)
#                res <- c(res, path)
#            }
#        }
#    }
#    return(res)
#}
#p <- Paths(test=c(EXACT, FTEST, RATIO), aligner=c(SAL, HIS), sample_set=c(AGE, MATCH))
#
## Function to read csv files as list of dfs
##Read_Tables <- function(paths){
##    res <- list()
#    for (p in paths){
#         lst <- Map(fread, header=FALSE, list.files(path = p, pattern = "*.csv", full.names=T))
#         res <- c(res, lst)
#    }
#    return(res)
#}
#Tables.Lst <- Read_Tables(paths=p) 

# Read in csv files as list of dfs;,pass list of file names with list.files
# Started writing function, but I need it to return a nested named list; a list of lists is too dangerous lol
# Exact test, salmon 
#Names <- c(ESA, ESM, EHA, EHM, FSA, FSM, FHA, FHM, RSA, RSM, RHA, RHM)
#
#ESA.lst <- Map(fread, header=FALSE, list.files(path = p1, pattern = "*.csv", full.names=T))
#ESM.lst <- Map(fread, header=FALSE, list.files(path = p2, pattern = "*.csv", full.names=T))
#EHA.lst <- Map(fread, header=FALSE, list.files(path = p3, pattern = "*.csv", full.names=T))
#EHM.lst <- Map(fread, header=FALSE, list.files(path = p4, pattern = "*.csv", full.names=T))
#FSA.lst <- Map(fread, header=FALSE, list.files(path = p5, pattern = "*.csv", full.names=T))
#FSM.lst <- Map(fread, header=FALSE, list.files(path = p6, pattern = "*.csv", full.names=T))
#FHA.lst <- Map(fread, header=FALSE, list.files(path = p7, pattern = "*.csv", full.names=T))
#FHM.lst <- Map(fread, header=FALSE, list.files(path = p8, pattern = "*.csv", full.names=T))
#RSA.lst <- Map(fread, header=FALSE, list.files(path = p9, pattern = "*.csv", full.names=T))
#RSM.lst <- Map(fread, header=FALSE, list.files(path = p10, pattern = "*.csv", full.names=T))
#RHA.lst <- Map(fread, header=FALSE, list.files(path = p11, pattern = "*.csv", full.names=T))
#RHM.lst <- Map(fread, header=FALSE, list.files(path = p12, pattern = "*.csv", full.names=T))
#
## Assign names to dfs in list
#p1 <- "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene/"
#ESA.lst <- Map(fread, header=FALSE, list.files(path = p1, pattern = "*.csv", full.names=T))
#names(ESA.lst) <- c("Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex",
#                 "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", 
#                 "Putamen", "Spinal_Cord", "Substantia_Nigra")
#
## Name items in list, name the df columns then drop last two cols by name
#Name_Drop <- function(x){
#    x <- data.frame(x)
#    colnames(x) <- c("gene_ID", "logFC", "logCPM", "PValue")
#    drp <- c("logFC", "PValue")
#    x <- x[,!colnames(x) %in% drp]
#    return(x)
#}    
#ESA.lst<- lapply(ESA.lst, Name_Drop)
#
## Left join list of dfs by gene_ID col; if row (gene ID) not present in a df, fill with zero
## After joining by gene_ID col, set it as rowname and drop col (this is how the example df looks) 
#ESA.df <- Reduce(function(df1, df2) merge(df1, df2, by="gene_ID", all=TRUE), ESA.lst)
#rownames(ESA.df) <- ESA.df[["gene_ID"]]
#ESA.df[["gene_ID"]] <- NULL
#colnames(ESA.df)[1:ncol(ESA.df)] <- names(ESA.lst)
#ESA.df[is.na(ESA.df)] <- 0
#
#print(head(ESA.df))
#
## Check for missing or non-finite values
#any(is.nan(ESA.df$logFC)) # No NANs
#all(is.finite(ESA.df$logFC)) # All values are finite
#
## Convert df to matrix (mash_set_data() only accepts matrices) 
#ESA.m <- as.matrix(ESA.df)



# ________________________________________________________________________________________________________  
# Convert test results (csv files) into effect matrices
# ________________________________________________________________________________________________________  
# Read in tables with logFC and make the effect matrix
# Cols: gene ID followed by tissue type, Rows: gene, Values: logFC
# Paths to csv files
# Salmon, Exact
p1 <- "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene"

# Read in csv files as list of dfs (they are actually tab delim...);,pass list of file names with list.files
ESA.lst <- Map(fread, header=FALSE, list.files(path = p1, pattern = "*.csv", full.names=T))

# Assign names to dfs in list
names(ESA.lst) <- c("Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex",
                 "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", 
                 "Putamen", "Spinal_Cord", "Substantia_Nigra")

# Name items in list, name the df columns then drop last two cols by name
Name_Drop <- function(x){
    x <- data.frame(x)
    colnames(x) <- c("gene_ID", "logFC", "logCPM", "PValue")
    drp <- c("logFC", "PValue")
    x <- x[,!colnames(x) %in% drp]
    return(x)
}    
ESA.lst<- lapply(ESA.lst, Name_Drop)

# Left join list of dfs by gene_ID col; if row (gene ID) not present in a df, fill with zero
# After joining by gene_ID col, set it as rowname and drop col (this is how the example df looks) 
ESA.df <- Reduce(function(df1, df2) merge(df1, df2, by="gene_ID", all=TRUE), ESA.lst)
rownames(ESA.df) <- ESA.df[["gene_ID"]]
ESA.df[["gene_ID"]] <- NULL
colnames(ESA.df)[1:ncol(ESA.df)] <- names(ESA.lst)
ESA.df[is.na(ESA.df)] <- 0

print(head(ESA.df))

# Check for missing or non-finite values
any(is.nan(ESA.df$logFC)) # No NANs
all(is.finite(ESA.df$logFC)) # All values are finite

# Convert df to matrix (mash_set_data() only accepts matrices) 
ESA.m <- as.matrix(ESA.df)

# ________________________________________________________________________________________________________  
# Mashr analysis 
# ________________________________________________________________________________________________________  
# Step 1: Select strong signals
# Create mashr object using option to automatically generate the standard error matrix (fill with 1s)
mash.ESA <- mash_set_data(Bhat=ESA.m, Shat=1)

# Run a condition-by-condition analysis on all the data
m.1by1 <- mash_1by1(mash.ESA)

# Selected strong signals; i.e. indices where the lfsr < 0.05
strong <- get_significant_results(m.1by1, 0.05)

# Step 2: Obtain initial data-driven covariance matrices
U.pca <- cov_pca(mash.ESA, 5, subset=strong)

# Step 3: Extreme deconvolution of principle components with strong signals only
U.ed <- cov_ed(mash.ESA, U.pca, subset=strong)

# Step 4: Run mash
# Crucial! Must fit mash to all tests, not just the strong subset
m.ed <- mash(mash.ESA, U.ed)
get_loglik(m.ed) # -433173.4

# Compare to results using cononical covariance matrix
U.c <- cov_canonical(mash.ESA)
m.c <- mash(mash.ESA, U.c)
get_loglik(m.c) #-404642.2

# Extract posterior summaries
Tissues <- c("Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex",
             "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", 
             "Putamen", "Spinal_Cord", "Substantia_Nigra")

# Barplot of estimated mixture proportions of the data-driven covariance and canonical covar
tiff("data_driven_barplot.tiff")
par(mar=c(10,3,2,1)) # bottom, left, top, right
barplot(get_estimated_pi(m.ed), las=2)
dev.off()

tiff("canonical_barplot.tiff")
par(mar=c(10,3,2,1)) # bottom, left, top, right
barplot(get_estimated_pi(m.c), las=2)
dev.off()

# Metaplot based on posterior means and posterior effects for both covar matrices
tiff("data_driven_metaplot.tiff")
mash_plot_meta(m.ed, get_significant_results(m.ed)[1], labels=Tissues)
dev.off()

tiff("canonical_metaplot.tiff")
mash_plot_meta(m.c, get_significant_results(m.c)[1], labels=Tissues)
dev.off()






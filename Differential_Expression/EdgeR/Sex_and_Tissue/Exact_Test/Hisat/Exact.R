# This script uses an exact test to tests for differential gene expression between males and females 
# within each brain tissue type using gene level counts on the age matched samples.
setwd("/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat")

# Constants
# Hisat/stringtie results are stored in seperate matrices because the same transcripts/genes reported are tissue-specific
# Using annotated gene count matrix
METADATA <- snakemake@input[[1]]
PATHS.1 <- snakemake@input[[2]]
PATHS.2 <- snakemake@input[[3]]
PATHS.3 <- snakemake@input[[4]]
PATHS.4 <- snakemake@input[[5]]
PATHS.5 <- snakemake@input[[6]]
PATHS.6 <- snakemake@input[[7]]
PATHS.7 <- snakemake@input[[8]]
PATHS.8 <- snakemake@input[[9]]
PATHS.9 <- snakemake@input[[10]]
PATHS.10 <- snakemake@input[[11]]
PATHS.11 <- snakemake@input[[12]]
PATHS.12 <- snakemake@input[[13]]
PATHS.13 <- snakemake@input[[14]]
# Combine to one vector
PATHS <- c(PATHS.1,PATHS.2,PATHS.3,PATHS.4,PATHS.5,PATHS.6,PATHS.7,PATHS.8,PATHS.9,PATHS.10,PATHS.11,PATHS.12,PATHS.13)

# Output
MD_PLOT <- snakemake@output[[1]]
VOLCANO_PLOT <- snakemake@output[[2]]
UP_JSON <- snakemake@output[[3]]
DOWN_JSON <- snakemake@output[[4]]

# Load packages                                                                 
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(rjson)
library(dplyr)
library(org.Hs.eg.db)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample   

# Count matrices
Tissues <- c('Amygdala', 'Anterior', 'Caudate', 'Cerbellar', 'Cerebellum', 'Cortex', 'Frontal_Cortex',
'Hippocampus', 'Hypothalamus', 'Nucleus_Accumbens', 'Putamen', 'Spinal Cord', 'Substantia_Nigra')

cts <- lapply(PATHS, function(x){
  t <- read.table(file=x, sep="\t")
})
names(cts) <- Tissues

# Replace . to - in colnames in each df
for (i in seq_along(cts)){
  colnames(cts[[i]]) <- str_replace_all(colnames(cts[[i]]), pattern = "\\.","-")
}

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(samples$Tissue))){
  Meta[[i]] <- samples[samples$Tissue == levels(samples$Tissue)[i],]
}
names(Meta) <- levels(samples$Tissue)

# Remove samples not present in metadata
Subset_Func <- function(df.1, df.2) {
  df.1 <- df.1[, intersect(colnames(df.1), as.character(df.2$Sample))]
  return(df.1)
}
cts <-  Map(Subset_Func, cts, Meta)

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
cts <- Map(Sort_Cols, x=cts, z=Meta)

# Create list of DGEList objects
DGE_Func <- function(df, lst){
  res <- DGEList(df, group=lst)
  return(res)
}
y <- Map(DGE_Func, cts, Groups)

# Create design matrix: Step 2; Part 1
Model_Func <- function(fc, df){
  res <- model.matrix(~0 + fc, data = df$samples)
  return(res)
}
Design <- Map(Model_Func, Groups, y)

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
Design <- Map(Set_Levels, x=Design, z=y)

# Keep only genes expressed in at least half the samples for each tissue type
Keep <- lapply(y, function(x){
  rowSums(cpm(x[['counts']])>1) >= ncol(x[['counts']]) 
})

Filter_Func <- function(x, k){
  x <- x[k, , keep.lib.sizes=FALSE]
}
y <- Map(Filter_Func, y, Keep)

# TMM Normalization
y <- lapply(y, calcNormFactors)

# Estimate common dispersion and tagwise dispersions in one run (recommended)
Dispersion_Func <- function(a, b){
  estimateDisp(a, b, robust=TRUE)
}
y <- Map(Dispersion_Func, y, Design)

#---------------------------------------------------------------------------------------------------------------------
# Test for DGX with Exact Test
#---------------------------------------------------------------------------------------------------------------------
# Comaprisons to test; set males as baseline
Pairs <- list(c("Amygdala.Female", "Amygdala.Male"), 
              c("Anterior.Female", "Anterior.Male"), 
              c("Caudate.Female", "Caudate.Male"),
              c("Cerebellar.Female", "Cerebellar.Male"),
              c("Cerebellum.Female", "Cerebellum.Male"),
              c("Cortex.Female", "Cortex.Male"),
              c("Frontal_Cortex.Female", "Frontal_Cortex.Male"),
              c("Hippocampus.Female", "Hippocampus.Male"), 
              c("Hypothalamus.Female", "Hypothalamus.Male"), 
              c("Nucleus_Accumbens.Female", "Nucleus_Accumbens.Male"), 
              c("Putamen.Female", "Putamen.Male"),  
              c("Spinal_Cord.Female", "Spinal_Cord.Male"), 
              c("Substantia_Nigra.Female", "Substantia_Nigra.Male"))

# Apply exact test
Exact_Func <- function(x, comp){
  exactTest(x, comp)
}
Exact_Res <- Map(Exact_Func, y, Pairs)

#---------------------------------------------------------------------------------------------------------------------
# Summarize results 
#---------------------------------------------------------------------------------------------------------------------
# Function to correct for multiple testing
Test_Correct <- function(x){
  x[['table']][['PValue']] <- p.adjust(x[['table']][['PValue']],method="BH")
  return(x)
}
Corrected_Exact <- lapply(Exact_Res, Test_Correct)

# Function to filter p-vals, and filter by logFC
Up_Reg <- function(x){
  res <- x[['table']][x[['table']][['PValue']] < 0.05, ] 
  res <- res[res[['logFC']] > 0, ]
  return(res)
}

Down_Reg <- function(x){
  res <- x[['table']][x[['table']][['PValue']] < 0.05, ] 
  res <- res[res[['logFC']] < 0, ]
  return(res)
}
Up_Top <- lapply(Corrected_Exact, Up_Reg)
Down_Top <- lapply(Corrected_Exact, Down_Reg)

# Make table of up and down regulated genes for each tissue
Get_Vec <- function(x){
  res <- rownames(x)
  return(res)
}
Up_Genes <- lapply(Up_Top, Get_Vec)
Down_Genes <- lapply(Down_Top, Get_Vec)

# Write to file
Up_Json <- toJSON(Up_Genes)
Down_Json <- toJSON(Down_Genes)

 write(Up_Json, UP_JSON)
 write(Down_Json, DOWN_JSON)

#---------------------------------------------------------------------------------------------------------------------
# Mean-Difference Plots
#---------------------------------------------------------------------------------------------------------------------
# Plot Mean-Difference  plots on one page
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Hisat: Gene Mean-Difference Plots; Exact Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}

#Write to file
pdf(MD_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE)  # margins: c(bottom, left, top, right)
Res_Plots <- Map(MD_Plot_Func, x=Exact_Res, w=Tissues)
legend(20.0,0.0, legend=c("Up","Not Sig", "Down"), pch = 16, col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Volcano Plots
#---------------------------------------------------------------------------------------------------------------------
# Make df of neg log p-vals and logFC results after correcting for multiple testing 
Volcano_Func <- function(x){
  cbind(x[["logFC"]], -log10(x[["PValue"]]))
}
Volcano_Up <- lapply(Up_Top, Volcano_Func)
Volcano_Down <- lapply(Down_Top, Volcano_Func) 

# Function to make df of log p-vals and logFC on untransformed exact test results
Untrans_Volcano <- function(x){
      cbind(x[["table"]][["logFC"]], -log10(x[["table"]][,"PValue"]))
}
Volcano_Res <- lapply(Exact_Res, Untrans_Volcano)

# Coerce to df from mtx
Volcano_Up <- lapply(Volcano_Up, as.data.frame)
Volcano_Down <- lapply(Volcano_Down, as.data.frame)
Volcano_Res <- lapply(Volcano_Res, as.data.frame)

# Rename columns
colnames <- c("logFC", "negLogPval")

Rename_Cols_Func <- function(x){
  setNames(x, colnames)
}
Volcano_Up <- lapply(Volcano_Up, Rename_Cols_Func)
Volcano_Down <- lapply(Volcano_Down, Rename_Cols_Func)
Volcano_Res <- lapply(Volcano_Res, Rename_Cols_Func)

# Set ylim and xlim
xmax <- ceiling(max(as.numeric(lapply(Volcano_Res, function(x) max(x[['logFC']])))))
xmin <- ceiling(min(as.numeric(lapply(Volcano_Res, function(x) min(x[['logFC']])))))
ymax <- ceiling(max(as.numeric(lapply(Volcano_Res, function(x) max(x[['negLogPval']])))))
ymin <- ceiling(min(as.numeric(lapply(Volcano_Res, function(x) min(x[['negLogPval']])))))

# Plot
Plot_Func <- function(RES, TISSUE, UP, DOWN){
  plot(RES, pch=19, main=TISSUE, xlab = '', ylab = '', las = 1, ylim=c(ymin, ymax), xlim=c(xmin,xmax))
  with(UP, points(logFC, negLogPval, pch=19, col="green"))
  with(DOWN, points(logFC, negLogPval, pch=19, col="blue"))
  abline(a=-log10(0.05), b=0, col="blue") 
  abline(v=c(2,-2), col="red")
  mtext('Hisat: Gene Volcano Plots; Exact Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, RES=Volcano_Res, TISSUE=Tissues, UP=Volcano_Up, DOWN=Volcano_Down)
legend(10.0, 8.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()


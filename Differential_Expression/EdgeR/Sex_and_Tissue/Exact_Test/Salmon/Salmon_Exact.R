# This script looks at differential gene expression between males and females within each brain tissue type.
setwd("/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon")

METADATA <- "/scratch/mjpete11/GTEx/Metadata/Matched_Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Count_Matrix.tsv"
MD_PLOT <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Plots/Matched_Salmon_Exact_MD.pdf'
VOLCANO_PLOT <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Plots/Matched_Salmon_Exact_Volcano.pdf'

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(grid)

# Read Metadata CSV.                                                            
Samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample   

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Read in count matrix
cts <- read.table(COUNT_MATRIX, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Remove Samples not present in metadata
cts <- cts[intersect(colnames(cts), as.character(Samples$Sample))]

# Function to split count matrix into list of dfs
Tissues <- names(Meta)

Split_Cols <- function(w, z){
  names <- which(colnames(w) %in% Samples[['Sample']] & Samples[['Tissue']] == z)
  res <- cts[, names]
  return(res)
}

# Split count matrix into list of dfs
res <- list()
for (i in Tissues){
  res[[i]] <- Split_Cols(w=cts, z=i)
}

# Check if columns were subset correctly
check <- list()
for (i in 1:13){
  check[[i]] <- colnames(res[[i]]) %in% subset(Samples$Sample, Samples$Tissue==Tissues[[i]])
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

# Keep only genes expressed in at least half the samples
Keep <- lapply(y, function(x){
  rowSums(cpm(x)>1)>=11
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
# Comaprisons to test
Pairs <- list(c("Amygdala.Male", "Amygdala.Female"), 
              c("Anterior.Male", "Anterior.Female"), 
              c("Caudate.Male", "Caudate.Female"),
              c("Cerebellar.Male", "Cerebellar.Female"),
              c("Cerebellum.Male", "Cerebellum.Female"),
              c("Cortex.Male", "Cortex.Female"),
              c("Frontal_Cortex.Male", "Frontal_Cortex.Female"),
              c("Hippocampus.Male", "Hippocampus.Female"), 
              c("Hypothalamus.Male", "Hypothalamus.Female"), 
              c("Nucleus_Accumbens.Male", "Nucleus_Accumbens.Female"), 
              c("Putamen.Male", "Putamen.Female"),  
              c("Spinal_Cord.Male", "Spinal_Cord.Female"), 
              c("Substantia_Nigra.Male", "Substantia_Nigra.Female"))

# Apply exact test
Exact_Func <- function(x, comp){
  exactTest(x, comp)
}
Exact_Res <- Map(Exact_Func, y, Pairs)

#---------------------------------------------------------------------------------------------------------------------
# Summarize results 
#---------------------------------------------------------------------------------------------------------------------
# Report sig DGX genes
Up_Reg <- function(x){
  res <- topTags(x, n=Inf, p=0.05)$table
  up <- res[res$logFC > 0, ]
  return(up)
}

Down_Reg <- function(x){
  res <- topTags(x, n=Inf, p=0.05)$table
  down <- res[res$logFC < 0, ]
  return(down)
}
Up_Top <- lapply(Exact_Res, Up_Reg)
Down_Top <- lapply(Exact_Res, Down_Reg)

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

write(Up_Json, "Salmon_Upreg_Exact.json")
write(Down_Json, "Salmon_Downreg_Exact.json")

# Get summary of results as table
Summary_Func <- function(x){
  res <- summary(decideTests(x))
  return(res)
}
Results_df <- lapply(Exact_Res, Summary_Func)

#---------------------------------------------------------------------------------------------------------------------
# Mean-Difference Plots
#---------------------------------------------------------------------------------------------------------------------
# Plot Mean-Difference  plots on one page
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Salmon: Mean-Difference Plots; Exact Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}

# Write to file
pdf(MD_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE)  # margins: c(bottom, left, top, right)
Res_Plots <- Map(MD_Plot_Func, x=Exact_Res, w=Tissues)
legend(50.0, 15.0, legend=c("Up","Not Sig", "Down"), pch = 16, col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Volcano Plots
#---------------------------------------------------------------------------------------------------------------------
# Make df of values for axis
Volcano_Func <- function(x){
  cbind(x$table$logFC, -log10(x$table[,"PValue"]))
}
Volcano_Res <- lapply(Exact_Res, Volcano_Func)

# Coerce to df
Volcano_Res <- lapply(Volcano_Res, as.data.frame)

# Rename columns
colnames <- c("logFC", "negLogPval")

Rename_Cols_Func <- function(x){
  setNames(x, colnames)
}
Volcano_Res <- lapply(Volcano_Res, Rename_Cols_Func)

# Plot
Plot_Func <- function(a, b, c, d){
  plot(a, pch=19, main=b, xlab = '', ylab = '', las = 1)
  with(inner_join(a, c), points(logFC, negLogPval, pch=19, col="green"))
  with(inner_join(a, d), points(logFC, negLogPval, pch=19, col="blue"))
  abline(a=-log10(0.05), b=0, col="blue") 
  abline(v=c(2,-2), col="red")
  mtext('Salmon: Volcano Plots; Exact Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, a=Volcano_Res, b=Tissues, c=Up_Top, d=Down_Top)
legend(25.0, 8.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()

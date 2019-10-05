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
samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Read in count matrix
cts <- read.table(COUNT_MATRIX, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Remove samples not present in metadata
cts <- cts[intersect(colnames(cts), as.character(samples$Sample))]

# Create design matrix: Step 1
# 26 combos of tissue and sex
Group <- factor(paste(samples$Tissue, samples$Sex, sep="."))
cbind(samples,group=Group)

# Create DGEList object: All combos
y <- DGEList(cts, group=Group)

#Create deisgn matrix: Step 2
design <- model.matrix(~0+group, data=y$samples) # No intercept
colnames(design) <- levels(y$samples$group)

# Filter out lowly expressed genes.
# Remove genes w/ <7 counts.
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# TMM Normalization
y <- calcNormFactors(y)

# Estimate common dispersion and tagwise dispersions in one run (recommended)
y <- estimateDisp(y, design, robust=TRUE)

#------------------------------------------------------------------------------------------------------------------
# Test for DGX with Exact Test
#------------------------------------------------------------------------------------------------------------------
# Comparisons to test
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
Exact_Func <- function(x){
  exactTest(y, x)
}
Exact_Res <- lapply(Pairs, Exact_Func)
names(Exact_Res) <- c('Amygdala', 'Anterior', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_Cortex',
                      'Hippocampus', 'Hypothalamus', 'Nucleus_Accumbens', 'Putamen', 'Spinal_Cord', 'Substantia_Nigra')

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

# Male table of up and down regulated genes for each tissue
Tissues <- list('Amygdala', 'Anterior', 'Caudate', 'Cerbellar', 'Cerebellum', 'Cortex', 'Frontal Cortex',
                'Hippocampus', 'Hypothalamus', 'Nucleus Accumbens', 'Putamen', 'Spinal Cord', 'Substantia Nigra')

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

# Get summary of results as tables
Summary_Func <- function(x){
  res <- summary(decideTests(x))
  return(res)
}
Results_df <- lapply(Exact_Res, Summary_Func)

#---------------------------------------------------------------------------------------------------------------------
# Mean-Difference plots
#---------------------------------------------------------------------------------------------------------------------
# To reset par
opar <- par(no.readonly = TRUE) 

# Plot Mean-Difference  plots on one page
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Salmon: Mean-Difference Plots; Exact  Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}
pdf(MD_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE) 
Res_Plots <- Map(MD_Plot_Func, x=Exact_Res, w=Tissues)
legend(26.0, 10.0, legend=c("Up","Not Sig", "Down"), pch = 16, col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
dev.off()
par(opar) 

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

# Add values of zero for tissues that did not return significant results
# None of the tissues returned significant DEGs
Write_Zero <- function(x){
  x <-  data.frame(logFC=c(0), PValue=c(0))
  return(x)
}
Up_Top <- lapply(Up_Top, Write_Zero)
Down_Top <- lapply(Down_Top, Write_Zero)

# Plot
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE) # margins: c(bottom, left, top, right) 
Plot_Func <- function(a, b, c, d){
  plot(a, pch=19, main=b, xlab = '', ylab = '', las = 1)
  with(inner_join(a, c), points(logFC, negLogPval, pch=19, col="green"))
  with(inner_join(a, d), points(logFC, negLogPval, pch=19, col="blue"))
  abline(a=-log10(0.05), b=0, col="blue") 
  abline(v=2, col="red")
  abline(v=-2, col="red")
  mtext('Salmon: Volcano Plots; Exact Ratio Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, a=Volcano_Res.TEST, b=Tissues, c=Up_Top, d=Down_Top)
legend(30.0, 2.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()
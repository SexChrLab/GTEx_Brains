# This script uses an exact test to tests for differential gene expression between males and females 
# within each brain tissue type using transcript level counts on the age matched samples.

# Constants
METADATA <- "/scratch/mjpete11/GTEx/Metadata/Age_Matched_Metadata.csv"
# Hisat/stringtie results are stored in seperate matrices because the same transcripts/genes reported are tissue-specific
PATHS <-c('/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Amygdala_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Anterior_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Caudate_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Cerebellar_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Cerebellum_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Cortex_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/FrontalCortex_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Hippocampus_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Hypothalamus_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/NucleusAccumbens_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/Putamen_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/SpinalCord_Transcript_Hisat_CountMatrix.tsv',
                     '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/SubstantiaNigra_Transcript_Hisat_CountMatrix.tsv')
# Plots/json files
UP_JSON <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Transcript/Hisat_Upreg_Exact.json'
DOWN_JSON <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Transcript/Hisat_Downreg_Exact.json'
MD_PLOT <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Transcript/Hisat_Exact_MD.pdf'
VOLCANO_PLOT <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Transcript/Hisat_Exact_Volcano.pdf'

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(grid)
library(rjson)
library(dplyr)
library(org.Hs.eg.db)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample   

# Count matrices
Tissues <- list('Amygdala', 'Anterior', 'Caudate', 'Cerbellar', 'Cerebellum', 'Cortex', 'Frontal Cortex',
                'Hippocampus', 'Hypothalamus', 'Nucleus Accumbens', 'Putamen', 'Spinal Cord', 'Substantia Nigra')

cts <- lapply(PATHS, function(x){
  t <- read.table(x, sep="\t")
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

write(Up_Json, UP_JSON)
write(Down_Json, DOWN_JSON)

# Get summary of results as table
# Summary_Func <- function(x){
#   res <- summary(decideTests(x))
#   return(res)
# }
# Results_df <- lapply(Exact_Res, Summary_Func)

#---------------------------------------------------------------------------------------------------------------------
# Mean-Difference Plots
#---------------------------------------------------------------------------------------------------------------------
# Plot Mean-Difference  plots on one page
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Hisat: Transcript Mean-Difference Plots; Exact Test', side = 3, outer = TRUE, cex=1.2, line=3)
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
  mtext('Hisat: Transcript Volcano Plots; Exact Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, a=Volcano_Res, b=Tissues, c=Up_Top, d=Down_Top)
legend(40.0, 80.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Gene ontology analysis
#---------------------------------------------------------------------------------------------------------------------
# Example
#qlf <- glmQLFTest(fit, coef=2)
# go <- goana(qlf, species="Hm")
# topGO(go, sort="up")

# # RD
# Gene_Ont <- function(x){
#   res <- goana(x, species="Hm")
#   return(res)
# }
# GO <- lapply(Exact_Res, Gene_Ont)
# 
# TOP_GO <- function(x){
#   res <- topGO(x, sort="up")
#   return(res)
# }
# GO_Res <- lapply(GO, TOP_GO)

# add step to write results



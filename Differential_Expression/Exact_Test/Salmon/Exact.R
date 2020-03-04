# This script looks at differential gene expression between males and females within each brain tissue type.
# gene, age matched
setwd("/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon")

# Input/Output
METADATA <- snakemake@input[[1]]
COUNTS <- snakemake@input[[2]] 
MD_PLOT <- snakemake@output[[1]]
VOLCANO_PLOT <- snakemake@output[[2]]
UP_JSON <- snakemake@output[[3]]
DOWN_JSON <-snakemake@output[[4]]

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(org.Hs.eg.db)
library(rjson)
library(dplyr)

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
cts <- read.table(COUNTS, sep="\t")

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
# Comaprisons to test
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
  plotMD(x, main=w, p.value=0.05, adjust.method="BH", legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Salmon: Gene Mean-Difference Plots; Exact Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}

# Write to file
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
  mtext('Salmon: Gene Volcano Plots; Exact Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, RES=Volcano_Res, TISSUE=Tissues, UP=Volcano_Up, DOWN=Volcano_Down)
legend(10.0, 8.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Gene ontology and pathway enrichment analysis
#---------------------------------------------------------------------------------------------------------------------
# Both use the NCBI RefSeq annotation.
# Convert ensemble annotation to NCBI RefSeq
#library(biomaRt)
#symbols <- mapIds(org.Hs.eg.db, keys = Up_Genes[[1]], keytype = "ENSEMBL", column="SYMBOL")
#
## GO
#Gene_Ont <- function(x){
#  res <- goana(rownames(x), species="Hs")
#  return(res)
#}
#GO <- lapply(Exact_Res, Gene_Ont)
#
## Ontology options: 'MF': molecular function, 'BP': biological process, 'CC': celular component
#TOP_GO <- function(x, w){
#  res <- topGO(x, ontology=c('BP'), sort=rownames(w), number=10)
#  return(res)
#}
#Up_GO_Res <- Map(TOP_GO, x=GO, w=Up_Genes)
#Down_GO_Res <- Map(TOP_GO, x=GO, w=Down_Genes)
#
## KEGG
#keg <- kegga(qlf, species="Mm")
#
#Kegg_Path <- function(x){
#  res <- kegga(rownames(x), species="Hs")
#  return(res)
#}
#KEGG <- lapply(Exact_Res, Kegg_Path)
#
#TOP_KEGG <- function(x, w){
#  res <- topKEGG(x, sort=rownames(w), number=10)
#  return(res)
#}
#Up_Kegg_Res <- Map(TOP_KEGG, x=KEGG, w=Up_Genes)
#Down_Kegg_Res <- Map(TOP_KEGG, x=KEGG, w=Down_Genes)





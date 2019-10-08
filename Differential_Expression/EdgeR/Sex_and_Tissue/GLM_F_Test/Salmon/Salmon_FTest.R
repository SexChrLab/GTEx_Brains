# This script looks at differential gene expression between males and females within each brain tissue type.
setwd("/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon")

METADATA <- "/scratch/mjpete11/GTEx/Metadata/Matched_Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Count_Matrix.tsv"
MD_PLOT <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Plots/Matched_Salmon_FTest_MD.pdf'
VOLCANO_PLOT <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Plots/Matched_Salmon_FTest_Volcano.pdf'

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(grid)
library(rjson)

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

# Read in counts 
cts <- read.csv(COUNT_MATRIX, sep = "\t")

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
  res <- model.matrix(~0 + fc, data = df$Samples)
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

# Filter out lowly expressed genes.
# Remove genes w/ <7 counts.
Keep <- lapply(y, function(x){
  rowSums(cpm(x)>1)>=2
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

#------------------------------------------------------------------------------------------------------------------
# Test for DGX with GLM F Test
#------------------------------------------------------------------------------------------------------------------
# Fit glm model
Fit_Func <- function(a, b){
  fit <- glmQLFit(a, b, robust=TRUE)
  return(fit)
}
Fit <- Map(Fit_Func, a=y, b=Design)

# Make contrasts: Sex by Tissue
Contrasts <- c('Amygdala.Female - Amygdala.Male',
               'Anterior.Female - Anterior.Male',
               'Caudate.Female - Caudate.Male',
               'Cerebellar.Female - Cerebellar.Male',
               'Cerebellum.Female - Cerebellum.Male',
               'Cortex.Female - Cortex.Male',
               'Frontal_Cortex.Female - Frontal_Cortex.Male',
               'Hippocampus.Female - Hippocampus.Male',
               'Hypothalamus.Female - Hypothalamus.Male',
               'Nucleus_Accumbens.Female - Nucleus_Accumbens.Male',
               'Putamen.Female - Putamen.Male',
               'Spinal_Cord.Female - Spinal_Cord.Male',
               'Substantia_Nigra.Female - Substantia_Nigra.Male')

Names <- c('Am.F.vs.M', 'At.F.vs.M', 'Ca.F.vs.M', 'Ce.F.vs.M', 'Co.F.vs.M', 'Fc.F.vs.M', 'Cm.F.vs.M', 
           'Hp.F.vs.M', 'Hy.F.vs.M', 'Na.F.vs.M', 'Pu.F.vs.M', 'Sp.F.vs.M', 'Sn.F.vs.M')

# Make contrasts
Contrast_Func <- function(a, b){
  res <- makeContrasts(contrasts = a, levels = colnames(b))
  return(res)
}
my.contrasts <- Map(Contrast_Func, a=Contrasts, b=Design)
names(my.contrasts) <- Names

# Apply GLM F test
GLM_Ratio_Func <- function(a, b){
  glmQLFTest(a, contrast=b)
}
GLM_Res <- Map(GLM_Ratio_Func, a=Fit, b=my.contrasts)
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
Up_Top <- lapply(GLM_Res, Up_Reg)
Down_Top <- lapply(GLM_Res, Down_Reg)

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

write(Up_Json, "Salmon_Upreg_FTest.json")
write(Down_Json, "Salmon_Downreg_FTest.json")

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
  mtext('Salmon: Mean-Difference Plots; GLM F Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}

# Write to file
pdf(MD_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE)  # margins: c(bottom, left, top, right)
Res_Plots <- Map(MD_Plot_Func, x=GLM_Res, w=Tissues)
legend(50.0, 15.0, legend=c("Up","Not Sig", "Down"), pch = 16, col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Volcano Plots
#---------------------------------------------------------------------------------------------------------------------
# Make df of values for axis
Volcano_Func <- function(x){
  cbind(x$table$logFC, -log10(x$table[,"PValue"]))
}
Volcano_Res <- lapply(GLM_Res, Volcano_Func)

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
  mtext('Salmon: Volcano Plots; GLM F Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, a=Volcano_Res, b=Tissues, c=Up_Top, d=Down_Top)
legend(25.0, 8.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()

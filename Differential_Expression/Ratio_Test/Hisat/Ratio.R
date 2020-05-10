# This script looks at differential gene expression between males and females within each brain tissue type.
# Gene, age matched

# Constants
METADATA <- snakemake@input[[1]]
# Hisat/stringtie results are stored in seperate matrices because the same transcripts/genes reported are tissue-specific
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
PATHS <- c(PATHS.1,PATHS.2,PATHS.3,PATHS.4,PATHS.5,PATHS.6,PATHS.7,
           PATHS.8,PATHS.9,PATHS.10,PATHS.11,PATHS.12,PATHS.13)

# Plots/json files
MD_PLOT <- snakemake@output[[1]]
VOLCANO_PLOT <- snakemake@output[[2]]
UP_JSON <- snakemake@output[[3]]
DOWN_JSON <- snakemake@output[[4]]
TABLE.1 <- snakemake@output[[5]]
TABLE.2 <- snakemake@output[[6]]
TABLE.3 <- snakemake@output[[7]]
TABLE.4 <- snakemake@output[[8]]
TABLE.5 <- snakemake@output[[9]]
TABLE.6 <- snakemake@output[[10]]
TABLE.7 <- snakemake@output[[11]]
TABLE.8 <- snakemake@output[[12]]
TABLE.9 <- snakemake@output[[13]]
TABLE.10 <- snakemake@output[[14]]
TABLE.11 <- snakemake@output[[15]]
TABLE.12 <- snakemake@output[[16]]
TABLE.13 <- snakemake@output[[17]]

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
Tissues <- list('Amygdala', 'Anterior', 'Caudate', 'Cerbellar', 'Cerebellum', 'Cortex', 'Frontal Cortex',
                'Hippocampus', 'Hypothalamus', 'Nucleus Accumbens', 'Putamen', 'Spinal Cord', 'Substantia Nigra')

cts <- lapply(PATHS, function(x){
  t <- read.csv(x, sep=",")
})
names(cts) <- Tissues

# Replace . to - in colnames in each df
for (i in seq_along(cts)){
  colnames(cts[[i]]) <- str_replace_all(colnames(cts[[i]]), pattern = "\\.","-")
}

# Remove "-" at end of caudate sample colnames
colnames(cts[[3]]) <- str_remove(colnames(cts[[3]]),"-$")

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

# Sort cols in cts in same order as rows in Meta$Samples
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

#------------------------------------------------------------------------------------------------------------------
# Test for DGX with GLM Ratio Test
#------------------------------------------------------------------------------------------------------------------
# Fit glm model
Fit_Func <- function(a, b){
  fit <- glmFit(a, b, robust=TRUE)
  return(fit)
}
Fit <- Map(Fit_Func, a=y, b=Design)

# Make contrasts: Sex by Tissue
Contrasts <- c('Amygdala.Male - Amygdala.Female',
               'Anterior.Male - Anterior.Female',
               'Caudate.Male - Caudate.Female',
               'Cerebellar.Male - Cerebellar.Female',
               'Cerebellum.Male - Cerebellum.Female',
               'Cortex.Male - Cortex.Female',
               'Frontal_Cortex.Male - Frontal_Cortex.Female',
               'Hippocampus.Male - Hippocampus.Female',
               'Hypothalamus.Male - Hypothalamus.Female',
               'Nucleus_Accumbens.Male - Nucleus_Accumbens.Female',
               'Putamen.Male - Putamen.Female',
               'Spinal_Cord.Male - Spinal_Cord.Female',
               'Substantia_Nigra.Male - Substantia_Nigra.Female')

Names <- c('Am.MvsF', 'At.MvsF', 'Ca.MvsF', 'Ce.MvsF', 'Co.MvsF', 'Fc.MvsF', 'Cm.MvsF', 
           'Hp.MvsF', 'Hy.MvsF', 'Na.MvsF', 'Pu.MvsF', 'Sp.MvsF', 'Sn.MvsF')

# Make contrasts
Contrast_Func <- function(a, b){
  res <- makeContrasts(contrasts = a, levels = colnames(b))
  return(res)
}
my.contrasts <- Map(Contrast_Func, a=Contrasts, b=Design)
names(my.contrasts) <- Names

# Apply GLM F test
# Resulting objects are of class DGELRT 
GLM_Ratio_Func <- function(a, b){
  glmLRT(a, contrast=b)
}
GLM_Res <- Map(GLM_Ratio_Func, a=Fit, b=my.contrasts)

# Write table from test results object to file for mashr analysis
Snake_Var <- c(TABLE.1, TABLE.2, TABLE.3, TABLE.4, TABLE.5, TABLE.6,
               TABLE.7, TABLE.8, TABLE.9, TABLE.10, TABLE.11, TABLE.12, TABLE.13)
Table_Func <- function(tabl, name){
    res <- write.table(tabl[['table']], 
                       file=name) 
    return(res)
}
Map(Table_Func, tabl=GLM_Res, name=Snake_Var)

#---------------------------------------------------------------------------------------------------------------------
# Summary stats
#---------------------------------------------------------------------------------------------------------------------
# Function to correct for multiple testing
Test_Correct <- function(x){
  x[['table']][['PValue']] <- p.adjust(x[['table']][['PValue']],method="BH")
  return(x)
}
Corrected_RGLM <- lapply(GLM_Res, Test_Correct)

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
Up_Top <- lapply(Corrected_RGLM, Up_Reg)
Down_Top <- lapply(Corrected_RGLM, Down_Reg)

# Male table of up and down regulated genes for each tissue
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
  mtext('Hisat: Isoform Mean-Difference Plots; GLM Ratio Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}

# Write to file
pdf(MD_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE)  # margins: c(bottom, left, top, right)
Res_Plots <- Map(MD_Plot_Func, x=GLM_Res, w=Tissues)
legend(20.0, 0.0, legend=c("Up","Not Sig", "Down"), pch = 16, col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
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
Volcano_Res <- lapply(GLM_Res, Untrans_Volcano)

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
  mtext('Hisat: Isoform Volcano Plots; GLM Ratioi Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, RES=Volcano_Res, TISSUE=Tissues, UP=Volcano_Up, DOWN=Volcano_Down)
legend(10.0, 8.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()


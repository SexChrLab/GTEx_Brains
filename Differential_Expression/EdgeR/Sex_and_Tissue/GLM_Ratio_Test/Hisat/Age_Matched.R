# This script looks at differential gene expression between males and females within each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata/Age_Matched_Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Count_Matrix.tsv"
PATHS <- c('/scratch/mjpete11/GTEx/Amygdala/Hisat_Stringtie/gene_count_matrix.csv', 
           '/scratch/mjpete11/GTEx/Anterior/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Caudate/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Cerebellar/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Cerebellum/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Cortex/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Frontal_Cortex/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Hippocampus/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Hypothalamus/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Nucleus_Accumbens/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Putamen/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Spinal_Cord/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Substantia_Nigra/Hisat_Stringtie/gene_count_matrix.csv')

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
samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Count matrices
cts <- lapply(PATHS, function(x) {
  t <- as.matrix(read.csv(x, row.names="gene_id"))
})
names(cts) <- c('Amygdala', 'Anterior', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_Cortex',
                'Hippocampus', 'Hypothalamus', 'Nucleus_Accumbens', 'Putamen', 'Spinal_Cord', 'Substantia_Nigra')

# Remove the . at the end of the Caudate sample names
colnames(cts$Caudate) <- substring(colnames(cts$Caudate), 1, nchar(colnames(cts$Caudate))-1)

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

# Create list of DGEList objects
DGE_Func <- function(df, lst){
  res <- DGEList(df, group=lst)
  return(res)
}
y <- Map(DGE_Func, cts, Groups)

# Create design matrix: Step 2
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

# Filter out lowly expressed genes.
# Remove genes w/ <7 counts.
Keep <- lapply(y, function(x){
  rowSums(cpm(x)>1)>=2
})
#summary(Keep$Amygdala) # example

Filter_Func <- function(x, k){
  x <- x[k, , keep.lib.sizes=FALSE]
}
y <- Map(Filter_Func, y, Keep)

# TMM Normalization
y <- lapply(y, calcNormFactors)
#View(y$Amygdala$samples)

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
  glmLRT(a, contrast=b)
}
GLM_Res <- Map(GLM_Ratio_Func, a=Fit, b=my.contrasts)

#---------------------------------------------------------------------------------------------------------------------
# Summary stats
#---------------------------------------------------------------------------------------------------------------------
# Get summary of results as tables
Summary_Func <- function(x){
  res <- summary(decideTests(x))
  return(res)
}
Results_df <- lapply(GLM_Res, Summary_Func)

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

write(Up_Json, "Hisat_Upreg_Ratio.json")
write(Down_Json, "Hisat_Downreg_Ratio.json")

#---------------------------------------------------------------------------------------------------------------------
# Plots
#---------------------------------------------------------------------------------------------------------------------
# To reset par
opar <- par(no.readonly = TRUE) 

# Plot Mean-Difference  plots on one page
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Hisat: Mean-Difference Plots; GLM Ratio Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}
Res_Plots <- Map(MD_Plot_Func, x=GLM_Res, w=Tissues)
legend(26.0, 10.0, legend=c("Up","Not Sig", "Down"), pch = 16, col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
#par(opar) 

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

# Volcano plots
# Subset pos and neg sig values
Subset_Func <- function(x){
  significant <- list()
  up <- subset(x,  negLogPval >= -log10(0.05) & logFC > 0)
  down <- subset(x, negLogPval >= -log10(0.05) & logFC < 0)
  significant <- append(significant, list(up))
  significant <- append(significant, list(down))
  return(significant)
}
Subset_Res <- lapply(Volcano_Res, Subset_Func)

# Add colored points for sig genes
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE) # margins: c(bottom, left, top, right) 
Plot_Func <- function(a, b, c){
  plot(a, pch=19, main=b, xlab = '', ylab = '', las = 1)
  with(inner_join(a, c[[1]]), points(logFC, negLogPval, pch=19, col="green"))
  with(inner_join(a, c[[2]]), points(logFC, negLogPval, pch=19, col="blue"))
  abline(a=-log10(0.05), b=0, col="blue") 
  abline(v=0, col="red")
  mtext('Hisat: Volcano Plots; GLM Ratio Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
Map(Plot_Func, a=Volcano_Res, b=Tissues, c=Subset_Res)
legend(16.0, 9.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)




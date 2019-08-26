# This script looks at differential gene expression between males and females within each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata/Age_Matched_Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Count_Matrix.tsv"

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

# Read in counts 
cts <- read.csv(COUNT_MATRIX, sep = "\t")

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
# Remove genes < 1 CPM.
keep <- rowSums(cpm(y)>1) >= 2 # Must be expressed in 2 or more libraries
y <- y[keep, , keep.lib.sizes=FALSE] # Recalculate lib sizes 

# TMM Normalization
y <- calcNormFactors(y)

# Estimate common dispersion and tagwise dispersions in one run (recommended)
y <- estimateDisp(y, design, robust=TRUE)

#------------------------------------------------------------------------------------------------------------------
# Test for DGX with GLM F Test
#------------------------------------------------------------------------------------------------------------------
# Make contrasts: Sex by Tissue
my.contrasts <- makeContrasts(Am.F.vs.M = Amygdala.Female - Amygdala.Male,
                              At.F.vs.M = Anterior.Female - Anterior.Male,
                              Ca.F.vs.M = Caudate.Female - Caudate.Male,
                              Ce.F.vs.M = Cerebellar.Female - Cerebellar.Male,
                              Co.F.vs.M = Cortex.Female - Cortex.Male,
                              Fc.F.vs.M = Frontal_Cortex.Female - Frontal_Cortex.Male,
                              Cm.F.vs.M = Cerebellum.Female - Cerebellum.Male,
                              Hp.F.vs.M = Hippocampus.Female - Hippocampus.Male,
                              Hy.F.vs.M = Hypothalamus.Female - Hypothalamus.Male,
                              Na.F.vs.M = Nucleus_Accumbens.Female - Nucleus_Accumbens.Male,
                              Pu.F.vs.M = Putamen.Female - Putamen.Male,
                              Sp.F.vs.M = Spinal_Cord.Female - Spinal_Cord.Male,
                              Sn.F.vs.M = Substantia_Nigra.Female - Substantia_Nigra.Male,
                              levels=design)

# Fit glm model
fit <- glmQLFit(y, design, robust=TRUE)

Contrasts <- c('Am.F.vs.M', 'At.F.vs.M', 'Ca.F.vs.M', 'Ce.F.vs.M', 'Co.F.vs.M', 'Fc.F.vs.M', 'Cm.F.vs.M', 
               'Hp.F.vs.M', 'Hy.F.vs.M', 'Na.F.vs.M', 'Pu.F.vs.M', 'Sp.F.vs.M', 'Sn.F.vs.M')

# Apply GLM F test
GLM_Func <- function(x){
  glmQLFTest(fit, contrast=my.contrasts[,x])
}
GLM_Res <- lapply(Contrasts, GLM_Func)
names(GLM_Res) <- c('Amygdala', 'Anterior', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_Cortex',
                      'Hippocampus', 'Hypothalamus', 'Nucleus_Accumbens', 'Putamen', 'Spinal_Cord', 'Substantia_Nigra')

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

# Make table of up and down regulated genes for each tissue
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

write(Up_Json, "Salmon_Upreg_FTest.json")
write(Down_Json, "Salmon_Downreg_FTest.json")

#---------------------------------------------------------------------------------------------------------------------
# Plots
#---------------------------------------------------------------------------------------------------------------------
# To reset par
opar <- par(no.readonly = TRUE) 

# Plot Mean-Difference  plots on one page
par(mfrow = c(3, 5), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.4)
  mtext('Salmon: Mean-Difference Plots; GLM F Test', side = 3, outer = TRUE, cex=1.2, line=3)
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
  mtext('Salmon: Volcano Plots; GLM F Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
Map(Plot_Func, a=Volcano_Res, b=Tissues, c=Subset_Res)
legend(16.0, 9.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)



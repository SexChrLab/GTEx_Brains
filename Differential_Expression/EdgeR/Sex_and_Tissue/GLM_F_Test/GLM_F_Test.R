# This script looks at differential gene expression between males and females within each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv"
PLOT_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Gene/Normalized/GLM_F_Test/Plots/"

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

# Read in counts 
cts <- read.csv(COUNT_MATRIX, sep = "\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Create design matrix: Step 1
# 26 combos of tissue and sex
Group <- factor(paste(samples$Tissue, samples$Sex, sep="."))
cbind(samples,group=Group)

# Create DGEList object: All combos
y <- DGEList(cts, group=Group)

#Create deisgn matrix: Step 2
design <- model.matrix(~0+group, data=y$samples) # No intercept
colnames(design) <- levels(y$samples$group)
design

# Filter out lowly expressed genes.
# Remove genes < 1 CPM.
keep <- rowSums(cpm(y)>1) >= 2 # Must be expressed in 2 or more libraries
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE] # Recalculate lib sizes 

# TMM Normalization
y <- calcNormFactors(y)

# Estimate common dispersion and tagwise dispersions in one run (recommended)
y <- estimateDisp(y, design, robust=TRUE)
sqrt(y$common.dispersion) 

# Data Exploration

# Plot
setwd(PLOT_DIR)

pdf('Normalized_GLM_F_Test.pdf')

# Barplot of libsizes
library_plot <- barplot(y$samples$lib.size*1e-6, names=1:1340, ylab="Library size (millions)", main="Barplot of sample library sizes", xlab="sample")

# BCV plot
plotBCV(y, main="Biological Coefficient of Variation Plot")

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

# Test for DGX using quasi-liklihood.

# Fit glm model
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit, main="Quasi-Likelihood Disperson")

# Amygdala Female vs Male F Test
qlf.Am.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Am.F.vs.M"])
df_Am <- summary(decideTests(qlf.Am.F.vs.M))
grid.newpage() # To keep summary table from plotting on top of other plots
grid.table(df_Am)
plotMD(qlf.Am.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Am.F.vs.M$table$logFC, -log10(qlf.Am.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala.Female-1*Amygdala.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(qlf.Am.F.vs.M$table[,"PValue"], breaks=50, main="Amygdala p-value frequency histogram")

# Anterior Female vs Male F Test
qlf.At.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"At.F.vs.M"])
df_At <- summary(decideTests(qlf.At.F.vs.M))
grid.newpage()
grid.table(df_At)
plotMD(qlf.At.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.At.F.vs.M$table$logFC, -log10(qlf.At.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior.Female-1*Anterior.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.F.vs.M$table[,"PValue"], breaks=50, main="Anterior p-value frequency histogram")

# Cortex Female vs Male F Test
qlf.Co.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Co.F.vs.M"])
df_Co <- summary(decideTests(qlf.Co.F.vs.M))
grid.newpage()
grid.table(df_Co)
plotMD(qlf.Co.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Co.F.vs.M$table$logFC, -log10(qlf.Co.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex.Female-1*Cortex.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.F.vs.M$table[,"PValue"], breaks=50, main="Cortex p-value frequency histogram")

# Cerebellum Female vs Male F Test
qlf.Cm.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Cm.F.vs.M"])
df_Cm <- summary(decideTests(qlf.Cm.F.vs.M))
grid.newpage()
grid.table(df_Cm)
plotMD(qlf.Cm.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Cm.F.vs.M$table$logFC, -log10(qlf.Cm.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellum.Female-1*Cerebellum.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Cm.F.vs.M$table[,"PValue"], breaks=50, main="Cerebellum p-value frequency histogram")

# Cerebellar Female vs Male
qlf.Ce.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.F.vs.M"])
df_Ce <- summary(decideTests(qlf.Ce.F.vs.M))
grid.newpage()
grid.table(df_Ce)
plotMD(qlf.Ce.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Ce.F.vs.M$table$logFC, -log10(qlf.Ce.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar.Female-1*Cerebellar.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.F.vs.M$table[,"PValue"], breaks=50, main="Cerebellar p-value frequency histogram")

# Hippocampus Female vs Male F Test
qlf.Hp.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Hp.F.vs.M"])
df_Hp <- summary(decideTests(qlf.Hp.F.vs.M))
grid.newpage()
grid.table(df_Hp)
plotMD(qlf.Hp.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Hp.F.vs.M$table$logFC, -log10(qlf.Hp.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus.Female-1*Hippocampus.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hp.F.vs.M$table[,"PValue"], breaks=50, main="Hippocampus p-value frequency histogram")

# Hypothalamus Female vs Male F Test
qlf.Hy.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Hy.F.vs.M"])
df_Hy <- summary(decideTests(qlf.Hy.F.vs.M))
grid.newpage()
grid.table(df_Hy)
plotMD(qlf.Hy.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Hy.F.vs.M$table$logFC, -log10(qlf.Hy.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus.Female-1*Hypothalamus.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hy.F.vs.M$table[,"PValue"], breaks=50, main="Hypothalamus p-value frequency histogram")

# Frontal Cortex Female vs Male F Test
qlf.Fc.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.F.vs.M"])
df_Fc <- summary(decideTests(qlf.Fc.F.vs.M))
grid.newpage()
grid.table(df_Fc)
plotMD(qlf.Fc.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Fc.F.vs.M$table$logFC, -log10(qlf.Fc.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex.Female-1*Frontal_Cortex.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.F.vs.M$table[,"PValue"], breaks=50, main="Frontal cortex p-value frequency histogram")

# Nucleus Accumbens Female vs Male F Test
qlf.Na.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Na.F.vs.M"])
df_Na <- summary(decideTests(qlf.Na.F.vs.M))
grid.newpage()
grid.table(df_Na)
plotMD(qlf.Na.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Na.F.vs.M$table$logFC, -log10(qlf.Na.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens.Female-1*Nucleus_Accumbens.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Na.F.vs.M$table[,"PValue"], breaks=50, main="Nucleus accumbens p-value frequency histogram")

# Putamen Female vs Male F Test
qlf.Pu.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Pu.F.vs.M"])
df_Pu <- summary(decideTests(qlf.Pu.F.vs.M))
grid.newpage()
grid.table(df_Pu)
plotMD(qlf.Pu.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Pu.F.vs.M$table$logFC, -log10(qlf.Pu.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen.Female-1*Putamen.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Pu.F.vs.M$table[,"PValue"], breaks=50, main="Putamen p-value frequency histogram")

# Spinal Cord Female vs Male F Test
qlf.Sp.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Sp.F.vs.M"])
df_Sp <- summary(decideTests(qlf.Sp.F.vs.M))
grid.newpage()
grid.table(df_Sp)
plotMD(qlf.Sp.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Sp.F.vs.M$table$logFC, -log10(qlf.Sp.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Spinal_Cord.Female-1*Spinal_Cord.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Sp.F.vs.M$table[,"PValue"], breaks=50, main="Spinal cord p-value frequency histogram")

# Substantiaa Nigra Female vs Male
qlf.Sn.F.vs.M <- glmQLFTest(fit, contrast=my.contrasts[,"Sn.F.vs.M"])
df_Sn <- summary(decideTests(qlf.Sn.F.vs.M))
grid.newpage()
grid.table(df_Sn)
plotMD(qlf.Sn.F.vs.M)

# Volcano plot
volcanoData <- cbind(qlf.Sn.F.vs.M$table$logFC, -log10(qlf.Sn.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Substantia_Nigra.Female-1*Substantia_Nigra.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Sn.F.vs.M$table[,"PValue"], breaks=50, main="Substantia nigra p-value frequency histogram")

dev.off()

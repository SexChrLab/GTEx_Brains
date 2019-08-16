# This script looks at differential gene expression between males and females within each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv"
PLOT_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Gene/Normalized/GLM_Ratio_Test/Plots/"

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
# Remove genes w/ <7 counts.
keep <- rowSums(cpm(y)>1) >= 2
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

# TMM Normalization
y <- calcNormFactors(y)
y$samples

# Estimate common dispersion and tagwise dispersions in one run (recommended)
y <- estimateDisp(y, design, robust=TRUE)

# Test for DGX using liklihood ratio test.
fit <- glmFit(y, design, robust=TRUE)
head(fit$coefficients)

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

# Plot
setwd(PLOT_DIR)

pdf('Normalized_GLM_Ratio_Test.pdf')

# Amygdala Female vs Male Ratio Test
lrt.Am.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Am.F.vs.M"])
df_Am <- summary(decideTests(lrt.Am.F.vs.M))
grid.newpage()
grid.table(df_Am)
plotMD(lrt.Am.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Am.F.vs.M$table$logFC, -log10(lrt.Am.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala.Female-1*Amygdala.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Am.F.vs.M$table[,"PValue"], breaks=50, main="Amygdala p-value frequency histogram")

# Anterior Female vs Male
lrt.At.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"At.F.vs.M"])
df_At <- summary(decideTests(lrt.At.F.vs.M))
grid.newpage()
grid.table(df_At)
plotMD(lrt.At.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.At.F.vs.M$table$logFC, -log10(lrt.At.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior.Female-1*Anterior.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.At.F.vs.M$table[,"PValue"], breaks=50, main="Anterior p-value frequency histogram")

# Cortex Female vs Male
lrt.Co.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Co.F.vs.M"])
df_Co <- summary(decideTests(lrt.Co.F.vs.M))
grid.newpage()
grid.table(df_Co)
plotMD(lrt.Co.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Co.F.vs.M$table$logFC, -log10(lrt.Co.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex.Female-1*Cortex.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Co.F.vs.M$table[,"PValue"], breaks=50, main="Cortex p-value frequency histogram")

# Cerebellum Female vs Male
lrt.Cm.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Cm.F.vs.M"])
df_Cm <- summary(decideTests(lrt.Cm.F.vs.M))
grid.newpage()
grid.table(df_Cm)
plotMD(lrt.Cm.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Cm.F.vs.M$table$logFC, -log10(lrt.Cm.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellum.Female-1*Cerebellum.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Cm.F.vs.M$table[,"PValue"], breaks=50, main="Cerebellum p-value frequency histogram")

# Cerebellar Female vs Male
lrt.Ce.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Ce.F.vs.M"])
df_Ce <- summary(decideTests(lrt.Ce.F.vs.M))
grid.newpage()
grid.table(df_Ce)
plotMD(lrt.Ce.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Ce.F.vs.M$table$logFC, -log10(lrt.Ce.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar.Female-1*Cerebellar.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Ce.F.vs.M$table[,"PValue"], breaks=50, main="Cerebellar p-value frequency histogram")

# Hippocampus Female vs Male
lrt.Hp.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Hp.F.vs.M"])
df_Hp <- summary(decideTests(lrt.Hp.F.vs.M))
grid.newpage()
grid.table(df_Hp)
plotMD(lrt.Hp.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Hp.F.vs.M$table$logFC, -log10(lrt.Hp.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus.Female-1*Hippocampus.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Hp.F.vs.M$table[,"PValue"], breaks=50, main="Hippocampus p-value frequency histogram")

# Hypothalamus Female vs Male
lrt.Hy.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Hy.F.vs.M"])
df_Hy <- summary(decideTests(lrt.Hy.F.vs.M))
grid.newpage()
grid.table(df_Hy)
plotMD(lrt.Hy.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Hy.F.vs.M$table$logFC, -log10(lrt.Hy.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus.Female-1*Hypothalamus.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Hy.F.vs.M$table[,"PValue"], breaks=50, main="Hypothalamus p-value frequency histogram")

# Frontal Cortex Female vs Male
lrt.Fc.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Fc.F.vs.M"])
df_Fc <- summary(decideTests(lrt.Fc.F.vs.M))
grid.newpage()
grid.table(df_Fc)
plotMD(lrt.Fc.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Fc.F.vs.M$table$logFC, -log10(lrt.Fc.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex.Female-1*Frontal_Cortex.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Fc.F.vs.M$table[,"PValue"], breaks=50, main="Frontal cortex p-value frequency histogram")

# Nucleus Accumbens Female vs Male
lrt.Na.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Na.F.vs.M"])
df_Na <- summary(decideTests(lrt.Na.F.vs.M))
grid.newpage()
grid.table(df_Na)
plotMD(lrt.Na.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Na.F.vs.M$table$logFC, -log10(lrt.Na.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens.Female-1*Nucleus_Accumbens.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Na.F.vs.M$table[,"PValue"], breaks=50, main="Nucleus accumbens p-value frequency histogram")

# Putamen Female vs Male
lrt.Pu.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Pu.F.vs.M"])
df_Pu <- summary(decideTests(lrt.Pu.F.vs.M))
grid.newpage()
grid.table(df_Pu)
plotMD(lrt.Pu.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Pu.F.vs.M$table$logFC, -log10(lrt.Pu.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen.Female-1*Putamen.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Pu.F.vs.M$table[,"PValue"], breaks=50, main="Putamen p-value frequency histogram")

# Spinal Cord Female vs Male
lrt.Sp.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Sp.F.vs.M"])
df_Sp <- summary(decideTests(lrt.Sp.F.vs.M))
grid.newpage()
grid.table(df_Sp)
plotMD(lrt.Sp.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Sp.F.vs.M$table$logFC, -log10(lrt.Sp.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Spinal_Cord.Female-1*Spinal_Cord.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Sp.F.vs.M$table[,"PValue"], breaks=50, main="Spinal cord p-value frequency histogram")

# Substantiaa Nigra Female vs Male
lrt.Sn.F.vs.M <- glmLRT(fit, contrast=my.contrasts[,"Sn.F.vs.M"])
df_Sn <- summary(decideTests(lrt.Sn.F.vs.M))
grid.newpage()
grid.table(df_Sn)
plotMD(lrt.Sn.F.vs.M)

# Volcano plot
volcanoData <- cbind(lrt.Sn.F.vs.M$table$logFC, -log10(lrt.Sn.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Substantia_Nigra.Female-1*Substantia_Nigra.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Histogram of unadjusted p-values
hist(lrt.Sn.F.vs.M$table[,"PValue"], breaks=50, main="Substantia nigra p-value frequency histogram")

dev.off()


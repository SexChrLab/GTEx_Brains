# This script looks at differential gene expression between males and females within each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv"
PLOT_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Gene/Normalized/Exact_Test/Plots/"

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
samples

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Read in count matrix
cts <- read.table(COUNT_MATRIX, sep="\t")

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
sqrt(y$common.dispersion) 

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

# Test for DGX with Exact Test

# Plot
setwd(PLOT_DIR)

pdf('Normalized_Exact_Test.pdf')

# Amygdala Female vs Male
et.Am.F.vs.M <- exactTest(y, c("Amygdala.Male", "Amygdala.Female"))
df_Am <- summary(decideTests(et.Am.F.vs.M))
grid.newpage() # To keep summary table from plotting on top of other plots
grid.table(df_Am)
plotMD(et.Am.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Am.F.vs.M$table$logFC, -log10(et.Am.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala.Female-1*Amygdala.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.F.vs.M$table[,"PValue"], breaks=50, main="Amygdala p-value frequency histogram")

# Anterior Female vs Male
et.At.F.vs.M <- exactTest(y, c("Anterior.Male", "Anterior.Female"))
df_At <- summary(decideTests(et.At.F.vs.M))
grid.newpage()
grid.table(df_At)
plotMD(et.At.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.At.F.vs.M$table$logFC, -log10(et.At.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior.Female-1*Anterior.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.F.vs.M$table[,"PValue"], breaks=50, main="Anterior p-value frequency histogram")

# Cortex Female vs Male
et.Co.F.vs.M <- exactTest(y, c("Cortex.Male", "Cortex.Female"))
df_Co <- summary(decideTests(et.Co.F.vs.M))
grid.newpage()
grid.table(df_Co)
plotMD(et.Co.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Co.F.vs.M$table$logFC, -log10(et.Co.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex.Female-1*Cortex.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.F.vs.M$table[,"PValue"], breaks=50, main="Cortex p-value frequency histogram")

# Cerebellum Female vs Male
et.Cm.F.vs.M <- exactTest(y, c("Cerebellum.Male", "Cerebellum.Female"))
df_Cm <- summary(decideTests(et.Cm.F.vs.M))
grid.newpage()
grid.table(df_Cm)
plotMD(et.Cm.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Cm.F.vs.M$table$logFC, -log10(et.Cm.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellum.Female-1*Cerebellum.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Cm.F.vs.M$table[,"PValue"], breaks=50, main="Cerebellum p-value frequency histogram")

# Cerebellar Female vs Male
et.Ce.F.vs.M <- exactTest(y, c("Cerebellar.Male", "Cerebellar.Female"))
df_Ce <- summary(decideTests(et.Ce.F.vs.M))
grid.newpage()
grid.table(df_Ce)
plotMD(et.Ce.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Ce.F.vs.M$table$logFC, -log10(et.Ce.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar.Female-1*Cerebellar.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.F.vs.M$table[,"PValue"], breaks=50, main="Cerebellar p-value frequency histogram")

# Hippocampus Female vs Male
et.Hp.F.vs.M <- exactTest(y, c("Hippocampus.Male", "Hippocampus.Female"))
df_Hp <- summary(decideTests(et.Hp.F.vs.M))
grid.newpage()
grid.table(df_Hp)
plotMD(et.Hp.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Hp.F.vs.M$table$logFC, -log10(et.Hp.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus.Female-1*Hippocampus.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hp.F.vs.M$table[,"PValue"], breaks=50, main="Hippocampus p-value frequency histogram")

# Hypothalamus Female vs Male
et.Hy.F.vs.M <- exactTest(y, c("Hypothalamus.Male", "Hypothalamus.Female"))
df_Hy <- summary(decideTests(et.Hy.F.vs.M))
grid.newpage()
grid.table(df_Hy)
plotMD(et.Hy.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Hy.F.vs.M$table$logFC, -log10(et.Hy.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus.Female-1*Hypothalamus.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hy.F.vs.M$table[,"PValue"], breaks=50, main="Hypothalamus p-value frequency histogram")

# Frontal Cortex Female vs Male
et.Fc.F.vs.M <- exactTest(y, c("Frontal_Cortex.Male", "Frontal_Cortex.Female"))
df_Fc <- summary(decideTests(et.Fc.F.vs.M))
grid.newpage()
grid.table(df_Fc)
plotMD(et.Fc.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Fc.F.vs.M$table$logFC, -log10(et.Fc.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex.Female-1*Frontal_Cortex.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.F.vs.M$table[,"PValue"], breaks=50, main="Frontal cortex p-value frequency histogram")

# Nucleus Accumbens Female vs Male
et.Na.F.vs.M <- exactTest(y, c("Nucleus_Accumbens.Male", "Nucleus_Accumbens.Female"))
df_Na <- summary(decideTests(et.Na.F.vs.M))
grid.newpage()
grid.table(df_Na)
plotMD(et.Na.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Na.F.vs.M$table$logFC, -log10(et.Na.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens.Female-1*Nucleus_Accumbens.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Na.F.vs.M$table[,"PValue"], breaks=50, main="Nucleus accumbens p-value frequency histogram")

# Putamen Female vs Male
et.Pu.F.vs.M <- exactTest(y, c("Putamen.Male", "Putamen.Female"))
df_Pu <- summary(decideTests(et.Pu.F.vs.M))
grid.newpage()
grid.table(df_Pu)
plotMD(et.Pu.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Pu.F.vs.M$table$logFC, -log10(et.Pu.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen.Female-1*Putamen.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Pu.F.vs.M$table[,"PValue"], breaks=50, main="Putamen p-value frequency histogram")

# Spinal Cord Female vs Male
et.Sp.F.vs.M <- exactTest(y, c("Spinal_Cord.Male", "Spinal_Cord.Female"))
df_Sp <- summary(decideTests(et.Sp.F.vs.M))
grid.newpage()
grid.table(df_Sp)
plotMD(et.Sp.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Sp.F.vs.M$table$logFC, -log10(et.Sp.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Spinal_Cord.Female-1*Spinal_Cord.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Sp.F.vs.M$table[,"PValue"], breaks=50, main="Spinal cord p-value frequency histogram")

# Substantiaa Nigra Female vs Male
et.Sn.F.vs.M <- exactTest(y, c("Substantia_Nigra.Male", "Substantia_Nigra.Female"))
df_Sn <- summary(decideTests(et.Sn.F.vs.M))
grid.newpage()
grid.table(df_Sn)
plotMD(et.Sn.F.vs.M)

# Volcano plot
volcanoData <- cbind(et.Sn.F.vs.M$table$logFC, -log10(et.Sn.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Substantia_Nigra.Female-1*Substantia_Nigra.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Sn.F.vs.M$table[,"PValue"], breaks=50, main="Substantia nigra p-value frequency histogram")

dev.off()



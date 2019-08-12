# This script looks at differential gene expression between males and females within each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata/Age_Matched_Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Count_Matrix.tsv"
FILE_NAME <-  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Gene/Normalized/Exact_Test/Stringtie/Tissue_Exact_AgeMatched.pdf"
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
# OG
# design <- model.matrix(~0+group, data=y$samples) # No intercept
# colnames(design) <- levels(y$samples$group)

Model_Func <- function(fc, df){
  res <- model.matrix(~0 + fc, data = df$samples)
  return(res)
}
Design <- Map(Model_Func, Groups, y)

### idk why this doesn't work; temp hardcode 
# Remove 'fc' from colnames
# Colname_Func <- function(x){
#   colnames(x) <- gsub("fc", "", colnames(x))
# }
# Design <- lapply(Design, Colname_Func)

colnames(Design$Amygdala) <- gsub("fc","",colnames(Design$Amygdala))
colnames(Design$Anterior) <- gsub("fc","",colnames(Design$Anterior))
colnames(Design$Caudate) <- gsub("fc","",colnames(Design$Caudate))
colnames(Design$Cerebellar) <- gsub("fc","",colnames(Design$Cerebellar))
colnames(Design$Cerebellum) <- gsub("fc","",colnames(Design$Cerebellum))
colnames(Design$Cortex) <- gsub("fc","",colnames(Design$Cortex))
colnames(Design$Frontal_Cortex) <- gsub("fc","",colnames(Design$Frontal_Cortex))
colnames(Design$Hippocampus) <- gsub("fc","",colnames(Design$Hippocampus))
colnames(Design$Hypothalamus) <- gsub("fc","",colnames(Design$Hypothalamus))
colnames(Design$Nucleus_Accumbens) <- gsub("fc","",colnames(Design$Nucleus_Accumbens))
colnames(Design$Putamen) <- gsub("fc","",colnames(Design$Putamen))
colnames(Design$Spinal_Cord) <- gsub("fc","",colnames(Design$Spinal_Cord))
colnames(Design$Substantia_Nigra) <- gsub("fc","",colnames(Design$Substantia_Nigra))

# Filter out lowly expressed genes.
# Remove genes w/ <7 counts.
Keep <- lapply(y, function(x){
  rowSums(cpm(x)>1)>=2
})
summary(Keep$Amygdala) # example

Filter_Func <- function(x, k){
  x <- x[k, , keep.lib.sizes=FALSE]
}
y <- Map(Filter_Func, y, Keep)

# TMM Normalization
y <- lapply(y, calcNormFactors)
y$Amygdala$samples

# Estimate common dispersion and tagwise dispersions in one run (recommended)
Dispersion_Func <- function(a, b){
  estimateDisp(a, b, robust=TRUE)
}
y <- Map(Dispersion_Func, y, Design)

# Test for DGX with Exact Test

# Plot
#pdf(FILE_NAME)

# Example of original method
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

# ____________________________________________________________________________
# Comaprisons to test
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
Exact_Func <- function(x, comp){
  exactTest(x, comp)
}
Exact_Res <- Map(Exact_Func, y, Pairs)

# Get summary of results
Summary_Func <- function(x){
  res <- summary(decideTests(x))
  return(res)
}
Results_df <- lapply(Exact_Res, Summary_Func)

Tables <- lapply(Results_df, grid.table)

# Plot summaries of results and MD plots
Plot_Func <- function(a, b){
  print(a)
  grid.newpage()
  grid.table(a)
  plotMD(b)
}
Res_Plots <- Map(Plot_Func, Results_df, Exact_Res) 

# Volcano plot
volcanoData <- cbind(et.Am.F.vs.M$table$logFC, -log10(et.Am.F.vs.M$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala.Female-1*Amygdala.Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

test <- cbind(Exact_Res$Amygdala$table$logFC, -log10(Exact_Res$Amygdala$table[,"PValue"]))

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


# ____________________________________________________________________________


#dev.off()



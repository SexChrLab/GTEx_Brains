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
summary(Keep$Amygdala) # example

Filter_Func <- function(x, k){
  x <- x[k, , keep.lib.sizes=FALSE]
}
y_Filtered <- Map(Filter_Func, y, Keep)

# TMM Normalization
y <- lapply(y, calcNormFactors)
View(y$Amygdala$samples)

y_Filtered <- lapply(y_Filtered, calcNormFactors)
View(y_Filtered$Amygdala$samples)

# Estimate common dispersion and tagwise dispersions in one run (recommended)
Dispersion_Func <- function(a, b){
  estimateDisp(a, b, robust=TRUE)
}
y <- Map(Dispersion_Func, y, Design)
y_Filtered <- Map(Dispersion_Func, y, Design)

# Test for DGX with Exact Test

# Plot
#pdf(FILE_NAME)

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
Exact_Res_Filtered <-  Map(Exact_Func, y_Filtered, Pairs)

# Get summary of results
Summary_Func <- function(x){
  res <- summary(decideTests(x))
  return(res)
}
Results_df <- lapply(Exact_Res, Summary_Func)
Results_df_Filtered <- lapply(Exact_Res_Filtered, Summary_Func)

# Add column to DGX results w/ and w/out filtering CPM < 1
# Does not make a difference...
# library(data.table)
# l <- list(Results_df, Results_df_Filtered)
# test <- rbindlist(l)

# Plots
Titles <- list('Amygdala', 'Anterior', 'Caudate', 'Cerbellum', 'Cerebellar', 'Cortex', 'Frontal Cortex',
               'Hippocampus', 'Hypothalamus', 'Nucleus Accumbens', 'Putamen', 'Spinal Cord', 'Substantia Nigra')

# Temporary: Prints plots followed by summary
Plot_Func <- function(a, b){
  print(a)
  grid.newpage()
  grid.table(a)
  plotMD(b)
}
TEMP_Res_Plots <- Map(Plot_Func, Results_df, Exact_Res)

# Plot MD plots on one page
opar <- par(no.readonly = TRUE) 
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
MD_Plot_Func <- function(x){
  plotMD(x)
  mtext('Mean-Difference Plots', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('log-fold-change', side = 2, outer = TRUE, line=2)
}
Res_Plots <- lapply(Exact_Res, MD_Plot_Func)
legend(legend=c("Up","Not Sig", "Down"), pch = 16, col = c("blue","black", "green"), bty = "n", xpd=NA)
#par(opar) # reset par

legend(1.75, 5.5, inset=0, legend=levels(Meta$Cortex$Sex), pch=16, cex=2.0, col=colors, xpd=NA)

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
  mtext('Hisat: Volcano Plots', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
Map(Plot_Func, Volcano_Res, Titles, Subset_Res)
legend(16.0, 9.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"), 
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)





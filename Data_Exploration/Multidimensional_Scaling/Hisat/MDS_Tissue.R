# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold changes between each pair of samples.
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Multidimensional_Scaling/Hisat")

METADATA = file.path('/scratch/mjpete11/GTEx/Metadata/', 'Matched_Metadata.csv')
TABLE_DIR =  '/scratch/mjpete11/GTEx/Data_Exploration/Sex_Tissue_Age/Multidimensional_Scaling/Hisat/Matched_Table'
MDS_K2 = '/scratch/mjpete11/GTEx/Data_Exploration/Multidimensional_Scaling/Hisat/Matched_Plots/Matched_Hisat_MDS_k2.pdf'
MDS_K4 = '/scratch/mjpete11/GTEx/Data_Exploration/Multidimensional_Scaling/Hisat/Matched_Plots/Matched_Hisat_MDS_k4.pdf'
REPS_K2 = '/scratch/mjpete11/GTEx/Data_Exploration/Multidimensional_Scaling/Hisat/Matched_Plots/Matched_Reps_Hisat_MDS_k2.pdf'
REPS_K4 = '/scratch/mjpete11/GTEx/Data_Exploration/Multidimensional_Scaling/Hisat/Matched_Plots/Matched_Reps_Hisat_MDS_k4.pdf'

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
library(readr)
library(stringr)
library(dplyr)
library(limma)
library(edgeR)

# Read Metadata CSV.                                                            
Samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample                                             

# Load count data 
cts <- lapply(PATHS, function(x) {
  t <- as.matrix(read.csv(x, row.names="gene_id"))
})
names(cts) <- c('Amygdala', 'Anterior', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_Cortex',
                'Hippocampus', 'Hypothalamus', 'Nucleus_Accumbens', 'Putamen', 'Spinal_Cord', 'Substantia_Nigra')

# Replace . to - in colnames in each df
for (i in seq_along(cts)){
  colnames(cts[[i]]) <- str_replace_all(colnames(cts[[i]]), pattern = "\\.","-")
}

# Remove the - at the end of the Caudate sample names
colnames(cts$Caudate) <- substring(colnames(cts$Caudate), 1, nchar(colnames(cts$Caudate))-1)

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Remove samples not present in metadata
Subset_Func <- function(df.1, df.2) {
  df.1 <- df.1[, intersect(colnames(df.1), as.character(df.2$Sample))]
  return(df.1)
}
Tissue_Count <-  Map(Subset_Func, cts, Meta)

# Check that rownames equals colnames 
Check_Func <- function(a, b){
  Check <- all(rownames(a) %in% colnames(b))
}
Res_1 <- Map(Check_Func, a=Meta, b=Tissue_Count)

# Match order of samples in metadata and count data
Match_Order <- function(a, b){ # Will throw an error if dfs are different lengths
  a <- a[, rownames(b)]
}
#Tissue_Count <- Map(Match_Order, a=Tissue_Count, b=Meta)
test <- Tissue_Count <- Map(Match_Order, a=Tissue_Count, b=Meta)

# Check that order matches
Match_Check <- function(a, b){
  Match <- all(rownames(a) == colnames(b))
}
Res_2 <- Map(Match_Check, a=Meta, b=Tissue_Count)

# Create DGEList object for each tissue count matrix
DGE_lst <- lapply(Tissue_Count, function(x){
  DGEList(x)
}) 

# Titles and counts without replicates
DGE_Without_Reps <- DGE_lst[-c(4,6)]
Meta_Without_Reps <- Meta[-c(4,6)]

# MDS plots on one page
colors <- c("blue", "darkgreen")
# Dim 1 & 2
MDS_FUN_k2 <- function(DGE, NAME, META, TOP) {
  k2_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = TOP, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(1,2), 
                    col = colors[META[['Sex']]],
                    main = NAME)
  mtext('Hisat MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes', side=3, outer=TRUE, line=3)
  mtext('Dimension 1', side = 1, outer = TRUE, line=1)
  mtext('Dimension 2', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
pdf(MDS_K2)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
k2_MDS <- Map(MDS_FUN_k2, DGE = DGE_Without_Reps, NAME = names(DGE_Without_Reps), META = Meta_Without_Reps, TOP = 100)
legend(5.0, 2.5, inset=0, legend=levels(Meta$Amygdala$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()

# Replicates only
DGE_Reps <- DGE_lst[c(4:7)]
Meta_Reps <- Meta[c(4:7)]

pdf(MDS_K2)
par(mfrow = c(2, 2), cex=0.4, mar = c(3, 2, 2, 6), oma =c(5, 5, 6, 6), xpd=TRUE)
k2_MDS <- Map(MDS_FUN_k2, DGE = DGE_Reps, NAME = names(DGE_Reps), META = Meta_Reps, TOP = 100)
legend(5.0, 2.5, inset=0, legend=levels(Meta$Amygdala$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()

# Dim 3 & 4
MDS_FUN_k4 <- function(DGE, NAME, META, TOP) {
  k4_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = TOP, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(3,4), 
                    col = colors[META$Sex],
                    main = NAME)
  mtext('Hisat MDS Plots: Dimensions 3 and 4; Top 100 Most Variable Genes', side=3, outer=TRUE, line=3)
  mtext('Dimension 3', side = 1, outer = TRUE, line=1)
  mtext('Dimension 4', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
pdf(MDS_K4)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
k4_MDS <- Map(MDS_FUN_k4, DGE = DGE_Without_Reps, NAME = names(DGE_Without_Reps), META = Meta_Without_Reps, TOP = 100)
legend(2.0, 1.5, inset=0, legend=levels(Meta$Amygdala$Sex), pch=16, cex=2.0, col=colors, xpd=NA) 
dev.off()

# Replicates only
pdf(REPS_K4)
par(mfrow = c(2, 2), cex=0.4, mar = c(3, 2, 2, 6), oma =c(5, 5, 6, 6), xpd=TRUE)
k4_MDS <- Map(MDS_FUN_k4, DGE = DGE_Reps, NAME = names(DGE_Reps), META = Meta_Reps, TOP = 100)
legend(1.5, 2.5, inset=0, legend=levels(Meta$Cortex$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()

# Eigenvector table
setwd(TABLE_DIR)
sapply(names(k2_MDS), 
       function (x) write.table(k2_MDS[[x]]['cmdscale.out'], file=paste(x, "txt", sep=".")))



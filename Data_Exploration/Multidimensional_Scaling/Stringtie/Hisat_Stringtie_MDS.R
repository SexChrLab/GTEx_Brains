# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-foldchanges between each pair of samples.

METADATA = file.path('/scratch/mjpete11/GTEx/Metadata/', 'Metadata.csv')
PLOT_DIR = '/scratch/mjpete11/GTEx/Data_Exploration/Sex_Tissue_Age/Multidimensional_Scaling/Stringtie_Plots/'
TABLE_DIR =  '/scratch/mjpete11/GTEx/Data_Exploration/Sex_Tissue_Age/Multidimensional_Scaling/Stringtie_Tables/'
PLOT_FILE = 'MDS_Plots.pdf'

# Load packages                                                                 
library(readr)
library(stringr)
library(dplyr)
library(limma)
library(edgeR)
library(DESeq2)

# Read Metadata CSV.                                                            
Samples <- read.csv(METADATA, sep=",", header=TRUE)

# Load gene(/transcript) count matrix and labels
# Function to load count data
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


cts <- lapply(PATHS, function(x) {
  t <- as.matrix(read.csv(x, row.names="gene_id"))
})
names(cts) <- c('Amygdala', 'Anterior', 'Caudate', 'Cerebellar', 'Cerebellum', 'Cortex', 'Frontal_Cortex',
                'Hippocampus', 'Hypothalamus', 'Nucleus_Accumbens', 'Putamen', 'Spinal_Cord', 'Substantia_Nigra')

# Remove the . at the end of the Caudate sample names
# No idea why its just the Caudate samples...
colnames(cts$Caudate) <- substring(colnames(cts$Caudate), 1, nchar(colnames(cts$Caudate))-1)

# Replace . to - in colnames in each df
for (i in seq_along(cts)){
  colnames(cts[[i]]) <- str_replace_all(colnames(cts[[i]]), pattern = "\\.","-")
}

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample            

# List of sample metadata for each tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

### TEMPORARY ###
# Drop sample from amygdala metadata that is missing count data.
Meta$Amygdala <- Meta$Amygdala[-c(53),]
Meta$Substantia_Nigra <- Meta$Substantia_Nigra[-c(16),]

# Check that samples in count data are listed in metadata
all(rownames(Meta$Amygdala) %in% colnames(cts$Amygdala))

Check_Func <- function(a, b){
  Check <- all(rownames(a) %in% colnames(b))
}
Res <- Map(Check_Func, a=Meta, b=cts)

# Match order of samples in metadata and count data
Match_Order <- function(a, b){ # Will throw an error if dfs are different lengths
  a <- a[, rownames(b)]
}

cts <- Map(Match_Order, a=cts, b=Meta)

# Check that order matches
Match_Check <- function(a, b){
  Match <- all(rownames(a) == colnames(b))
}
Res_2 <- Map(Match_Check, a=Meta, b=cts)

# Create DGEList object for each tissue count matrix
DGE_lst <- lapply(cts, function(x){
  DGEList(x)
})

# MDS plots
#setwd(PLOT_DIR)
#pdf(PLOT_FILE)

# Plot constants 
colors <- c("blue", "darkgreen")
par(mar=c(8, 4.1, 4.1, 8), xpd=TRUE) # margins: bottom, left, top, and right

# MDS plot function for dimensions 1 and 2
MDS_FUN <- function(DGE, NAME, META, TOP, NUM_GENES) {
  k2_MDS <- plotMDS(DGE,
                   gene.selection = "common",
                   top = TOP, 
                   pch = 16, 
                   cex = 1, 
                   dim.plot = c(1,2), 
                   col = colors[META$Sex],
                   main = paste('MDS Plot: Dim 1 and 2; ', NUM_GENES,'; ', NAME))
  legend("topright", inset=c(-0.3,0), legend=levels(META$Sex), pch=16, cex=1,col=colors, xpd=TRUE) # legend will be out of range in viewer
  return(k2_MDS)
}
k2_MDS <- Map(MDS_FUN, DGE = DGE_lst, NAME = names(DGE_lst), META = Meta, TOP = 100, NUM_GENES='Top 100')

# MDS plot function for dimensions 3 and 4
MDS_FUN <- function(DGE, NAME, META, TOP, NUM_GENES) {
  k4_MDS <- plotMDS(DGE,
                   gene.selection = "common",
                   top = TOP, 
                   pch = 16, 
                   cex = 1, 
                   dim.plot = c(3,4), 
                   col = colors[META$Sex],
                   main = paste('MDS Plot: Dim 3 and 4; ', NUM_GENES,'; ', NAME))
  legend("topright", inset=c(-0.3,0), legend=levels(Samples$Sex), pch=16, cex=1,col=colors, xpd=TRUE)
  return(k4_MDS)
}
k4_MDS <- Map(MDS_FUN, DGE = DGE_lst, NAME = names(DGE_lst), META = Meta, TOP = 100, NUM_GENES='Top 100')


# Plot MDS plots on one page

# Dim 1 & 2
# Titles and counts without replicates
DGE_Without_Reps <- DGE_lst[-c(4,6)]
Meta_Without_Reps <- Meta[-c(4,6)]

pdf('Hisat_MDS_k2_One_Page.pdf')
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
MDS_FUN_k2 <- function(DGE, NAME, META, TOP) {
  k2_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = TOP, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(1,2), 
                    col = colors[META$Sex],
                    main = NAME)
  mtext('Hisat MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes', side=3, outer=TRUE, line=3)
  mtext('Principle Component 1', side = 1, outer = TRUE, line=1)
  mtext('Principle Component 2', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
k2_MDS <- Map(MDS_FUN_k2, DGE = DGE_Without_Reps, NAME = names(DGE_Without_Reps), META = Meta_Without_Reps, TOP = 100)
legend(3.0, 2.5, inset=0, legend=levels(Meta$Amygdala$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()

# Dim 3 & 4
pdf('Hisat_MDS_k4_One_Page.pdf')
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
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
  mtext('Principle Component 3', side = 1, outer = TRUE, line=1)
  mtext('Principle Component 4', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
k4_MDS <- Map(MDS_FUN_k4, DGE = DGE_Without_Reps, NAME = names(DGE_Without_Reps), META = Meta_Without_Reps, TOP = 100)
legend(3.5, 3.0, inset=0, legend=levels(Meta$Amygdala$Sex), pch=16, cex=2.0, col=colors, xpd=NA) 
dev.off()

# Plot MDS plots on one page

# Dim 1 & 2
# Repplicates only
DGE_Reps <- DGE_lst[c(4:7)]
Meta_Reps <- Meta[c(4:7)]

pdf('Reps_Hisat_MDS_k2_One_Page.pdf')
par(mfrow = c(2, 2), cex=0.4, mar = c(3, 2, 2, 6), oma =c(5, 5, 6, 6), xpd=TRUE) # margins: c(bottom, left, top, right)
MDS_FUN_k2 <- function(DGE, NAME, META, TOP) {
  k2_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = TOP, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(1,2), 
                    col = colors[META$Sex],
                    main = NAME)
  mtext('Hisat MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes', side=3, outer=TRUE, line=3)
  mtext('Principle Component 1', side = 1, outer = TRUE, line=1)
  mtext('Principle Component 2', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
k2_MDS <- Map(MDS_FUN_k2, DGE = DGE_Reps, NAME = names(DGE_Reps), META = Meta_Reps, TOP = 100)
legend(1.75, 5.5, inset=0, legend=levels(Meta$Cortex$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()

# Dim 3 & 4
pdf('Reps_Hisat_MDS_k4_One_Page.pdf')
par(mfrow = c(2, 2), cex=0.4, mar = c(3, 2, 2, 6), oma =c(5, 5, 6, 6), xpd=TRUE) # margins: c(bottom, left, top, right)
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
  mtext('Principle Component 3', side = 1, outer = TRUE, line=1)
  mtext('Principle Component 4', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
k4_MDS <- Map(MDS_FUN_k4, DGE = DGE_Reps, NAME = names(DGE_Reps), META = Meta_Reps, TOP = 100)
legend(2.5, 3.0, inset=0, legend=levels(Meta$Cortex$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()


# Eigenvector table
setwd(TABLE_DIR)
sapply(names(k2_MDS), 
       function (x) write.table(k2_MDS[[x]]['cmdscale.out'], file=paste(x, "txt", sep=".")))

k2_table <- lapply(k2_MDS, function(x){
  k2_MDS[[x]]['cmdscale.out']
})
       
write.csv(rbind(df1, d32, df3), "filename.csv")





#!/usr/bin/env Rscript

# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold changes between each pair of samples.

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain/data"
METADATA <- file.path(BASE, "metadata.csv")
COUNTS <- file.path(BASE, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
MDS <- file.path(BASE, "MDS.pdf"

# Load packages
library(data.table)
library(readr)
library(stringr)
library(dplyr)
library(limma)
library(edgeR)

# Read in files                                                            
Meta <- read.table(METADATA, header = TRUE, sep=",")
Counts <- fread(COUNTS) 

# Quick view
head(Meta)
head(colnames(Counts))

# Set rownames of metadata object equal to sample names
rownames(Meta) <- Meta$Sample_ID

# Make list of lists of samples for each tissue
Tissue_Lst <- list()
for(i in 1:length(levels(Meta$Tissue))){
		  Tissue_Lst[[i]] <- as.vector(Meta[Meta$Tissue == levels(Meta$Tissue)[i], "Sample_ID"])
}

# Rename lists in list as tissue names
names(Tissue_Lst) <- levels(Meta$Tissue)
str(Tissue_Lst)

# Split counts into list of dfs by tissue
Tissue_Count <- list()
for(i in 1:length(Tissue_Lst)){
		  Tissue_Count[[i]] <- Counts[,which(colnames(Counts) %in% Tissue_Lst[[i]])]
}
names(Tissue_Count) <- levels(Meta$Tissue)
str(Tissue_Count)

# Metadata split into list of dfs by tissue
Meta_Lst <- list()
for(i in 1:length(levels(Meta$Tissue))){
		  Meta_Lst[[i]] <- Meta[Meta$Tissue == levels(Meta$Tissue)[i],]
}
names(Meta_Lst) <- levels(Meta$Tissue)

# Check that rownames equals colnames 
Check <- function(a, b){
	all(rownames(a) %in% colnames(b))
}
Res_1 <- Map(Check, a=Meta, b=Tissue_Count)
all(Res_1==TRUE)

# Check that order matches
Match_Check <- function(a, b){
	 Match <- all(rownames(a) == colnames(b))
}
Res_2 <- Map(Match_Check, a=Meta, b=Tissue_Count)
all(Res_2==TRUE)

# Create DGEList object for each tissue count matrix
DGE_lst <- lapply(Tissue_Count, function(x){ DGEList(x)}) 

# MDS plots on one page
colors <- c("blue", "darkgreen")

MDS_FUN_k2 <- function(DGE, NAME, META, TITLE) {
  k2_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = 100, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(1,2), 
                    col = colors[META[['Sex']]],
                    main = NAME)
  mtext(TITLE, side=3, outer=TRUE, line=3)
  mtext('Dimension 1', side = 1, outer = TRUE, line=1)
  mtext('Dimension 2', side = 2, outer = TRUE, line=2)
  return(k2_MDS)
}
pdf(MDS)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) # margins: c(bottom, left, top, right)
k2_MDS <- Map(MDS_FUN_k2, 
			  DGE = DGE_Lst, 
			  NAME = names(DGE_Lst), 
			  META = Meta_Lst, 
			  TITLE='GTEx v8 MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes')
legend(5.0, 2.5, inset=0, legend=levels(Meta$Amygdala$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()
